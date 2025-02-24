program WSME_genTden_loopy
  !using the hamoltonian:
  !     H=- T  \sum_i q + eps*ASAunfoprotsa + sum_{i<=j} (-eps Delta_ij x_{i,j}) +\sum_i (T q+alfa c) x_{i,i})
  ! where Delta_ij represent the contribution of contact ij to the DeltaASA(U-N).
  ! Notice that i<=j also contribution from an isolated native residue
  use defreal
  use protdep_par
  use globalpar
  use phys_const
  use calc_interact
  use thermo
  use thermoab
  use disulfide_module
  implicit none
  real(kind=db)::parv(nparmax) !vector of npar parameters on which the energy depends (i.e. q, eps, ...)
  real(kind=db):: Phi(3),natbase(3), Phi_garbage(3), natbase_garbage(3) !garbages are used to get the flipped map w/out breaking anything. Useless output.
  real(kind=db), allocatable :: e(:,:,:),F(:),S(:,:),m(:,:),nu(:,:),auxe(:,:,:)
  real(kind=db):: FreeonRT,EonRT,EnthonRT,logZeta,M0,Minf,ConR,Mavg,sigmaavg !thermodinamic variables
  real(kind=db):: EonRTsquaredaux,ConR_fixedconfaux !PIER: added this 27/8/24

  real(kind=db):: T,cden
  real(kind=db):: fracfold
  real(kind=db), allocatable:: logZetaab(:,:), EonRTab(:,:), ConRab(:,:), sigmaab(:,:),sigmaiab(:,:,:) !matrix to allocate loop influence to thermodinamic variables of folded islands (a,b)
  real(kind=db), allocatable::  EonRTabsquared(:,:),ConRab_fixedconf(:,:),EonRTabtot(:,:),&
       & EonRTabsquaredtot(:,:),ConRab_fixedconftot(:,:),Unsquared(:,:) !PIER: added this 27/8/24  

  real(kind=db), allocatable:: Hn(:,:), Un(:,:),Cvn(:,:) !matrix to allocate the termodinamic variables of islands (a,b) of the all-native case                 
               
  real(kind=db), allocatable::Mi(:),sigmai(:),Mij(:,:),sigmaij(:,:) !arrays to allocate the magnetizations per residue and per island					        								 
  real(kind=db):: logZetaaux, EonRTaux, ConRaux, sigmaaux,Maux,sig,aux !auxilliar variables 		                 					 
  real(kind=db), allocatable ::              sigmaiaux(:)
  real(kind=db), allocatable :: sigma_st_ab_matrix(:,:,:), sigma_st_ab_aux(:), sigma_st(:)
  real(kind=db), allocatable :: sigma_st_ab_all_matrix(:,:,:,:), sigma_st_ab_all_aux(:,:), sigma_st_all(:,:) !first two dimensions are all possibilities of (s,t) intervals
  real(kind=db), allocatable :: L(:),fun_L(:),sigma_nu(:,:) !summatory on all (s,t) islands in function of L=1:N
  real(kind=db):: ConRtest,ConRfixave !PIER: checks 27/8/24

  real(kind=db),allocatable :: deltaT_array(:), t_exp_aux(:), c_exp_aux(:) !used while optimization

  real(kind=db), allocatable :: folded_ab_ij_matrix(:,:,:,:), folded_ab_matrix(:,:), fold_nu(:,:),&
                                log_ZETA_aux(:),  log_ZETA_abj(:,:,:) !log_ZETA_abj is log of Z_{a,j}^{a,b}/Z_{a,j-1}^{a,b}
                                                         
  integer:: i,j,k,iM,l_index,q,s_index,t_index,g,aux2,aux3									 
  integer:: aaux,baux !variables to go through islands								 
  !  real(kind=db):: nu(1:N,1:N)

  integer, allocatable :: SS_matrix(:,:) ! matrix with the disulfide bonds
  integer :: num_rows, num_cols, ss_index, ss_a, ss_b! auxiliar variables to build the matrix and exclude states from Z
  logical :: ss_ok_do_calculation

  !we give dimension to matrix
  read (*,*) N
  read (*,*) N1
  allocate(delta(5,1:N,1:N)) 
  allocate(rCalpha(1:N,1:N))
  allocate(v(N1,5)) 
  allocate(e(4,N,0:N),F(0:N),S(0:N+1,0:N+1),m(0:N+1,0:N+1),nu(1:N,1:N),logZetaab(1:N,1:N), EonRTab(1:N,1:N),&
       & ConRab(1:N,1:N), sigmaab(1:N,1:N),Hn(1:N,1:N), Un(1:N,1:N),Cvn(1:N,1:N),Mi(1:N),sigmai(1:N),sigmaij(1:N,1:N),Mij(1:N,1:N),&
       & sigmaiaux(1:N),sigmaiab(1:N,1:N,1:N)) 		
  allocate(EonRTabsquared(1:N,1:N),ConRab_fixedconf(1:N,1:N)) !PIER: added this 27/8/24
  allocate(EonRTabtot(1:N,1:N),EonRTabsquaredtot(1:N,1:N),ConRab_fixedconftot(1:N,1:N),Unsquared(1:N,1:N)) !PIER: added this 27/8/24

  call read_init(&              ! No imput because we pass the file with the console in format ____.exe < input.in
                                !     O:
       &     parv)

  allocate(sigma_st_ab_matrix(1:ST_length,1:N,1:N),sigma_st_ab_aux(1:ST_length),sigma_st(1:ST_length))
  allocate(sigma_st_ab_all_matrix(1:N,1:N,1:N,1:N), sigma_st_ab_all_aux(1:N,1:N), sigma_st_all(1:N,1:N))
  allocate(L(1:N),fun_L(1:N),sigma_nu(N,N),fold_nu(N,N))
  allocate(t_exp_aux(1:nexp),c_exp_aux(1:nexp))
  allocate(folded_ab_ij_matrix(N,N,N,N),folded_ab_matrix(N,N),log_ZETA_aux(N),log_ZETA_abj(N,N,N))

  !if only the specific heat is requested, the input flags are redefined accordingly:
  if(onlyC) then
     wEave=.false.
     wC=.true.
     wfoldfr=.false.
     wFprof=.false.
     wMave=.false.
     wMres=.false.
     wMisland=.false.
     wmprof=.false.
     wstr=.false.
     wProd_ms=.false.
     f_L=.false.
     fold_profile=.false.
  endif
  ! Mavg(T=0K)=1:
  M0=1._db
  ! calculate Minf=Mavg(T=481K)
  Minf=0._db

  open(20,file='Output/profthermo.dat')
  !write(20,*) "#cden        T        (Mavg-Minf)/(M0-Minf)        sigmaavg       Free energy        Enthalpy        R*(EnthonRT-FreeonRT)      Cp     R*Phi(3)      R*natbase(3)" !error with line length of compilator
  open(50,file='Output/magnet.dat')
  write(50,*) "#T       Residue #i       m_i         sigma_i         m_i - sigma_i "
  open(60,file='Output/ener.dat')
  open(94,file="Output/simCp.txt")
  if(wFprof) open(30,file='Output/Fprof.dat') !profiles
  if(wmprof) open(35,file='Output/mprofile.dat') !profiles
  if(wstr) open(40,file='Output/strings.dat') !native strings
  if(wProd_ms) then
      open(70,file="Output/sigma_st.dat") !s t <prod_k=s^t m_k sigma_k>
      write(70,*) "#[Cden]    T   <prod_k=s^t m_k sigma_k> for all (s,t) intervals"
  end if
  if(f_L) then
   open(73,file="Output/f_L.txt") !f(L) = 1/(N-L+1) * sum_i^N-L+1 [ <prod_k=i^i+L-1 m_k sigma_k> ] 
   open(74,file="Output/sigma_profile.txt")
   open(75,file="Output/sigma_profile cont.txt")
  endif
  if(fold_profile) then
   open(47,file="Output/fold_profile.txt") !<prod 1-sigma><prod m>
   open(48,file="Output/f_L fold.txt")
   open(49,file="Output/fold_profile cont.txt")
  endif

  ! Disulfide bridges
  if(SS_breakable .or. SS_flag)  call get_disulfide_bonds_matrix(pdb_code, SS_matrix, num_rows, num_cols) ! SS_matrix: each row is a bond, the two columns have the two residues which conformates it
  
  allocate(auxe(4,1:N,0:N))

  !This section of code if used if we want to only simulate temperatures we have experimental data for
  if(constant_deltaT) then
   allocate(deltaT_array(1:1000))
      deltaT_array=deltaT
  else
      allocate(deltaT_array(1:nexp))
      open(56,file="Input/"//trim(expfile))
      do i=1,nexp
         read(56,*) t_exp_aux(i),c_exp_aux(i)
      enddo
      do i=1,nexp-1
         deltaT_array(i)=t_exp_aux(i+1)-t_exp_aux(i)
      enddo
      deltaT_array(nexp)=1 !just to get out of the while bucle
      Tmin=t_exp_aux(1)
      Tmax=t_exp_aux(size(t_exp_aux))
  endif 

  cden=cdenmin ! concentration of denaturant. We only use 0
  do while (cden.le.cdenmax)
     q=1 
     T=Tmin
     do while (T <= Tmax)
        FreeonRT=0._db
        logZeta=0._db
        EonRT=0._db
        ConR=0._db

        call calc_e_Phi(& ! calculates interaction between residues without taking into account their folding state (conceptually: h_ij)
                                !     I:
             &          T,cden,parv,SS_matrix,&
                                !     O:
             &          e,Phi,natbase)
        !     natbase,phi (1,2,3)= completely native and unfolded baselines for FRT,H/RT,C/R
        auxe=e !probably useless

        !loop to callculate the partition fuction and observables inside a folded island (loops influence)
        logZetaab=0._db
        EonRTab=0._db
        ConRab=0._db
        sigmaab=0._db
        sigmaiab=0._db
        sigma_st_ab_matrix=0._db
        sigma_st=0._db
        sigma_st_ab_all_aux=0._db
        sigma_st_ab_all_matrix=0._db
        sigma_st_all=0._db
        EonRTabsquared=0._db !PIER: added this and below (4 lines) 27/8/24
        ConRab_fixedconf=0._db
        EonRTabsquaredtot=0._db
        ConRab_fixedconftot=0._db
        log_ZETA_aux=0._db
        log_ZETA_abj=0._db
        folded_ab_ij_matrix=0._db
        folded_ab_matrix=0._db
         do aaux=1,N
            do baux=aaux,N

                  call calc_thermoab(aaux,baux-aaux+1,auxe,& !calculates the contributions of each (a->b) island of m=1,s=0 in a sea of m=1,s=1
                        logZetaaux,EonRTaux,EonRTsquaredaux,ConR_fixedconfaux,ConRaux,sigmaaux,sigmaiaux,fracfold,&
                        sigma_st_ab_aux,& !sigmaaux is <s>_(a,b). sigmai(aux) is <s_i>_ab = sum_(a,b) f^i_(a,b)*Z^(a,b) where f=1 if a in (a,b) and 0 if not
                        sigma_st_ab_all_aux,parv(8),log_ZETA_aux,SS_matrix) ! <prod_k=s_t m_k sigma_k> for each (s,t) island (only a->b contribution)
                        
      !!$           !PIER: To build the pure WSME limit, without loops, comment lines below: from here...
                  logZetaab(aaux,baux)=logZetaaux
                  EonRTab(aaux,baux)=EonRTaux
                  EonRTabsquared(aaux,baux)=EonRTsquaredaux !PIER: added this  27/8/24
                  ConRab_fixedconf(aaux,baux)=ConR_fixedconfaux !PIER: added this  27/8/24
                  sigmaab(aaux,baux)=sigmaaux
                  do k=1,N
                     sigmaiab(k,aaux,baux)=sigmaiaux(k)
                  enddo

                  do k=1,ST_length
                        sigma_st_ab_matrix(k,aaux,baux)=sigma_st_ab_aux(k) ! _aux is an array with ST_interval columns, and a value for one (a,b) interval for each (S,T)
                  end do
                  
                  do s_index=1,N
                     do t_index=s_index,N
                           sigma_st_ab_all_matrix(s_index,t_index,aaux,baux)=sigma_st_ab_all_aux(s_index,t_index)
                     enddo 
                  enddo

                  if(fold_profile) then
                     !save ZETA
                     do k=1,N
                        log_ZETA_abj(aaux,baux,k)=log_ZETA_aux(k) !imagine log_ZETA_abj as a matrix, each row is an (a,b) island and j-th column is ZETA_j (Z_1^(ab), Z2, Z3... 0 0 0 0...0). There are a-b numbers /= 0, and index numbers are relative to "a" (Z1=Za, Z2=Z_a+1...)
                     enddo
                     do j=aaux,baux
                        do i=j,aaux,-1 !i goes from j to a (j, j-1, j-2...a)
                           if (i==j) then
                              folded_ab_ij_matrix(aaux,baux,i,j)=1-sigmaiab(i,aaux,baux) !remember sigmaiab is <sigma_k>^(a,b)
                           else
                              folded_ab_ij_matrix(aaux,baux,i,j)=exp(log_ZETA_abj(aaux,baux,i-aaux+1))*&
                                       folded_ab_ij_matrix(aaux,baux,i+1,j) ! in log_ZETA_abj(aaux,baux,i-a+1), the "-a+1" is an offset
                           endif
                           !if ( folded_ab_ij_matrix(aaux,baux,i,j) < 0 .or.  folded_ab_ij_matrix(aaux,baux,i,j)>1) &
                           !write(*,*) "F(",aaux,",",baux,",",i,",",j,")=",folded_ab_ij_matrix(aaux,baux,i,j)
                        enddo
                     enddo
                  endif
   !!$        !...to here.
            enddo
         enddo


         !we have to add the off-set of the all-native case:
         !Hn is betaH^nn
         Hn=0.0_db
         Un=0.0_db
         Unsquared=0.0_db !PIER: added this  27/8/24
         Cvn=0.0_db
         do j=1,N
               do i=1,j
               do k=i,j
                  do g=k,j
                     !auxe is just e(), the matrix calculated with calc_e_Phi_genTden_loopy
                     Hn(i,j)=Hn(i,j)+auxe(1,k,g) ! Hn(i,j)=(effective energy of the i,j native island)/RT (being 0 the completely unfolded energy as reference)
                     Un(i,j)=Un(i,j)+auxe(2,k,g)
                     Cvn(i,j)=Cvn(i,j)+auxe(3,k,g)
                     if (k.eq.g) then
                        Hn(i,j)=Hn(i,j)- auxe(4,k,k-1) !PIER: added this 27/08/2024
                     endif
                  enddo
               enddo
               Unsquared(i,j)=Un(i,j)**2

            enddo
         enddo 

        !PIER: added this part for Esquared and C_fixedconf
        do j=1,N !PIER: added this do loop  27/8/24
           do i=1,j
              EonRTabsquaredtot(i,j)=EonRTabsquared(i,j)+2*EonRTab(i,j)*Un(i,j)+Unsquared(i,j)
           enddo
        enddo
        logZetaab=logZetaab- Hn !PIER: changed this sign 27/08/2024
        EonRTabtot=EonRTab+Un
        ConRab_fixedconftot=ConRab_fixedconf+Cvn !PIER: added this  27/8/24

        !finally we calculate partition function and observables:
        !notice calc_thermoab was in bucle of (a->b) while here the argument is just the number of residues N bc the bucle is within the subroutine
        call calc_thermo2(N, logZetaab, EonRTabtot,EonRTabsquaredtot, ConRab_fixedconftot, sigmaab,sigmaiab,sigma_st_ab_matrix,&
        &      sigma_st_ab_all_matrix, logZeta, EonRT, ConR,ConRfixave, Mavg, sigmaavg,fracfold,Mi,sigmai,Mij,sigmaij,&
               sigma_st,sigma_st_all,folded_ab_ij_matrix,folded_ab_matrix,SS_matrix)
        !RESULTS:

        FreeonRT=-logZeta+Phi(1)
        EnthonRT=EonRT+Phi(2)
        ConR=ConR+Phi(3)

        write(20,*) cden,T,(Mavg-Minf)/(M0-Minf),sigmaavg,&
             &       R*T*FreeonRT,R*T*EnthonRT,R*(EnthonRT-FreeonRT),R*ConR,R*Phi(3),R*natbase(3)

        write(94,*) R*ConR ! specific heat
        if (wstr) then
           call stringhe(m,S,e,N)
           do i=1,N
              do j=i,N


                 write(40,*) cden,T,i,j,m(i,j),S(i-1,j+1)
                 !     m(i,j): island of  1s between  i and j; S(i-1,j+1): island of 1s between  i and j, with 0 at i-1, j+1,
              enddo
           enddo
        endif

        if(wProd_ms) then
            write(70,*) Cden, T, sigma_st ! remember sigma_st is an array with the value of <prod_k=s^t m_k sigma_k> for all (s,t) intervals
        end if

        if(f_L) then
            fun_L=0._db

            do l_index=1,N
               do i=1,N-l_index+1  !f(L)~ sum_i^N-L+1 [ <prod_k=i^i+L-1 m_k sigma_k> ] where <prod_k=i^i+L-1 m_k sigma_k> is sigma_st_all(s,t). No bucle for k is needed.
                  j=i+l_index-1
                  if(i/=1 .and. j/=N) then
                     fun_L(l_index)=fun_L(l_index)+&
                     & ( sigma_st_all(i,j) - sigma_st_all(i,j+1) - &
                     &   sigma_st_all(i-1,j) + sigma_st_all(i-1,j+1))
                  endif
                  if(i/=1 .and. j==N) then
                     fun_L(l_index)=fun_L(l_index)+&
                      & ( sigma_st_all(i,j) -  sigma_st_all(i-1,j)) 
                  endif
                  if(i==1 .and. j/=N) then
                     fun_L(l_index)=fun_L(l_index)+&! max(,0) to avoid case L=i=j=1, <ij> = 0 and <i,j+1> = 1
                     & max(sigma_st_all(i,j) - sigma_st_all(i,j+1),0.0_db)
                  endif
                  if(i==1 .and. j==N) fun_L(l_index)=fun_L(l_index)+sigma_st_all(i,j)
               enddo !i
               fun_L(l_index)=fun_L(l_index)/(N-l_index+1) !normalization
            enddo

            do i=1,N
               do j=i,N
                  if(i/=1 .and. j/=N) then
                     sigma_nu(i,j)= &
                     & ( sigma_st_all(i,j) - sigma_st_all(i,j+1) - &
                     &   sigma_st_all(i-1,j) + sigma_st_all(i-1,j+1))
                  endif
                  if(i/=1 .and. j==N) then
                     sigma_nu(i,j)=&
                      & max( sigma_st_all(i,j) -  sigma_st_all(i-1,j),0.0_db) 
                  endif
                  if(i==1 .and. j/=N) then
                     sigma_nu(i,j)=&
                     & max(sigma_st_all(i,j) - sigma_st_all(i,j+1),0.0_db)
                  endif
                  if(i==1 .and. j==N) sigma_nu(i,j)=sigma_st_all(i,j)

               enddo
            enddo
            write(73,*) T, fun_L
            do i=1,N
               write(74,*) sigma_nu(i,:)
               write(75,*) sigma_st_all(i,:)
            enddo
        endif

        if(fold_profile) then ! <prod_k=s^t (1-sigma_k)*m_k>
         fun_L=0._db

         do l_index=1,N
            do i=1,N-l_index+1
               j=i+l_index-1
               if(i/=1 .and. j/=N) then
                  fun_L(l_index)=fun_L(l_index)+&
                  & ( folded_ab_matrix(i,j) - folded_ab_matrix(i,j+1) - &
                  &   folded_ab_matrix(i-1,j) + folded_ab_matrix(i-1,j+1))
               endif
               if(i/=1 .and. j==N) then
                  fun_L(l_index)=fun_L(l_index)+&
                   & ( folded_ab_matrix(i,j) - folded_ab_matrix(i-1,j)) 
               endif
               if(i==1 .and. j/=N) then
                  fun_L(l_index)=fun_L(l_index)+&! max(,0) to avoid case L=i=j=1, <ij> = 0 and <i,j+1> = 1
                  & max(folded_ab_matrix(i,j) - folded_ab_matrix(i,j+1),0.0_db)
               endif
               if(i==1 .and. j==N) fun_L(l_index)=fun_L(l_index)+folded_ab_matrix(i,j)
            enddo !i
            fun_L(l_index)=fun_L(l_index)/(N-l_index+1) !normalization
         enddo

         do i=1,N
            do j=i,N
               if(i/=1 .and. j/=N) then
                  fold_nu(i,j)= &
                  & ( folded_ab_matrix(i,j) - folded_ab_matrix(i,j+1) - &
                  &   folded_ab_matrix(i-1,j) + folded_ab_matrix(i-1,j+1))
               endif
               if(i/=1 .and. j==N) then
                  fold_nu(i,j)=&
                   folded_ab_matrix(i,j) -  folded_ab_matrix(i-1,j)
               endif
               if(i==1 .and. j/=N) then
                  fold_nu(i,j)=&
                  folded_ab_matrix(i,j) - folded_ab_matrix(i,j+1)
               endif
               if(i==1 .and. j==N) fold_nu(i,j)=folded_ab_matrix(i,j)

            enddo
         enddo
         write(48,*) T, fun_L
         do i=1,N
            write(47,*) fold_nu(i,:)
            write(49,*) folded_ab_matrix(i,:)
         enddo
        endif

        do i=1,N
           write(50,*) T, i, Mi(i),sigmai(i),(Mi(i)-sigmai(i)) ! magnetization for each residue
        enddo

       if(show_cmd_output) write(*,*) 'T, Cp, <m>, <s>:', T, R*ConR, Mavg, sigmaavg ! the magnetizations are the average of the whole system, not just of one residue or island

        T=T+deltaT_array(q)
        q=q+1
     enddo
     cden=cden+deltacden
  enddo

  if(wFprof) close (30)
  if(wmprof) close (35)
  if(wstr)   close (40)
  if(wProd_ms) close(70)
  if(f_L) then
   close(73)
   close(74)
   close(75)
  endif
  if(fold_profile) then
   close(47)
   close(48)
   close(49)
  endif
  close (20)
  close(50) 
  close(60)
  close(94)


if(show_cmd_output)  write(*,*) '***** THE END *****'

end program WSME_genTden_loopy


!************************************************
!************************************************
!************************************************
!************************************************



subroutine read_init(& !we call this routine at the very start of the main
                                !    O:
     &     parv)
  !More output hidden in module globalpar, protdeppar

  use defreal
  use phys_const
  ! use phenom_const
  use globalpar
  !  use protdep_phenom_par
  use protdep_par
  implicit none
  real(kind=db),intent(out):: parv(nparmax)
  character(len=80):: mapaCalpha
  character(len=80):: cmapfile
  character(len=80):: elecmapfile
  character(len=80):: solvmapfile
  character(len=80):: MapasContacto
  integer::  i,j,l,p,io_status
  real(kind=db):: difftot,nASA,nct,valor
  double precision s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12
  ! calculated with the logarithm to better manage big numbers. 27*log(10) is the conversion from liters to A^3
  prefac_kappa = 0.5*(log(2.)+log(Navo)+2*log(qe)-27*log(10.)-log(vac_eps0)-log(kB))
  prefac_kappa = exp(prefac_kappa)
  
  read(*,*) Mw
  read(*,*) (parv(i),i=1,nparmax)
  read(*,*) Tmin,Tmax,deltaT,T_ref
  read(*,*) nexp
  read(*,*) cdenmin,cdenmax,deltacden

  read(*,*) ST_length
  allocate(S_interval(ST_length), T_interval(ST_length))
  read(*,*) (S_interval(p), p=1,ST_length), (T_interval(p), p=1,ST_length)

  read(*,*) mapaCalpha     !rCalpha
  read(*,*) cmapfile       !PBI.map
  read(*,*) elecmapfile    !PBI_elec.map
  read(*,*) solvmapfile    !cmapASA_PBI_uf_HP.dat
  read(*,*) expfile        !file with experimental data for T and Cp (only used while optimization)
  read(*,*) MapasContacto  !ct (this is a flag to use contact map instead of ASA)
  read(*,*) wEave
  read(*,*) wC
  read(*,*) wMave
  read(*,*) wfoldfr
  read(*,*) wstr
  read(*,*) wFprof
  read(*,*) wmprof
  read(*,*) wMres
  read(*,*) wMisland
  read(*,*) wProd_ms
  read(*,*) onlyC
  read(*,*) SS_flag
  read(*,*) SS_breakable
  read(*,*) show_cmd_output
  read(*,*) constant_deltaT
  read(*,*) f_L
  read(*,*) fold_profile

  pdb_code = cmapfile(1:4) !.map file is in format PDB.map, for example 1DPX.map. PDB has always 4 characters. Also, pdb_code is defined in moudle protdep_par

if (show_cmd_output) then
  write(*,*) '************************************************************************************'
  write(*,*) "INPUT:"
  write(*,*) "PDB:", pdb_code
  write(*,*) "N=",N
  write(*,*) "N1=",N1
  write(*,*) "Mw=",Mw
  write(*,*) "eps=",parv(1),"q=",parv(2),"eps_eff=",parv(3)
  write(*,*) "Isolv=",parv(4),"deltaCp=",parv(5)
  write(*,*) "a=",parv(6),"b=",parv(7)
  write(*,*) "Tmin,Tmax=",Tmin,Tmax
  write(*,*) "step,T_ref=",deltaT,T_ref
  write(*,*) "cdenmin,cdenmax,deltacden=",cdenmin,cdenmax,deltacden

  write(*,'(A)', advance="no") " ST intervals: "
  do i = 1, ST_length
      write(*,'(A,I0,A,I0,A)', advance="no") "(", S_interval(i), ",", T_interval(i), ")"
      if (i < ST_length) then
          write(*,'(A)', advance="no") ", "
      end if
  end do
  print * 

  write (*,*) "rCalpha map: ",mapaCalpha
  write (*,*) "Contact map: ",cmapfile
  write (*,*) "Elecmap: ",elecmapfile
  write(*,*) "Solvmap: ",solvmapfile
  write(*,*) "Experimental data: ", expfile
  write (*,*) "Case for map using: ",MapasContacto
  write(*,*) "wEave,wC,wMave,wfoldfr,wstr,wFprof,wmprof,wMres,wMisland",wEave,wC,wMave,wfoldfr,wstr,wFprof,wmprof,wMres,wMisland
  write(*,*) "onlyC=",onlyC
  write(*,*) "Use disulfide bridges at Z level: ",SS_flag
  write(*,*) "Calculate f(L) i.e. sumatory of all (s,t) islands: ", f_L
  write(*,*) "Calculate fold profile: ", fold_profile
  write(*,*) '************************************************************************************'
endif
  nct=0 !nº contactos ct
  !VdW contact map
  delta = 0._db
  open(1,file="Input/" // trim(cmapfile))
1 read(1,*,end=2) i,j,difftot 
  !  if (i.le.j-2) then
  if (i.le.j) then
     delta(1,i,j) = difftot
     nct = nct + difftot
  endif
  go to 1
2 close(1)

  rCalpha=0.0_db
  open(2, file="Input/" // trim(mapaCalpha), status='old', action='read')
  do
     read(2, *, iostat=io_status) i, j, valor
     if (io_status /= 0) exit
     rCalpha(i,j)=valor

  end do
  close(2)

  !Contribuciones electricas
  v = 0._db
  open(31,file="Input/" // trim(elecmapfile))
  do l=1,N1
     read(31,*) v(l,1),v(l,2),v(l,3),v(l,4),v(l,5)
     !v(k,1)=residue to which the first atom involved in a-a contact k belongs; 
     !v(k,2)=residue to which the second atom involved in a-a contact k belongs; 
     !v(k,3)=charge of the first atom involved in a-a contact k; 
     !v(k,4)=charge of the first atom involved in a-a contact k; 
     !v(k,5)= a-a distance for contact k 
  enddo
  close(31)

  !!nASA=0 !nº contactos ASA
  !!!Mapa solvatación
  !!open(12,file="Input/" // trim(solvmapfile))
!!3 read(12,*,end=4) i,j,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12
  !!!write(*,*) 's1, s9 =',s1,s9
  !!!   if (i.le.j-2) then
  !!if (i.le.j) then
  !!   delta(2,i,j) = s9
   !!  nASA = nASA + s9
  !!endif
  !!go to 3
!!4 close(12)
close(12)

  select case (MapasContacto) !casos para los mapas de contacto. ct=ambos cutoff, ASA=ambos ASA, mix=VDW ct, solvatacion ASA
  case ('ct') 
     delta(2,:,:)  = delta(1,:,:) 
  case ('ASA')
     delta(1,:,:)  = delta(2,:,:) 
  case ('mix') 
  end select



  return
end subroutine read_init


!************************************************
!************************************************
!************************************************
!************************************************
!************************************************


subroutine calc_fracfold(&
                                !I:	 
     &F,&
                                !O:
     &fracfold)
  use defreal
  use protdep_par


  implicit none
  integer:: i,j,iM
  real(kind=db),intent(in):: F(0:N)
  real(kind=db),intent(out):: fracfold
  real(kind=db):: dright,dleft,Fmax,Ftot,Fext(-1:N+1)
  integer:: listamin(N),posmax

  Fext=0.
  Ftot=0.
  !         if (T.eq.T_min) Fref=Phi(1)
  Fext(-1)=10000000.
  Fext(N+1)=10000000.
  do iM=0,N
     Fext(iM)=F(iM)      !+Phi(1)-Fref
     Ftot=Ftot+dexp(-F(iM))
  enddo
  Ftot=-log(Ftot)
  j=0
  do i=0,N
     dright=Fext(i)-Fext(i+1)
     dleft=Fext(i-1)-Fext(i)
     if((dright.lt.0).and.(dleft.gt.0)) then
        j=j+1
        listamin(j)=i
     endif
  enddo
  Fmax=-1000000000.
  if(j.le.1) then 
     posmax=-1
     fracfold=-1.
  else
     do i=listamin(1),listamin(j)
        dright=Fext(i)-Fext(i+1)
        dleft=Fext(i-1)-Fext(i)
        if((dright.gt.0).and.(dleft.lt.0).and.(Fext(i).gt.Fmax)) then
           posmax=i
           Fmax=Fext(i)
        endif
     enddo
     fracfold=0.
     !            do i=1,posmax-1
     !               Zu=Zu+dexp(-Fext(i))
     !            enddo
     do i=posmax+1,N
        fracfold=fracfold+exp(-Fext(i)+Ftot)
        !               Zf=Zf+dexp(-Fext(i))
     enddo
     !            fracfold=Zf/(Zu+Zf+dexp(-Fmax))
  endif

  return
end subroutine calc_fracfold
