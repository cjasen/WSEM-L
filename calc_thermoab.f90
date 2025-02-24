module thermoab
    use defreal
    use phys_const
    use globalpar
    use protdep_par, only: N
    implicit none
    private  ! All entities are now module-private by default
    public ::  calc_thermoab,datiab

    contains

      subroutine calc_thermoab(& !calculates the contributions of each (a->b) island of m=1,s=0 in a sea of m=1,s=1
           
                                ! I:
           & a_start,leng,econtrib,& ! econtrib is auxe which is also e, the matrix with the ~h_ij calculated in calc_e_Phi
                                ! O:
           & logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai,fracfold,&
           sigma_st_ab,sigma_st_ab_all,gamma,log_ZETA, SS_matrix)
        !     use defreal
        !     use phys_const
        !     use globalpar , only : wfoldfr
        !     use protdep_par, only : N


        implicit none
        integer,intent(in) :: a_start,leng ! for region (a->b), a_start=a, leng=b-a+1

        real(kind=db),intent(in):: econtrib(:,:,:) !econtrib(4,leng,leng)
        real(kind=db),intent(out):: logZeta,EonRT,ConR,sigma,sigmai(1:N),fracfold,sigma_st_ab(1:ST_length),sigma_st_ab_all(1:N,1:N)
        real(kind=db),intent(out):: EonRTsquared,ConR_fixedconf    !PIER: added this, also in the argument 27/08/24
        real(kind=db):: nu(1:leng,1:leng),F(0:leng),gamma
        real(kind=db):: log_ZETA(N) !log of Z_{a,j}^{a,b}/Z_{a,j-1}^{a,b} . No relation with logZeta. We only use the first "leng" slots of the array
        integer :: SS_matrix(:,:)

        logical,parameter :: withmprof=.false.

        logZeta=0._db
        EonRT=0._db
        ConR=0._db
        sigma=0._db
        sigma_st_ab=0._db
        sigma_st_ab_all=0._db
        nu=0._db
        log_ZETA=0._db
        call datiab(log_ZETA,logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai,sigma_st_ab,sigma_st_ab_all,&
        & econtrib,a_start,leng,gamma, SS_matrix) !dati is the plural for data, as this calculates thermo data. Terrible name.

        return
      end subroutine calc_thermoab


   !***********************************************
   !***********************************************
   !***********************************************
   !***********************************************

   subroutine datiab(&     
                                !     O:
         & log_ZETA,logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai,sigma_st_ab,sigma_st_ab_all, &
                                !     I:
         &  econtrib,a_start,leng,gamma, SS_matrix)
      !   use defreal
      !   use protdep_par, only: N
      !   use globalpar

      implicit none
      integer,intent(in) :: a_start,leng ! for region (a,b), a_start=a, leng=b-a+1
      real(kind=db),intent(in):: econtrib(:,:,0:) !econtrib(3,leng,leng)
      real(kind=db),intent(out):: logZeta,EonRT, ConR, sigma, sigmai(1:N)
      real(kind=db) :: sigma_st_ab(1:ST_length), sigma_st_ab_all(1:N,1:N),gamma 
      real(kind=db),intent(out):: EonRTsquared,ConR_fixedconf    !PIER: added this, also in the argument 27/08/24
      integer :: i,j,k,offset,s,t,p,i_start
      real(kind=db):: H(leng,leng),O(leng,leng),A(leng+1),B(leng+1),C(leng+1),D(leng+1),E(leng+1),Z,X,Y,zaux
      real(kind=db):: Zeta,O2(leng,leng),B2(leng+1),X2
      real(kind=db):: aux,aux1,auxmin,auxmax,aux1min,aux1max,auxHmin,auxHmax !FOR CHECKS
      real(kind=db):: log_ZETA(N)
      integer :: SS_matrix(:,:)


      offset=a_start-1 !to know in which part of the whole residue chain our island starts

      H=0.0_db
      do j=1,leng
         do i=1,j
            
            aux=0.;aux1=0.;
            !PIER: here H(i,j)=beta (DeltaH(i->j)-DeltaH(i->j-1))
            !PIER: Here we are assuming that V^{l,l} =0 the interactions within a loop have no energy, like those in the unfolded regions, so that the contdition 8, 14 in the notes become DeltaV^(n,l)=DeltaV^(l,l)/2=-V^nn/2
            do k=1,j-1
               H(i,j)=H(i,j)-(1+gamma)*econtrib(1,k+offset,j+offset)/2 !PIER: changed  27/08/24, CARLOS: added gamma 26/11/24
            enddo
            do k=j+1,leng
               H(i,j)=H(i,j)-(1+gamma)*econtrib(1,j+offset,k+offset)/2 !PIER: changed  27/08/24
            enddo
            H(i,j)=H(i,j)-(1+gamma)*econtrib(1,j+offset,j+offset)
            aux=H(i,j)
            
            if (j==i) then
               H(i, i) = H(i, i)+econtrib(4,i+offset,i+offset)
               aux1=econtrib(4,i+offset,i+offset)             
            endif
            
            if (j-i>=1) then
               H(i, j) = H(i, j)+econtrib(4,i+offset,j+offset)-econtrib(4,i+offset,j+offset-1)
               aux1=econtrib(4,i+offset,j+offset)-econtrib(4,i+offset,j+offset-1)               
            endif
            H(i,j)=H(i,j)+ econtrib(4,j+offset,j-1+offset) !PIER: changed sign and econtrib(4) 27/08/24
            aux1=aux1+ econtrib(4,j+offset,j-1+offset)

            H(i,j)=exp(-H(i,j))
            if(SS_flag .and. j==SS_matrix(1,2)-offset) H(i,j)=0.0_db !Restriction of disulfide bridges
         enddo
      enddo

      ! O are, I think the Delta_Phi where Phi is the contribution of changing the state of one residue to a certain observable. O is the energy.
      if(wEave.or.wC) then
         !     energia e(2,i,j)=vij/RT
         O=0.0_db
         do j=1,leng
            do i=1,j
               !PIER: aquí estaba un error importante, en el cálculo de las energías 27/8/24
                  !   energy contribution from loop stretch a+i-1 to a+j-1
               do k=1,j-1
                  O(i,j)=O(i,j)-econtrib(2,k+offset,j+offset)/2
               enddo
               do k=j+1,leng
                  O(i,j)=O(i,j)-econtrib(2,j+offset,k+offset)/2
               enddo
               O(i,j)=O(i,j)-econtrib(2,j+offset,j+offset)
            enddo
         enddo
      endif

      ! O2 is the contribution for the specific heat (yes, another terrible name bc specific heat is e(3) and its output is O2)
      if(wC) then
         !     cal spec della configurazione e(3,i,j)=dij/R
         O2=0.0_db
         do j=1,leng
            do i=1,j
               do k=1,j-1
                  O2(i,j)=O2(i,j)-econtrib(3,k+offset,j+offset)/2
               enddo
               do k=j+1,leng
                  O2(i,j)=O2(i,j)-econtrib(3,j+offset,k+offset)/2
               enddo
               O2(i,j)=O2(i,j)-econtrib(3,j+offset,j+offset)

            enddo
         enddo
      endif

      A=0.0_db
      A(1)=1.0_db
      B=0.0_db
      B2=0.0_db
      C=0.0_db
      D=0.0_db
      E=0._db
      X=0.0_db
      X2=0.0_db
      Y=0.0_db
      sigma=0.0_db
      Zeta=1.0
      logZeta=0.0_db
      log_ZETA=0._db


      do j=1,leng
         i_start=1
         if(SS_flag .and. (SS_matrix(1,2)<j+offset)) i_start=max(SS_matrix(1,2)-offset+1,1) !if there's a Cysteine with SS-bridge between a and j, then start the calculation after the Cysteine, which breaks the loopy chain
         Z=1.0_db
         do i=i_start,j
            Z=Z+H(i,j)*A(i) ! Notes' nomenclature: xi_j = 1 + sum_i A_i^(j-1)*w_ij, here Z=xi and H(i,j)=w_ij
         enddo

         logZeta=logZeta+log(Z)
         Z=1.0_db/Z
         log_ZETA(j)=log(Z) !save Z. We only use this if we want to calculate <prod_k=i^j (1-sigma_k)><prod_k m_k>
         do i=1,j !notice here we start on i=1, not i=i_start
            A(i)=Z*H(i,j)*A(i) ! A_i^j = A_i^(j-1)*w_ij/xi_j for i<=j
         enddo
         A(j+1)=Z ! for A_(j+1)^j = 1/xi_j (now, Z=1/xi_j)

         !     energia media 
         if(wEave.or.wC) then
            do i=1,j
               B(i)=Z*H(i,j)*B(i)+O(i,j)*A(i) ! remember that O is Phi_i^j - Phi_i^(j-1)
            enddo
            B(j+1)=Z*X

            X=0.0_db
            do i=i_start,j+1
               X=X+B(i)
            enddo
         endif

         if (wC) then
            !     added by pier
            !     spec. heat of conformation {sigma}
            do i=1,j
               B2(i)=Z*H(i,j)*B2(i)+O2(i,j)*A(i)
            enddo
            B2(j+1)=Z*X2

            X2=0.0_db
            do i=i_start,j+1
               X2=X2+B2(i)
            enddo
            !     endadded by pier

            !    <E^2>
            do i=1,j
               C(i)=Z*H(i,j)*C(i)+2.0_db*O(i,j)*B(i)-(O(i,j)**2.0_db)*A(i)
            enddo
            C(j+1)=Z*Y

            Y=0.0_db
           do i=i_start,j+1
               Y=Y+C(i)
            enddo
         endif

         !     magnetization per island -> sigma(a,b), but not the one needed to calculate <prod m_k sigma_k>, just an average of sigmas
         if(wMave .or. wMres .or. wMisland) then
            do i=1,j
               D(i)=Z*H(i,j)*D(i)+A(i)
            enddo
            D(j+1)=Z*sigma

            sigma=0.0_db
            do i=i_start,j+1
               sigma=sigma+D(i)
            enddo
         endif

      enddo ! j=1,leng

      if (wProd_ms) then ! <prod_k=s^t m_k  sigma_k>
         do p=1,ST_length
            s=S_interval(p)
            t=T_interval(p)
            E=0._db
            A=0._db
            A(1)=1.0_db

            ! This section is quite inefficient because we do the same calculations in the previous and in the next section.
            do j=1,leng
               i_start=1
               if(SS_flag .and. (SS_matrix(1,2)<j+offset)) i_start=max(SS_matrix(1,2)-offset+1,1) 

               Zaux=1.0_db
               do i=i_start,j
                  Zaux=Zaux+H(i,j)*A(i)
               enddo
               Zaux=1.0_db/Zaux
               do i=1,j
                  A(i)=Zaux*H(i,j)*A(i)
               enddo
               A(j+1)=Zaux  
            
               do i=1,j 
                  if (s >= (i+offset) .and. t == (j+offset)) then ! <prod m sigma>^(ab) = O(s-i)*d_tj where O is the Heavyside function and d the Kroneker's delta
                     E(i)=Zaux*H(i,j)*E(i)+A(i)
                  else
                     E(i)=Zaux*H(i,j)*E(i)
                  end if
               end do
               E(j+1)=Zaux*sigma_st_ab(p)

               sigma_st_ab(p)=0.0_db
               do i=i_start,j+1
                  sigma_st_ab(p)=sigma_st_ab(p)+E(i)
               enddo
            
            end do ! j

         end do ! p in ST_length
      end if

      !all (s,t) islands, used to calculate f(L) which is a summatory to see the lengths of the islands in the whole protein
      if(f_L) then
         do s=1,N
            do t=s,N

               E=0._db
               A=0._db
               A(1)=1.0_db
   
               do j=1,leng
                  i_start=1
                  if(SS_flag .and. (SS_matrix(1,2)<j+offset)) i_start=max(SS_matrix(1,2)-offset+1,1) 

                  Zaux=1.0_db
                  do i=i_start,j
                     Zaux=Zaux+H(i,j)*A(i)
                  enddo
                  Zaux=1.0_db/Zaux
                  do i=1,j
                     A(i)=Zaux*H(i,j)*A(i)
                  enddo
                  A(j+1)=Zaux  
               
                  do i=1,j 
                     if (s >= (i+offset) .and. t == (j+offset)) then ! <prod m sigma>^(ab) = O(s-i)*d_tj where O is the Heavyside function and d the Kroneker's delta
                        E(i)=Zaux*H(i,j)*E(i)+A(i)
                     else
                        E(i)=Zaux*H(i,j)*E(i)
                     end if
                  end do
                  E(j+1)=Zaux*sigma_st_ab_all(s,t)
   
                  sigma_st_ab_all(s,t)=0.0_db
                  do i=i_start,j+1
                     sigma_st_ab_all(s,t)=sigma_st_ab_all(s,t)+E(i)
                  enddo
               end do ! j

            enddo
         enddo
      endif

      !   magnetization of residue k per island -> sigma_k(a,b). That it's <m_i*s_i>^(a,b) 
      if (wMave .or. wMres .or. wMisland) then
         sigmai=0.0_db
         do k=1,N !PIER: changed leng ->N . This was an error! 5/9/24
            D=0.0_db
            A=0.0_db
            A(1)=1.0_db
            do j=1,leng
               i_start=1
               if(SS_flag .and. (SS_matrix(1,2)<j+offset)) i_start=max(SS_matrix(1,2)-offset+1,1) 

               Zaux=1.0_db
               do i=i_start,j
                  Zaux=Zaux+H(i,j)*A(i)
               enddo
               Zaux=1.0_db/Zaux

               do i=1,j
                  A(i)=Zaux*H(i,j)*A(i)
               enddo
               A(j+1)=Zaux  


               do i=1,j
                  if ((j+offset) .eq. k) then
                     D(i)=Zaux*H(i,j)*D(i)+A(i)
                  else
                     D(i)=Zaux*H(i,j)*D(i)
                  endif
               enddo
               D(j+1)=Zaux*sigmai(k)

               sigmai(k)=0.0_db 
               do i=i_start,j+1
                  sigmai(k)=sigmai(k)+D(i) 
               enddo

            enddo
         enddo
      endif
      
      EonRT=X
      ConR=Y-X**2 + X2
      !PIER:ADDED THIS TO OUTPUT <E^2> and Cv(fixed_m_and_sigma):
      ConR_fixedconf=X2
      EonRTsquared=Y


      return
    end subroutine datiab
end module thermoab

