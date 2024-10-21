
 program WSME_loopy_optimizable
 !using the hamoltonian:
 !     H=- T  \sum_i q + eps*ASAunfoprotsa + sum_{i<=j} (-eps Delta_ij x_{i,j}) +\sum_i (T q+alfa c) x_{i,i})
 ! where Delta_ij represent the contribution of contact ij to the DeltaASA(U-N).
 ! Notice that i<=j also contribution from an isolated native residue
   use  defreal
   use protdep_par
   use globalpar
   use phys_const
   include 'nlopt.f'
!  include '/usr/local/include/nlopt.f'
   
   external dist
   real (kind=db)::parv(nparmax),y(nparmax) !vector of npar parameters on which the energy depends (i.e. q, eps, ...)
   real(kind=db):: Phi(3),natbase(3),Mavg
   real(kind=db), allocatable :: e(:,:,:),F(:),S(:,:),m(:,:),nu(:,:)
   real(kind=db):: FreeonRT,EonRT,EnthonRT,logZeta,M0,Minf,ConR 
   real(kind=db):: T,cden
   real(kind=db):: fracfold
   double precision d, params3(3), grad3(3),params5(5),grad5(5),minf_opt,eps,ds,dC
   double precision eps_i,ds_i,dC_i
   integer :: npar, flaggrad, ires, maxeval
   integer*8 opt,l
   REAL time_begin, time_end, tol1
   
   
   logical,parameter :: onlyC=.true.
   !if only the specific heat is requested, the input flags are redefined accordingly:
   if(onlyC) then
      wEave=.false.
      wC=.true.
      wfoldfr=.false.
      wFprof=.false.
      wMave=.false.
      wmprof=.false.
      wstr=.false.
   endif  
   
   read (*,*) N !nº residuos
   read (*,*) N1 !nº contactos eléctricos
   read (*,*) N2 !nº datos experimentales
   read (*,*) flagpar !nº de parámetros a minimizar

   allocate(delta(5,1:N,1:N)) 
   allocate(rCalpha(1:N,1:N)) 
   allocate(v(N1,5)) 
   allocate(e(4,N,N),F(0:N),S(0:N+1,0:N+1),m(0:N+1,0:N+1),nu(1:N,1:N))

   call read_init(parv)
   open(20,file='HEWL_min_ct15.dat')
   eps=parv(1)
   ds=parv(2)
   dC=parv(5)
   l=0
         parv(1)=eps
         parv(2)=ds
         parv(5)=dC
         write(*,*)'eps_i,ds_i,dC_i=',eps,ds,dC
            !********************************
            !Inicio Optimizacion
            flaggrad = 0 !flaggrad=0 (para no usar gradiente)
            grad = 0._db
            maxeval = 1000 !nº evaluaciones permitidas
            tol1=1E-1
            opt = 0 
            eval=0 !nº evaluaciones

            CALL CPU_TIME ( time_begin )!calcula el tiempo de calculo

            select case (flagpar)
            case (3) 
               npar = 3 !nº parámetros a optimizar

               !parámetros a optimizar
               params3(1) = parv(1) !epsilon
               params3(2) = parv(2) !s
               params3(3) = parv(5) !deltaCp

               !parámetros que no se optimizan
               y(1)=parv(3) !epsilon_eff
               y(2)=parv(4) !I_solv
               y(3)=parv(6) !a
               y(4)=parv(7) !b
            write(*,*) 'Optimizando eps,s,dC' 

            call dist(d,npar,params3,grad3,flaggrad,y) !calcula distancia entre curva experimental y teorica
            write(*,*)  'd_min, dolddef iniciales = ',d,d*sqrt(real(N2)-1)/real(N2) 
            
            call nlo_create(opt, NLOPT_LN_NELDERMEAD, npar)
            call nlo_set_min_objective(ires, opt, dist, y)
            !call nlo_set_xtol_abs1(ires, opt, tol1)
            !call nlo_get_ftol_abs(tol1, opt)
            call nlo_set_maxeval(ires, opt, maxeval)
            call nlo_get_maxeval(maxeval, opt)
            call nlo_optimize(ires, opt, params3, minf_opt)

            if (ires.lt.0) then
               write(*,*) 'nlopt failed!'
            else
               write(*,*) 'found min at eps,s,dC=', params3(1), params3(2), params3(3)
               write(*,*) 'd_min final, dolddef final = ',minf_opt, minf_opt*sqrt(real(N2)-1)/real(N2) 
               write(*,*) 'nº iteraciones = ', eval-1
               write(20,*) minf_opt,minf_opt*sqrt(real(N2)-1)/real(N2) , eval-1, params3(1), params3(2), params3(3)
               write(20,*) "allpars=", params3(1), params3(2),parv(3),parv(4), params3(3),parv(6),parv(7),"\n"
               flush(20)
            endif
      
            call nlo_destroy(opt)

            case (5)
               npar = 5 !nº parámetros a optimizar

               !parámetros a optimizar
               params5(1) = parv(1) !epsilon
               params5(2) = parv(2) !s
               params5(3) = parv(5) !deltaCp
               params5(4) = parv(6) !a
               params5(5) = parv(7) !b

               !parámetros que no se optimizan
               y(1)=parv(3) !epsilon_eff
               y(2)=parv(4) !I_solv
            write(*,*) 'Optimizando eps,s,dC,a,b' 
            !write(20,*) params5(1), params5(2), params5(3),params5(4),params5(5)
            !flush(20)

            call dist(d,npar,params5,grad5,flaggrad,y) !calcula distancia entre curva experimental y teorica
            write(*,*) 'd_min, dolddef iniciales = ',d,d*sqrt(real(N2)-1)/real(N2) 

            call nlo_create(opt, NLOPT_LN_NELDERMEAD, npar)
            call nlo_set_min_objective(ires, opt, dist, y)
            call nlo_set_maxeval(ires, opt, maxeval)
            call nlo_get_maxeval(maxeval, opt)
            call nlo_optimize(ires, opt, params5, minf_opt)

            if (ires.lt.0) then
               write(*,*) 'nlopt failed!'
            else
               write(*,*) 'found min at eps,s,dC=', params5(1), params5(2), params5(3),minf_opt
               write(*,*) 'a,b=',params5(4),params5(5)
               write(*,*) 'd_min final = ', minf_opt
               write(*,*) 'nº iteraciones = ', eval-1
               if(minf_opt .lt. 1.0) then
                  if(params5(1) .lt. 0) then
                     if(params5(2) .lt. 0) then
                        if(params5(3) .lt. 0) then
                           write(20,*)l,params5(1), params5(2), params5(3),params5(4),params5(5),minf_opt
                           flush(20)
                        endif
                     endif
                  endif
               endif
            endif
            write(*,*) 'd_min final, dolddef final = ',minf_opt, minf_opt*sqrt(real(N2)-1)/real(N2) 
      
            call nlo_destroy(opt)
 
         end select

         CALL CPU_TIME ( time_end )
         WRITE (*,*) 'Time of operation was ', time_end - time_begin, ' seconds'
         l=l+1
    
     
   close(20)
   write(*,*) '***** THE END *****'
   
 end program WSME_loopy_optimizable

 !*******************************************************
 !*******************************************************
 !*******************************************************
 !*******************************************************
 !*******************************************************

 subroutine dist(d,npar,params,grad,flaggrad,y)
 
 use  defreal
  use protdep_par
  use globalpar
  use phys_const
  use calc_interact
  use thermo
  use thermoab
  implicit none
  
 
 
   integer :: npar, flaggrad,contadorT !nº parametros a optimizar, flaggrad=0 (para no usar gradiente)
   double precision d, params(npar), grad(npar) !distancia entre curvas, parametros a optimizar, gradiente=(0,0,0)




  real (kind=db)::parv(nparmax),y(nparmax) !vector of npar parameters on which the energy depends (i.e. q, eps, ...)
  real(kind=db):: Phi(3),natbase(3)
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
  real(kind=db):: ConRtest,ConRfixave !PIER: checks 27/8/24
  integer:: i,j,k,iM,l,p								 
  integer:: aaux,baux !variables to go through islands								 
  !  real(kind=db):: nu(1:N,1:N)

  logical,parameter :: onlyC=.false.

  !we give dimension to matrix
  

  allocate(e(4,N,N),F(0:N),S(0:N+1,0:N+1),m(0:N+1,0:N+1),nu(1:N,1:N),logZetaab(1:N,1:N), EonRTab(1:N,1:N),&
       & ConRab(1:N,1:N), sigmaab(1:N,1:N),Hn(1:N,1:N), Un(1:N,1:N),Cvn(1:N,1:N),Mi(1:N),sigmai(1:N),sigmaij(1:N,1:N),Mij(1:N,1:N),&
       & sigmaiaux(1:N),sigmaiab(1:N,1:N,1:N)) 		
  allocate(EonRTabsquared(1:N,1:N),ConRab_fixedconf(1:N,1:N)) !PIER: added this 27/8/24
  allocate(EonRTabtot(1:N,1:N),EonRTabsquaredtot(1:N,1:N),ConRab_fixedconftot(1:N,1:N),Unsquared(1:N,1:N)) !PIER: added this 27/8/24


 aaux=1
  baux=N
  allocate(auxe(4,aaux:baux,aaux:baux))


   eval=eval+1 !nº de iteraciones
   
   if(flagpar==3) then
      parv(1) = params(1) !eps
      parv(2) = params(2) !s
      parv(3)= y(1) !epsilon_eff
      parv(4) = y(2) !I_solv
      parv(5)= params(3) !deltaCp
      parv(6)= y(3) !a
      parv(7)= y(4) !b
   endif

   if(flagpar==5) then
      parv(1) = params(1) !eps
      parv(2) = params(2) !s
      parv(3)= y(1) !epsilon_eff
      parv(4) = y(2) !I_solv
      parv(5)= params(3) !deltaCp
      parv(6)= params(4) !a
      parv(7)= params(5) !b
   endif

   d = 0._db
   e=0._db
   contadorT=0
       do p=1,N2 !N2 = nº datos exp. Declarado en el modulo protdep_par
         T = T_exp(p) !T_exp declarado en el modulo protdep_par
!          FreeonRT=0.
	  contadorT=contadorT+1
	  write(*,*) "temperatura numero:", contadorT
          logZeta=0.
          EonRT=0.
          ConR=0.
          Mavg=0.
          ConR=0.
          fracfold=0.
       
          call calc_e_Phi(&
                !     I:
             &          T,cden,parv,&
             !     O:
             &          e,Phi,natbase)
             
        auxe=e(1:4,1:N,1:N)
             
        logZetaab=0._db
        EonRTab=0._db
        ConRab=0._db
        sigmaab=0._db
        sigmaiab=0._db
        EonRTabsquared=0._db !PIER: added this and below (4 lines) 27/8/24
        ConRab_fixedconf=0._db
        EonRTabsquaredtot=0._db
        ConRab_fixedconftot=0._db
            
         
          
         do aaux=1,N
           do baux=aaux,N
             call calc_thermoab(aaux,baux-aaux+1,auxe,&
                   &  logZetaaux,EonRTaux,EonRTsquaredaux,ConR_fixedconfaux,ConRaux,sigmaaux,sigmaiaux,fracfold)
!!$ !PIER: To build the pure WSME limit, without loops, comment lines below: from here...
              logZetaab(aaux,baux)=logZetaaux
!              write (*,*) "check: a,b,logZetaab= ",aaux,baux,logZetaaux
              EonRTab(aaux,baux)=EonRTaux
              EonRTabsquared(aaux,baux)=EonRTsquaredaux !PIER: added this  27/8/24
              !              ConRab(aaux,baux)=ConRaux
              ConRab_fixedconf(aaux,baux)=ConR_fixedconfaux !PIER: added this  27/8/24
              sigmaab(aaux,baux)=sigmaaux
              do k=1,N
                 sigmaiab(k,aaux,baux)=sigmaiaux(k)
              enddo
!!$              !...to here.
           enddo
        enddo
                
           
   
            Hn=0.0_db
   	do j=1,N
      	   do i=1,j

              do k=i,j
                 do l=k,j
                    Hn(i,j)=Hn(i,j)+auxe(1,k,l)
                    if (k.eq.l) then
                       Hn(i,j)=Hn(i,j)- auxe(4,k,k) !PIER: added this 27/08/2024
                    endif
                 enddo
              enddo
              ! Hn(i,j)=(effective energy of the i,j native island)/RT (being 0 the completely unfolded energy as reference)
!              write (*,*) "check: a,b, betaHnn(a,b)= ",i,j,Hn(i,j)
           enddo
  	enddo

        logZetaab=logZetaab- Hn !PIER: changed this sign 27/08/2024

!!!!!!!!!        ARREGLAR AQUI ABAJO la contribución a E^2(a,b) que sera EonRTsquared(a,b)+2*Hn(a,b)EonRT(a,b)+Unsquared (este hay que calcularlo). El calor specifico nativo CVn debería estar bien, pero ConRab=ConRab+Cvn debería ser ConRab_fixedconf(a,b)=ConRab_fixedconf(a,b)+Cvn
        !PIER: added this part for Esquared and C_fixedconf.


        Un=0.0_db
        Unsquared=0.0_db !PIER: added this  27/8/24
        do j=1,N
      	   do i=1,j
              do k=i,j
                 do l=k,j
                    Un(i,j)=Un(i,j)+auxe(2,k,l)
                 enddo
              enddo
              Unsquared(i,j)=Un(i,j)**2   !PIER: added this  27/8/24
           enddo
  	enddo
        do j=1,N !PIER: added this do loop  27/8/24
           do i=1,j
              EonRTabsquaredtot(i,j)=EonRTabsquared(i,j)+2*EonRTab(i,j)*Un(i,j)+Unsquared(i,j)
           enddo
        enddo
        EonRTabtot=EonRTab+Un



        Cvn=0.0_db
        do j=1,N
      	   do i=1,j
              do k=i,j
                 do l=k,j
                    Cvn(i,j)=Cvn(i,j)+auxe(3,k,l)
                 enddo
              enddo
           enddo
  	enddo
  	ConRab_fixedconftot=ConRab_fixedconf+Cvn !PIER: added this  27/8/24


  	
  	
          
  call calc_thermo2(N, logZetaab, EonRTabtot,EonRTabsquaredtot, ConRab_fixedconftot, sigmaab,sigmaiab,&
             &logZeta, EonRT, ConR,ConRfixave, Mavg, sigmaavg,fracfold,Mi,sigmai,Mij,sigmaij)


          !write(*,*)'Temperatura y logZeta=',T,logZeta
          ConR=ConR+Phi(3)
          write(*,*)'T y ConR=',T,R*ConR
          d = d + ABS(R*ConR - C_exp(p))*ABS(R*ConR - C_exp(p)) !C_exp declarado en el modulo protdep_par
            
       enddo
       

       !d = sqrt(d)/(real(N2)) 
       d=sqrt(d/(real(N2)-1.))
       write(*,*) "allpars=",(parv(i),i=1,7)
   write(*,*)'distancias=',d,d*sqrt(real(N2)-1)/real(N2)

 end subroutine dist

 !************************************************************
 !************************************************************
 !************************************************************
 !************************************************************
 !************************************************************

 
 subroutine read_init(&
      !    O:
      &     parv)

   use defreal
   use phys_const
  ! use phenom_const
   use globalpar
 !  use protdep_phenom_par
   use protdep_par
   implicit none
   real(kind=db),intent(out):: parv(nparmax)
   character(len=80):: cmapfile
   character(len=80):: elecmapfile
   character(len=80):: expfile
   character(len=80):: solvmapfile
   integer::  i,j,l,io_status
   real(kind=db):: difftot,valor
   double precision s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12
   prefac_kappa = 0.5*(log(2.)+log(Navo)+2*log(qe)-27*log(10.)-log(vac_eps0)-log(kB))
   prefac_kappa = exp(prefac_kappa)
 ! calculated with the logarithm to better manage big numbers. 27*log(10) is the conversion from liters to A^3
 
 !  input file format:
 !  N  !n of residues
 !  N1 !number of electrostatic contacts
 !  Mw !protein mass
 !  csi DeltaS eps_eff I DeltaC a b  !as in Naganathan 2012
 !  Tmin, Tmax, deltaT, Tref 
 !  concdenamin, concdenatmax, deltaconcdenat !(en nuestro caso, 0,0,1 ; el valor diferente de 1 sirve para salir del bucle)
 !  cmapfile
 !  elecmapfile
 !  solvmapfile
 !  case flag: ct=mapas de contacto a partir de cut off, ASA= mapas a partir de las ASAs, mix= VDW con cutfoff, solv con ASA
 !  wEave   ! =.true., .false. triggers the calculation of the average energy
 !  wC      ! =.true., .false. triggers the calculation of the specific heat
 !  wMave   ! =.true., .false. triggers the calculation of the average native fraction
 !  wfoldfr ! =.true., .false. triggers the calculation of the folded fraction (Zf/Z, folding constant)
 !  wstr    ! =.true., .false. triggers the calculation of the native strings
 !  wFprof  ! =.true., .false. triggers the calculation of the free energy profiles 
 !  wmprof  ! =.true., .false. triggers the calculation of the profile for folding probability of any residue <m_i>_M 
   read(*,*) Mw
   read(*,*) (parv(i),i=1,nparmax)
   read(*,*) Tmin,Tmax,deltaT,T_ref
   read(*,*) cdenmin,cdenmax,deltacden
   read(*,*) cmapfile
   read(*,*) elecmapfile
   read(*,*) solvmapfile
   read(*,*) expfile
   read (*,*) MapasContacto
   read(*,*) wEave
   read(*,*) wC
   read(*,*) wMave
   read(*,*) wfoldfr
   read(*,*) wstr
   read(*,*) wFprof
   read(*,*) wmprof
   read(*,*) wMres
   read(*,*) wMisland

   write(*,*) '************************************************************************************'
   write(*,*) "INPUT:"
   write(*,*) "N=",N
   write(*,*) "N1=",N1
   write(*,*) "N2=",N2
   write(*,*) "Mw=",Mw
   write(*,*) "eps=",parv(1),"q=",parv(2),"eps_eff=",parv(3)
   write(*,*) "Isolv=",parv(4),"deltaCp=",parv(5)
   write(*,*) "a=",parv(6),"b=",parv(7)
   write(*,*) "Tmin,Tmax=",Tmin,Tmax
   write(*,*) "step,T_ref=",deltaT,T_ref
   write(*,*) "cdenmin,cdenmax,deltacden=",cdenmin,cdenmax,deltacden
   write (*,*) "contact map: ",cmapfile
   write (*,*) "elecmap: ",elecmapfile
   write(*,*) "solvmap: ",solvmapfile
   write (*,*) "expfile: ",expfile
   write (*,*) "Caso: ",MapasContacto
   write(*,*) "wEave,wC,wMave,wfoldfr,wstr,wFprof,wmprof= ",wEave,wC,wMave,wfoldfr,wstr,wFprof,wmprof,wMres,wMisland
   write(*,*) '************************************************************************************'
 
   !VdW contact map
   delta = 0._db
   open(1,file=cmapfile)
 1 read(1,*,end=2) i,j,difftot 
   if (i.le.j-2) then
    delta(1,i,j) = difftot
   endif
   go to 1
 2 close(1)
 
 rCalpha=0.0_db
 open(2, file="rCalpha.txt", status='old', action='read')
 do
         read(2, *, iostat=io_status) i, j, valor
        if (io_status /= 0) exit
        rCalpha(i,j)=valor
              
 end do
 close(2)

   !Contribuciones electricas
   v = 0._db
    open(31,file=elecmapfile)
    do l=1,N1
       read(31,*) v(l,1),v(l,2),v(l,3),v(l,4),v(l,5)
       !v(k,1)=residue to which the first atom involved in a-a contact k belongs; 
       !v(k,2)=residue to which the second atom involved in a-a contact k belongs; 
       !v(k,3)=charge of the first atom involved in a-a contact k; 
       !v(k,4)=charge of the first atom involved in a-a contact k; 
       !v(k,5)= a-a distance for contact k 
    enddo
    close(31)

  !Mapa solvatación
    open(12,file=solvmapfile)
    3 read(12,*,end=4) i,j,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12
      if (i.le.j) then
       delta(2,i,j) = s9
      endif
      go to 3
    4 close(12)

   write(*,*) "OPTIMIZACIÓN:"

   select case (MapasContacto)
    case ('ct') 
      delta(2,:,:) = delta(1,:,:)
      write(*,*) 'ct case' 
    case ('ASA')
      delta(1,:,:) = delta(2,:,:) 
      write(*,*) 'ASA case' 
    case ('mix') 
      write(*,*) 'mix case' 
   end select

   !exp data
   do l=1,N2
      T_exp(l)=0._db
      C_exp(l)=0._db
   enddo
   open(10,file=expfile)
   do l=1,N2
    read(10,*) T_exp(l),C_exp(l)
   enddo
   close(10) 
   T_exp=T_exp + 273.15

   return
 end subroutine read_init
 
 
 !***************************************************
 !***************************************************
 !***************************************************
 !***************************************************
 !***************************************************

 
 subroutine profiles(&
 !     I:
      &     e,withmprofile, &
      !     O:
      &     beF,nu)
   !     CALCULATE FREE ENERGY PROFILES F(M)/RT and average nu(i,M)=<m_i(M)>
   use defreal
   use protdep_par,only:N
   implicit none
 
 !  integer:: N
   logical,intent(in) :: withmprofile
   real(kind=db),intent(in):: e(3,N,N)
   real(kind=db),intent(out):: beF(0:N),nu(1:N,1:N)
 
   integer:: i,j,k,M
   real(kind=db):: w(1:N,1:N)
   real(kind=db):: logZ,zeta,A(1:N+1),eta(0:N),alfa(1:N+1,0:N) !alfa(k,M)
   real(kind=db):: R(1:N,1:N+1,1:N) !nu(i,M), R(i,k,M)
   w=1.
   nu=0.
   zeta=1.
   A=0.
   A(1)=1./zeta
   alfa=0.
   alfa(1,0)=1.
   eta=0.
   R=0.
   eta(0)=1
   logZ=0.
 
   do j=1,N
      do i=1,j
         do k=i,j
            w(i,j)=w(i,j)*exp(e(1,k,j))
            !              e(1,i,j)=-hij/RT
         enddo
      enddo
   enddo
 
 
   do j=1,N
      zeta=1.
      do k=1,j
         zeta=zeta+A(k)*w(k,j)
      enddo
      logZ=logZ+log(zeta)
      
      do k=1,j
         A(k)=A(k)*w(k,j)/zeta
      enddo
      A(j+1)=1./zeta
 
      alfa(1,j)=alfa(1,j-1)*w(1,j)/zeta
 
      do k=2,j
         do M=j-1,j-k+1,-1
                  ! NOTICE: M from top down, to avoid using the new values (at j) instead of the old one (at j-1) when adjourning alfa(..,M) from alfa(..,M-1)!
            alfa(k,M)=alfa(k,M-1)*w(k,j)/zeta
         enddo
      enddo
 
      do M=0,j-1
         alfa(j+1,M)=eta(M)/zeta
         eta(M)=alfa(j+1,M)
      enddo
 
      do k=2,j
         do M=j-k+1,j-1
            eta(M)=eta(M)+alfa(k,M)
         enddo
      enddo
      eta(j)=alfa(1,j)
 
      if(withmprofile) then
      if(j.gt.1) then
         do M=1,j-1
            do i=1,j-1
               R(i,j+1,M)=nu(i,M)/zeta
            enddo
         enddo
      endif
 
      if(j.gt.2) then
         do i=1,j-2
            do k=i+2,j
               do M=j-1,j-k+2,-1
                  ! NOTICE: M from top down, to avoid using the new values (at j) instead of the old one (at j-1) when adjourning R(..,..,M) from R(..,..,M-1)!
                  R(i,k,M)=R(i,k,M-1)*w(k,j)/zeta
               enddo
            enddo
         enddo
      endif
 ! i=j case:
      do k=2,j
         do M=j-k+1,j-1
            nu(j,M)=nu(j,M)+alfa(k,M)
         enddo
      enddo
      nu(j,j)=alfa(1,j)
 
 ! j>i case:
      if(j.gt.1) then
         do i=1,j-1
            do M=1,j-1
               nu(i,M)=R(i,j+1,M)
            enddo
 
            do k=2,i
               do M=j-k+1,j-1
                  nu(i,M)=nu(i,M)+alfa(k,M)
               enddo
            enddo
            nu(i,j)=alfa(1,j)
 
            if(j.ge.i+2) then
               do k=i+2,j
                  do M=j-k+2,j-1
                     nu(i,M)=nu(i,M)+R(i,k,M)
                  enddo
               enddo
            endif           
         enddo
      endif
     endif
     
   enddo
 
   do M=0,N
      beF(M)=-log(eta(M))-logZ
   enddo
   
   return
 end subroutine profiles
 
 !***************************************************
 !***************************************************
 !***************************************************
 !***************************************************
 !***************************************************

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
       Ftot=Ftot+F(iM)
    enddo
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
          fracfold=fracfold+dexp(-Fext(i)+Ftot)
          !               Zf=Zf+dexp(-Fext(i))
       enddo
       !            fracfold=Zf/(Zu+Zf+dexp(-Fmax))
    endif
  
    return
 end subroutine calc_fracfold
