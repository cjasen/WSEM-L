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
           & logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai,fracfold)
        !     use defreal
        !     use phys_const
        !     use globalpar , only : wfoldfr
        !     use protdep_par, only : N


        implicit none
        integer,intent(in) :: a_start,leng ! for region (a->b), a_start=a, leng=b-a+1

        real(kind=db),intent(in):: econtrib(:,:,:) !econtrib(4,leng,leng)
        real(kind=db),intent(out):: logZeta,EonRT,ConR,sigma,sigmai(1:N),fracfold
        real(kind=db),intent(out):: EonRTsquared,ConR_fixedconf    !PIER: added this, also in the argument 27/08/24
        real(kind=db):: nu(1:leng,1:leng),F(0:leng)

        logical,parameter :: withmprof=.false.


        logZeta=0._db
        EonRT=0._db
        ConR=0._db
        sigma=0._db
        nu=0._db
        call datiab(logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai,econtrib,a_start,leng) !dati is the plural for data, as this calculates thermo data. Terrible name.
        ! write(*,*) 'main:: structF/RT=' ,-logZeta


        if(wfoldfr) then
           call profiles(econtrib,withmprof,F,nu)
           call calc_fracfold(F,fracfold)
        endif
        return
      end subroutine calc_thermoab


   !***********************************************
   !***********************************************
   !***********************************************
   !***********************************************

   subroutine datiab(&     
                                !     O:
         & logZeta,EonRT,EonRTsquared,ConR_fixedconf,ConR,sigma,sigmai, &
                                !     I:
         &  econtrib,a_start,leng)
      !   use defreal
      !   use protdep_par, only: N
      !   use globalpar

      implicit none
      integer,intent(in) :: a_start,leng ! for region (a,b), a_start=a, leng=b-a+1
      real(kind=db),intent(in):: econtrib(:,:,:) !econtrib(3,leng,leng)
      real(kind=db),intent(out):: logZeta,EonRT, ConR, sigma, sigmai(1:N)
      real(kind=db),intent(out):: EonRTsquared,ConR_fixedconf    !PIER: added this, also in the argument 27/08/24
      integer :: i,j,k,offset
      real(kind=db):: H(leng,leng),O(leng,leng),A(leng+1),B(leng+1),C(leng+1),D(leng+1),Z,X,Y,zaux
      real(kind=db):: Zeta,O2(leng,leng),B2(leng+1),X2
      real(kind=db):: aux,aux1,auxmin,auxmax,aux1min,aux1max,auxHmin,auxHmax !FOR CHECKS
      integer:: auximin,auximax,auxjmin,auxjmax !FOR CHECKS

      offset=a_start-1 !to know in which part of the whole residue chain our island starts

      H=0.0_db
      auxHmax=-100000000000000.
      auxHmin=100000000000000.
      auximin=0;auxjmin=0;auximax=0;auxjmax=0;
      do j=1,leng
         do i=1,j
            aux=0.;aux1=0.;
            !PIER: here H(i,j)=beta (DeltaH(i->j)-DeltaH(i->j-1))
            !PIER: Here we are assuming that V^{l,l} =0 the interactions within a loop have no energy, like those in the unfolded regions, so that the contdition 8, 14 in the notes become DeltaV^(n,l)=DeltaV^(l,l)/2=-V^nn/2
            do k=1,j-1
               H(i,j)=H(i,j)-econtrib(1,k+offset,j+offset)/2 !PIER: changed  27/08/24
            enddo
            do k=j+1,leng
               H(i,j)=H(i,j)-econtrib(1,j+offset,k+offset)/2 !PIER: changed  27/08/24
            enddo
            H(i,j)=H(i,j)-econtrib(1,j+offset,j+offset)
            aux=H(i,j)
            
            if (j-i>1) then
               H(i, j) = H(i, j)+econtrib(4,i+offset,j+offset)-econtrib(4,i+offset,j+offset-1)
               aux1=econtrib(4,i+offset,j+offset)-econtrib(4,i+offset,j+offset-1)               
            endif
            H(i,j)=H(i,j)+ econtrib(4,j+offset,j+offset) !PIER: changed sign and econtrib(4) 27/08/24
            aux1=aux1+ econtrib(4,j+offset,j+offset)

            if(H(i,j)<auxHmin) then
               auxmin=aux
               aux1min=aux1
               auxHmin=H(i,j)
               auximin=i
               auxjmin=j
            endif
            if(H(i,j)>auxHmax) then
               auxmax=aux
               aux1max=aux1
               auxHmax=H(i,j)
               auximax=i
               auxjmax=j
            endif


!            aux=0
!            aux1=0
            H(i,j)=exp(-H(i,j))
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
      X=0.0_db
      X2=0.0_db
      Y=0.0_db
      sigma=0.0_db
      Zeta=1.0
      logZeta=0.0_db
      do j=1,leng

         Z=1.0_db
         do i=1,j
            Z=Z+H(i,j)*A(i) ! Pablo's nomenclature: xi_j = 1 + sum_i A_i^(j-1)*w_ij, here Z=xi and H(i,j)=w_ij
         enddo
         logZeta=logZeta+log(Z)
         Z=1.0_db/Z

         do i=1,j
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
            do i=1,j+1
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
            do i=1,j+1
               X2=X2+B2(i)
            enddo
            !     endadded by pier

            !    <E^2>
            do i=1,j
               C(i)=Z*H(i,j)*C(i)+2.0_db*O(i,j)*B(i)-(O(i,j)**2.0_db)*A(i)
            enddo
            C(j+1)=Z*Y

            Y=0.0_db
           do i=1,j+1
               Y=Y+C(i)
            enddo
         endif

         !     magnetization per island -> sigma(a,b)
         if(wMave .or. wMres .or. wMisland) then
            do i=1,j
               D(i)=Z*H(i,j)*D(i)+A(i)
            enddo
            D(j+1)=Z*sigma

            sigma=0.0_db
            do i=1,j+1
               sigma=sigma+D(i)
            enddo
         endif
      enddo

      !   magnetization of residue k per island -> sigma_k(a,b). That it's <m_i*s_i>^(a,b) 

      sigmai=0.0_db
      do k=1,N !PIER: changed leng ->N . This was an error! 5/9/24
         D=0.0_db
         A=0.0_db
         A(1)=1.0_db
         do j=1,leng

            Zaux=1.0_db
            do i=1,j
               Zaux=Zaux+H(i,j)*A(i)
               !write(*,*) i,j,Z
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
            do i=1,j+1
               sigmai(k)=sigmai(k)+D(i) 
            enddo

         enddo
      enddo

      EonRT=X
      ConR=Y-X**2 + X2
      !PIER:ADDED THIS TO OUTPUT <E^2> and Cv(fixed_m_and_sigma):
      ConR_fixedconf=X2
      EonRTsquared=Y


      return
    end subroutine datiab
end module thermoab

