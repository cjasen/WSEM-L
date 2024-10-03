module thermo
  use defreal
  use phys_const
  use globalpar
  !    use protdep_par, only: N
  implicit none
  private  ! All entities are now module-private by default
  public ::  calc_thermo2,dati2 !the number is bc the routine for calc_thermo_ab is called calc_thermo

contains

  subroutine calc_thermo2(& 
       & leng,logZetaab,EonRTab,EonRTabsquared,ConRabfixed,sigmaab,sigmaiab, & !PIER: changed arguments 27/8/24
                                ! O:
       & logZeta,EonRT,ConR,ConRfixed,Mavg,sigmaavg,fracfold,mi,sigmai,mij,sigmaij)   
    !     use defreal
    !     use phys_const
    !     use globalpar , only : wfoldfr
    !     use protdep_par, only : N

  
    implicit none
    integer,intent(in) :: leng ! PIER: leng will be always N

    real(kind=db),intent(in)::logZetaab(:,:),EonRTab(:,:),EonRTabsquared(:,:),ConRabfixed(:,:),sigmaab(:,:),sigmaiab(:,:,:) !PIER: changed this 27/8/24
    real(kind=db),intent(out):: logZeta,EonRT,ConR,ConRfixed,Mavg,sigmaavg,fracfold,mi(:),sigmai(:),mij(:,:),sigmaij(:,:)
    real(kind=db):: nu(1:leng,1:leng),F(0:leng)

    logical,parameter :: withmprof=.false.


    logZeta=0.
    EonRT=0.
    ConR=0.
    Mavg=0.
    sigmaavg=0.
    nu=0.
    Mi=0.0_db
    sigmai=0.0_db
    call dati2(logZeta,EonRT,ConR,ConRfixed,Mavg,sigmaavg,mi,sigmai,mij,sigmaij,&
         leng,logZetaab,EonRTab,EonRTabsquared,ConRabfixed,sigmaab,sigmaiab)
    ! write(*,*) 'main:: structF/RT=' ,-logZeta

!THE PROFILE PART BELOW WILL BE FIXED IN THE FUTURE:
!!$    if(wfoldfr) then
!!$       call profiles(econtrib,withmprof,F,nu)
!!$       call calc_fracfold(F,fracfold)
!!$    endif
    return
  end subroutine calc_thermo2


  !***********************************************
  !***********************************************
  !***********************************************
  !***********************************************


  subroutine dati2(& !PIER: changed arguments 27/8/24
                                !     O:
       & logZeta,EonRT,ConR,ConRfixed,M,sigma,mi,sigmai,mij,sigmaij, & 
                                !     I:
       &  leng,logZetaab,EonRTab,EonRTabsquared,ConRabfixed,sigmaab,sigmaiab)
    !   use defreal
    !   use protdep_par, only: N
    !   use globalpar

    implicit none
    integer,intent(in) :: leng ! PIER: leng =N in the relevant case
    real(kind=db),intent(in):: logZetaab(:,:),EonRTab(:,:),EonRTabsquared(:,:),ConRabfixed(:,:),sigmaab(:,:),sigmaiab(:,:,:) !PIER: changed this 27/8/24
    real(kind=db),intent(out):: logZeta,EonRT, ConR,ConRfixed, M, sigma,mi(:), sigmai(:),mij(:,:),sigmaij(:,:)

   !*********************************************************************************************************************************
   !sigmaab is sigmaaux or sigma, calculated in calc_thermo_ab. It has the value of <m*s>^(a->b) for each (a->b) island
   !sigmaiab has the sigma per residue <m_i*s_i>_(a->b) for each possible (a->b) island. It's an hypermatrix of 3 dimensions (residue, and two coordenates of the island)

   !sigmai is <m_i*s_i>, i.e. the probability of residue "i" to have s_i=1. Idem with mi
   !sigmaij is the probability that all residues from i->j have s=1 (and i-1 and j+1 = 0?) ? Probably <prod_k m_k*s_k>

   !*********************************************************************************************************************************

    integer :: i,j,k,l,offset
    real(kind=db):: H(leng,leng),O(leng,leng),A(leng+1),B(leng+1),C(leng+1),D(leng+1),E(leng+1),Z,X,Y
    real(kind=db):: Zeta,O2(leng,leng),O3(leng,leng),B2(leng+1),X2,Zaux,logZetaaux

    H=1.0_db
    do j=1,leng
       do i=1,j
          if (j == 1) then
             H(i, j) = exp(logZetaab(i, 1))
          else
             H(i, j) = exp(logZetaab(i, j) - logZetaab(i, j - 1)) !OK contribution of all folded string ij to the interaction with j
          end if
       enddo
    enddo


    if(wEave.or.wC) then !Now, these O are Delta_Theta, which depend on Delta_Phi (the same concept, applied to loop island on native sea istead of native sea on unfolded sea)
       !     energia e(2,i,j)=vij/RT
       O=0.0_db
       do j=1,leng
          do i=1,j
             if (j== 1) then
                O(i, j) = EonRTab(i, 1) !This things were calculated on calc_thermo_ab
             else
                O(i, j) = EonRTab(i, j) - EonRTab(i, j - 1)
             end if
          enddo
       enddo
    endif

    if(wC) then
       !     to calculate the average on m of EonRTsquared
       O2=0.0_db
       do j=1,leng
          do i=1,j 
             if (j== 1) then !PIER: changed this if 27/8/24
                O2(i, j) = EonRTabsquared(i, 1)
             else
                O2(i, j) = EonRTabsquared(i, j) - EonRTabsquared(i, j - 1)
             end if
          enddo
       enddo


       !Now O3 has the output for the specific heat, but on calc_thermo_ab O3 didn't exist and that contribution was in O2
       O3=0.0_db
       do j=1,leng
          do i=1,j 
             if (j== 1) then !PIER: changed this if 27/8/24
                O3(i, j) = ConRabfixed(i, 1)
             else
                O3(i,j)=ConRabfixed(i,j)-ConRabfixed(i,j-1)
             end if
          enddo
       enddo

    endif


    !if(wMave) then
    !     cal spec della configurazione e(3,i,j)=dij/R
    !  O3=0.0_db
    !  do j=1,leng
    !     do i=1,j
    !        if (j == 1) then
    !              O3(i, j) = sigmaab(i, i)
    !          else
    !              O3(i, j) = sigmaab(i, j) - sigmaab(i, j - 1)
    !          end if
    !                  

    !          enddo
    !      enddo
    !  endif

    A=0.0_db
    A(1)=1.0_db
    B=0.0_db
    B2=0.0_db
    C=0.0_db
    D=0.0_db
    E=0.0_db
    X=0.0_db
    X2=0.0_db
    Y=0.0_db
    M=0.0_db
    sigma=0.0_db
    Zeta=1.0
    logZeta=0.0_db
    do j=1,leng

       Z=1.0_db
       do i=1,j
          Z=Z+H(i,j)*A(i)
       enddo
       logZeta=logZeta+log(Z)
       Z=1.0_db/Z

       do i=1,j
          A(i)=Z*H(i,j)*A(i)
       enddo
       A(j+1)=Z

       !     energia media 
       if(wEave.or.wC) then
          do i=1,j
             B(i)=Z*H(i,j)*B(i)+O(i,j)*A(i)
          enddo
          B(j+1)=Z*X

          X=0.0_db
          do i=1,j+1
             X=X+B(i)
          enddo
       endif

       if (wC) then
          !     added by pier
          !     spec. heat of conformation {m}
          do i=1,j
             B2(i)=Z*H(i,j)*B2(i)+O3(i,j)*A(i)
          enddo
          B2(j+1)=Z*X2

          X2=0.0_db
          do i=1,j+1
             X2=X2+B2(i)
          enddo
          !     endadded by pier

          !    <EonRTabsquared>
          do i=1,j
             C(i)=Z*H(i,j)*C(i)+2.0_db*O(i,j)*B(i) +A(i)* ( O2(i,j) - 2.*O (i,j)* EonRTab(i,j)   )
!             C(i)=Z*H(i,j)*C(i)+O2(i,j)*A(i) !PIER: changed this  27/8/24
          enddo
          C(j+1)=Z*Y

          Y=0.0_db
          do i=1,j+1
             Y=Y+C(i)
          enddo
       endif

       !      m and sigma  (nativefraction). This should be <m> and <m*s>, for the whole chain
       if(wMave) then

          do i=1,j
             D(i)=Z*H(i,j)*D(i)+A(i)
          enddo
          D(j+1)=Z*M

          M=0.0_db
          do i=1,j+1
             M=M+D(i)
          enddo

          do i=1,j
             E(i)=Z*H(i,j)*E(i)+A(i)*(sigmaab(i,j)-sigmaab(i,j-1))

          enddo
          E(j+1)=Z*sigma

          sigma=0.0_db
          do i=1,j+1
             sigma=sigma+E(i)
          enddo

       endif
    enddo



    !    m y sigma de cada resido
    if(wMres) then
       do k=1,leng ! remember leng=N always
          D=0.0_db
          E=0.0_db
          A=0.0_db
          A(1)=1
          do j=1,leng

             Zaux=1.0_db
             do i=1,j
                Zaux=Zaux+H(i,j)*A(i)
                !write(*,*) i,j,Z
             enddo
             logZetaaux=logZetaaux+log(Zaux)
             Zaux=1.0_db/Zaux

             do i=1,j
                A(i)=Zaux*H(i,j)*A(i)
             enddo
             A(j+1)=Zaux  


             do i=1,j
                if (j == k) then
                   D(i)=Zaux*H(i,j)*D(i)+A(i) ! The Delta_Theta for <m_i> are 1 if k==j and 0 if not
                else
                   D(i)=Zaux*H(i,j)*D(i)
                endif
             enddo
             D(j+1)=Zaux*Mi(k)

             Mi(k)=0.0_db
             do i=1,j+1
                Mi(k)=Mi(k)+D(i) ! <m_i>. Notice this calculation doesn't need any data from calc_thermo_ab, as it's independent of the existence of loops
             enddo

             do i=1,j
                if (j > k .and. i .le. k) then
                   E(i)=Zaux*H(i,j)*E(i)+A(i)*(sigmaiab(k,i,j)-sigmaiab(k,i,j-1))
                elseif (j == k .and. i .le. k) then
                   E(i)=Zaux*H(i,j)*E(i)+A(i)*(sigmaiab(k,i,j))
                else
                   E(i)=Zaux*H(i,j)*E(i)
                endif !PIER: esto de arriba se podría simplificar  con solo la primera expresión, ya que sigmaiab(k,i,j) debería ser 0 si k<i o k>j

             enddo
             E(j+1)=Zaux*sigmai(k)

             sigmai(k)=0.0_db
             do i=1,j+1
                sigmai(k)=sigmai(k)+E(i)
             enddo
          enddo
       enddo
    endif


    !magnetizaciones por isla
    if (wMisland) then
       do l=1,leng
          do k=1,l
             D=0.0_db
             E=0.0_db
             A=0.0_db
             A(1)=1
             do j=1,leng

                Zaux=1.0_db
                do i=1,j
                   Zaux=Zaux+H(i,j)*A(i)
                   !write(*,*) i,j,Z
                enddo
                logZetaaux=logZetaaux+log(Zaux)
                Zaux=1.0_db/Zaux

                do i=1,j
                   A(i)=Zaux*H(i,j)*A(i)
                enddo
                A(j+1)=Zaux  


                do i=1,j
                   if (k .le. i .and. j .le. l ) then
                      D(i)=Zaux*H(i,j)*D(i)+A(i)
                   else
                      D(i)=Zaux*H(i,j)*D(i)
                   endif
                enddo
                D(j+1)=Zaux*Mij(k,l)

                Mij(k,l)=0.0_db
                do i=1,j+1
                   Mij(k,l)=Mij(k,l)+D(i)
                enddo

!CREO QUE LO DE ABAJO ESTA MAL, COMPROBAR...


                do i=1,j
                   if (k .le. i .and. j .le. l ) then                
                      E(i)=Zaux*H(i,j)*E(i)+A(i)*(sigmaab(i,j)-sigmaab(i,j-1))


                   else
                      E(i)=Zaux*H(i,j)*E(i)
                   endif

                enddo
                E(j+1)=Zaux*sigmaij(k,l)

                sigmaij(k,l)=0.0_db
                do i=1,j+1
                   sigmaij(k,l)=sigmaij(k,l)+E(i)
                enddo

             enddo
          enddo
       enddo



       do k = 1, leng
          do l = k, leng
             sigmaij(k, l) = sigmaij(k, l) / (l - k + 1)
             Mij(k,l)=Mij(k,l)/(l - k + 1)

          end do
       end do

    endif



    EonRT=X
    ConR=Y-X**2 + X2
    ConRfixed=X2
    M=M/leng 
    sigma=sigma/leng

    return

  end subroutine dati2
end module thermo

