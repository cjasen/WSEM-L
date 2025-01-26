module calc_interact
  use defreal  
  use phys_const
  use globalpar
  use protdep_par
  implicit none
  private  ! All entities are now module-private by default
  public ::  calc_e_Phi 
  
  
contains 
  !Some adjustements to include i,i interactions, and to consider correctly the enthalpy (not used in gankyrin work) 13/11/2015

  subroutine calc_e_Phi(&
                                !     I:
       &  T,c_den,y,SS_matrix, &
                                !hidden input: in modules
                                !     O:
       &  e,Phi,natbase)


    real(kind=db),intent(in):: y(nparmax),T,c_den !c_den= denaturant concentration
    integer, intent(in) :: SS_matrix(:,:) 
    real(kind=db),intent(out):: e(4,N,0:N),Phi(3),natbase(3) !native  baselines (the unfolded one are the Phi: 1->F,2->H,3->C)
    !     e(1,i,j)= h(i,j)/RT= eq (37b,c)
    !     e(2,i,j)=v(i,j)/RT  
    !     e(3,i,j)=D(i,j)/R   
    !     e(4,i,j)=-1.5*log (real(j - i+1))       coste entropico de tener un loop entre i y j. OJO! hasta ahora todo escrito en términos de interacció entre i y j ahora este termino da cuenta de la modificación a la energía libre de una región no de un enlace
    !The e(1,i,j), as in previous versions, account for the (free-)energy change upon forming the contact i,j. 
    !So, the sum of the e(1,i,j) for the completely unfolded conformation is 0. 
    integer :: i,j,l
    real(kind=db)::  ee(3,N,N)
    !     ee(1,i,j) = (Kco*q(i)q(j)/r(i,j))*exp(-kappa*r(i,j))=> Contribución eléctrica a e(1,i,j)
    !     ee(2,i,j) = Kco*q(i)q(j)*(1/r(i,j)-kappa/2)*exp(-kappa*r(i,j))=> Contribución eléctrica a e(2,i,j)
    !     ee(3,i,j) = (Kco*q(i)q(j)*kappa*(3-kappa*r(i,j))/4T)*exp(-kappa*r(i,j))=> Contribución eléctrica a e(3,i,j)
    real(kind=db):: csi_nag,DeltaC,cmapd(N,N),ASAmap(N,N),DeltaS,epseff,Ionicforce,aonR,bonR,varCp,kappa,Kco 
    real(kind=db) :: cmapd_copy(N,N), ASAmap_copy(N,N)
    real(kind=db):: w, lc, lp, d,pi
    logical :: now_flip
    pi=3.141592653589793
    !Inicialización de arrays a 0
    e=0._db
    ee = 0._db 
    Phi=0._db
    natbase=0._db

    !Mapas de contactos a utilizar
    cmapd=delta(1,:,:) !cmapd is the contact map calculated with cutoff distance
    ASAmap=delta(2,:,:) !ASAtotal map
    cmapd_copy=cmapd
    ASAmap_copy=ASAmap


    !Parámetros del modelo
    csi_nag=y(1) !eps in J/mol
    DeltaS=y(2) !s in J/molK
    epseff=y(3) !eps_eff 
    Ionicforce=y(4) !I_solv in M
    DeltaC=y(5) !DeltaC in J/molK
    aonR=y(6)*Mw/R  !a in J/(g K), Mw in g/mol
    bonR=y(7)*Mw/(1000*R) !b in J/(g K^2)
    ! according to Naganathan 2012, the specific heat's unfolded  baseline is Mw(a+(b/1000)(T-T0C)) 
    ! here we divide by R, to "dimensionally" agree with the rest of the output for specific heat

    !Parámetros no dependientes de los átomos, se precalculan fuera del bucle en átomos
    Kco=Kcoul/epseff
    varcp=((T-T_ref)-T*log(T/T_ref))
    kappa = prefac_kappa*sqrt(Ionicforce/(T*epseff)) !prefac_kappa is calculated on read_init
    !Contribuciones eléctricas. Bucle en átomos
    do l=1,N1
      !v(i,j) are the electric contributions coded in elecmapfile (PDB_elec.map). v(i,1) and v(i,2) are the residues of interaction number "i", v(i,3) and v(i,4) have charges and v(i,5) has distances
       i=int(v(l,1))
       j=int(v(l,2))
       !write(*,*) "line=",l, "i=",i, "v1=",v(l,1), "j=", j, "v2=", v(l,2)
       ee(1,i,j)= ee(1,i,j) + Kco*v(l,3)*v(l,4)*exp(-v(l,5)*kappa)/v(l,5)
       ee(2,i,j)= ee(2,i,j) + Kco*v(l,3)*v(l,4)*exp(-v(l,5)*kappa)*(1/v(l,5)-kappa/2)
       ee(3,i,j)= ee(3,i,j) + Kco*v(l,3)*v(l,4)*exp(-v(l,5)*kappa)*kappa*(3-kappa*v(l,5))/(4*T)
    enddo
    
    !Hamiltoniano completo+contribuciones C dependiente de T. Bucle en residuos
    do i=1,N
       do j=i,N
          if(i.le.j-1) then
             e(1,i,j)= e(1,i,j)+ ee(1,i,j)
             e(2,i,j)= e(2,i,j) + ee(2,i,j)
             e(3,i,j)= e(3,i,j) + ee(3,i,j)
          endif
          if(i.le.j-2) then
             e(1,i,j)=e(1,i,j)+ csi_nag*cmapd(i,j) + varcp*DeltaC*ASAmap(i,j)
             e(2,i,j)= e(2,i,j) + csi_nag*cmapd(i,j)+ DeltaC*(T-T_ref)*ASAmap(i,j)
             e(3,i,j)=e(3,i,j) + DeltaC*ASAmap(i,j)
          endif
          e(1,i,j)= e(1,i,j)/(R*T) !PIER: CHANGED SIGN 27/8/24
          e(2,i,j)= e(2,i,j)/(R*T)
          e(3,i,j)= e(3,i,j)/(R)
          natbase(1)=natbase(1)+e(1,i,j) !PIER: CHANGED SIGN 27/8/24
          natbase(2)=natbase(2)+e(2,i,j)
          natbase(3)=natbase(3)+e(3,i,j)
          if (j.eq.i) then !PIER: CHANGED THE STRUCTURE OF THE e, putting all the entropy in e(4,,)  27/8/24
             e(4,i,i-1)=DeltaS/R !PIER: notice that this represents the entropy cost of native residue i (so, this is negative)
             natbase(1)=natbase(1)-deltaS/R !PIER: CHANGED  27/8/24
          endif

          if (j-i>=0) then
             if((i==1).or.(j==N)) then
               e(4,i,j)=0._db
             else
               !Ooka's entropy:   
               !e(4,i,j)=1.5*log (real(j - i+2))+1.5*(rCalpha(i-1,j+1)**2-3.8**2)/((j-i+2)*2*20*3.8) 

               !Zhou's entropy:
               lc=(j-i+2)*3.8 !total lenght of the chain
               lp=2.04 !persistence length in Amstrongs (Pablo used 20, but i have NaN)
               d=rCalpha(i-1,j+1) !I use "d" to mantain the notation of Zhou et all
               w = (5.0_db * lp / (4.0_db * lc)) - (2.0_db * d**2 / lc**2) &
                     & + (33.0_db * d**4 / (80.0_db * lp * lc**3)) &
                     & + (79.0_db * lp**2 / (160.0_db * lc**2)) &
                     & + (329.0_db * d**2 * lp / (120.0_db * lc**3)) &
                     & - (6799.0_db * d**4 / (1600.0_db * lc**4)) &
                     & + (3441.0_db * d**6 / (2800.0_db * lp * lc**5)) &
                     & - (1089.0_db * d**8 / (12800.0_db * lp**2 * lc**6))

               !e(4,i,j) = 1.5*log(4*pi*lp*lc/3) + 3*d**2/(4*lp*lc) - log(1.6-w) +1 ! +1 is for gyration radio

               lc=(j-i+2)*3.8
               lp=1.4 ! 0.04 for 1dpx with 1 ss bond
               d=rCalpha(i-1,j+1)

               call compute_Q_r(lp/lc, d/lc, e(4,i,j))
               e(4,i,j)=-log(e(4,i,j)) 
               !if(.not. ieee_is_finite(e(4, i, j))) e(4,i,j)=200._db
               endif
            else
            e(4,i,j)=60.0_db
          endif
          if(isnan(e(4,i,j))) e(4,i,j)=e(4,i,j-1)
       end do
    enddo

    !Disulfide bridge. As a covalent link, we overwrite the h_ij of two bridged residues with a very high value
    if (SS_breakable) then
      do i=1,size(SS_matrix,1)
            e(1,SS_matrix(i,1),SS_matrix(i,2)) = -10.0_db ! Without bridge, the orgiginal value is -0.5
      end do
    endif 

    !Cálculo de las Phis
    Phi(1)=aonR*((T-T0C)/T-log(T/T0C))+ bonR*((T0C**2-T**2)/(2*T)+T0C*log(T/T0C)) !Free energy
    natbase(1)=natbase(1)+ Phi(1)

    Phi(2)=aonR*(T-T0C)/T+bonR*((T-T0C)**2)/(2.*T) !Enthalpy
    natbase(2)=natbase(2)+ Phi(2)

    Phi(3)=aonR+bonR*(T-T0C) !Specific heat
    natbase(3)=natbase(3)+ Phi(3)

    return
  end subroutine calc_e_Phi

  subroutine compute_Q_r(kappa, r, Q_r) ! Becker's entropy
   use defreal
   implicit none
 
   ! Argumentos de entrada y salida
   real(kind=db), intent(in) :: kappa, r
   real(kind=db), intent(out) :: Q_r
 
   ! Declaración de variables
   real(kind=db) :: a, b, c, d, J_SYD, term1, term2, term3, term4
   real(kind=db), dimension(2, 3) :: cij
   integer :: i, j
   real(kind=db), dimension(2) :: i_values
   real(kind=db), dimension(3) :: j_values
 
   ! Inicialización de constantes
   a = 14.054_db
   b = 0.473_db
 
   c = 1.0_db - (1.0_db + (0.38_db * kappa**(-0.95_db))**(-5))**(-1.0_db / 5.0_db)
 
   if (kappa < 1.0_db / 8.0_db) then
     d = 0.0_db
   else
     d = 1.0_db / (0.177_db / (kappa - 0.111_db) + 6.40_db * (kappa - 0.111_db)**0.783_db)
   end if
 
   cij = reshape([ &
     -3.0_db / 4.0_db, 23.0_db / 64.0_db, -7.0_db / 64.0_db, &
     -1.0_db / 2.0_db, 17.0_db / 16.0_db, -9.0_db / 16.0_db], &
     shape(cij), order=[1, 2])
 
   if (kappa <= 1.0_db / 8.0_db) then
     J_SYD = (3.0_db / 4.0_db * 3.141592653589793_db * kappa)**(3.0_db / 2.0_db) * (1.0_db - 5.0_db * kappa / 4.0_db)
   else
     J_SYD = 112.04_db * kappa**2 * exp(0.246_db / kappa - a * kappa)
   end if
 
   term1 = ((1.0_db - c * r**2) / (1.0_db - r**2))**(5.0_db / 2.0_db)
 
   term2 = 0.0_db
   i_values = [-1.0_db, 0.0_db]
   j_values = [1.0_db, 2.0_db, 3.0_db]
   
   do i = 1, size(i_values)
     do j = 1, size(j_values)
       term2 = term2 + cij(i, j) * kappa**i_values(i) * r**(2.0_db * j_values(j))
     end do
   end do
 
   term2 = exp(term2 / (1.0_db - r**2))
 
   term3 = exp(-d * kappa * a * b * (1.0_db + b) * r**2 / (1.0_db - b**2 * r**2))
   call bessel_i0(-d * kappa * a * (1.0_db + b) * r / (1.0_db - b**2 * r**2), term4)
 
   Q_r = J_SYD * term1 * term2 * term3 * term4
 end subroutine compute_Q_r
 
 subroutine bessel_i0(x, result)
   implicit none
 
   ! Argumentos de entrada y salida
   real(8), intent(in) :: x
   real(8), intent(out) :: result
 
   ! Constantes
   real(8), parameter :: eps = 1.0d-10  ! Tolerancia para la serie
   real(8), parameter :: x_threshold = 3.75d0  ! Umbral para usar la serie o la aproximación asintótica
   real(8) :: abs_x, y, t, sum, term
   integer :: k
 
   ! Valor absoluto de x
   abs_x = abs(x)
 
   if (abs_x <= x_threshold) then
     ! Serie de potencias para valores pequeños de x
     t = abs_x / x_threshold
     t = t * t
     sum = 1.0d0
     term = 1.0d0
     k = 1
     do
       term = term * (t / (4.0d0 * k * k))
       sum = sum + term
       if (term < eps) exit
       k = k + 1
     end do
     result = sum
   else
     ! Aproximación asintótica para valores grandes de x
     y = x_threshold / abs_x
     result = exp(abs_x) / sqrt(abs_x) * (1.0d0 + y * (3.5156229d0 + y * (3.0899424d0 + y * (1.2067492d0 + &
               y * (0.2659732d0 + y * (0.0360768d0 + y * 0.0045813d0))))))
   end if
 
 end subroutine bessel_i0

! !*********************************************
! function pol(v,x,deg)
! use defreal 
! implicit none
! real(kind=db):: pol
! integer :: deg,i      
! real(kind=db):: v(0:deg),x
! pol=v(deg)
! do i=deg-1,0,-1
!   pol=pol*x+v(i)
! enddo
! return
! end function pol

! !**********************************************
end module calc_interact

