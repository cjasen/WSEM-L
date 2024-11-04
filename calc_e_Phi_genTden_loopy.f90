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
    real(kind=db):: w, lc, lp, d,pi
    pi=3.141592653589793
    !Inicialización de arrays a 0
    e=0._db
    ee = 0._db 
    Phi=0._db
    natbase=0._db

    !Mapas de contactos a utilizar
    cmapd=delta(1,:,:) !cmapd is the contact map calculated with cutoff distance
    ASAmap=delta(2,:,:) !ASAtotal map

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

          if (j.eq.i) then !not sure this block shouldn't be after the following if-block
            e(4,i,i-1)=DeltaS/R !PIER: notice that this represents the entropy cost of native residue i (so, this is negative)
            natbase(1)=natbase(1)-deltaS/R
         endif

         if (j-i>=4) then !where the number is lmin, the minimum lenght allowed. Previous versions had condition: j-i > 1
           !Ooka's entropy: 
           e(4,i,j)=1.5*log (real(j+1 - (i-1)))+1.5*(rCalpha(i-1,j+1)**2-3.8**2)/((j-i+2)*2*20*3.8)
           
           !Zhou's entropy:
           lc=(j-i)*3.8 !total lenght of the chain
           lp=6 !persistence length in Amstrongs (Pablo used 20, but i have NaN)
           d=rCalpha(i,j) !I use "d" to mantain the notation of Zhou et all
           w = (5.0_db * lp / (4.0_db * lc)) - (2.0_db * d**2 / lc**2) &
           & + (33.0_db * d**4 / (80.0_db * lp * lc**3)) &
           & + (79.0_db * lp**2 / (160.0_db * lc**2)) &
           & + (329.0_db * d**2 * lp / (120.0_db * lc**3)) &
           & - (6799.0_db * d**4 / (1600.0_db * lc**4)) &
           & + (3441.0_db * d**6 / (2800.0_db * lp * lc**5)) &
           & - (1089.0_db * d**8 / (12800.0_db * lp**2 * lc**6))
           !e(4,i,j) = -1.5*log(4*pi*lp*lc/3) - 3*d**2/(4*lp*lc) + log(1-w)

         else
            e(4,i,j)=500 !the number is AVERYBIGNUMBER, a penalty for forbiden loops (if >1000 it produces NaN)
         endif
       end do
    enddo

    !Disulfide bridge. As a covalent link, we overwrite the h_ij of two bridged residues with a very high value
    if (SS_flag) then
      do i=1,size(SS_matrix,1)
            e(1,SS_matrix(i,1),SS_matrix(i,2)) = -20.0_db ! Without bridge, the orgiginal value is -0.5
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

