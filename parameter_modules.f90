module defreal
  implicit none
  save
  !choice of precision
  integer, parameter :: db = selected_real_kind(p=15, r=307) !64-bits reals)

end module defreal
!********************************************
module phys_const
  use defreal
  implicit none
  save
  real(kind=db), parameter :: Kcoul=1389.0_db !Coulomb constant:
  !     in KJ A/mol
  real(kind=db), parameter :: qe = 1.602176634d-19 !elementary charge in C
  real(kind=db), parameter :: pi=acos(-1.D0)
  real(kind=db), parameter :: kB=1.380649d-23 !Boltzmann constant in J/K
  real(kind=db), parameter :: vac_eps0=8.8541878128d-22 !vacuum permittivity constant in C^2/(J A)
  real(kind=db), parameter :: Navo=6.02214076d23
  real(kind=db),parameter:: T0C=273.0_db
  real(kind=db) :: prefac_kappa
  real(kind=db), parameter :: R=0.008314462_db !in KJ/(K mol). setting R=1 is the choice by Alessandro and Marco
  real(kind=db),parameter:: cKj=0.0041868_db !1 cal = 0.0041868 KJ
end module phys_const

!!********************************************

 module globalpar
   use defreal
   implicit none
   save
   integer,parameter:: nparmax=8 !!number of model parameters
   real(kind=db) :: Tmin, Tmax,deltaT,T_ref
   integer :: ST_length, nexp
   integer, allocatable :: S_interval(:), T_interval(:) ! intervals to calculate <prod_k=S^T m_k sigma_k>
   real(kind=db) :: cdenmin, cdenmax,deltacden
   logical :: wEave,wC,wMave,wfoldfr,wstr,wFprof,wmprof, wMres, wMisland, wProd_ms, constant_deltaT,fold_profile
   logical :: onlyC, SS_flag, show_cmd_output, f_L, SS_breakable
   character(len=80)::expfile
   ! wEave= calculates the average energy; wC: calculate specific heat, wMave=calculate average "magnetization", wfoldfr: with fraction folded; wst=with strings, wFprof=with F profiles, wmprof=with m profiles
   ! character(len=80):: outtherm,outprof,outstring,outmprofile
   ! real(kind=db):: parv(nparmax)
 end module globalpar


!********************************************

 module opt_aux
  use defreal
  implicit none
  save

  integer:: nexp, eval, N_res, N_elec, Mw_opt !some of this are parameters already defined, but I do it again to avoid problems
  real(kind=db), allocatable :: T_exp(:),C_exp(:) !an array to carry the experimental datapoints
  real(kind=db):: y(8)
  character(len=80):: expfile, simfile !the name of the .txt with the experimental datapoints and name of the .txt with the simulated data
  character(len=4):: pdb_code_opt

end module opt_aux

!*******************************************


module protdep_par
  use defreal
  implicit none
  save
 !distancia entre carbonos alfa calculado con el programa radPDB.py
  integer :: N,N1 !Number of residues, number of elec contacts
  character(len=4) :: pdb_code ! PDB of the protein. For example, 1DPX
  real(kind=db) :: Mw ! protein mass [grams/mol]
  real(kind=db), allocatable :: delta(:,:,:),rCalpha(:,:) !contact maps
  real, allocatable :: v(:,:)!list of atom-atom contacts, with their charge and distance; namely:
  !v(k,1)=residue to which the first atom involved in a-a contact k belongs;
  !v(k,2)=residue to which the second atom involved in a-a contact k belongs;
  !v(k,3)=charge of the first atom involved in a-a contact k;
  !v(k,4)=charge of the first atom involved in a-a contact k;
  !v(k,5)= a-a distance for contact k


  !  integer,allocatable :: fluo_res(:) !list of fluorescent residues' indices (1..nfluo)
!  integer :: nfluo,nfluocts
!  real(kind=db), allocatable :: fluocts(:,:,:)
!! fluocts(i,j,k) is the fraction of ASA exposed of fluorescent residue i in the isolated  native string j,k. i is reset to 1,2,etc, to save space.j,k run from 1 to N, but fluocts(i,j,k)=0 for irrelevant configurations

  !  character(len=1),allocatable:: secstr(:) SI EN EL FUTURO VAMOS A PONER ENTROPÍA DEPENDIENTE DE LA ESTRUCTURA SECUNDARIA, ESTO ES RELEVANTE
!  integer :: nsecstr

!!$ secstr(i)= one letter DSSP code:
!!$    * H = alpha helix
!!$    * B = residue in isolated beta-bridge
!!$    * E = extended strand, participates in beta ladder
!!$    * G = 3-helix (3/10 helix)
!!$    * I = 5 helix (pi helix)
!!$    * T = hydrogen bonded turn
!!$    * S = bend

end module protdep_par

!********************************************
!module phenom_const
!  use defreal
!  use phys_const
!  implicit none
!  save


!  real(kind=db),parameter:: T0=298.15_db
!  !     from Gomez..Freire, Proteins 22,404(1995), table 5:
!  real(kind=db),parameter::  vfr=0.28_db*cKj/R
!  real(kind=db),parameter::  wfr=9.75_db*0.0001_db*cKj/R
!  real(kind=db),parameter::  pfr=8.7_db*0.001_db*cKj/R
!  real(kind=db),parameter::  qfr=6.43_db*0.0001_db*cKj/R
!  real(kind=db),parameter:: a1fr=0.45_db*cKj/R
!  real(kind=db),parameter::  a2fr=2.63_db*0.0001_db*cKj/R
!  real(kind=db),parameter::  a3fr=-4.2_db*0.00001_db*cKj/R
!  real(kind=db),parameter::  b1fr=-0.265_db*cKj/R
!  real(kind=db),parameter::  b2fr=2.85_db*0.0001_db*cKj/R
!  real(kind=db),parameter:: b3fr=4.31_db*0.00001_db*cKj/R
!  !     above:values converted to use KJ/K mol as a measure of spec.heat.
!  !     Areas in A^2

!!  real(kind=db),parameter:: facCys=20._db
!  real(kind=db):: facCys


!  real(kind=db),dimension(3),parameter:: cS_A=(/ 0._db, 4.1_db, 4.1_db/)
!  real(kind=db),dimension(3),parameter:: cS_R=(/7.11_db, 2.56_db, 9.67_db/)
!  real(kind=db),dimension(3),parameter:: cS_N=(/3.29_db, 5.64_db, 8.93_db/)
!  real(kind=db),dimension(3),parameter:: cS_D=(/2._db,  5.56_db, 7.56_db/)
!  real(kind=db),dimension(3),parameter:: cS_C=(/3.55_db ,4.01_db ,7.56_db/)
!  real(kind=db),dimension(3),parameter:: cS_Q=(/ 5.02_db, 5.52_db, 10.54_db/)
!  real(kind=db),dimension(3),parameter:: cS_E=(/ 3.53_db, 5.67_db, 9.2_db/)
!  real(kind=db),dimension(3),parameter:: cS_G=(/ 0._db, 6.5_db,6.5_db /)
!  real(kind=db),dimension(3),parameter:: cS_H=(/ 3.44_db,4.19_db ,7.63_db /)
!  real(kind=db),dimension(3),parameter:: cS_I=(/ 1.74_db,2.85_db,4.59_db/)
!  real(kind=db),dimension(3),parameter:: cS_L=(/1.63_db,3.65_db,5.28_db/)
!  real(kind=db),dimension(3),parameter:: cS_K=(/5.86_db,4.42_db,10.28_db/)
!  real(kind=db),dimension(3),parameter:: cS_M=(/4.55_db,3.98_db,8.53_db/)
!  real(kind=db),dimension(3),parameter:: cS_F=(/1.40_db,6.29_db,7.69_db/)
!  real(kind=db),dimension(3),parameter:: cS_S=(/3.68_db,3.95_db,7.63_db/)
!  real(kind=db),dimension(3),parameter:: cS_T=(/3.31_db,3.88_db,7.19_db/)
!  real(kind=db),dimension(3),parameter:: cS_W=(/2.74_db,4.55_db,7.29_db/)
!  real(kind=db),dimension(3),parameter:: cS_Y=(/2.78_db,6.52_db,9.3_db/)
!  real(kind=db),dimension(3),parameter:: cS_V=(/0.12_db,3.47_db,3.59_db/)
!  real(kind=db),dimension(3),parameter:: cS_P=(/0.0_db,0.0_db,0.0_db/)

!end module phenom_const
!********************************************

!module protdep_phenom_par
!  use defreal
!  implicit none
!  save
!
!  real(kind=db):: enpar(9,0:4)
!  real(kind=db),allocatable :: qfrei(:,:)
!  !     enpar(9,0:4) coefficients appearing in the different theermodyn functions
!  !     enpar(1,0:2) -->B(0:2)
!  !     enpar(2,0:3) -->u(0:3)
!  !     enpar(3,0:4) -->phi(0:4)
!  !     enpar(4,0:2) -->D_NP(0:2)
!  !     enpar(5,0:2) -->D_P(0:2)
!  !     enpar(6,0:3) -->v_NP(0:3)
!  !     enpar(7,0:3) -->v_P(0:3)
!  !     enpar(8,0:4) -->h_NP(0:4)
!  !     enpar(9,0:4) -->h_P(0:4)
!
!end module protdep_phenom_par


!********************************************

!!$module expdata
!!$  use defreal
!!$  implicit none
!!$  save
!!$  integer,parameter:: nfiles=3 !!NUMERO MASSIMO DI TECNICHE SPERIMENTALI
!!$  integer :: expflag(nfiles,2)
!!$  ! expflag(i,1)=0,1 depending on whether there is an experimental file for technique i (Calorimetry, FLuorescence, CD, ...)
!!$  !expflag(i,2)=1 if the denaturating agent for technique i is T
!!$  !expflag(i,2)=2 if the denaturating agent for technique i is [D]
!!$  !NOTICE: if more denaturants are used in different experiments, the dimension 2 can be increased.
!!$  integer,parameter :: ndatamax=1000 !max num of exp. points
!!$  integer :: ndata(nfiles)
!!$  real(kind=db)::TorDenFix(nfiles) !value of fixed T in chemical denaturation experiments, of of fixed [D] in thermal denaturation experiments.
!!$  real(kind=db),allocatable :: xexp(:,:),yexp(:,:),erry(:,:)
!!$end module expdata

!***************************************************************************
