module energy_profiles
    implicit none

contains

    subroutine profiles_ij(xi_j,N_interval,e,&
                           nu_ijm)
        use defreal
        use protdep_par,only:N
        implicit none
        ! xi_j are the same as the calculated in calc_thermo_ab
        ! creo que no necesito "e" si saco las zetas de calc_thermo_ab, ¿o mejor calcularlo todo aquí? ¿Las A's son las mismas? Para el cálculo sólo necesito las w's y las xi's





    end subroutine profiles_ij

    !********************************
    !********************************
    !********************************
    !********************************
    !********************************

    subroutine profiles(nu_ijm,xi_j,N_inverval,&
                        Z_M)
        use defreal
        use protdep_par,only:N
        implicit none
        ! nu^ij_m is what I have called psi in the notes and it's calculated in profiles_ij. Probably, we would have to change it a bit to work because of the normalization
        ! xi_j are the same as the calculated in calc_thermo (not in calc_thermo_ab)
        ! N_interval is the arbitrary interval of residues we choose as variable
        ! Z_M := Z(sum m_k*sigma_k=M) which is an array from 0 to N. We could also calculate Z_L. F/RT=-log(Z_M) )

        real(kind=db),intent(in):: nu_ijm(:,:,:), chi_j(:), N_interval(:)
        real(kind=db),intent(out):: Z_M(:)
        real(kind=db) :: A_ijMm(1:N+1,1:N,0:N,0:M), sum_im ! don't know if A_ijMm(1:N+1,...) that N+1 is ok
        integer(kind=db) :: delta_K, i, j, k, m, M_ind

        A_ijMm=0._db 
        !A(1,1,0,0) = 1._db
        Z_M=1._db
        sum_im=0

        do M_ind=0,N
         do j=1,N 
            if(j==M) then !give value to Kronecker's delta
               delta_K=1
            else
               delta_K=0
            endif

            sum_im=
            do i=1,j
               do m=max(M_ind-i+2,0),min(j-i+1,M_ind) !here, we have to sum A_ijMm
                 ! A_ijMm(i,j,M_ind,m)=
               enddo !m
            enddo !i

            Z_M(M_ind)=(1-delta_k)*Z_M(M_ind)/xi_j(j)+sum_im !calculate Z_j(M) with Z_{j-1}(M)

         enddo !j
      end do !M  

    end subroutine profiles

    !********************************
    !********************************
    !********************************
    !********************************
    !********************************

    subroutine profiles_old(&
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
   end subroutine profiles_old

end module energy_profiles