module profiles
    use defreal
    use protdep_par,only:N
    implicit none

contains

    subroutine profiles_ab(e, withmprofile, beF, nu)
        use defreal
        use protdep_par,only:N
        implicit none
       
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


    
    end subroutine profiles_ab


    !*************************
    !*************************
    !*************************
    !*************************

    subroutine profiles_2()
        ! Código de la segunda subrutina
    end subroutine profiles_2

    !*************************
    !*************************
    !*************************
    !*************************

    subroutine profiles_old(e, withmprofile, beF, nu) 
         use defreal
         use protdep_par,only:N
         implicit none
        
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




end module profiles