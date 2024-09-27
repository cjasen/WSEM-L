subroutine stringhe(&
     !O:
     &m,S,&
     !I:
     &e,N)
  !modified by pierpaolo for speedup, eliminating T (12/8/2011) 

  !     m(i,j): stretch of 1 between i and j (included); S(i-1,j+1): stretch of 1 between i and j,
  !     with 0 in i-1 and j+1  
  use defreal
  implicit none
  integer:: N,i,j,k
 
  real (kind=db) S(0:N+1,0:N+1),e(3,N,N)
  real (kind=db) aux(N,0:N),U(N,0:N),V(N,0:N)! ,T(N-1,0:N-1,0:N)
  real (kind=db) Z
  real (kind=db) m(0:N+1,0:N+1)

  do j=1,N
     do i=0,j
        aux(j,i)=0.0_db
     enddo
  enddo
!!$  do i=1,N
!!$     do j=0,N
!!$        aux(i,j)=0.0_db
!!$     enddo
!!$  enddo
  
  do j=1,N
     if (j.gt.0) then
        do k=1,j-1
           aux(j,0)=aux(j,0)-e(1,k,j)
        enddo
     endif
     do i=1,j
        aux(j,i)=e(1,j,j)
        if (i.lt.j) then
           do k=1,j-i
              aux(j,i)=aux(j,i)-e(1,k,j)
           enddo
        endif
     enddo
     do i=0,j
        aux(j,i)=dexp(aux(j,i))
     enddo
  enddo
!  aux(j,i)=exp(\sum_{k=1}^{j-i} beta h_{k,j})-beta h_{j,j} 
  
!!$  do i=1,N
!!$     do j=0,N-1
!!$        do k=0,N
!!$           T(i,j,k)=0.0d0
!!$        enddo
!!$     enddo
!!$  enddo
!!$  
!!$  do j=1,N-1
!!$     do i=0,j
!!$        T(j,i,0)=aux(j,i)
!!$        T(j,i,i+1)=aux(j,i)
!!$     enddo
!!$  enddo

!!$  do j=1,N
!!$     do i=0,N
!!$        U(j,i)=0.0_db
!!$     enddo
!!$  enddo
!!$  do j=0,N
!!$     U(1,j)=1
!!$  enddo

  do j=1,N
     do i=0,j
        U(j,i)=0.0_db
     enddo
  enddo
  do j=0,1
     U(1,j)=1
  enddo
  
  do j=2,N
     do k=0,j-1
        U(j,0)=U(j,0)+U(j-1,k)*aux(j-1,k)
     enddo
     do i=1,j
        U(j,i)=U(j-1,i-1)*aux(j-1,i-1)
     enddo
  enddo
!!$  do j=2,N
!!$     do i=0,j
!!$        do k=0,j-1
!!$           U(j,i)=U(j,i)+U(j-1,k)*T(j-1,k,i)
!!$        enddo
!!$     enddo
!!$  enddo
!!$

  do j=1,N
     do i=0,j
        V(j,i)=aux(j,i)
     enddo
  enddo
   
  do j=N-1,1,-1
     do i=0,j
! OCCHIO! per via di approssimazioni numeriche, non e' esattamente uguale fare:
!      V(j,i)=aux(j,i)*V(j+1,0)+aux(j,i)*V(j+1,i+1)
! (che coincide con il modo di procedere della routine originale di Zamparo) oppure fare come qui sotto:
!            V(j,i)=V(j+1,0)+V(j+1,i+1)
!            V(j,i)=aux(j,i)*V(j,i)
!     oppure ancora :
!            V(j,i)=aux(j,i)*(V(j+1,0)+V(j+1,i+1))
        V(j,i)=V(j,i)*(V(j+1,0)+V(j+1,i+1))

!     questi ultimi due danno risultati uguali, ma differenti dal primo.
!     NB: da prove con mathematica, ne' il primo ne' il secondo-terzo danno riultati esatti, ma il secondo-terzo si avvicina di piu'. Pero' lo lascio come sta.
!c$$$            write(*,*) "j,i,0,addend=",j,i,0,aux(j,i)*V(j+1,0)
!c$$$            write(*,*) "j,i,i+1,addend,V(j,i)=",j,i,i+1,
!c$$$     >           aux(j,i)*V(j+1,i+1),V(j,i),aux(j,i),
!c$$$     >           V(j+1,0),V(j+1,i+1),aux(j,i)*(V(j+1,0)+V(j+1,i+1))     
     enddo
  enddo

!!$  do i=1,N
!!$     do j=0,N
!!$        V(i,j)=0.0d0
!!$     enddo
!!$  enddo
!!$  do j=0,N
!!$     V(N,j)=aux(N,j)
!!$  enddo
!!$
!!$  do j=N-1,1,-1
!!$     do i=0,j
!!$        do k=0,j+1
!!$           V(j,i)=V(j,i)+T(j,i,k)*V(j+1,k)
!!$        enddo
!!$     enddo
!!$  enddo
 
  Z=V(1,0)+V(1,1)
 ! write (*,*) "           stringhe:: Z=",Z
  do j=0,N+1
     do i=0,j
!     do j=0,N+1
        m(i,j)=0.0_db
     enddo
  enddo
  
  do j=1,N
     do i=1,j
        if(Z.ne.0.) then
           !AGGIUNTO QUESTO IF (16/4/2013)
           do k=j-i+1,j   
              m(i,j)=m(i,j)+U(j,k)*V(j,k)/Z
           enddo
        else
           m(i,j)=0.0_db
        endif
     enddo
  enddo
  
  do j=0,N+1
     do i=0,j
        if (j.eq.i) then
           S(i,j)=1.0d0-m(i,j)
        elseif (j.eq.i+1) then
           S(i,j)=1.0d0-m(i,j-1)-m(i+1,j)+m(i,j)
        else
           S(i,j)=m(i+1,j-1)-m(i,j-1)-m(i+1,j)+m(i,j)
        endif
     enddo
  enddo
  
  return
end subroutine stringhe


