	subroutine spec(xm,Aw,W,kmin,kmax)
        include 'dmft.dat'
	complex*16 T(Nume,Nume),work(lwork)
        complex*16 GLT(nume,-nw:nw)
        real*8     Aw(nume,-nw:nw)
	complex*16 sum,omega
	dimension ipiv(nume)

	beta_pi=pi/Beta
        omax=(2*Im+1)*beta_pi
c-------GREEN FUNCTION - OMEGA LOOP
        dw=w/real(nw)
	do n=-nw,nw
	 do j=1,M
	  aw(j,n)=zero
         enddo
        enddo
c********k-point summation
	do 100 k=kmin,kmax
         do 200 n=-nw,nw
          omega=n*dw+(5.d-2)*cim
	  do j=1,M
	   do i=1,M
	    T(i,j)=-H(i,j,k)  
	   enddo
           T(j,j)=T(j,j)+omega+xm
	  enddo
c.........inversion
          call ZGETRF(M,M,T,Nume,ipiv,info)
          call ZGETRI(M,T,Nume,ipiv,work,lwork,info)
	  do i=1,M
	   aw(i,n)=aw(i,n)-imag(T(i,i))*weight(k)/Pi
	  enddo
 200     continue
 100    continue
        end 


