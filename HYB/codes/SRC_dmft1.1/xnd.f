        subroutine xnd(xm,a0,da0,kmin,kmax)
c       reads hamilonian, overlap and self-energy matrices
	include 'dmft.dat'

	complex*16 E,T,work,omega
        real*8 tmp
        dimension T(nume,nume),Tinv(nume,nume)
	dimension ipiv(Nume),work(Lwork)

	beta_pi=pi/Beta
	om=(2*Im+1)*pi/Beta
	tol=0.1*emax/(om*om)

c.....asymptotics 
        wn=zero
	sig_trace=zero
        hs_trace=zero
	do k=kmin,kmax
        wn=wn+weight(k)
	 do i=1,M
	  hs_trace=hs_trace+dreal(h(i,i,k))*weight(k)
	 enddo
	 do jatom=1,natom_int
          Nint=n_int(jatom)
	  do i=1,Nint
	  sig_trace=sig_trace+
     &    dreal(etainf(i,i,jatom))
     &    *weight(k)
	  enddo
         enddo
	enddo

        tmp=0.d0
        do n=0,Im
         tmp=tmp+1.d0/(Two*n+One)**2
        enddo
        tail=Beta*(0.25d0-Two/(Pi*Pi)*tmp)
c.....end asymptotics

	a0=zero
        da0=zero
        do 400 k=kmin,kmax
         if (weight(k).gt.1d-6) then
          sum=zero
          sumd=zero
          do 300 n=0,Im
	   omega=(n+n+1)*beta_pi*cim
	   do j=1,M
	    do i=1,M
	     T(i,j)=-H(i,j,k)
	    enddo
	    T(j,j)=T(j,j)+xm+omega
	   enddo
c......add selfenergy
           do 200 jatom=1,natom_int
	    index0=n_min(jatom)-1
	    Nint=n_int(jatom)
	    do j=1,Nint
             do i=1,Nint
              T(index0+i,index0+j)=T(index0+i,index0+j)
     &        -eta(i,j,n,jatom)
             enddo
            enddo
 200       continue	
c.......end add selfenergy
           call ZGETRF(M,M,T,Nume,ipiv,info)
	   call ZGETRI(M,T,Nume,ipiv,work,lwork,info)
	   do i=1,M
	    sum=sum+T(i,i)
            do j=1,M
             sumd=sumd+dreal(T(i,j)*T(j,i))
            enddo
           enddo
 300      continue
	  a0=a0+Two*sum*weight(k)
          da0=da0-Two*sumd*weight(k)
         endif
 400    continue
	a0=a0/Beta+real(M)*wn/Two
     &  +(M*xm*wn-hs_trace-sig_trace)*tail
        da0=da0/Beta+M*tail*wn
 	a0=spinfac*a0
        da0=spinfac*da0
	return
	end
