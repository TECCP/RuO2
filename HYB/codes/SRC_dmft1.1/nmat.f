        subroutine nmat(xm,a0,kmin,kmax)
c    Orthogonal basis only !!!!!!!!!!!!!!!!!!!!!!!!!
	include 'dmft.dat'

	complex*16 HS,E,T,work,omega,sum,a0(nume,nume)
        complex*16 Hk,Sk
        real*8 tmp,delta,trace
        dimension T(nume,nume),Tinv(nume,nume)
        dimension Sk(nume,nume),Hk(nume,nume)
	dimension ipiv(Nume),work(Lwork)
        dimension delta(nume,nume)

        do j=1,M
         do i=1,M
          delta(i,j)=0.d0
         enddo
         delta(j,j)=One
        enddo

	beta_pi=pi/Beta
	om=(2*Im+1)*pi/Beta
c.....asymptotics for sigma
c..... Tr [S ReSigma]
        do j=1,M
         do i=1,M
          Hk(i,j)=czero
          Sk(i,j)=czero
         enddo
        enddo
        do jatom=1,natom_int
         index0=n_min(jatom)-1
         Nint=n_int(jatom)
         do i=1,Nint
          do j=1,Nint
           Sk(i+index0,j+index0)=etainf(i,j,jatom)
          enddo
         enddo
        enddo
	do k=kmin,kmax
	 do j=1,M
          do i=1,M
	   Hk(i,j)=Hk(i,j)+H(i,j,k)*weight(k)
          enddo
	 enddo
	enddo
        tail=0.d0
        do iw=0,Im
         tail=tail+1.d0/(2*iw+One)**2
        enddo
        tail=Beta*(0.25d0-Two/(Pi*Pi)*tail)
c.....end asymptotics
        do j=1,M
         do i=1,M
	  a0(i,j)=czero
         enddo
        enddo
        wn=zero
        do 400 k=kmin,kmax
        wn=wn+weight(k)
         do 300 iw=0,Im
	  omega=(2*iw+1)*beta_pi*cim
	  do j=1,M
	   do i=1,M
	    T(i,j)=-H(i,j,k)
	   enddo
	   T(j,j)=T(j,j)+xm+omega
	  enddo
c......add selfenergy
          do jatom=1,natom_int
	   index0=n_min(jatom)-1
	   Nint=n_int(jatom)
	   do j=1,Nint
            do i=1,Nint
             T(index0+i,index0+j)=T(index0+i,index0+j)-eta(i,j,iw,jatom)  
            enddo
           enddo
          enddo
c.......end add selfenergy
          call ZGETRF(M,M,T,Nume,ipiv,info)
	  call ZGETRI(M,T,Nume,ipiv,work,lwork,info)
	  do j=1,M
           do i=1,M
	    a0(i,j)=a0(i,j)+(T(i,j)+dconjg(T(j,i)))*weight(k)
           enddo
          enddo
 300     continue
 400    continue
        do j=1,M
         do i=1,M
	  a0(i,j)=a0(i,j)/Beta
     &    +(xm*delta(i,j)*wn-Hk(i,j)-Sk(i,j)*wn)*tail
         enddo
         a0(j,j)=a0(j,j)+0.5d0*wn
        enddo
	return
	end
