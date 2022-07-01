        subroutine energy(xm,a,EE,EU,kmin,kmax)
c    Orthogonal basis only !!!!!!!!!!!!!!!!!!!!!!!!!
	include 'dmft.dat'

	complex*16 E,T,work,omega,sum,a0(nume,nume)
        complex*16 a(nume,nume)
        complex*16 Sk
        dimension T(nume,nume),Sk(nume,nume)
	dimension ipiv(Nume),work(Lwork)
        dimension delta(nume,nume)

        do i=1,M
         do j=1,M
          Sk(i,j)=czero
          a(i,j)=czero
         enddo
        enddo
        EE=0.d0
        EU=0.d0
	beta_pi=pi/Beta
	om=(2*Im+1)*beta_pi

c.....unwrap Sigma(w_max) 
        do jatom=1,natom_int
         index0=n_min(jatom)-1
         Nint=n_int(jatom)
         do i=1,Nint
          do j=1,Nint
           Sk(i+index0,j+index0)=etainf(i,j,jatom)
          enddo
         enddo
        enddo
c.....asymptotic sum
        tail=0.d0
        do iw=0,Im
         tail=tail+1.d0/(2*iw+One)**2
        enddo
        tail=Beta*(0.25d0-Two/(Pi*Pi)*tail)
c.....k-point sum
        do 400 k=kmin,kmax
         do j=1,M
          do i=1,M
           a0(i,j)=czero
          enddo
         enddo
         eku=0.d0
c.....w_n sum
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
             T(index0+i,index0+j)=
     &       T(index0+i,index0+j)-eta(i,j,iw,jatom)
            enddo
           enddo
          enddo
c.......matrix inversion
          call ZGETRF(M,M,T,Nume,ipiv,info)
	  call ZGETRI(M,T,Nume,ipiv,work,lwork,info)
	  do j=1,M
           do i=1,M
	    a0(i,j)=a0(i,j)+(T(i,j)+dconjg(T(j,i)))
           enddo
          enddo
          do jatom=1,natom_int
           index0=n_min(jatom)-1
           Nint=n_int(jatom)
           do j=1,Nint
            do i=1,Nint
             eku=eku+dreal(T(index0+i,index0+j)*eta(j,i,iw,jatom))
            enddo
           enddo
          enddo
 300      continue
c....end w_n sum
          eku=eku/Beta
          do j=1,M
           do i=1,M
	    a0(i,j)=a0(i,j)/Beta-(H(i,j,k)+Sk(i,j))*tail
            eku=eku-dreal((H(i,j,k)+Sk(i,j))*Sk(j,i))*tail/2
           enddo
           a0(j,j)=a0(j,j)+xm*tail+0.5d0
           eku=eku+(xm*tail+0.5d0)*dreal(Sk(j,j))/2
          enddo
          do j=1,M
           do i=1,M
            EE=EE+dreal(H(j,i,k)*a0(i,j))*weight(k)
            a(i,j)=a(i,j)+a0(i,j)*weight(k)
           enddo
          enddo
          EU=EU+eku*weight(k)
 400    continue
c...end k-sum
        EE=spinfac*EE
        EU=spinfac*EU
	return
	end
