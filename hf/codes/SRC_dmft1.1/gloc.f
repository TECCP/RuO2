	subroutine gloc(xmu)
        include 'dmft.dat'
        include 'mpif.h'
        
	complex*16 T(Nume,Nume),work(lwork),Hk(nume,nume)
        complex*16 Hk2(nume,nume),T2(nume,nume)
        complex*16 Hkt(nume,nume),Hk2t(nume,nume)
        complex*16 Glo(nume,nume,0:momega),Glt(nume,nume,0:momega)
	complex*16 GL(Niter,Niter),Gw(Niter,Niter,0:momega,nat)
        complex*16 xm1(Niter,Niter,nat),xm2(Niter,Niter,nat)
        complex*16 xm1t(Niter,Niter,nat),xm2t(Niter,Niter,nat)
        complex*16 tmp2(Niter,Niter),xm0(Niter,Niter)
	complex*16 sum,omega
        complex*16 ym1(Niter,Niter),ym2(Niter,Niter)
        integer size
	dimension ipiv(nume),shift(niter,nat),shiftt(niter,nat)

        if (myid.eq.0) then
c	 open(20,file='mom1.tmp',status='unknown',form='formatted')
c	 open(21,file='mom2.tmp',status='unknown',form='formatted')
c        open(22,file='shift.tmp',status='unknown',form='formatted')
         write(6,*)'***********************************************'
         write(6,*)'GREEN FUNCTION CALCULATION'
         write(6,*)'***********************************************'
         write(6,*)
        endif

	beta_pi=pi/Beta
        omax=(2*Im+1)*beta_pi
        do i=1,M
         do j=1,M
          Hk(i,j)=czero
          Hk2(i,j)=czero
         enddo
        enddo
        do j=1,Niter
         do i=1,Niter
          xm0(i,j)=zero
         enddo
         xm0(j,j)=one
        enddo

        kmin=nk1(myid+1)
        kmax=nk2(myid+1)

c.........preparation for moment calculation
        do 111 k=kmin,kmax
         do j=1,M
          do i=1,M
           T(i,j)=H(i,j,k)
          enddo
          T(j,j)=T(j,j)-xmu
         enddo
         do jatom=1,natom_int
          Nint=n_int(jatom)
          index0=n_min(jatom)-1
          do j=1,Nint
           do i=1,Nint
            T(index0+i,index0+j)=T(index0+i,index0+j)+
     &      etainf(i,j,jatom)
           enddo
          enddo
         enddo
         do j=1,M
          do i=1,M
           sum=czero
           do l=1,M
            sum=sum+T(i,l)*T(l,j)
           enddo
           T2(i,j)=sum
          enddo
         enddo
c........Hk=sum_k (H(k)-mu+eta(inf))
c........Hk2=sum_k (H(k)-mu+eta(inf))^2
         do j=1,M
          do i=1,M
           Hk(i,j)=Hk(i,j)+T(i,j)*weight(k)
           Hk2(i,j)=Hk2(i,j)+T2(i,j)*weight(k)
          enddo
         enddo
 111     continue  
         size=nume*nume
         call MPI_REDUCE(Hk,Hkt,size,MPI_DOUBLE_COMPLEX,
     &   MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(Hk2,Hk2t,size,MPI_DOUBLE_COMPLEX,
     &   MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(Hkt,size,MPI_DOUBLE_COMPLEX,0,
     &   MPI_COMM_WORLD,ierr)
         call MPI_BCAST(Hk2t,size,MPI_DOUBLE_COMPLEX,0,
     &   MPI_COMM_WORLD,ierr)
         
c.......moments of local & bath green functions
c.......ym1, ym2 1st and 2nd moment of the local spectal function
c.......xm1, xm2 1st and 2nd moment of the local bath function
       do jatom=1,natom_int
        Nint=n_int(jatom)
        index0=n_min(jatom)-1
         do j=1,Nint
          do i=1,Nint
           ym1(i,j)=Hkt(index0+i,index0+j)
           ym2(i,j)=eta1(i,j,jatom)+Hk2t(index0+i,index0+j)
           xm1(i,j,jatom)=ym1(i,j)-etainf(i,j,jatom)
          enddo
          shift(j,jatom)=dreal(xm1(j,j,jatom))
         enddo
         do j=1,Nint
          do i=1,Nint
           sum=czero
            do k=1,Nint
             sum=sum+xm1(i,k,jatom)*xm1(k,j,jatom)-
     &               ym1(i,k)*ym1(k,j)
            enddo
           xm2(i,j,jatom)=sum+ym2(i,j)-eta1(i,j,jatom)
          enddo
         enddo
        enddo

c********k-point summation
        do iw=0,Im
         do j=1,M
          do i=1,M
           glo(i,j,iw)=czero
          enddo
         enddo
        enddo

	do 100 k=kmin,kmax
         do 200 iw=0,Im
          omega=(2*iw+1)*beta_pi*cim
	  do j=1,M
	   do i=1,M
	    T(i,j)=-H(i,j,k)
	   enddo
           T(j,j)=omega+xmu+T(j,j)
	  enddo
c.........add -self-energy 
          do jatom=1,natom_int
           index0=n_min(jatom)-1
           Nint=n_int(jatom)
	   do j=1,Nint
	    do i=1,Nint
	     T(index0+i,index0+j)=T(index0+i,index0+j)
     &       -eta(i,j,iw,jatom)
	    enddo
	   enddo
          enddo
c.........invertion
          call ZGETRF(M,M,T,Nume,ipiv,info)
          call ZGETRI(M,T,Nume,ipiv,work,lwork,info)
	  do j=1,M
	   do i=1,M
	    Glo(i,j,iw)=Glo(i,j,iw)+T(i,j)*weight(k)
           enddo
	  enddo
 200     continue
 100    continue
c**********
        size=nume*nume*(momega+1)
        call MPI_REDUCE(Glo,Glt,size,MPI_DOUBLE_COMPLEX,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)

c-------rest is done on single CPU
        if (myid.eq.0) then
c........local (small size) GF for each interacting site
         do 250 jatom=1,natom_int
          Nint=n_int(jatom)
          index0=n_min(jatom)-1
          do 260 iw=0,Im
           do j=1,Nint
            do i=1,Nint
             GL(i,j)=GLT(index0+i,index0+j,iw)
            enddo
           enddo
           call ZGETRF(Nint,Nint,GL,Niter,ipiv,info)
           call ZGETRI(Nint,GL,Niter,ipiv,work,lwork,info)	 
c.........add self-energy to G^-1 -> G_bath^-1 
           do j=1,Nint
	    do i=1,Nint
	     GL(i,j)=GL(i,j)+eta(i,j,iw,jatom)
             if (i.eq.j) GL(i,j)=GL(i,j)+shift(i,jatom)
	    enddo
	   enddo
c.........inversion
           call ZGETRF(Nint,Nint,GL,Niter,ipiv,info)
           call ZGETRI(Nint,GL,Niter,ipiv,work,lwork,info)	
c........Gw - bath GF for each omega and interacting site
           do j=1,Nint
            do i=1,Nint
             Gw(i,j,iw,jatom)=GL(i,j)
            enddo
           enddo
 260      continue
 250     continue
c---------
         do 400 jatom=1,natom_int
          Nint=n_int(jatom)
c........calculate moments after shift
          do j=1,nint
           do i=1,nint
            sum=czero
            do k=1,nint
             sum=sum+xm1(i,k,jatom)*xm1(k,j,jatom)
            enddo
            xm2(i,j,jatom)=xm2(i,j,jatom)-sum
           enddo
          enddo
          do i=1,nint
           xm1(i,i,jatom)=xm1(i,i,jatom)-shift(i,jatom)
          enddo
          do i=1,nint
           do j=1,nint
            sum=czero
            do k=1,nint
             sum=sum+xm1(i,k,jatom)*xm1(k,j,jatom)
            enddo
            xm2(i,j,jatom)=xm2(i,j,jatom)+sum
           enddo
          enddo
c.......diag_only filter
          if (diag_only) then
           do i=1,nint
            do j=1,nint
             if (.not.i.eq.j) then
              xm1(i,j,jatom)=0.d0
              xm2(i,j,jatom)=0.d0
              do iw=0,Im
               Gw(i,j,iw,jatom)=0.d0
              enddo
             endif
            enddo
           enddo
          endif
c.........print moments
c	  do i=1,nint
c          write(22,*)shift(i,jatom)
c	   do j=1,nint
c	    write(20,*)dreal(xm1(i,j,jatom)),imag(xm1(i,j,jatom))
c	    write(21,*)dreal(xm2(i,j,jatom)),imag(xm2(i,j,jatom))
c...........print Gw
c	    do iw=0,Im
c	     w=(2*iw+1)*beta_pi
c	     write(15,*)w,dreal(gw(i,j,iw,jatom)),imag(gw(i,j,iw,jatom))
c	    enddo
c          enddo
c	  enddo
c.........Fourier transform w->tau
          iout=16
	  call w2t1(Gw(1,1,0,jatom),jatom,Nint,xm0,
     &    xm1(1,1,jatom),xm2(1,1,jatom),iout)
 400     continue         
        endif
        end 


