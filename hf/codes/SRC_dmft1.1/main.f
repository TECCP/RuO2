	program main
	include 'dmft.dat'
        include 'mpif.h'
        external cputim
        double precision my_x,xmu0,dxmu
	character*80 aline,tmp,input_data
        complex*16 xocc(nume,nume),xocct(nume,nume)     
        dimension Aw(nume,-nw:nw),Awt(nume,-nw:nw)
        integer size
        logical fix,licht,newt


        CALL MPI_Init( ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
        Pi=acos(-1.d0)
        diag_only=.false.
        fix=.false.
        licht=.false.
        newt=.false.

        iposition=0
        size=0
        in=5
        call MPI_PACK_SIZE(in,MPI_LOGICAL,MPI_COMM_WORLD,
     &  iout,ierr)
        size=size+iout
        in=6+2*NAT
        call MPI_PACK_SIZE(in,MPI_INTEGER,MPI_COMM_WORLD,
     &  iout,ierr)
        size=size+iout
        in=5
        call MPI_PACK_SIZE(in,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,
     &  iout,ierr)
        size=size+iout

c-----------read myid=0		
        if (myid.eq.0) then
	 open(unit=10,file='hamilt',form='formatted',status='old')
         open(unit=5,file='dmft.in',form='formatted',status='old')
         open(unit=6,file='dmft.out',form='formatted',status='unknown')
	 open(unit=12,file='sigma',form='formatted',status='unknown')
c        open(unit=15,file='gw',form='formatted',status='unknown')
c        open(unit=16,file='gtau',form='formatted',status='unknown')
         open(unit=23,file='n.dmft',form='formatted',status='unknown')
c        open(unit=17,file='delta.w',form='formatted',status='unknown')
c        open(unit=18,file='delta.tau',form='formatted',
c    &   status='unknown')
!******Band input*******************************
	 read(10,*)kpt,M
         if (kpt.gt.nkpt) stop 'NKPT too small'
         if (m.gt.nume) stop 'NUME too small'
	 wsum=Zero
	 do k=1,kpt
	  read(10,*)weight(k)
 	  wsum=wsum+weight(k)
	  do i=1,M
	   do j=1,M
	    read(10,*)HR,HI
	    H(i,j,k)=HR+cim*HI
	   enddo
	  enddo
          do i=1,M
           do j=i,M
            diff=cdabs(h(i,j,k)-dconjg(h(j,i,k)))
            if (diff.gt.1d-5) then
             write(*,*)'NON_HERMITEAN:',i,j,k,h(i,j,k),h(j,i,k)
             STOP 'NON HERMITEAN HAMILTONIAN'
            endif
           enddo
          enddo
	 enddo
	 do k=1,kpt
	  weight(k)=weight(k)/wsum
	 enddo
        
!******Band input - end ***********************
!******Local input*****************************
	 read(5,*)beta,Ecut,LL
         read(5,'(a80)')aline
         if(index(aline(1:72),'fix').ne.0) fix=.true.
         if(index(aline(1:72),'diag').ne.0) diag_only=.true.
         if(index(aline(1:72),'licht').ne.0) licht=.true.
          if(index(aline(1:72),'newt').ne.0) newt=.true.
         read(aline,*)x_el,xmu,emax,dum,spinfac 
         if (.not.fix) then
c.......looks for coverged xmu from previous iteration
          dxmu=1.
 110      read(6,'(a80)',end=100)aline
          write(tmp,'(a80)')aline
          if (aline(1:3).eq.':MU') then
           read(tmp,'(3x,2f12.7)')xmu0,dxmu
           xmu=xmu0
           dxmu=4*dxmu
           dxmu=max(0.001,dxmu)
           dxmu=min(1.,dxmu)
          endif
          goto 110
 100      rewind(6)
         endif
c...................................................
         read(5,*)MM,natom_int
         if (natom_int.gt.nat) stop 'NAT too small'
         if (.not.MM.eq.M) stop 'MM.neq.M'
         do jatom=1,natom_int
	  if (licht) then
           read(5,*)n_min(jatom),n_int(jatom)
	   edc(jatom)=0.
	  else
	   read(5,*)n_min(jatom),n_int(jatom),edc(jatom)
	  endif
          if (n_int(jatom).gt.niter) stop 'NITER too small'
          do i=1,spinfac*n_int(jatom)
           read(5,*)
          enddo
         enddo
         write(6,*)'HAMILTONIAN STRUCTURE'      
         write(6,*)'MPI-buffer size',size
         write(6,*)'number of k-points:',kpt
         write(6,*)'size:',M,' x',M
         write(6,*)'interacting blocks:'
         write(6,*)(n_min(i),'-',n_min(i)+n_int(i)-1,' ',i=1,natom_int)
         Im=int(Ecut*Beta/(2.d0*Pi))
         write(6,*)'GREEN FUNCTION PARAMETERS'
	 write(6,*)'Beta=',beta
         write(6,*)'Energy cut-off:',Ecut
	 write(6,*)'Number of Matsubara frequencies=',
     &             Im
         write(6,*)
         if (fix) write(6,*)'chem. potential is fixed at',xmu
         if (diag_only) write(6,*) 'ONLY DIAGONAL ELEMENTS OF 
     &     SIGMA CONSIDERED'
         if (Im.gt.momega) 
     &   stop 'Number of Matsubara frequencies too large'
!********Self-energy input**********************
	 read(12,*,iostat=ios)
         if (ios.eq.0) then
	  rewind(12)
          do 10 jatom=1,natom_int
	  do 10 i=1,n_int(jatom)
	  do 10 j=1,n_int(jatom)
	  do 10 n=0,Im
 	   read(12,*)dum,SR,SI
 	   eta(i,j,n,jatom)=SR+cim*SI 
 10	  continue
         else
          write(6,*)'ZERO SELF-ENERGY ASSUMED'
          write(6,*)
          do 13 jatom=1,natom_int
          do 13 i=1,n_int(jatom)
          do 13 j=1,n_int(jatom)
          do 13 n=0,Im
           eta(i,j,n,jatom)=(0.d0,0.d0)
 13       continue
         endif
!**********Self-energy assymptotics*************
c------end read myid=0
c----------MPI broadcast------------------------------------
        call MPI_PACK(fix,1,MPI_LOGICAL,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(diag_only,1,MPI_LOGICAL,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(licht,1,MPI_LOGICAL,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(newt,1,MPI_LOGICAL,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
c    
        call MPI_PACK(M,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(spinfac,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(natom_int,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(n_min,NAT,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(n_int,NAT,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(Im,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(LL,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
        call MPI_PACK(Kpt,1,MPI_INTEGER,input_data,size,iposition,
     &  MPI_COMM_WORLD,ierr)
      
        call MPI_PACK(Beta,1,MPI_DOUBLE_PRECISION,input_data,size,
     &  iposition,MPI_COMM_WORLD,ierr)
        call MPI_PACK(xmu,1,MPI_DOUBLE_PRECISION,input_data,size,
     &  iposition,MPI_COMM_WORLD,ierr)
        call MPI_PACK(dxmu,1,MPI_DOUBLE_PRECISION,input_data,size,
     &  iposition,MPI_COMM_WORLD,ierr)
        call MPI_PACK(x_el,1,MPI_DOUBLE_PRECISION,input_data,size,
     &  iposition,MPI_COMM_WORLD,ierr)
        call MPI_PACK(emax,1,MPI_DOUBLE_PRECISION,input_data,size,
     &  iposition,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(input_data,size,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
     
        call MPI_BCAST(weight,Nkpt,MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edc,Nat,MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr)
        
        nsize=Niter*Niter*Nat*(momega+1)
        call MPI_BCAST(eta,nsize,MPI_DOUBLE_COMPLEX,0,
     &  MPI_COMM_WORLD,ierr)
        
        nsize=Nume*Nume*Nkpt
        call MPI_BCAST(H,nsize,MPI_DOUBLE_COMPLEX,0,
     &  MPI_COMM_WORLD,ierr)
c-----------end myid=0
       else
        call MPI_BCAST(input_data,size,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    
        call MPI_UNPACK(input_data,size,iposition,fix,1,MPI_LOGICAL,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,diag_only,1,
     &  MPI_LOGICAL,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,licht,1,MPI_LOGICAL,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,newt,1,MPI_LOGICAL,
     &  MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(input_data,size,iposition,M,1,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,spinfac,1,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,natom_int,1,
     &  MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,n_min,NAT,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,n_int,NAT,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,Im,1,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,LL,1,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,Kpt,1,MPI_INTEGER,
     &  MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(input_data,size,iposition,Beta,1,
     &  MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,xmu,1,
     &  MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,dxmu,1,
     &  MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,x_el,1,
     &  MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(input_data,size,iposition,emax,1,
     &  MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(weight,Nkpt,MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr)
        call MPI_BCAST(edc,Nat,MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr)

        nsize=Niter*Niter*Nat*(momega+1)
        call MPI_BCAST(eta,nsize,MPI_DOUBLE_COMPLEX,0,
     &  MPI_COMM_WORLD,ierr)       

        nsize=Nume*Nume*Nkpt
        call MPI_BCAST(H,nsize,MPI_DOUBLE_COMPLEX,0,
     &  MPI_COMM_WORLD,ierr)        
       endif
c----------end MPI braodcast--------------------------------
        call asym(licht)
        kp=kpt/numprocs
        ir=kpt-kp*numprocs
        kt=1
        do ip=1,ir
         nk1(ip)=kt
         nk2(ip)=kt+kp
         kt=kt+kp+1
        enddo
        do ip=ir+1,numprocs
         nk1(ip)=kt
         nk2(ip)=kt+kp-1
         kt=kt+kp
        enddo
          kmin=nk1(myid+1)
          kmax=nk2(myid+1)

!******Chem. potential calculation**************
         call cputim(t0)
  	 if (fix) then
c-----I/O myid=0
          if (myid.eq.0) then
	   write(6,*)'CHEM. POTENTIAL FIXED'
           write(6,*)
          endif
c----------------
	  call xn(xmu,my_x,kmin,kmax)
          call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &    MPI_SUM,0,MPI_COMM_WORLD,ierr)
c-------I/O myid=0
          if (myid.eq.0) then
           write(6,*)
           write(6,'(a3,2f12.7)')':MU',xmu,abs(xmu-xmu0)
            write(6,'(a7,f12.7)')':NTOTAL',x
          endif
c------------
         else
          if (myid.eq.0) then
           write(6,*)'OCCUPATION FIXED  ',x_el
           write(6,*)
           if (newt) then 
            write(6,*)'mu search: newton'
           else
            write(6,*)'mu search: bisection'
           endif
          endif
          
          if (newt) then
           call mu_calc_newton(xmu,x,iout)
           if (.not.(iout.eq.0)) then
            if (myid.eq.0) write(6,*)'newton failed'
            newt=.false.
           endif
          endif
          if (.not.newt) then
           call mu_calc(xmu,dxmu,x)
          endif
          if (myid.eq.0) then
           write(6,'(a3,2f12.7)')':MU',xmu,abs(xmu-xmu0)
           write(6,'(a7,f12.7)')':NTOTAL',x
          endif
         endif
        write(6,*)
        call cputim(t1)
!*********Bath GF for each interacting site     
c       reset weights
         do k=1,kpt
c         weight(k)=One/real(kpt)
         enddo
  	call gloc(xmu)
        call delta(xmu)
        call cputim(t2)
!*********Occupation matrix************************
c       call nmat(xmu,xocc,kmin,kmax)
c       size=nume*nume
c       call MPI_REDUCE(xocc,xocct,size,MPI_DOUBLE_COMPLEX,
c    &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call spec(xmu,Aw,10.d0,kmin,kmax)
        size=nume*(2*nw+1)
        call MPI_REDUCE(Aw,Awt,size,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (myid.eq.0) then
         open(25,file='aw.dmft',status='unknown',form='formatted')
         write(6,*)'***********************************************'
         write(6,*)'SPECTRAL FUNCTION CALCULATION'
         write(6,*)'***********************************************'
         write(6,*)
         dw=10.d0/nw
         do i=1,M
          write(25,*)'#',i
          do n=-nw,nw
           write(25,*)n*dw,Awt(i,n)
          enddo
          write(25,*)
         enddo
        endif
        
        call energy(xmu,xocc,ee,eu,kmin,kmax)
        call MPI_REDUCE(ee,eet,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(eu,eut,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        size=nume*nume
        call MPI_REDUCE(xocc,xocct,size,MPI_DOUBLE_COMPLEX,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (myid.eq.0) then
         write(6,*)'OCCUPATION MATRIX:'
         if (M.le.10) then
          write(6,*)'real:'
          do i=1,M
           write(6,123)(dreal(xocct(i,j)),j=1,M)
          enddo
          write(6,*)
          write(6,*)'imaginary:'
          do i=1,M
           write(6,123)(imag(xocct(i,j)),j=1,M)
          enddo
         else
          do jatom=1,natom_int
           nint=n_int(jatom)
           ind=n_min(jatom)
           write(*,*)'atom #',jatom
           write(6,*)'real:'
           do i=ind,ind+nint-1
            write(6,123)(dreal(xocct(i,j)),
     &      j=ind,ind+nint-1)
           enddo
           write(6,*)
           write(6,*)'imaginary:'
           do i=ind,ind+nint-1
            write(6,123)(imag(xocct(i,j)),
     &      j=ind,ind+nint-1)
           enddo
           write(6,*)
          enddo
          write(6,*)'DIAGONAL ELEMENTS'
          do i=1,M
           write(*,*)i,dreal(xocct(i,i))
          enddo
         endif
         trace=0.d0
         do i=1,M
          do j=1,M
           write(23,*)dreal(xocct(i,j)),imag(xocct(i,j))
          enddo
          trace=trace+dreal(xocct(i,i))
         enddo
 123     format(10f14.8)
         write(6,*)
         write(6,*)'spinfac*trace=',spinfac*trace
         write(6,*)
         write(6,*)'<E_kin>=',EEt
         write(6,*)'<U>=',EUt
         write(6,*)'E_tot=',EEt+EUt
        endif
        call cputim(t3)
        if (myid.eq.0) then
         write(6,*)
         write(6,*)'CPU-summary:'
         write(6,*)'mu_calc:',t1-t0
         write(6,*)'gloc:',t2-t1
         write(6,*)'ene:',t3-t2
         write(6,*)'total:',t3-t0
        endif
        call MPI_FINALIZE(ierr)        

	end

	
	   
