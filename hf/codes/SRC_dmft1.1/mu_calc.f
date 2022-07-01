	subroutine mu_calc(xmu,tinc,x)
        include 'dmft.dat'
        include 'mpif.h'
        double precision my_x
        logical cont

        if (myid.eq.0) then
         write(6,*)'****************************************'
         write(6,*)'CALULATION OF CHEM. POTENTTIAL'
         write(6,*)'****************************************'
         write(6,*)'SEARCH FOR UPPER and LOWER BOUNDS'
        endif
c--------initial occupation
        kmin=nk1(myid+1)
        kmax=nk2(myid+1)
        
        call xn(xmu,my_x,kmin,kmax)
        call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
c
        if (myid.eq.0) then
         if (x.gt.x_el) tinc=-tinc
         x0=x
         t=xmu+tinc
        endif
c
        iter=1
 11     continue
        iter=iter+1
        call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call xn(t,my_x,kmin,kmax)
        call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,err)
        
        if (myid.eq.0) then
         if ((x-x_el)*(x0-x_el).gt.0.d0) then
          t=t+tinc
          x0=x
          cont=.true.
         else
          cont=.false.
          t0=t-tinc
         endif
        endif
        call MPI_BCAST(cont,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (cont) then
         goto 11
        else
         if (myid.eq.0) then
          t1=min(t0,t)
          t2=max(t0,t)
          write(6,*)
          write(6,*)'XMU BOUNDS:',t1,t2,' with ',iter,' xn calls'
          write(6,*)
         endif
        endif
         
	iter=0
        om=(2*Im+1)*pi/beta
        Xtol=1.d-1*(emax/om)**3
        Xtol=max(Xtol,1.d-6)
        if (myid.eq.0) write(6,*)'n_tot toletance =',xtol

 10     continue
        if (myid.eq.0) then
         t=(t1+t2)/Two
        endif
        call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call xn(t,my_x,kmin,kmax)
        call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,err)
        if (myid.eq.0) then
         iter=iter+1
	 if (abs(x_el-x).le.Xtol) then
	  cont=.false.
	 else
	  if (x.gt.x_el) then
	   t2=t
	   x2=x
	  else
           t1=t
           x1=x
	  endif
          cont=.true.
         endif
         if (iter.gt.Maxiter) then
          write(6,*)'chem. pot. too many iterations'
          cont=.false.
          stop
         endif
        endif
        call MPI_BCAST(cont,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (cont) goto 10
        xmu=t
        if (myid.eq.0) write(6,*)iter,' iterations' 
        return
	end 
	
