	subroutine mu_calc_newton(xmu,x,iret)
        include 'dmft.dat'
        include 'mpif.h'
        double precision my_x,my_dx
        logical cont

        if (myid.eq.0) then
         write(6,*)'****************************************'
         write(6,*)'CALULATION OF CHEM. POTENTTIAL'
         write(6,*)'****************************************'
        endif
c--------initial occupation
        kmin=nk1(myid+1)
        kmax=nk2(myid+1)
        
        call xnd(xmu,my_x,my_dx,kmin,kmax)
        call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(my_dx,dx,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,ierr)
c
	iter=1
        om=(2*Im+1)*pi/beta
        Xtol=1.d-1*(emax/om)**3
        Xtol=max(Xtol,1.d-6)
        if (myid.eq.0) then
         write(6,*)'n_tot toletance =',xtol
         write(6,*)'convergence:'
        endif
        
        if (myid.eq.0) then
         t0=xmu
         diff=x_el-x
         dt=diff/dx
         t=t0+dt
         write(6,*)iter,t0,x
        endif 
        

 10     call MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call xnd(t,my_x,my_dx,kmin,kmax)
        call MPI_REDUCE(my_x,x,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,err)
        call MPI_REDUCE(my_dx,dx,1,MPI_DOUBLE_PRECISION,
     &  MPI_SUM,0,MPI_COMM_WORLD,err)
        iter=iter+1
        if (myid.eq.0) then
         write(6,*)iter,t,x
        endif
        if (myid.eq.0) then
         if (abs(x_el-x).le.Xtol) then
          iret=0
          cont=.false.
         else
          if (abs(x_el-x).gt.abs(diff)) then
           dt=0.5*dt
	  else
           diff=x_el-x
           t0=t
           dt=diff/dx
           t=t0+dt
          endif
          cont=.true.
         endif
         if (iter.gt.10) then
          write(6,*)'chem. pot. too many iterations'
          iret=1
          cont=.false.
         endif
        endif
        call MPI_BCAST(cont,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iret,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (cont) goto 10
        if (iret.eq.0) xmu=t
        return
	end 
	
