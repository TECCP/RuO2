	subroutine w2t1(Gw,jatom,Nint,xm0,xm1,xm2,iout)
        include 'dmft.dat'
	complex*16 cdummy,Gw,Gt,ga
        complex*16 xm0,xm1,xm2
	dimension Gw(niter,niter,0:momega)
 	dimension Gt(niter,niter,mtau)
	dimension xm0(Niter,Niter),xm1(Niter,Niter)
	dimension xm2(Niter,Niter)

         beta_pi=Pi/Beta
         do 1 i=1,LL
          tau=(i-1)*beta/real(LL)
c.........tail functions
          tail1=0.d0
          tail2=0.d0
          tb=tau/beta
          do n=0,Im
           om=(2*n+1)*beta_pi
           tail1=tail1+cos(om*tau)/(om*om)
           tail2=tail2+sin(om*tau)/(om*om*om)
          enddo
           tail1=beta*(One-Two*tb)/4.d0-Two*tail1/beta
           tail2=beta**2*tb*(One-tb)/4.d0-Two*tail2/beta
c..........numerical summation to the cut-off          
          do 3 k=1,nint
           do 3 l=1,nint
            sum=czero
            do 2 j=0,Im
             om=(2*j+1)*beta_pi
             cdummy=Cim*mod(om*tau,2*Pi)
             dummy=(Gw(k,l,j)+xm0(k,l)*cim/om)*exp(-cdummy)
     &	          +dconjg(Gw(l,k,j)+xm0(k,l)*cim/om)*exp(cdummy)
            sum=sum+dummy
2         continue
c........add assymptotic sums
	 Gt(k,l,i)=sum/Beta-xm1(k,l)*tail1+xm2(k,l)*tail2
     &            -xm0(k,l)*One/Two
3      continue
1     continue	
      if (iout.eq.18.and.diag_only) then
c.....ensure that delta.tau > 0
       do i=1,nint
        do l=1,ll
         if (dreal(gt(i,i,l)).gt.0.d0) gt(i,i,l)=-1.d-10
        enddo
       enddo
c      do k=1,LL
c       if (spinfac.eq.2) then
c        write(18,111)(k-1)*beta/real(LL),
c    &   (-dreal(gt(i,i,k)),i=1,nint),(-dreal(gt(i,i,k)),i=1,nint)
c       else
c        write(18,111)(k-1)*beta/real(LL),
c    &   (-dreal(gt(i,i,k)),i=1,nint)     
c       endif
c      enddo
c      if (spinfac.eq.2) then
c       write(18,111)beta,(dreal(xm0(i,i)+gt(i,i,1)),i=1,nint),
c    &  (dreal(xm0(i,i)+gt(i,i,1)),i=1,nint)
c      else
c       write(18,111)beta,(dreal(xm0(i,i)+gt(i,i,1)),i=1,nint)
c      endif
c     else
c      do i=1,nint
c       do j=1,nint
c        do k=1,LL
c.....sign change for QMC definition of G(t)
c         write(iout,*)(k-1)*beta/real(LL),-dreal(gt(i,j,k)),
c    &    -imag(gt(i,j,k))
c        enddo
c       enddo
c      enddo
      endif
 111  format(f8.4,1x,14f18.10)
 
      end
