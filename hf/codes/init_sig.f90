program init_sig
! ifort -o init_sig init_sig.f90

implicit none

integer*4  i,j,norb,im,iw
parameter (norb=6)
integer*4  beta,ecut
real*8     pi

open(unit=20,file='sigma',form='formatted',status='unknown')

beta=40
ecut=200
write(*,*)'beta and ecut? (default: 40 and 200)'
read(*,*)beta,ecut
pi=dacos(-1.d0)
im=int(ecut*beta/(2.d0*pi))
write(*,*)'beta, ecut, im:',beta,ecut,im

do i=1,norb
   do j=1,norb
      do iw=0,im
         write(20,40)(2.d0*dble(iw)+1.d0)*pi/beta
      enddo
   enddo
enddo

40 format(e20.10,2(4x,'0.0'))
stop
end program init_sig
