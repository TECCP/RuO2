program mk_klist

implicit none

integer*4 i,j,k,nkx,nky,nkz
parameter (nkx=100,nky=100,nkz=150)
real*8    kx,ky,kz

open(101,file='klist.in',status='unknown',form='formatted')

write(101,*)'#',nkx+1,'x',nky+1,'x',nkz+1,'=',(nkx+1)*(nky+1)*(nkz+1)

do i=0,nkx
   do j=0,nky
      do k=0,nkz
         kx=dble(i)/dble(nkx)-0.5d0
         ky=dble(j)/dble(nky)-0.5d0
         kz=dble(k)/dble(nkz)-0.5d0
         write(101,9991)kx,ky,kz
      enddo
   enddo
enddo


9991 format(3f20.8)
stop
end program mk_klist
