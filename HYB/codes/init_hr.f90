program init_hr
! ifort -o init_hr init_hr.f90

implicit none

integer*4  i,j,k,nr,nh,nh2,is
integer*4  i1,i2,kzero,count,count2
parameter (nr=2057)
parameter (nh=6,nh2=nh*2)
integer*8  new(nh2)
integer*4  avec(3,nr)
real*8     rvec(3,nr),a,b,hr_cutoff
parameter (hr_cutoff=0.d0)
complex*16 hamr0(2,nh,nh,nr)
complex*16 hamr(nh2,nh2,nr),hamr2(nh2,nh2,nr)

open(91,file='wannier90_hr.datup',status='old',    form='formatted')
open(92,file='wannier90_hr.datdn',status='old',    form='formatted')
open(99,file='wannier90_hr.dat',  status='unknown',form='formatted')



!-----[begin] read H(R)
do is=1,2

   read(90+is,*)
   read(90+is,*)
   read(90+is,*)
   count=0
   count2=0
   avec(:,:)=0
   rvec(:,:)=0.d0
   hamr0(is,:,:,:)=0.d0
   do k=1,nr
      do i=1,nh
         do j=1,nh
            read(90+is,*)avec(1,k),avec(2,k),avec(3,k),i1,i2,a,b
            if (dsqrt(a**2.d0+b**2.d0)>hr_cutoff) then
               count=count+1
               hamr0(is,i1,i2,k)=dcmplx(a,b)
            else
               count2=count2+1
               hamr0(is,i1,i2,k)=0.d0
            endif
            if (avec(1,k)**2+avec(2,k)**2+avec(3,k)**2==0) then
               kzero=k
            endif
         enddo
      enddo
   enddo
   rvec(:,:)=dble(avec(:,:))
   write(*,*)'nonzero elements:',count
   write(*,*)'   zero elements:',count2
   write(*,*)'  total elements:',count+count2
   write(*,*)'nr,kzero:',nr,kzero
   write(*,*)rvec(1,kzero),rvec(2,kzero),rvec(3,kzero)

   write(*,*)'mat.re[H(R0)]:'
   do j=1,nh
      write(*,9992)(dble(hamr0(is,i,j,kzero)),i=1,nh)
   enddo
   write(*,*)'mat.im[H(R0)]:'
   do j=1,nh
      write(*,9992)(dimag(hamr0(is,i,j,kzero)),i=1,nh)
   enddo

enddo
!-----[end]



!-----[begin] merge the up and down
hamr(:,:,:)=0.d0
do k=1,nr
   do i=1,nh
      do j=1,nh
         hamr(i,j,k)=hamr0(1,i,j,k)
         hamr(i+nh,j+nh,k)=hamr0(2,i,j,k)
      enddo
   enddo
enddo
!-----[end]



!-----[begin] rearrange H(R)
do i=1,nh2
   new(i)=i
enddo
hamr2(:,:,:)=hamr(:,:,:)
new(4)=7
new(5)=8
new(6)=9
new(7)=4
new(8)=5
new(9)=6
do i=1,nh2
   do j=1,nh2
      hamr(i,j,:)=hamr2(new(i),new(j),:)
   enddo
enddo

write(*,8002)(new(i),i=1,nh2)
write(*,*)'form.order.mat.re.diag[H(R0)]:'
do j=1,nh2
   write(*,9992)(dble(hamr(i,j,kzero)),i=1,nh2)
enddo
write(*,*)'form.order.mat.im.diag[H(R0)]:'
do j=1,nh2
   write(*,9992)(dimag(hamr(i,j,kzero)),i=1,nh2)
enddo
!-----[end]



!-----[begin] write the diagonalized H(R)
write(99,*)
write(99,*)nh2
write(99,*)nr
do k=1,nr
   do j=1,nh2
      do i=1,nh2
         write(99,9993)avec(1,k),avec(2,k),avec(3,k),i,j,dble(hamr(i,j,k)),dimag(hamr(i,j,k))
      enddo
   enddo
enddo
!-----[end]



8002 format(12i5)
9992 format(12f12.6)
9993 format(5i5,2f12.6)

stop
end program init_hr
