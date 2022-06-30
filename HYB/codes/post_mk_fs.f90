program post_mk_fs
! ifort -O -o post_mk_fs post_mk_fs.f90 -L/opt/intel-14.0.3/composer_xe_2013_sp1.3.174/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

implicit none

integer*4  i,j,k,l,nr,ii,jj,kpt,nh,bh
integer*4  i1,i2,kzero,count,count2
parameter (nr=2057,kpt=1540351,nh=12,bh=nh/2)
integer*4  avec(3,nr)
real*8     rvec(3,nr),klist(3),phase,a,b,pi2,hr_cutoff
parameter (pi2=2.d0*dacos(-1.d0))
parameter (hr_cutoff=0.d0)
complex*16 hamr(nh,nh,nr)
complex*16 hk(nh,nh),hk2(nh,nh)

! for zheev
integer*4  lwork,info
parameter (lwork=300)
real*8     ene(bh),rwork(300)
complex*16 work(lwork)

! for sum
integer*4  new(nh)
!real*8    esum(2),diff_cut
real*8     kcoa,ef
!parameter(diff_cut=-0.2d0)
parameter (kcoa=2.022846d0/1.398758d0)
complex*16 bandup(bh,bh),banddn(bh,bh)

open( 98,file='wannier90_hr_mf.dat',status='old',    form='formatted')
open(100,file='klist.in',           status='old',    form='formatted')
open(101,file='ef.in',              status='old',    form='formatted')

open(111,file='fs_3d_up_band1.out', status='unknown',form='formatted')
open(112,file='fs_3d_up_band2.out', status='unknown',form='formatted')
open(113,file='fs_3d_up_band3.out', status='unknown',form='formatted')
open(114,file='fs_3d_up_band4.out', status='unknown',form='formatted')
open(115,file='fs_3d_up_band5.out', status='unknown',form='formatted')
open(116,file='fs_3d_up_band6.out', status='unknown',form='formatted')

open(121,file='fs_3d_dn_band1.out', status='unknown',form='formatted')
open(122,file='fs_3d_dn_band2.out', status='unknown',form='formatted')
open(123,file='fs_3d_dn_band3.out', status='unknown',form='formatted')
open(124,file='fs_3d_dn_band4.out', status='unknown',form='formatted')
open(125,file='fs_3d_dn_band5.out', status='unknown',form='formatted')
open(126,file='fs_3d_dn_band6.out', status='unknown',form='formatted')



!-----[begin] read H(R)
read(98,*)
read(98,*)
read(98,*)
count=0
count2=0
avec(:,:)=0
rvec(:,:)=0.d0
hamr(:,:,:)=0.d0
do k=1,nr
   do i=1,nh
      do j=1,nh
         read(98,*)avec(1,k),avec(2,k),avec(3,k),i1,i2,a,b
         if (dsqrt(a**2.d0+b**2.d0)>hr_cutoff) then
            count=count+1
            hamr(i1,i2,k)=dcmplx(a,b)
         else
            count2=count2+1
            hamr(i1,i2,k)=0.d0
         endif
         if (avec(1,k)**2+avec(2,k)**2+avec(3,k)**2==0) then
            kzero=k
         endif
      enddo
   enddo
enddo
rvec(:,:)=dble(avec(:,:))
write(*,*)count,count2,count+count2
write(*,*)rvec(1,kzero),rvec(2,kzero),rvec(3,kzero)
write(*,*)'re[H(R0)]:'
do j=1,nh
   write(*,9990)(dble(hamr(i,j,kzero)),i=1,nh)
enddo
write(*,*)'im[H(R0)]:'
do j=1,nh
   write(*,9990)(dimag(hamr(i,j,kzero)),i=1,nh)
enddo
!-----[end]



!-----[begin] make H(k)
do i=1,nh
   new(i)=i
enddo
new(4)=7
new(5)=8
new(6)=9
new(7)=4
new(8)=5
new(9)=6
ef=0.d0
read(101,*)ef
read(100,*)

do l=1,kpt!l

   hk(:,:)=0.d0
   read(100,*)klist(1),klist(2),klist(3)

   do j=1,nr
      phase=0.d0
      do i=1,3
         phase=phase+pi2*klist(i)*rvec(i,j)
      enddo
      do ii=1,nh
         do jj=1,nh
            hk(ii,jj)=hk(ii,jj) &
   &                 +hamr(ii,jj,j)*dcmplx(dcos(phase),dsin(phase))
         enddo
      enddo
   enddo

   ! rearrange up and down
   hk2(:,:)=hk(:,:)
   do i=1,nh
      do j=1,nh
         hk(i,j)=hk2(new(i),new(j))
      enddo
   enddo

   do i=1,bh
      do j=1,bh
         bandup(i,j)=hk(i,j)
      enddo
   enddo
   call zheev('V','U',bh,bandup,bh,ene,work,lwork,rwork,info)
   do k=1,bh
      write(110+k,9994)klist(1),klist(2),klist(3)*kcoa,ene(k)-ef
   enddo
   do i=1,bh
      do j=1,bh
         banddn(i,j)=hk(i+bh,j+bh)
      enddo
   enddo
   call zheev('V','U',bh,banddn,bh,ene,work,lwork,rwork,info)
   do k=1,bh
      write(120+k,9994)klist(1),klist(2),klist(3)*kcoa,ene(k)-ef
   enddo

enddo!l
!-----[end]



9990 format(12f12.6)
!9992 format(25e14.6)
9994 format(4f12.8)

stop
end program post_mk_fs
