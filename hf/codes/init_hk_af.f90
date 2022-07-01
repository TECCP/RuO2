program init_hk_af
!
! ifort -o init_hk_af init_hk_af.f90
!
! ifort -O -o init_hk_af init_hk_af.f90 -L/opt/intel/mkl/10.2.5.035/lib/em64t -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

implicit none

integer*4  i,j,k,l,nr,ii,jj,kpt,nh
integer*4  i1,i2,kzero,count,count2
parameter (nr=2057,kpt=1600,nh=12)
integer*4  avec(3,nr)
real*8     rvec(3,nr),klist(3,kpt),phase,a,b,pi2,hr_cutoff
parameter (pi2=2.d0*dacos(-1.d0))
parameter (hr_cutoff=0.d0)
complex*16 hamr(nh,nh,nr)
complex*16 U(nh,nh)
complex*16 hk(nh,nh,kpt)

! for gfwork
!integer*4  mainit,ipiv(nh),info
!real*8     ef
!parameter (ef=0.d0)
!complex*16 z,gfwork(nh),gf(nh,nh)

! for af
real*8     es

open( 99,file='wannier90_hr.dat',status='old',    form='formatted')
open(100,file='klist.in',        status='old',    form='formatted')
open( 90,file='hamilt',          status='unknown',form='formatted')
!open( 92,file='dos.out',         status='unknown',form='formatted')



!-----[begin] read exchange energy
es=0.d0
write(*,*)'es? >>> E(up)-es/2 and E(dn)+es/2'
read(*,*)es
!-----[end]



!-----[begin] read k-points for hamilt
klist(:,:)=0.d0
read(100,*)
do i=1,kpt
   read(100,*)klist(1,i),klist(2,i),klist(3,i)
enddo
!-----[end]



!-----[begin] read H(R)
read(99,*)
read(99,*)
read(99,*)
count=0
count2=0
avec(:,:)=0
rvec(:,:)=0.d0
hamr(:,:,:)=0.d0
do k=1,nr
   do i=1,nh
      do j=1,nh
         read(99,*)avec(1,k),avec(2,k),avec(3,k),i1,i2,a,b
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
!
hamr(1,1,kzero)=hamr(1,1,kzero)-es/2.d0
hamr(2,2,kzero)=hamr(2,2,kzero)-es/2.d0
hamr(3,3,kzero)=hamr(3,3,kzero)-es/2.d0
hamr(4,4,kzero)=hamr(4,4,kzero)+es/2.d0
hamr(5,5,kzero)=hamr(5,5,kzero)+es/2.d0
hamr(6,6,kzero)=hamr(6,6,kzero)+es/2.d0
hamr(7,7,kzero)=hamr(7,7,kzero)+es/2.d0
hamr(8,8,kzero)=hamr(8,8,kzero)+es/2.d0
hamr(9,9,kzero)=hamr(9,9,kzero)+es/2.d0
hamr(10,10,kzero)=hamr(10,10,kzero)-es/2.d0
hamr(11,11,kzero)=hamr(11,11,kzero)-es/2.d0
hamr(12,12,kzero)=hamr(12,12,kzero)-es/2.d0
!
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
hk(:,:,:)=0.d0
do l=1,kpt
   do j=1,nr
      phase=0.d0
      do i=1,3
         phase=phase+pi2*klist(i,l)*rvec(i,j)
      enddo
      do ii=1,nh
         do jj=1,nh
            hk(ii,jj,l)=hk(ii,jj,l) &
&                      +hamr(ii,jj,j)*dcmplx(dcos(phase),dsin(phase))
         enddo
      enddo
   enddo
enddo
write(90,*)kpt,nh
do l=1,kpt
   write(90,9991)1,klist(1,l),klist(2,l),klist(3,l)
   do i=1,nh
      do j=1,nh
         write(90,*)dble(hk(i,j,l)),dimag(hk(i,j,l))
      enddo
   enddo
enddo
!-----[end]



!-----[begin] make DOS
!do mainit=-200,200
!   z=dcmplx(dble(mainit)/dble(100),0.05d0)
!   U(:,:)=0.d0
!   gf(:,:)=0.d0
!   do k=1,kpt
!      do i=1,nh
!         do j=1,nh
!            U(i,j)=-hk(i,j,k)
!         enddo
!      enddo
!      do i=1,nh
!         U(i,i)=z+ef-hk(i,i,k)
!      enddo
!      call zgetrf(nh,nh,U,nh,ipiv,info)
!      if (info/=0) write(*,*)'info1=',info
!      call zgetri(nh,U,nh,ipiv,gfwork,nh,info)
!      if (info/=0) write(*,*)'info2=',info
!      do i=1,nh
!         do j=1,nh
!            gf(i,j)=gf(i,j)+U(i,j)/dble(kpt)
!         enddo
!      enddo
!   enddo
!   write(92,9993)dble(z),(dimag(-gf(i,i)),i=1,nh)
!enddo
!-----[end]



9990 format(12f12.6)
9991 format(i2,3f10.6)
!9993 format(25e13.4)

stop
end program init_hk_af
