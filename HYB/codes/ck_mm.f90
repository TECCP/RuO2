program ck_mm
! ifort -o ck_mm ck_mm.f90

implicit none

integer*4  i,j,nh
parameter (nh=12)
real*8     a,b
real*8     ru(nh)
real*8     mm1,mm2,mmavg

! for occupation matrix
complex*16 occ_full(nh,nh)

open( 81,file='n.dmft',    status='old',    form='formatted')



occ_full(:,:)=0.d0
ru(:)=0.d0
do i=1,nh
   do j=1,nh
      read(81,*)a,b
      occ_full(i,j)=dcmplx(a,b)
   enddo
enddo

do i=1,nh
   ru(i)=occ_full(i,i)
enddo

mm1=ru(1)+ru(2)+ru(3)-ru(4)-ru(5)-ru(6)
mm2=ru(7)+ru(8)+ru(9)-ru(10)-ru(11)-ru(12)
mmavg=(mm1-mm2)*0.5d0

write(*,*)'MM_Ru1=',mm1
write(*,*)'MM_Ru2=',mm2
write(*,*)'MM_avg=',mmavg



stop
end program ck_mm
