program mx_occ
! ifort -o mx_occ mx_occ.f90

implicit none

integer*4  i,j,k,nh,bh,natom,ia
parameter (nh=12,natom=2,bh=nh/natom)
real*8     a,b

! for occupation matrix
complex*16 occ_full_old(nh,nh)
complex*16 occ_old(natom,bh,bh)
complex*16 occ_full(nh,nh)
complex*16 occ(natom,bh,bh)
complex*16 occ_full_mix(nh,nh)
complex*16 occ_mix(natom,bh,bh)

open( 80,file='n.dmft_old',status='old',    form='formatted')
open( 81,file='n.dmft',    status='old',    form='formatted')
open( 82,file='n.dmft_mix',status='unknown',form='formatted')



!-----[begin] read occupation matrix
occ_full_old(:,:)=0.d0
occ_old(:,:,:)=0.d0
occ_full(:,:)=0.d0
occ(:,:,:)=0.d0
occ_full_mix(:,:)=0.d0
occ_mix(:,:,:)=0.d0
do i=1,nh
   do j=1,nh
      read(80,*)a,b
      occ_full_old(i,j)=dcmplx(a,b)
      read(81,*)a,b
      occ_full(i,j)=dcmplx(a,b)
   enddo
enddo

do i=1,nh
   do j=1,nh
      occ_full_mix(i,j)=(occ_full_old(i,j)+occ_full(i,j))*0.5d0
   enddo
enddo

do ia=1,natom
   do i=1,bh
      do j=1,bh
         occ_old(ia,i,j)=occ_full_old(i+(ia-1)*bh,j+(ia-1)*bh)
         occ(ia,i,j)=occ_full(i+(ia-1)*bh,j+(ia-1)*bh)
         occ_mix(ia,i,j)=occ_full_mix(i+(ia-1)*bh,j+(ia-1)*bh)
      enddo
   enddo
enddo

do ia=1,natom
   write(*,*)'re[occ_new-occ_old] for atom',ia,':'
   do j=1,bh
      write(*,9990)(dble(occ(ia,i,j)-occ_old(ia,i,j)),i=1,bh)
   enddo
enddo
do ia=1,natom
   write(*,*)'re[occ_mix-occ_old] for atom',ia,':'
   do j=1,bh
      write(*,9990)(dble(occ_mix(ia,i,j)-occ_old(ia,i,j)),i=1,bh)
   enddo
enddo

do i=1,nh
   do j=1,nh
      write(82,*)dble(occ_full_mix(i,j)),dimag(occ_full_mix(i,j))
   enddo
enddo
!-----[end]



9990 format(12f8.4)

stop
end program mx_occ
