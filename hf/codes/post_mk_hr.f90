program post_mk_hr
!====================================================================================================
! ifort -o mk_kanamori_mf mk_kanamori_mf.f90
!
! ifort -O -o mk_kanamori_mf mk_kanamori_mf.f90 -L/opt/intel/mkl/10.2.5.035/lib/em64t -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
!
! date; mpirun -np 2 /home/asurada/program/dmft/SRC_dmft1.1/dmft; date
!====================================================================================================

implicit none

integer*4  i,j,k,l,nr,ii,jj,kpt,nh,bh,natom,ia
integer*4  i1,i2,kzero,count,count2
parameter (nr=2057,kpt=1600,nh=12,natom=2,bh=nh/natom)
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

! for occupation matrix
complex*16 occ_full(nh,nh)
complex*16 occ(natom,bh,bh)

! for Hartree-Fock
real*8     u1,u2,jh,i6(bh,bh)
complex*16 hmf_full(nh,nh)
complex*16 hmf(natom,bh,bh)
complex*16 mat1(natom,bh,bh),shift1(natom)
complex*16 mat2(natom,bh,bh),shift2(natom)
complex*16 mat3(natom,bh,bh),shift3(natom)
complex*16 mat4(natom,bh,bh),shift4(natom)
complex*16 mat5(natom,bh,bh),shift5(natom)
complex*16 shift(natom),dum

open( 99,file='wannier90_hr.dat',   status='old',    form='formatted')
!open(100,file='klist.in',          status='old',    form='formatted')
open( 81,file='n.dmft',             status='old',    form='formatted')
open( 82,file='u.in',               status='old',    form='formatted')
open( 90,file='wannier90_hr_mf.dat',status='unknown',form='formatted')




!-----[begin] read u and j
read(82,*)u1
read(82,*)jh
u2=u1-2.d0*jh
write(*,*)"U =",u1
write(*,*)"U'=",u2,"= U-2J"
write(*,*)"J =",jh
!-----[end]



!-----[begin] read occupation matrix
occ_full(:,:)=0.d0
occ(:,:,:)=0.d0
do i=1,nh
   do j=1,nh
      read(81,*)a,b
      occ_full(i,j)=dcmplx(a,b)
   enddo
enddo
write(*,*)'re[occ]:'
do j=1,nh
   write(*,9990)(dble(occ_full(i,j)),i=1,nh)
enddo
write(*,*)'im[occ]:'
do j=1,nh
   write(*,9990)(dimag(occ_full(i,j)),i=1,nh)
enddo
do ia=1,natom
   do i=1,bh
      do j=1,bh
         occ(ia,i,j)=occ_full(i+(ia-1)*bh,j+(ia-1)*bh)
      enddo
   enddo
enddo
!do ia=1,natom
!   write(*,*)'re[occ] for atom',ia,':'
!   do j=1,bh
!      write(*,9990)(dble(occ(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'im[occ] for atom',ia,':'
!   do j=1,bh
!      write(*,9990)(dimag(occ(ia,i,j)),i=1,bh)
!   enddo
!enddo
!-----[end]



!-----[begin] make Hartree-Fock matrix
hmf_full(:,:)=0.d0
hmf(:,:,:)=0.d0
mat1(:,:,:)=0.d0
mat2(:,:,:)=0.d0
mat3(:,:,:)=0.d0
mat4(:,:,:)=0.d0
mat5(:,:,:)=0.d0
shift1(:)=0.d0
shift2(:)=0.d0
shift3(:)=0.d0
shift4(:)=0.d0
shift5(:)=0.d0
i6(:,:)=0.d0
do i=1,bh
   i6(i,i)=1.d0
enddo

! part (1): U
do ia=1,natom
   mat1(ia,1,1)=occ(ia,4,4)
   mat1(ia,2,2)=occ(ia,5,5)
   mat1(ia,3,3)=occ(ia,6,6)
   mat1(ia,4,4)=occ(ia,1,1)
   mat1(ia,5,5)=occ(ia,2,2)
   mat1(ia,6,6)=occ(ia,3,3)

   mat1(ia,1,4)=-occ(ia,4,1)
   mat1(ia,2,5)=-occ(ia,5,2)
   mat1(ia,3,6)=-occ(ia,6,3)
   mat1(ia,4,1)=-occ(ia,1,4)
   mat1(ia,5,2)=-occ(ia,2,5)
   mat1(ia,6,3)=-occ(ia,3,6)

   shift1(ia)= occ(ia,1,4)*occ(ia,4,1)+occ(ia,2,5)*occ(ia,5,2)+occ(ia,3,6)*occ(ia,6,3) &
   &          -occ(ia,1,1)*occ(ia,4,4)-occ(ia,2,2)*occ(ia,5,5)-occ(ia,3,3)*occ(ia,6,6)
enddo

! part (2): U'
do ia=1,natom
   mat2(ia,1,1)=occ(ia,5,5)+occ(ia,6,6)
   mat2(ia,2,2)=occ(ia,4,4)+occ(ia,6,6)
   mat2(ia,3,3)=occ(ia,4,4)+occ(ia,5,5)
   mat2(ia,4,4)=occ(ia,2,2)+occ(ia,3,3)
   mat2(ia,5,5)=occ(ia,1,1)+occ(ia,3,3)
   mat2(ia,6,6)=occ(ia,1,1)+occ(ia,2,2)

   mat2(ia,1,5)=-occ(ia,5,1)
   mat2(ia,1,6)=-occ(ia,6,1)
   mat2(ia,2,4)=-occ(ia,4,2)
   mat2(ia,2,6)=-occ(ia,6,2)
   mat2(ia,3,4)=-occ(ia,4,3)
   mat2(ia,3,5)=-occ(ia,5,3)

   mat2(ia,4,2)=-occ(ia,2,4)
   mat2(ia,4,3)=-occ(ia,3,4)
   mat2(ia,5,1)=-occ(ia,1,5)
   mat2(ia,5,3)=-occ(ia,3,5)
   mat2(ia,6,1)=-occ(ia,1,6)
   mat2(ia,6,2)=-occ(ia,2,6)

   shift2(ia)= occ(ia,1,5)*occ(ia,5,1)+occ(ia,1,6)*occ(ia,6,1) &
   &          +occ(ia,2,4)*occ(ia,4,2)+occ(ia,2,6)*occ(ia,6,2) &
   &          +occ(ia,3,4)*occ(ia,4,3)+occ(ia,3,5)*occ(ia,5,3) &
   &          -occ(ia,1,1)*occ(ia,5,5)-occ(ia,1,1)*occ(ia,6,6) &
   &          -occ(ia,2,2)*occ(ia,4,4)-occ(ia,2,2)*occ(ia,6,6) &
   &          -occ(ia,3,3)*occ(ia,4,4)-occ(ia,3,3)*occ(ia,5,5)
enddo

! part (3): U'-J
do ia=1,natom
   mat3(ia,1,1)=occ(ia,2,2)+occ(ia,3,3)
   mat3(ia,2,2)=occ(ia,1,1)+occ(ia,3,3)
   mat3(ia,3,3)=occ(ia,1,1)+occ(ia,2,2)
   mat3(ia,4,4)=occ(ia,5,5)+occ(ia,6,6)
   mat3(ia,5,5)=occ(ia,4,4)+occ(ia,6,6)
   mat3(ia,6,6)=occ(ia,4,4)+occ(ia,5,5)

   mat3(ia,1,2)=-occ(ia,2,1)
   mat3(ia,1,3)=-occ(ia,3,1)
   mat3(ia,2,1)=-occ(ia,1,2)
   mat3(ia,2,3)=-occ(ia,3,2)
   mat3(ia,3,1)=-occ(ia,1,3)
   mat3(ia,3,2)=-occ(ia,2,3)

   mat3(ia,4,5)=-occ(ia,5,4)
   mat3(ia,4,6)=-occ(ia,6,4)
   mat3(ia,5,4)=-occ(ia,4,5)
   mat3(ia,5,6)=-occ(ia,6,5)
   mat3(ia,6,4)=-occ(ia,4,6)
   mat3(ia,6,5)=-occ(ia,5,6)

   shift3(ia)= occ(ia,1,2)*occ(ia,2,1)+occ(ia,4,5)*occ(ia,5,4) &
   &          +occ(ia,1,3)*occ(ia,3,1)+occ(ia,4,6)*occ(ia,6,4) &
   &          +occ(ia,2,3)*occ(ia,3,2)+occ(ia,5,6)*occ(ia,6,5) &
   &          -occ(ia,1,1)*occ(ia,2,2)-occ(ia,4,4)*occ(ia,5,5) &
   &          -occ(ia,1,1)*occ(ia,3,3)-occ(ia,4,4)*occ(ia,6,6) &
   &          -occ(ia,2,2)*occ(ia,3,3)-occ(ia,5,5)*occ(ia,6,6)
enddo

! part (4): J
do ia=1,natom
   mat4(ia,1,2)=occ(ia,5,4)
   mat4(ia,1,3)=occ(ia,6,4)
   mat4(ia,2,1)=occ(ia,4,5)
   mat4(ia,2,3)=occ(ia,6,5)
   mat4(ia,3,1)=occ(ia,4,6)
   mat4(ia,3,2)=occ(ia,5,6)

   mat4(ia,4,5)=occ(ia,2,1)
   mat4(ia,4,6)=occ(ia,3,1)
   mat4(ia,5,4)=occ(ia,1,2)
   mat4(ia,5,6)=occ(ia,3,2)
   mat4(ia,6,4)=occ(ia,1,3)
   mat4(ia,6,5)=occ(ia,2,3)

   mat4(ia,1,4)=-occ(ia,5,2)-occ(ia,6,3)
   mat4(ia,2,5)=-occ(ia,4,1)-occ(ia,6,3)
   mat4(ia,3,6)=-occ(ia,4,1)-occ(ia,5,2)
   mat4(ia,4,1)=-occ(ia,2,5)-occ(ia,3,6)
   mat4(ia,5,2)=-occ(ia,1,4)-occ(ia,3,6)
   mat4(ia,6,3)=-occ(ia,1,4)-occ(ia,2,5)

   shift4(ia)=-occ(ia,1,2)*occ(ia,5,4)-occ(ia,1,3)*occ(ia,6,4) &
   &          -occ(ia,2,1)*occ(ia,4,5)-occ(ia,2,3)*occ(ia,6,5) &
   &          -occ(ia,3,1)*occ(ia,4,6)-occ(ia,3,2)*occ(ia,5,6) &
   &          +occ(ia,1,4)*occ(ia,5,2)+occ(ia,1,4)*occ(ia,6,3) &
   &          +occ(ia,2,5)*occ(ia,4,1)+occ(ia,2,5)*occ(ia,6,3) &
   &          +occ(ia,3,6)*occ(ia,4,1)+occ(ia,3,6)*occ(ia,5,2)
enddo

! part (5): J
do ia=1,natom
   mat5(ia,1,2)=occ(ia,4,5)
   mat5(ia,1,3)=occ(ia,4,6)
   mat5(ia,2,1)=occ(ia,5,4)
   mat5(ia,2,3)=occ(ia,5,6)
   mat5(ia,3,1)=occ(ia,6,4)
   mat5(ia,3,2)=occ(ia,6,5)

   mat5(ia,4,5)=occ(ia,1,2)
   mat5(ia,4,6)=occ(ia,1,3)
   mat5(ia,5,4)=occ(ia,2,1)
   mat5(ia,5,6)=occ(ia,2,3)
   mat5(ia,6,4)=occ(ia,3,1)
   mat5(ia,6,5)=occ(ia,3,2)

   mat5(ia,1,5)=-occ(ia,4,2)
   mat5(ia,1,6)=-occ(ia,4,3)
   mat5(ia,2,4)=-occ(ia,5,1)
   mat5(ia,2,6)=-occ(ia,5,3)
   mat5(ia,3,4)=-occ(ia,6,1)
   mat5(ia,3,5)=-occ(ia,6,2)

   mat5(ia,4,2)=-occ(ia,1,5)
   mat5(ia,4,3)=-occ(ia,1,6)
   mat5(ia,5,1)=-occ(ia,2,4)
   mat5(ia,5,3)=-occ(ia,2,6)
   mat5(ia,6,1)=-occ(ia,3,4)
   mat5(ia,6,2)=-occ(ia,3,5)

   shift5(ia)=-occ(ia,1,2)*occ(ia,4,5)-occ(ia,1,3)*occ(ia,4,6) &
   &          -occ(ia,2,1)*occ(ia,5,4)-occ(ia,2,3)*occ(ia,5,6) &
   &          -occ(ia,3,1)*occ(ia,6,4)-occ(ia,3,2)*occ(ia,6,5) &
   &          +occ(ia,1,5)*occ(ia,4,2)+occ(ia,1,6)*occ(ia,4,3) &
   &          +occ(ia,2,4)*occ(ia,5,1)+occ(ia,2,6)*occ(ia,5,3) &
   &          +occ(ia,3,4)*occ(ia,6,1)+occ(ia,3,5)*occ(ia,6,2)
enddo

!do ia=1,natom
!   write(*,*)'re[hmf] for atom',ia,':'
!   write(*,*)'mat1:'
!   do j=1,bh
!      write(*,9990)(dble(mat1(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat2:'
!   do j=1,bh
!      write(*,9990)(dble(mat2(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat3:'
!   do j=1,bh
!      write(*,9990)(dble(mat3(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat4:'
!   do j=1,bh
!      write(*,9990)(dble(mat4(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat5:'
!   do j=1,bh
!      write(*,9990)(dble(mat5(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat6:'
!   do j=1,bh
!      write(*,9990)(dble(mat6(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'im[hmf] for atom',ia,':'
!   write(*,*)'mat1:'
!   do j=1,bh
!      write(*,9990)(dimag(mat1(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat2:'
!   do j=1,bh
!      write(*,9990)(dimag(mat2(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat3:'
!   do j=1,bh
!      write(*,9990)(dimag(mat3(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat4:'
!   do j=1,bh
!      write(*,9990)(dimag(mat4(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat5:'
!   do j=1,bh
!      write(*,9990)(dimag(mat5(ia,i,j)),i=1,bh)
!   enddo
!   write(*,*)'mat6:'
!   do j=1,bh
!      write(*,9990)(dimag(mat6(ia,i,j)),i=1,bh)
!   enddo
!enddo

shift(:)=     u1*shift1(:) &
&       +     u2*shift2(:) &
&       +(u2-jh)*shift3(:) &
&       +     jh*shift4(:) &
&       +     jh*shift5(:)
write(*,*)'shift_atom1:',shift(1)
write(*,*)'shift_atom2:',shift(2)

! sym_shift: symmetrize the chemical potentials of Ru1 and Ru2
!             if not, one of them goes upward, and another one goes downward!
!             why this happens, OTL ...
dum=0.d0
do ia=1,natom
   dum=dum+shift(ia)
enddo
dum=dum/dble(natom)
do ia=1,natom
   shift(ia)=dum
enddo
!write(*,*)'shift_atom1:',shift(1),'(symmetrized)'
!write(*,*)'shift_atom2:',shift(2),'(symmetrized)'
!
!write(*,*)'shift1_atom1:',shift1(1)
!write(*,*)'shift1_atom2:',shift1(2)
!
!write(*,*)'shift2_atom1:',shift2(1)
!write(*,*)'shift2_atom2:',shift2(2)
!
!write(*,*)'shift3_atom1:',shift3(1)
!write(*,*)'shift3_atom2:',shift3(2)
!
!write(*,*)'shift4_atom1:',shift4(1)
!write(*,*)'shift4_atom2:',shift4(2)
!
!write(*,*)'shift5_atom1:',shift5(1)
!write(*,*)'shift5_atom2:',shift5(2)

hmf(:,:,:)=     u1*mat1(:,:,:) &
&         +     u2*mat2(:,:,:) &
&         +(u2-jh)*mat3(:,:,:) &
&         +     jh*mat4(:,:,:) &
&         +     jh*mat5(:,:,:)

! IMPORTANT NOTE ON 11JAN2019:
! CURRENTLY, I TURNED OFF THE SHIFT.
! WHEN I UNDERSTAND THE EFFECTS OF THE SHIFT,
! I WILL RECOVER HERE.
!
!do ia=1,natom
!   do i=1,bh
!      hmf(ia,i,i)=hmf(ia,i,i)+shift(ia)
!   enddo
!enddo

do ia=1,natom
   do i=1,bh
      do j=1,bh
         hmf_full(i+(ia-1)*bh,j+(ia-1)*bh)=hmf(ia,i,j)
      enddo
   enddo
enddo

write(*,*)'re[hmf]:'
do j=1,nh
   write(*,9990)(dble(hmf_full(i,j)),i=1,nh)
enddo
write(*,*)'im[hmf]:'
do j=1,nh
   write(*,9990)(dimag(hmf_full(i,j)),i=1,nh)
enddo
!-----[end]



!-----[begin] read k-points for hamilt
!klist(:,:)=0.d0
!read(100,*)
!do i=1,kpt
!   read(100,*)klist(1,i),klist(2,i),klist(3,i)
!enddo
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
write(*,*)'re[H(R0)]:'
do j=1,nh
   write(*,9990)(dble(hamr(i,j,kzero)),i=1,nh)
enddo
write(*,*)'im[H(R0)]:'
do j=1,nh
   write(*,9990)(dimag(hamr(i,j,kzero)),i=1,nh)
enddo
!
hamr(:,:,kzero)=hamr(:,:,kzero)+hmf_full(:,:)
write(*,*)'new.mf.re[H(R0)]:'
do j=1,nh
   write(*,9990)(dble(hamr(i,j,kzero)),i=1,nh)
enddo
write(*,*)'new.mf.im[H(R0)]:'
do j=1,nh
   write(*,9990)(dimag(hamr(i,j,kzero)),i=1,nh)
enddo
!-----[end]



!-----[begin] write the new H(R)
write(90,*)
write(90,*)nh
write(90,*)nr
do k=1,nr
   do j=1,nh
      do i=1,nh
         write(90,9993)avec(1,k),avec(2,k),avec(3,k),i,j,dble(hamr(i,j,k)),dimag(hamr(i,j,k))
      enddo
   enddo
enddo
!-----[end]



!-----[begin] make H(k)
!hk(:,:,:)=0.d0
!do l=1,kpt
!   do j=1,nr
!      phase=0.d0
!      do i=1,3
!         phase=phase+pi2*klist(i,l)*rvec(i,j)
!      enddo
!      do ii=1,nh
!         do jj=1,nh
!            hk(ii,jj,l)=hk(ii,jj,l) &
!&                      +hamr(ii,jj,j)*dcmplx(dcos(phase),dsin(phase))
!         enddo
!      enddo
!   enddo
!enddo
!do l=1,kpt
!   hk(:,:,l)=hk(:,:,l)+hmf_full(:,:)
!enddo
!write(90,*)kpt,nh
!do l=1,kpt
!   write(90,9991)1,klist(1,l),klist(2,l),klist(3,l)
!   do i=1,nh
!      do j=1,nh
!         write(90,*)dble(hk(i,j,l)),dimag(hk(i,j,l))
!      enddo
!   enddo
!enddo
!-----[end]



!-----[begin] make DOS
!do mainit=-200,1200
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
9993 format(5i5,2f12.6)

stop
end program post_mk_hr
