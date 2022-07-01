        subroutine asym(licht)
        include 'dmft.dat'

        logical licht
        complex*16 her(niter,niter),aher(niter,niter)
        complex*16 y,omega

c......calculates the constant and 1/iw parts
c......of self-energy from num+1 largest frequencies
        
        do 100 jatom=1,natom_int
        Nint=n_int(jatom)
        index0=n_min(jatom)-1

        beta_pi=Pi/Beta
        do i=1,nint
         do j=1,nint
          her(i,j)=(0.d0,0.d0)
          aher(i,j)=(0.d0,0.d0)
         enddo
        enddo
        ia=10
        do n=Im-ia+1,Im
         om=(2*n+1)*beta_pi
         omega=cim*om
         do i=1,Nint
          do j=1,Nint
           y=(eta(i,j,n,jatom)+
     &        dconjg(eta(j,i,n,jatom)))/Two
           her(i,j)=her(i,j)+y
           y=omega*(eta(i,j,n,jatom)-
     &        dconjg(eta(j,i,n,jatom)))/Two
           aher(i,j)=aher(i,j)+y
           enddo
          enddo
         enddo

         if (myid.eq.0) then
          write(6,*)
          write(6,*)'SELF-ENERGY ASSYMPTOSICS:',jatom
          write(6,*)
         endif


         do i=1,Nint
          do j=1,Nint
           etainf(i,j,jatom)=her(i,j)/ia
           eta1(i,j,jatom)=aher(i,j)/ia
          enddo
         enddo

         trS=0.d0
         do i=1,Nint
          trS=trS+dreal(etainf(i,i,jatom))
         enddo
        
        trS=trS/Nint
        if (licht) then
	write(6,*)':eDC:  ',-trS,' atom # ',jatom
         do k=1,kpt
          do i=1,Nint
           H(index0+i,index0+i,k)=H(index0+i,index0+i,k)-trS
          enddo
         enddo
	else
         do k=1,kpt
          do i=1,Nint
           H(index0+i,index0+i,k)=H(index0+i,index0+i,k)-edc(jatom)
          enddo
         enddo
        endif

        if (myid.eq.0) then
         write(6,*)'ETA(infinity):'
         do i=1,Nint
          write(6,5)(etainf(i,j,jatom),j=1,Nint)
         enddo
         write(6,*)
         write(6,*)'ETA(1):'
         do i=1,Nint
          write(6,5)(eta1(i,j,jatom),j=1,Nint)
         enddo
         write(6,*)
        endif
 100     continue
 5       format(5(2f10.5,4x))

         end
           

         


