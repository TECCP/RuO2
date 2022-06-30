        subroutine delta_as(xm2,xm3,xm4,Gw,jatom)
        include 'dmft.dat'

        complex*16 xm2(Niter,Niter,nat),xm3(Niter,Niter,nat),
     &             xm4(Niter,NIter,nat),Gw(niter,niter,0:momega),
     &             y,omega

        ia=10
        beta_pi=Pi/Beta
        Nint=n_int(jatom)
        do i=1,Nint
         do j=1,Nint
          xm3(i,j,jatom)=zero
          xm4(i,j,jatom)=zero
          do n=Im-ia+1,Im
           om=(2*n+1)*beta_pi
           omega=cim*om
           y=(Gw(i,j,n)+dconjg(Gw(j,i,n)))/Two
           y=y*omega**2
           xm3(i,j,jatom)=xm3(i,j,jatom)+y
           y=(Gw(i,j,n)-dconjg(Gw(j,i,n)))/Two
           y=(y*omega-xm2(i,j,jatom))*(omega**2)
           xm4(i,j,jatom)=xm4(i,j,jatom)+y
          enddo
          xm3(i,j,jatom)=xm3(i,j,jatom)/ia
          xm4(i,j,jatom)=xm4(i,j,jatom)/ia
         enddo
         write(6,*)'i:',xm3(i,i,jatom),xm4(i,i,jatom)
        enddo
        end
           
          
          
       
         
