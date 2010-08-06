c***********************************************************************
c Initializing particles.
      subroutine pinit(subcycle)
c Common data:
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include '3dcom.f'
      include 'meshcom.f'
      external linregion
      logical linregion
c Local dummy variables for partlocate.
      real xfrac(ndims_mesh)
      integer ixp(ndims_mesh)
      logical linmesh

      ntries=0
c      ntrapped=0
c      rmax=0.99999*rs
c      write(*,*) 'Pinit. rmax=',rmax,'xmeshend=',xmeshend
c      rmax2=rmax*rmax
c The maximum used slot is the same as the number of particles initially
      ioc_part=n_part
c     We initialize the 'true' particles'
      do i=1,n_part
         if_part(i)=1
 1       continue
         ntries=ntries+1
         x_part(1,i)=xmeshstart(1)+
     $        ran1(myid)*(xmeshend(1)-xmeshstart(1))
         x_part(2,i)=xmeshstart(2)+
     $        ran1(myid)*(xmeshend(2)-xmeshstart(2))
         x_part(3,i)=xmeshstart(3)+
     $        ran1(myid)*(xmeshend(3)-xmeshstart(3))
c     If we are not in the plasma region, try again.
         if(.not.linregion(ibool_part,npdim,x_part(1,i)))then
c            write(*,*)'Injection of',i,' wrong region',inewregion
c     $           ,(x_part(kk,i),kk=1,3)
            goto 1
         endif
         Ti0=Ti
         tisq=sqrt(Ti0)
         x_part(4,i)=tisq*gasdev(myid)
         x_part(5,i)=tisq*gasdev(myid)
         x_part(6,i)=tisq*gasdev(myid) + vd
         dtprec(i)=0.
c Initialize the mesh fraction data in x_part.
         call partlocate(i,ixp,xfrac,iregion,linmesh)
c This test rejects particles exactly on mesh boundary:
         if(.not.linmesh)goto 1
      enddo
c Set flag of unused slots to 0
      do i=n_part+1,n_partmax
         if_part(i)=0
      enddo
      write(*,*)'Initialized ','id=',myid,
     $     '  n=',n_part,'  ntries=',ntries
c Initialize rhoinf:
      if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rs**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo
c-----------------------------------------------------------------
c A special orbit.
c Pinit resets x_part. So set it for the special first particle.
      x_part(1,1)=2.
      x_part(2,1)=0.
      x_part(3,1)=0.
c Prior half-step radial velocity
      x_part(4,1)=0.5*subcycle*(abs(phip)/x_part(1,1)**2)
c      x_part(4,1)=0.5*dt*(abs(phip)/x_part(1,1)**2)
c      x_part(4,1)=0.25*dt*(abs(phip)/x_part(1,1)**2)
c Tangential velocity of circular orbit at r=?.
      x_part(5,1)=sqrt(abs(phip/x_part(1,1))-x_part(4,1)**2)
      x_part(6,1)=0.
      call partlocate(1,ixp,xfrac,iregion,linmesh)
      write(*,*)'PINIT',(x_part(i,1),i=1,9),phip
      end
c***********************************************************************
