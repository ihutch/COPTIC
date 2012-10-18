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
      if_part(1)=1
      i1=1
      call partlocate(1,ixp,xfrac,iregion,linmesh,nreloc)
      if(.not.linmesh.and.myid.eq.0)then
         write(*,*)'WARNING Special particle-1 outside region.'
         write(*,*)'Coordinates',(x_part(i,1),i=1,9),' Resetting.'
      endif
      if(linmesh)i1=2
c-----------------------------------------------------------------
      ntries=0
c      ntrapped=0
c The maximum used slot is the same as the number of particles initially
      ioc_part=n_part
c     We initialize the 'true' particles'
      do i=i1,n_part
         if_part(i)=1
 1       continue
c         write(*,'(i8,$)')i
         ntries=ntries+1
         do j=1,ndims_mesh
            x_part(j,i)=xmeshstart(j)+
     $        ran1(myid)*(xmeshend(j)-xmeshstart(j))
         enddo
c     If we are not in the plasma region, try again.
         if(.not.linregion(ibool_part,npdim,x_part(1,i)))then
c            write(*,*)'Initialization of',i,' wrong region',inewregion
c     $           ,(x_part(kk,i),kk=1,3)
            goto 1
         endif
         Ti0=Ti
         tisq=sqrt(Ti0)
         x_part(4,i)=tisq*gasdev(myid)
         x_part(5,i)=tisq*gasdev(myid)
         x_part(6,i)=tisq*gasdev(myid) + vd

c One idea that might improve the initialization would be to add the 
c (electric) potential energy here, by scaling up the total velocity,
c but keeping the direction the same. 

c However, to make this accurate, one would have to have solved the
c initial potential better than the simple Laplace solve currently
c used. A solution of the electron-shielded case might be nearer the
c mark, but it would not be correct when Ti<Te. It is perhaps
c advantageous because it constitutes a minimal estimate of the
c shielding, which is probably what one wants.

c One advantage of using the Laplace or electron-shielded case is that
c it enhances the kinetic energy. Thus if the potential rises towards
c zero during initial evolution, one has indeed avoided spuriously
c trapped ions. 

c Asymmetries are of course not accounted for. The particles in the
c close wake are in fact going to have mostly negative velocities. They
c are therefore completely unmodelled in this initialization, which
c fills in positive drifting ions even where they could not be present.

         dtprec(i)=0.
c Initialize the mesh fraction data in x_part.
         call partlocate(i,ixp,xfrac,iregion,linmesh,nreloc)
c This test rejects particles exactly on mesh boundary:
         if(.not.linmesh)goto 1
      enddo
c Set flag of unused slots to 0
      do i=n_part+1,n_partmax
         if_part(i)=0
      enddo
      if(myid.eq.0)then
       write(*,'('' Initialized '',i3,''x n='',i7,'' ntries='',i7,$)')
     $     nprocs,n_part,ntries
      endif
c Initialize rhoinf:
      if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rs**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo
      end
c***********************************************************************
