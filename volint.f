
c************************************************************************
c Routine to be passed to mditerarg, to store the volumes calculated
c from meshcom for unintersected points or by volintegrate.
c Can't be called for edge nodes.
c Does not work for objects that do not set the intersections,
c e.g. point-charge objects. 
      subroutine volnode(inc,ipoint,indi,ndims,iused,
     $     volumes,cij)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real volumes(*)
      real cij(*)

      include '3dcom.f'
      external linregion
      logical linregion
      include 'meshcom.f'
      real xm(ndims_mesh),xi(ndims_mesh),xp(ndims_mesh)
      parameter (npoints=10000)
c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs
      include 'myidcom.f'
      integer icall
      save icall
      data icall/0/

c Silence warnings with spurious iused acces.
      icb=iused(1)
      icall=icall+1
      if(mod(icall,300).eq.0.and.myid.eq.0)write(*,'(''.'',$)')
c The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
c Object pointer
      icij=icb*(ipoint+1)
      if(cij(icij).ne.0.) goto 1
c See if there is a chance the volume is intersected by a boundary.
      do id=1,ndims
         icij=icb*(ipoint+iLs(id)+1)
         if(cij(icij).ne.0.) goto 1
         icij=icb*(ipoint-iLs(id)+1)
         if(cij(icij).ne.0.) goto 1         
      enddo
c This is an unintersected case. Calculate simply.
      vol=1.
      do id=1,ndims_mesh
         index=ixnp(id)+indi(id)+1
         xi(id)=xn(index)
         vol=vol*(xn(index+1)-xn(index-1))*0.5
      enddo
      if(.not.linregion(ibool_part,ndims_mesh,xi)) goto 3
      volumes(ipoint+1)=vol
      inc=1
      return

c Intersected case.
 1    continue
      do id=1,ndims_mesh
         index=ixnp(id)+indi(id)+1
         xp(id)=xn(index+1)
         xi(id)=xn(index)
         xm(id)=xn(index-1)
      enddo
c If we are outside the active region use unintersected case
      if(.not.linregion(ibool_part,ndims_mesh,xi)) goto 3

c      write(*,*)'Volintegrate call:',indi,xi,cij(icij)
c Use volintegrate function.
      volumes(ipoint+1)=volintegrate(ndims,xm,xi,xp,npoints)
      inc=1
      return

c Outside region. Set a large volume
 3    volumes(ipoint+1)=1.e30

      end
c********************************************************************
c Volume integrations by monte-carlo.  For a node i in the mesh the
c cic-volume is the integral of 1-|f| from x_{i-1}=xm to x_{i+1}=xp,
c where f is the mesh-fractional distance from x_i.  This integration is
c performed in each dimension.  To do this integration using monte-carlo
c techniques, for each direction, we choose a random position uniformly
c between x_{i-1} and x_{i+1}; (not uniform in mesh-fraction, although a
c scheme could be constructed uniform in mesh-fraction, which would
c probably be slightly less efficient); we weight by
c (1-|f|)(x_{i+1}-x_{i-1}). Then if it is in the region we add it on. If
c not we throw away. Eventually we divide by the total number of points
c examined.

c      real function volintegrate(ndims,xm,xi,xp,iregion,npoints)
      real function volintegrate(ndims,xm,xi,xp,npoints)

      real xm(ndims),xi(ndims),xp(ndims)
      integer npoints

      parameter (mdims=10)
      real x(mdims)
      include '3dcom.f'
      external linregion
      logical linregion

      wtot=0.
      do i=1,npoints
         w=1.
         do id=1,ndims
c Using the improved random number generator seems a significant hit on
c time here. So use the direct C rand which is perhaps quicker.
            p=rand()
            x(id)=xp(id)*p+xm(id)*(1-p)
            f=x(id)-xi(id)
            if(f.lt.0)then
               f=f/(xm(id)-xi(id))
            else
               f=f/(xp(id)-xi(id))
            endif
            w=w*(1.-f)*(xp(id)-xm(id))
         enddo
c         irg=insideall(ndims,x)
         if(linregion(ibool_part,ndims,x))wtot=wtot+w
      enddo
      volintegrate=wtot/npoints

      end
c**********************************************************************
c Read and/or write geometric information.
c On entry 
c   volumes, iuds, ifull are the current values
c   istat indicates 0: Write out data. 1: Try to read data.
c On exit
c   volumes contains read data if successful
c   istat indicates 0: No file exists (we failed if trying to read)
c                   1: File exists (we succeeded if trying to read).
c Thus two subsequent calls to this routine make sense: the first
c with istat=1, and the second with istat unchanged. That will read in
c the data if it exists and write it out if it doesn't. Between the 
c two calls the data should be generated if istat=0 after the first
c read.
      subroutine stored3geometry(volumes,iuds,ifull,istat)
      parameter (ndims=3)
      integer iuds(ndims),ifull(ndims)
      real volumes(ifull(1),ifull(2),ifull(3))
      integer iuds1(ndims)
      include '3dcom.f'
      include 'meshcom.f'
      real obj1(odata,ngeomobjmax)
      real xn1(1000)
      parameter (iunit=14)

c Use istat to decide action.
      if(istat.eq.0)goto 10
      if(istat.eq.1)goto 11
c Else do nothing
      return

c Write out the current geometric data.
 10   open(file='storedgeom.dat',unit=iunit,status='unknown',err=101)
      close(iunit,status='delete')
      open(file='storedgeom.dat',unit=iunit,status='new'
     $     ,form='unformatted',err=101)
      write(iunit)ndims
      write(iunit)ngeomobj
c      write(iunit)((obj_geom(i,j),i=1,oabc-1),j=1,ngeomobj)
      write(iunit)((obj_geom(i,j),i=1,odata),j=1,ngeomobj)
      write(iunit)(iuds(i),i=1,ndims)
      write(iunit)(xn(i),i=1,ixnp(ndims_mesh+1))
      write(iunit)
     $     (((volumes(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(iunit)
      istat=1
      write(*,*)'Successful storedgeom write completed.'
      return
c Read in the old geometric data if it exists.
c 
 11   open(file='storedgeom.dat',unit=iunit,status='old'
     $     ,form='unformatted',err=101)
      read(iunit,err=101,end=102)ndims1
      if(ndims1.ne.ndims) goto 101
      read(iunit,err=101,end=102)ngeomobj1
      if(ngeomobj1.ne.ngeomobj) goto 101
c      read(iunit,err=101,end=102)((obj1(i,j),i=1,oabc-1),j=1,ngeomobj)
      read(iunit,err=101,end=102)((obj1(i,j),i=1,odata),j=1,ngeomobj)
      do j=1,ngeomobj
         do i=1,odata
            if(
     $           (i.le.otype .or. i.ge.ocenter) .and.
     $           (obj1(i,j).ne.obj_geom(i,j))
     $        ) goto 101
         enddo
      enddo
      read(iunit,err=101,end=102)(iuds1(i),i=1,ndims)
      if(iuds1(1).ne.iuds(1)) goto 101
      if(iuds1(2).ne.iuds(2)) goto 101
      if(iuds1(3).ne.iuds(3)) goto 101
      read(iunit)(xn1(i),i=1,ixnp(ndims_mesh+1))
      do i=1,ixnp(ndims_mesh+1)
         if(xn(i).ne.xn1(i))goto 101
      enddo
      read(iunit,err=101,end=102)
     $     (((volumes(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(iunit)
      istat=1
c      write(*,*)'Successful storedgeom read completed.'
c Successful read.
      return
      
 102  write(*,*)'End-file error, reading storedgeom.dat data.'
c No existing file.
 101  istat=0
      close(iunit,status='delete')
      return

      end
