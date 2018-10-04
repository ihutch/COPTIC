
!************************************************************************
! Routine to be passed to mditerarg, to store the volumes calculated
! from meshcom for unintersected points or by volintegrate.
! Can't be called for edge nodes.
! Does not work for objects that do not set the intersections,
! e.g. point-charge objects. 
      subroutine volnode(inc,ipoint,indi,mdims,iLs,iused,
     $     volumes,cij)
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real volumes(*)
      real cij(*)

      include 'ndimsdecl.f'
      include 'partcom.f'
      include '3dcom.f'
      external linregion
      logical linregion
      include 'meshcom.f'
      include 'ptchcom.f'
      real xm(ndims),xi(ndims),xp(ndims)
!      parameter (npoints=10000)
! Sobel numbers are more "favorable" if the number of them is 2^k
! where k>6 (for 3 dimensions). So make it a power of 2.
      parameter (npoints=2**13)
! Structure vector needed for finding adjacent u values.
      integer iLs(ndims+1)
      include 'myidcom.f'
      integer icall
      external insideall
      save icall
      data icall/0/

! Silence warnings with spurious iused acces.
      icb=iused(1)
      icall=icall+1
      if(mod(icall,1000).eq.0.and.myid.eq.0)write(*,'(''.'',$)')
! The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
! Object pointer
      icij=icb*(ipoint+1)
      if(cij(icij).ne.0.) goto 1
! See if there is a chance the volume is intersected by a boundary.
      do id=1,ndims
         icij=icb*(ipoint+iLs(id)+1)
         if(cij(icij).ne.0.) goto 1
         icij=icb*(ipoint-iLs(id)+1)
         if(cij(icij).ne.0.) goto 1         
      enddo
! This is an unintersected case. Calculate simply.
      vol=1.
      do id=1,ndims
         index=ixnp(id)+indi(id)+1
         xi(id)=xn(index)
         vol=vol*(xn(index+1)-xn(index-1))*0.5
      enddo
      if(.not.linregion(ibool_part,ndims,xi)) goto 3
      volumes(ipoint+1)=vol
      inc=1
      return

! Intersected case.
 1    continue
      do id=1,ndims
         index=ixnp(id)+indi(id)+1
         xp(id)=xn(index+1)
         xi(id)=xn(index)
         xm(id)=xn(index-1)
      enddo
! If we are outside the active region use unintersected case
      if(.not.linregion(ibool_part,ndims,xi)) goto 3

!      write(*,*)'Volintegrate call:',indi,xi,cij(icij)
! Use volintegrate function.
      volumes(ipoint+1)=volintegrate(ndims,xm,xi,xp,npoints)
      inc=1
      return

! Outside region. Set a large volume. 10^30 if we are outside
! everything.  10^20 if we are inside a pointcharge region, in which
! case we (later) set the density to 0 rather than faddu alone.
 3    continue
      iregion=insideall(ndims,xi)
      if(IAND(iregion,iptch_mask).ne.0)then
! We are inside a point-charge region.
         volumes(ipoint+1)=1.e20
      else
         volumes(ipoint+1)=1.e30
      endif

      end
!********************************************************************
! Volume integrations by monte-carlo.  
! New version using Sobel quasi-random numbers from toms659.
! Using intrinsic IEOR it is faster than using rand. And less noisy.
! For a node i in the mesh the
! cic-volume is the integral of 1-|f| from x_{i-1}=xm to x_{i+1}=xp,
! where f is the mesh-fractional distance from x_i.  This integration is
! performed in each dimension.  To do this integration using monte-carlo
! techniques, for each direction, we choose a random position uniformly
! between x_{i-1} and x_{i+1}; (not uniform in mesh-fraction, although a
! scheme could be constructed uniform in mesh-fraction, which would
! probably be slightly less efficient); we weight by
! (1-|f|)(x_{i+1}-x_{i-1}). Then if it is in the region we add it on. If
! not we throw away. Eventually we divide by the total number of points
! examined.
! Old version just uncomment rand() [and comment out sobel calls].
! But randc.c would have to be reinstated. See end of file.

      real function volintegrate(mdims,xm,xi,xp,npoints)

      real xm(mdims),xi(mdims),xp(mdims)
      integer npoints
      double precision qrn(mdims)
      include 'ndimsdecl.f'
      include 'partcom.f'
      include '3dcom.f'
      real x(ndims)
      external linregion
      logical linregion
!     .. Scalar Arguments ..
      INTEGER ATMOST,DIMEN,TAUS
!     ..
!     .. Array Arguments ..
      LOGICAL FLAG(2)
!     ..

      ATMOST=npoints
      DIMEN=ndims
      wtot=0.
      vmax=0.
      vmin=1.e30
      call INSOBL(FLAG,DIMEN,ATMOST,TAUS)
      do i=1,npoints
         call GOSOBL(qrn)
         w=1.
         do id=1,ndims
            p=real(qrn(id))
!            p=rand()
            x(id)=xp(id)*p+xm(id)*(1-p)
            f=x(id)-xi(id)
            if(f.lt.0)then
               f=f/(xm(id)-xi(id))
            else
               f=f/(xp(id)-xi(id))
            endif
            w=w*(1.-f)*(xp(id)-xm(id))
         enddo
         if(linregion(ibool_part,ndims,x))wtot=wtot+w
!         if(npoints-i.lt.100)then
!            v=wtot/i
!            if(v.gt.vmax)vmax=v
!            if(v.lt.vmin)vmin=v
!         endif
      enddo
      volintegrate=wtot/npoints
!      write(*,'(i6,f8.5,'' +-fraction'',f7.4,f7.4)')npoints,volintegrate
!     $     ,(vmax-volintegrate)/vmax,(volintegrate-vmin)/vmax
      end
!**********************************************************************
! Read and/or write geometric information.
! On entry 
!   volumes, iuds, ifull are the current values
!   istat indicates 0: Write out data. 1: Try to read data.
! On exit
!   volumes contains read data if successful
!   istat indicates 0: No file exists (we failed if trying to read)
!                   1: File exists (we succeeded if trying to read).
! Thus two subsequent calls to this routine make sense: the first
! with istat=1, and the second with istat unchanged. That will read in
! the data if it exists and write it out if it doesn't. Between the 
! two calls the data should be generated if istat=0 after the first
! read.
! Normally this is called with lstrict=.true. but if not, then
! the checks for object identities are skipped.
      subroutine stored3geometry(volumes,iuds,ifull,istat,lstrict)
      include 'ndimsdecl.f'
      integer iuds(ndimsmax),ifull(ndimsmax)
      real volumes(ifull(1),ifull(2),ifull(3))
      logical lstrict
      integer iuds1(ndimsmax)
      include '3dcom.f'
      include 'meshcom.f'
      real obj1(odata,ngeomobjmax)
      real xn1(ixnlength)
      parameter (iunit=14)

! Use istat to decide action.
      if(istat.eq.0)goto 10
      if(istat.eq.1)goto 11
! Else do nothing
      return

! Write out the current geometric data.
 10   open(file='storedgeom.dat',unit=iunit,status='unknown',err=101)
      close(iunit,status='delete')
      open(file='storedgeom.dat',unit=iunit,status='new'
     $     ,form='unformatted',err=101)
      write(iunit)ndims
      write(iunit)ngeomobj
!      write(iunit)((obj_geom(i,j),i=1,oabc-1),j=1,ngeomobj)
      write(iunit)((obj_geom(i,j),i=1,odata),j=1,ngeomobj)
      write(iunit)(iuds(i),i=1,ndims)
      write(iunit)(xn(i),i=1,ixnp(ndims+1))
      write(iunit)
     $     (((volumes(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(iunit)
      istat=1
      write(*,*)'Successful storedgeom write completed.'
      return
! Read in the old geometric data if it exists.
! 
 11   open(file='storedgeom.dat',unit=iunit,status='old'
     $     ,form='unformatted',err=101)
      read(iunit,err=101,end=102)ndims1
      if(lstrict .and. ndims1.ne.ndims) goto 101
      read(iunit,err=101,end=102)ngeomobj1
      if(lstrict.and. ngeomobj1.ne.ngeomobj) goto 101
      read(iunit,err=101,end=102)((obj1(i,j),i=1,odata),j=1,ngeomobj)
      do j=1,ngeomobj1
         do i=1,odata
            if(  lstrict .and.
     $           (i.le.otype .or. i.ge.ocenter) .and.
     $           (obj1(i,j).ne.obj_geom(i,j))
     $        ) goto 101
         enddo
      enddo
      read(iunit,err=101,end=102)(iuds1(i),i=1,ndims)
      if(iuds1(1).ne.iuds(1)) goto 101
      if(iuds1(2).ne.iuds(2)) goto 101
      if(iuds1(3).ne.iuds(3)) goto 101
      read(iunit)(xn1(i),i=1,ixnp(ndims+1))
      do i=1,ixnp(ndims+1)
         if(xn(i).ne.xn1(i))goto 101
      enddo
      read(iunit,err=101,end=102)
     $     (((volumes(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(iunit)
      istat=1
!      write(*,*)'Successful storedgeom read completed.'
! Successful read.
      return
      
 102  write(*,*)'End-file error, reading storedgeom.dat data.'
! No existing file.
 101  istat=0
      if(lstrict)then
         close(iunit,status='delete')
      else
         close(iunit)
      endif
      return

      end
!*********************************************************************
! Obsolete c function:
!$$$#include <stdlib.h>
!$$$#include <stdio.h>
!$$$
!$$$float rand_()
!$$${
!$$$  long i;
!$$$  float x;
!$$$  static double xfac;
!$$$  static long il=0;
!$$$  if(!il){
!$$$    il=(long)RAND_MAX;
!$$$    xfac=1./( ((double) il) + (((double) il)/1000000.));
!$$$  }
!$$$  i = random(); 
!$$$  x = ((double) i)*xfac ;
!$$$  if(x >= 1. || x < 0.) {
!$$$    printf("RAND Error: x=%f, i=%ld, il=%ld\n",x,i,il);
!$$$  }
!$$$  return x;
!$$$}
!$$$
!$$$void srand_(long *iseed)
!$$${
!$$$  srandom( (unsigned int) *iseed); 
!$$$}
!$$$/*******************************************************************/
