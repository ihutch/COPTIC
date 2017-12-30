      program diagexamine

      include 'examdecl.f'

      parameter (ndiagmax=8)
! diagmax here must be 7+1 to accommodate potential possibly.
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
! Volumes are stored in ndiagmax+1

! Extra work array for arrowplotting in sliceGweb.
      real vp(na_m,na_m2,3,3)
! 1-d plotting arrays.
      real z1d(na_m),u1d(na_m),dene1d(na_m),deni1d(na_m)
      real zminmax(2)
      character*20 mname(ndiagmax+1)

      integer mcell,icontour
      logical lvtk
      character*70 xtitle,ytitle,label
      integer iuphi(3),iurs(3),isrs(3)
      integer iunp,i1d,isingle,i1,iwr,istd
      data iunp/0/i1d/0/iwr/0/zminmax/0.,0./icontour/0/istd/1/
      data ierr/0/

      isingle=0
      lvtk=.false.
      fluxfilename=' '
      mname(1)='Cellcount'
      mname(2)='v!d1!d'
      mname(3)='v!d2!d'
      mname(4)='v!d3!d'
      mname(5)='T!d1!d'
      mname(6)='T!d2!d'
      mname(7)='T!d3!d'
      mname(ndiagmax)='Potential'
      mname(ndiagmax+1)='Density'
      xleg=.75
      ipp=0
      iworking=0
!      pscale=3.
! 
 1    continue
      call diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle,lvtk
     $     ,mcell,zminmax,icontour,iworking)
      if(iworking.lt.0)call exit(0)
      if(zminmax(1).gt.zminmax(2))then
         istd=0
      endif
      i1=1
      ied=ndiagmax
! Create label
      label=diagfilename(lentrim(diagfilename)-3:lentrim(diagfilename))
      call array3read(diagfilename,ifull,iuds,ied,diagsum,ierr)
      if(ierr.eq.1)stop 'Error reading diag file'

      ndiags=ied

! Attempt to read the volumes data.
      istat=1
      call stored3geometry(diagsum(1,1,1,ndiags+1),iuds,ifull,istat
     $     ,.false.)
      if(istd.gt.0)then
      if(istat.eq.1)then
         write(*,*)'Read volumes successfully'
     $     ,diagsum(2,2,2,ndiags+1)
      else
         write(*,*)'Storedgeom failed; returned',istat,iuds,ifull
         write(*,*)'***** No volume-corrected density is available.'
      endif
      endif
!-------------------------------------
! These are up to two additional quantities \phi and n for separate plotting
      if(lentrim(phifilename).gt.1)then
! Read in a separate potential as well.
         ied=1
         ierr=1
         call array3read(phifilename,ifull,iuphi,ied,u,ierr)
         do j=1,3
            if(iuphi(j).ne.iuds(j))then 
               write(*,*)'Potential array dimensions',iuphi
     $              ,'  incompatible with diagnostics',iuds
               stop
            endif
         enddo
      endif
      if(lentrim(denfilename).gt.1)then
! Read in a density as well.
         ied=1
         ierr=1
         call array3read(denfilename,ifull,iuphi,ied,q,ierr)
         do j=1,3
            if(iuphi(j).ne.iuds(j))then 
               write(*,*)'Density array dimensions',iuphi
     $              ,'  incompatible with diagnostics',iuds
               stop
            endif
         enddo
      endif
!-------------------------------------

      if(iunp.ne.0)then
! Unnormalized diagnostic plotting.
         do k=i1,ndiags
! Suppress help.
            zp(1,1,1)=99
            ifix=2
            write(fluxfilename,'(''diagsum('',i1,'')'',a1)')k,char(0)
            call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $           ixnp,xn,ifix,fluxfilename,dum,dum)
         enddo
      endif

! Normalize by dividing by density, which is first diagnostic.
! But as of 4 Dec 2013 the ndiags is phi, so don't normalize that.
      do id=2,ndiags-1
         do k=1,iuds(3)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  if(abs(diagsum(i,j,k,1)).gt.mcell+1.e-5)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)/diagsum(i,j,k
     $                    ,1)
                  else
                     diagsum(i,j,k,id)=0.
                  endif
! Subtract v^2 from the second moment to give temperature.             
                  if(id.gt.ndims+1)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)
     $                    -diagsum(i,j,k,id-3)**2
                  endif
               enddo
            enddo
         enddo
      enddo
      if(istd.gt.0)write(*,*)'rs,debyelen,vd,Ti',rs,debyelen,vd,Ti

! Normalize the density= diagsum(1)/volumes into the volumes 
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               if(diagsum(i,j,k,ndiags+1).ne.0.)then
                  diagsum(i,j,k,ndiagmax+1)
     $                 =diagsum(i,j,k,1)/diagsum(i,j,k,ndiags+1)
               else
                  diagsum(i,j,k,ndiagmax+1)=0.
!                  write(*,*)'Warning zero volumes',i,j,k
               endif
            enddo
         enddo
      enddo
!---------------------------------------------------------------
! Decide the actual ends and beginnings of the relevant data.
! If the middle of a face has zero density, that is a dummy face.
! e.g. because of periodicity.
      isrs(1)=1
      if(diagsum(1,iuds(2)/2,iuds(3)/2,1).eq.0)isrs(1)=2
      iurs(1)=iuds(1)+1-isrs(1)
      if(diagsum(iuds(1),iuds(2)/2,iuds(3)/2,1).eq.0)iurs(1)=iurs(1)-1
      isrs(2)=1
      if(diagsum(iuds(1)/2,1,iuds(3)/2,1).eq.0)isrs(2)=2
      iurs(2)=iuds(2)+1-isrs(2)
      if(diagsum(iuds(1)/2,iuds(2),iuds(3)/2,1).eq.0)iurs(2)=iurs(2)-1
      isrs(3)=1
      if(diagsum(iuds(1)/2,iuds(2)/2,1,1).eq.0)isrs(3)=2
      iurs(3)=iuds(3)+1-isrs(3)
      if(diagsum(iuds(1)/2,iuds(2)/2,iuds(3),1).eq.0)iurs(3)=iurs(3)-1

!      write(*,*)'isrs,iurs',isrs,iurs
!      write(*,*)'Normalized diagnostics.'

!      lvtk=.false.
      if(lvtk)then
! Write Visit-readable vtk file of density, velocity, and potential.
         if(fluxfilename(1:1).eq.' ')then
            fluxfilename='Velocity'//char(0)
         else
            write(*,*)'Variable name: '
     $           ,fluxfilename(1:lentrim(fluxfilename))
! Spaces are not allowed in visit data names. Fix:
            do i=1,lentrim(fluxfilename)
               if(fluxfilename(i:i).eq.' ')fluxfilename(i:i)='_'
            enddo
            call termchar(fluxfilename)
         endif
         ibinary=0
         call vtkwritescalar(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),1)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary
     $        ,'N'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,'Cell_Particles'//char(0))
         if(istat.eq.1)call vtkwritescalar(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),ndiagmax+1)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary
     $        ,'n'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,'Density'//char(0))
!         write(*,*)'Finished density vtkwrite',ifull,iuds,u(1,1,1)
         if(ndiags.ge.4)then
            call vtkwritevector(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),2)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary,3
     $        ,'V'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,fluxfilename)
            write(*,*)'Finished velocity vtkwrite',ifull,iuds,u(1,1,1)
! This is supposed to be the potential. It is the last array ndiags.
            call vtkwritescalar(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),ndiags)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary
     $        ,'u'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,'Potential'//char(0))
!         write(*,*)'Finished potential vtkwrite',ifull,iuds,u(1,1,1)
         endif
        stop
      endif


!-------------------------------------------------------------
      if(i1d.ne.0)then
! Calculate average profiles in direction i1d.
         do i=1,iuds(i1d)
            z1d(i)=0.
            u1d(i)=0.
            deni1d(i)=0.
         enddo
         write(*,*)'i1d,iuds',i1d,iuds
         do k=2,iuds(3)-1
            do j=2,iuds(2)-1
               do i=2,iuds(1)-1
                  ii=1
                  if(i1d.eq.2)ii=j
                  if(i1d.eq.3)ii=k
                  z1d(ii)=z1d(ii)+diagsum(i,j,k,1+i1d)
                  u1d(ii)=u1d(ii)+u(i,j,k)
                  deni1d(ii)=deni1d(ii)+q(i,j,k)
               enddo
            enddo
         enddo
         znorm=(iuds(1)-2)*(iuds(2)-2)*(iuds(3)-2)/(iuds(i1d)-2)
!         write(*,*)phifilename
         if(phifilename(1:1).ne.' ')write(*,*)'   z           v       '
     $        ,'    phi         n_i        z_E     4pi grad(phi)'
         if(iwr.ne.0)then
            if(phifilename(1:1).ne.' ')then
               write(*,*)iuds(i1d)-2,5
            else
               write(*,*)iuds(i1d)-2
            endif
         endif
! Fix up first and last value for gradient calculation
         u1d(1)=(2.*u1d(2)-u1d(3))/znorm
         u1d(iuds(i1d))=(2.*u1d(iuds(i1d)-1)-u1d(iuds(i1d)-2))/znorm
         do i=2,iuds(i1d)-1
            z1d(i)=z1d(i)/znorm
            u1d(i)=u1d(i)/znorm
            deni1d(i)=deni1d(i)/znorm
            dene1d(i)=exp(u1d(i))
!            write(*,*)i,u1d(i),dene1d(i),deni1d(i)
            if(iwr.ne.0)then
               if(phifilename(1:1).ne.' ')then
                  write(*,'(6f12.5)')xn(ixnp(i1d)+i),z1d(i),u1d(i)
     $                 ,deni1d(i)
     $                 ,(xn(ixnp(i1d)+i)+xn(ixnp(i1d)+i-1))/2.
     $                 ,(u1d(i)-u1d(i-1))*4.*3.14159
     $                 /(xn(ixnp(i1d)+i)-xn(ixnp(i1d)+i-1))
               else
                  write(*,*)xn(ixnp(i1d)+i),z1d(i)
               endif
            endif
         enddo
! Line-out plot.
         call fitinit(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),z1d(2)
     $        ,max(z1d(iuds(i1d)-1),0.))
         call polyline(xn(ixnp(i1d)+2),z1d(2),ixnp(i1d+1)-ixnp(i1d)-2)
         call axptset(1.,0.)
! xticoff reverses the tics.
         call xaxis(0.,0.)
         call ticset(-.015,.015,-.03,.007,0,0,0,0)
         call yaxis(0.,0.)
         call axlabels(xtitle(1:lentrim(xtitle))
     $        ,ytitle(1:lentrim(ytitle)))
         call axlabels('','v!di!d')
         call ticset(-.015,.015,-.03,.007,-1,0,0,0)
         call axptset(0.,1.)
         call altxaxis(1.,1.)
         call axptset(0.,1.)
         call ticset(0.,0.,0.,0.,0,0,0,0)


         if(phifilename(1:1).ne.' ')then
            call legendline(xleg,.15,0,'v!di!d')
            if(ipp.eq.1)then
               call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),u1d(2)
     $              ,u1d(iuds(i1d)-1),.false.,.false.)
               call axptset(-.2,0.)
               call color(4)
               call dashset(1)
               call altyaxis(1.,1.)
               call polyline(xn(ixnp(i1d)+2),u1d(2),ixnp(i1d+1)
     $              -ixnp(i1d)-2)
               call ticset(.015,.015,-.03,-.007,0,0,0,0)
               call axlabels('','!Af!@')
               call legendline(xleg,.2,0,'!Af!@')
            endif
            call axptset(0.,0.)
            call color(5)
            call dashset(2)
!            write(*,*)i1d,iuds(i1d)-1,u1d(2),dene1d(2),u1d(iuds(i1d)-1)
!     $           ,dene1d(iuds(i1d)-1)
            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),dene1d(2)
     $           ,dene1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),dene1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
!            call ticrev()
            call altyaxis(1.,1.)
            call ticset(.015,.015,-.03,-.007,0,0,0,0)
!            call ticset(0.,0.,0.,0.,0,0,0,0)
            call axlabels('','Density')
            call legendline(xleg,.05,0,'n!de!d')
            call ticset(0.,0.,0.,0.,0,0,0,0)
         endif
         if(denfilename(1:1).ne.' ')then
            call color(6)
            call dashset(6)
!            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),deni1d(2)
!     $           ,deni1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),deni1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
            call legendline(xleg,.1,0,'n!di!d')
         endif
         call pltend()
      elseif(lentrim(phifilename).gt.1)then
! Arrow plotting of velocity:
         fluxfilename=ytitle
         ifix=2+4
         call sliceGweb(ifull,iuds,u,na_m,zp,
     $        ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $        ,diagsum(1,1,1,2),vp)
      else
         if(istat.eq.1.and.isingle.eq.0)then
            k=ndiagmax+1
            zp(1,1,1)=99
            ifix=2
!         write(fluxfilename,'(''diagnorm('',i1,'')'')')k
            fluxfilename=mname(k)
            write(*,*)k,mname(k),fluxfilename(1:lentrim(fluxfilename))
            call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $           ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $           ,dum,dum)
         endif
!----------------------------------------------------------------
! Default examination of all diagnostics.
            call examineall(diagsum,zminmax,mname,i1,isingle,istd
     $     ,ndiags,icontour,label)
      endif
!----------------------------------------------------------------
! Arrow plotting of velocity:
!      fluxfilename=mname(1)
!      ifix=2+4
!      if(isingle.eq.0)call sliceGweb(ifull,iuds,diagsum(1,1,1,1),na_m,zp
!     $     ,ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
!     $     ,diagsum(1,1,1,2),vp)

      goto 1
      end


!*************************************************************
      subroutine diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle
     $     ,lvtk,mcell,zminmax,icontour,iworking)
      integer iunp,isingle,i1d
      integer mcell,icontour
      real zminmax(2)
      character*70 xtitle,ytitle
      logical lvtk
      include 'examdecl.f'

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

! silence warnings:
      zp(1,1,1)=0.
! Defaults
      diagfilename=' '
      phifilename=' '
      denfilename=' '
      ytitle=' '
      xtitle=' '
      mcell=5.
      ipfs=3

      write(*,*)'diagexamargs',iworking,iargc()
! Deal with arguments
      if(iargc().eq.0) goto 201
      do i=iworking+1,iargc()
         iworking=i
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:13).eq.'--objfilename')
     $           read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-p')then
               if(lentrim(phifilename).gt.1)then
                  read(argument(3:),'(a)',err=201)denfilename
               else
                  read(argument(3:),'(a)',err=201)phifilename
               endif
            endif
            if(argument(1:3).eq.'-ly')
     $           read(argument(4:),'(a)',err=201)ytitle
            if(argument(1:3).eq.'-lx')
     $           read(argument(4:),'(a)',err=201)xtitle
            if(argument(1:2).eq.'-u')iunp=1
            if(argument(1:2).eq.'-d')
     $           read(argument(3:),*,err=201)isingle
            if(argument(1:2).eq.'-o')iwr=1
            if(argument(1:2).eq.'-f')ipp=1
            if(argument(1:2).eq.'-g')then
               read(argument(3:),*,err=204,end=204)ipfs
 204           call pfset(ipfs)
            endif
            if(argument(1:2).eq.'-r')call noeye3d(0)
            if(argument(1:2).eq.'-c')
     $           read(argument(3:),*,err=201)icontour
            if(argument(1:2).eq.'-a')
     $           read(argument(3:),*,err=201)i1d
            if(argument(1:2).eq.'-m')
     $           read(argument(3:),*,err=201)mcell
            if(argument(1:2).eq.'-z')then
                read(argument(3:),*,err=201)zminmax
! Turn on trucation if the z-range is set
                call togi3trunc()
             endif
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
            if(argument(1:2).eq.'-w')then
               lvtk=.not.lvtk
               read(argument(3:),'(a)',err=201)fluxfilename
            endif
         else
            read(argument(1:),'(a)',err=201)diagfilename
            return
         endif         
      enddo
      iworking=-1
      return
      goto 202
!------------------------------------------------------------
! Help text
 201  continue
      write(*,*)'=====Error reading command line argument',argument
 203  continue
 301  format(a,i5)
 302  format(a,2f8.3)
      write(*,301)'Usage: diagexamine [switches] <diagfile>'
      write(*,*)' Plot diagnostics from file'
      write(*,*)' If additional parameter file is set,'
     $     ,' do arrow plot of velocity on it.'
      write(*,301)' -p<name>   set name of additional parameter file.'
     $     //' (Given <=twice: phi, den.)'
      write(*,301)' -ly  set label of parameter. -lx label of xaxis'
      write(*,301)' -d   set single diagnostic to be examined [',isingle
      write(*,301)' -a   set dimension number for ave profile [',i1d
      write(*,301)' -o   write out the profiles               [',iwr
      write(*,301)' -f   plot potential profile (if -p given) [',ipp
      write(*,301)' -u   plot un-normalized diagnostics       [',iunp
      write(*,301)' -c<i>  set contour plotting               ['
     $     ,icontour
      write(*,301)' -m<f>  set minimum non-zero cell count    [',mcell
      write(*,302)' -z<min><max> set range of values plotted  [',zminmax
      write(*,*)' if zmin>zmax, print zmin and zmax of first diag'
      write(*,*)'-g<i>   print ps-graphics files to unit i [3,-3]'
      write(*,*)'-r   run without pausing'
      write(*,*)'-w<name> toggle vtk file writing,'
     $     ,' optionally name the vector[',lvtk,' '
     $     ,fluxfilename(1:lentrim(fluxfilename))
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -h -?  Print usage and current switch values.'
      call exit(0)
 202  continue
      if(lentrim(diagfilename).lt.5)goto 203
      end
!*****************************************************************
      subroutine minmax3(ifull,iuds,u,umin,umax)
! Find the minimum and maximum of a 3d array. Allocated dimensions ifull,
! used dimensions iuds.
      parameter (ndims=3)
      integer ifull(ndims),iuds(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      umin=1.e30
      umax=-1.e30
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               if(u(i,j,k).gt.umax)umax=u(i,j,k)
               if(u(i,j,k).lt.umin)umin=u(i,j,k)
            enddo
         enddo
      enddo
      end
!*****************************************************************
!********************************************************************
      subroutine minmaxN(ndims,ifull,iuds,u,umin,umax)
! Find the minimum and maximum of a ndims array. Allocated dimensions ifull,
! used dimensions iuds. Using the general mditerator. 
      implicit none
      integer ndims
      integer ifull(ndims),iuds(ndims)
      real u(*),umin,umax
! Need local storage.
      integer ndimsmax,indexcontract,mditerator
      parameter (ndimsmax=5)
      integer indi(ndimsmax),iview(3,ndimsmax),i
      external mditerator,indexcontract
      umin=1.e30
      umax=-1.e30
      i=mditerator(ndims,iview,indi,4,iuds)
 1    i=1+indexcontract(ndims,ifull,indi)
         if(u(i).gt.umax)umax=u(i)
         if(u(i).lt.umin)umin=u(i)
      if(mditerator(ndims,iview,indi,0,iuds).eq.0)goto 1
      end
!*************************************************************
! This demonstrates the minimalistic iterator.
      subroutine minmaxM(ndims,ifull,iuds,u,umin,umax)
      integer ndims,ifull(ndims),iuds(ndims)
      real u(*),umin,umax
      integer mditer
      external mditer
      integer ii,i
      umin=1.e30
      umax=-1.e30
      ii=mditer(ndims,ifull,iuds,i)
 1       if(u(i).gt.umax)umax=u(i)
         if(u(i).lt.umin)umin=u(i)
      if(mditer(ndims,ifull,iuds,i).eq.0)goto 1
      end
!**************************************************************
      subroutine examineall(diagsum,zminmax,mname,i1,isingle,istd
     $     ,ndiags,icontour,label)
! Default examination of all diagnostics.
      implicit none
      include 'examdecl.f'
      integer ndiagmax
      parameter (ndiagmax=8)
! diagmax here must be 7+1 to accommodate potential possibly.
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
! Volumes are stored in ndiagmax+1
      real zminmax(2)
      integer i1,isingle,istd,ndiags,icontour
      character*20 mname(ndiagmax+1)
      character*70 label
      integer lentrim
      external lentrim

! Local variables
      real dum
      integer k,ifix
      do k=i1,ndiags
         zp(1,1,1)=99
         ifix=3
!            write(fluxfilename,'(''diagnorm('',i1,'')'')')k
         fluxfilename=mname(k)(1:lentrim(mname(k)))
     $        //'('//label(1:lentrim(label))//')'
!            write(*,*)k,isingle
         if(k.eq.ndiags)
     $        fluxfilename='!Af!@('//label(1:lentrim(label))//')'
         if(istd.gt.0.and.(k.eq.isingle.or.isingle.eq.0))write(*,*)k,
     $        fluxfilename(1:lentrim(fluxfilename))
         if(isingle.eq.0.or.isingle.eq.k)then
            if(zminmax(1).lt.zminmax(2))then
               write(*,*)'zminmax',zminmax
               call scale3(0.,1.,0.,1.,zminmax(1),zminmax(2))
               call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $              ixnp,xn,ifix+256+icontour*16+512
     $              ,fluxfilename(1:lentrim(fluxfilename)+2) ,dum
     $              ,dum)   
            elseif(zminmax(1).gt.zminmax(2))then
! If limits are crossed on entry. Find minimum and maximum.
               call minmaxM(ndims,ifull,iuds,diagsum(1,1,1,k),
     $              zminmax(1),zminmax(2))
               write(*,*)zminmax,'=zmin/max'
               call exit(0)
            else
               call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $              ixnp,xn,ifix+icontour*16+512
     $              ,fluxfilename(1:lentrim(fluxfilename)+2) ,dum
     $              ,dum)
            endif
         endif
      enddo
      end
