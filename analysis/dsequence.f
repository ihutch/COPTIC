      program dsequence
! Stripped down diagexamine to facilitate plotting of chosen slices.

      include 'examdecl.f'
      parameter (ndiagmax=8)
! diagmax here must be 7+1 to accommodate potential possibly.
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
! Volumes are stored in ndiagmax+1

! Extra work array for arrowplotting in sliceGweb.
      real vp(na_m,na_m2,3,3)
! 1-d plotting arrays.
      real z1d(na_m),u1d(na_m),dene1d(na_m),deni1d(na_m)
      real zbase(na_m)
      real zminmax(2),xminmax(2)
      character*20 mname(ndiagmax+1)

      integer mcell,icontour,nskip
      logical lvtk
      character*70 xtitle,ytitle,label
      integer iuphi(3),iurs(3),isrs(3)
      integer iunp,i1d,isingle,i1,iwr,istd
      data iunp/0/i1d/0/iwr/0/zminmax/0.,0./icontour/0/istd/1/
      data ierr/0/xminmax/0.,0./nskip/1/jbase/0/


      lvtk=.false.
      fluxfilename=' '
      mname(1)='Cellcount'
      mname(2)='v!d1!d'
      mname(3)='v!d2!d'
      mname(4)='v!d3!d'
      mname(5)='T!d1!d'
      mname(6)='T!d2!d'
      mname(7)='T!d3!d'
      xleg=.75
! 
      call dsequenceargs(iunp,isingle,idfix,iwr,ipp,xtitle,ytitle,lvtk
     $     ,mcell,zminmax,xminmax,icontour,nskip,jbase)
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
      mname(ndiags)='Potential'
      mname(ndiags+1)='Density'
!----------------------------------------------------------------
      k=isingle
      if(isingle.eq.0)k=ndiags
      if(ipp.eq.0)ipp=(iuds(idfix)+1)/2
      idp1=mod(idfix,3)+1
      idp2=mod(idfix+1,3)+1
!----------------------------------------------------------------
! Get a slice at fixed dimension node n1, default middle. 
      n1=ipp
      nf1=iuds(idp1)-1
      nf2=iuds(idp2)-1
      nff=iuds(idfix)-1
      if1=1+1
      if2=1+1
      iff=1+1
      write(*,*)idfix,n1,nf1,nf2,idp1,idp2
! We might wish to swap the indices order for plotting here.
      id1=1
      do j=if2,nf2
         do i=if1,nf1
            if(idfix.eq.1)then
               zp(i,j,id1)=diagsum(n1,i,j,k)
            elseif(idfix.eq.2)then
               zp(i,j,id1)=diagsum(j,n1,i,k)
            elseif(idfix.eq.3)then
               zp(i,j,id1)=diagsum(i,j,n1,k)
            endif
         enddo
      enddo
      if(jbase.gt.0)then
! Subtract a base case from all
         do i=if1,nf1
            zbase(i)=zp(i,jbase,id1)
         enddo
         do j=if2,nf2
            do i=if1,nf1
               zp(i,j,id1)=zp(i,j,id1)-zbase(i)
            enddo
         enddo
      endif
!-------------------------------------------------------------
! Start of plotting.

! Possibly find vertical range.
      if(zminmax(1).eq.zminmax(2))call minmax2(zp(if1,if2,id1)
     $     ,na_m,nf1+1-if1,nf2+1-if2,zminmax(1),zminmax(2))
      if(xminmax(1).eq.xminmax(2))then
         xminmax(1)=xn(ixnp(idp1)+if1)
         xminmax(2)=xn(ixnp(idp1)+nf1)
      endif
! Stacked plot of all profiles.
      call pltinit(xminmax(1),xminmax(2),zminmax(1),zminmax(2))
      call axis()
      call axis2()
      call axlabels('x',mname(k)(1:lentrim(mname(k))))
      call winset(.true.)
      do j=if2,nf2,nskip
         jc=1+(j-if2)/nskip
!               write(*,*)j,jc
         if((nf2-if2+1)/nskip.lt.16)call color(jc)
         call polyline(xn(ixnp(idp1)+if1),zp(if1,j,id1),nf1+1
     $        -if1)
      enddo
      call pltend()

      call pltinit(0.,1.,0.,1.)
      if(.false.)then
         call puteye(0.0,.1,0.)
         call axon(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),zp(if1
     $        ,if2,id1),na_m,nf1+1-if1,nf2+1-if2)
      endif
!       Plot the surface. With axes (1). Web color 10, axis color 7.
      j=1 + 256*10 + 256*256*7
!            call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),zp(if1,if2
!     $           ,id1),na_m,nf1+1-if1,nf2+1-if2,j)
!            call pltend()

      end


!*************************************************************
      subroutine dsequenceargs(iunp,isingle,idfix,iwr,ipp,xtitle,ytitle
     $     ,lvtk,mcell,zminmax,xminmax,icontour,nskip,jbase)
      integer iunp,isingle,idfix
      integer mcell,icontour
      real zminmax(2),xminmax(2)
      character*70 xtitle,ytitle
      logical lvtk
      include 'examdecl.f'

         ifull(1)=na_i
         ifull(2)=na_j
         ifull(3)=na_k

! silence warnings:
      zp(1,1,1)=0.
! Defaults
      isingle=0
      idfix=3
      ipp=0
      diagfilename=' '
      phifilename=' '
      denfilename=' '
      ytitle=' '
      xtitle=' '
      mcell=5.
      ipfs=3

! Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
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
            if(argument(1:2).eq.'-f')read(argument(3:),*,err=201)ipp
            if(argument(1:2).eq.'-g')then
               read(argument(3:),*,err=204,end=204)ipfs
 204           call pfset(ipfs)
            endif
            if(argument(1:2).eq.'-r')call noeye3d(0)
            if(argument(1:2).eq.'-c')
     $           read(argument(3:),*,err=201)icontour
            if(argument(1:2).eq.'-a')read(argument(3:),*,err=201)idfix
            if(argument(1:2).eq.'-m')read(argument(3:),*,err=201)mcell
            if(argument(1:2).eq.'-z')then
               read(argument(3:),*,err=201)zminmax
! Turn on trucation if the z-range is set
               call togi3trunc()
            endif
            if(argument(1:2).eq.'-x')read(argument(3:),*,err=201)xminmax
            if(argument(1:2).eq.'-s')read(argument(3:),*,err=201)nskip
            if(argument(1:2).eq.'-b')read(argument(3:),*,err=201)jbase
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
            if(argument(1:2).eq.'-w')then
               lvtk=.not.lvtk
               read(argument(3:),'(a)',err=201)fluxfilename
            endif
         else
            read(argument(1:),'(a)',err=201)diagfilename
         endif
         
      enddo
      goto 202
!------------------------------------------------------------
! Help text
 201  continue
      write(*,*)'=====Error reading command line argument',argument
 203  continue
 301  format(a,i5)
 302  format(a,2f8.3)
      write(*,301)'Usage: dsequence [switches] <diagfile>'
      write(*,*)' Plot diagnostics from file'
      write(*,301)' -ly  set label of parameter. -lx label of xaxis'
      write(*,301)' -d   set single diagnostic to be examined [',isingle
      write(*,301)' -a   set dimension number of plane normal [',idfix
!      write(*,301)' -o   write out the profiles               [',iwr
      write(*,301)' -f   set fixed plane number               [',ipp
      write(*,301)' -s   set skip size in plot sequence       [',nskip
      write(*,301)' -b   set base case to subtract            [',jbase
!      write(*,301)' -u   plot un-normalized diagnostics       [',iunp
!      write(*,301)' -m<f>  set minimum non-zero cell count    [',mcell
      write(*,302)' -z<min><max> set range of values plotted  [',zminmax
      write(*,302)' -x<min><max> set range of values plotted  [',xminmax
      write(*,*)' if zmin>zmax, print zmin and zmax of first diag'
      write(*,*)'-g<i>   print ps-graphics files to unit i [3,-3]'
      write(*,*)'-r   run without pausing'
      write(*,*)'-w<name> toggle vtk file writing,'
     $     ,' optionally name the vector[',lvtk,' '
     $     ,fluxfilename(1:lentrim(fluxfilename))
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
