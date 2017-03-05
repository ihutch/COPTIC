!****************************************************************
! Routine for setting cij, which is used by the general mditerarg.
      subroutine cijroutine(inc,ipoint,indi,mdms,iLs,iused,cij,debyelen
     $     ,error)
! This routine sets the object pointer in cij if there is an object
! crossing next to the the point at indi(ndims) in the 
! dimension id (plus or minus), adjusts the cij values,
! and sets the object data into the place pointed to in obj. 
! mdms is dummy here because ndimsdecl is used.
      include 'ndimsdecl.f'
! Effective index in dimension, c-style (zero based)
!      integer indm1(ndims)
! Shifted to correct for face omission.
      integer indi(ndims),iused(ndims)
! Error indicator
      real error
! Not used in this routine
      real cij(*)
      integer iLs(ndims+1)
!--------------------------------------------------
! Object-data storage.
      include 'objcom.f'
! Object information
      include '3dcom.f'
!-----------------------------------------------------
      real dpm(2),fraction(2),dplus(2),deff(2)
      real conditions(3,2)
      real tiny
      parameter (tiny=1.e-15)
      real fn(ndims)
      integer ipt2,oihere
!      integer ipa(ndims)
!      logical newbox
      integer istart,ireport
      integer ijbin(ovlen)
      integer ipb(ndims)
      data ipb/ndims*-1/
      data istart/0/ireport/0/

! Default no chain.
      ichain=0
! ipoint here is the offset (zero-based pointer)
! The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
      icij=icb*(ipoint+1)
! Object pointer defaults to zero, unset initially.
! We are setting only trailing adjacent pointers in the loop below.
! Therefore we know that this cij(icij) has not yet been set.
      if(cij(icij).ne.0)then
! When we are sure this is working, we can simplify to cij(icij)=0.
         if(ireport.eq.0)write(*,*)'cij(icij) needed initialization'
     $        ,icij,cij(icij)
         ireport=ireport+1
         cij(icij)=0
      endif
! Prevent divide by zero issues with debyelen
      debyehere=debyelen
      if(debyehere.lt.tiny)debyehere=tiny
! New boxedge management section (newbox)
      idiag=0
! Look at the trailing box relative to point (since ipb is -1).
! Report any intersections as npoints. 
      call boxedge(ndims,ipb,indi,fn,npoints,idiag)
! There was doubt about what value of npoints to require, >0 is tight.
      if(npoints.gt.0)then
! Do over all the corners of this box
         do ic=0,2**ndims-1
            ipt2=ipoint
            do ib=0,ndims-1
! Avoiding the lower edge regions, which are boundaries.
! Upper boundaries are not set because they aren't trailing.
               if(indi(ib+1).le.1)goto 31
               if(btest(ic,ib))ipt2=ipt2-iLs(ib+1)
            enddo
            icij2=icb*(ipt2+1)
! Conditionally start the object: only if it does not already exist.
            call objstart(cij(icij2),istart,ipt2)
 31        continue
         enddo
      else
         cij(icij)=0.
      endif
      
! It is impossible now to assume this object is oi_cij because we might
! have started an addtional obj. So don't
      if(cij(icij).eq.0.)then
         oihere=oi_cij
      else
         oihere=int(cij(icij))
      endif
! Iterate over dimensions.
      do id=1,ndims
         ifound=0
! Find the intersection, if any, of each mesh branch with the objects.
! Return the fraction of the mesh distance of the intersection, and 
! the value of a, b, c to be attached to it. Fraction=1. is the 
! default case, where no intersection occurs. However, we need to know
! the fractions for opposite directions before we can calculate the
! coefficients; so we have to loop twice over the opposite directions.
! For each direction in this dimension,
         do i=1,2
            ioad=ndata_cij*(2*(id-1)+(i-1))+1
            ipm=1-2*(i-1)
! Determine whether this is a boundary point: adjacent a fraction ne 1.
! This potlsect call might not need ijbin.
            call potlsect(id,ipm,ndims,indi,
     $           fraction(i),conditions(1,i),dpm(i),iobjno,ijbin(1))
            if(fraction(i).lt.1. .and. fraction(i).ge.0.)then
               ifound=ifound+1
! Start object data for this point if not already started.
               call objstart(cij(icij),istart,ipoint)
!---------------
               oihere=int(cij(icij))
! If this object's a,b,c depends on flux, start extra info.
               if(int(obj_geom(otype,iobjno))/256.eq.4 .or.
     $              int(obj_geom(otype,iobjno))/256.eq.5)then
                  if(idob_cij(iextra_cij,oihere).eq.0)then
! Start chained extra data. This really is in position oi_cij+ichain.
                     ichain=1
!                     write(*,*)'Chaining',ipoint,(indi(kw),kw=1,ndims)
!     $                    ,id,i,iobjno,ijbin(1),cij(icij)
! First time, initialize the extra pointer
                     idob_cij(iextra_cij,oihere)=oi_cij+ichain
! and initialize the data for all directions
                     do ie=1,nobj_cij
                        dob_cij(ie,oi_cij+ichain)=0.
                     enddo
                  endif
! Set the data for this direction.          
                  idob_cij(ioad,oi_cij+ichain)=iobjno
                  idob_cij(ioad+1,oi_cij+ichain)=ijbin(1)
! ioad+2 is set below from coef/a                 
               endif
!---------------
! Calculate dplus and deff for this direction.
! dplus becomes dminus for the other direction.
               a=conditions(1,i)
               b=conditions(2,i)/debyehere
               c=conditions(3,i)
               dxp1=fraction(i)*dpm(i)/debyehere
               if(a.eq.0.) a=tiny
               if(b.ge.0.)then
! Active
                  apb=a*dxp1+b
                  dplus(i)=dxp1*(apb+b)/apb
                  deff(i)=apb/a
               else
! Inactive
                  if(dxp1.eq.0.) dxp1=tiny
                  dxp2=(1.-fraction(i))*dpm(i)/debyehere
                  apb=a*dxp2+b
                  boapb=b/apb
!  Without \psi'', i.e. setting it to zero and boundary.
!                  dplus(i)=dxp1
!  With \psi'', i.e. extrapolating it uniformly over the boundary
                  dplus(i)=dxp1+boapb*dxp2**2/dxp1
                  deff(i)=dxp1
               endif
               if(deff(i).eq.0)deff(i)=tiny
! The object data consists of data enumerated as
! ndims*(2=forward/backward from point)*(3=fraction,b/a,c/a)
! + diagonal + potential terms.
! Prevent subsequent divide by zero danger.
               dob_cij(ioad,oihere)=max(fraction(i),tiny)
               if(a.eq.0.) a=tiny
               dob_cij(ioad+1,oihere)=b/a*debyehere
               dob_cij(ioad+2,oihere)=c/a
               idob_cij(iinter_cij,oihere)=iobjno
            else
! No intersection.
               dplus(i)=dpm(i)/debyehere
               deff(i)=dpm(i)/debyehere
            endif
         enddo
! Now the dplus and deff are set correctly for each direction.
! Calculate the coefficients for each direction, using the opposite
! dplus for the dminus.
         do i=1,2
            im=mod(i,2)+1
! coefficient Cd same for all cases
            coef=2./(deff(i)*(dplus(i)+dplus(im)))
            if(ifound.gt.0)then
! This is a boundary point.
               ioad=ndata_cij*(2*(id-1)+(i-1))+1
               if(dob_cij(ioad,oihere).lt.1
     $              .and.dob_cij(ioad,oihere).ge.0.)then
! We intersected an object in this direction. Adjust Cij and B_y
                  a=conditions(1,i)
                  b=conditions(2,i)/debyehere
                  c=conditions(3,i)
                  if(a.eq.0.) a=tiny
                  if(b.ge.0.)then
! Active side.
                     cij(icb*ipoint+2*(id-1)+i)=0.
! Diagonal term (denominator) difference Cd-Cij stored
                     dob_cij(idgs_cij,oihere)=dob_cij(idgs_cij,oihere)
     $                 + coef
! Adjust potential sum (B_y)
                     dob_cij(ibdy_cij,oihere)=dob_cij(ibdy_cij,oihere)
     $                    - coef*c/a
                     if(ichain.gt.0)dob_cij(ioad+2,oihere+ichain)=coef/a
                  else
! Inactive side. Continuity.
                     dxp2=(1.-fraction(i))*dpm(i)/debyehere
                     apb=a*dxp2+b
                     boapb=b/apb
                     cij(icb*ipoint+2*(id-1)+i)=coef*boapb
! Diagonal term (denominator) difference Cd-Cij stored.
                     dob_cij(idgs_cij,oihere)=dob_cij(idgs_cij,oihere)
     $                 + coef*(1.-boapb)
! Adjust potential sum (B_y)
                     dob_cij(ibdy_cij,oihere)=dob_cij(ibdy_cij,oihere)
     $                    - coef*c*dxp2/apb
                     if(ichain.gt.0)dob_cij(ioad+2,oi_cij+ichain)=coef
     $                    *dxp2/apb
                  endif
               else
! We did not intersect an object in this direction.
                  cij(icb*ipoint+2*(id-1)+i)=coef
               endif
            else
! Standard non-boundary setting: not a boundary point.
               cij(icb*ipoint+2*(id-1)+i)=coef            
            endif
         enddo
! End of dimension iteration. cij coefficients now set.
      enddo
! Here is where the obsolete code for newbox went.

! Diagnostics sample -----------------------------------
      if(.false. .and. cij(ipoint*icb+ndims*2+1).ne.0)then
         cnon=0.
         do k=1,2*ndims
            cnon=cnon+abs(dob_cij(3*k-1,oihere))
     $           +abs(dob_cij(3*k,oihere))
         enddo
         if(cnon.eq.0.)then
            write(*,201)(indi(k),k=1,ndims),
     $           (cij(ipoint*icb+k),k=1,ndims*2+1) 
            write(*,202)(dob_cij(3*k-2,oihere),k=1,nobj_cij/3)
!            write(*,203)(dob_cij(3*k-1,oihere),k=1,nobj_cij/3)
!            write(*,204)(dob_cij(3*k,oihere),k=1,nobj_cij/3)
         elseif(.false. .and. mod(oihere,100).eq.0 )then
            write(*,201)(indi(k),k=1,ndims),
     $           (cij(ipoint*icb+k),k=1,ndims*2+1) 
 201        format('cij(',i2,',',i2,',',i2,')=',12f8.1)
            write(*,202)(dob_cij(3*k-2,oihere),k=1,nobj_cij/3)
 202        format('fract  (3k-2)=',14f8.3)
            write(*,203)(dob_cij(3*k-1,oihere),k=1,nobj_cij/3)
 203        format('b/a    (3k-1)=',14f8.2)
            write(*,204)(dob_cij(3*k,oihere),k=1,nobj_cij/3)
 204        format('c/a    (3k  )=',14f8.2)
         endif
      endif
! Skip extra data created above, if any.
      oi_cij=oi_cij+ichain
!----------------------------------------
!     Return increment of 1
      inc=1
      end
!********************************************************************
!********************************************************************
      subroutine ddn_cij(ip,dden,dnum)
! Routine to do the adjustment to dden and dnum for this point (ip)
      include 'ndimsdecl.f'
      include 'objcom.f'
      dden=dden+dob_cij(idgs_cij,ip)
      dnum=dnum+dob_cij(ibdy_cij,ip)
      end
!********************************************************************
! Initialize a specific object. 1s for frac, 0 for diag,potterm.
! Reverse pointer.
      subroutine objinit(dob,idob,ipoint)
      include 'ndimsdecl.f'
      include 'objcom.f'
      real dob(nobj_cij)
      integer idob(nobj_cij)
      do j=1,2*ndims
         dob(3*j-2)=1.
         dob(3*j-1)=0.
         dob(3*j)=0.
      enddo
      dob(idgs_cij)=0.
      dob(ibdy_cij)=0.
      idob(iflag_cij)=0
! Set the reverse pointer to the u/c arrays:
      idob(ipoint_cij)=ipoint
! Zero the chained pointer.
      idob(iextra_cij)=0
      end
!******************************************************************
      subroutine objstart(cijp,istart,ipoint)
      real cijp
      include 'ndimsdecl.f'
      include 'objcom.f'
      logical lfirst/.true./
      save lfirst
! Initialization to save block-data cost.
      if(lfirst)then
         oi_cij=0
         lfirst=.false.
      endif

! Start object data for this point if not already started.
      if(cijp.eq.0)then
         oi_cij=oi_cij+1
         if(oi_cij.gt.Lobjmax) then
            write(*,*)'cijobj: oi_cij overflow',oi_cij,
     $           ' Increase Lobjmax in objcom.'
            stop 
         endif
! Initialize the object data: 1s for frac, 0 for b/a,c/a,diag,...
         call objinit(dob_cij(1,oi_cij),idob_cij(1,oi_cij),ipoint)
! Set the pointer
         cijp=oi_cij
         istart=istart+1
      endif
      end
!******************************************************************
!*******************************************************************
! Spagetti code for general number of dimensions.
      subroutine boxedge(ndims,ipm,indi,fn,npoints,idin)
! Iterate over the edges of a box-cell in ndims dimensions.
! For the box with corners at indices indi(ndims) and indi+ipm(ndims),
! return the fractional intersections fn(ndims), and number of points.
! Diagnostics are controlled by idin. 
! The number of dimensions
      integer ndims
! The coordinates of this point
      integer indi(ndims)
! The directions (+-1) along each dimension.
      integer ipm(ndims)
! The fractions for each dimension.
      real fn(ndims)
      real conditions(3)
!-----------------------------------------------
! Iteration parameters. Leave these alone.
      integer mdims,mpoints
      parameter (mdims=3,mpoints=20)
! Array of counters at the different levels.
      integer jd(mdims)
! Array of coefficients in the different dimensions.
      integer icp(mdims)
!------------------------------------------------
! Local Index array
      integer indl(mdims)
! Storage/workspace for intersection points
      real xf(mpoints,mdims)
!      real af(mdims,mdims),
      real afs(mpoints,mdims)
      real bs(mpoints),V(mdims,mdims),W(mdims)
! Number of points found
      integer npoints
! Unfortunately ovlen is not available here. So as a kludge, using ndims.
      integer ijbin(ndims)

      idiag=idin
      npoints=0

! Iteration parameter and other initialization.
      do i=1,ndims
         icp(i)=0
         jd(i)=0
         fn(i)=1.
      enddo
! The edges of a box can be considered to be arrived at by
! First choose 1 direction of ndims to increment. 
! Each is an edge from  000.... We call these the level 1 edges.
! Then from each of those points (e.g. 000100) choose a second increment
! from ndims-1. These are the edges at level 2, to points at level 2
! Then from each of the level-2 points of the form (010100) (001010)...
! choose a third increment from ndims-2, giving edges of level 3. 
! Carry on like this till we reach level ndims,
! which is the opposite corner of the box.
! When selecting a starting point we exclude multiple paths.
! Let icp(mdims) contain 0 or 1 corresponding to the coordinates above.
! iw is the level of edge we are constructing
! il is the level we are iterating, il<=iw.
! i is the dimension within that level.

      ibreak=0
      iw=1
! The level we are iterating 
 100  il=1
 101  i=0
! Do over dimensions i in this level il
! Skipping dimensions that are already chosen.
 102  i=i+1
      if(i.gt.ndims)goto 104
      if(icp(i).eq.1)goto 102
      jd(il)=i
      if(il.eq.iw)then
! If we are at desired level: MAIN ACTION. 
! ---------------------------- Diagnostics ------------------
         if(idiag.gt.5)then
            icp(i)=9
            write(*,*)il,(icp(j),j=1,ndims),' jd=',(jd(j),j=1,ndims)
            icp(i)=0
         endif
! ------------------------- End Diagnostics ------------------
! Commands iterated at each active level.
! ===================================================
! Local indices and fractions of this edge start.
         do j=1,ndims
            indl(j)=indi(j)+ipm(j)*icp(j)
         enddo
! Look for intersection along this edge.
!      if(idiag.ge.5)write(*,'(a,i2,$)')'Calling potlsect '
         call potlsect(i,ipm(i),ndims,indl,fraction,conditions,dpm
     $        ,iobjno,ijbin(1))
!      if(idiag.ge.5)write(*,'(a,$)')'Returned'
         if(fraction.ne.1. .and. npoints.lt.mpoints)then
!            idiag=idiag+1
            npoints=npoints+1
            do j=1,ndims
               xf(npoints,j)=icp(j)
            enddo
            xf(npoints,i)=fraction
!            if(idiag.gt.0)write(*,*)'Intersection fraction=',fraction
!     $           ,i,il,iw,(indl(j),j=1,ndims),' x='
!     $           ,(xf(npoints,j),j=1,ndims)
         endif
! =================================================
! End of active level iteration.
         if(i.lt.ndims)goto 102
      else
! Else iterate at the next level up with this dimension excluded:
         icp(i)=1
         il=il+1
!         write(*,*)'Incremented to',il,i
! For levels below iw, don't allow multiple routes. Use higher i's.
         if(il.lt.iw)goto 102
         goto 101
      endif
 104  continue
!      write(*,*)'Finished level',il
      if(il.eq.1)then
!------------------------------------------
! We've finished an entire iw working-level
! Commands executed at the end of the active level:
! =================================================
         if(npoints.ge.ndims)then
            ibreak=ibreak+1
            do j=1,npoints
               bs(j)=1.
               do k=1,ndims
                  afs(j,k)=xf(j,k)
               enddo
            enddo
! return solution of af.as=bs=1, which is fn=1/fractions.
            call SVDsol(fn,ndims,bs,npoints,afs,V,W,mpoints,mdims)
            if(idiag.gt.0) then
!               write(*,'(12f6.3)')((xf(k1,k2),k2=1,ndims),k1=1,npoints)
!               write(*,*)
!     $              npoints,' svdsol ',
!     $              ((afs(j,k1),j=1,npoints),k1=1,ndims)
!               write(*,*)'1/fn=',(1./fn(k1),k1=1,ndims)
            endif
         endif
! =================================================
! End of commands executed
! If a break has been set by the routine, exit. 
         if(ibreak.gt.0)goto 200
         if(iw.lt.ndims)then
            iw=iw+1
            goto 100
         else
! We've finished all iw levels.
            goto 200
         endif
      else
! Working above il=1; decrement back to prior level
         il=il-1
         i=jd(il)
         icp(i)=0
!         write(*,*)'decremented to',il,i
         if(i.ge.ndims)goto 104
         goto 102
      endif
! Finished.
 200  continue
      if(idiag.gt.5)write(*,*)'Leaving boxedge',npoints

      end
!**********************************************************************
!************************************************************************
      subroutine cijedge(inc,ipoint,indi,mdims,iLs,iused,cij)
! Edge routine which sets the iregion for all boundary points.
! mdims argument is unused here because we use ndimsdecl to give ndims.
      integer ipoint,inc
      include 'ndimsdecl.f'
      integer indi(ndims),iused(ndims)
      real cij(*)
      include 'objcom.f'
      include 'partcom.f'
! Silence unused warnings
      icb=iLs

      icb=2*ndims+1
      icij=icb*ipoint+icb
! Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
! Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
! Initialize the cij because cijroutine now does not do so for boundary.
      cij(icij)=0
! Tests to see if we do the pointer setting (revert by uncommenting):
!      goto 102
      do id=1,ndims
         if(ipartperiod(id)/64-(ipartperiod(id)/128)*2.ne.1
     $        .and.ipartperiod(id).ne.4)then
            if(indi(id).eq.0)goto 102
         endif
         if(ipartperiod(id)/128-(ipartperiod(id)/256)*2.ne.1
     $        .and.ipartperiod(id).ne.4)then
            if(indi(id).eq.iused(id)-1)goto 102
         endif
      enddo
      goto 103
 102  continue
! Start object data for this point. 
      call objstart(cij(icij),ist,ipoint)
      idob_cij(iregion_cij,int(cij(icij)))=-1
 103   continue
! Calculate the increment:
      do n=ndims,2,-1
         if(indi(n).eq.0)then
! On boundary face 0 of dimension n>1. Break.
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
! On boundary face iused(n) of dimension n>1. Break.
            goto 101
         endif
      enddo
! This is where the boundary setting is done for n=1
      if(indi(n).eq.0)then
         inc=(iused(1)-1)
      elseif(indi(n).eq.iused(n)-1)then
         inc=1
      else
         write(*,*)'CIJEDGE Error. We should not be here',
     $        n,ipoint,indi,inc
         stop
      endif
 101  continue
      end
!****************************************************************
!*********************************************************************
! Print a text graphic of slices through the regions 
      subroutine text3vgraph(mdims,iuds,ifull,volumes)
      integer iuds(mdims),ifull(mdims)
!      real cij(*)
      real volumes(*)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      character*40 form1
! Standard volume for uniform mesh:
         vs=.99999
         do id=1,ndims
            vs=vs*(xn(ixnp(id)+2)-xn(ixnp(id)+1))
         enddo
         vs=vs
! Text graphic of slice through volumes
         write(*,*)'Volumes percentages:'
         write(form1,'(''('',i3,''i4)'')')iuds(1)
         write(*,form1)((nint(100.*
     $        volumes(j-1+(iuds(2)/2-1)*ifull(1)
     $                +(k-1)*ifull(1)*ifull(2)+1)/vs),
     $        j=1,iuds(1)),k=1,iuds(3))
         end
!*********************************************************************
      subroutine text3rgraph(mdims,iuds,ifull,cij)
      integer iuds(mdims),ifull(mdims)
      real cij(*)
!      real volumes(*)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      character*40 form1
! Text graphic of slice through cij
!      write(*,*)'iregion:'
!         write(form1,'(''('',i3,''i1)'')')iuds(1)
!         write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),
!     $        j=1,iuds(1)),k=1,iuds(3))
! Slightly more flexible ascii form:
      if(ndims.ne.3.or.mdims.ne.3)then
         write(*,*)'****text3rgraph expects 3-D. Trying anyway.'
      endif
      write(form1,'(''('',i3,''a1)'')')iuds(1)
      write(*,form1)
     $        ((char(min(126,48+ireg3(j,iuds(2)/2,k,ifull,cij))),
     $        j=1,iuds(1)),k=1,iuds(3))
      end
!*******************************************************************
      function ireg3(i,j,k,ifull,cij)
! Return the region value for point i,j,k specified in 3d.
      integer i,j,k
      include 'ndimsdecl.f'
      include 'objcom.f'
      integer ifull(3)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))

      ipoint=int(cij(ndims*2+1,i,j,k))
      if(ipoint.ne.0)then
         ireg3=idob_cij(iregion_cij,ipoint)
      else
         ireg3=99
      endif
      end
!*********************************************************************
      subroutine text3graphs(ndims,iuds,ifull,cij,volumes)
      integer iuds(ndims),ifull(ndims)
      real cij(*),volumes(*)
      write(*,*)'Edge volumes sample at iuds(2)/2, iuds(3)/2'
      write(*,*)(((volumes(j-1+(iuds(2)/2-1)*ifull(1)
     $                +(k-1)*ifull(1)*ifull(2)+1)),
     $        j=1,10),k=iuds(3)/2,iuds(3)/2)
      call text3vgraph(ndims,iuds,ifull,volumes)
      call text3rgraph(ndims,iuds,ifull,cij)
      end
! Do least squares fit using svd, based on numerical recipes routines.
!********************************************************************
! Given matrix U(nd/mp,ma) solve in a least squares sense the
! equation U.a = b, to find a(ma). (b is b(nd))
! U V W are returned as the SVD of the original U. (U is changed!)
! np>=ma is the leading storage dimension of V (and trailing of U).
!********************************************************************
      SUBROUTINE SVDsol(A,MA,b,nd,U,V,W,MP,NP)
      PARAMETER(TOL=1.E-5)
      real A(MA),U(MP,NP),V(NP,NP),W(NP),B(nd)
      CALL SVDCMP(U,ND,MA,MP,NP,W,V)
! Adjust singular values.
      WMAX=0.
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.
14    CONTINUE
      CALL SVBKSB(U,W,V,ND,MA,MP,NP,B,A)
      END
!********************************************************************
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
      PARAMETER (NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END
!********************************************************************
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
! Silence warnings.
      l=0
      nm=0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            RV1(I)=C*RV1(I)
            IF ((ABS(F)+ANORM).EQ.ANORM) GO TO 2
            G=W(I)
            H=SQRT(F*F+G*G)
            W(I)=H
            H=1.0/H
            C= (G*H)
            S=-(F*H)
            DO 42 J=1,M
              Y=A(J,NM)
              Z=A(J,I)
              A(J,NM)=(Y*C)+(Z*S)
              A(J,I)=-(Y*S)+(Z*C)
42          CONTINUE
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30)
     $         write(*,*) 'SVDCMP No convergence in 30 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 JJ=1,N
              X=V(JJ,J)
              Z=V(JJ,I)
              V(JJ,J)= (X*C)+(Z*S)
              V(JJ,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 JJ=1,M
              Y=A(JJ,J)
              Z=A(JJ,I)
              A(JJ,J)= (Y*C)+(Z*S)
              A(JJ,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
!********************************************************************
      real function phiofcount(count,area)
! Turn the particle count for the current timestep over an area 
! into a potential at which current density is zero.
      include 'ndimsdecl.f'
!ontains eoverms
      include 'plascom.f'
!ontains dt and rhoinf:
      include 'partcom.f'
! For preventing logarithm infinities when count is zero.
      parameter (small=1.e-2)
      
      if(count.lt.0.)then
         write(*,*)'Negative count in phiofcount!',count
         phiofcount=-2.
         return
      endif
      flogfac=0.5*alog(2.*3.1415926*eoverm/1836.)
      fluxdensity=(count+small)/(area*rhoinf*dt)
      phiofcount=alog(fluxdensity)+flogfac

! Trap for singularities etc.
      if(.not.(abs(phiofcount).lt.1.e12))then
         write(*,*)'phiofcount error',count,area,phiofcount
         stop
      endif
!      write(*,*)area,rhoinf,dt,count,fluxdensity,phiofcount
      end
!****************************************************************
! Routine for directly updating cij, which needs no iteration.
! The geometry: fractions etc, are assumed not to have changed.
! It just treats the auxiliary data that is present, rather than
! looking at every mesh node
      subroutine cijdirect(debyelen,error)
! If there is an object crossing next to the point at indi(ndims) whose
! pointer is ipoint, in the dimension id (plus or minus), this routine
! adjusts the cij values for situations where the cij are variables
! changing from step to step.  At present is is assumed only C changes
! and it is updated in a smoothed manner toward the potential that
! corresponds to floating.
! Error indicator
      real error,debyelen
!--------------------------------------------------
! Object-data storage.
      include 'ndimsdecl.f'
      include 'objcom.f'
! Object information
      include '3dcom.f'
!-----------------------------------------------------
      real tiny
      parameter (tiny=1.e-15)
      integer oisor,oiextra
! The total areas and flux counts of all flux-tracked objects
      real totarea(nf_obj),totflux(nf_obj)
! Probably this ought to be set up in a common. But for now:
      integer iavemax
      integer ijbin(ovlen)
      data iavemax/50/      

! Calculate the current step's totals for use in floating cases.
      do ifobj=1,mf_obj
         i2type=int(obj_geom(otype,nf_geommap(ifobj)))/256
         if(i2type.eq.5)then
            call objfluxtotal(ifobj,totflux(ifobj),totarea(ifobj))
            if(totflux(ifobj).lt.0.)then
               write(*,*)'Negative totflux',ifobj,totflux(ifobj)
     $              ,totarea(ifobj)
            endif
         endif
      enddo
!      write(*,*)'Flux and area totals',(totflux(k),totarea(k),k=1
!     $     ,mf_obj)

! Prevent divide by zero issues with debyelen
      debyehere=debyelen
      if(debyehere.lt.tiny)debyehere=tiny

! Instead of the above iteration over mesh, we simply iterate over 
! the auxiliary data. That is, oisor. oi_cij is the maximum number
! we've reached.
      oisor=1
      do ioi=1,oi_cij
      oiextra=idob_cij(iextra_cij,oisor)
      if(oiextra.eq.0)then
! Do nothing but advance to next.
         oisor=oisor+1
      else
! Byte 2 of the type: If it's 4 insulating, 5 floating.
         iobj=idob_cij(iinter_cij,oisor)
         i2type=int(obj_geom(otype,iobj))/256
!         write(*,*)'i2type=',i2type,iobj
!         ichain=1
         dibdy=dob_cij(ibdy_cij,oisor)
         dob_cij(ibdy_cij,oisor)=0.
! Iterate over dimensions.
         do id=1,ndims
! For each direction in this dimension,
            do i=1,2
               ipm=1-2*(i-1)
               im=mod(i,2)+1
               ioad=ndata_cij*(2*(id-1)+(i-1))+1
               if(dob_cij(ioad,oisor).lt.1
     $              .and.dob_cij(ioad,oisor).ge.0.)then
! We intersected an object in this direction. Adjust Cij and B_y
! Get back coef
                  iobj=idob_cij(ioad,oiextra)
                  coefoa=dob_cij(ioad+2,oiextra)
! These must get the information from somewhere.
! Assume that a and b are unchanged by the variability:
                  a=obj_geom(oabc,iobj)
                  b=obj_geom(oabc+1,iobj)/debyehere
! c is the thing that depends on flux, so it's going to be different.
                  ifobj=nf_map(iobj)
!-------------
                  if(i2type.eq.4)then
!                     write(*,*)'Insulating',ifobj
! Address the flux data
                     ijbin(1)=idob_cij(ioad+1,oiextra)
! Pull the area of this facet into here
                     iaddress=ijbin(1)+nf_address(nf_flux,ifobj,nf_pa)
                     area=ff_data(iaddress)
! Pull the flux for this timestep and this facet.
                     flux=ff_data(ijbin(1)
     $                    +nf_address(nf_flux,ifobj,nf_step))
! Calculate new potential
                     if(flux.lt.0.)then
                        write(*,*)'phiofcount flux negative',flux
     $                       ,ijbin(1),area,ifobj
                        flux=0.
                     endif
                     cnew=-a*phiofcount(flux,area)
                  elseif(i2type.eq.5)then
! Floating. Use totals
!                     write(*,*)'Floating',ifobj
                     cnew=-a*phiofcount(totflux(ifobj),totarea(ifobj))
                  else
! This should not happen.
                     cnew=0.
                  endif
! Smooth over nave steps. 
                  nave=min(nf_step,iavemax)
                  cold=a*dob_cij(ioad+2,oisor)
                  c=(cold*(nave-1)+cnew)/nave
!         if(oisor.lt.20)write(*,*)'cijupdate',nf_step,iobj,ijbin(1)
!     $                 ,cold,cnew,c
! or test that it gives the same answer from the old value.
!                  c=cold
!-------------

                  if(a.eq.0.) a=tiny
! Diagonal term (denominator) difference Cd-Cij stored
! This does not change
!                     dob_cij(idgs_cij,oisor)=dob_cij(idgs_cij
!     $                    ,oisor)+ coef
! Adjust potential sum (B_y or tau in new reference)
                  dob_cij(ibdy_cij,oisor)=dob_cij(ibdy_cij
     $                    ,oisor)- coefoa*c
! Now we need to update coa in the top data set. Not boa or fraction 
! since they haven't changed.                  
                  if(c .ne. a*dob_cij(ioad+2,oisor))then
!                     write(*,*)'Updating coa from',dob_cij(ioad+2,oisor)
!     $                    ,' to',c/a
                     dob_cij(ioad+2,oisor)=c/a
                  endif
               endif
            enddo
! End of dimension iteration. cij coefficients now set.
         enddo
! This was using the incorrect assumption that chained data follows.
!         oisor=oisor+ichain+1
         oisor=oisor+1
         error=0.
      endif
      if(oisor.gt.oi_cij)return
      enddo
      
      end
!********************************************************************
      subroutine objfluxtotal(ifobj,flux,area)
! Calculate and return the total flux and total area of a
! flux-collecting object ifobj.
      integer ifobj
      real flux,area
      include 'ndimsdecl.f'
      include '3dcom.f'
      area=0.
      flux=0.
! Starting addresses of area and flux of this step.
      iada=nf_address(nf_flux,ifobj,nf_pa)
      iadf=nf_address(nf_flux,ifobj,nf_step)
!      write(*,*)
!      write(*,*)nf_step,ifobj,', 3 steps forward from where we are:'
!     $     ,nf_address(1,ifobj,nf_step)-1
!      write(*,'(10f8.4)')(ff_data(k),k=nf_address(1,1,nf_step)
!     $     ,nf_address(1,ifobj,nf_step+3)-1)

      do i=1,nf_posno(nf_flux,ifobj)
         area=area+ff_data(iada+i-1)
         flux=flux+ff_data(iadf+i-1)
         if(ff_data(iadf+i-1).lt.0)then
            write(*,*)'Negative Flux in objfluxtotal',ifobj,nf_step,i
     $           ,ff_data(iadf+i-1),iada,iadf
            write(*,*)'Are you excluding the interior'
     $           ,' of floating objects?'
         endif
      enddo

      
      end
!**************************************************************
