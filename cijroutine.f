c********************************************************************
c This initialization of oi_sor costs a big size hit on object.
c Now done in objstart.
c      block data objcomset
c      include 'objcom.f'
c      data oi_sor/0/
c      end
c****************************************************************
c Routine for setting cij, which is used by the general mditerate.
c      subroutine cijroutine(inc,ipoint,indm1,ndims,iused,cij)
      subroutine cijroutine(inc,ipoint,indi,ndims,iused,cij,debyelen
     $     ,error)
c This routine sets the object pointer in cij if there is an object
c crossing next to the the point at indi(ndims) in the 
c dimension id (plus or minus), adjusts the cij values,
c and sets the object data into the place pointed to in obj. 
      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
c      integer indm1(mdims)
c Shifted to correct for face omission.
      integer indi(mdims)
      integer ndims
c Error indicator
      real error
c Not used in this routine
      integer iused(ndims)
      real cij(*)
c--------------------------------------------------
c Object-data storage.
      include 'objcom.f'
c-----------------------------------------------------
      real dpm(2),fraction(2),dplus(2),deff(2)
      real conditions(3,2)
      real tiny
      parameter (tiny=1.e-15)
      integer ipa(mdims)
      real fn(mdims)
      integer istart
      data istart/0/

c iused is not used here. But silence the warning by spurious use:
      icb=iused(1)
c ipoint here is the offset (zero-based pointer)
c The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
c Object pointer defaults to zero, unset initially..
c      icij=icb*ipoint+2*ndims+1
      icij=icb*(ipoint+1)
      cij(icij)=0.
c Prevent divide by zero issues with debyelen
      debyehere=debyelen
      if(debyehere.lt.tiny)debyehere=tiny

c Iterate over dimensions.
      do id=1,ndims
         ifound=0
c Find the intersection, if any, of each mesh branch with the objects.
c Return the fraction of the mesh distance of the intersection, and 
c the value of a, b, c to be attached to it. Fraction=1. is the 
c default case, where no intersection occurs. However, we need to know
c the fractions for opposite directions before we can calculate the
c coefficients; so we have to loop twice over the opposite directions.
c not used         icb2=2*icb
c For each direction in this dimension,
         do i=1,2
            iobj=ndata_sor*(2*(id-1)+(i-1))+1
            ipm=1-2*(i-1)
c Determine whether this is a boundary point: adjacent a fraction ne 1.
c            write(*,'(a,$)')'potlsect'
            call potlsect(id,ipm,ndims,indi,
     $           fraction(i),conditions(1,i),dpm(i),iobjno)
            if(fraction(i).lt.1. .and. fraction(i).ge.0.)then
               ifound=ifound+1
c Start object data for this point if not already started.
c Here on 1 Sep 09 istart was ist. Which seemed incorrect.
               call objstart(cij(icij),istart,ipoint)
c Calculate dplus and deff for this direction.
c dplus becomes dminus for the other direction.
               a=conditions(1,i)
               b=conditions(2,i)/debyehere
               c=conditions(3,i)
               dxp1=fraction(i)*dpm(i)/debyehere
               if(a.eq.0.) a=tiny
               if(b.ge.0.)then
c Active
                  apb=a*dxp1+b
                  dplus(i)=dxp1*(apb+b)/apb
                  deff(i)=apb/a
               else
c Inactive
                  if(dxp1.eq.0.) dxp1=tiny
                  dxp2=(1.-fraction(i))*dpm(i)/debyehere
                  apb=a*dxp2+b
                  boapb=b/apb
c  Without \psi'', i.e. setting it to zero and boundary.
c                  dplus(i)=dxp1
c  With \psi'', i.e. extrapolating it uniformly over the boundary
                  dplus(i)=dxp1+boapb*dxp2**2/dxp1
                  deff(i)=dxp1
               endif
               if(deff(i).eq.0)deff(i)=tiny
c The object data consists of data enumerated as
c ndims*(2=forward/backward from point)*(3=fraction,b/a,c/a)
c + diagonal + potential terms.
c Prevent subsequent divide by zero danger.
               dob_sor(iobj,oi_sor)=max(fraction(i),tiny)
               if(a.eq.0.) a=tiny
               dob_sor(iobj+1,oi_sor)=b/a*debyehere
               dob_sor(iobj+2,oi_sor)=c/a
               idob_sor(iinter_sor,oi_sor)=iobjno
            else
c No intersection.
               dplus(i)=dpm(i)/debyehere
               deff(i)=dpm(i)/debyehere
            endif
         enddo
c Now the dplus and deff are set correctly for each direction.
c Calculate the coefficients for each direction, using the opposite
c dplus for the dminus.
         do i=1,2
            im=mod(i,2)+1
c coefficient Cd same for all cases
            coef=2./(deff(i)*(dplus(i)+dplus(im)))
            if(ifound.gt.0)then
c This is a boundary point.
               iobj=ndata_sor*(2*(id-1)+(i-1))+1
               if(dob_sor(iobj,oi_sor).lt.1
     $              .and.dob_sor(iobj,oi_sor).ge.0.)then
c We intersected an object in this direction. Adjust Cij and B_y
                  a=conditions(1,i)
                  b=conditions(2,i)/debyehere
                  c=conditions(3,i)
                  if(a.eq.0.) a=tiny
                  if(b.ge.0.)then
c Active side.
                     cij(icb*ipoint+2*(id-1)+i)=0.
c Diagonal term (denominator) difference Cd-Cij stored
                     dob_sor(idgs_sor,oi_sor)=dob_sor(idgs_sor,oi_sor)
     $                 + coef
c Adjust potential sum (B_y)
                     dob_sor(ibdy_sor,oi_sor)=dob_sor(ibdy_sor,oi_sor)
     $                    - coef*c/a
                  else
c Inactive side. Continuity.
                     dxp2=(1.-fraction(i))*dpm(i)/debyehere
                     apb=a*dxp2+b
                     boapb=b/apb
                     cij(icb*ipoint+2*(id-1)+i)=coef*boapb
c Diagonal term (denominator) difference Cd-Cij stored.
                     dob_sor(idgs_sor,oi_sor)=dob_sor(idgs_sor,oi_sor)
     $                 + coef*(1.-boapb)
c Adjust potential sum (B_y)
                     dob_sor(ibdy_sor,oi_sor)=dob_sor(ibdy_sor,oi_sor)
     $                    - coef*c*dxp2/apb
                  endif
               else
c We did not intersect an object in this direction.
                  cij(icb*ipoint+2*(id-1)+i)=coef
               endif
            else
c Standard non-boundary setting: not a boundary point.
               cij(icb*ipoint+2*(id-1)+i)=coef            
            endif
         enddo
c End of dimension iteration. cij coefficients now set.
      enddo

c----------------------------------------
c The following sets the object pointer in cij if there is an object
c crossing within the one of the 2**ndims boxes which have
c ijk...=indi(ndims) as a vertex, if it is not already set.
c A "box" is an ndims-dimensional cube with opposite vertices 
c indi(ndims) and indi(ndims)+ipa(ndims) with ipa(i)=+-1, defining
c the box of interest. If the dimension i is under consideration, then
c ipa(i) defines the direction under consideration. 
c
c It then calculates the extended intersection fraction for all the
c directions from indi(ndims) and enters it into the object data.
c The routine boxedge does the calculations.
c If a fraction is already set <1 (in a direction). Use that.  If not,
c the fraction in each direction represents the axis minimal crossing of
c bounding planes for that direction. Minimal means the closest to +1,
c in the order 1->\infty - \infty -> -0. Bounding planes are planes
c through ndims-lets of box-edge intersections. The boxes relevant to a
c direction are the 2**(ndims-1) that share this edge.
c If the returned fraction is 0<f<1, but the prior calculated was not,
c so that the b/a and c/a are both zero. Don't shift the fraction.
c If there are no intersections of a box, do no setting for that box (which
c may mean that some edges remain unset=1, but some might already be set).
c
c----------------------------------------
c Do the box examination and fraction update.
c Section for using the boxedge call. For each ndims-let of directions
c call boxedge with the appropriate ipmarray. 2**ndims in all.
      do j=1,2**ndims
         k=j-1
         do i=1,ndims
c This could be done with bit manipulation probably more efficiently:
c direction(j)=bit(j)of(i). 0->+1, 1->-1.
            ipa(i)=1-2*(k-2*(k/2))
            k=k/2
         enddo
c Now ipa is set. Call boxedge, returning the inverse of fractions
c in fn, and the number of intersections found in npoints.
         idiag=0
c         if(oi_sor.eq.3276)idiag=5
c         if(oi_sor.eq.288.or.oi_sor.eq.289)idiag=5
c         if(ipoint.gt.101210)idiag=5
         call boxedge(ndims,ipa,indi,fn,npoints,idiag)
c         if(idiag.ne.0.and.npoints.ne.0)then
c            write(*,*)oi_sor,npoints,' fn=',(fn(iw),iw=1,ndims)
c         endif
c
         if(npoints.ge.ndims)then
            ftot=0
            do i=1,ndims
               ftot=ftot+fn(i)
            enddo
c See if this plane actually cuts the cell. If not, don't add it.
c This drastically reduces the number of boxes counted.
c Also in 3-D it limits additions to 6-intersection cases.
c            if(ftot.gt.2.)then
c However, then it is impossible to assume, when dealing with a box
c containing a point in known region, that any box vertex with no pointer
c is in that region. This breaks getpotential fillinlin, so disable
c for now.
            if(.true.)then
c Conditionally start the object: only if it does not already exist.
               call objstart(cij(icij),istart,ipoint)
c Set the flag
               ifl=2**(j-1)
               idob_sor(iflag_sor,oi_sor)=idob_sor(iflag_sor,oi_sor)+ifl
c               idob_sor(iflag_sor,oi_sor)=
c     $              ibset(idob_sor(iflag_sor,oi_sor),j-1)
c-------- Diagnostics
c               itemp=istart
c               if(istart.ne.itemp)then 
c                  write(*,*)'Added start',istart,npoints
c     $              ,(indi(iw),iw=1,ndims)
c     $              ,' fractions=',(1./fn(iw),iw=1,ndims)
c     $              ,(fn(iw),iw=1,ndims)
c     $                 ,(ipa(iw),iw=1,ndims)
c                  write(*,*)ndims,ipa,(indi(iw),iw=1,ndims)
c     $                 ,(fn(iw),iw=1,ndims),npts
c Redo the intersect with diagnostics on:
c                  call boxedge(ndims,ipa,indi,fn,npts,1)
c                  write(*,*)'Diag:',npts,(fn(iw),iw=1,ndims)
c               endif
c 211        format(6i4,3f10.4,i4)
c            write(*,211)(indi(kk),ipa(kk),kk=1,ndims),
c     $           (1./fn(kk),kk=1,ndims),npoints
c----------Diagnostics end.
c Now use these ndims fractions to update the fractions already inserted,
c [not] if new ones are "smaller" (closer to zero on the +ve side),
c or if the fraction is 1, implying not set. 
               do i=1,ndims
                  iobj=ndata_sor*(2*(i-1)+(1-ipa(i))/2)+1
                  f0=dob_sor(iobj,oi_sor)
c Only if this is the first entry this direction 
                  if(f0.eq.1)then 
c                  if( f0.eq.1 
c     $                 .or. abs(f0*fn(i)-1.).gt.0.1
c     $                 .or. 1./(f0+sign(tiny,f0)).lt.fn(i)-1.e-3
c     $                 )then
                     f1=1./(sign(max(abs(fn(i)),tiny),fn(i)))
c Warn if any strange crossings found.
                     if(f1.lt.1. .and. f1.ge.0)then 
c This should never happen. If it does, it's a bug or mesh clash.
                        write(*,*)'Warning: Box Recut ',npoints
     $                       ,i,ipa(i),oi_sor
     $                       ,(indi(kk),kk=1,ndims),(1./fn(i))
     $                       ,' Adjust mesh!'
c     $                       ,' Try changing the grid'
c     $                       ,' dimensions a little.'
                        error=error+1
c     $                       ,(ipa(kk),kk=1,ndims)
c     $                       ,f0,f1,ftot
c     $                       ,(1./fn(kk),kk=1,ndims)
c     $                      ,idob_sor(iflag_sor,oi_sor)
c But set it not equal to 1, so we know it was set.
                        dob_sor(iobj,oi_sor)=1.001
c                        dob_sor(iobj,oi_sor)=f1
                     else
                        dob_sor(iobj,oi_sor)=f1
                     endif
                  endif
               enddo
            endif
         endif
      enddo


c Diagnostics sample -----------------------------------

      if(.false. .and. cij(ipoint*icb+ndims*2+1).ne.0)then
         cnon=0.
         do k=1,2*ndims
            cnon=cnon+abs(dob_sor(3*k-1,oi_sor))
     $           +abs(dob_sor(3*k,oi_sor))
         enddo
         if(cnon.eq.0.)then
            write(*,201)(indi(k),k=1,ndims),
     $           (cij(ipoint*icb+k),k=1,ndims*2+1) 
            write(*,202)(dob_sor(3*k-2,oi_sor),k=1,nobj_sor/3)
c            write(*,203)(dob_sor(3*k-1,oi_sor),k=1,nobj_sor/3)
c            write(*,204)(dob_sor(3*k,oi_sor),k=1,nobj_sor/3)
         elseif(.false. .and. mod(oi_sor,100).eq.0 )then
            write(*,201)(indi(k),k=1,ndims),
     $           (cij(ipoint*icb+k),k=1,ndims*2+1) 
 201        format('cij(',i2,',',i2,',',i2,')=',12f8.1)
            write(*,202)(dob_sor(3*k-2,oi_sor),k=1,nobj_sor/3)
 202        format('fract  (3k-2)=',14f8.3)
            write(*,203)(dob_sor(3*k-1,oi_sor),k=1,nobj_sor/3)
 203        format('b/a    (3k-1)=',14f8.2)
            write(*,204)(dob_sor(3*k,oi_sor),k=1,nobj_sor/3)
 204        format('c/a    (3k  )=',14f8.2)
         endif
      endif
c----------------------------------------
c     Return increment of 1
      inc=1
      end
c********************************************************************
c********************************************************************
      subroutine ddn_sor(ip,dden,dnum)
c Routine to do the adjustment to dden and dnum for this point (ip)
      include 'objcom.f'
      dden=dden+dob_sor(idgs_sor,ip)
      dnum=dnum+dob_sor(ibdy_sor,ip)
      end
c********************************************************************
c Initialize a specific object. 1s for frac, 0 for diag,potterm.
c Reverse pointer.
      subroutine objinit(dob,idob,ipoint)
      include 'objcom.f'
      real dob(nobj_sor)
      integer idob(nobj_sor)
      do j=1,2*ndims_sor
         dob(3*j-2)=1.
         dob(3*j-1)=0.
         dob(3*j)=0.
      enddo
      dob(idgs_sor)=0.
      dob(ibdy_sor)=0.
      idob(iflag_sor)=0
c Set the reverse pointer to the u/c arrays:
      idob(ipoint_sor)=ipoint
      end
c******************************************************************
      subroutine objstart(cijp,istart,ipoint)
      real cijp
      include 'objcom.f'
      logical lfirst/.true./
      save lfirst
c Initialization to save block-data cost.
      if(lfirst)then
         oi_sor=0
         lfirst=.false.
      endif

c Start object data for this point if not already started.
      if(cijp.eq.0)then
         oi_sor=oi_sor+1
         if(oi_sor.gt.lobjmax) then
            write(*,*)'cijobj: oi_sor overflow',oi_sor,
     $           ' Increase in objcom.'
            stop 
         endif
c Initialize the object data: 1s for frac, 0 for b/a,c/a,diag,...
         call objinit(dob_sor(1,oi_sor),idob_sor(1,oi_sor),ipoint)
c Set the pointer
         cijp=oi_sor
         istart=istart+1
      endif
      end
c******************************************************************
c*******************************************************************
c Spagetti code for general number of dimensions.
      subroutine boxedge(ndims,ipm,indi,fn,npoints,idin)
c Iterate over the edges of a box-cell in ndims dimensions.
c For the box with corners at indices indi(ndims) and indi+ipm(ndims),
c return the fractional intersections fn(ndims), and number of points.
c Diagnostics are controlled by idin. 
c The number of dimensions
      integer ndims
c The coordinates of this point
      integer indi(ndims)
c The directions (+-1) along each dimension.
      integer ipm(ndims)
c The fractions for each dimension.
      real fn(ndims)
      real conditions(3)
c-----------------------------------------------
c Iteration parameters. Leave these alone.
      integer mdims,mpoints
      parameter (mdims=10,mpoints=20)
c Array of counters at the different levels.
      integer jd(mdims)
c Array of coefficients in the different dimensions.
      integer icp(mdims)
c------------------------------------------------
c Local Index array
      integer indl(mdims)
c Storage/workspace for intersection points
      real xf(mpoints,mdims)
c      real af(mdims,mdims),
      real afs(mpoints,mdims)
      real bs(mpoints),V(mdims,mdims),W(mdims)
c Number of points found
      integer npoints

      idiag=idin
      npoints=0

c Iteration parameter and other initialization.
      do i=1,ndims
         icp(i)=0
         jd(i)=0
         fn(i)=1.
      enddo
c The edges of a box can be considered to be arrived at by
c First choose 1 direction of ndims to increment. 
c Each is an edge from  000.... We call these the level 1 edges.
c Then from each of those points (e.g. 000100) choose a second increment
c from ndims-1. These are the edges at level 2, to points at level 2
c Then from each of the level-2 points of the form (010100) (001010)...
c choose a third increment from ndims-2, giving edges of level 3. 
c Carry on like this till we reach level ndims,
c which is the opposite corner of the box.
c When selecting a starting point we exclude multiple paths.
c Let icp(mdims) contain 0 or 1 corresponding to the coordinates above.
c iw is the level of edge we are constructing
c il is the level we are iterating, il<=iw.
c i is the dimension within that level.

      ibreak=0
      iw=1
c The level we are iterating 
 100  il=1
 101  i=0
c Do over dimensions i in this level il
c Skipping dimensions that are already chosen.
 102  i=i+1
      if(i.gt.ndims)goto 104
      if(icp(i).eq.1)goto 102
      jd(il)=i
      if(il.eq.iw)then
c If we are at desired level: MAIN ACTION. 
c ---------------------------- Diagnostics ------------------
         if(idiag.gt.5)then
            icp(i)=9
            write(*,*)il,(icp(j),j=1,ndims),' jd=',(jd(j),j=1,ndims)
            icp(i)=0
         endif
c ------------------------- End Diagnostics ------------------
c Commands iterated at each active level.
c ===================================================
c Local indices and fractions of this edge start.
         do j=1,ndims
            indl(j)=indi(j)+ipm(j)*icp(j)
         enddo
c Look for intersection along this edge.
c      if(idiag.ge.5)write(*,'(a,i2,$)')'Calling potlsect '
         call potlsect(i,ipm(i),ndims,indl,fraction,conditions,dpm
     $        ,iobjno)
c      if(idiag.ge.5)write(*,'(a,$)')'Returned '
         if(fraction.ne.1. .and. npoints.lt.mpoints)then
c            idiag=idiag+1
            npoints=npoints+1
            do j=1,ndims
               xf(npoints,j)=icp(j)
            enddo
            xf(npoints,i)=fraction
            if(idiag.gt.0)write(*,*)'Intersection fraction=',fraction
     $           ,i,il,iw,(indl(j),j=1,ndims),' x='
     $           ,(xf(npoints,j),j=1,ndims)
         endif
c =================================================
c End of active level iteration.
         if(i.lt.ndims)goto 102
      else
c Else iterate at the next level up with this dimension excluded:
         icp(i)=1
         il=il+1
c         write(*,*)'Incremented to',il,i
c For levels below iw, don't allow multiple routes. Use higher i's.
         if(il.lt.iw)goto 102
         goto 101
      endif
 104  continue
c      write(*,*)'Finished level',il
      if(il.eq.1)then
c------------------------------------------
c We've finished an entire iw working-level
c Commands executed at the end of the active level:
c =================================================
         if(npoints.ge.ndims)then
            ibreak=ibreak+1
            do j=1,npoints
               bs(j)=1.
               do k=1,ndims
                  afs(j,k)=xf(j,k)
               enddo
            enddo
c return solution of af.as=bs=1, which is fn=1/fractions.
            call SVDsol(fn,ndims,bs,npoints,afs,V,W,mpoints,mdims)
            if(idiag.gt.0) then
               write(*,'(12f6.3)')((xf(k1,k2),k2=1,ndims),k1=1,npoints)
               write(*,*)
     $              npoints,' svdsol ',
     $              ((afs(j,k1),j=1,ndims),k1=1,npoints)
            endif
         endif
c =================================================
c End of commands executed
c If a break has been set by the routine, exit. 
         if(ibreak.gt.0)goto 200
         if(iw.lt.ndims)then
            iw=iw+1
            goto 100
         else
c We've finished all iw levels.
            goto 200
         endif
      else
c Working above il=1; decrement back to prior level
         il=il-1
         i=jd(il)
         icp(i)=0
c         write(*,*)'decremented to',il,i
         if(i.ge.ndims)goto 104
         goto 102
      endif
c Finished.
 200  continue
      if(idiag.gt.5)write(*,*)'Leaving boxedge',npoints

      end
c**********************************************************************
c************************************************************************
      subroutine cijedge(inc,ipoint,indi,ndims,iused,cij)
c Edge routine which sets the iregion for all boundary points.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real cij(*)
      parameter (mdims=10)
c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      integer iLs(mdims+1)
      common /iLscom/iLs

      include 'objcom.f'

      icb=2*ndims+1
      icij=icb*ipoint+icb
c Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
c Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
c Start object data for this point if not already started.
      call objstart(cij(icij),ist,ipoint)
      idob_sor(iregion_sor,int(cij(icij)))=-1
c Calculate the increment:
      do n=ndims,2,-1
         if(indi(n).eq.0)then
c On boundary face 0 of dimension n>1. Break.
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On boundary face iused(n) of dimension n>1. Break.
            goto 101
         endif
      enddo
c This is where the boundary setting is done for n=1
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
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
c****************************************************************
c*********************************************************************
c Print a text graphic of slices through the regions 
      subroutine text3vgraph(ndims,iuds,ifull,cij,volumes)
      integer iuds(ndims),ifull(ndims)
      real cij(*)
      real volumes(*)
      include 'meshcom.f'
      character*40 form1
c Standard volume for uniform mesh:
         vs=.99999
         do id=1,ndims
            vs=vs*(xn(ixnp(id)+2)-xn(ixnp(id)+1))
         enddo
         vs=vs
c Text graphic of slice through volumes
         write(*,*)'Volumes percentages:'
         write(form1,'(''('',i3,''i4)'')')iuds(1)
         write(*,form1)((nint(100.*
     $        volumes(j-1+(iuds(2)/2-1)*ifull(1)
     $                +(k-1)*ifull(1)*ifull(2)+1)/vs),
     $        j=1,iuds(1)),k=1,iuds(3))
         end
c*********************************************************************
      subroutine text3rgraph(ndims,iuds,ifull,cij,volumes)
      integer iuds(ndims),ifull(ndims)
      real cij(*)
      real volumes(*)
      include 'meshcom.f'
      character*40 form1
c Text graphic of slice through cij
c      write(*,*)'iregion:'
c         write(form1,'(''('',i3,''i1)'')')iuds(1)
c         write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),
c     $        j=1,iuds(1)),k=1,iuds(3))
c Slightly more flexible ascii form:
         write(form1,'(''('',i3,''a1)'')')iuds(1)
         write(*,form1)
     $        ((char(min(126,48+ireg3(j,iuds(2)/2,k,ifull,cij))),
     $        j=1,iuds(1)),k=1,iuds(3))
      end
c*********************************************************************
      subroutine text3graphs(ndims,iuds,ifull,cij,volumes)
      call text3vgraph(ndims,iuds,ifull,cij,volumes)
      call text3rgraph(ndims,iuds,ifull,cij,volumes)
      end
c Do least squares fit using svd, based on numerical recipes routines.
c********************************************************************
c Given matrix U(nd/mp,ma) solve in a least squares sense the
c equation U.a = b, to find a(ma). (b is b(nd))
c U V W are returned as the SVD of the original U. (U is changed!)
c np>=ma is the leading storage dimension of V (and trailing of U).
c********************************************************************
      SUBROUTINE SVDsol(A,MA,b,nd,U,V,W,MP,NP)
      PARAMETER(TOL=1.E-5)
      real A(MA),U(MP,NP),V(NP,NP),W(NP),B(nd)
      CALL SVDCMP(U,ND,MA,MP,NP,W,V)
c Adjust singular values.
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
c********************************************************************
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
c********************************************************************
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
c Silence warnings.
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
