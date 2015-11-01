c****************************************************************
c Routine for setting cij, which is used by the general mditerarg.
      subroutine cijroutine(inc,ipoint,indi,mdms,iLs,iused,cij,debyelen
     $     ,error)
c This routine sets the object pointer in cij if there is an object
c crossing next to the the point at indi(ndims) in the 
c dimension id (plus or minus), adjusts the cij values,
c and sets the object data into the place pointed to in obj. 
c mdms is dummy here because ndimsdecl is used.
      include 'ndimsdecl.f'
c Effective index in dimension, c-style (zero based)
c      integer indm1(ndims)
c Shifted to correct for face omission.
      integer indi(ndims),iused(ndims)
c Error indicator
      real error
c Not used in this routine
      real cij(*)
      integer iLs(ndims+1)
c--------------------------------------------------
c Object-data storage.
      include 'objcom.f'
c Object information
      include '3dcom.f'
c-----------------------------------------------------
      real dpm(2),fraction(2),dplus(2),deff(2)
      real conditions(3,2)
      real tiny
      parameter (tiny=1.e-15)
      real fn(ndims)
      integer ipt2,oihere
c      integer ipa(ndims)
c      logical newbox
      integer istart,ireport
      integer ijbin(ovlen)
      integer ipb(ndims)
      data ipb/ndims*-1/
      data istart/0/ireport/0/

c Default no chain.
      ichain=0
c ipoint here is the offset (zero-based pointer)
c The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
      icij=icb*(ipoint+1)
c Object pointer defaults to zero, unset initially.
c We are setting only trailing adjacent pointers in the loop below.
c Therefore we know that this cij(icij) has not yet been set.
      if(cij(icij).ne.0)then
c When we are sure this is working, we can simplify to cij(icij)=0.
         if(ireport.eq.0)write(*,*)'cij(icij) needed initialization'
     $        ,icij,cij(icij)
         ireport=ireport+1
         cij(icij)=0
      endif
c Prevent divide by zero issues with debyelen
      debyehere=debyelen
      if(debyehere.lt.tiny)debyehere=tiny
c New boxedge management section (newbox)
      idiag=0
c Look at the trailing box relative to point (since ipb is -1).
c Report any intersections as npoints. 
      call boxedge(ndims,ipb,indi,fn,npoints,idiag)
c There was doubt about what value of npoints to require, >0 is tight.
      if(npoints.gt.0)then
c Do over all the corners of this box
         do ic=0,2**ndims-1
            ipt2=ipoint
            do ib=0,ndims-1
c Avoiding the lower edge regions, which are boundaries.
c Upper boundaries are not set because they aren't trailing.
               if(indi(ib+1).le.1)goto 31
               if(btest(ic,ib))ipt2=ipt2-iLs(ib+1)
            enddo
            icij2=icb*(ipt2+1)
c Conditionally start the object: only if it does not already exist.
            call objstart(cij(icij2),istart,ipt2)
 31        continue
         enddo
      else
         cij(icij)=0.
      endif
      
c It is impossible now to assume this object is oi_cij because we might
c have started an addtional obj. So don't
      if(cij(icij).eq.0.)then
         oihere=oi_cij
      else
         oihere=int(cij(icij))
      endif
c Iterate over dimensions.
      do id=1,ndims
         ifound=0
c Find the intersection, if any, of each mesh branch with the objects.
c Return the fraction of the mesh distance of the intersection, and 
c the value of a, b, c to be attached to it. Fraction=1. is the 
c default case, where no intersection occurs. However, we need to know
c the fractions for opposite directions before we can calculate the
c coefficients; so we have to loop twice over the opposite directions.
c For each direction in this dimension,
         do i=1,2
            ioad=ndata_cij*(2*(id-1)+(i-1))+1
            ipm=1-2*(i-1)
c Determine whether this is a boundary point: adjacent a fraction ne 1.
c This potlsect call might not need ijbin.
            call potlsect(id,ipm,ndims,indi,
     $           fraction(i),conditions(1,i),dpm(i),iobjno,ijbin(1))
            if(fraction(i).lt.1. .and. fraction(i).ge.0.)then
               ifound=ifound+1
c Start object data for this point if not already started.
               call objstart(cij(icij),istart,ipoint)
c---------------
               oihere=int(cij(icij))
c If this object's a,b,c depends on flux, start extra info.
               if(int(obj_geom(otype,iobjno))/256.eq.4 .or.
     $              int(obj_geom(otype,iobjno))/256.eq.5)then
                  if(idob_cij(iextra_cij,oihere).eq.0)then
c Start chained extra data. This really is in position oi_cij+ichain.
                     ichain=1
c                     write(*,*)'Chaining',ipoint,(indi(kw),kw=1,ndims)
c     $                    ,id,i,iobjno,ijbin(1),cij(icij)
c First time, initialize the extra pointer
                     idob_cij(iextra_cij,oihere)=oi_cij+ichain
c and initialize the data for all directions
                     do ie=1,nobj_cij
                        dob_cij(ie,oi_cij+ichain)=0.
                     enddo
                  endif
c Set the data for this direction.          
                  idob_cij(ioad,oi_cij+ichain)=iobjno
                  idob_cij(ioad+1,oi_cij+ichain)=ijbin(1)
c ioad+2 is set below from coef/a                 
               endif
c---------------
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
               dob_cij(ioad,oihere)=max(fraction(i),tiny)
               if(a.eq.0.) a=tiny
               dob_cij(ioad+1,oihere)=b/a*debyehere
               dob_cij(ioad+2,oihere)=c/a
               idob_cij(iinter_cij,oihere)=iobjno
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
               ioad=ndata_cij*(2*(id-1)+(i-1))+1
               if(dob_cij(ioad,oihere).lt.1
     $              .and.dob_cij(ioad,oihere).ge.0.)then
c We intersected an object in this direction. Adjust Cij and B_y
                  a=conditions(1,i)
                  b=conditions(2,i)/debyehere
                  c=conditions(3,i)
                  if(a.eq.0.) a=tiny
                  if(b.ge.0.)then
c Active side.
                     cij(icb*ipoint+2*(id-1)+i)=0.
c Diagonal term (denominator) difference Cd-Cij stored
                     dob_cij(idgs_cij,oihere)=dob_cij(idgs_cij,oihere)
     $                 + coef
c Adjust potential sum (B_y)
                     dob_cij(ibdy_cij,oihere)=dob_cij(ibdy_cij,oihere)
     $                    - coef*c/a
                     if(ichain.gt.0)dob_cij(ioad+2,oihere+ichain)=coef/a
                  else
c Inactive side. Continuity.
                     dxp2=(1.-fraction(i))*dpm(i)/debyehere
                     apb=a*dxp2+b
                     boapb=b/apb
                     cij(icb*ipoint+2*(id-1)+i)=coef*boapb
c Diagonal term (denominator) difference Cd-Cij stored.
                     dob_cij(idgs_cij,oihere)=dob_cij(idgs_cij,oihere)
     $                 + coef*(1.-boapb)
c Adjust potential sum (B_y)
                     dob_cij(ibdy_cij,oihere)=dob_cij(ibdy_cij,oihere)
     $                    - coef*c*dxp2/apb
                     if(ichain.gt.0)dob_cij(ioad+2,oi_cij+ichain)=coef
     $                    *dxp2/apb
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
c Here is where the obsolete code for newbox went.

c Diagnostics sample -----------------------------------
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
c            write(*,203)(dob_cij(3*k-1,oihere),k=1,nobj_cij/3)
c            write(*,204)(dob_cij(3*k,oihere),k=1,nobj_cij/3)
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
c Skip extra data created above, if any.
      oi_cij=oi_cij+ichain
c----------------------------------------
c     Return increment of 1
      inc=1
      end
c********************************************************************
c********************************************************************
      subroutine ddn_cij(ip,dden,dnum)
c Routine to do the adjustment to dden and dnum for this point (ip)
      include 'ndimsdecl.f'
      include 'objcom.f'
      dden=dden+dob_cij(idgs_cij,ip)
      dnum=dnum+dob_cij(ibdy_cij,ip)
      end
c********************************************************************
c Initialize a specific object. 1s for frac, 0 for diag,potterm.
c Reverse pointer.
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
c Set the reverse pointer to the u/c arrays:
      idob(ipoint_cij)=ipoint
c Zero the chained pointer.
      idob(iextra_cij)=0
      end
c******************************************************************
      subroutine objstart(cijp,istart,ipoint)
      real cijp
      include 'ndimsdecl.f'
      include 'objcom.f'
      logical lfirst/.true./
      save lfirst
c Initialization to save block-data cost.
      if(lfirst)then
         oi_cij=0
         lfirst=.false.
      endif

c Start object data for this point if not already started.
      if(cijp.eq.0)then
         oi_cij=oi_cij+1
         if(oi_cij.gt.Lobjmax) then
            write(*,*)'cijobj: oi_cij overflow',oi_cij,
     $           ' Increase Lobjmax in objcom.'
            stop 
         endif
c Initialize the object data: 1s for frac, 0 for b/a,c/a,diag,...
         call objinit(dob_cij(1,oi_cij),idob_cij(1,oi_cij),ipoint)
c Set the pointer
         cijp=oi_cij
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
      parameter (mdims=3,mpoints=20)
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
c Unfortunately ovlen is not available here. So as a kludge, using ndims.
      integer ijbin(ndims)

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
     $        ,iobjno,ijbin(1))
c      if(idiag.ge.5)write(*,'(a,$)')'Returned'
         if(fraction.ne.1. .and. npoints.lt.mpoints)then
c            idiag=idiag+1
            npoints=npoints+1
            do j=1,ndims
               xf(npoints,j)=icp(j)
            enddo
            xf(npoints,i)=fraction
c            if(idiag.gt.0)write(*,*)'Intersection fraction=',fraction
c     $           ,i,il,iw,(indl(j),j=1,ndims),' x='
c     $           ,(xf(npoints,j),j=1,ndims)
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
c               write(*,'(12f6.3)')((xf(k1,k2),k2=1,ndims),k1=1,npoints)
c               write(*,*)
c     $              npoints,' svdsol ',
c     $              ((afs(j,k1),j=1,npoints),k1=1,ndims)
c               write(*,*)'1/fn=',(1./fn(k1),k1=1,ndims)
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
      subroutine cijedge(inc,ipoint,indi,mdims,iLs,iused,cij)
c Edge routine which sets the iregion for all boundary points.
c mdims argument is unused here because we use ndimsdecl to give ndims.
      integer ipoint,inc
      include 'ndimsdecl.f'
      integer indi(ndims),iused(ndims)
      real cij(*)
      include 'objcom.f'
      include 'partcom.f'
c Silence unused warnings
      icb=iLs

      icb=2*ndims+1
      icij=icb*ipoint+icb
c Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
c Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
c Initialize the cij because cijroutine now does not do so for boundary.
      cij(icij)=0
c Tests to see if we do the pointer setting (revert by uncommenting):
c      goto 102
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
c Start object data for this point. 
      call objstart(cij(icij),ist,ipoint)
      idob_cij(iregion_cij,int(cij(icij)))=-1
 103   continue
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
      end
c****************************************************************
c*********************************************************************
c Print a text graphic of slices through the regions 
      subroutine text3vgraph(mdims,iuds,ifull,volumes)
      integer iuds(mdims),ifull(mdims)
c      real cij(*)
      real volumes(*)
      include 'ndimsdecl.f'
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
      subroutine text3rgraph(mdims,iuds,ifull,cij)
      integer iuds(mdims),ifull(mdims)
      real cij(*)
c      real volumes(*)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      character*40 form1
c Text graphic of slice through cij
c      write(*,*)'iregion:'
c         write(form1,'(''('',i3,''i1)'')')iuds(1)
c         write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),
c     $        j=1,iuds(1)),k=1,iuds(3))
c Slightly more flexible ascii form:
      if(ndims.ne.3.or.mdims.ne.3)then
         write(*,*)'****text3rgraph expects 3-D. Trying anyway.'
      endif
      write(form1,'(''('',i3,''a1)'')')iuds(1)
      write(*,form1)
     $        ((char(min(126,48+ireg3(j,iuds(2)/2,k,ifull,cij))),
     $        j=1,iuds(1)),k=1,iuds(3))
      end
c*******************************************************************
      function ireg3(i,j,k,ifull,cij)
c Return the region value for point i,j,k specified in 3d.
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
c*********************************************************************
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
c********************************************************************
      real function phiofcount(count,area)
c Turn the particle count for the current timestep over an area 
c into a potential at which current density is zero.
      include 'ndimsdecl.f'
Contains eoverms
      include 'plascom.f'
Contains dt and rhoinf:
      include 'partcom.f'
c For preventing logarithm infinities when count is zero.
      parameter (small=1.e-2)
      
      if(count.lt.0.)then
         write(*,*)'Negative count in phiofcount!',count
         phiofcount=-2.
         return
      endif
      flogfac=0.5*alog(2.*3.1415926*eoverm/1836.)
      fluxdensity=(count+small)/(area*rhoinf*dt)
      phiofcount=alog(fluxdensity)+flogfac

c Trap for singularities etc.
      if(.not.(abs(phiofcount).lt.1.e12))then
         write(*,*)'phiofcount error',count,area,phiofcount
         stop
      endif
c      write(*,*)area,rhoinf,dt,count,fluxdensity,phiofcount
      end
c****************************************************************
c Routine for directly updating cij, which needs no iteration.
c The geometry: fractions etc, are assumed not to have changed.
c It just treats the auxiliary data that is present, rather than
c looking at every mesh node
      subroutine cijdirect(debyelen,error)
c If there is an object crossing next to the point at indi(ndims) whose
c pointer is ipoint, in the dimension id (plus or minus), this routine
c adjusts the cij values for situations where the cij are variables
c changing from step to step.  At present is is assumed only C changes
c and it is updated in a smoothed manner toward the potential that
c corresponds to floating.
c Error indicator
      real error,debyelen
c--------------------------------------------------
c Object-data storage.
      include 'ndimsdecl.f'
      include 'objcom.f'
c Object information
      include '3dcom.f'
c-----------------------------------------------------
      real tiny
      parameter (tiny=1.e-15)
      integer oisor,oiextra
c The total areas and flux counts of all flux-tracked objects
      real totarea(nf_obj),totflux(nf_obj)
c Probably this ought to be set up in a common. But for now:
      integer iavemax
      integer ijbin(ovlen)
      data iavemax/50/      

c Calculate the current step's totals for use in floating cases.
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
c      write(*,*)'Flux and area totals',(totflux(k),totarea(k),k=1
c     $     ,mf_obj)

c Prevent divide by zero issues with debyelen
      debyehere=debyelen
      if(debyehere.lt.tiny)debyehere=tiny

c Instead of the above iteration over mesh, we simply iterate over 
c the auxiliary data. That is, oisor. oi_cij is the maximum number
c we've reached.
      oisor=1
      do ioi=1,oi_cij
      oiextra=idob_cij(iextra_cij,oisor)
      if(oiextra.eq.0)then
c Do nothing but advance to next.
         oisor=oisor+1
      else
c Byte 2 of the type: If it's 4 insulating, 5 floating.
         iobj=idob_cij(iinter_cij,oisor)
         i2type=int(obj_geom(otype,iobj))/256
c         write(*,*)'i2type=',i2type,iobj
c         ichain=1
         dibdy=dob_cij(ibdy_cij,oisor)
         dob_cij(ibdy_cij,oisor)=0.
c Iterate over dimensions.
         do id=1,ndims
c For each direction in this dimension,
            do i=1,2
               ipm=1-2*(i-1)
               im=mod(i,2)+1
               ioad=ndata_cij*(2*(id-1)+(i-1))+1
               if(dob_cij(ioad,oisor).lt.1
     $              .and.dob_cij(ioad,oisor).ge.0.)then
c We intersected an object in this direction. Adjust Cij and B_y
c Get back coef
                  iobj=idob_cij(ioad,oiextra)
                  coefoa=dob_cij(ioad+2,oiextra)
c These must get the information from somewhere.
c Assume that a and b are unchanged by the variability:
                  a=obj_geom(oabc,iobj)
                  b=obj_geom(oabc+1,iobj)/debyehere
c c is the thing that depends on flux, so it's going to be different.
                  ifobj=nf_map(iobj)
c-------------
                  if(i2type.eq.4)then
c                     write(*,*)'Insulating',ifobj
c Address the flux data
                     ijbin(1)=idob_cij(ioad+1,oiextra)
c Pull the area of this facet into here
                     iaddress=ijbin(1)+nf_address(nf_flux,ifobj,nf_pa)
                     area=ff_data(iaddress)
c Pull the flux for this timestep and this facet.
                     flux=ff_data(ijbin(1)
     $                    +nf_address(nf_flux,ifobj,nf_step))
c Calculate new potential
                     if(flux.lt.0.)then
                        write(*,*)'phiofcount flux negative',flux
     $                       ,ijbin(1),area,ifobj
                        flux=0.
                     endif
                     cnew=-a*phiofcount(flux,area)
                  elseif(i2type.eq.5)then
c Floating. Use totals
c                     write(*,*)'Floating',ifobj
                     cnew=-a*phiofcount(totflux(ifobj),totarea(ifobj))
                  else
c This should not happen.
                     cnew=0.
                  endif
c Smooth over nave steps. 
                  nave=min(nf_step,iavemax)
                  cold=a*dob_cij(ioad+2,oisor)
                  c=(cold*(nave-1)+cnew)/nave
c         if(oisor.lt.20)write(*,*)'cijupdate',nf_step,iobj,ijbin(1)
c     $                 ,cold,cnew,c
c or test that it gives the same answer from the old value.
c                  c=cold
c-------------

                  if(a.eq.0.) a=tiny
c Diagonal term (denominator) difference Cd-Cij stored
c This does not change
c                     dob_cij(idgs_cij,oisor)=dob_cij(idgs_cij
c     $                    ,oisor)+ coef
c Adjust potential sum (B_y or tau in new reference)
                  dob_cij(ibdy_cij,oisor)=dob_cij(ibdy_cij
     $                    ,oisor)- coefoa*c
c Now we need to update coa in the top data set. Not boa or fraction 
c since they haven't changed.                  
                  if(c .ne. a*dob_cij(ioad+2,oisor))then
c                     write(*,*)'Updating coa from',dob_cij(ioad+2,oisor)
c     $                    ,' to',c/a
                     dob_cij(ioad+2,oisor)=c/a
                  endif
               endif
            enddo
c End of dimension iteration. cij coefficients now set.
         enddo
c This was using the incorrect assumption that chained data follows.
c         oisor=oisor+ichain+1
         oisor=oisor+1
         error=0.
      endif
      if(oisor.gt.oi_cij)return
      enddo
      
      end
c********************************************************************
      subroutine objfluxtotal(ifobj,flux,area)
c Calculate and return the total flux and total area of a
c flux-collecting object ifobj.
      integer ifobj
      real flux,area
      include 'ndimsdecl.f'
      include '3dcom.f'
      area=0.
      flux=0.
c Starting addresses of area and flux of this step.
      iada=nf_address(nf_flux,ifobj,nf_pa)
      iadf=nf_address(nf_flux,ifobj,nf_step)
c      write(*,*)
c      write(*,*)nf_step,ifobj,', 3 steps forward from where we are:'
c     $     ,nf_address(1,ifobj,nf_step)-1
c      write(*,'(10f8.4)')(ff_data(k),k=nf_address(1,1,nf_step)
c     $     ,nf_address(1,ifobj,nf_step+3)-1)

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
c**************************************************************
