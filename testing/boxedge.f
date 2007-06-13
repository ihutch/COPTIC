c Iterate over the edges of a box-cell in ndims dimensions.
c      subroutine boxedge(ndims,ipm,iLs,indi,boxroutine,arguments)
c The number of dimensions
c      integer ndims
c Testing main
      parameter (ndims=3)
c
c The coordinates of this point
c      integer indi(ndims)
c The directions (+-1) along each dimension.
      integer ipm(ndims)
c      real arguments(*)
c The fractions for each dimension.
      real fn(ndims)
      real conditions(3)

c-----------------------------------------------
c Iteration parameters. Leave these alone.
      integer mdims
      parameter (mdims=10)
c Array of counters at the different levels.
      integer jd(mdims)
c Array of coefficients in the different dimensions.
      integer icp(mdims)
c------------------------------------------------
c Local Index array
      integer indl(mdims)
c Storage/workspace for intersection points
      real xf(mdims,mdims)
c Number of points found
      integer npoints

      npoints=0

c Iteration parameter initialization.
      do i=1,ndims
         icp(i)=0
         jd(i)=0
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
         icp(i)=1
         write(*,*)il,(icp(j),j=1,ndims),' jd=',(jd(j),j=1,ndims)
         icp(i)=0
c ------------------------- End Diagnostics ------------------
c Don't mess around outside the following box.
c ------------------------------------------ 
c Local indices of this edge start.
         do j=1,ndims
            indl(j)=indi(j)+ipm(j)*icp(j)
         enddo
c Look for intersection along this edge.
         call potlsect(i,ipm(i),ndims,indl,fraction,conditions,dpm)
         if(fraction.ne.1.)then
            if(npoints.eq.ndims)then
c We already have ndims points. Decide whether to replace one.
            else
               npoints=npoints+1
            endif
         endif
c ------------------------------------------
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
c      write(*,*)'Finished level',il
c If a break has been set by the routine, exit. 
 104  if(ibreak.gt.0)goto 200
      if(il.eq.1)then
c We've finished an entire iw working-level
         if(iw.lt.ndims)then
            iw=iw+1
            goto 100
         else
c We've finished all iw levels.
            goto 200
         endif
      else
c Working above il=1; decrement back to prior level
 103     il=il-1
         i=jd(il)
         icp(i)=0
c         write(*,*)'decremented to',il,i
         if(i.ge.ndims)goto 104
         goto 102
      endif
c Finished.
 200  continue
      end
c**********************************************************************
      subroutine boxroutine(arguments)
c For this edge, find if there's an intersection. 
c If so, then if we have not yet found ndims intersections, add it to the
c list, return. If we have found ndims, calculate the intersection's position
c relative to the previously identified plane.
c If it is interior, replace some point (which?) with the new one.
c If we are progressing away from the most interior point, probably we
c should replace only points at the same level, never lower. So there
c have to be replaceable and non-replaceable points.

      end
