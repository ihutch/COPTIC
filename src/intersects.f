!**********************************************************************
      function insideall(mdims,x)
! For an ndims-dimensional point x, return the integer insideall
! consisting of bits i=0-30 that are zero or one according to whether
! the point x is outside or inside object i.
      integer mdims
      real x(mdims)
! Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
!
      insideall=0
      do i=ngeomobj,1,-1
         insideall=2*insideall+inside_geom(ndims,x,i)
      enddo
!      if(insideall.ne.0) write(*,*)
!     $     'ngeomobj=',ngeomobj,' insideall=',insideall,x

      end
!****************************************************************
      function insidemask(mdims,x)
! For an ndims-dimensional point x, return the integer insidemask
! consisting of bits i=1-31 that are zero or one according to whether
! the point x is outside or inside object i. But bits are masked by
! the regionmask.
      integer mdims
      real x(mdims)
      intrinsic btest
! Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
!
      insidemask=0
      ibit=ngeomobj
      do i=ngeomobj,1,-1
         insidemask=2*insidemask
         if(btest(ifield_mask,i-1))
     $        insidemask=insidemask+inside_geom(ndims,x,i)
      enddo

      end
!*****************************************************************
! Return the masked iregion. Used only in fieldatpoint now.
      function imaskregion(iregion)
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iregion
      imaskregion=IAND(iregion,ifield_mask)
      end
!*****************************************************************
      function inside_geom(mdims,x,i)
! Return integer 0 or 1 according to whether ndims-dimensional point x
! is outside or inside object number i. Return 0 for non-existent object.
      integer mdims,i
      real x(mdims)
! Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
      external interp

      inside_geom=0
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
! Use only bottom 8 bits:
      itype=itype-256*(itype/256)
! Start of inside type choices
      if(itype.eq.1)then
!------------------------------------------------
! Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims) 
         r2=0
         do k=1,ndims
            r2=r2+((x(k)-obj_geom(ocenter-1+k,i))/
     $           obj_geom(oradius-1+k,i))**2
         enddo
         if(r2.lt.1.)inside_geom=1
      elseif(itype.eq.2)then
!------------------------------------------------
! Coordinate-Aligned Cuboid data:
         do k=1,ndims
            xk=x(k)-obj_geom(ocenter-1+k,i)
            if(abs(xk).ge.abs(obj_geom(oradius-1+k,i)))return
         enddo
         inside_geom=1
      elseif(itype.eq.3)then
!------------------------------------------------
! Coordinate-Aligned Cylinder data:  Center(ndims), Semi-axes(ndims), 
! Axial coordinate. (Signed Axial 1/2 length 
! = semi-axis of the axial coordinate).
         ic=int(obj_geom(ocylaxis,i))
         xa=(x(ic)-obj_geom(ocenter+ic-1,i))
         if(abs(xa).ge.abs(obj_geom(oradius+ic-1,i))) return
         r2=0.
         do k=1,ndims
            if(k.ne.ic)
     $         r2=r2+((x(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i))**2
         enddo
!         write(*,*)'Cyl. ic=',ic,' r2=',r2,' x=',x
         if(r2.lt.1.)inside_geom=1
      elseif(itype.eq.4)then
!------------------------------------------------
! General Cuboid is equivalent to General Parallelopiped
! "center(ndims)", vectors(ndims,ndims), contravariants.
         do k=1,ndims
            xcj=0.
            do j=1,ndims
               xcj=xcj+(x(j)-obj_geom(ocenter+j-1,i))
     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
            enddo
            if(abs(xcj).gt.1.)return
         enddo
         inside_geom=1
      elseif(itype.eq.5)then
!------------------------------------------------
! Non-aligned cylinder
         r2=0.
         do k=1,ndims
            xcj=0.
            do j=1,ndims
               xcj=xcj+(x(j)-obj_geom(ocenter+j-1,i))
     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
!               write(*,'(2i2,3f10.4)')k,j,xcj,obj_geom(ocenter+j-1,i)
!     $              ,obj_geom(ocontra+ndims*(k-1)+j-1,i)
            enddo
            if(abs(xcj).gt.1.)return
! radius:
            if(k.lt.ndims)r2=r2+xcj**2
         enddo
!         write(*,*)r2
         if(abs(r2).lt.1.)inside_geom=1
      elseif(itype.eq.6)then
!------------------------------------------------
! Surface of revolution. Monotonic z.
! This proves to be about 10% faster than the general version 7,
! on a surface with one pair or 15% with 4 pairs, with a run whose
! time is dominated by cijroutine.
         r2=0.
         z=0
         do j=1,ndims
            z=z+obj_geom(ovax+j-1,i)*(x(j)-obj_geom(obase+j-1,i))
         enddo
         if(z.lt.0.or.z.gt.1)return
! Subtract from position z times axial vector-> cylindrical radius
         do j=1,ndims
            r2=r2+(x(j)-obj_geom(obase+j-1,i)-z*
     $           (obj_geom(ovec+2*ndims+j-1,i)))**2
         enddo
! Find the z-real-index 
         nq=int(obj_geom(onpair,i))
         iz=interp(obj_geom(opz,i),nq,z,zind)
         if(iz.le.0)then
            write(*,*)'inside_geom SoR interp failure'
            write(*,*)nq,z,zind,iz
            stop
         endif
! Then find the r-value for this z-index scaled to world coordinates.
         r=obj_geom(opr+iz-1,i)+(zind-iz)
     $        *(obj_geom(opr+iz,i)-obj_geom(opr+iz-1,i))
     $        *obj_geom(orscale,i)
         if(r2.lt.r**2)inside_geom=1
      elseif(itype.eq.7)then
!------------------------------------------------
! Surface of revolution. General.
         r=0.
         z=0.
         do j=1,ndims
            z=z+obj_geom(ovax+j-1,i)*(x(j)-obj_geom(obase+j-1,i))
         enddo
         do j=1,ndims
            r=r+(x(j)-obj_geom(obase+j-1,i)-z*
     $           (obj_geom(ovec+2*ndims+j-1,i)))**2
         enddo
         r=sqrt(r)/obj_geom(orscale,i)
         isect=int(w2sect(r,z,1.e5,0.,obj_geom(opr,i),obj_geom(opz,i)
     $        ,int(obj_geom(onpair,i)),fsect,psect))
         if(mod(isect,2).eq.1)inside_geom=1
      endif

      end
!************************************************************
! Return whether position x is inside the region specified by ibool.
      logical function linregion(ibool,ndims,x)
! ibool structure: n1, n1*values, n2, n2*values, ... ,0
      integer ibool(*)
      integer ndims
      real x(ndims)
      logical ltemp,lt1,lt2

!  linregion = Prod_1^nj Sum_1^ni inside(bool(ni,nj))
! where inside(n) is true if n is +/-ve and x is inside/outside |n|. 
!      write(*,'(11i4,3f10.4)')ibool(1:10),ndims,x
      lt2=.true.
! Special case for zero particle boolean.
      if(ibool(1).eq.0)goto 10
      lt1=.false.
      i=1
      n1=ibool(i)+i
 1    if(i.lt.n1)then
        i=i+1
         ib=ibool(i)
         if(ib.ne.0)then
            ltemp=inside_geom(ndims,x,abs(ib)).eq.1
            if(ib.lt.0)ltemp=.not.ltemp
            lt1=lt1.or.ltemp
            goto 1
         endif
      else
         i=i+1
         lt2=lt1.and.lt2
         lt1=.false.
         n1=i+ibool(i)
         if(i.lt.n1)goto 1
      endif
 10   linregion=lt2
      end
!************************************************************
      integer function leaveregion(ibool,ndims,x1,x2,icross,f)
! Return non-zero if we left the region in passing from x1 to x2.
! leaveregion is equal to the object whose intersection was detected.
!
! On entry 
!    icross has bits nonzero for objects possibly crossed, so zero bits
!      mask out tests. 
!    ibool contains the boolean defining the region.
!    f(*) must have sufficient length to accept multiple crossings.
!
! On exit 
!    f(1) contains the minimum leaving intersection fraction detected
!       or 1 if no intersections.
!    f contains in increasing order all the intersection fractions.
!    x2 is changed to be the position just past fraction f(1).

! A boolean block is violated iff there is a crossing of one object
! that occurs outside all the other booleans of the block. If any 
! boolean block is violated, then the ray leaves the region. Thus 
! for each crossing, for each block that contains the object crossed,
! check intersection for being inside all the other objects in the block.
! If there is none that it is inside, then it is a leaving intersection.
! We need to find the lowest crossing fraction (nearest to x1).
      integer ndims,ibool(*),icross
      real x1(ndims),x2(ndims)
      real f(*)
! The following ought to be consistent with 3dcom.f
      integer ibmax
      parameter (ibmax=100) 
      real fobs(ibmax)
      integer iobs(ibmax)
      real xx(3)
      logical linregion
      
! it is the object currently whose crossing is being examined
! ic represents the crossed objects by its bits
! ibl points to the boolean slot containing current block length 
! nbl is current block length
! io is increment position being examined within block 
! inblock says the crossing object is in this block 
! icit is the number of crossings of it
! ic is the current crossing being examined
      idebug=0
 99   it=0
      ic=icross
      leaveregion=0
      nleave=0
!      iobs(1)=0
      f(1)=1.
 1    if(ic.eq.0)goto 5
      icb2=ic/2
      it=it+1
      if((ic-2*icb2).ne.0)then
! Check this non-zero crossing against the entire boolean.
         ibl=1
 2       nbl=ibool(ibl)
         if(nbl.eq.0)goto 3
         inblock=0
         do io=1,nbl
! Search for this object within the block.
            iobj=ibool(ibl+io)
            if(abs(iobj).eq.it)then
               inblock=1               
               if(idebug.ne.0)write(*,*)'inblock',ibl,nbl,iobj,ic,it
            endif
         enddo
         if(inblock.ne.0)then
! This object is in this block. Test whether intersection is in others.
! Get back the position of the intersection(s) in fobs
            icit=icross_geom(x1,x2,it,fobs)
            if(idebug.ne.0)then
               write(*,*)'it,icit,fobs',it,icit,(fobs(j),j=1,icit)
               write(*,'(a,6f7.3)')'x1,x2',(x1(j),x2(j),j=1,ndims)
            endif
            do ic=1,icit
! Find the position of the intersection.
               do k=1,ndims
                  xx(k)=x1(k)+(x2(k)-x1(k))*fobs(ic)
               enddo
               do io=1,nbl
                  iobj=ibool(ibl+io)
                  if(abs(iobj).ne.it)then
                     is=inside_geom(ndims,xx,abs(iobj))
                     if(iobj.gt.0 .and. is.eq.1 .or.
     $                    iobj.lt.0 .and. is.eq.0)then
! Intersection is inside one of the other objects of block. 
! It therefore doesn't count as a leaving. Break to next intersection.
                        goto 4
                     endif
                  endif
               enddo
! Here when we have an intersection that is a real crossing.
! We crossed the particle region boundary here.
! Insert the fraction and the object number into stack ordered by fraction.
               call insertsorted2(nleave,f,fobs(ic),iobs,it)
! In order to find all leavings and choose the minimum fraction we mustnt
! return here. The next three statements revert to prior algorithm which
! agrees with earlier calculations for make test.
!               f(1)=fobs(ic)
!               leaveregion=it
!               return
!               write(*,*)'Would have left',f(1),it,icit
 4             continue
            enddo
         endif
! We found this block to be clean of leaving. Go to next bool block.
         ibl=ibl+nbl+1
!         write(*,*)it,ic,icb2,icit,ibl,n
!         write(*,*)(ibool(k),k=1,10)
         goto 2
      endif
 3    ic=icb2
! Iterate to next object crossed.
      goto 1
! Finished. f(i) 
 5    continue
      if(f(1).ne.1.)then
         leaveregion=iobs(1)
      endif
! Test for leakage. Not normally used. 
      if(.false.)then
      if(.not.linregion(ibool,ndims,x2).and.leaveregion.eq.0)then
         if(idebug.ne.1)then
            write(*,*)'Leaveregion error',icross,ibl,nbl
            write(*,*)'ibool',(ibool(i),i=1,5)
            write(*,*)nleave,(f(i),i=1,nleave)
! Repeat with diagnostics.
            idebug=1
            goto 99
         endif
      endif
      if(idebug.eq.1)stop
      endif

      end

!************************************************************
! Return whether the particle region specified by ibool has any
! "enclosed regions". An enclosed region is one which consists only of
! the inside of its objects ORed together. If it does, then cartesian
! reinjection is not correct.
      logical function lbounded(ibool,ifmask)
! ibool structure: n1, n1*values, n2, n2*values, ... ,0
      integer ibool(*)
      integer ifmask
      logical ltemp
!  linregion = Prod_1^nj Sum_1^ni inside(bool(ni,nj))
! where inside(n) is true if n is +/-ve and x is inside/outside |n|. 
      ltemp=.false.
      i=1
 1    nb=ibool(i)
      n1=nb+i
!      write(*,*)'lbounded outer',i,nb,n1
      if(nb.eq.0)goto 10
      ltemp=.true.
 2    continue
      if(i.lt.n1)then
! Reading objects for group ending at n1-1
         i=i+1
         ib=ibool(i)
!         write(*,*)'lbounded',i,ib,ltemp,ibits(ifmask,abs(ib)-1,1)
! Trap subtle object error.
! It's not clear that this is logically correct.
! Nonphysically we could banish particles from a non-field-boundary.
         if(.false..and.ibits(ifmask,abs(ib)-1,1).eq.0)then
            if(ib.gt.0)then
               write(*,'(a,6i4,a)')' Particle region boolean error.'
     $              ,(ibool(k),k=1,6),' ...'
               write(*,*)'Specified non-field-boundary object',i,ib
               call reportfieldmask()
               write(*,*)'You must fix the object file. Aborting ...'
               stop
            else
               write(*,*)'WARNING: particle region excluding a '
     $              ,'non-field-object',ib
            endif
         endif
         if(ib.lt.0)then
! A positive ib value defines inside object. If there are any outsides
! in this compound, it is not enclosed.
            ltemp=.false.
         endif
         goto 2
      else
         i=i+1
         goto 1
      endif
 10   lbounded=ltemp
!(ibool(k),k=1,16),
!      write(*,*)' lbounded=',lbounded
      end

!****************************************************************
! Convert from contravariant coefficients to world cartesian.
! Covariant and Contravariant vectors are stored in obj_geom(...,iobj)
! xcontra and xw can be the same storage positions if desired.
      subroutine contra3world(mdims,xcontra,xw,iobj)
! Number of dimensions. 
      integer mdims
! Object number
      integer iobj
! World coordinates
      real xw(mdims)
! Contravariant coordinates
      real xcontra(mdims)
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xwl(ndims)
! Cartesian obtained as sum of covariant vectors times contra coeffs.
      do i=1,ndims
         xwl(i)=obj_geom(ocenter+i-1,iobj)
! Contra 
         do j=1,ndims
            xwl(i)=xwl(i)
     $           +xcontra(j)*obj_geom(ovec+ndims*(j-1)+i-1,iobj)
!            write(*,*)i,j,obj_geom(ovec+ndims*(j-1)+i-1,iobj),xcontra(j)
         enddo
      enddo
      do i=1,ndims
         xw(i)=xwl(i)
      enddo
      end
!****************************************************************
      subroutine world3contra(mdims,xw,xcontra,iobj)
! Convert from world cartesian to contravariant coefficients 
! Covariant and Contravariant vectors are stored in obj_geom(...,iobj)
! xcontra and xw can overlap, if aligned.
      integer mdims,iobj
      real xw(mdims),xcontra(mdims)
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xd(ndims),xcl(ndims)

      do j=1,ndims
! Cartesian world relative to center
         xd(j)=xw(j)-obj_geom(ocenter+j-1,iobj)
      enddo
      do i=1,ndims
         xcl(i)=0.
! Contra coefficients are obtained by dotting with contra vectors
         do j=1,ndims
            xcl(i)=xcl(i)
     $           +xd(j)*obj_geom(ocontra+ndims*(i-1)+j-1,iobj)
         enddo
      enddo
      do i=1,ndims
         xcontra(i)=xcl(i)
      enddo
      end
!**********************************************************************     
      real function w2sect(r1,z1,r2,z2,rwall,zwall,nwall,fsect,psect)
! Find intersection between line joining two points and a wall. 2-D.
!     
! On input    
!     (r1,z1) and (r2,z2) are two points
!     rwall, zwall a possibly open piecewise-linear contour with length
!     nwall

! On output 
!     w2sect is the total number of intersections with segment 1->2
!     fsect the smallest fractional distance along 1->2 of an intersection.
!     psect the wall fractional index of this intersection.
!     fsect, psect =-1 if w2sect=0. 

!  Find the intersection of the line joining the two points and the
!  wall. Return the number of intersections w2sect. If w2sect.ne.0
!  find the fractional distance between rz1, rz2 of the intersection
!  closest to r1,z1: fsect, and the wall index to which the intersection
!  closest to r1,z1 corresponds, including the fractional distance along
!  the intersecting section, psect. (Otherwise they are both -1).

      real rwall(nwall),zwall(nwall)
!     
!     Coefficients of equation of line connecting two points:
      a=+(z1-z2)
      b=-(r1-r2)
      c= r1*z2-r2*z1
!     
!     For each wall segment, determine if it intersects the
!     line. The test for intersection is based on the fact that the
!     endpoints of each segment must lie on opposite sides of the equation
!     describing the line of the other segment.
!     
      icount=0
      psect=-1
      fsect=-1
! a large number
      valmin=1.e30
      val2=a*rwall(1)+b*zwall(1)+c
      do 10 i=1,nwall-1
         val1=val2
         val2=a*rwall(i+1)+b*zwall(i+1)+c
         if (((val1.le.0.).and.(val2.gt.0)).or.
     &        ((val1.gt.0).and.(val2.le.0))) then
            aprime=+(zwall(i+1)-zwall(i))
            bprime=-(rwall(i+1)-rwall(i))
            cprime= rwall(i+1)*zwall(i)-rwall(i)*zwall(i+1)
            val3=aprime*r1+bprime*z1+cprime
            val4=aprime*r2+bprime*z2+cprime
!     Prescribe that the second end of the line R2,Z2 does not count if it is
!     exactly on the boundary but the first end always does. 
            if (((val3.le.0.).and.(val4.gt.0)).or.
     &           ((val3.ge.0).and.(val4.lt.0))) then
!     section for intersecting wall sections
               icount=icount+1
               f=val3/(val3-val4)
               if(abs(f).lt.valmin)then
                  valmin=abs(f)
                  fsect=f
                  psect=i+ val1/(val1-val2)
               endif
            endif
         endif
 10   continue
      w2sect=icount
      return
      end
!
!**********************************************************************
      function icrossall(x1,x2)
! For two ndims-dimensional points x1,x2, return the integer icrossall
! consisting of bits i=0-30 that are zero or one according to whether
! the line x1-x2 crosses object i.
! Common object geometric data.
      include 'ndimsdecl.f'
      real x1(ndims),x2(ndims)
      include '3dcom.f'
      real f(10)
!
      icrossall=0
      do i=ngeomobj,1,-1
         icrossall=2*icrossall
         if(icross_geom(x1,x2,i,f).ne.0)icrossall=icrossall+1
      enddo
      end
!*****************************************************************
      integer function icross_geom(x1,x2,i,f)
! Return the number of intersections that the straight line from x1 to x2
! makes with the object i. Return 0 for a non-existent object.
! Place the intersect fractions in f.
      integer i
! Common object geometric data.
      include 'ndimsdecl.f'
      real x1(ndims),x2(ndims)
      real f(*)
      include '3dcom.f'
      real xn1(ndims),xn2(ndims)
      integer ids(ovlen)
      external interp

      icross_geom=0
!      write(*,*)'icross_geom',i,ngeomobj
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
! Use only bottom 8 bits:
      itype=itype-256*(itype/256)
! Start of type choices
      if(itype.eq.1)then
!------------------------------------------------
! Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims)
         call spheresect(ndims,0,x1,x2,obj_geom(ocenter,i)
     $        ,obj_geom(oradius,i),f(1),f(2),sd,C,D)
         if(sd.eq.0)return
! No intersection.
         if(f(1).gt.1. .or. f(1).lt.0.)return
! Nearest intersection is beyond the segment.
         if(f(2).gt.1. .or. f(2).lt.0.)then
            icross_geom=1
         else
            icross_geom=2
         endif
      elseif(itype.eq.2)then
!------------------------------------------------
! Coordinate-Aligned Cuboid:
! Normalize
         do k=1,ndims
            xn1(k)=(x1(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i)
            xn2(k)=(x2(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i)
         enddo
!         call cubenormsect(xn1,xn2,icross_geom,f)
         call cubeusect(xn1,xn2,icross_geom,f,ids)
      elseif(itype.eq.3)then
!------------------------------------------------
! Coordinate-Aligned Cylinder data:  Center(ndims), Semi-axes(ndims), 
! Axial coordinate. (Signed Axial 1/2 length 
! = semi-axis of the axial coordinate).
! Normalize
         ic=int(obj_geom(ocylaxis,i))
         do k=1,ndims
            kn=mod(2+k-ic,ndims)+1
            xn1(kn)=(x1(k)-obj_geom(ocenter+k-1,i))
     $           /obj_geom(oradius-1+k,i)
            xn2(kn)=(x2(k)-obj_geom(ocenter+k-1,i))
     $           /obj_geom(oradius-1+k,i)
         enddo
         call cylusect(xn1,xn2,i,icross_geom,f,ids)
      elseif(itype.eq.4)then
!------------------------------------------------
! General Cuboid is equivalent to General Parallelopiped
! "center(ndims)", vectors(ndims,ndims), contravariants.
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
! This may be faster than the above calls.
!         do k=1,ndims
!            xn1(k)=0.
!            xn2(k)=0.
!            do j=1,ndims
!               xn1(k)=xn1(k)+(x1(j)-obj_geom(ocenter+j-1,i))
!     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
!               xn2(k)=xn2(k)+(x2(j)-obj_geom(ocenter+j-1,i))
!     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
!            enddo
!         enddo
!         call cubenormsect(xn1,xn2,icross_geom,f)
         call cubeusect(xn1,xn2,icross_geom,f,ids)
      elseif(itype.eq.5)then
!------------------------------------------------
! Non-aligned cylinder
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
         call cylusect(xn1,xn2,i,icross_geom,f,ids)
!         call cylnormsect(xn1,xn2,icross_geom,f)
      elseif(itype.eq.6.or.itype.eq.7)then
!------------------------------------------------
! Surface of revolution. General. Scale to contravariant components.
! Then use svrsect. It is normalized to reference dimensions, but
! the SoR parameters call for renormalizations internal to srvsect.
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
         call srvsect(xn1,xn2,i,icross_geom,f,ids)
      endif

      end
!*********************************************************************
      subroutine insertsorted(n,f,fnew)
! Increment n and insert the value fnew into the ordered list f(*) of n
! increasing values. 
      integer n
      real fnew,f(*)
      n=n+1
      do i=n,2,-1
         if(f(i-1).gt.fnew)then
            f(i)=f(i-1)
         else
            f(i)=fnew
            goto 1
         endif
      enddo
      f(1)=fnew
 1    continue
      end
!*********************************************************************
      subroutine srvsect(xp1,xp2,iobj,icross_geom,f,ids)
! Given a general SoR object iobj. Count the number of intersections of
! the line joining points xn1,xn2 (contravariant components) with
! it. This version uses explicit cone intersection, not iteration.

      integer iobj
      include 'ndimsdecl.f'
      real xp1(ndims),xp2(ndims)
      real f(*)
      integer ids(*)
      include '3dcom.f'
      include 'dbgcom.f'
! Local storage`
      real f1,f2
      real xn1(ndims),xn2(ndims),xr(ndims)
      equivalence (xn1(ndims),z1),(xn2(ndims),z2)
      real zero(ndims)
      data zero/ndims*0./

! Example of choosing a particular particle case
      icross_geom=0
      f(1)=1.
      do i=1,ndims
         xn1(i)=xp1(i)
         xn2(i)=xp2(i)
      enddo
!      write(*,*)'srvsect',xp1,xp2
!      write(*,*)z1,z2
! For each SoR-segment find if the line xp1-xp2 intersects it.
! Order the two wall positions increasing.
      do i=1,int(obj_geom(onpair,iobj))-1
         if(obj_geom(opz+i,iobj).ge.obj_geom(opz+i-1,iobj))then
            rw1=obj_geom(opr+i-1,iobj)
            zw1=obj_geom(opz+i-1,iobj)
            rw2=obj_geom(opr+i,iobj)
            zw2=obj_geom(opz+i,iobj)
         else
            rw2=obj_geom(opr+i-1,iobj)
            zw2=obj_geom(opz+i-1,iobj)
            rw1=obj_geom(opr+i,iobj)
            zw1=obj_geom(opz+i,iobj)
         endif
         zd=zw2-zw1
         rd=rw2-rw1
         rc=(rw2+rw1)/2.
         if((z1.lt.zw1).and.(z2.lt.zw1).or.(z1.gt.zw2).and.(z2.gt.zw2))
     $        then
! Non intersecting: points both off the same end of segment.
! Hopefully this case is the dominant one and quick. Continue to next.
!               write(*,*)'Both off end',z1,z2
         elseif(abs(zd).lt. 1.e-6*rc)then
! Near Degenerate: disc identical z. Alternate treatment, unscaled.
            if(idbug.gt.0)write(*,*)'small zd',zd
            fr=(zw1-z1)/(z2-z1)
            if(fr.ge.0. .and. fr.lt.1.)then
               r=sqrt((xp1(1)+fr*(xp2(1)-xp1(1)))**2
     $              +(xp1(2)+fr*(xp2(2)-xp1(2)))**2)
               if(r.gt.min(rw1,rw2) .and. r.le.max(rw1,rw2))then
                  call insertsorted2(icross_geom,f,fr,ids,i)
               endif
            endif
         else
            sd=0.
            if(abs(rd).lt.1.e-8*rc)then
! Degenerate: cylinder. Use cylinder code.
!               write(*,*)'Cylinder approximation'
               if(abs(zd).ge.2.e-8*rc)then
                  do k=1,ndims-1
                     xr(k)=rc
                  enddo
                  xr(ndims)=1.
! Find intersections with curved surface by projection. In direction 3.
                  call spheresect(ndims,3,xn1,xn2,zero,xr,f1,f2,sd,C
     $                 ,D)
               endif
            else
! Standard cone. zscale=dz/dr.
               dzdr=zd/rd
! Vertex position
               zx=zw1-rw1*dzdr
! Scale z to the unit cone. z1 and z2 equivalence xn?(3)
               z1=(z1-zx)/dzdr
               z2=(z2-zx)/dzdr
               call conesect(xn1,xn2,f1,f2,sd)
               if(idbug.gt.0)write(*,'(a,i3,6f10.2)'
     $              )'i,zx,z1,z2,f1,f2,sd',i,zx,z1,z2,f1,f2,sd
! Restore zs for next cone, and world calculation.
               z1=xp1(ndims)
               z2=xp2(ndims)
            endif
            if(sd.ne.0.)then
            if(f1.ge.0. .and. f1.lt.1.)then
! Here we need to tell if the intersection is within the face or beyond
! it on the cone. The criterion needs to be robust to rounding errors.
! Take it to be that the 2d position of the intersection observes the
! wall element as an oblique angle.  
               zx=z1+f1*(z2-z1)
               if(idbug.gt.0)write(*,*)'zx1,zw1,zw2',zx,zw1,zw2
               rx=sqrt((xp1(1)+f1*(xp2(1)-xp1(1)))**2
     $              +(xp1(2)+f1*(xp2(2)-xp1(2)))**2)
               if((zw1-zx)*(zw2-zx)+(rw1-rx)*(rw2-rx).lt.0.)then
                  if(idbug.ne.0)write(*,*)'Insert f1',f1
                  call insertsorted2(icross_geom,f,f1,ids,i)
               endif
            endif
            if(f2.ge.0. .and. f2.lt.1.)then
               zx=z1+f2*(z2-z1)
               if(idbug.gt.0)write(*,*)'zx1,zw1,zw2',zx,zw1,zw2
               rx=sqrt((xp1(1)+f2*(xp2(1)-xp1(1)))**2
     $              +(xp1(2)+f2*(xp2(2)-xp1(2)))**2)
               if((zw1-zx)*(zw2-zx)+(rw1-rx)*(rw2-rx).lt.0.)then
                  if(idbug.ne.0)write(*,*)'Insert f2',f2
                  call insertsorted2(icross_geom,f,f2,ids,i)
               endif
            endif
            endif
         endif
      enddo
      end
!*********************************************************************
      subroutine conesect(xp1,xp2,f1,f2,sd)
! Find the fractional intersection points of the line joining xp1 to xp2
! with the unit cone x^2+y^2=z^2. If xd=xp2-xp1, then the solution of
! the intersection is x=xp1+f*xd, with Af^2+2Bf+C=0. Where
! A=xd^2+yd^2-zd^2, B=xd.xp1+yd.yp1-zd.zp1, C=xp1^2+yp1^2-zp1^2 So f=
! (-B +- sqrt(B^2-AC) )/A Return f=1 for no intersection. f1 should
! be always positive if possible, and closest to zero. sd=1 if point 1
! is outside, meaning r1>z1 (C>0), and the first-f crossing is
! inward. sd=-1 if the first crossing is outward. sd=0 if there are no
! crossings even of the extrapolated line.
      include 'dbgcom.f'
      integer ndims
      parameter (ndims=3)
      real xp1(ndims),xp2(ndims),f1,f2,sd
      real xd(ndims)
      real A,B,C,D,DS,One
      data One/1.0/

      if(idbug.ne.0)write(*,*)'Conesect xp1,xp2',xp1,xp2
      do i=1,ndims
         xd(i)=xp2(i)-xp1(i)
      enddo
      A=xd(1)**2+xd(2)**2-xd(3)**2
      B=xd(1)*xp1(1)+xd(2)*xp1(2)-xd(3)*xp1(3)
      C=xp1(1)**2+xp1(2)**2-xp1(3)**2
      D=xd(3)**2*(xp1(1)**2+xp1(2)**2)+xp1(3)**2*(xd(1)**2+xd(2)**2)
     $     -2.*xd(3)*xp1(3)*(xd(1)*xp1(1)+xd(2)*xp1(2))
     $     -(xd(1)*xp1(2)-xd(2)*xp1(1))**2
! D =B**2-A*C evaluated carefully to avoid rounding.
!      D=B**2-A*C   ! This gave particle leakage.
      if(idbug.ne.0)write(*,*)'A,B,C,D,B^2-AC',A,B,C,D,B**2-A*C
      if(D.lt.0)then
! No intersections.
         f1=1.
         f2=1.
         sd=0.
      elseif(A.eq.0)then
! One intersection
         f1=-0.5*C/B
         f2=1.
         sd=sign(One,C)
      else
! Two intersections.
          DS=sqrt(D)
         if( (A.gt.0. .and. C.gt.0. .and. B.lt.0) .or.
     $        (A.lt.0. .and.(C.ge.0. .or. B.le.0)) )then
! Can take minus sign and still get positive for A>0.
            f1=(-B-DS)/A
            f2=(-B+DS)/A
         else
! Can take plus sign and still get positive for A<0.
            f1=(-B+DS)/A
            f2=(-B-DS)/A
         endif
         sd=sign(One,C)
      endif
      end
!********************************************************************
      subroutine conesectz(xp1,xp2,z1,z2,sd)
! Find the z-positions of the intersection of the line from xp1 to xp2
! with the unit cone x^2+y^2=z^2. 

! They are the solution of Az^2+Bz+C=0, with
! A=x'^2+y'^2-1;  B=2[x_0x'+y_0y'];  C= x_0^2+y_0^2,
! where ' means d/dz, and x_0=x1-z1.x' and y_0=y1-z1.y'.  

! Logic of the sorted order needs work. Not done. 
      integer ndims
      parameter (ndims=3)
      real xp1(ndims),xp2(ndims),z1,z2,sd
      real A,B,C,disc,One
      data One/1.0/
      zd=xp2(3)-xp1(3)
      xp=(xp2(1)-xp1(1))/zd
      yp=(xp2(2)-xp1(2))/zd
      x0=xp1(1)-xp1(3)*xp
      y0=xp1(2)-xp1(3)*yp
      A=xp**2+yp**2-1
      B=2.*(x0*xp+y0*yp)
      C=x0**2+y0**2
      disc=B**2-A*C
      if(disc.lt.0)then
! No intersections.
         sd=0.
         z1=0.
         Z2=0.
      elseif(A.eq.0)then
! One intersection
         z1=-0.5*C/B
         z2=0.
         sd=sign(One,C)
      else
! Two intersections.
         disc=sqrt(disc)
         if( (A.gt.0. .and. C.gt.0. .and. B.lt.0) .or.
     $        (A.lt.0. .and.(C.ge.0. .or. B.le.0)) )then
! Can take minus sign and still get positive for A>0.
            z1=(-B-disc)/A
            z2=(-B+disc)/A
         else
! Can take plus sign and still get positive for A<0.
            z1=(-B+disc)/A
            z2=(-B-disc)/A
         endif
         sd=sign(One,C)
      endif

      end
! Obsolete Below here.
!*********************************************************************
      subroutine cylnormsect(xn1,xn2,icross_geom,f)
! Find the number of intersections of the line segment joining xn1,xn2,
! with the unit (half-length and radius) cylinder.

! This code is logically incorrect and now not used anywhere. The problem
! is that the use of f() for ends is overwritten by the spheresect
! search.

      include 'ndimsdecl.f'
      integer icross_geom
      real xn1(ndims),xn2(ndims)
      real zero(ndims),one(ndims)
      real f(*)
      data zero/ndims*0./one/ndims*1/

      icross_geom=0
!      write(*,*)'Cylnormsect entry',xn1,xn2
! On entry xn1,2 are the coordinates relative to the unit cylinder.
      z1=xn1(ndims)
      z2=xn2(ndims)
! Test if both off same end, to avoid extra work.
      if(z1.gt.1. .and. z2.gt.1. .or. z1.lt.-1..and.z2.lt.-1)then
         return
      endif
! The square roots here are not necessary. Nor the abs(r..).
      r1=sqrt(xn1(1)**2+xn1(2)**2)
      r2=sqrt(xn2(1)**2+xn2(2)**2)
      if(abs(z1).lt.1. .and. abs(z2).lt.1.
     $     .and.abs(r1).lt.1. .and. abs(r2).lt.1.)then
! 1,2 both inside
         return
      endif
! Not both inside. Need more interpolation.
! z-faces:
      ida=3
      iend=0
      xd=z2-z1
      if(abs(xd).lt.1.e-20)xd=1.e-20
      do j=-1,1,2
         fr=(j-z1)/xd
         write(*,*)'xd, fr', xd,fr,iend+1
         if(fr.ge.0.and.fr.le.1.)then
! Intersects this plane. Test whether r at intersect is inside
            r=0.
            do i=mod(ida,ndims)+1,mod(ida+ndims-2,ndims)+1
               r=r+(xn1(i)+fr*(xn2(i)-xn1(i)))**2
            enddo
!            write(*,*)'Plane',j,' r^2=',r
            if(r.gt.1.)goto 32
! Intersection is inside circle. It counts.
            iend=iend+1
            f(iend)=fr
         endif
 32      continue
      enddo
      icross_geom=iend
      if(icross_geom.ne.2)then
! Still looking. Find the intersection (if any) with the
! Curved surface by projecting along ida and testing a circle
! of center 0 and radius 1.
!         write(*,*)'Calling spheresect'
         call spheresect(ndims,ida,xn1,xn2,zero,one,f(1),f(2),sd,C,D)
         do j=1,2
            if(f(j).ge.0.and.f(j).lt.1.)then
! The intersection j is on the circle. At relevant z?
               zx=xn1(ndims)+f(j)*(xn2(ndims)-xn1(ndims))
!               write(*,*)'j,zx,f(j)',j,zx,f(j)
               if(zx.lt.1. .and. zx.gt.-1.)then
                  icross_geom=icross_geom+1
               endif
            endif
         enddo
         if(iend.eq.1)then
! Put back the end intersection in the correct order.
            write(*,*)'fr,f(1)',fr,f(1)
            if(fr.lt.f(1))then
               f(2)=f(1)
               f(1)=fr
            else
               f(2)=fr
            endif
         endif
      endif
      end
!*******************************************************************
      subroutine cubenormsect(xn1,xn2,icross_geom,f)
! Find the number of intersections of the line segment from xn1 to xn2
! with the centered normalized cube of half-side length 1. This code
! ought to be equivalent to cubeusect except that it does not require 
! the plane on which the intersections occur to be tracked and so lacks
! cubeusect's last parameter. It also searches planes in different order.
      include 'ndimsdecl.f'
      real xn1(ndims),xn2(ndims)
      real f(*)
      icross_geom=0
      f(1)=1.
! Find intersection(s) with bounding planes. If any lies inside
! the square defined by the orthogonal planes, there's an intersection.
      do k=1,ndims
         xd=xn2(k)-xn1(k)
         if(xd.ne.0)then
! It is possible the line intersects planes of constant k-coordinate.
            do j=-1,1,2
               fr=(j-xn1(k))/xd
               if(fr.ge.0.and.fr.le.1.)then
! Intersects this plane. Test whether other dimensions are inside
                  do ii=1,ndims-1
                     i=mod(k+ii-1,ndims)+1
                     xi=xn1(i)+fr*(xn2(i)-xn1(i))
                     if(abs(xi).gt.1)goto 22
! This prevents duplicates when the intersection is at an edge:
                     if(icross_geom.ge.1.and.fr.eq.f(icross_geom))goto
     $                    22
                  enddo
! Intersection is inside all other planes. It counts. Put it into the
! right place in the list of fs, i.e. with the lowest positive first.
                  call insertsorted(icross_geom,f,fr)
               endif
! Intersection is not inside all others. Invalid.
 22            continue
            enddo
         endif
      enddo
      end
