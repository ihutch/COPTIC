c**********************************************************************
      function insideall(mdims,x)
c For an ndims-dimensional point x, return the integer insideall
c consisting of bits i=0-30 that are zero or one according to whether
c the point x is outside or inside object i.
      integer mdims
      real x(mdims)
c Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
c
      insideall=0
      do i=ngeomobj,1,-1
         insideall=2*insideall+inside_geom(ndims,x,i)
      enddo
c      if(insideall.ne.0) write(*,*)
c     $     'ngeomobj=',ngeomobj,' insideall=',insideall,x

      end
c****************************************************************
      function insidemask(mdims,x)
c For an ndims-dimensional point x, return the integer insidemask
c consisting of bits i=1-31 that are zero or one according to whether
c the point x is outside or inside object i. But bits are masked by
c the regionmask.
      integer mdims
      real x(mdims)
      intrinsic btest
c Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
c
      insidemask=0
      ibit=ngeomobj
      do i=ngeomobj,1,-1
         insidemask=2*insidemask
         if(btest(ifield_mask,i-1))
     $        insidemask=insidemask+inside_geom(ndims,x,i)
      enddo

      end
c*****************************************************************
c Return the masked iregion. Used only in fieldatpoint now.
      function imaskregion(iregion)
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iregion
      imaskregion=IAND(iregion,ifield_mask)
      end
c*****************************************************************
      function inside_geom(mdims,x,i)
c Return integer 0 or 1 according to whether ndims-dimensional point x
c is outside or inside object number i. Return 0 for non-existent object.
      integer mdims,i
      real x(mdims)
c Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
      external interp

      inside_geom=0
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)
c Start of inside type choices
      if(itype.eq.1)then
c------------------------------------------------
c Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims) 
         r2=0
         do k=1,ndims
            r2=r2+((x(k)-obj_geom(ocenter-1+k,i))/
     $           obj_geom(oradius-1+k,i))**2
         enddo
         if(r2.lt.1.)inside_geom=1
      elseif(itype.eq.2)then
c------------------------------------------------
c Coordinate-Aligned Cuboid data:
         do k=1,ndims
            xk=x(k)-obj_geom(ocenter-1+k,i)
            if(abs(xk).ge.abs(obj_geom(oradius-1+k,i)))return
         enddo
         inside_geom=1
      elseif(itype.eq.3)then
c------------------------------------------------
c Coordinate-Aligned Cylinder data:  Center(ndims), Semi-axes(ndims), 
c Axial coordinate. (Signed Axial 1/2 length 
c = semi-axis of the axial coordinate).
         ic=int(obj_geom(ocylaxis,i))
         xa=(x(ic)-obj_geom(ocenter+ic-1,i))
         if(abs(xa).ge.abs(obj_geom(oradius+ic-1,i))) return
         r2=0.
         do k=1,ndims
            if(k.ne.ic)
     $         r2=r2+((x(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i))**2
         enddo
c         write(*,*)'Cyl. ic=',ic,' r2=',r2,' x=',x
         if(r2.lt.1.)inside_geom=1
      elseif(itype.eq.4)then
c------------------------------------------------
c General Cuboid is equivalent to General Parallelopiped
c "center(ndims)", vectors(ndims,ndims), contravariants.
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
c------------------------------------------------
c Non-aligned cylinder
         r2=0.
         do k=1,ndims
            xcj=0.
            do j=1,ndims
               xcj=xcj+(x(j)-obj_geom(ocenter+j-1,i))
     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
c               write(*,'(2i2,3f10.4)')k,j,xcj,obj_geom(ocenter+j-1,i)
c     $              ,obj_geom(ocontra+ndims*(k-1)+j-1,i)
            enddo
            if(abs(xcj).gt.1.)return
c radius:
            if(k.lt.ndims)r2=r2+xcj**2
         enddo
c         write(*,*)r2
         if(abs(r2).lt.1.)inside_geom=1
      elseif(itype.eq.6)then
c------------------------------------------------
c Surface of revolution. Monotonic z.
c This proves to be about 10% faster than the general version 7,
c on a surface with one pair or 15% with 4 pairs, with a run whose
c time is dominated by cijroutine.
         r2=0.
         z=0
         do j=1,ndims
            z=z+obj_geom(ovax+j-1,i)*(x(j)-obj_geom(obase+j-1,i))
         enddo
         if(z.lt.0.or.z.gt.1)return
c Subtract from position z times axial vector-> cylindrical radius
         do j=1,ndims
            r2=r2+(x(j)-obj_geom(obase+j-1,i)-z*
     $           (obj_geom(ovec+2*ndims+j-1,i)))**2
         enddo
c Find the z-real-index 
         nq=int(obj_geom(onpair,i))
         iz=interp(obj_geom(opz,i),nq,z,zind)
         if(iz.le.0)then
            write(*,*)'inside_geom SoR interp failure'
            write(*,*)nq,z,zind,iz
            stop
         endif
c Then find the r-value for this z-index scaled to world coordinates.
         r=obj_geom(opr+iz-1,i)+(zind-iz)
     $        *(obj_geom(opr+iz,i)-obj_geom(opr+iz-1,i))
     $        *obj_geom(orscale,i)
         if(r2.lt.r**2)inside_geom=1
      elseif(itype.eq.7)then
c------------------------------------------------
c Surface of revolution. General.
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
         isect=w2sect(r,z,1.e5,0.,obj_geom(opr,i),obj_geom(opz,i)
     $        ,int(obj_geom(onpair,i)),fsect,psect)
         if(mod(isect,2).eq.1)inside_geom=1
      endif

      end
c************************************************************
c Return whether position x is inside the region specified by ibool.
      logical function linregion(ibool,ndims,x)
c ibool structure: n1, n1*values, n2, n2*values, ... ,0
      integer ibool(*)
      integer ndims
      real x(ndims)
      logical ltemp,lt1,lt2
c The following ought to be consistent with 3dcom.f
c      integer ibmax
c      parameter (ibmax=100) 

c  linregion = Prod_1^nj Sum_1^ni inside(bool(ni,nj))
c where inside(n) is true if n is +/-ve and x is inside/outside |n|. 
c      write(*,'(11i4,3f10.4)')ibool(1:10),ndims,x
      lt2=.true.
c Special case for zero particle boolean.
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
c************************************************************
c Return whether the particle region specified by ibool has any
c enclosed regions (regions inside an object. If it does, then
c cartesian reinjection is not correct. 
      logical function lbounded(ibool,ifmask)
c ibool structure: n1, n1*values, n2, n2*values, ... ,0
      integer ibool(*)
      integer ifmask
      logical ltemp
c  linregion = Prod_1^nj Sum_1^ni inside(bool(ni,nj))
c where inside(n) is true if n is +/-ve and x is inside/outside |n|. 
      ltemp=.false.
c Special case for zero particle boolean.
      if(ibool(1).eq.0)goto 10
      i=1
      n1=ibool(i)+i
 1    if(i.lt.n1)then
c Reading objects for group ending at n1-1
        i=i+1
         ib=ibool(i)
c Trap subtle object error.
         if(ibits(ifmask,abs(ib)-1,1).eq.0)then
            if(ib.gt.0)then
               write(*,*)'Particle region boolean error.'
     $              ,(ibool(k),k=1,6),' ...'
               write(*,*)'Specified non-field-boundary object',ib
               call reportfieldmask()
               write(*,*)'You must fix the object file. Aborting ...'
               stop
            else
               write(*,*)'WARNING: particle region excluding a '
     $              ,'non-field-object',ib
            endif
         endif
         if(ib.gt.0)then
c A positive ib value defines inside object.
            ltemp=.true.
         endif
         goto 1
      else
         i=i+1
         n1=i+ibool(i)
         if(i.lt.n1)goto 1
      endif
 10   lbounded=ltemp
c(ibool(k),k=1,16),
c      write(*,*)' lbounded=',lbounded
      end

c************************************************************
      subroutine potlsect(id,ipm,mdims,indi,fraction,conditions,dp,
     $     iobjno,ijbin)

c In dimension id, direction ipm, 
c from mesh point at indi(ndims) (zero-based indices, C-style),
c find any intersection of the mesh leg from this point to its neighbor
c with a bounding surface. Return the "fraction" of the leg at which
c the intersection occurs (1 if no intersection), the "conditions" at
c the intersection (irrelevant if fraction=1), the +ve length
c in computational units of the full leg in dp, and the object number
c in iobjno
      integer id,ipm,mdims,iobjno,ijbin
      integer indi(mdims)
      real fraction,dp
      real conditions(3)
c Equivalence: c1,c2,c3==a,b,c
c The boundary conditions are in the form c1\psi + c2\psi^\prime + c3
c =0.  Where ci=conditions(i).  Thus a fixed potential is (e.g.) c1=1
c c3=-\psi.  A fixed gradient is (e.g.) c2=1 c3=-\psi^\prime. 
c If c1=c2=c3=0 then no BC is applied and this object is skipped.
c
c In multiple dimensions, if the condition is placed on the gradient
c normal to the surface, in direction n, then for axis direction i
c (assuming an orthogonal set) the value of c2 should become c2_i =
c B/(n.e_i) where e_i is the unit vector in the axis direction.
c
c [If n.e_i=0 this is a pathological intersection.  The mesh leg is
c tangential to the surface, and a normal-gradient BC is irrelevant.
c Therefore, a fraction of 1 should be returned.]
c
c If the surface conditions include non-zero c2, then it is not possible
c to apply the same BC to both sides, without a discontinuity arising.
c Continuity alone should be applied on the other (inactive) side of the
c surface.  A negative value of c2 is used to indicate that just the
c appropriate continuity condition is to be applied, but c1,c2,c3,
c should all be returned with the magnitudes that are applied to the
c active side of the surface, with consistently reversed
c signs. (e.g. a=1, b=1, c=-2 -> a=-1, b=-1, c=2)
c
c A fraction of 1 causes all the bounding conditions to be ignored.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include '3dcom.f'
c      integer debug
      save numreports
      data numreports/0/

      real xx(ndims),xd(ndims)
      real xp1(ndims),xp2(ndims)
c Default no intersection.
      fraction=1
c      debug=0
c      if(iobjno.gt.100)debug=iobjno-100
c Silence warnings.
      istype=0
      iobjno=0

c-------------------------------------------------------------
c Store the xp1 and xp2 positions of the point and its neighbor
c so we are in a position to call the flux routines.
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         xp1(i)=xn(ix1)
         xp2(i)=xn(ix1)
         if(i.eq.id)then
            ix2=ix1+ipm
            xp2(i)=xn(ix2)
            dp=abs(xp2(i)-xp1(i))
         endif
      enddo
c-------------------------------------------------------------
c Process data stored in obj_geom.
      do i=1,ngeomobj
         if(obj_geom(oabc,i).ne.0. .or. obj_geom(oabc+1,i).ne.0. .or.
     $        obj_geom(oabc+2,i).ne.0.)then
c Only for non-null BCs
c Find the fractional intersection point if any.
c            if(debug.gt.0)fraction=101
            itype=int(obj_geom(otype,i))
            istype=itype/256
            itype=itype-256*istype
            if(itype.eq.1)then
c First implemented just for spheres.
               call spherefsect(ndims,xp1,xp2,i,ijbin,sd,fraction
     $              ,ijbin2)
            elseif(itype.eq.2)then
c Coordinate-aligned cuboid.               
               call cubefsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.3)then
c Coordinate aligned cylinder.
               call cylfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.4)then
c Parallelopiped.
               call pllelofsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
c     I don't think this ever actually happens, but we'll see.
               if(fraction.gt.1.or.fraction.lt.0.)then
                  write(*,*)'***********pllelofsect fraction',fraction
                  fraction=1.
               endif
            elseif(itype.eq.5)then
c Non-aligned cylinder.
               call cylgfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.6)then
               call srvfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.7)then
               call srvfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            else
               write(*,*)"Unknown object type",obj_geom(otype,i),
     $              " in potlsect"
            endif
         else
c Null BC. Ignore.
            fraction=1.
         endif

c Now we have decided the fraction for this object. Test if there's
c some subtraction for this object
         if(fraction.ne.1. .and. normv(i).ne.0)then
            if(numreports.lt.1)then
               write(*,*)'Subtraction from obj',i,'; total',normv(i)
     $           ,' objects subtracted:',(ormv(k,i),k=1,normv(i))
c     $              ,' indi',indi
c     $              ,' ipm',ipm
               numreports=numreports+1
            endif
c If this intersection is subtracted, exclude by setting fraction=1.
c Construct the intersection.
            do k=1,ndims
               xx(k)=xp1(k)+fraction*(xp2(k)-xp1(k))
            enddo
c Test if +-inside objects.
            do k=1,normv(i)
               isubtr=abs(ormv(k,i))
               inout=sign(1,ormv(k,i))
               inout=inout*(2*inside_geom(ndims,xx,isubtr)-1)
               if(inout.eq.1)then
c                  write(*,*)'Removed',indi,i,k,xx,fraction
                  fraction=1.
               endif
            enddo
c Ad hoc diagnostics
c            if(xx(1)**2+xx(2)**2.lt.1.)then
c               ii=inside_geom(ndims,xx,2)
c               write(*,*)'Unremoved',i,(ormv(k,i),ii,k=1,normv(i)),xx
c     $              ,fraction
c            endif
         endif

         if(fraction.ne.1.)then
c Default:
            projection=0.
            conditions(1)=obj_geom(oabc,i)
            conditions(2)=obj_geom(oabc+1,i)
            conditions(3)=obj_geom(oabc+2,i)
c            write(*,'(a,$)')'Setting conditions. '
            if(istype.eq.3)then
c            if(.false.)then
c Gradient of conditions exists. Calculate the local conditions.
               do k=1,ndims
                  xi=(xp1(k)-obj_geom(ocenter+k-1,i))*(1.-fraction)
     $                 +(xp2(k)-obj_geom(ocenter+k-1,i))*fraction
                  do ic=1,3
c gradients are in reverse order
                     conditions(ic)=conditions(ic)
     $                    +xi*obj_geom(oagrad+k-1-(ic-1)*ndims,i)
                  enddo
               enddo
c               write(*,*)'Gradient. cond=',conditions,i
            endif
            if(conditions(2).ne.0)then
c Derivative term present. Calculate the projection:
c    Calculate Coordinates of crossing
c Unnecessarily complicated:
c               do idi=1,ndims
c    Address of mesh point. 
c                  ix=indi(idi)+ixnp(idi)+1
c                  xx(idi)=xn(ix)
c This needs to be generalized for other objects.
c                  if(idi.eq.id)
c     $                 xx(idi)=xx(idi)+fraction*(xn(ix+ipm)-xn(ix))
c    Component of outward vector = Sphere outward normal.
c                  xd(idi)=xx(idi)-obj_geom(ocenter+idi-1,i)
c               enddo
               do k=1,ndims
                  xx(k)=xp1(k)+fraction*xp2(k)
c    Component of outward vector = Sphere outward normal.
                  xd(k)=xx(k)-obj_geom(ocenter+k-1,i)
               enddo
c    Signed projection cosine:
               projection=ipm*(xd(id)/obj_geom(oradius,i))
               if(projection.eq.0.)then
                  fraction=1.
               else
c I'm a bit nervous about when conditions are negative. Is this right?
c    Continuity works for both directions.
                  conditions(1)=sign(conditions(1),projection)
                  conditions(2)=conditions(2)/projection
                  conditions(3)=sign(conditions(3),projection)
c Special Zero outside, rather than continuity alternative
                  if(int(obj_geom(otype,i))/256.eq.1
     $                 .and. projection.lt.0.)then
                     conditions(1)=1.
                     conditions(2)=0.
                     conditions(3)=0.
                  endif
               endif
            endif
c            write(*,*)'ABC,projection',(obj_geom(oabc+k,i),k=0,2)
c     $           ,projection,conditions
            iobjno=i
            return
         endif
      enddo
c----------------------------------------------------------
c Default return for no intersection.
c Address of mesh point. 
      ix=indi(id)+ixnp(id)+1
      dp=abs(xn(ix+ipm)-xn(ix))
      end
c****************************************************************
c****************************************************************
c Convert from contravariant coefficients to world cartesian.
c Covariant and Contravariant vectors are stored in obj_geom(...,iobj)
c xcontra and xw can be the same storage positions if desired.
      subroutine contra3world(mdims,xcontra,xw,iobj)
c Number of dimensions. 
      integer mdims
c Object number
      integer iobj
c World coordinates
      real xw(mdims)
c Contravariant coordinates
      real xcontra(mdims)
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xwl(ndims)
c Cartesian obtained as sum of covariant vectors times contra coeffs.
      do i=1,ndims
         xwl(i)=obj_geom(ocenter+i-1,iobj)
c Contra 
         do j=1,ndims
            xwl(i)=xwl(i)
     $           +xcontra(j)*obj_geom(ovec+ndims*(j-1)+i-1,iobj)
c            write(*,*)i,j,obj_geom(ovec+ndims*(j-1)+i-1,iobj),xcontra(j)
         enddo
      enddo
      do i=1,ndims
         xw(i)=xwl(i)
      enddo
      end
c****************************************************************
      subroutine world3contra(mdims,xw,xcontra,iobj)
c Convert from world cartesian to contravariant coefficients 
c Covariant and Contravariant vectors are stored in obj_geom(...,iobj)
c xcontra and xw can overlap, if aligned.
      integer mdims,iobj
      real xw(mdims),xcontra(mdims)
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xd(ndims),xcl(ndims)

      do j=1,ndims
c Cartesian world relative to center
         xd(j)=xw(j)-obj_geom(ocenter+j-1,iobj)
      enddo
      do i=1,ndims
         xcl(i)=0.
c Contra coefficients are obtained by dotting with contra vectors
         do j=1,ndims
            xcl(i)=xcl(i)
     $           +xd(j)*obj_geom(ocontra+ndims*(i-1)+j-1,iobj)
c            write(*,*)i,j,obj_geom(ocontra+ndims*(j-1)+i-1,iobj)
c     $           ,xd(i)
         enddo
      enddo
      do i=1,ndims
         xcontra(i)=xcl(i)
      enddo
      end
C**********************************************************************     
      real function w2sect(r1,z1,r2,z2,rwall,zwall,nwall,fsect,psect)
c Find intersection between line joining two points and a wall. 2-D.
c     
c On input    
c     (r1,z1) and (r2,z2) are two points
c     rwall, zwall a possibly open piecewise-linear contour with length
c     nwall

c On output 
c     w2sect is the total number of intersections with segment 1->2
c     fsect the smallest fractional distance along 1->2 of an intersection.
c     psect the wall fractional index of this intersection.
c     fsect, psect =-1 if w2sect=0. 

c  Find the intersection of the line joining the two points and the
c  wall. Return the number of intersections w2sect. If w2sect.ne.0
c  find the fractional distance between rz1, rz2 of the intersection
c  closest to r1,z1: fsect, and the wall index to which the intersection
c  closest to r1,z1 corresponds, including the fractional distance along
c  the intersecting section, psect. (Otherwise they are both -1).

      real rwall(nwall),zwall(nwall)
c     
c     Coefficients of equation of line connecting two points:
      a=+(z1-z2)
      b=-(r1-r2)
      c= r1*z2-r2*z1
c     
c     For each wall segment, determine if it intersects the
c     line. The test for intersection is based on the fact that the
c     endpoints of each segment must lie on opposite sides of the equation
c     describing the line of the other segment.
c     
      icount=0
      psect=-1
      fsect=-1
c a large number
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
c     Prescribe that the second end of the line R2,Z2 does not count if it is
c     exactly on the boundary but the first end always does. 
            if (((val3.le.0.).and.(val4.gt.0)).or.
     &           ((val3.ge.0).and.(val4.lt.0))) then
c     section for intersecting wall sections
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
c
c*****************************************************************
      integer function icross_geom(x1,x2,i)
c Return the number of intersections that the straight line from x1 to x2
c makes with the object i. Return 0 for a non-existent object.
      integer i
c Common object geometric data.
      include 'ndimsdecl.f'
      real x1(ndims),x2(ndims)

      include '3dcom.f'
      real xn1(ndims),xn2(ndims)
      external interp

      icross_geom=0
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)
c Start of type choices
      if(itype.eq.1)then
c------------------------------------------------
c Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims)
         call sphereinterp(ndims,0,x1,x2,obj_geom(ocenter,i)
     $        ,obj_geom(oradius,i),f1,f2,sd,C,D)
         if(sd.eq.0)return
c No intersection.
         if(f1.gt.1. .or. f1.lt.0.)return
c Nearest intersection is beyond the segment.
         if(f2.gt.1. .or. f2.lt.0.)then
            icross_geom=1
         else
            icross_geom=2
         endif
      elseif(itype.eq.2)then
c------------------------------------------------
c Coordinate-Aligned Cuboid:
c Normalize
         do k=1,ndims
            xn1(k)=(x1(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i)
            xn2(k)=(x2(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i)
         enddo
         call cubenormsect(xn1,xn2,icross_geom)
      elseif(itype.eq.3)then
c------------------------------------------------
c Coordinate-Aligned Cylinder data:  Center(ndims), Semi-axes(ndims), 
c Axial coordinate. (Signed Axial 1/2 length 
c = semi-axis of the axial coordinate).
c Normalize
         ic=int(obj_geom(ocylaxis,i))
         do k=1,ndims
            kn=mod(2+k-ic,3)+1
            xn1(kn)=(x1(k)-obj_geom(ocenter+k-1,i))
     $           /obj_geom(oradius-1+k,i)
            xn2(kn)=(x2(k)-obj_geom(ocenter+k-1,i))
     $           /obj_geom(oradius-1+k,i)
         enddo
         call cylnormsect(xn1,xn2,icross_geom)
      elseif(itype.eq.4)then
c------------------------------------------------
c General Cuboid is equivalent to General Parallelopiped
c "center(ndims)", vectors(ndims,ndims), contravariants.
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
c This may be faster than the above calls.
c         do k=1,ndims
c            xn1(k)=0.
c            xn2(k)=0.
c            do j=1,ndims
c               xn1(k)=xn1(k)+(x1(j)-obj_geom(ocenter+j-1,i))
c     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
c               xn2(k)=xn2(k)+(x2(j)-obj_geom(ocenter+j-1,i))
c     $              *obj_geom(ocontra+ndims*(k-1)+j-1,i)
c            enddo
c         enddo
         call cubenormsect(xn1,xn2,icross_geom)
      elseif(itype.eq.5)then
c------------------------------------------------
c Non-aligned cylinder
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
c If need be, use the explicit code from previous case for speed.
         call cylnormsect(xn1,xn2,icross_geom)
      elseif(itype.eq.6.or.itype.eq.7)then
c------------------------------------------------
c Surface of revolution. General. Scale to contravariant components.
c Then use svrnormsect.
         call world3contra(ndims,x1,xn1,i)
         call world3contra(ndims,x2,xn2,i)
         call srvnormsect(xn1,xn2,i,icross_geom)
      endif

      end
C*******************************************************************
      subroutine cubenormsect(xn1,xn2,icross_geom)
c Find the number of intersections of the line segment from xn1 to xn2
c with the centered normalized cube of half-side length 1.
      include 'ndimsdecl.f'
c Count the number of cube ranges each point is outside.
      k1=0
      k2=0
      do k=1,ndims
         if(abs(xn1(k)).ge.1.)k1=k1+1
         if(abs(xn2(k)).ge.1.)k2=k2+1
c Break if we've found two outside.
         if(k1.ne.0.and.k2.ne.0)goto 21
      enddo
c Not both are outside
      if(k1.eq.0)then
         if(k2.eq.0)then
c 1,2 both inside
            icross_geom=0
         else
c 1 inside 2 outside
            icross_geom=1
         endif
      else
c 2 inside 1 outside
         icross_geom=1
      endif
      return
 21   continue
c 1,2 both outside. Zero or 2 crossings. Further tests needed.
c Find intersection(s) with bounding planes. If any lies inside
c the square defined by the orthogonal planes, there's an intersection.
c In that case, there are two intersections total. 
      do k=1,ndims
         xd=xn2(k)-xn1(k)
         if(xd.ne.0)then
c It is possible the line intersects planes of constant k-coordinate.
            do j=-1,1,2
               f=(j-xn1(k))/xd
               if(f.ge.0.and.f.le.1.)then
c Intersects this plane. Test whether other dimensions are inside
                  do i=mod(k,ndims)+1,mod(k+ndims-2,ndims)+1
                     xi=xn1(i)+f*(xn2(i)-xn1(i))
                     if(abs(xi).gt.1)goto 22
                  enddo
c Intersection is inside all other planes. So there is >=1, hence 2.
                  icross_geom=2
                  return
               endif
c Intersection is not inside all others. Invalid.
 22            continue
            enddo
         endif
      enddo
c Failed to find an intersection of two outside points.
      icross_geom=0

      end
c*********************************************************************
      subroutine cylnormsect(xn1,xn2,icross_geom)
c Find the number of intersections of the line segment joining xn1,xn2,
c with the unit (half-length and radius) cylinder.
      include 'ndimsdecl.f'
      integer icross_geom
      real xn1(ndims),xn2(ndims)
      real zero(ndims),one(ndims)
      data zero/ndims*0./one/ndims*1/

c On entry xn1,2 are the coordinates relative to the unit cylinder.
      z1=xn1(ndims)
      z2=xn2(ndims)
      r1=sqrt(xn1(1)**2+xn1(2)**2)
      r2=sqrt(xn2(1)**2+xn2(2)**2)
      k1=0
      if(abs(z1).gt.1.)k1=k1+1
      if(abs(r1).gt.1.)k1=k1+1
      k2=0
      if(abs(z2).gt.1.)k2=k2+1
      if(abs(r2).gt.1.)k2=k2+1
      if(k1.eq.0.or.k2.eq.0)then
c Not both are outside
         if(k1.eq.0)then
            if(k2.eq.0)then
c 1,2 both inside
               icross_geom=0
            else
c 1 inside 2 outside
               icross_geom=1
            endif
         else
c 2 inside 1 outside
            icross_geom=1
         endif
      else
c Both outside. Need more interpolation.
c z-faces:         
         ida=3
         xd=xn2(ndims)-xn1(ndims)
         if(xd.ne.0)xd=1.e-20
         do j=-1,1,2
            f=(j-xn1(ndims))/xd
            if(f.ge.0.and.f.le.1.)then
c Intersects this plane. Test whether r at intersect is inside
               r2=0.
               do i=mod(ida,ndims)+1,mod(ida+ndims-2,ndims)+1
                  r2=r2+(xn1(i)+f*(xn2(i)-xn1(i)))**2
               enddo
               if(r2.gt.1.)goto 32
c Intersection is inside circle. So there is >=1, hence 2.
               icross_geom=2
               return
            endif
 32         continue
         enddo
c No z-face intersections. Find the intersection (if any) with the
c curved surface by projecting along ida and testing a circle
c of center 0 and radius 1.
         call sphereinterp(ndims,ida,xn1,xn2,zero,one,f1,f2,sd,C,D)
         if(sd.ne.0)then
c The projected chord intersects the circle. At relevant z?
c Only need to test one intersection, because 1=>2, since no ends.
            zx1=xn1(ndims)+f1*(xn2(ndims)-xn1(ndims))
            if(zx1.lt.1. .and. zx1.gt.0.)then
               icross_geom=2
            else
c Intersections with circle are off the z-ends. 
               icross_geom=0
            endif
         else
            icross_geom=0
         endif
      endif
      end
c*********************************************************************
      subroutine srvnormsect(xp1,xp2,iobj,icross_geom)
c Given a general SoR object iobj. Count the number of intersections of
c the line joining points xn1,xn2 (contravariant components) with
c it. This version uses r-scaling, not iteration.

      integer iobj
c Local storage
      include 'ndimsdecl.f'
      real xp1(ndims),xp2(ndims)
      include '3dcom.f'
      real xn1(ndims),xn2(ndims)
      equivalence (xn1(ndims),z1),(xn2(ndims),z2)
      real zero(ndims),one(ndims)
      data zero/ndims*0./one/ndims*1/

      icross_geom=0
      zp1=xp1(ndims)
      rp1=sqrt(xn1(1)**2+xn1(2)**2)
      zp2=xp2(ndims)
      rp2=sqrt(xn2(1)**2+xn2(2)**2)

c For each SoR-segment find if the line xp1-xp2 intersects it.
c We scale the radius proportional to the end radii of the SoR-segment.
c Then the scaled test is essentially a straight cylinder intersection.

      do i=1,int(obj_geom(onpair,iobj))-1
         rw1=obj_geom(opr+i-1,iobj)
         zw1=obj_geom(opz+i-1,iobj)
         rw2=obj_geom(opr+i,iobj)
         zw2=obj_geom(opz+i,iobj)
         zd=(zw2-zw1)/2.
         zc=(zw2+zw1)/2.
         rc=(rw2+rw1)/2.
         if(abs(zd).gt.1.e-8*rc)then
c Standard surface can be scaled to unit cylinder. 
c Get point positions relative to unit cylinder.
            z1=(zp1-zc)/zd
            z2=(zp2-zc)/zd
            if(z1.gt.1. .and. z2.gt.1. .or.
     $           z1.lt.-1. .and. z2.lt.-1.)then
c Non intersecting: points both off the same end of segment.
c Hopefully this case is the dominant one and quick. Continue to next.
               goto 1
            else
c Scale points to linear cone segment ends
               rd=(rw2-rw1)/2.
               rs1=rc+z1*rd
               if(abs(rs1).lt.1.e20)rs1=1.e20
               xn1(1)=xp1(1)/rs1
               xn1(2)=xp1(2)/rs1
               rs2=rc+z2*rd
               if(abs(rs2).lt.1.e20)rs2=1.e20
               xn2(1)=xp2(1)/rs2
               xn2(2)=xp2(2)/rs2
c Now find intersections with straight unit radius cylinder by projection. 
c Only the round surface matters.
               ida=3
               call sphereinterp(ndims,ida,xn1,xn2,zero,one,f1,f2,sd,C
     $              ,D)
               if(sd.ne.0)then
                  zx1=z1+f1*(z2-z1)
                  if(zx1.le.1. .and. zx1.gt.0.)then
                     icross_geom=icross_geom+1
                  endif
                  zx2=z1+f2*(z2-z1)
                  if(zx2.le.1. .and. zx2.gt.0.)then
                     icross_geom=icross_geom+1
                  endif
               endif
            endif
         else
c Degenerate Plane surface with identical z. Alternate treatment, unscaled.
            f=(zw1-zp1)/(zp2-zp1)
            if(f.gt.0. .and. f.lt.1.)then
               r=sqrt((xp1(1)+f*(xp2(1)-xp1(1)))**2
     $              +(xp1(2)+f*(xp2(2)-xp1(2)))**2)
               if(r.gt.min(rw1,rw2) .and. r.le.max(rw1,rw2))then
                  icross_geom=icross_geom+1
               endif
            endif
         endif
 1       continue
      enddo
      end
C*********************************************************************
