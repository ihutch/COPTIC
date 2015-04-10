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
      include 'myidcom.f'
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
      real xp1(ndims),xp2(ndims),xn1(ndims),xn2(ndims)
      real fmin(ovlen)
      integer ids(ovlen)
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
c First implemented just for spheres
               ida=0
               call spheresect(ndims,ida,xp1,xp2,obj_geom(ocenter,i)
     $              ,obj_geom(oradius,i),fraction ,f2,sd,C,D)
               if(sd.eq.0.or.fraction-1..gt.0. .or. fraction.lt.0.)then
                  fraction=1.
               else
                  call ijbinsphere(i,fraction,xp1,xp2,ijbin)
               endif
            elseif(itype.eq.2)then
c Coordinate-aligned cuboid.               
c Convert into normalized position for object i.
               do j=1,ndims
                  xn1(j)=(xp1(j)-obj_geom(ocenter+j-1,i))
     $                 /obj_geom(oradius+j-1,i)
                  xn2(j)=(xp2(j)-obj_geom(ocenter+j-1,i))
     $                 /obj_geom(oradius+j-1,i)
               enddo
               call cubeusect(xn1,xn2,nsect,fmin,ids)
               call ijbincube(i,ids(1),fmin(1),xn1,xn2,ijbin,idebug)
               fraction=fmin(1)
            elseif(itype.eq.3)then
c Coordinate aligned cylinder.
c Convert to normalized.
               ida=int(obj_geom(ocylaxis,i))
               do j=1,ndims
                  ii=mod(j-ida+2,ndims)+1
                  xn1(ii)=(xp1(j)-obj_geom(ocenter+j-1,i))
     $                 /obj_geom(oradius+j-1,i)
                  xn2(ii)=(xp2(j)-obj_geom(ocenter+j-1,i))
     $                 /obj_geom(oradius+j-1,i)
               enddo
               call cylusect(xn1,xn2,i,nsect,fmin,ids)
               fraction=fmin(1)
               if(nsect.ge.1)then
                  call ijbincyl(i,ids(1),fmin(1),xn1,xn2,ijbin)
               endif
            elseif(itype.eq.4)then
c Parallelopiped.
               call xp2contra(i,xp1,xp2,xn1,xn2,ins1,ins2)
               call cubeusect(xn1,xn2,nsect,fmin,ids)
               fraction=fmin(1)
               if(nsect.ge.1)then
                  call ijbincube(i,ids(1),fmin(1),xn1,xn2,ijbin,idebug)
               endif
            elseif(itype.eq.5)then
c Non-aligned cylinder.
               call xp2contra(i,xp1,xp2,xn1,xn2,ins1,ins2)         
               call cylusect(xn1,xn2,i,nsect,fmin,ids)
               fraction=fmin(1)
               if(nsect.ge.1)then
                  call ijbincyl(i,ids(1),fmin(1),xn1,xn2,ijbin)
               endif
            elseif(itype.eq.6.or.itype.eq.7)then
c Surface of revolution.
               call xp2contra(i,xp1,xp2,xn1,xn2,ins1,ins2)
               call srvsect(xn1,xn2,i,nsect,fmin,ids)
               if(nsect.ge.1)then
                  call ijbinsrv(i,ids(1),fmin(1),xp1,xp2,ijbin)
               endif
               fraction=fmin(1)
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
            if(myid.eq.0.and.numreports.lt.1)then
               write(*,'(a,i2,a,i3,a,20i2)')
     $              ' Subtraction from obj',i,'; total',normv(i)
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
c*********************************************************************
c*********************************************************************
      subroutine xp2contra(iobj,xp1,xp2,xn1,xn2,ins1,ins2)
c Transform from world (xp1,xp2) to normalized contravariant (xn1,xn2)
c wrt object iobj. Assumed storage not overlapping.
c Return ins1/2 whether points 1,2 are inside (1) unit cube or outside (0).

c The object is specified by center and three vectors pqr. Each of nsdim
c pairs of parallel planes consists of the points: +-p + c_q q + c_r r.
c Where p is one of the three base (covariant) vectors and qr the
c others, and c_q are real coefficients.  Inside corresponds to between
c these two planes, i.e. contravariant coefficients <1. We define the
c facets of the cube to be the faces (planes) in the following order:
c +v_1,+v_2,+v_3,-v_1,-v_2,-v_3. 
c Then within each face the facet indices are in cyclic order. But that
c is determined by the cubeusect code.
      implicit none
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iobj
      real xp1(3),xp2(3),xn1(3),xn2(3)
      integer ins1,ins2

      integer j,i,ii,ji
      real xc
      ins1=1
      ins2=1
      do j=1,ndims
         xn1(j)=0.
         xn2(j)=0.
c Contravariant projections.
         do i=1,ndims
c Cartesian coordinates.
            ii=(ocenter+i-1)
            xc=obj_geom(ii,iobj)
c xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(ocontra+ndims*(j-1)+i-1)
c            write(*,*)'ji',ji
            xn1(j)=xn1(j)+(xp1(i)-xc)*obj_geom(ji,iobj)
            xn2(j)=xn2(j)+(xp2(i)-xc)*obj_geom(ji,iobj)
         enddo
         if(abs(xn1(j)).ge.1.)ins1=0
         if(abs(xn2(j)).ge.1.)ins2=0
      enddo
      end
c*********************************************************************
      subroutine srvfsect(nsdim,xp1,xp2,iobj,ijbin,sd,fraction)
c Given a general SoR object iobj. Find the point of intersection of the
c line joining xp1,xp2, with it, and determine the ijbin to which it is
c therefore assigned, and the direction it is crossed (sd=+1 means
c inward from 1 to 2). Fractional distance xp1->xp2 is returned.
c xp1,xp2 here are world, not contra. Have to be transformed for comparison.

c This version uses only inside_geom and bisection to find the point
c of intersection. It is not now used.

      integer nsdim,iobj,ijbin
      real xp1(nsdim),xp2(nsdim),sd,fraction
c Local storage
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xm(ndims),xc1(ndims),xc2(ndims),rm

      ins1=inside_geom(ndims,xp1,iobj)
      ins2=inside_geom(ndims,xp2,iobj)
      fu=1.
      fl=0.
      if(ins1.eq.ins2)then
c Endpoints on same side.
c         write(*,*)'srvfsect error: root not bracketed.',iobj,xp1,xp2
c We need some sort of check on double intersection.
c Perhaps this could be iterated:
         fraction=0.5
         do i=1,ndims
            xm(i)=xp1(i)+(xp2(i)-xp1(i))*fraction
         enddo
         ins=inside_geom(ndims,xm,iobj)
         if(ins.ne.ins1)then 
c Double crossing. Choose solution closest to fraction=0.
c            write(*,*)'Double crossing fixed.'
            fu=fraction
            fl=0.
         else
            sd=0.
            fraction=1.
            return
         endif
      endif
c Bisection, bipolar on value of inside_geom, to find fraction.
c Iterate to accuracy of 2^16. Fairly costly.
      niter=16
      do k=1,niter
         fraction=(fu+fl)*0.5
         do i=1,ndims
            xm(i)=xp1(i)+(xp2(i)-xp1(i))*fraction
         enddo
         ins=inside_geom(ndims,xm,iobj)
         if(ins.ne.ins1)then
            fu=fraction
         else
            fl=fraction
         endif
      enddo
c Now fraction should be the intersection point.

c Its theta value is based on contravariant coordinates:
      call world3contra(ndims,xm,xm,iobj)
      theta=atan2(xm(2),xm(1))
      zm=xm(ndims)
      rm=sqrt(xm(1)**2+xm(2)**2)

c Find the wall position by calling w2sect with adjacent positions.
c Construct a short chord in the same order as xp1,xp2, that brackets
c the solution already found.
      delta=.01
      do i=1,ndims
         xc1(i)=xp1(i)+(xp2(i)-xp1(i))*(fraction-delta)
         xc2(i)=xp1(i)+(xp2(i)-xp1(i))*(fraction+delta)
      enddo
c Transform to contravariant components.
      call world3contra(ndims,xc1,xc1,iobj)
      call world3contra(ndims,xc2,xc2,iobj)
c Find the r,z normalized coordinates of the positions.
      r1=sqrt(xc1(1)**2+xc1(2)**2)
      r2=sqrt(xc2(1)**2+xc2(2)**2)
      z1=xc1(ndims)
      z2=xc2(ndims)
      isect=int(w2sect(r1,z1,r2,z2,obj_geom(opr,iobj)
     $     ,obj_geom(opz,iobj),int(obj_geom(onpair,iobj)),fsect,psect))
      if(isect.eq.0)then
c This should not happen. And diagnostics should be cleared up eventually.
         write(*,*)'No intersection in srvsect',isect,psect
         write(*,*)'r1,r2,z1,z2',r1,r2,z1,z2
c         write(*,*)'onpair',int(obj_geom(onpair,iobj))
         write(*,*)fraction,fu,fl,ins1,ins2,ins,rm,zm
         write(*,*)xc1,xc2
         write(*,*)xp1,xp2
c         write(*,*)(obj_geom(opr+i,iobj),obj_geom(opz+i,iobj),i=0
c     $        ,int(obj_geom(onpair,iobj))-1)
c         write(*,*)inside_geom(ndims,xp1,iobj)
c     $        ,inside_geom(ndims,xp2,iobj)
         call rztell(xp1,iobj)
         call rztell(xp2,iobj)
         do i=1,ndims
            xm(i)=xp1(i)+(xp2(i)-xp1(i))*fraction
         enddo
         call rztell(xm,iobj)
         call world3contra(ndims,xm,xm,iobj)
         write(*,*)sqrt(xm(1)**2+xm(2)**2),xm(3)
         write(*,'(3f8.4)')(obj_geom(ovec+k-1,iobj),k=1,18)
         do i=1,ndims
            xm(i)=xp1(i)+(xp2(i)-xp1(i))*fl
         enddo
         call rztell(xm,iobj)
         stop
      endif

c Calculate ijbin.
c Decide the irz index based upon wall position. The face is the integer
c part of psect. The facet is based upon equal r^2 division of face.
c See objplot for the principles.
c The face:
      irz=int(psect)
c Encapsulated world units.
      call ijbinsrv(iobj,irz,fraction,xp1,xp2,ijbin)

      if(ins1.eq.0)then
         sd=1.
      else
         sd=-1
      endif
      end
c*********************************************************************
      subroutine ijbinsrv(iobj,ids,fraction,xp1,xp2,ijbin)
c Given the two world (non-normalized) points, xp1,xp2, the fraction
c between them and the face index ids, of object iobj, get the facet
c index ijbin.
      implicit none
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iobj,ijbin,ids
      real fraction,xp1(ndims),xp2(ndims)

      real xm(ndims),theta,zm,rm,rb,rt,zb,zt,dr2
      integer i,nz,ifobj,itc,k,ifct

      do i=1,ndims
         xm(i)=xp1(i)+(xp2(i)-xp1(i))*fraction
      enddo
      call world3contra(ndims,xm,xm,iobj)
      theta=atan2(xm(2),xm(1))
      zm=xm(ndims)
      rm=sqrt(xm(1)**2+xm(2)**2)
      itc=int(obj_geom(ofn1,iobj)*(theta/(2.*3.141593)+0.5))

c The ends of the segment chosen.
      rb=obj_geom(opr+ids-1,iobj)
      rt=obj_geom(opr+ids,iobj)
      zb=obj_geom(opz+ids-1,iobj)
      zt=obj_geom(opz+ids,iobj)
      nz=int(obj_geom(opdiv+ids-1,iobj))
c This test sometimes fails very near the end of a segment, leading to
c benign application of a crossing to an adjacent segment. If it fails
c by a larger amount. Worry!
      if((rb-rm)*(rm-rt)*abs(rt-rb).lt.-1.e-6
     $     .or. (zb-zm)*(zm-zt)*abs(zt-zb).lt.-1.e-6)then
         write(*,*)
         write(*,*)'Puzzling Values for',ids,iobj,fraction
         write(*,*)'rb,rm,rt,zb,zm,zt',rb,rm,rt,zb,zm,zt
         write(*,'(20f7.3)')(obj_geom(opr+k,iobj),k=0
     $        ,int(obj_geom(onpair,iobj))-1)
         write(*,'(20f7.3)')(obj_geom(opz+k,iobj),k=0
     $        ,int(obj_geom(onpair,iobj))-1)
         stop
      endif
      dr2=(rt**2-rb**2)/nz
c The facet number: k3-1
      if(abs(dr2).ne.0.)then
         ifct=int((rm**2-rb**2)/dr2)
      else
         ifct=int((nz)*(zm-zb)/(zt-zb))
      endif
      if(ifct.lt.0)ifct=0
      if(ifct.gt.nz-1)ifct=nz-1
      ifobj=nf_map(iobj)
      ijbin=(itc)+nf_dimlens(nf_flux,ifobj,1)*(ifct)
     $        +nf_faceind(nf_flux,ifobj,ids)
      
      end
c*********************************************************************
      subroutine spheresect(nsdim,ida,xp1,xp2,xc,rc,f1,f2,sd,C,D)
c Given two different nsdim dimensioned vectors xp1,xp2,and a sphere
c center xc radii rc, find the intersection of the line joining x1,x2,
c with the sphere and return it as the value of the fraction f1 of
c x1->x2 to which this corresponds, chosen always positive if possible,
c and closest to 0. The other intersection fraction in f2.  Also return
c the direction of crossing in sd, and the fractional radial distance^2
c outside the sphere of the two points in C and D.  (positive means
c inward from x1 to x2). If there is no intersection even of the
c extrapolated line, return f1=1., sd=0.  If ida is non-zero then
c form the radius only over the other dimensions.  In that case the
c subsurface (circle) is the figure whose intersection is sought.
      integer nsdim,ida
      real xp1(nsdim),xp2(nsdim),xc(nsdim),rc(nsdim)
      real sd,f1,f2,C,D
c
      real x1,x2,A,B,E,disc,f3
c Prevent a singularity if the vectors coincide.
      A=1.e-25
      B=0.
      E=-4.
      C=-1.
      D=-1.
c x1 and x2 are the coordinates in system in which sphere 
c has center 0 and radius 1.
      ni=nsdim
      if(ida.ne.0)ni=nsdim-1
      do ii=1,ni
         i=mod(ida+ii-1,nsdim)+1
         xci=xc(i)
         rci=rc(i)
         x1=(xp1(i)-xci)/rci
         x2=(xp2(i)-xci)/rci
         A=A+(x2-x1)**2
         B=B+(x1+x2)*(x2-x1)
         E=E+(x1+x2)**2
         C=C+x1**2
         D=D+x2**2
      enddo
      B=B*0.5
      E=E*0.25
c Crossing direction from x1 to intersection
      sd=sign(1.,C)
c Symmetric discriminant calculation:
      disc=B*B-A*E
      if(disc.ge.0.)then
         disc=sqrt(disc)
c With different discriminant we need a different test for sign.
c Most negative f
         f1=(-B-disc)/A+0.5
         if(f1.ge.0.)then
            f2=(-B+disc)/A+0.5
         else
c f1 negative
            f3=(-B+disc)/A+0.5
            if(f3.ge.0..or.f3.le.f2)then
               f2=f1
               f1=f3
            else
               f2=f3
            endif
         endif
c Fix cases where xp1 exactly on sphere to have correct direction.
         if(f1.eq.0.and.C.lt.1.e-5.and.f2.lt.0.)sd=-sign(1.,D)
      else
c         write(*,*)'Sphere-crossing discriminant negative.',A,B,C
         sd=0.
         f1=1.
         f2=1.
c         return
      endif
      end
c***********************************************************************
      subroutine cubeusect(xp1,xp2,nsect,fmin,ids)
c For a unit cube, center 0 radii 1, find the intersection(s) of line
c xp1,xp2 with it, and the face that it intersects.
      include 'ndimsdecl.f'
      real xp1(ndims),xp2(ndims)
      include '3dcom.f'
      real fmin(ovlen)
      integer ids(ovlen)

      fmin(1)=1.
      ids(1)=0
      nsect=0
      fn=1.
      do i=1,2*ndims
         im=mod(i-1,ndims)+1
c First half of the i's are negative. Second half positive.
         xc1=float(((i-1)/ndims)*2-1)
         xd1=(xp1(im)-xc1)
         xd2=(xp2(im)-xc1)
         if(xd1.lt.0. .neqv. xd2.lt.0)then
c Crossed this plane at fraction:
            fn=xd1/(xd1-xd2)
c Intersects this plane. Test whether other dimensions are inside
            do jj=1,ndims-1
               j=mod(i+jj-1,ndims)+1
               xi=xp1(j)+fn*(xp2(j)-xp1(j))
               if(abs(xi).gt.1)goto 1
c This prevents duplicates when the intersection is at an edge:
               if(nsect.ge.1.and.fn.eq.fmin(nsect))goto 1
            enddo
c This is sorted:
            call insertsorted2(nsect,fmin,fn,ids,i)
 1          continue
         endif
c This escape prevents (unobserved) pathological third intersections.
         if(nsect.ge.2)return
      enddo
      end
c*********************************************************************
      subroutine cylusect(xp1,xp2,iobj,nsect,fmin,ids)
c Find the points of intersection of the line joining xp1,xp2, with the
c UNIT cylinder The 1-axis is where theta is measured from and the 3-axis
c is the axial direction. Return the number of intersections nsect 
c ordered in fmin and the corresponding faces ids.
c The facets of the cylinder are the end faces -xr +xr, and the curved
c side boundary. 3 altogether.  The order of faces is bottom, side, top.
      
      integer iobj
      include 'ndimsdecl.f'
      real xp1(ndims),xp2(ndims)
      include '3dcom.f'
      real fmin(ovlen)
      integer ids(ovlen)
c 3D here.
      parameter (nds=3)
c      real x12(nds)
      real fn(4),zrf(4)
      integer isdf(4)
      real xc(nds),rc(nds)
      data xc/0.,0.,0./rc/1.,1.,1./

      ida=3
      nsect=0
      fmin(1)=1.
c First shortcut, return if both points are beyond the same axial end.
      z1=xp1(ida)
      z2=xp2(ida)
      if((z1.gt.1. .and. z2.gt.1.).or.
     $     (-z1.gt.1. .and. -z2.gt.1))return
      r1=(xp1(1)**2+xp1(2)**2)+1.e-5
      r2=(xp2(1)**2+xp2(2)**2)+1.e-5
c Second shortcut, if both points are inside cylinder, return.
      if(r1.lt.1..and.r2.lt.1..and.abs(z1).lt.1..and.abs(z2).lt.1.)
     $     return
      xd=z2-z1
      sds=0.
c Find the intersection (if any) with the circular surface in the plane
c perpendicular to direction ida (projected along ida).
      call spheresect(ndims,ida,xp1,xp2, xc,rc,fn(1),fn(2),sds,d1,d2)
      if(sds.ne.0)then
c There is at least one intersection with the circle. Find where.
         isdf(1)=0
         isdf(2)=0
         zrf(1)=(1.-fn(1))*xp1(ida)+fn(1)*xp2(ida)
         zrf(2)=(1.-fn(2))*xp1(ida)+fn(2)*xp2(ida)
      else
c No radial intersections
         zrf(1)=2.
         zrf(2)=2.
c I think we should return here because if the chord does not intersect
c the cylinder the only way it can intersect the ends is if it is
c exactly parallel to the cylinder axis. spheresect prevents that happening.
         return
      endif
c Find the axial intersection fractions with the end planes
c and their radial positions.
      if(xd.ne.0)then
         fn(3)=(1.-z1)/xd
         isdf(3)=1
         fn(4)=(-1.-z1)/xd
         isdf(4)=-1
         zrf(3)=0.
         zrf(4)=0.
         do k=1,ndims
            if(k.ne.ida)then
               xkg1=(1.-fn(3))*xp1(k)+fn(3)*xp2(k)
               zrf(3)=zrf(3)+xkg1**2
               xkg2=(1.-fn(4))*xp1(k)+fn(4)*xp2(k)
               zrf(4)=zrf(4)+xkg2**2
            endif
         enddo
      else
c Pure radial difference. No end-intersections anywhere.
         zrf(3)=2.
         zrf(4)=2.
      endif
c Now we have 4 possible fractions fn(4). Up to 2 of those are
c true. Truth is determined by abs(zrf(k))<=1. and 0.<=fn<1.  The
c corresponding face numbers are in isdf.
      nsect=0
      do k=1,4
         if(abs(zrf(k)).le.1. .and. fn(k).ge.0..and.fn(k).lt.1.)then
c This prevents duplicates when the intersection is at an edge I hope:
            if(nsect.ge.1.and.fn(k).eq.fmin(nsect))then
            else
               call insertsorted2(nsect,fmin,fn(k),ids,isdf(k))
            endif
         endif
      enddo
      if(nsect.gt.2)then
         write(*,*)'cylusect too many intersections!',nsect
         write(*,*)fn
         write(*,*)zrf
         stop
      endif
c         if(fmin(1).ne.1.)write(*,*)(fmin(k),k=1,4)
      end
c*********************************************************************
      subroutine ijbinsphere(iobj,fraction,xp1,xp2,ijbin)
c Given the fractional distance between xp1 and xp2, fraction, which is the
c crossing point. Index the face position using information in
c obj_geom(ofn,iobj):
c This code decides which of the nf_posno for this object to update
c corresponding to this crossing.
      implicit none
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iobj,ijbin
      real fraction,xp1(ndims),xp2(ndims)

      real tiny,onemtiny
      parameter (tiny=1.e-5,onemtiny=1.-2.*tiny)
      real x12(ndims),psi
      integer i
      integer infobj,ibin,jbin

      infobj=nf_map(iobj)
c Calculate normalized intersection coordinates.
      do i=1,ndims
         x12(i)=((1.-fraction)*xp1(i)+fraction*xp2(i)
     $        -obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
      enddo

c Bin by cos(theta)=x12(3) uniform grid in first nf_dimension. 
c ibin runs from 0 to N-1 cos = -1 to 1.
      ibin=int(nf_dimlens(nf_flux,infobj,1)*(onemtiny*x12(3)+1.)*0.5)
      psi=atan2(x12(2),x12(1))
c jbin runs from 0 to N-1 psi = -pi to pi.
      jbin=int(nf_dimlens(nf_flux,infobj,2)
     $     *(0.999999*psi/3.1415926+1.)*0.5)
      ijbin=ibin+jbin*nf_dimlens(nf_flux,infobj,1)
      if(ijbin.gt.nf_posno(1,infobj))then
         write(*,*)'ijbin error in ijbinsphere'
         write(*,*)infobj,ijbin,nf_posno(1,infobj),ibin,jbin
     $        ,nf_dimlens(nf_flux,infobj,1),nf_dimlens(nf_flux,infobj,2)
     $        ,x12(3),obj_geom(ocenter+2,iobj),xp1(3),xp2(3)
     $        ,fraction
      endif

      end
c*******************************************************************
      subroutine ijbincube(iobj,ids,fmin,xp1,xp2,ijbin,idebug)
c Given the fractional distance between xp1 and xp2, fmin, which is the
c crossing point and ids the face-index of this crossing. Index
c within the face on equal spaced grid whose numbers have been read into
c obj_geom(ofn.,iobj):
      implicit none
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iobj,idebug,ids,ijbin
      real fmin,xp1(ndims),xp2(ndims)
      integer infobj,ibstep,ibin,k,i
      real xk,xcr
      infobj=nf_map(iobj)
      ibstep=1
      ibin=0
      if(idebug.eq.1)write(*,'(''Position'',$)')
      do i=1,ndims-1
c The following defines the order of indexation. ofn1 is the next highest
c cyclic index following the face index. So on face 1 or 4 the other
c two indices on the face are 2,3. But on face 2,5 they are 3,1.
         k=mod(mod(ids-1,3)+1+i-1,ndims)+1
         xk=(1.-fmin)*xp1(k)+fmin*xp2(k)
c This has xcr run from 0 to 1 as xk goes from -1. to +1. :
         xcr=(1.+xk)*.5
         if(idebug.eq.1)write(*,'(i2,f7.2,i5,i5,'',''$)')k,xcr
     $        ,nf_dimlens(nf_flux,infobj,k)
     $        ,int(nf_dimlens(nf_flux,infobj,k)*(0.999999*xcr))
         ibin=ibin+ibstep*
     $        int(nf_dimlens(nf_flux,infobj,k)*(0.999999*xcr))
         ibstep=ibstep*nf_dimlens(nf_flux,infobj,k)
      enddo
c Now we have ibin equal to the face-position index, and ibstep equal
c to the face-position size. Add the face-offset for ids. This is
c tricky for unequal sized faces. So we need to have stored it. 
      ijbin=ibin+nf_faceind(nf_flux,infobj,ids)
      if(idebug.eq.1)write(*,'(a,3i4)')'Ending ijbincube',ijbin
     $     ,nf_faceind(nf_flux,infobj,ids),ids
c That's it.
      end
c*********************************************************************
      subroutine ijbincyl(iobj,ids,fmin,xp1,xp2,ijbin)
c Given the fractional distance between xp1 and xp2, fmin, which is the
c crossing point and ids the face-index ids of this crossing. Index
c within the face on equal spaced grid whose numbers have been read into
c obj_geom(ofn.,iobj). Return ijbin index.
      implicit none
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer iobj,ids,ijbin
      real fmin,xp1(ndims),xp2(ndims)
      integer infobj,k,i,ir,it,iz,ida
      real x12(ndims),z,theta,r2

      ida=3
c Calculate normalized intersection coordinates.
      do i=1,ndims
         x12(i)=(1.-fmin)*xp1(i)+fmin*xp2(i)
      enddo
c Calculate r,theta,z (normalized) relative to the ida direction as z.
      z=x12(ida)
      theta=atan2(x12(mod(ida+1,ndims)+1),x12(mod(ida,ndims)+1))
      r2=0.
      do i=1,ndims-1
         k=mod(ida+i-1,ndims)+1
         r2=r2+x12(k)**2
      enddo
c      write(*,'(a,7f7.4,3i3)')'r2,theta,z,x12,fmin,ids(1)'
c     $     ,r2,theta,z,x12,fmin,ids(1)
c End blocks are of size nr x nt, and the curved is nt x nz.

      infobj=nf_map(iobj)
      ijbin=0
      if(ids.ne.0)then
c Ends
         if(ids.eq.1)then
c offset by (nr+nz)*nt
            ijbin=(nf_dimlens(nf_flux,infobj,1)
     $           +nf_dimlens(nf_flux,infobj,3))
     $           *nf_dimlens(nf_flux,infobj,2)
         endif
c Uniform mesh in r^2 normalized. 
         ir=int(nf_dimlens(nf_flux,infobj,1)*(0.999999*r2))
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         ijbin=ijbin+ir+it*nf_dimlens(nf_flux,infobj,1)
      else
c Side. Offset to this facet nr*nt:
         ijbin=nf_dimlens(nf_flux,infobj,1)*nf_dimlens(nf_flux,infobj
     $        ,2)
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         iz=int(nf_dimlens(nf_flux,infobj,3)*(0.999999*z+1.)*0.5)
c Index in order theta,z
         ijbin=ijbin+it+iz*nf_dimlens(nf_flux,infobj,2)
      endif
c      write(*,'(6f8.4,3i3)')xp1,xp2,ir,it,iz
      end
c********************************************************************
      subroutine segmentcross(x11,x12,x21,x22,f1,f2)
c Given two pairs of points in 2-D, find the intersection of the 
c line joining the first pair with that joining the second pair. 
c Return as fractions of the distances along the lines
      real x11(2),x12(2),x21(2),x22(2),f1,f2
c Distances of each point from the opposite line      
      pd21=(x12(2)-x11(2))*(x21(1)-x11(1))
     $     -(x12(1)-x11(1))*(x21(2)-x11(2))
      pd22=(x12(2)-x11(2))*(x22(1)-x11(1))
     $     -(x12(1)-x11(1))*(x22(2)-x11(2))
      pd11=(x22(2)-x21(2))*(x11(1)-x21(1))
     $     -(x22(1)-x21(1))*(x11(2)-x21(2))
      pd12=(x22(2)-x21(2))*(x12(1)-x21(1))
     $     -(x22(1)-x21(1))*(x12(2)-x21(2))
c Hence fractional intersect positions:
      f2=(0.-pd21)/(pd22-pd21)
      f1=(0.-pd11)/(pd12-pd11)
      end
c********************************************************************
      subroutine rztell(x,i)
c Common object geometric data.
      include 'ndimsdecl.f'
      include '3dcom.f'
      real x(ndims)
c Tell the fractional position in a 
c Surface of revolution. General code from inside_geom
      r=0.
      z=0.
      do j=1,ndims
         z=z+obj_geom(ovax+j-1,i)*(x(j)-obj_geom(obase+j-1,i))
      enddo
      do j=1,ndims
         r=r+( x(j)-obj_geom(obase+j-1,i)
     $        -z*(obj_geom(ovec+2*ndims+j-1,i)) )**2
      enddo
      r=sqrt(r)/obj_geom(orscale,i)
      write(*,*)'rztell',r,z,w2sect(r,z,1.e5,0.,obj_geom(opr,i)
     $     ,obj_geom(opz,i),int(obj_geom(onpair,i)),fsect,psect)
c     $     ,obj_geom(ovax+3-1,i),obj_geom(ovec+2*ndims+3-1,i)
c      write(*,'(3f8.4)')(obj_geom(ovec+k-1,i),k=1,18)
      end
c*********************************************************************
      subroutine insertsorted2(n,f,fnew,j,jnew)
c Increment n and insert the value fnew into the ordered list f(*) of n
c increasing values, and also jnew in j-list with same sort. 
      integer n
      real fnew,f(*)
      integer jnew,j(*)
      n=n+1
      do i=n,2,-1
         if(f(i-1).gt.fnew)then
            f(i)=f(i-1)
            j(i)=j(i-1)
         else
            f(i)=fnew
            j(i)=jnew
            goto 1
         endif
      enddo
      f(1)=fnew
      j(1)=jnew
 1    continue
c      write(*,*)'insertsorted',n,(f(k),j(k),k=1,n)
      end
