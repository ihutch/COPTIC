c*****************************************************************
c Initialize with zero 3d objects.
      block data com3dset
      include '3dcom.f'
      include 'meshcom.f'
      data ngeomobj/0/
c Default track no objects.
      data nf_map/ngeomobjmax*0/
c And no reverse-map pointers.
      data nf_geommap/nf_obj*0/
      integer ibm
c Default particle region: outside object 1, inside 2.
c      parameter (ibm=ibtotal_part-4)
c      data ibool_part/1,-1,1,2,ibm*0/
c Default particle region: zero boolean. Particles everywhere.
      parameter (ibm=ibtotal_part)
      data ibool_part/ibm*0/
      integer ifm1
c Set all the bits of ifield_mask: =2**31-1=2*(2**30-1)+1 avoid overflow.
      parameter (ifm1=2*(2**30-1)+1)
c      parameter (ifm1=2**10-1)
c      parameter (ifm1=31)
      data ifield_mask/ifm1/

c Normally there's no external field.
      data lextfield/.false./extfield/ns_ndims*0./

c Mesh default initialization (meshcom.f)
      parameter (imsr=ndims_mesh*(nspec_mesh-2))
      data imeshstep/ndims_mesh*1,ndims_mesh*32,imsr*0/
      data xmeshpos/ndims_mesh*-5.,ndims_mesh*5.,imsr*0./

c We don't do flux initialization in a block data. Too big.
      end
c**********************************************************************
      subroutine geomdocument()
      write(*,*)'Format and meaning of the object geometry file'
     $     ,' default [ ccpicgeom.dat'
      write(*,*)'First line number of dimensions: 3'
      write(*,*)'Thereafter ignored comment lines start with #'
      write(*,*)'Object lines have the format: type, a,b,c, center(3),'
     $     ,' radii(3), [extra data]'
      write(*,*)' if a=b=c=0, no potential BC is applied and object'
     $     ,' is masked.'
      write(*,*)'type indicates how to use the line. Higher bytes'
     $     ,' of type indicate specials.'
      write(*,*)'byte-1: 1 Spheroid, 2 Cuboid, 3 Cylinder, 4 Parallel'
     $     ,'opiped,'
      write(*,*)' 99 Boolean region,  91-3 Set mesh in dimension 1-3. '
      write(*,*)'byte-2: 1(x256) Special boundary phi=0 instead of '
     $     ,'continuity.'
      write(*,*)'byte-2: 2(x256) Tally exit?'
      write(*,'(a)')
     $     'Boolean particle region 99, n1, n1*values, n2, values,.. 0:'
     $     ,'The nk values are ored together. The groups are anded.'
     $     ,'value n means inside object n, -n outside, 0 do nothing.'
      write(*,'(2a)')
     $     ' E.g. -1: outside object 1, and 2 inside object 2:'
     $     ,' 99, 1, -1, 1, 2, 0'
     $     ,' or: (inside 3 OR inside 4) AND outside 5 [ors first]:'
     $     ,' 99, 2, 3, 4, 1, -5, 0'
      write(*,'(a)')'Mesh setting: 9d,is1,is2,...,isN,0,xm1,xm2...xmN'
      write(*,'(2a)')'Set structure for dimension d. N-1 blocks.'
     $     ,' Steps between is1,is2 '
      write(*,'(2a)')'equally spaced between xm1,xm2: x(is1)=xm1;'
     $     ,'x(is2)=xm2, etc.'
      write(*,*)'E.g. 92,1,12,20,32,0;-5.,-1.,1.,5. y has 3 domains: '
     $     ,'1-12 covers -5. to -1.;',' 12-20 covers -1. to 1.;'
     $     ,'20-32 covers 1. to 5.'
      write(*,*)'Default equivalent to 9d,1,32,0,-5.,5.'


      end

c**********************************************************************
      subroutine readgeom(filename,myid)
c Read the geometric data about objects from the file filename
      character*(*) filename
      character*128 cline
c      character*128 fstring
      include '3dcom.f'
      include 'meshcom.f'
c Common data containing the object geometric information. 
c Each object, i < 64 has: type, data(odata).
c      integer ngeomobjmax,odata,ngeomobj
c      parameter (ngeomobjmax=31,odata=16)
c      real obj_geom(odata,ngeomobjmax)
c      common /objgeomcom/ngeomobj,obj_geom
      intrinsic IBCLR
      real valread(2*nspec_mesh)

c Zero the obj_geom data.
      do j=1,odata
         do i=1,ngeomobjmax 
            obj_geom(odata,ngeomobjmax)=0.
         enddo
      enddo
c Read
      open(1,file=filename,status='old',err=101)
      iline=1
c First line must be the number of dimensions.
      read(1,'(a)',end=902)cline
      read(cline,*,err=901)nd

      if(myid.eq.0)write(*,'(a,a50)')'Reading objects from file: '
     $     ,filename
      if(myid.eq.0)write(*,*)'Object Descr   type',
     $     '       (BCs)              (center)              (radii)'
c Loop over lines of the input file.
 1    iline=iline+1
      read(1,'(a)',end=902)cline
      if(cline(1:1).eq.'#') goto 1
      if(cline(1:6).eq.'      ') goto 1

      read(cline,*,err=901)type
c Use only lower byte.
      itype=int(type)
      type=int(itype - 256*(itype/256))
      ngeomobj=ngeomobj+1
      if(type.eq.1.)then
         read(cline,*,err=901,end=801)
     $        (obj_geom(k,ngeomobj),k=1,odata)
 801     if(myid.eq.0)write(*,820)ngeomobj,' Spheroid '
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
      elseif(type.eq.2.)then
         read(cline,*,err=901,end=802)
     $        (obj_geom(k,ngeomobj),k=1,odata)
 802     if(myid.eq.0)write(*,820)ngeomobj,' Cuboid '
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
      elseif(type.eq.3.)then
         read(cline,*,err=901,end=803)
     $        (obj_geom(k,ngeomobj),k=1,odata)
 803     if(myid.eq.0)write(*,820)ngeomobj,' Cylinder '
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
      elseif(type.eq.4.)then
         read(cline,*,err=901,end=804)
     $        (obj_geom(k,ngeomobj),k=1,odata)
 804     if(myid.eq.0)write(*,820)ngeomobj,
     $        ' General Cuboid/Parallelepiped '
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),
     $        k=1,1+nd*(1+nd)+3)
      elseif(type.eq.99)then
c Specify the particle region.
         read(cline,*,err=901,end=899)idumtype,ibool_part
 899     if(myid.eq.0)write(*,898)
     $        ngeomobj,idumtype,(ibool_part(i),i=1,16)
 898     format(i3,' Boolean ',17i4)
c Don't count this as an object.
         ngeomobj=ngeomobj-1
         goto 1
      elseif(type.gt.90.and.type.le.90+ndims_mesh)then
c Mesh specification for a particular dimension. 
         id=int(type-90)
         ist=0
         read(cline,*,err=901,end=880)idumtype,valread
 880     do i=1,nspec_mesh
            ist=i
            imeshstep(id,i)=int(valread(i))
            if(valread(i).eq.0.)goto 881
         enddo
 881     do j=1,ist-1
            xmeshpos(id,j)=valread(ist+j)
            if(j.gt.1)then
               if(xmeshpos(id,j).le.xmeshpos(id,j-1))then
                  write(*,*)'Readgeom: Meshpos not monotonic'
     $                 ,j,xmeshpos
                  stop
               endif
            endif
         enddo
c Don't count this as an object.
         ngeomobj=ngeomobj-1
         goto 1
      endif
c If this is a null boundary condition clear the relevant bit.
      if(obj_geom(oabc,ngeomobj).eq.0.
     $     .and. obj_geom(oabc+1,ngeomobj).eq.0.
     $     .and. obj_geom(oabc+2,ngeomobj).eq.0.)
     $     ifield_mask=IBCLR(ifield_mask,ngeomobj-1)
 820  format(i3,a,$)
 821  format(f4.0,9f7.3)
      goto 1

 901  write(*,*)'Readgeom error reading line',iline,':'
      write(*,*)cline
      
 902  continue

      return

 101  write(*,*) 'Readgeom File ',filename,' could not be opened.'
      stop

      end
c****************************************************************
      subroutine spheresect(id,ipm,ndims,indi,rc,xc,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c the ndims-dimensional spheroid of semi-radii rc(ndims), center xc.
c
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real xc(ndims),rc(ndims)
      real fraction,dp
      include 'meshcom.f'
      A=0.
      B=0.
      C=-1.
      D=-1.
c x1 and x2 are the coordinates in system in which sphere has radius 1.
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=(xn(ix1)-xc(i))/rc(i)
         x2=x1
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=(xn(ix2)-xc(i))/rc(i)
            A=A+(x2-x1)**2
            B=B+x1*(x2-x1)
            dp=abs(x2-x1)*rc(i)
         endif
         C=C+x1**2
         D=D+x2**2
      enddo
      fraction=1.
c This condition tests for a sphere crossing.
      if(D.ne.0. .and. D*C.le.0.)then
         if(B.ge.0. .and. A*C.le.0) then
            fraction=(-B+sqrt(B*B-A*C))/A
         elseif(B.lt.0. .and. A*C.ge.0)then
            fraction=(-B-sqrt(B*B-A*C))/A
         endif
c That should exhaust the possibilities.
      endif
      end
c****************************************************************
      function insideall(ndims,x)
c For an ndims-dimensional point x, return the integer insideall
c consisting of bits i=0-30 that are zero or one according to whether
c the point x is outside or inside object i.
      integer ndims
      real x(ndims)
c Common object geometric data.
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
      function insidemask(ndims,x)
c For an ndims-dimensional point x, return the integer insidemask
c consisting of bits i=1-31 that are zero or one according to whether
c the point x is outside or inside object i. But bits are masked by
c the regionmask.
      integer ndims
      real x(ndims)
      intrinsic btest
c Common object geometric data.
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
c Return the masked iregion.
      function imaskregion(iregion)
      include '3dcom.f'
      imaskregion=IAND(iregion,ifield_mask)
      end
c*****************************************************************
      function inside_geom(ndims,x,i)
c Return integer 0 or 1 according to whether ndims-dimensional point x
c is outside or inside object number i. Return 0 for non-existent object.
      integer ndims,i
      real x(ndims)

c Common object geometric data.
      include '3dcom.f'
      
      inside_geom=0
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)
      if(itype.eq.0)then
         return

      elseif(itype.eq.1)then
c Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims) 
         r2=0
         do k=1,ndims
            r2=r2+((x(k)-obj_geom(ocenter-1+k,i))/
     $           obj_geom(oradius-1+k,i))**2
         enddo
         if(r2.lt.1.) inside_geom=1

      elseif(itype.eq.2)then
c Coordinate-Aligned Cuboid data: low-corner(ndims), high-corner(ndims)
         do k=1,ndims
            xk=x(k)
            xl=obj_geom(ocenter-1+k,i)
            xh=obj_geom(oradius-1+k,i)
            if((xk-xl)*(xh-xk).lt.0) return
         enddo
         inside_geom=1

      elseif(itype.eq.3)then
c Coordinate-Aligned Cylinder data:  Face center(ndims), 
c Semi-axes(ndims), Axial coordinate, Signed Axial length.
         ic=int(obj_geom(ocylaxis+1,i))
         xa=(x(ic)-obj_geom(1+ic,i))
         if(xa*(obj_geom(ocylaxis+2,i)-xa).lt.0.) return
         r2=0.
         do k=1,ndims
            if(k.ne.ic)
     $         r2=r2+((x(k)-obj_geom(ocenter-1+k,i))
     $           /obj_geom(oradius-1+k,i))**2
         enddo
c         write(*,*)'Cyl. ic=',ic,' r2=',r2,' x=',x
         if(r2.lt.1.) inside_geom=1
      elseif(itype.eq.4)then
c General Cuboid data: Origin corner(ndims), vectors(ndims,ndims) 
c to adjacent corners. Is equivalent to 
c General Parallelepiped, where vectors are the face normals of
c length equal to the distance to the opposite face.
         do k=1,ndims
            proj=0.
            plen=0
            do j=1,ndims
               proj=proj+(x(j)-obj_geom(ocenter-1+j,i))
     $              *obj_geom(oradius-1+ndims*(k-1)+j,i)
               plen=plen+obj_geom(oradius+ndims*(k-1)+j,i)**2
            enddo
            if(proj.gt.plen)return
            if(proj.lt.0.)return
         enddo
         inside_geom=1
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
c Specific routine for this problem.
      subroutine potlsect(id,ipm,ndims,indi,fraction,conditions,dp,
     $     iobjno)

c In dimension id, direction ipm, 
c from mesh point at indi(ndims) (zero-based indices, C-style),
c find any intersection of the mesh leg from this point to its neighbor
c with a bounding surface. Return the "fraction" of the leg at which
c the intersection occurs (1 if no intersection), the "conditions" at
c the intersection (irrelevant if fraction=1), the +ve length
c in computational units of the full leg in dp, and the object number
c in iobjno
      integer id,ipm,ndims,iobjno
      integer indi(ndims)
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
c Continuity alone should be applied
c on the other (inactive) side of the surface.  A negative value of c2
c is used to indicate that just the appropriate continuity condition is
c to be applied, but c1,c2,c3, should all be returned with the magnitudes
c that are applied to the active side of the surface, with consistently
c reversed signs. (e.g. a=1, b=1, c=-2 -> a=-1, b=-1, c=2)
c
c A fraction of 1 causes all the bounding conditions to be ignored.

      include 'meshcom.f'
      include '3dcom.f'

      real xx(10),xd(10)

c Default no intersection.
      fraction=1
      iobjno=0

c-------------------------------------------------------------
c Process data stored in obj_geom.
      do i=1,ngeomobj
         if(obj_geom(oabc,i).ne.0. .or. obj_geom(oabc+1,i).ne.0. .or.
     $        obj_geom(oabc+2,i).ne.0.)then
c Only for non-null BCs
c Find the fractional intersection point if any.
c Currently implemented just for spheres.
            call spheresect(id,ipm,ndims,indi,
     $        obj_geom(oradius,i),obj_geom(ocenter,i)
     $        ,fraction,dp)
         else
c Null BC. Ignore.
            fraction=1.
         endif
         if(fraction.ne.1.)then
            if(obj_geom(oabc+1,i).eq.0)then
c No derivative term. Fixed Potential. No projection needed.
               conditions(1)=obj_geom(oabc,i)
               conditions(2)=obj_geom(oabc+1,i)
               conditions(3)=obj_geom(oabc+2,i)
            else
c Derivative term present. Calculate the projection:
c    Calculate Coordinates of crossing
               do idi=1,ndims
c    Address of mesh point. 
                  ix=indi(idi)+ixnp(idi)+1
                  xx(idi)=xn(ix)
                  if(idi.eq.idi)
     $                 xx(idi)=xx(idi)+fraction*(xn(ix+ipm)-xn(ix))
c    Component of outward vector = Sphere outward normal.
                  xd(idi)=xx(idi)-obj_geom(ocenter+idi-1,i)
               enddo
c    Signed projection cosine:
               projection=ipm*(xd(id)/obj_geom(oradius,i))
               if(projection.eq.0.)then
                  fraction=1.
               else
c    Continuity works for both directions.
                  conditions(1)=sign(obj_geom(oabc,i),projection)
                  conditions(2)=obj_geom(oabc+1,i)/projection
                  conditions(3)=sign(obj_geom(oabc+2,i),projection)
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
c*************************************************************
c Initialize the iregion flags of the existing nodes with boundary
c object data.
      subroutine iregioninit(ndims,ifull)
      integer ifull(ndims)

      include 'objcom.f'
      include 'meshcom.f'
      include '3dcom.f'

      integer ix(ndims_sor)
      real x(ndims_sor)

      if(ndims.ne.ndims_sor)then 
         write(*,*)'iregioninit error; incorrect dimensions:',
     $        ndims,ndims_sor
         call exit(0)
      endif

      do i=1,oi_sor
         ipoint=idob_sor(ipoint_sor,i)
c Convert index to multidimensional indices.
         call indexexpand(ndims,ifull,ipoint,ix)
         do k=1,ndims
            x(k)=xn(ixnp(k)+ix(k))
         enddo
c Store in object-data.
c         idob_sor(iregion_sor,i)=insideall(ndims,x)
         idob_sor(iregion_sor,i)=insidemask(ndims,x)

c         write(*,101)i,ipoint,idob_sor(iregion_sor,i),x
 101     format(3i8,5f10.4)
      enddo

      end
c*****************************************************************
      subroutine reportfieldmask()
      include '3dcom.f'
      integer ip(32)
c Calculate the bits of the field mask.
      ifd=ifield_mask
      do i=1,32
         ip(32-i+1)=ifd - 2*(ifd/2)
         ifd=ifd/2
      enddo
c      write(*,*)'Initializing Object Regions:No, pointer, region, posn.'
c This is an unportable extension. Hence the calculation above.
c      write(*,'('' Mask='',i11,'' ='',b32.32)')ifield_mask,ifield_mask
      write(*,'('' Field Mask='',i11,'' ='',32i1)')ifield_mask,ip
      end
c*******************************************************************
      function ireg3(i,j,k,ifull,cij)
      include 'objcom.f'
      integer ifull(3)
      real cij(ndims_sor*2+1,ifull(1),ifull(2),ifull(3))

      ipoint=int(cij(ndims_sor*2+1,i,j,k))
      if(ipoint.ne.0)then
         ireg3=idob_sor(iregion_sor,ipoint)
      else
         ireg3=99
      endif

      end
c*******************************************************************
      subroutine objsetabc(iobject,a,b,c)
c Set/Reset the boundary conditions for an object that has already been
c read in. 
      integer iobject
      real a,b,c

      include '3dcom.f'

      if(iobject.gt.ngeomobj)then
         write(*,'(a,i4,a,i4)')
     $        '**** objsetabc error. Attempt to set object',iobject
     $        ,' greater than existing',ngeomobj
         return
      endif
      obj_geom(oabc,iobject)=a
      obj_geom(oabc+1,iobject)=b
      obj_geom(oabc+2,iobject)=c
c      write(*,'(''Set obj_geom(oabc,'',i2,'')='',3f8.4)')
c     $     iobject,(obj_geom((oabc+i),iobject),i=0,2)
      end
