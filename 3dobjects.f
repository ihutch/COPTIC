c*****************************************************************
c Initialize with zero 3d objects.
      block data com3dset
      include 'ndimsdecl.f'
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
      data lextfield/.false./extfield/ndims*0./
c Mesh default initialization (meshcom.f)
      parameter (imsr=ndims*(nspec_mesh-2))
      data imeshstep/ndims*1,ndims*32,imsr*0/
      data xmeshpos/ndims*-5.,ndims*5.,imsr*0./
c Default no point charges:
      data iptch_mask/0/
c Default no subtractive objects
      data normv/ngeomobjmax*0/
c We don't do flux initialization in a block data. Too big.
      end
c**********************************************************************
      subroutine geomdocument()
      write(*,*)'######################################################'
      write(*,*)'Format and meaning of the object geometry file'
     $     ,' default [ copticgeom.dat'
      write(*,*)'First line is just a comment. Nothing is used there.'
      write(*,*)'Thereafter ignored comment lines start with #'
      write(*,*)'Arguments: -v1. -l3.5 (e.g.) adds additional command'
     $     ,' line arguments'
      write(*,*)'######################################################'
      write(*,'(2a)')'Object lines have the format: otype, a,b,c,'
     $     ,' center(3), radii(3), [extra data]'
      write(*,*)'a,b,c, set BC in Robin form aU+bU''+c=0' 
      write(*,*)'if a=b=c=0, no potential BC is applied and object'
     $     ,' is masked.'
      write(*,*)'otype indicates how to use the line. Higher bytes'
     $     ,' of otype indicate specials.'
      write(*,*)'byte-1: 1 Spheroid, 2 Cuboid, 3 Cylinder, 4 Parallelopi
     $ped, ',' 5 General cylinder,',' 6 Surface of Revolution'
      write(*,'(2a)')'Extra data for more complex objects and'
     $     ,' to specify flux accumulation.'
      write(*,*)'For Cylinder, center(3), radii(3), iaxial=number'
     $     ,' of axial axis.'
      write(*,*)'Parallelopiped, center(3)->origin. radii(3)->vec1.'
     $     ,' Then vec2(3), vec3(3).'
      write(*,*)'Non-aligned Cylinder, center(3), axial-vector(3),'
     $     ,' reference-vector(3), radius.'
      write(*,*)'S-of-R, Base(3), Apex(3), ref-vec(3), radius,'
     $     ,' n-pairs(1), zr-pairs(2n).' 
      write(*,*)'Flux accumulation (follows geometry): ofluxtype,'
     $     ,'ofn1,ofn2[,ofn3].'
      write(*,*)' Indicates number-of-flux-types, sizes of uniform '
     $     ,'index arrays [direction].'
      write(*,*)' Sphere needs 2 arrays (cos(th),psi).'
     $     ,' Cylinder 3 (r,th,z) arrays.'
      write(*,*)' Cuboid (x,y,z) arrays. '
     $     ,'Pllelopiped (1,2,3) arrays. SoR not done.'
      write(*,*)'byte-1: 99 Boolean particle region,  91-3 Set mesh in'
     $     ,' dimension 1-3.'
      write(*,*)'byte-2: 1(x256) Special boundary phi=0 instead of '
     $     ,'continuity.'
      write(*,*)'byte-2: 2(x256) Point-charge (spherical) object.'
     $     ,' Special treatment'
      write(*,*)' for which 2nd of extra (ofn1) is charge magnitude='
     $     ,' coulomb-phi at radius.'
      write(*,*)'byte-2: 3(x256) Variable-coefficients.'
     $     ,' Gradient vectors of c,b,a follow'
     $     ,' flux data [which must be specified].'
      write(*,*)'byte-2: 4(x256) [a,b,]c float by local fluxes:'
     $     ,' insulating.'
      write(*,*)'byte-2: 5(x256) [a,b,]c floating total fluxes.'
      write(*,*)
      write(*,'(a)')
     $     'Boolean particle region 99, n1, n1*values, n2, values,.. 0:'
     $     ,'The nk values are ored together. The groups are anded.'
     $     ,'value n means inside object n, -n outside, 0 do nothing.'
      write(*,'(2a)')
     $     ' E.g. -1: outside object 1, and 2 inside object 2:'
     $     ,' 99, 1, -1, 1, 2, 0'
     $     ,' or: (inside 3 OR inside 4) AND outside 5 [ors first]:'
     $     ,' 99, 2, 3, 4, 1, -5, 0'
      write(*,*)
      write(*,'(a)')'Mesh setting: 9d,is1,is2,...,isN,0,xm1,xm2...xmN'
      write(*,'(2a)')'Set structure for dimension d. N-1 blocks.'
     $     ,' Steps between is1,is2 '
      write(*,'(2a)')'equally spaced between xm1,xm2: x(is1)=xm1;'
     $     ,'x(is2)=xm2, etc.'
      write(*,*)'E.g. 92,1,12,20,32,0;-5.,-1.,1.,5. y has 3 domains: '
     $     ,'1-12 covers -5. to -1.;',' 12-20 covers -1. to 1.;'
     $     ,'20-32 covers 1. to 5.'
      write(*,*)'Default equivalent to 9d,1,32,0,-5.,5.'
8     format(9a)
      write(*,*)
      write(*,8)'Face Boundary Condition Setting: 10d,A,B,C0[,Cx,Cy,Cz'
     &         ,'] where d is face index,'
      write(*,8)'  [1-6] equivalent to -bfd,A,B,C0[,Cx,Cy,Cz]'
      write(*,8)'Periodic Face Potential Condition Setting: 11d toggles'
     &         ,' periodity dimension d.'
      write(*,8)'Periodic Particle Region Setting: 98,ix,iy,iz equivale'
     &         ,'nt to -ppix,iy,iz.'
      write(*,8)'Subtractive object association; 89, nassoc, nrmv,'
     $     ,' irmv1,irmv2,...,irmvnrmv.'
      write(*,*)'Set objects irmv1-irmvnrmv as subtractive for nassoc.'
      end

c**********************************************************************
      subroutine readgeom(filename,myid,ifull,CFin,iCFcount,LPF,ierr
     $     ,argline)
c Read the geometric data about objects from the file filename
      character*(*) filename
      integer myid
      integer ifull(*)
      character*256 cline
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'meshcom.f'
      include 'partcom.f'
      real CFin(3+ndims,2*ndims)
      logical LPF(ndims)
      character*256 argline
c Common data containing the object geometric information. 
c Each object, i < 64 has: type, data(odata).
c      integer ngeomobjmax,odata,ngeomobj
c      parameter (ngeomobjmax=...,odata=...)
c      real obj_geom(odata,ngeomobjmax)
      intrinsic IBCLR
      real valread(2*nspec_mesh)
      logical lbounded
      external lbounded

c silence spurious warning
      ist=0

      ierr=0
c Zero the obj_geom data.
      do j=1,odata
         do i=1,ngeomobjmax 
            obj_geom(j,i)=0.
         enddo
      enddo
c Read
      open(1,file=filename,status='old',err=101)
      iline=1
      read(1,'(a)',end=902)cline
c First line must be the number of dimensions. Pointless. Remove.
c      read(cline,*,err=901)nd
      nd=3

      if(myid.eq.0)write(*,'(a,a50)')'Reading objects from file: '
     $     ,filename
      if(myid.eq.0)write(*,*)'Object Descr   type',
     $     '       (BCs)              (center)            (radii)'
c Loop over lines of the input file.
 1    iline=iline+1
      read(1,'(a)',end=902)cline
      if(cline(1:1).eq.'#') goto 1
      if(cline(1:6).eq.'      ') goto 1
      if(cline(1:10).eq.'Arguments:')then
         argline(lentrim(argline)+2:)=cline(11:)
         goto 1
      endif

      read(cline,*,err=901)type
c Use only lower byte.
      itype=int(type)
      type=int(itype - 256*(itype/256))
      ngeomobj=ngeomobj+1
      if(ngeomobj.gt.ngeomobjmax)then
         write(*,*)'More objects than can be managed in ',ngeomobjmax
         goto 901
      endif
c Start of type choices
      if(type.eq.1.)then
c------------------------------------------------
c Sphere Read the geometry definition variables and flux counts.
         read(cline,*,err=901,end=801)
     $        (obj_geom(k,ngeomobj),k=1,oradius+nd-1)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn2)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
 801     if(myid.eq.0)write(*,820)ngeomobj,' Spheroid '
c Sphere has just one facet and we must make n3=1 too:
         obj_geom(ofn3,ngeomobj)=1
         obj_geom(offc,ngeomobj)=1
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
      elseif(type.eq.2.)then
c------------------------------------------------
c Cuboid has 6 facets and uses numbering in three coordinates:
         read(cline,*,err=901,end=802)
     $        (obj_geom(k,ngeomobj),k=1,oradius+nd-1)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn3)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
 802     if(myid.eq.0)write(*,820)ngeomobj,' Cuboid '
         obj_geom(offc,ngeomobj)=2*nd
         if(myid.eq.0)then 
            write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
            do k=1,ndims
               if(obj_geom(oradius+k-1,ngeomobj).eq.0.)
     $              stop 'Zero volume cube not allowed'
            enddo
         endif
      elseif(type.eq.3.)then
c------------------------------------------------
c Cylinder
         read(cline,*,err=901,end=803)
     $        (obj_geom(k,ngeomobj),k=1,ocylaxis)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn3)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
 803     if(myid.eq.0)write(*,820)ngeomobj,' Cylinder '
c 3 facets.
         obj_geom(offc,ngeomobj)=3
         if(obj_geom(ocylaxis,ngeomobj).le.0. .or.
     $        obj_geom(ocylaxis,ngeomobj).ge.4. )then
            write(*,*)'Geometry ERROR. No cylinder axial direction in:'
            write(*,*)cline
            stop
         endif
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
      elseif(type.eq.4.)then
c------------------------------------------------
c Parallelopiped also serves as general cuboid.
         read(cline,*,err=901,end=804)
     $        (obj_geom(k,ngeomobj),k=1,oradius+nd*nd-1)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn3)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
c     $        (obj_geom(k,ngeomobj),k=1,odata)
 804     if(myid.eq.0)write(*,820)ngeomobj,
     $        ' Pllelopiped '
         obj_geom(offc,ngeomobj)=2*nd
         call plleloinit(ngeomobj)
         if(myid.eq.0)write(*,822)(obj_geom(k,ngeomobj),
     $        k=1,1+nd*(1+nd)+3)
      elseif(type.eq.5)then
c------------------------------------------------
c Non-aligned cylinder
         read(cline,*,err=901,end=805)
     $        (obj_geom(k,ngeomobj),k=1,ocylrad)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn3)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
 805     if(myid.eq.0)write(*,820)ngeomobj,' General Cyl '
         obj_geom(offc,ngeomobj)=3
         if(obj_geom(ocylrad,ngeomobj).le.0.)then
            write(*,*)'Geometry ERROR. No positive cylinder radius in:'
            write(*,*)cline
            stop
         endif
         if(myid.eq.0)write(*,822)(obj_geom(k,ngeomobj),k=1,1+4*nd+1)
         call cylinit(obj_geom(1,ngeomobj))
      elseif(type.eq.6.or.type.eq.7)then
c------------------------------------------------
c Surface of revolution. Now like cylinder except for Base/Apex
         do k=1,ovlen
            obj_geom(opdiv+k-1,ngeomobj)=0
         enddo
         obj_geom(ofluxtype,ngeomobj)=0
         isrnpair=0
c Try to read the base,apex,npair,r,z;fluxtype,ntheta,div
c For type 6 there is one more div than pair, for type 7 one less.
         read(cline,*,err=901,end=806)
     $        (obj_geom(k,ngeomobj),k=1,ocylrad)
     $        ,isrnpair
     $        ,(obj_geom(opr+k-1,ngeomobj),
     $        obj_geom(opz+k-1,ngeomobj),k=2,isrnpair+1)
     $        ,obj_geom(ofluxtype,ngeomobj),obj_geom(ofn1,ngeomobj)
     $        ,(obj_geom(opdiv+k-1,ngeomobj),k=1,isrnpair+1)
 806     if(myid.eq.0)then
            write(*,820)ngeomobj,' Surf-of-Revln'
            write(*,822)(obj_geom(k,ngeomobj),k=1,oapex+2)
         endif
         obj_geom(onpair,ngeomobj)=isrnpair
c Make the onpair equal to the number of segment ends.
         if(type.eq.6)obj_geom(onpair,ngeomobj)=isrnpair+2
c         write(*,*)'npair',obj_geom(onpair,ngeomobj)
         if(isrnpair.le.0.or.isrnpair.gt.ovlen-2.and.myid.eq.0)then
            write(*,*)' Impossible number of pairs specified',isrnpair
            stop
         endif
         call srvinit(ngeomobj,type,myid)         
         if(myid.eq.0)then
            np=obj_geom(onpair,ngeomobj)
            write(*,'(a,i3,i3,a)')' rz-pairs',isrnpair,np,':'
            write(*,'(10f7.2)')(obj_geom(opr+k-1,ngeomobj)
     $           ,obj_geom(opz+k-1,ngeomobj),k=1,np)
            if(obj_geom(ofluxtype,ngeomobj).gt.0)then
               write(*,'(a,20i3)')' Face divisions'
     $              ,(int(obj_geom(opdiv+k-1,ngeomobj)),k=1,np-1)
            endif
         endif
      elseif(type.eq.99)then
c------------------------------------------------
c Specify the particle region.
         read(cline,*,err=901,end=899)idumtype,ibool_part
 899     if(myid.eq.0)write(*,898)
     $        ngeomobj,idumtype,(ibool_part(i),i=1,16)
 898     format(i2,' Boolean ',17i4)
c Don't count this as an object.
         ngeomobj=ngeomobj-1
         goto 1
      elseif(type.eq.98)then
         read(cline,*,err=901,end=897)idumtype,ipartperiod
 897     if(myid.eq.0)write(*,'(''PartPeriod'',5i6)')
     $        ngeomobj,idumtype,ipartperiod
         ngeomobj=ngeomobj-1
         goto 1
      elseif(type.gt.90.and.type.le.90+ndims)then
c------------------------------------------------
c Mesh specification for a particular dimension. 
         id=int(type-90)
         ist=0
         read(cline,*,err=901,end=880)idumtype,valread
 880     do i=1,nspec_mesh
            ist=i
            imeshstep(id,i)=int(valread(i))
            if(imeshstep(id,i).gt.ifull(id))then
               write(*,*)cline
               write(*,*)'Readgeom: Meshpos',imeshstep(id,i)
     $              ,'  too large for ifull',ifull(id),id
               stop
            endif
            if(valread(i).eq.0.)goto 881
         enddo
 881     do j=1,ist-1
            xmeshpos(id,j)=valread(ist+j)
c            write(*,*)id,j,xmeshpos(id,j)
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
      elseif(type.ge.100.and.type.lt.100+2*ndims+1)then
c------------------------------------------------
c Face boundary conditions.
         id=int(type)-100
         if(iCFcount.eq.0)then
c Reset all.
            do ii=1,2*ndims
               do idc=1,3+ndims
                  CFin(idc,ii)=0.
               enddo
            enddo
         endif
         read(cline,*,err=901,end=882)idumtype,(CFin(ii,id),ii=1,3
     $        +ndims)
 882     continue
         LPF(mod(id-1,ndims)+1)=.false.
         iCFcount=iCFcount+1
c Don't count this as an object.
         ngeomobj=ngeomobj-1
         goto 1         
      elseif(type.gt.110.and.type.lt.110+ndims)then
c------------------------------------------------
c Periodic boundary condition on this dimension
         if(iCFcount.eq.0)then
c Reset all if this is the first face call.
            do ii=1,2*ndims
               do idc=1,3+ndims
                  CFin(idc,ii)=0.
               enddo
            enddo
         endif
         id=int(type)-110
         LPF(id)=.not.LPF(id)
         do ii=id,id+ndims,ndims
            CFin(1,ii)=0.
            CFin(2,ii)=1.
         enddo
         iCFcount=iCFcount+1
c Don't count this as an object.
         ngeomobj=ngeomobj-1
      elseif(type.eq.89.or.type.eq.88)then
c------------------------------------------------
c Subtraction
         read(cline,*,err=901,end=881)idumtype,iassoc,normv(iassoc)
     $        ,(ormv(k,iassoc),k=1,normv(iassoc))
         if(myid.eq.0)write(*,*)ngeomobj,' Subtracting association',type
     $        ,iassoc,normv(iassoc),(ormv(k,iassoc),k=1,normv(iassoc))
c Don't count this line as an object.
         ngeomobj=ngeomobj-1
         if(iassoc.gt.ngeomobj)stop
     $        'Input error subtracting from non-existent object'
         if(type.eq.88)then
c Implement the implied subtractive effects of each subtracted object.
c Make the complement of object iassoc subtractive for each subtracted.
         do k=1,normv(iassoc)
c isoc is the subtractive object
            isoc=abs(ormv(k,iassoc))
            if(isoc.gt.ngeomobj)stop
     $           'Input error subtracting non-existent object' 
c Increment its list of subtractives.
            normv(isoc)=normv(isoc)+1
c Make its last entry the complement of the additive.
            ormv(normv(isoc),isoc)=-iassoc
            if(obj_geom(oabc,isoc).eq.0.and.
     $           obj_geom(oabc+1,isoc).eq.0.and.
     $           obj_geom(oabc+2,isoc).eq.0.)then
c Set subtractive boundary conditions equal to additive
               do ii=1,ndims
                  obj_geom(oabc+ii-1,isoc)=obj_geom(oabc+ii-1,iassoc)
               enddo
               ifield_mask=IBSET(ifield_mask,isoc-1)
            endif
            if(myid.eq.0)write(*,*)'  Set subtractive surfaces',isoc
     $           ,normv(isoc),(ormv(ii,isoc),ii=1,normv(isoc))
     $           ,(obj_geom(oabc+ii-1,isoc),ii=1,3*ndims)
         enddo
         endif
      endif
c End of type choices
c------------------------------------------------
c If this is a null boundary condition clear the relevant bit.
      if(obj_geom(oabc,ngeomobj).eq.0.
     $     .and. obj_geom(oabc+1,ngeomobj).eq.0.
     $     .and. obj_geom(oabc+2,ngeomobj).eq.0.)
     $     ifield_mask=IBCLR(ifield_mask,ngeomobj-1)
c If this is a point-charge object, set the relevant mask bit.
      if(itype/256.eq.2)then
         if(
     $        obj_geom(oradius,ngeomobj).ne.obj_geom(oradius+1,ngeomobj)
     $        .or.
     $        obj_geom(oradius,ngeomobj).ne.obj_geom(oradius+2,ngeomobj)
     $        )then
            write(*,*)'Unequal radii not allowed for ptch'
            stop
         endif
         iptch_mask=IBSET(iptch_mask,ngeomobj-1)
      elseif(itype/256.eq.3)then
         do k=ocgrad,oagrad+2
            if(obj_geom(k,ngeomobj).ne.0.)goto 110
         enddo
         write(*,*)'ERROR in Readgeom.'
     $        ,' No non-zero gradient cooefficient'
         stop
 110     continue
         write(*,*)'   Gradients:'
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
      endif
 820  format(i2,a,$)
 821  format(f5.0,15f7.3)
 822  format(f3.0,24f6.2)
      goto 1

 901  write(*,*)'Readgeom error reading line',iline,':'
      write(*,*)cline
      
 902  continue
c Set whether particle region has a part inside an object.
      lboundp=lbounded(ibool_part,ifield_mask)
      if(lboundp.and.rjscheme(1:4).eq.'cart')then
         write(*,'(3a)')'ERROR: using cartesian injection'
     $        ,' with bounded region.'
         stop
      endif
      return

 101  write(*,*) 'Readgeom File ',filename(1:lentrim(filename))
     $     ,' could not be opened.'
      ierr=1

      end
c****************************************************************
      subroutine cylinit(objg)
c Initialize the non-aligned cylinder contravariant vectors
c They are such that they yield the covariant coefficients relative
c to a unit cylinder.
c There are two input vectors (that are overwritten on return): 
c    v the axis (ovec), and 
c    u the ref-vec (ovec+ndims), plus
c    r the single scalar radius (ovec+2*ndims=ocylrad)
c The contravariant vectors are in the order vb,vg,va, (ocontra...)
c where vb gives the component of u perpendicular to va, vg is the 
c cross product of va and u, and va is the axial vector v/|v|^2.
c Consequently, vacontra.x=projected length in axial direction, in units
c of the length of v.
c Also the other two contra vector are normalized so that they yield 
c the length in units of the radius, when dotted into a vector.
c Finally the ovec storage positions are overwritten with
c covariant vectors equal to contra/|contra|^2.
      include 'ndimsdecl.f'
      include '3dcom.f'
      real objg(odata)

      radius=objg(ovec+2*ndims)
      vamag=0.
      umag=0.
      uv=0.
      do i=1,ndims
c v^2, u^2
         vamag=vamag+objg(ovec+i-1)**2
         umag=umag+objg(ovec+ndims+i-1)**2
c u.v
         uv=uv+objg(ovec+i-1)*objg(ovec+ndims+i-1)
c v x u
         ip=mod(i,ndims)
         im=mod(i+1,ndims)
         objg(ocontra+ndims+i-1)=
     $        objg(ovec+ip)*objg(ovec+ndims+im)
     $        -objg(ovec+im)*objg(ovec+ndims+ip)
      enddo
      if(vamag.eq.0.)stop 'cylinit error axial vector zero'
      vbmag=0.
      vgmag=0.
      do i=1,ndims
         vacontra=objg(ovec+i-1)/vamag
         objg(ocontra+2*ndims+i-1)=vacontra
         vbcontra=(objg(ovec+ndims+i-1)-uv*vacontra)
         objg(ocontra+i-1)=vbcontra
         vbmag=vbmag+vbcontra**2
         vgmag=vgmag+objg(ocontra+ndims+i-1)**2
      enddo
      if(vbmag.eq.0. .or. vgmag.eq.0)
     $     stop 'cylinit error perp vector zero'
      vbmag=sqrt(vbmag)
      vgmag=sqrt(vgmag)
      c1mag=0.
      c2mag=0.
      c3mag=0.
      do i=1,ndims
         objg(ocontra+i-1)=objg(ocontra+i-1)/(vbmag*radius)
         objg(ocontra+ndims+i-1)=objg(ocontra+ndims+i-1)
     $        /(vgmag*radius)
         c1mag=c1mag+objg(ocontra+i-1)**2
         c2mag=c2mag+objg(ocontra+ndims+i-1)**2
         c3mag=c3mag+objg(ocontra+2*ndims+i-1)**2
      enddo
c Now ocontra=vb/(|vb|*r), vg/(|vg|r), va/|va^2|
c The covariant vectors are needed for reconstruction in fluxexamine.
c They are such that \sum_j (x.contra^j) covar_j = x, which is true
c for orthogonal vectors when covar = contra/|contra|^2.
      do i=1,ndims
         objg(ovec+i-1)=objg(ocontra+i-1)/c1mag
         objg(ovec+ndims+i-1)=objg(ocontra+ndims+i-1)/c2mag
         objg(ovec+2*ndims+i-1)=objg(ocontra+2*ndims+i-1)/c3mag
      enddo
c      write(*,*)'Covariant and contravariant:'
c      write(*,'(9f8.4)')(objg(ovec+i-1),i=1,18)
c      stop

      end
c****************************************************************
      subroutine plleloinit(iobj)
c Initialize the iobj by calculating the contravariant 
c vectors from the covariant vectors.
      include 'ndimsdecl.f'
      include '3dcom.f'

      triple=0.
      do j=1,ndims
         jpv=ovec+(j-1)*ndims-1
c Other vectors:
         jpv2=ovec+mod(j,ndims)*ndims-1
         jpv3=ovec+mod(j+1,ndims)*ndims-1
         jpc=ocontra+(j-1)*ndims-1
c Set obj_geom(jpc..,iobj) equal to the cross product between the other
c vectors.
         do i=1,ndims
            i2=mod(i,ndims)+1
            i3=mod(i+1,ndims)+1
            obj_geom(jpc+i,iobj)=(obj_geom(jpv2+i3,iobj)*obj_geom(jpv3
     $           +i2,iobj)-obj_geom(jpv2+i2,iobj)*obj_geom(jpv3+i3
     $           ,iobj))
         enddo
c calculate the scalar triple product the first time:
         if(j.eq.1)then
            do i=1,ndims
               triple=triple+obj_geom(jpc+i,iobj)*obj_geom(jpv+i,iobj)
            enddo
         endif
         if(triple.eq.0.)then
            write(*,*)'Parallelopiped of zero volume ERROR.'
            write(*,*)(obj_geom(k,iobj),k=1,21)
            write(*,*)jpv,jpc
            stop
         endif
c normalize
         do i=1,ndims
            obj_geom(jpc+i,iobj)=obj_geom(jpc+i,iobj)/triple
         enddo
      enddo

      end
c****************************************************************
      subroutine srvinit(iobj,type,myid)
      include 'ndimsdecl.f'
      include '3dcom.f'
      if(type.eq.6.)then
c Insert the base and apex into the z,r pair arrays.
         obj_geom(opz,iobj)=0.
         obj_geom(opr,iobj)=0.
         j=int(obj_geom(onpair,iobj))
         obj_geom(opz+j-1,iobj)=1.
         obj_geom(opr+j-1,iobj)=0.
c Check other z values for consistency.
         z=0.
         do i=2,j-1
            z1=obj_geom(opz+i-1,iobj)
            if(z1.le.z.or.z1.ge.1.)then
               write(*,*)'z values must be monotonic 0.<z<1.'
     $              ,(obj_geom(opz+k-1,iobj),k=1,j)
               stop
            endif
            z=z1
         enddo
      else
c type 7 reads the entire wall. Shuffle in.
         do i=1,obj_geom(onpair,iobj)
            obj_geom(opz+i-1,iobj)=obj_geom(opz+i,iobj)
            obj_geom(opr+i-1,iobj)=obj_geom(opr+i,iobj)
         enddo
c Inconsistent if not closed and either end not on axis
         if((obj_geom(opz,iobj).ne.obj_geom(opz+i-2,iobj).or.
     $        obj_geom(opr,iobj).ne.obj_geom(opr+i-2,iobj)).and.
     $        (obj_geom(opr,iobj).ne.0..or.obj_geom(opr+i-2,iobj).ne.0.)
     $        )then
            write(*,*)'rz-contour not closed and ends not on axis'
            stop
         endif
      endif

      do j=1,ndims
c Convert the apex into apex-base and call it ovec
c ovec is actually the same as oapex
         obj_geom(ovec+j-1,iobj)=obj_geom(oapex+j-1,iobj)
     $        -obj_geom(obase+j-1,iobj)
      enddo

c Initialize Surface of Revolution contravariant vectors E.g. the axial
c ovec whose length is such that axial.(apex-base)=1. So it is
c (apex-base)/|apex-base|^2. See cylinit for fuller description.
      call cylinit(obj_geom(1,iobj))

c Flux accounting requires all relevant facets have non-zero div 
c number. By this point there is one less facet than npair.
      if(obj_geom(ofluxtype,iobj).ne.0)then
         if(obj_geom(ofn1,iobj).le.0)goto 1
         do i=1,obj_geom(onpair,iobj)-1
            if(obj_geom(opdiv+i-1,iobj).le.0)goto 1
         enddo
c Here if all is well.
         return
 1       continue
c Else There's an inadequacy in the flux collection specification
         if(myid.eq.0)then
            write(*,*)'Flux specification for object' ,iobj
     $           ,' is in ERROR at division',i
            write(*,*)int(obj_geom(ofluxtype,iobj))
     $           ,int(obj_geom(ofn1,iobj)),',',
     $         (obj_geom(opdiv+i-1,iobj),i=1,obj_geom(onpair,iobj)-1)
     $           ,' Flux turned off.'
         endif
         obj_geom(ofluxtype,iobj)=0
      endif
      end
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
c     $           (obj_geom(ovec+j-1,i)))**2
c Using cylinit it has been moved to 
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
c Then find the r value for this z-index.
         r=obj_geom(opr+iz-1,i)+(zind-iz)
     $        *(obj_geom(opr+iz,i)-obj_geom(opr+iz-1,i))
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
         r=sqrt(r)
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
c*************************************************************
c Initialize the iregion flags of the existing nodes with boundary
c object data.
      subroutine iregioninit(ifull)
c      integer ndims
      include 'ndimsdecl.f'
      integer ifull(ndims)

      include 'objcom.f'
      include 'meshcom.f'
      include '3dcom.f'

      integer ix(ndims)
      real x(ndims)

      do i=1,oi_cij
         ipoint=idob_cij(ipoint_cij,i)
c Convert index to multidimensional indices.
         call indexexpand(ndims,ifull,ipoint,ix)
         do k=1,ndims
            x(k)=xn(ixnp(k)+ix(k))
         enddo
c Store in object-data.
         idob_cij(iregion_cij,i)=insidemask(ndims,x)

c         write(*,101)i,ipoint,idob_cij(iregion_cij,i),x
c 101     format(3i8,5f10.4)
      enddo

      end
c*****************************************************************
      subroutine reportfieldmask()
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer ipb(32),ifb(32)
c Calculate the bits of the field mask and iptch_mask.
      ifd=ifield_mask
      ipp=iptch_mask
      do i=1,32
         ipb(32-i+1)=ipp - 2*(ipp/2)
         ipp=ipp/2
         ifb(32-i+1)=ifd - 2*(ifd/2)
         ifd=ifd/2
      enddo
c      write(*,*)'Initializing Object Regions:No, pointer, region, posn.'
c This is an unportable extension. Hence the calculation above.
c      write(*,'('' Mask='',i11,'' ='',b32.32)')ifield_mask,ifield_mask
      write(*,'('' Field Mask='',i11,'' ='',32i1)')ifield_mask,ifb
c      if(ipp.ne.0)
      write(*,'('' Ptch Mask= '',i11,'' ='',32i1)')iptch_mask,ipb
      end
c*******************************************************************
      subroutine objsetabc(iobject,a,b,c)
c Set/Reset the boundary conditions for an object that has already been
c read in. 
      integer iobject
      real a,b,c

      include 'ndimsdecl.f'
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
C
C*******************************************************************
