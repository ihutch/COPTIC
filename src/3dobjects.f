c This replaces the block data program com3dset which causes giant objects.
      subroutine blockdatainit()
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptchcom.f'
      include '3dcom.f'
      include 'partcom.f' 
      ngeomobj=0
c Default track no objects.
c Default no subtractive objects
      do i=1,ngeomobjmax
         nf_map(i)=0
         normv(i)=0
      enddo
c And no reverse-map pointers.
      do i=1,nf_obj
         nf_geommap(i)=0
c      data nf_geommap/nf_obj*0/
      enddo
c Default particle region: zero boolean. Particles everywhere.
      do i=1,ibtotal_part
         ibool_part(i)=0
      enddo
c Set all the bits of ifield_mask: =2**31-1=2*(2**30-1)+1 avoid overflow.
      ifield_mask=2*(2**30-1)+1
c Normally there's no external field.
      lextfield=.false.
c Mesh default initialization (meshcom.f)
      do i=1,ndims
         extfield(i)=0
         do j=1,nspec_mesh
            imeshstep(i,j)=0
            xmeshpos(i,j)=0
            if(j.le.1)then
               imeshstep(i,j)=1
               xmeshpos(i,j)=-5.
            elseif(j.le.2)then
               imeshstep(i,j)=32
               xmeshpos(i,j)=5.
            endif
         enddo
      enddo
c Default no point charges:
      iptch_mask=0
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
     $ped, ',' 5 General cylinder,',' 6,7 Surface of Revolution'
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
     $     ,' Cylinder 3 (r,th,z). Cuboid 3 (x,y,z). '
      write(*,*)' Pllelopiped 3 (1,2,3).'
     $     ,' SoR 1+nfs (thdivs,n1,n2,...,nnfs) nfs=npairs[+1type6].'
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
      character*512 cline
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptchcom.f'
      include '3dcom.f'
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
 2    continue
      iline=1
      read(1,'(a)',end=902)cline
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
            np=int(obj_geom(onpair,ngeomobj))
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
      elseif(type.gt.110.and.type.le.110+ndims)then
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
         read(cline,*,err=901,end=883)idumtype,iassoc,normv(iassoc)
     $        ,(ormv(k,iassoc),k=1,normv(iassoc))
         goto 884
 883     write(*,*)'Error reading subtraction line',cline
         stop
 884     continue
         if(myid.eq.0)write(*,'(i2,a,f5.0,20i3)')
     $        ngeomobj,
     $        ' Subtract',type ,iassoc,normv(iassoc)
     $        ,(ormv(k,iassoc),k=1,normv(iassoc))
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
c Write only for debugging.
            if(.false. .and. myid.eq.0)then
               write(*,'(a,i3,i2,10i3)')
     $              '  Setting subtractive surfaces',isoc,normv(isoc),
     $              (ormv(ii,isoc),ii=1,normv(isoc))
c     $           ,(obj_geom(oabc+ii-1,isoc),ii=1,3*ndims)
            endif
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
      if(filename(1:lentrim(filename)).ne.'copticgeom.dat')then
         filename='copticgeom.dat'
         open(1,file=filename,status='old',err=102)  
         goto 2
 102     write(*,*) 'And copticgeom.dat could not be opened either.'
      endif
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
c Save the radial scale.
      objg(orscale)=radius
c      write(*,*)'Cylinit radius=',radius,' opz=',opz,' odata=',odata
c     $     ,'ofluxtype=',ofluxtype
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
         do i=1,int(obj_geom(onpair,iobj))
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
         do i=1,int(obj_geom(onpair,iobj)-1)
            if(obj_geom(opdiv+i-1,iobj).le.0)goto 1
         enddo
c Here if all is well.
         return
 1       continue
c Else There's an inadequacy in the flux collection specification
         if(myid.eq.0)then
            write(*,*)'Flux specification for object' ,iobj
     $           ,' is in ERROR at division',i
            write(*,*)int(obj_geom(ofluxtype,iobj)) ,int(obj_geom(ofn1
     $           ,iobj)),',', (obj_geom(opdiv+i-1,iobj),i=1
     $           ,int(obj_geom(onpair,iobj)-1)) ,' Flux turned off.'
         endif
         obj_geom(ofluxtype,iobj)=0
      endif
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
      include 'griddecl.f'
      include 'ptchcom.f'
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
c*******************************************************************
c*****************************************************************
      subroutine phipset(myid)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include '3dcom.f'

      if(obj_geom(oabc,1).ne.0)then
         phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
         if(myid.eq.0)write(*,*)'Object 1 potential=',phip
      elseif(obj_geom(oradius,1).ne.0.)then
         phip=obj_geom(omag,1)*obj_geom(oradius,1)
         if(myid.eq.0)write(*,*)'Potential from point charge'
     $        ,obj_geom(omag,1),' at radius ',obj_geom(oradius,1)
     $        ,' Charge:',phip
      else
         phip=0.
         if(myid.eq.0)write(*,*)'Potential phip not set from objects.'
      endif

      end
c******************************************************************
