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
      write(*,*)'byte-1: 1 Spheroid, 2 Cuboid, 3 Cylinder, 4 Parallel'
     $     ,'opiped, 5 General cylinder'
      write(*,'(2a)')'Extra data for more complex objects and'
     $     ,' to specify flux accumulation.'
      write(*,*)'For Cylinder, center(3), radii(3), iaxial=number'
     $     ,' of axial axis.'
      write(*,*)'Parallelopiped, center(3)->origin. radii(3)->vec1.'
     $     ,' Then vec2(3), vec3(3).'
      write(*,*)'Non-aligned Cylinder, center(3), axial-vector(3),'
     $     ,' reference-vector(3), radius.'
      write(*,*)'Flux accumulation (follows geometry): ofluxtype,'
     $     ,'ofn1,ofn2[,ofn3].'
      write(*,*)' Indicates number-of-flux-types, sizes of uniform '
     $     ,'index arrays [direction].'
      write(*,*)' Sphere needs 2 arrays (cos(th),psi).'
     $     ,' Cylinder 3 (r,th,z) arrays.'
      write(*,*)' Cuboid (x,y,z) arrays. '
     $     ,'Pllelopiped (1,2,3) arrays.'
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
     $     '       (BCs)              (center)              (radii)'
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
      if(type.eq.1.)then
c Read the geometry definition variables and then the flux counts.
         read(cline,*,err=901,end=801)
     $        (obj_geom(k,ngeomobj),k=1,oradius+nd-1)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn2)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
 801     if(myid.eq.0)write(*,820)ngeomobj,' Spheroid '
c Sphere has just one facet and we must make n3=1 too:
         obj_geom(ofn3,ngeomobj)=1
         obj_geom(offc,ngeomobj)=1
         if(myid.eq.0)write(*,821)(obj_geom(k,ngeomobj),k=1,1+2*nd+3)
c         if(myid.eq.0 .and. obj_geom(ofluxtype,ngeomobj).ne.0)
c     $        write(*,821)(obj_geom(k,ngeomobj),k=ofluxtype,ofn2)
      elseif(type.eq.2.)then
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
               if(obj_geom(ocenter+k-1,ngeomobj).eq.obj_geom(oradius+k-1
     $              ,ngeomobj)) stop 'Zero volume cube not allowed'
            enddo
         endif
      elseif(type.eq.3.)then
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
c Parallelopiped also serves as general cuboid.
         read(cline,*,err=901,end=804)
     $        (obj_geom(k,ngeomobj),k=1,oradius+nd*nd-1)
     $        ,(obj_geom(k,ngeomobj),k=ofluxtype,ofn3)
     $        ,(obj_geom(k,ngeomobj),k=ocgrad,oagrad+2)
c     $        (obj_geom(k,ngeomobj),k=1,odata)
 804     if(myid.eq.0)write(*,820)ngeomobj,
     $        ' Pllelopiped '
         obj_geom(offc,ngeomobj)=2*nd
c         call plleloinit(obj_geom(1,ngeomobj))
         call plleloinit(ngeomobj)
         if(myid.eq.0)write(*,822)(obj_geom(k,ngeomobj),
     $        k=1,1+nd*(1+nd)+3)
      elseif(type.eq.5)then
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
      elseif(type.eq.99)then
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
      endif
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
c to a unit cylinder. They are in the order vb,vg,va.
c Where vb is the component of u perpendicular to va, and vg is the 
c cross product of u and va.
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
c Initialize the iobj pp_structure by calculating the contravariant 
c vectors from the covariant vectors.
      include 'ndimsdecl.f'
      include '3dcom.f'

      triple=0.
      do j=1,ndims
         jpv=pp_vec+(j-1)*ndims-1
c Other vectors:
         jpv2=pp_vec+mod(j,ndims)*ndims-1
         jpv3=pp_vec+mod(j+1,ndims)*ndims-1
         jpc=pp_contra+(j-1)*ndims-1
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
      
      inside_geom=0
      if(i.gt.ngeomobj) return

      itype=int(obj_geom(otype,i))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)
      if(itype.eq.1)then
c Coordinate-Aligned Spheroid data : center(ndims), semi-axes(ndims) 
         r2=0
         do k=1,ndims
            r2=r2+((x(k)-obj_geom(ocenter-1+k,i))/
     $           obj_geom(oradius-1+k,i))**2
         enddo
         if(r2.lt.1.)inside_geom=1
      elseif(itype.eq.2)then
c Coordinate-Aligned Cuboid data:
         do k=1,ndims
            xk=x(k)-obj_geom(ocenter-1+k,i)
            if(abs(xk).ge.abs(obj_geom(oradius-1+k,i)))return
         enddo
         inside_geom=1
      elseif(itype.eq.3)then
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
c Continuity alone should be applied
c on the other (inactive) side of the surface.  A negative value of c2
c is used to indicate that just the appropriate continuity condition is
c to be applied, but c1,c2,c3, should all be returned with the magnitudes
c that are applied to the active side of the surface, with consistently
c reversed signs. (e.g. a=1, b=1, c=-2 -> a=-1, b=-1, c=2)
c
c A fraction of 1 causes all the bounding conditions to be ignored.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include '3dcom.f'
c      integer debug

      real xx(10),xd(10)
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
c                  call sphereinterp(ndims,0,xp1,xp2
c     $                 ,obj_geom(ocenter,i),obj_geom(oradius,i)
c     $                 ,fraction,f2,sd,C,D)
            elseif(itype.eq.2)then
c Coordinate-aligned cuboid.               
               call cubefsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.3)then
c Coordinate aligned cylinder.
               call cylfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            elseif(itype.eq.4)then
c Parallelopiped.
               call pllelofsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
               if(fraction.gt.1.)fraction=1.
               if(fraction.lt.0.)fraction=1.
            elseif(itype.eq.5)then
c Non-aligned cylinder needed
               call cylgfsect(ndims,xp1,xp2,i,ijbin,sd,fraction)
            else
               write(*,*)"Unknown object type",obj_geom(otype,i),
     $              " in potlsect"
            endif
         else
c Null BC. Ignore.
            fraction=1.
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
               do idi=1,ndims
c    Address of mesh point. 
                  ix=indi(idi)+ixnp(idi)+1
                  xx(idi)=xn(ix)
c This needs to be generalized for other objects.
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
      function ireg3(i,j,k,ifull,cij)
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
c> Convert from contravariant coefficients to world cartesian.
c> Covariant and Contravariant vectors are stored in obj_geom(...,iobj)
c> xcontra and xw can be the same storage positions if desired.
      subroutine contra3world(mdims,xcontra,xw,iobj)
c> Number of dimensions. 
      integer mdims
c> Object number
      integer iobj
c> World coordinates
      real xw(mdims)
c> Contravariant coordinates
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
         xcl(j)=0.
      enddo
c Cartesian world relative to center
      do i=1,ndims
         xd(i)=xw(i)-obj_geom(ocenter+i-1,iobj)
c Contra coefficients are obtained by dotting with contra vectors
         do j=1,ndims
            xcl(j)=xcl(j)
     $           +xd(i)*obj_geom(ocontra+ndims*(j-1)+i-1,iobj)
         enddo
         xcontra(i)=xcl(i)
      enddo
      end
