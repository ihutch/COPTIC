c*****************************************************************************
      subroutine surf3d(x,y,z,iLx,nx,ny,isw,work)
c Draw a 3-d surface of z(x,y), dim(iLx\nx,ny), viewed from (xe,ye,ze) scaled.
c This version designed for non-monotonic x,y, does sophisticated web order.
c work(0:iLx+1,0:ny+1) is a work array that is trashed.
c Switch: isw.lt.0 draw no axes.
c Second byte of isw gives web color. Third byte gives axis color.
c Lowest nibble levd (absolute value): 
c   =0 perform no scaling and use last perspective 2-D x,y
c   =1 scale to fit the region, 1-d x,y
c   =2 scale to fit the region, 2-d x,y
c   =4 perform no scaling, 1-d x,y.
c Second nibble (16+): isw for the surfdr3 call. bit0=1 => smooth tri.
c bit1=1 use the first 5 elements of work for the directional vector.
c Eye obtained from file eye.dat, or default if none exists.
      integer iLx, nx,ny,isw
      real x(iLx,*),y(iLx,*),z(iLx,ny),work((iLx+2)*(ny+2))
      integer colw,cola
      real x2,y2,z2
      real zmin,zmax
      real d(5)
      data x2/0./y2/0./z2/0./
      save

      cola=isw/256
      levd=abs(isw-(isw/16)*16)
c levd is lowest 4 bits.
      colw=abs(cola-(cola/256)*256)
c color of web is next 8 bits
      cola=abs(cola/256 -(cola/65536)*256)
c color of axes is next 8 bits.
      isw=isw/16 - (isw/256)*16
c isw is second nibble (bits 4-7).
      if(colw.ge.0)then
         call color(colw)
      endif
c      write(*,*)'levd,cola,colw,isw',levd,cola,colw,isw
      if((isw/2-(isw/4)*2).ne.0)then
c Set the d-vector from work
         do k=1,5
            d(k)=work(k)
         enddo
         write(*,*)'Set d-vector:',d
      endif
c Set the top and bottom horizons.
      call hidinit(0.,1.)
      if(levd.eq.0)then
         call surfdr3(x,y,z,iLx,nx,ny,work,isw,d)
      elseif(levd.eq.4)then
         call surf1dr3(x,y,z,iLx,nx,ny,work,isw,d)
      else
c     Set to non-hiding 3-D calls.
         call hdprset(0,0.)
c The first time, set the perspective transform.
c But not thereafter to allow rotatezoom control.
         if(x2.eq.0..and.y2.eq.0..and.z2.eq.0.)then
            call geteye(x2,y2,z2)
            call trn32(0.,0.,0.,x2,y2,z2,1)
         endif
c     Set the scaling.
         call minmax2(z,iLx,nx,ny,zmin,zmax)
         itics=5
         call fitrange(zmin,zmax,itics,ipow,fac10,delta,first,xlast)
         zmin=first
         zmax=xlast
         if(levd.eq.2)then
            call minmax2(x,iLx,nx,ny,xmin,xmax)
            call minmax2(y,iLx,nx,ny,ymin,ymax)
            call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
            call surfdr3(x,y,z,iLx,nx,ny,work,isw,d)
         elseif(levd.eq.1)then
            call minmax(x,nx,xmin,xmax)
            call minmax(y,ny,ymin,ymax)            
            call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
c            write(*,*)'Calling surf1dr3'
            call surf1dr3(x,y,z,iLx,nx,ny,work,isw,d)
         endif
      endif
      if(cola.ne.0) call color(cola)
      if(isw.ge.0)then
c Draw axes.
         call axproj(igetcorner())
      endif
      end
c********************************************************************
c Attempt was made to create a webdrw that drew in the order of closest
c to furthest. However, hiding does not in the end work correctly for 
c that. So the surf3d was done instead.
c********************************************************************
      subroutine surfdr3(x,y,z,iLx,nx,ny,work,isw,d)
c Draw a 3-d surface of z(x,y), dim(iLx\nx,ny), using current scaling.
c Version using filled quadrilaterals.
c isw bit 0 set => triangular fills, else chunks (uniform quadrilaterals). 
c If not chunks, then
c isw bit 1 set => directional distance shading, 
c    using optional argument d(5)
c    d(1-3) gives direction d(4-5) gives distance limits.
c    Otherwise shade by z height.
c isw bit 2 set => face orientation shading d(1-3) ref direction.
c    Overrules bit 1 if set, otherwise according to bit 1.
c isw bit 3, 4 set => mesh is periodic in direct i,j
c Array work is used to calculate the distance of the facet from the eye.
c Facets are drawn from furthest to nearest to give correct projection.
      integer nx,ny,iLx
      real x(iLx,ny),y(iLx,ny),z(iLx,ny),work(0:iLx+1,0:ny+1)
      integer isw
      real d(5)

      real x2,y2,z2
      real dout,ddone
      parameter (dout=2.e32,ddone=1.e32)
      real xp(5),yp(5),zp(5),zz(5)
      integer icx(0:4), icy(0:4), ic
      logical lchunk,ldir,lorient,liper,ljper
      data icx/0,1,1,0,0/
      data icy/0,0,1,1,0/
      data liper/.false./ljper/.false./
      save
c
c      write(*,*)'isw',isw
      lchunk=.true.
      if(isw-2*(isw/2) .ne.0)lchunk=.false.
      ldir=.false.
      if(isw/2-2*(isw/4) .ne.0)ldir=.true.
      if(isw/4-2*(isw/8) .ne.0)then
         lorient=.true.
         ldir=.false.
      else
         lorient=.false.
      endif
      if(isw/8-2*(isw/16) .ne.0)liper=.true.
      if(isw/16-2*(isw/32) .ne.0)ljper=.true.
      incolor=igetcolor()
c Get eye position x0,y0,z0 are unused here.
      call trn32(x0,y0,z0,x2,y2,z2,-1)
c      write(*,*)'Eye=',x2,y2,z2
c Set up the work array.
      zmin=ddone
      zmax=-ddone
      do j=1,ny-1
         do i=1,nx-1
            if(z(i,j).gt.zmax)zmax=z(i,j)
            if(z(i,j).lt.zmin)zmin=z(i,j)
            work(i,j)=0.
            do ic=0,3
               work(i,j)=work(i,j)
     $              +x(i+icx(ic),j+icy(ic))*x2
     $              +y(i+icx(ic),j+icy(ic))*y2
     $              +z(i+icx(ic),j+icy(ic))*z2

            enddo
         enddo
      enddo
      if(ldir.or.lorient)then
         zmin=d(4)
         zmax=d(5)
         dmag=zmax-zmin
         if(.not.abs(dmag).gt.0)stop 'ERROR: surf3d d-limits'
      endif

      do k=1,(nx-1)*(ny-1)
c Search for the face not yet treated, furthest from eye
         dmin=1.e35
         do j=1,ny-1
            do i=1,nx-1
               if(work(i,j).lt.dmin)then
                  imin=i
                  jmin=j
                  dmin=work(i,j)
               endif
            enddo
         enddo
c         write(*,*)'imin,jmin,work,x,y',imin,jmin,work(imin,jmin),
c     $        x(imin,jmin),y(imin,jmin)
         do ic=0,4
            xp(1+ic)=x(imin+icx(ic),jmin+icy(ic))
            yp(1+ic)=y(imin+icx(ic),jmin+icy(ic))
            zp(1+ic)=z(imin+icx(ic),jmin+icy(ic))
            if(lorient)then
c Calculate the vertex direction.
c               write(*,*)'Calling vertexorient',ic,imin,icx(ic),jmin
c     $              ,icy(ic),imin+icx(ic),jmin+icy(ic)
               call vertexorient(x,y,z,iLx,nx,ny,imin+icx(ic),jmin
     $              +icy(ic),d,liper,ljper,forient)
c Use a separate vector, zz, for it.
               zz(1+ic)=forient
            endif
         enddo
c Calculate height of centroid.
         if(lchunk)then
            zcol=0.
            do ic=0,3
               zcol=zcol+z(imin+icx(ic),jmin+icy(ic))
            enddo
            zcol=zcol/4.
            icolor=int((240*(zcol-zmin))/(zmax-zmin))
            call gradcolor(icolor)
            call poly3line(xp,yp,zp,5)
            call pathfill()
         else
            if(ldir)then
               call gradquad(xp,yp,zp,d,zmin,zmax,0,239,3)
            elseif(lorient)then
               call gradquad(xp,yp,zp,zz,zmin,zmax,0,239,1)
            else
               call gradquad(xp,yp,zp,zp,zmin,zmax,0,239,1)
            endif
         endif
c Here there's a problem for the very last quadrilateral. 
c All the previous ones in postscript inherit the 0 setlinewidth
c that is called from the gradquad immediately afterwards. But the 
c last one does not have it, so it has a thicker line. One option
c would be to do a stroke before 0 setlinewidth, but that makes all the 
c lines thick (as usual). 
c I'm not sure if that's what I want. Trying it; seems better.
         call color(incolor)
         if(incolor.ne.0)call poly3line(xp,yp,zp,5)
c Finish at the last filled position. Necessary for glx driver.
         call vec3w(xp,yp,zp,0)
c Set face as done
         work(imin,jmin)=ddone
      enddo
      
      end
c********************************************************************
      subroutine surf1dr3(x,y,z,iLx,nx,ny,work,isw,d)
c Draw a 3-d surface of z(x,y), dim(iLx\nx,ny), using current scaling.
c This version is on a regular mesh x(nx),y(ny),
c given with 1-d vectors rather than 2D arrays.
c isw bit 0 set => triangular fills, else chunks. If in addition,
c isw bit 1 set => directional shading, optional argument d(5)
c    d(1-3) gives direction d(4-5) gives distance limits.
      integer nx,ny,iLx
      real x(nx),y(ny),z(iLx,ny),work(0:iLx+1,0:ny+1)
      integer isw
      real d(5)

      real x2,y2,z2
      real dout,ddone
      parameter (dout=2.e32,ddone=1.e32)
      real xp(5),yp(5),zp(5)
      integer icx(0:4), icy(0:4), ic
      logical lchunk,ldir
      data icx/0,1,1,0,0/
      data icy/0,0,1,1,0/
      save
c
c      write(*,*)'isw',isw
      lchunk=.true.
      if(isw-isw/2 .ne.0)lchunk=.false.
      ldir=.false.
      if(isw/2-isw/4 .ne.0)ldir=.true.
      incolor=igetcolor()
c Get eye position x0,y0,z0 are unused here.
      call trn32(x0,y0,z0,x2,y2,z2,-1)
c      write(*,*)'Eye=',x2,y2,z2
c Set up the work array.
      zmin=ddone
      zmax=-ddone
      do j=1,ny-1
         do i=1,nx-1
            if(z(i,j).gt.zmax)zmax=z(i,j)
            if(z(i,j).lt.zmin)zmin=z(i,j)
            work(i,j)=0.
            do ic=0,3
               work(i,j)=work(i,j)
     $              +x(i+icx(ic))*x2
     $              +y(j+icy(ic))*y2
     $              +z(i+icx(ic),j+icy(ic))*z2
            enddo
         enddo
      enddo
      if(ldir)then
         zmin=d(4)
         zmax=d(5)
         if(.not.abs(zmax-zmin).gt.0)stop 'ERROR: surf3d d-limits'
      endif

      do k=1,(nx-1)*(ny-1)
c Search for the face not yet treated, furthest from eye
         dmin=1.e35
         do j=1,ny-1
            do i=1,nx-1
               if(work(i,j).lt.dmin)then
                  imin=i
                  jmin=j
                  dmin=work(i,j)
               endif
            enddo
         enddo
         do ic=0,4
            xp(1+ic)=x(imin+icx(ic))
            yp(1+ic)=y(jmin+icy(ic))
            zp(1+ic)=z(imin+icx(ic),jmin+icy(ic))
         enddo
c Calculate height of centroid.
         if(lchunk)then
            zcol=0.
            do ic=0,3
               zcol=zcol+z(imin+icx(ic),jmin+icy(ic))
            enddo
            zcol=zcol/4.
            icolor=int((240*(zcol-zmin))/(zmax-zmin))
            call gradcolor(icolor)
            call poly3line(xp,yp,zp,5)
            call pathfill()
         else
            if(ldir)then
               call gradquad(xp,yp,zp,d,zmin,zmax,0,239,3)
            else
               call gradquad(xp,yp,zp,zp,zmin,zmax,0,239,1)
            endif
         endif
c Here there's a problem for the very last quadrilateral. 
c All the previous ones in postscript inherit the 0 setlinewidth
c that is called from the gradquad immediately afterwards. But the 
c last one does not have it, so it has a thicker line. One option
c would be to do a stroke before 0 setlinewidth, but that makes all the 
c lines thick (as usual). 
c I'm not sure if that's what I want. Trying it; seems better.
c         write(*,*)'incolor',incolor
         call color(incolor)
c         write(*,*)'incolor2',incolor
         if(incolor.ne.0)call poly3line(xp,yp,zp,5)
c Set face as done
         work(imin,jmin)=ddone
      enddo
      
      end

c ********************************************************************

c Surface direction at a vertex. (Actually dot product with direction).
c A vertex is surrounded by 4 adjacent points. The cross product between
c the vectors from the vertex to two adjacent adjacent points defines
c the normal to the facet they enclose and its area.  The dot product of
c that normal with the ref-direction defines the orientation relative to
c reference.  We regard the vertex face orientation as the average of
c the four facets. If some facets are bigger than others, it is sensible
c to accord them correspondingly greater weight in the average. So they
c remain weighted by their area
c
c Edges of the mesh need some special treatment. In some cases, periodic
c assumptions make sense (e.g. for a sphere or a torus). We take the 
c first and last mesh positions to be the same place. The whole surface
c must be covered between 1 and ni/j
c
c In other cases, an edge vertex has 3 adjacent vertexes defining 2
c facets. The best one can do linearly is average the two. A corner
c vertex has 1 facet.  This is accomplished by effectively setting
c adjacent vector length to zero if they project beyond the mesh.
c
c On input, x,y,z are 2-D (ni,nj) arrays of 3-D coordinates
c           iv,jv reference the vertex whose orientation is to be found
c           d     is a 3-D direction vector
c           lpi,lpi are two logicals true if mesh is periodic 
c On exit,  returned value forient is equal to the dot product of d
c           with the sum of the cross-products divided by the sum of
c           their magnitudes. (So if d is unit, forient is a cosine.)

      subroutine vertexorient(x,y,z,iLi,ni,nj,iv,jv,d,lpi,lpj,forient)
      real x(iLi,nj),y(iLi,nj),z(iLi,nj)
      real d(3)
      logical lpi,lpj
c Local variables
      real p(3,4),cross(3)
      logical lprint
      data lprint/.false./

c Get the adjacent point data into the p array +x, +y, -x, -y
 101  do k=1,4
         i=iv+mod(k,2)*(2-k)
         j=jv+mod(k+1,2)*(3-k)
         if(lpi)then
c Periodic x
            if(i.gt.ni)i=i-ni+1
            if(i.lt.1)i=i+ni-1
         else
c Open
            if(i.gt.ni)i=ni
            if(i.lt.1)i=1
         endif
         if(lpj)then
c Periodic y
            if(j.gt.nj)j=j-nj+1
            if(j.lt.1)j=j+nj-1
         else
            if(j.gt.nj)j=nj
            if(j.lt.1)j=1
         endif
         p(1,k)=x(i,j)-x(iv,jv)
         p(2,k)=y(i,j)-y(iv,jv)
         p(3,k)=z(i,j)-z(iv,jv)
         if(lprint)write(*,*)k,i,j,(p(m,k),m=1,3)
      enddo
c There's a problem with degenerate quadrilaterals, having coincident
c vertices. If such a vertex occurs at a mesh boundary, then there may
c only be one non-zero vector to adjacent vertices. 
c There ought to be a way to revert to the direction of the quadrilateral
c that is being analysed.

c Form facet areas and accumulate their projections.
      tcum=0.
      dcum=0.
      do k=1,4
         k1=mod(k,4)+1
c Cross product of adjacent vectors
         call crossprod(p(1,k),p(1,k1),cross,cmag)
c Dotted into the direction vector.
         ddot=0.
         do m=1,3
            ddot=ddot+cross(m)*d(m)
         enddo
c Accumulated for each facet
         dcum=dcum+ddot
         tcum=tcum+cmag
      enddo
      if(tcum.le.0.)then
         if(.not.lprint)then
            lprint=.true.
            goto 101
         endif
         write(*,*)'vertexorient error. tcum,dcum=',tcum,dcum
         write(*,*)'iv,jv,ni,nj',iv,jv,ni,nj
         write(*,'(3f8.4)')p
         write(*,*)'Probably a mesh/coordinate singularity'
         stop
      else
c And normalized
         forient=dcum/tcum
      endif
      if(lprint)lprint=.false.
      end
