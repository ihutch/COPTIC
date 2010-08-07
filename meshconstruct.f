c**********************************************************************
c Version constructs a mesh from xmstart(id) to xmend(id)
c Convert from previously set imeshstep, xmeshpos to ixnp.
c There must be at least 2 non-zero and monotonically increasing
c values starting imeshstep.
c Return values of iuds(ndims) per constructed mesh.
c Also initialize the value of rs in plascom. Put equal to
c half the largest mesh box side length.
      subroutine meshconstruct(ndims,iuds)
      integer iuds(ndims)
      include 'meshcom.f'
      include 'plascom.f'

      iof=0
      rs=0.
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
         xmeshstart(id)=xmeshpos(id,1)
         xn(iof+1)=xmeshpos(id,1)
         iof=iof+1
         do iblk=1,nspec_mesh-1
            msteps=imeshstep(id,iblk+1)-imeshstep(id,iblk)
            if(msteps.le.0) goto 11
c Mesh data.
            do i=1,msteps
               xn(i+iof)=((msteps-i)*xmeshpos(id,iblk)+
     $              (i)*xmeshpos(id,iblk+1))/(msteps)
            enddo
            iof=iof+msteps
         enddo
         iblk=nspec_mesh-1
 11      continue
         if(iblk.le.1) write(*,*)'Too few mesh steps. Dimension',id
         xmeshend(id)=xmeshpos(id,iblk)
         rsi=0.5*abs(xmeshend(id)-xmeshstart(id))
         if(rsi.gt.rs)rs=rsi
c Set iuds according to specified mesh
         iuds(id)=imeshstep(id,iblk)
c         write(*,'(a,i3,10f8.3)')
c     $        ' Meshspec',id,(xmeshpos(id,kk),kk=1,iblk)
      enddo
c      write(*,*)'Meshcontructed',xmeshstart,xmeshend,' rs=',rs
      ixnp(ndims+1)=iof
c      write(*,*)(ixnp(k),k=1,ndims+1)
c      write(*,*)(xn(k),k=1,iof)
      end
c***************************************************************
c Version constructs a mesh from xmstart(id) to xmend(id)
      subroutine meshconstructold(ndims,iuds,xmstart,xmend)
      integer iuds(ndims)
      real xmstart(ndims),xmend(ndims)
      include 'meshcom.f'

      iof=0
      do id=1,ndims
         xmeshstart(id)=xmstart(id)
         xmeshend(id)=xmend(id)
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data.
         do i=1,iuds(id)
            xn(i+iof)=((iuds(id)-i)*xmstart(id)+
     $           (i-1.)*xmend(id))/(iuds(id)-1.)
c            write(*,*)xmstart(id),xmend(id),xn(i+iof)
c Non-uniform needs specific programming.
         enddo
         iof=iof+iuds(id)
      enddo
      ixnp(ndims+1)=iof

      end
