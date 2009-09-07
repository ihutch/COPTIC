c**********************************************************************
c General version constructs a mesh from xmstart(id) to xmend(id)
      subroutine meshconstruct(ndims,iuds,xmstart,xmend)
      integer iuds(ndims)
      real xmstart(ndims),xmend(ndims)
      include 'meshcom.f'

      iof=0
      do id=1,ndims 
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data. Make it extend from -rs to +rs
         do i=1,iuds(id)
c            xn(i+iof)=rs*(2*(i-1.)/(iuds(id)-1.) - 1.)
c            xn(i+iof)=abs(xmstart(1))*(2*(i-1.)/(iuds(id)-1.) - 1.)
            xn(i+iof)=((iuds(id)-i)*xmstart(id)+
     $           (i-1.)*xmend(id))/(iuds(id)-1.)
c            write(*,*)xmstart(id),xmend(id),xn(i+iof)
c Two region Non-uniform assuming iuds even.
c            iu2=iuds(id)/2
c Uniform -1 to 1:
c            xi=(i-iu2-0.5)/(iu2-0.5)
c Half of mesh at each spacing
c            xmid=0.5
c First half of mesh extends to this fraction of the total rs.
c            xs=0.25
c            if(abs(xi).lt.xmid)then
c               xn(i+iof)=rs*xs*xi/xmid
c            else
c               xn(i+iof)=rs*sign(((abs(xi)-xmid)+xs*(1.-abs(xi))),xi)
c     $              /(1-xmid)
c            endif
         enddo
         iof=iof+iuds(id)
      enddo
      ixnp(ndims+1)=iof

      end
