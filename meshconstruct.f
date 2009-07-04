c**********************************************************************
      subroutine meshconstruct(ndims,iuds,ium2,rs)
      integer iuds(ndims),ium2(ndims)
      include 'meshcom.f'

      iof=0
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data. Make it extend from -rs to +rs
         do i=1,iuds(id)
            xn(i+iof)=rs*(2*(i-1.)/(iuds(id)-1.) - 1.)
            if(.false.)then
c Two region Non-uniform assuming iuds even.
            iu2=iuds(id)/2
c Uniform -1 to 1:
            xi=(i-iu2-0.5)/(iu2-0.5)
c Half of mesh at each spacing
            xmid=0.5
c First half of mesh extends to this fraction of the total rs.
            xs=0.25
            if(abs(xi).lt.xmid)then
               xn(i+iof)=rs*xs*xi/xmid
            else
               xn(i+iof)=rs*sign(((abs(xi)-xmid)+xs*(1.-abs(xi))),xi)
     $              /(1-xmid)
            endif
            endif
         enddo
         iof=iof+iuds(id)
c Mesh size with the faces removed:
         ium2(id)=iuds(id)-2
      enddo
      ixnp(ndims+1)=iof

      end
c*********************************************************************
      subroutine text3graphs(ndims,iuds,ifull,cij,volumes)
      integer iuds(ndims),ifull(ndims)
      real cij(*)
      real volumes(*)
      include 'meshcom.f'
      character*40 form1
c Standard volume for uniform mesh:
         vs=.99999
         do id=1,ndims
            vs=vs*(xn(ixnp(id)+2)-xn(ixnp(id)+1))
         enddo
         vs=vs
c Text graphic of slice through volumes
         write(*,*)'Volumes percentages:'
         write(form1,'(''('',i2,''i4)'')')iuds(1)
         write(*,form1)((nint(100.*
     $        volumes(j-1+(iuds(2)/2-1)*ifull(1)
     $                +(k-1)*ifull(1)*ifull(2)+1)/vs),
     $        j=1,iuds(1)),k=1,iuds(3))

c Text graphic of slice through cij
c      write(*,*)'iregion:'
         write(form1,'(''('',i2,''i1)'')')iuds(1)
         write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),
     $        j=1,iuds(1)),k=1,iuds(3))
      end
