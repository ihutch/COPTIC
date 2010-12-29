c Vectorial Arrow-Plot
      subroutine arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,iswin,is1,is2)
c Draw arrows representing the field magnitude and direction of
c E1/2(imax/of/Li,jmax) on a domain x(imax[,jmax]),y([imax,]jmax).
c On entry the domain should have been set up so that x, y span it.
c Erange is the scaling parameter for arrow length, relative to 
c 1 unit of normalized position.
c isw: 0: x,y unused; just index. 1: x,y vectors, 2: x,y matrices
c isw: if 4 is added to isw, use is1, is2 as the steps in x,y position

      integer Li,imax,jmax,isw
      real E1(Li,jmax),E2(Li,jmax),x(*),y(*)

      isw=iswin
      if(isw/4-2*(isw/8).ne.0)then
         isl1=is1
         isl2=is2
      else
         isl1=1
         isl2=1
      endif
      isw=isw-4*(isw/4)

      do j=1,jmax,isl1
         do i=1,imax,isl2
            Ex=E1(i,j)
            Ey=E2(i,j)
            angle=atan2(Ey,Ex)
            Emag=sqrt(Ex**2 + Ey**2)
            csize=Emag/Erange+1.e-6
            if(isw.eq.0)then
               xc=i
               yc=j
            elseif(isw.eq.1)then
               xc=x(i)
               yc=y(j)
            elseif(isw.eq.2)then
               xc=x(i+Li*(j-1))
               yc=y(i+Li*(j-1))
            endif
            call charangl(180.*angle/3.1415926)
            call charsize(csize,0.3*csize)
            call jdrwstr(wx2nx(xc),wy2ny(yc),'!A_',0.)
c            write(*,*)Ex,Ey,Emag,angle
         enddo
      enddo

c Restore defaults.
      call charsize(0.,0.)
      call charangl(0.)
      call drcstr('!@')
      end
c*******************************************************************
