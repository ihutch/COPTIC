c Vectorial Arrow-Plot
      subroutine arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,iswin,is1,is2)
c Draw arrows representing the field magnitude and direction of
c E1/2(imax/of/Li,jmax) on a domain x(imax[,jmax]),y([imax,]jmax).
c On entry the domain should have been set up so that x, y span it.
c Erange is the scaling parameter for arrow length, relative to 
c 1 unit of normalized position.
c If on entry Erange=0, calculate it automatically from the arrays.
c isw: 0: x,y unused; just index. 1: x,y vectors, 2: x,y matrices
c isw: if 4 is added to isw, use is1, is2 as the steps in x,y position
c      if not, then is1, and is2 arguments may be omitted in call.
c isw: if 8 is added to isw, limit the arrow arrays to narm in each direction
c      i.e. calculate is1 and is2 algorithmically internally.


      include 'plotcom.h'
      integer Li,imax,jmax,isw
      real E1(Li,jmax),E2(Li,jmax),x(*),y(*)
      character*20 string
      integer narm
      parameter (narm=20)

      isw=iswin
      if(isw/4-2*(isw/8).ne.0)then
         isl1=is1
         isl2=is2
      else
         isl1=1
         isl2=1
      endif
      if(isw/8-2*(isw/16).ne.0)then
         isl1=(imax-1)/narm+1
         isl2=(jmax-1)/narm+1
      endif
      isw=isw-4*(isw/4)

      if(Erange.eq.0)then
c Find Erange by scaling.
         call minmax2(E1,Li,imax,jmax,e1min,e1max)
         call minmax2(E2,Li,imax,jmax,e2min,e2max)
         e2max=max(abs(e2min),abs(e2max))*(jmax/float(isl2))
         e1max=max(abs(e1min),abs(e1max))*(imax/float(isl1))
c         write(*,*)'e1max=',e1max,'  e2max=',e2max
         Erange=max(e1max,e2max)
         if(Erange.eq.0.)then
            write(*,*)'Arrowplot scaling error. Null vector field.'
            return
         endif
         call fitrange(0.,Erange*.15,2,ipow,fac10,delta,first,xlast) 
c         write(*,*)Erange,delta,first,xlast,fac10
         csize=delta/Erange
         call charangl(90.)
         call fwrite(delta,iw,1,string)
         call jdrwstr(naxmax+chrshght+.01,naymin+0.1,string(1:iw),1.)
         call charsize(csize,0.3*csize)
         call drcstr('!A_!@')
      endif

      do j=1,jmax,isl2
         do i=1,imax,isl1
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
c            call jdrwstr(wx2nx(xc),wy2ny(yc),'!A!!',0.)
c            write(*,*)Ex,Ey,Emag,angle
         enddo
      enddo

c Restore defaults.
      call charsize(0.,0.)
      call charangl(0.)
      call drcstr('!@')
      end
c*******************************************************************
