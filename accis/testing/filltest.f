c Demo of filling and line simultaneously
      include 'plotcom.h'
      integer nx, ny
      parameter (nx=16,ny=16)
      real x(nx),y(ny)
      real xp(5),yp(5)
      integer idx(5),idy(5)
      data idx/0,1,1,0,0/
      data idy/0,0,1,1,0/

      do i=1,nx
         x(i)=i
      enddo
      do j=1,ny
         y(j)=j
      enddo

      call pfset(3)
      call pltinit(x(1),x(nx),y(1),y(ny))
      call accisgradinit(-40000,12000,-40000,64000,64000,64000)
      do j=1,ny-1
         do i=1,nx-1
            do id=1,5
               xp(id)=x(i+idx(id))
               yp(id)=y(j+idy(id))
            enddo
            icol=(i+(j-1)*nx)
            if(icol.le.15) then
               call color(icol)
            else
               call gradcolor(icol)
            endif
            call polyline(xp,yp,5)
            call pathfill()
c Second write is necessary to get the line.
            call color(0)
            call polyline(xp,yp,5)
         enddo
      enddo
      call pltend()
      end
