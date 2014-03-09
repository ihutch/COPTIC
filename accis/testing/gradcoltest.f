c***********************************************************************
      program gradcoltest
      integer nx,ny
      parameter (nx=100,ny=100)
      real z(nx,ny)

      do j=1,ny
         do i=1,nx
            z(i,j)=(j-1.)/(ny-1.)
         enddo
      enddo

      call pfset(3)

      call blueredgreenwhite()
      call autocolcont(z,nx,nx,ny)

      call pltend()

      end
