      program threadtest
      parameter (nx=100)
      real x(nx),y(nx)
      character*(100) input

      do i=1,nx
         x(i)=i
         y(i)=cos(.2*i)
      enddo

      do i=1,100
         call autoplot(x,y,nx)
         call iwrite(i,iwidth,input)
         call jdrwstr(.5,.5,'Instance'//input(1:iwidth),0.)
c This is still needed:         
         call accisflush()
         write(*,*)'Now waiting on input'

         read(*,*)input

      enddo

      end
