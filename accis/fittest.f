c  Test fitrange.
      parameter (nbmax=20,ntmax=20)


      do itics=1,10
      do i=1,nbmax
         do j=1,ntmax
            xmin=i*10./nbmax
            xmax=j*10./ntmax
            if(xmin.ne.xmax)then
               call fitrange(xmin,xmax,itics,ipow,fac10,delta,first
     $              ,xlast) 
               write(*,'(2f6.3,i3,i4,4f7.3)')xmin,xmax,itics,ipow,fac10
     $              ,delta,first
            endif
         enddo
      enddo
      enddo

      end
