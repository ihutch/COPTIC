c*************************************************************
      program test3d
      character*100 argument,filename
      parameter (imax=100,jmax=100,kmax=10)
      real x(3,imax,jmax,kmax)
      real u(imax,jmax,kmax)
      parameter (x1s=1.,x2s=1.,x3s=1.)
      real xp(3)
      
      filename='geominput.dat'
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:2).eq.'-f') read(argument(3:),'(a)')filename
         if(argument(1:2).eq.'-p') read(argument(3:),*)xp
      enddo

      call readgeom(filename)

      is=inside_geom(3,xp,1)
      inside=insideall(3,xp)
      write(*,*)'xp=',xp,' insideall=',inside,' is=',is

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               x(1,i,j,k)=x1s*(i-1)/(imax-1.)
               x(2,i,j,k)=x2s*(j-1)/(jmax-1.)
               x(3,i,j,k)=x3s*(k-1)/(kmax-1.)
               u(i,j,k)=insideall(3,x(1,i,j,k))
c               if(u(i,j,k).ne.0) write(*,*)i,j,k
            enddo
         enddo
         if(imax.le.10)write(*,'(10f5.0)')((u(i,j,k),i=1,imax),j=1,jmax)
         call autocolcont(u(1,1,k),imax,imax,imax) 
c         call gradlegend(0.,1.,-.1,0.,-.1,1.,.02,.false.)
         call pltend()
      enddo

      

      end
