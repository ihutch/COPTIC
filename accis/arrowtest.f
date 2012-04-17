      program arrowtest
      parameter (Li=200,imax=20,jmax=50)
      real E1(Li,Li),E2(Li,Li),x(Li),y(Li)

      do i=1,imax
         x(i)=i
      enddo
      do j=1,jmax
         y(j)=j
      enddo
      do i=1,imax
         do j=1,jmax
            r=sqrt(x(i)**2+y(j)**2)
            theta=atan2(y(j),x(i))
            Emag= sin(2*theta)**2*r/(imax+(r**2/(.2*imax)))
            E1(i,j)=-Emag*sin(theta)
            E2(i,j)= Emag*cos(theta)
         enddo
      enddo

      call pfset(3)
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,imax/2.,Li,imax,jmax,x,y,1,idum,idum) 
      call pltend()

      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,imax/2.,Li,imax,jmax,x,y,5,4,4) 
      call pltend()

      Erange=0.
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,5,3,6) 
      write(*,*)'Scaling Erange=',Erange
      call pltend()

      Erange=0.
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,9,idum,idum) 
      call pltend()

      end
