c*********************************************************************
c This is just stuff for testing, with replacement getfield/potentials.
c********************************************************************
c test program
      program sphtest

      parameter (ntmax=100,npmax=100,ndims=3)
      real surfobj(2*ndims*ntmax*npmax)
      real force(ndims)


      nt=80
      np=80

c the call      
      call spheremesh(nt,np,surfobj)

c the tests:
 1    call pltinit(0.,1.,0.,1.)
      call setcube(0.2,0.2,0.2,.5,.4)
c      call scale3(-1.,1.,-1.,1.,-1.,1.)
c      call geteye(x2,y2,z2)
c      call trn32(0.,0.,0.,x2,y2,z2,1)
      call axproj(igetcorner())
      call ax3labels('x','y','z')
      call cubeproj()

      scale=nt*np/100.
      do i=1,nt*np
         ipt=6*(i-1)
c         write(*,'(6f8.4)')(surfobj(ipt+k),k=1,6)
         call vec3w(surfobj(ipt+1),surfobj(ipt+2),surfobj(ipt+3),0)
         call vec3w(surfobj(ipt+1)+surfobj(ipt+4)*scale
     $        ,surfobj(ipt+2)+surfobj(ipt+5)*scale
     $        ,surfobj(ipt+3)+surfobj(ipt+6)*scale
     $        ,1)
c         write(*,'(3f8.4)')(surfobj(ipt+k),k=1,3)
      enddo
      if(ieye3d().ne.0)goto 1

c Add up areas
      areax=0.
      areay=0.
      areaz=0.
      area=0.
      do i=1,nt*np
         ipt=6*(i-1)
         if(surfobj(ipt+4).gt.0.)areax=areax+surfobj(ipt+4)
         if(surfobj(ipt+5).gt.0.)areay=areay+surfobj(ipt+5)
         if(surfobj(ipt+6).gt.0.)areaz=areaz+surfobj(ipt+6)
         area=area+sqrt(
     $        surfobj(ipt+4)**2+
     $        surfobj(ipt+5)**2+
     $        surfobj(ipt+6)**2)
      enddo

      write(*,'(a,f10.6,a,f10.6,a,f10.6,a,f10.5)')
     $     'areax=',areax,' areay=',areay,' areaz=',areaz,' area=',area

      write(*,*)'4pi=',4.*3.141593,' ratio',area/(4.*3.1415926)

c It is noticeable that these sum tests give good results only if the
c number of cells in the mesh is even. 

      km=nt*np
      call maxwellforce(ndims,km,surfobj,force,  u,cij,iLs)
      write(*,'(a,3f8.4)')'Maxwellforce=',(force(k),k=1,ndims)

      call pressureforce(ndims,km,surfobj,force,  u,cij,iLs)
      write(*,'(a,3f8.4)')'Pressureforce=',(force(k),k=1,ndims)

      call totalcharge(ndims,km,surfobj,charge,  u,cij,iLs)
      write(*,'(a,3f8.4)')'Totalcharge=',charge

      end
c********************************************************************
c Version of field at point for testing.
c Just use analytic formula for field at point x from charge,dipole. 
      subroutine fieldatpoint(x,u,cij,iLs,field)

      real x(3),x0(3),field(3),u,cij
      integer iLs
      real dipole(3)

c charge at x0
      charge=1.
      x0(1)=0.
      x0(2)=.7
      x0(3)=0.
c dipole at origin
      dipole(1)=0.
      dipole(2)=0.
      dipole(3)=0.
c External field in x direction
      extfield=1.

      field(1)=extfield
      field(2)=0.
      field(3)=0.
c Add on the rest.
      r=0.
      rd=0
      do i=1,3
         r=r+(x(i)-x0(i))**2
         rd=rd+x(i)**2
      enddo
      r=sqrt(r)
      rd=sqrt(rd)
      do i=1,3
         field(i)=field(i)+charge*(x(i)-x0(i))/(r*r**2)/(4.*3.141596)
         sprod=0.
         do k=1,3
            sprod=sprod+dipole(k)*x(k)
         enddo
         field(i)=field(i)
     $        +(3.*sprod*x(i)/rd**2 - dipole(i))/(rd*rd**2)
     $        /(4.*3.141596)
      enddo

      end
c********************************************************************
c Version of potential at point for testing.
c Just use analytic formula for potl at point x from charge,dipole. 
      real function potentialatpoint(x,u,cij,iLs)

      real x(3),x0(3),u,cij
      integer iLs
      real dipole(3)


c charge at x0
      charge=1.
      x0(1)=0.
      x0(2)=0.
      x0(3)=0.
c dipole at origin
      dipole(1)=0.
      dipole(2)=0.
      dipole(3)=0.
c External field in x direction
      extfield=0.

      potential=extfield*x(1)
c Add on the rest.
      r=0.
      rd=0
      do i=1,3
         r=r+(x(i)-x0(i))**2
         rd=rd+x(i)**2
      enddo
      r=sqrt(r)
      rd=sqrt(rd)
      do i=1,3
         potential=potential+charge/(4.*3.141596*r)
         sprod=0.
         do k=1,3
            sprod=sprod+dipole(k)*x(k)
         enddo
         potential=potential+sprod/(2.*rd*rd**2)
      enddo

      potentialatpoint=potential

      end
