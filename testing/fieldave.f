      program fieldave
c Calculates the Maxwell stress force by averaging over a series 
c of spheres round an object at 0,0 for testing.

      include 'examdecl.f'

      parameter (ntheta=4,nr=20)
      character*1 ppath(0:nr,na_k)
      character*20 string
      integer iLs(ndims+1)
      parameter (na_ij=na_i*na_j,na_ijk=na_ij*na_k)
      
      integer ns_nt,ns_np
c the size of the stress-calculating mesh in theta and psi directions
      parameter (ns_nt=10,ns_np=10)
      real surfobj(2*ndims,ns_nt,ns_np)
c
      real xc(ndims),rp(nr)
      parameter (ndiagmax=7)
c      real diagsum(na_i,na_j,na_k,ndiagmax)
      parameter (nplot=200)
      real xp(nplot),yp(nplot),xpt(ndims),xfield(ndims)
      real fieldforce(ndims)
      real charge

      real phimax
      data phimax/0./
      data rp/0.5,1.,(nr-2)*0./xc/0.,0.,0./
      data iLs/1,na_i,na_ij,na_ijk/

c 
      string=''
      ppath(0,1)=''
      diagfilename=''
      isw=1
      charge=1.
c      write(*,*)rp
      call examargs(xc,rp,phimax,isw,charge)
      charge=charge*4.*3.14159

      ied=1
      call array3read(phifilename,ifull,iuds,ied,u,ierr)
      if(ierr.eq.1)stop

      write(*,*)'iLs',iLs
      write(*,'(a,3i4,$)')'On grid',iuds
      write(*,'(a,2f7.2,a,2f7.2,a,2f7.2)')
     $     (',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)

      ifix=1
c      call noeye3d(0)
      call sliceGweb(ifull,iuds,u,na_m,zp,
     $     ixnp,xn,ifix,'potential:'//'!Af!@')
c      call sliceGcont(ifull,iuds,u,na_m,zp,
c     $     ixnp,xn,ifix,'potential:'//'!Af!@')

c Initialize the force tracking for one object with center xc and 
c radius r. With angular mesh ns_nt, ns_np.
      do i=1,nr
         if(rp(i).eq.0)goto 102
         call forceinitone(ns_nt,ns_np,ndims,rp(i),xc,surfobj)
         call chargeforce(ndims,ns_nt*ns_np
     $        ,surfobj,fieldforce,u,iLs)
c Convert to units of nTr^2:
         do k=1,ndims
            fieldforce(k)=charge*fieldforce(k)*debyelen**2
         enddo
         write(*,'(a,f10.4,a,3f12.6)')'r=',rp(i)
     $        ,' Fieldforce',fieldforce
      enddo
 102  continue

c Simple plot of the field along a z-chord.
      xrange=1.
      xpt(1)=0.
      xpt(2)=0.
      do i=1,nplot
         xp(i)=-xrange+i*2.*xrange/nplot
         xpt(3)=xp(i)
         call fieldsimple3atpoint(xpt,u,iLs,xfield)
c         write(*,*)xfield,xpt
         yp(i)=xfield(3)
      enddo
      call autoplot(xp,yp,nplot)
      call axlabels('z','z-Field')
      call pltend()

      call exit(0)
 101  write(*,*)'Error writing output'
      end

c*************************************************************
      subroutine examargs(xc,rp,phimax,isw,charge)
      include 'examdecl.f'
      real xc(ndims)
      parameter (nr=20)
      real rp(nr)
      real charge

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

c Defaults and silence warnings.
      phifilename=' '
      zp(1,1,1)=0.

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
         if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)',err=201)phifilename
         if(argument(1:2).eq.'-d')
     $        read(argument(3:),'(a)',err=201)diagfilename
         if(argument(1:2).eq.'-q')
     $        read(argument(3:),*,err=201)charge
         if(argument(1:2).eq.'-c')
     $        read(argument(3:),*,err=201)xc
         if(argument(1:2).eq.'-r')then
            read(argument(3:),*,end=210,err=201)rp
 210        write(*,*)'Assigned radii',rp
         endif
         if(argument(1:2).eq.'-p')
     $        read(argument(3:),'(f8.4)',err=201)phimax
         if(argument(1:2).eq.'-l')
     $        read(argument(3:),*,err=201)isw
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)phifilename
         endif
         
      enddo
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5)
      write(*,301)'Usage: fieldave [switches] <phifile>'
c      write(*,301)' --objfile<filename>  set name of object data file.'
c     $     //' [copticgeom.dat'
      write(*,301)
     $     ' -d<diagfilename>  -r<rp1,rp2,...> -p<phimax> -l<isw>'
      write(*,301)' isw: Byte 1: gradlegend(1) 2: ...'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
      end
