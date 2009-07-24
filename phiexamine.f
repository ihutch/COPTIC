      program phiexamine

      include 'examdecl.f'
c 
      call examargs

c      partfilename=phifilename(1:lentrim(phifilename)-4)
c      nb=nbcat(partfilename,'.flx')
c      myid=0
c      nb=nameappendint(partfilename,'.',myid,3)
c         call readfluxfile(fluxfilename,ierr)
c         if(ierr.ne.0)goto 401
c         call partread(partfilename,ierr)
c         if(ierr.ne.0)goto 401

      call phiread(phifilename,ifull,iuds,u,ierr)
      call sliceGweb(ifull,iuds,u,Li,zp,
     $        ixnp,xn,ifix,'potential:'//'!Ay!@')

c plot potential versus radius.

      write(*,*)rs
      call pltinit(0.,rs,-2.,0.)
      call axis()
      call axlabels('radius','potential')
      call charsize(.001,.001)
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               x=xn(ixnp(1)+i)
               y=xn(ixnp(2)+j)
               z=xn(ixnp(3)+k)
               r=sqrt(x**2+y**2+z**2)
               call polymark(r,u(i,j,k),1,10)
               if(r.gt.rs .and. u(i,j,k).ne.0)then
                  write(*,'(4f12.6,3i3)')x,y,z,u(i,j,k),i,j,k
               endif
            enddo
         enddo
      enddo
      call charsize(0.,0.)
      call pltend()

      end


c*************************************************************
      subroutine examargs()
      include 'examdecl.f'

      do i=1,ndims
         ifull(i)=Li
      enddo

c Defaults
      phifilename=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
         if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)',err=201)phifilename
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
 302  format(a,f8.3)
      write(*,301)'Usage: phiexamine [switches] <phifile>'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
c      write(*,301)' -f   set name of phifile.'
c      write(*,301)'Debugging switches for testing'
c      write(*,301)' -gc   set wireframe/stencils(-) mask.'//
c     $     ' objects<->bits. [',iobpl
c      write(*,301)' -gt   Plot solution tests.'
c      write(*,301)' -gs   Plot slices of solution potential. '
c      write(*,301)' -go   set No of orbits'
c     $     //'(to plot on objects set by -gc). [',norbits
c      write(*,301)' -at   set test angle.'
c     $     //' -an   set No of angles. '
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
c*****************************************************************
