      program partexamine
c Examine the particle data, showing distribution function(s) for 
c selected cell ranges.
      include 'examdecl.f'
      parameter (nfilemax=999)
      include 'partcom.f'
 
      character*10 chartemp
      character*100 name
      character*50 string

      include 'ptaccom.f'
c Distributions
c      parameter (ndiag=100,mdims=3)
c      real xr(3*mdims)
c      real fv(ndiag,mdims)
c      real px(ndiag,mdims)
c      real diagv(ndiag)
c      real diagx(ndiag,mdims)
c      common /cartdiag/fv,px,diagv,diagx

c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims)

      do id=1,mdims
         xlimit(1,id)=-5.
         xlimit(2,id)=5.
      enddo
c      xlimit(2,3)=-4.5

      call partexamargs(xlimit)
c Now the base filename is in partfilename.

      do i=0,nfilemax
         write(chartemp,'(''.'',i3.3)')i
         name=partfilename(1:lentrim(partfilename))//chartemp
         write(*,*)'Reading file ',name(1:lentrim(name))
         call partread(name,ierr)
         if(ierr.ne.0)goto 1
c Do the accumulation for this file up to maximum relevant slot. 
         do j=1,ioc_part
c Only for filled slots
            if(if_part(j).eq.1)then
               do id=1,mdims
                  x=x_part(id,j)
                  if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 2
               enddo
               call partaccum(x_part(1,j))
            endif
 2          continue
         enddo
      enddo
 1    continue

      call multiframe(2,1,2)
      do id=1,mdims
         call ticnumset(10)
         call autoplot(diagv,fv(1,id),ndiag)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('velocity',string(1:lentrim(string)))
c         call pltend()
         call ticnumset(10)
         call autoplot(diagx(1,id),px(1,id),ndiag)
         call axlabels('position',string(1:lentrim(string)))
         call pltend()
      enddo

      end
c*************************************************************
      subroutine partexamargs(xlimit)
      include 'examdecl.f'
      real xlimit(2,3)

      do i=1,ndims
         ifull(i)=Li
      enddo
c Defaults
      partfilename=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,err=201)
     $              xlimit(1,1),xlimit(2,1)
            endif
            if(argument(1:2).eq.'-y')then
               read(argument(3:),*,err=201)
     $              xlimit(1,2),xlimit(2,2)
            endif
            if(argument(1:2).eq.'-z')then
               read(argument(3:),*,err=201)
     $              xlimit(1,3),xlimit(2,3)
            endif
            if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-f')
     $           read(argument(3:),'(a)',err=201)partfilename
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)partfilename
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
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' -x -y -z<fff,fff>  set position range. [ -5,5'
      write(*,301)' -f   set name of partfile.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
c*****************************************************************
      subroutine partaccum(xr)
c Test the reinjection scheme by forming cartesian distributions.
      include 'ptaccom.f'
      include 'plascom.f'
      include 'meshcom.f'
c      include 'creincom.f'
      character*100 string
      logical lfirst
      data lfirst/.true./
      save lfirst


      vrange=5.
c Initialization
      if(lfirst)then
c Default mesh data.
         do id=1,mdims
            xmeshstart(id)=-5.
            xmeshend(id)=5.
         enddo
         
         do i=1,ndiag
            diagv(i)=vrange*(-1.+2.*(i-0.5)/ndiag)
            do id=1,mdims
               fv(i,id)=0.
               px(i,id)=0.
               diagx(i,id)=xmeshstart(id)+(i-0.5)*
     $              (xmeshend(id)-xmeshstart(id))/(ndiag)
            enddo
         enddo
         lfirst=.false.
      endif

c Assign velocities to bins.
      do id=1,mdims
         v=xr(mdims+id)
         v=sign(min(vrange,abs(v)),v)
         ibin=nint(0.5*(ndiag)*(1.+0.99999*v/vrange)+0.5)
         if(ibin.lt.1.or.ibin.gt.ndiag)
     $        write(*,*)k,nin,id,' ibin',ibin,v
         fv(ibin,id)=fv(ibin,id)+1
c Assign positions to bins
         x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
         ibin=nint(0.50000+x*(ndiag-.00000))
         if(ibin.lt.1 .or. ibin.gt.ndiag)then
            write(*,*)'ibin=',ibin,id,x,xr(id)
         else
            px(ibin,id)=px(ibin,id)+1.
         endif
      enddo


      end
c**********************************************************************
