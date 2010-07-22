      program philineread

      include 'examdecl.f'

      character*100 filenames(Li)
      character*50 string
      integer ild,ilinechoice(ndims_mesh,20),ip(ndims_mesh)
      real philine(Li),xline(Li)
      real darray(20),pmax(20),punscale(20),rp(20)

      logical lrange,lwrite
      common /linecom/xmin,xmax,ymin,ymax,lrange,lwrite
c 
c silence warnings:
      zp(1,1,1)=0.
      fluxfilename=' '
      xmin=0.
      ymin=0.
      xmax=0.
      ymax=0.
      lwrite=.false.

c      write(*,*)nf,idl
      call lineargs(filenames,nf,ild,ilinechoice,rp)
c      write(*,*)'Filenames',nf

      nplot=0
      do inm=1,nf
         phifilename=filenames(inm)
         call array3read(phifilename,ifull,iuds,u,ierr)
         if(ierr.eq.1)stop
         do i=1,ndims_mesh
            if(iuds(i).gt.Li) then
               write(*,*)'Data too large for Li',Li,iuds(i)
               stop
            endif
         enddo

c Select the lineout into the plotting arrays.      
         if(ild.ne.0)then
            do k=1,ndims_mesh-1
               ip(mod(ild+k-1,ndims_mesh)+1)=
     $              ilinechoice(mod(ild+k-1,ndims_mesh)+1,inm)
c               write(*,*)k,ilinechoice(mod(ild+k-1,ndims_mesh)+1,inm)
            enddo
c            write(*,*)'ild,ip(ild)',ild,ip(ild)
            do i=1,iuds(ild)
               ip(ild)=i
               xline(i)=xn(ixnp(ild)+i)/debyelen
c               philine(i)=u(ip(1),ip(2),ip(3))
c Assuming the probe to be of radius 1 unit, normalized phi:
c               philine(i)=u(ip(1),ip(2),ip(3))/abs(phip/debyelen)
               philine(i)=u(ip(1),ip(2),ip(3))
     $              /abs(rp(inm)*phip/debyelen)
               if(lwrite)write(*,'(2f10.5)')xline(i),philine(i)
            enddo
            call minmax(philine,iuds(ild),pmin,pa)
            pmax(inm)=pa
            punscale(inm)=pmax(inm)*abs(rp(inm)*phip/debyelen)
            darray(inm)=abs(rp(inm)*phip/debyelen)

            call winset(.true.)
            call pfset(3)
            if(inm.eq.1)then
               if(lrange)then
                  if(xmin-xmax.eq.0.)then
                     call minmax(xline,iuds(ild),xmin,xmax)
                  endif
                  if(ymin-ymax.eq.0.)then
                     call minmax(philine,iuds(ild),ymin,ymax)
                  endif
                  call pltinit(xmin,xmax,ymin,ymax)
                  call axis()
                  call axlabels('z/!Al!@',
     $                 '!Af!@/(|!Af!@!dp!d|r!dp!d/!Al!@)')
                  call winset(.true.)
                  call polyline(xline,philine,iuds(ild))
               else
                  call autoplot(xline,philine,iuds(ild))
               endif
               call axis2()
            else
               call color(inm)
c               call iwrite(inm,iwd,string)
c               call labeline(xline,philine,iuds(ild),string,iwd)
               call dashset(inm)
               call polyline(xline,philine,iuds(ild))
            endif
            string=' !Af!@!dp!d='
            call fwrite(phip,iwd,2,string(lentrim(string)+1:))
            string(lentrim(string):)='@'
            call fwrite(rp(inm)/debyelen,iwd,2,
     $           string(lentrim(string)+1:))
            call legendline(.5,(.01+inm*.05),0,
     $           string(1:lentrim(string)))
            nplot=nplot+1
         endif

         write(*,'(a,3i4,$)')'On grid',iuds
         write(*,*)(',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)
     $        ,ip
c         write(*,*)ild,ilinechoice
      enddo
      if(nplot.gt.0)then
         call pltend()
         call dashset(0)
         call lautomark(darray,pmax,nplot,.true.,.true.,0)
         do k=1,nplot
c            write(*,*)nint(rp(k)/.2)
            call polymark(darray(k),pmax(k),1,nint(rp(k)/.2))
         enddo
         call axlabels('|!Af!@!dp!d|r!dp!d/!Al!@'
     $        //'!A ~ !@Q/4!Ape!@!d0!d!Al!@',
     $        '!Af!@!dmax!d/(|!Af!@!dp!d|r!dp!d/!Al!@)')
         call vecw(0.04,3.,0)
         call vecw(1.,.12,1)
         call vecw(0.01,2.,0)
         call vecw(0.1,2.,1)
         call pltend()
         call lautomark(darray,punscale,nplot,.true.,.true.,0)
         do k=1,nplot
c            write(*,*)nint(rp(k)/.2)
            call polymark(darray(k),punscale(k),1,nint(rp(k)/.2))
         enddo
         write(*,*)'punscale',(punscale(k),k=1,nplot)
         call axlabels('|!Af!@!dp!d|r!dp!d/!Al!@'
     $        //'!A ~ !@Q/4!Ape!@!d0!d!Al!@',
     $        '!Af!@!dmax!d')
         call pltend()
      endif


      end

c*************************************************************
      subroutine lineargs(filenames,nf,ild,ilinechoice,rp)
      include 'examdecl.f'
      character*100 filenames(Li)
      real rp(20)
      integer nf,ild,ilinechoice(ndims_mesh,20)
      integer idj(ndims_mesh)

      logical lrange,lwrite
      common /linecom/xmin,xmax,ymin,ymax,lrange,lwrite

      do i=1,ndims
         ifull(i)=Li
      enddo

c Defaults and silence warnings.
      phifilename=' '
      fluxfilename=' '
      nf=1
      zp(1,1,1)=0.
      ild=0
      lrange=.false.
      do id=1,ndims_mesh
         idj(id)=0
         ilinechoice(id,1)=0
      enddo
      rread=1.
      

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:2).eq.'-l')then
               read(argument(3:),*,end=11)ild,(idj(k),k=1,ndims_mesh-1)
 11            continue
c               write(*,*)ild, idj
               do k=1,ndims_mesh-1
                  ilinechoice(mod(ild+k-1,ndims_mesh)+1,nf)=idj(k)
c                  write(*,*)'nf,k,idj(k)',nf,k,idj(k)
               enddo
            endif
            if(argument(1:2).eq.'-y')then
               read(argument(3:),*,end=12)ymin,ymax
               lrange=.true.
 12            continue
            endif
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,end=13)xmin,xmax
               lrange=.true.
 13            continue
            endif
            if(argument(1:2).eq.'-w')lwrite=.true.

            if(argument(1:13).eq.'--objfilename')
     $           read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-r')
     $           read(argument(3:),'(f8.4)',err=201)rread
            if(argument(1:2).eq.'-f')
     $           read(argument(3:),'(a)',err=201)phifilename
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)phifilename
            filenames(nf)(1:)=phifilename
            rp(nf)=rread
            nf=nf+1
         endif
      enddo
      nf=nf-1

c      write(*,*)ild,
c     $     (ilinechoice(2,i),' ',rp(i),' ',filenames(i)(1:30),i=1,nf)
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5)
      write(*,301)'Usage: philineread [switches] <phifile>'//
     $     ' [<phifile2> ...]'
c      write(*,301)' --objfile<filename>  set name of object data file.'
c     $     //' [ccpicgeom.dat'
      write(*,301)' -l<idim>,<irow1>,<irow2> set fixed dimension [',ild
      write(*,301)'    and row position in other dimensions.'
      write(*,301)' -y<min>,<max>   -x<min>,<max>  set plot ranges'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
c*****************************************************************

