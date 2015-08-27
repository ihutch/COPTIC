      program phasespace
c Read particle data and construct phase space plots.
      include 'examdecl.f' 
c (Examdecl itself includes meshcom.f plascom.f, objcom.f)
      parameter (nfilemax=999)
      include '../partcom.f'
      include '../ptaccom.f'
 
      character*10 chartemp
      character*100 name
      character*100 string
      logical ldoc
      real extra(nptdiagmax,ndimsmax),diff(nptdiagmax,ndimsmax)
c Spatial limits bottom-top, dimensions
      real xlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      real Bdirs(ndimsmax+1)
c Velocity limits
      real vlimit(2,ndimsmax),wicell,wjcell,wkcell
      
      integer nx,nv
      parameter (nx=50,nv=50)
      real hist(nx,nv),cworka(nx,nv),zclv(nx)

      integer iuin(ndimsmax)
      integer ivtk
      integer ndfirst,ndlast
      character*200 ivarnames(2*ndimsmax)

      nfmax=nfilemax

c Defaults
      ndfirst=1
      ndlast=3
      do id=1,ndimsmax
         xlimit(1,id)=-5.
         xlimit(2,id)=5.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
         vlimit(1,id)=-5.
         vlimit(2,id)=5.
      enddo

      call partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc,ivtk
     $     ,ispecies,ndfirst,ndlast)

      do id=1,ndimsmax
c Needed? initialization removed from partacinit.
         xmeshstart(id)=min(-5.,xlimit(1,id))
         xmeshend(id)=max(5.,xlimit(2,id))
      enddo
c Now the base filename is in partfilename.

      ip=lentrim(partfilename)-3
      if(partfilename(ip:ip).eq.'.')then
c If filename has a 3-character extension or is a numbered pex
c file, assume it is complete and that we are reading just one
c file. Should do this only the first time.
         nfmax=0
         name=partfilename
         write(*,*)'Reading single file ',name(1:lentrim(name))
         if(partfilename(ip:ip+3).eq.'.pex')then
            write(*,*)'Using stored distribution file'
            nfmax=-1
         endif
      elseif(partfilename(ip-4:ip-1).eq.'.pex')then
         nfmax=-1
         name=partfilename
         write(*,*)'Reading numbered pex file ',name(1:lentrim(name))
      endif
c Possible multiple files.
      nfvaccum=0
      nhist=0
      do i=0,nfmax
         if(nfmax.ne.0)then
            write(chartemp,'(''.'',i3.3)')i
            name=partfilename(1:lentrim(partfilename))//chartemp
            write(*,*)'Reading file ',name(1:lentrim(name))
         endif
         Bt=0.
         call partread(name,ierr)
         if(ispecies.gt.nspecies)then
            ispecies=nspecies
            write(*,*)'nspecies=',nspecies,'  Reset ispecies',ispecies
         endif
         if(ldoc)then
            write(*,*)'debyelen,Ti,vd,rs,phip=',debyelen,Ti,vd,rs,phip
            write(*,*)'iregion_part,n_part,dt,ldiags,rhoinf,nrein,',
     $           'phirein,numprocs='
            write(*,*)iregion_part,n_part,dt,ldiags,rhoinf,nrein,phirein
     $           ,numprocs
            write(*,*)'eoverm,Bt,Bfield,vpar,vperp=',eoverm,Bt,Bfield
     $           ,vpar,vperp
            write(*,*)'nspecies',nspecies
c            stop
         endif
         if(ierr-4*(ierr/4).ne.0)goto 11
         if(Bdirs(4).gt.0. .or. Bt.eq.0)then
c All directions were set by commandline. Or none were read from file.
            do k=1,ndims
               Bfield(k)=Bdirs(k)
            enddo
         endif
         if(cellvol.eq.-1)write(*,*)'Bfield (projection)',Bfield
c Now act on the particle data we have read in:
c Choose active axis based on b-settings
         do k=1,ndims
            if(Bfield(k).ne.0)iaxis=k
         enddo
         call plotlocalphasepoints(xlimit,vlimit,iaxis,ispecies
     $        ,nfvaccum)
         call fhistinc(xlimit,vlimit,iaxis,ispecies,nhist,hist,nx
     $        ,nv)
      enddo
 11   continue
      call pltend()
      
c      write(*,*)((i,j,hist(i,j),i=nx/2,nx/2+2),j=nv/2,nv/2+2)
c Color Contour of Histogram.
      call minmax2(hist,nx,nx,nv,hmin,hmax)
      ncont=10
      do i=1,ncont
         zclv(i)=hmin+(i-1)*(hmax-hmin)/(ncont-1)
      enddo
      icl=ncont
      icsw=0+16
      call pltinit(1.,float(nx),1.,float(nv))
      call  blueredgreenwhite
      call contourl(hist,cworka,nx,nx,nx,zclv,icl,dummy,dummy,icsw) 
      call gradlegend(zclv(1),zclv(icl),-.2,.1,-.2,1.,-.05,.false.) 
      call scalewn(xlimit(1,iaxis),xlimit(2,iaxis),vlimit(1,iaxis)
     $     ,vlimit(2,iaxis),.false.,.false.)
      call axis
      call pltend()
      end
c*************************************************************
      subroutine plotlocalphasepoints(xlimit,vlimit,iaxis,ispecies
     $     ,nfvaccum)
c Plot points in phase space in direction iaxis for particles that lie
c in the space and velocity ranges given by xlimit, vlimit.
      include '../ndimsdecl.f'
      include '../partcom.f'
      real xlimit(2,ndims),vlimit(2,ndims)
      character*1 axnames(3)
      data axnames/'x','y','z'/ 

c Initialize the plot
      if(nfvaccum.eq.0)then
c         call pfset(3)
         call pltinit(xlimit(1,iaxis),xlimit(2,iaxis),vlimit(1,iaxis)
     $     ,vlimit(2,iaxis))
         call axis()
         call axlabels(axnames(iaxis),'v!d'//axnames(iaxis)//'!d')
         call color(ibrickred())
      endif
c      write(*,*)iicparta(ispecies),iocparta(ispecies)
      do j=iicparta(ispecies),iocparta(ispecies)
c Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
               x=x_part(id,j)
               v=x_part(id+ndims,j)
c               write(*,*)'x',id,x,xlimit(1,id),xlimit(2,id)
c               write(*,*)'v',id,v,vlimit(1,id),vlimit(2,id)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(v.lt.vlimit(1,id).or.v.gt.vlimit(2,id))goto 12
c               if(limadj.eq.1)then
c Adjust the xnewlim.
c                  if(x.lt.xnewlim(1,id))xnewlim(1,id)=x
c                  if(x.gt.xnewlim(2,id))xnewlim(2,id)=x
c               endif
            enddo
            nfvaccum=nfvaccum+1
c Here for a particle lying within the limits.
c Plot a point in phase space
c            write(*,'(2f8.4)')x_part(iaxis,j),x_part(iaxis+ndims,j)

            call vecw(x_part(iaxis,j),x_part(iaxis+ndims,j),-1)

         endif
 12      continue
      enddo
      write(*,*)'Accumulated',nfvaccum,' particles'
      end
c*************************************************************
c*************************************************************
      subroutine fhistinc(xlimit,vlimit,iaxis,ispecies
     $     ,nhist,hist,nx,nv)
c Increment histogram of phase space in direction iaxis for particles
c that lie in the space and velocity ranges given by xlimit, vlimit.
      include '../ndimsdecl.f'
      include '../partcom.f'
      real xlimit(2,ndims),vlimit(2,ndims)
      real hist(nx,nv)
      character*1 axnames(3)
      data axnames/'x','y','z'/ 
      save x1,x2,xdiff,v1,v3,vdiff

      if(nhist.eq.0)then
c Initialize
         x1=xlimit(1,iaxis)
         x2=xlimit(2,iaxis)
         xdiff=x2*1.00001-x1
         v1=vlimit(1,iaxis)
         v2=vlimit(2,iaxis)
         vdiff=v2*1.00001-v1
         write(*,*)'x1,x2,xdiff,v1,v2,vdiff',x1,x2,xdiff,v1,v2,vdiff
         do j=1,nv
            do i=1,nx
               hist(i,j)=0.
            enddo
         enddo
      endif
c      write(*,*)iicparta(ispecies),iocparta(ispecies)
      do j=iicparta(ispecies),iocparta(ispecies)
c Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
               x=x_part(id,j)
               v=x_part(id+ndims,j)
c               write(*,*)'x',id,x,xlimit(1,id),xlimit(2,id)
c               write(*,*)'v',id,v,vlimit(1,id),vlimit(2,id)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(v.lt.vlimit(1,id).or.v.gt.vlimit(2,id))goto 12
               ix=nx*(x-x1)/xdiff+1
               iv=nv*(v-v1)/vdiff+1
               if(id.eq.iaxis)then
c Here for a particle lying within the limits.
                  nhist=nhist+1
                  hist(ix,iv)=hist(ix,iv)+1.
               endif
            enddo
         endif
 12      continue
      enddo
      write(*,*)'Accumulated',nhist,' particles'
      end
c*************************************************************
c*************************************************************
      subroutine partexamargs(xlimit,vlimit
     $           ,iuin,cellvol,Bdirs,ldoc,ivtk,ispecies,ndfirst,ndlast)
      include 'examdecl.f'
      include '../ptaccom.f'
      real xlimit(2,3),vlimit(2,3),Bdirs(4)
      integer iuin(3)
      logical ldoc
      integer ivtk
c I think unused here 26 May 12. But I'm not sure.
      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k
      ivtk=0
      ispecies=1

      do i=1,3
         iuin(i)=9
c Default field/projection direction:
         Bdirs(i)=0
         if(i.eq.3)Bdirs(i)=1
      enddo
c Use cellvol=0. by default.
      cellvol=0.
      ldoc=.false.

c silence warnings:
      fluxfilename=' '
      zp(1,1,1)=0.
c Defaults
      partfilename=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,err=201) xlimit(1,1),xlimit(2,1)
            elseif(argument(1:2).eq.'-y')then
               read(argument(3:),*,err=201) xlimit(1,2),xlimit(2,2)
            elseif(argument(1:2).eq.'-z')then
               read(argument(3:),*,err=201) xlimit(1,3),xlimit(2,3)
            elseif(argument(1:2).eq.'-u')then
               read(argument(3:),*,err=201) vlimit(1,1),vlimit(2,1)
            elseif(argument(1:4).eq.'-vtk')then
               ivtk=1
            elseif(argument(1:2).eq.'-v')then
               read(argument(3:),*,err=201) vlimit(1,2),vlimit(2,2)
            elseif(argument(1:2).eq.'-w')then
               read(argument(3:),*,err=201) vlimit(1,3),vlimit(2,3)
            elseif(argument(1:3).eq.'-pu')then
               nptdiag=nsbins
            elseif(argument(1:2).eq.'-p')then
               cellvol=-1
               bread=0.
               read(argument(3:),*,err=102)(Bdirs(j),j=1,3)
               bread=1.
 102           continue
               Bdirs(4)=0.
               do j=1,3
                  Bdirs(4)=Bdirs(4)+Bdirs(j)**2
               enddo
               if(Bdirs(4).gt.0)then
                  Bdirs(4)=sqrt(Bdirs(4))
                  do j=1,3
                     Bdirs(j)=Bdirs(j)/Bdirs(4)
                  enddo
               endif
               Bdirs(4)=Bdirs(4)*bread
            elseif(argument(1:2).eq.'-b')then
               read(argument(3:),*,err=201)iuin
            elseif(argument(1:2).eq.'-q')then
               ldoc=.true.
            elseif(argument(1:3).eq.'-sp')then
               ispecies=ispecies+1
            elseif(argument(1:2).eq.'-d')then
               read(argument(3:),*,err=201)ndfirst,ndlast
            endif
            if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-f')
     $           read(argument(3:),'(a)',err=201)partfilename
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)partfilename
c            write(*,*)partfilename
         endif
         
      enddo
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5,a)
 302  format(a,4f8.3)
 303  format(a,5i5)
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,302)' -x -y -z<fff,fff>  set position range. [',
     $     xlimit(1,1),xlimit(2,1)
      write(*,302)' -u -v -w<fff,fff>  set velocity range. ['
     $     ,vlimit(1,1),vlimit(2,1)
      write(*,303)' -b<nx,ny,nz>  set spatial block range. [',iuin
      write(*,302)' -p[bx,by,bz]  project [in direction]   [',Bdirs
      write(*,'(a,$)')' -d<i,i> Set first,last v-index'
      write(*,303)' to plot [',ndfirst,ndlast
      write(*,301)' -sp     Increment species number to examine'
      write(*,301)' -pu     Set nptdiag for uniform bins   [',nptdiag
     $     ,'  (Hence 2-D f plots)'
      write(*,301)' -q      Output diagnostics of file reading'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -vtk    Output distribution function vtk files.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)then
         write(*,*)'Short filename, length<5 not allowed'
         goto 203
      endif
      end
