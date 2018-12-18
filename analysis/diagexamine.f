      program diagexamine

      include 'examdecl.f'
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
! Volumes are stored in ndiagmax+1
! Extra work array for arrowplotting in sliceGweb.
      real vp(na_m,na_m2,3,3)
! 1-d plotting arrays.
!      real z1d(na_m),u1d(na_m),dene1d(na_m),deni1d(na_m)
      real zminmax(2)
      character*20 mname(ndiagmax+1)
      integer nfiles,nmodes
      parameter (nfiles=2000,nmodes=11)
      real xmean(nfiles),xvar(nfiles),xmax(nfiles),xmin(nfiles)
      real xms(nfiles),xmcent(nfiles)
      real time(nfiles),xamp(nfiles),work(nfiles)
      real xcentroids(na_m,nfiles),xmodes(nmodes,nfiles)
      real xmodenums(nmodes)
      data xmodenums/0,1,2,3,4,5,6,7,8,9,10/
      integer ifullphi(3),iudsphi(3)
      complex phimodes(nfiles,na_m,nmodes)
      real rphimodes(nfiles,na_m,nmodes),iphimodes(nfiles,na_m,nmodes)
      real theta(nfiles,nmodes)
      real dt

      parameter (nxns=na_m+nfiles+nmodes+1)
      integer ixnps(4)   ! Local ixnp for 3-D plotting.
      real xns(nxns)
      integer mcell,icontour
      logical lvtk,ldebug,nopause
      character*70 xtitle,ytitle,label
      integer iuphi(ndims),iurs(ndims),isrs(ndims)
      integer iunp,i1d,isingle,i1,iwr,istd,linevec(3)
      data iunp/0/i1d/0/iwr/0/zminmax/0.,0./icontour/0/istd/1/
      data ierr/0/ldebug/.false./nopause/.false./
      data mname(1)/'Cellcount'/mname(2)/'v!d1!d'/mname(3)/'v!d2!d'/
      data mname(4)/'v!d3!d'/mname(5)/'T!d1!d'/mname(6)/'T!d2!d'/
      data mname(7)/'T!d3!d'/mname(ndiagmax)/'Potential'/
      data mname(ndiagmax+1)/'Density'/fluxfilename/' '/
      data linevec/-1,0,-1/    ! Default line along y.
      data ifullphi/nfiles,na_m,nmodes/
      isingle=0
      lvtk=.false.
      xleg=.75
      ipp=0
      iworking=0
      iyp=0
      dt=0.
      k=0

! Loop over a maximum of nfiles files: 
      do ifile=1,nfiles
         call diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle,lvtk
     $        ,mcell,zminmax,icontour,iworking,iyp,dt,ldebug,nopause
     $        ,linevec)
! diagexamargs returns here when it reads a diagnostic file.
         if(iworking.ge.0)then
            if(zminmax(1).gt.zminmax(2))istd=0
!-------------------------------------
! The first data file only, get dt and volumes data.
            if(ifile.le.1) call getrundata(dt,istd,ndiags,ldebug)
!-------------------------------------
            i1=1
            ied=ndiagmax
! Create label
            label=diagfilename(istrstr(diagfilename,'dia')+3
     $           :lentrim(diagfilename))
            read(label,*)istep
            time(ifile)=dt*istep
! Read the file whose name we have found in the arguments.
            call array3read(diagfilename,ifull,iuds,ied,diagsum,ierr)
            if(ierr.eq.1)stop 'Error reading diag file'
            ndiags=ied
!-------------------------------------
! Maybe Get up to two additional quantities \phi and n for separate plotting
            call readextra(iuphi)
            if(iunp.ne.0) call unnormplot(diagsum,i1,ndiags)
! Normalize the diagnostics.
            call donormalize(ndiags,diagsum,mcell,isrs,iurs)
            if(ldebug.and.istd.gt.0.and.ifile.eq.1)write(*,'(a,4f8.4)'
     $           )' rs,debyelen,vd,Ti',rs,debyelen,vd,Ti
! Maybe write vtk files and stop.
            if(lvtk)call diavtkwrite(iurs,isrs,diagsum,ndiags,istat)
!-------------------------------------------------------------
! Process and plot this file in various optional ways:
            if(i1d.ne.0)then
! Calculate average profiles in direction i1d.
               call lineout(i1d,ipp,iwr,xtitle,ytitle,xleg,diagsum)
            elseif(lentrim(phifilename).gt.1)then
! Arrow plotting of velocity:
               fluxfilename=ytitle
               ifix=2+4
               write(*,*)'Arrow plotting of ',phifilename
               call sliceGweb(ifull,iuds,u,na_m,zp, ixnp,xn,ifix
     $              ,fluxfilename(1:lentrim(fluxfilename)+2) ,diagsum(1
     $              ,1,1,2),vp)
            elseif(iyp.eq.1)then
! Hole position analysis.
               call hole2position(diagsum,xcentroids(1,ifile),ndiags)
               call getstats(xcentroids(1,ifile),iuds(2),xmean(ifile)
     $              ,xvar(ifile),xmax(ifile),xmin(ifile))
               xamp(ifile)=xmax(ifile)-xmin(ifile)
               if(ldebug)write(*,*)'xmean.max.min',xmean(ifile)
     $              ,xmax(ifile),xmin(ifile)
            elseif(linevec(1).ne.-1)then
! Save a lineout of the chosen diagnostic along specified dimension.
               k=ndiags
               if(isingle.gt.0)k=isingle
               call linewrite(diagsum(1,1,1,k),linevec,dt)
            elseif(iyp.eq.2)then
! Fourier modes in the y direction
               k=ndiags
               if(isingle.gt.0)k=isingle
               call fftdata(diagsum(1,1,1,k),phimodes,nfiles,ifile
     $              ,nmodes,linevec)
!               write(*,*)'Returned from fftdata',nfiles,ifile,nmodes
            else
!----------------------------------------------------------------
! Default case. Default examination of all diagnostics.
               if(istat.eq.1.and.isingle.eq.0)then
! First plot the density if volumes found successfully.
                  k=ndiagmax+1
                  ifix=2
                  fluxfilename=mname(k)
                  write(*,*)k,mname(k)
     $                 ,fluxfilename(1:lentrim(fluxfilename))
                  call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $                 ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)
     $                 +2) ,dum,dum)
               endif
               call examineall(diagsum,zminmax,mname,i1,isingle,istd
     $              ,ndiags,icontour,label)
            endif
!----------------------------------------------------------------
         else   ! iworking.lt.0
! Finished all the arguments: Either do final processing or quit
            if(iyp.gt.0)goto 101
            call exit(0)      
         endif
      enddo
!---------  End of loop over filenames ----------------------------
 101  ifile=ifile-1             ! Final processing
      if(iyp.eq.1)then
! Fitting of xpos vs y kink amplitudes as fn of time.
! Fourier transform the xcentroids and sum low order 10 modes into xms.
         call getmodes(iuds(2),ifile,xcentroids,nmodes,xmodes,10,xms
     $        ,xmcent)
! Plot the hole position and excursion.         
         if(ldebug)call plotexcursion(ifile,time,xmean,xmax,xmin,dt)
! Show the wiggles and modes
         call webinteract(xn(ixnp(2)+1),time,xcentroids,na_m,iuds(2)
     $        ,ifile,'y','time','x!dh!d      ')
         call webinteract(xmodenums,time,xmodes,nmodes,nmodes,ifile,
     $        'mode','time','Amplitude    ')
! Fit and plot the exponential growth.
         call fitexp(time,xamp,ifile,work,it0,imax,A,tau,ldebug,nopause
     $        ,'p-p amplitude')
! Fit to modesum gives nearly the same answer but might be quieter.
         call fitexp(time,xms,ifile,work,it0,imax,A,tau,ldebug,nopause
     $        ,'Mode amplitude')
! Get the mode centroid, averaged over the fit time.
         call modefortime(time,xmcent,ifile,it0,imax,.not.nopause)
      elseif(iyp.eq.2)then
! Fourier analysis in y-direction and finding perturbations.
         write(*,*)'fftdata:',ifile,iuds(2),nmodes
!         write(*,'(10f8.3)')((1.e4*phimodes(j,iuds(2)/2,i),i=1,5),j=1,5)
!         write(*,'(10f6.1)')(time(i),i=1,ifile)
         label=mname(k)(1:lentrim(mname(k)))
         rphimodes=real(phimodes)
         iphimodes=imag(phimodes)
         if(k.eq.ndiags)label='!Af!@ '
         iclipped=10
         ix1=max((iuds(1))/2+1-iclipped,1)
         ix2=min((iuds(1))/2+iclipped,iuds(1))

! Slicing call
         iudsphi(1)=ifile
         iudsphi(2)=ix2-ix1+1
         iudsphi(3)=nmodes
         call ixnpcreate(iudsphi(1),iudsphi(2),iudsphi(3)
     $        ,time,xn(ixnp(1)+ix1),xmodenums,ixnps,xns)
         call setax3chars('time','x','mode')

         if(.false.)then
         call webinteract3(time,xn(ixnp(1)+ix1),rphimodes(1,ix1,1)
     $        ,nfiles,ifile,na_m,ix2-ix1+1,nmodes,'time','x'
     $        ,label(1:lentrim(label)+1) //'r')
         call webinteract3(time,xn(ixnp(1)+ix1),iphimodes(1,ix1,1)
     $        ,nfiles,ifile,na_m,ix2-ix1+1,nmodes,'time','x'
     $        ,label(1:lentrim(label)+1) //'i')
         endif

         call projectphi(nfiles,na_m,nmodes,ifile,ix1,ix2
     $     ,phimodes,rphimodes,theta)
         call sliceGweb(ifullphi,iudsphi,rphimodes(1,ix1,1),na_m,zp,
     $              ixnps,xns,3+64,'Amplitude' ,dum,dum)   

         call pltinit(0.,time(ifile),-1.6,1.6)
         call axis
         call axlabels('time','theta for m=1-4')
         do m=1,4
            call color(m) 
            call polyline(time,theta(1,m),ifile)
         enddo
         call pltend()

      endif

      end
! End of main program.
!********************************************************************
!#####################################################################
! Operational subroutines:
!*************************************************************
      subroutine diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle
     $     ,lvtk,mcell,zminmax,icontour,iworking,iyp,dt,ldebug,nopause
     $     ,linevec)
! Read command line arguments, until a diag name is found.
! If we reach the end of them, return iworking=-1, otherwise return
! iworking= the argument we are working on.
      integer iunp,isingle,i1d
      integer mcell,icontour
      real zminmax(2)
      integer linevec(3)
      character*70 xtitle,ytitle
      logical lvtk,ldebug,nopause
      include 'examdecl.f'

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

! silence warnings:
      zp(1,1,1)=0.
! Defaults
      diagfilename=' '
      phifilename=' '
      denfilename=' '
      ytitle=' '
      xtitle=' '
      mcell=5
      ipfs=3

!      write(*,*)'diagexamargs',iworking,iargc()
! Deal with arguments
      if(iargc().eq.0) goto 201
      do i=iworking+1,iargc()
         iworking=i
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:13).eq.'--objfilename')
     $           read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-p')then
               if(lentrim(phifilename).gt.1)then
                  read(argument(3:),'(a)',err=201)denfilename
               else
                  read(argument(3:),'(a)',err=201)phifilename
               endif
            endif
            if(argument(1:3).eq.'-ly')then
               read(argument(4:),'(a)',err=201)ytitle
            elseif(argument(1:3).eq.'-lx')then
               read(argument(4:),'(a)',err=201)xtitle
            elseif(argument(1:2).eq.'-l')then
               read(argument(3:),*,err=201,end=201)linevec
            endif
            if(argument(1:2).eq.'-u')iunp=1
            if(argument(1:3).eq.'-dt')then
               read(argument(4:),*,err=201)dt
            elseif(argument(1:2).eq.'-d')then
               read(argument(3:),*,err=201)isingle
            endif
            if(argument(1:2).eq.'-o')iwr=1
            if(argument(1:2).eq.'-f')ipp=1
            if(argument(1:2).eq.'-g')then
               read(argument(3:),*,err=204,end=204)ipfs
 204           call pfset(ipfs)
            endif
            if(argument(1:2).eq.'-r')then 
               call noeye3d(0)
               nopause=.true.
            endif
            if(argument(1:2).eq.'-i')ldebug=.not.ldebug
            if(argument(1:2).eq.'-c')
     $           read(argument(3:),*,err=201)icontour
            if(argument(1:2).eq.'-a')
     $           read(argument(3:),*,err=201)i1d
            if(argument(1:2).eq.'-m')
     $           read(argument(3:),*,err=201)mcell
            if(argument(1:3).eq.'-yp')read(argument(4:),*,err=201,end
     $           =201)iyp
            if(argument(1:2).eq.'-z')then
                read(argument(3:),*,err=201)zminmax
! Turn on trucation if the z-range is set
                call togi3trunc()
             endif
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
            if(argument(1:2).eq.'-w')then
               lvtk=.not.lvtk
               read(argument(3:),'(a)',err=201)fluxfilename
            endif
         else
            read(argument(1:),'(a)',err=201)diagfilename
            if(ldebug)write(*,*)diagfilename(1:lentrim(diagfilename))
            return
         endif         
      enddo
      iworking=-1
      return
      goto 202
!------------------------------------------------------------
! Help text
 201  continue
      write(*,*)'=====Error reading command line argument',argument
 203  continue
 301  format(a,i5)
 302  format(a,4f8.3)
 303  format(a,4i5)
      write(*,301)'Usage: diagexamine [switches] <diagfile>'
      write(*,*)' Plot, Analyse, or Output diagnostics from file'
      write(*,*)' If additional parameter file is set,'
     $     ,' do arrow plot of velocity on it.'
      write(*,301)' -p<name>   set name of additional parameter file.'
     $     //' (Given <=twice: phi, den.)'
      write(*,301)' -ly  set label of parameter. -lx label of xaxis'
      write(*,301)' -d   set single diagnostic to be examined [',isingle
      write(*,301)' -a   set dimension number for ave profile [',i1d
      write(*,301)' -yp  if=1 get posn of hole as fn of y     [',iyp
      write(*,301)'      if=2 fourier analyse in y direction   '
      write(*,301)' -o   write out the profiles               [',iwr
      write(*,301)' -f   plot potential profile (if -p given) [',ipp
      write(*,301)' -u   plot un-normalized diagnostics       [',iunp
      write(*,301)' -c<i>  set contour plotting               ['
     $     ,icontour
      write(*,301)' -m<f>  set minimum non-zero cell count    [',mcell
      write(*,302)' -z<min><max> set range of values plotted  [',zminmax
      write(*,302)' -dt  set time-step for multiple steps     [',dt
      write(*,303)' -l<x,y,z> set linevec (+ve=choice,0=direc)[',linevec
      write(*,*)' if zmin>zmax, print zmin and zmax of first diag'
      write(*,*)'-g<i>   print ps-graphics files to unit i [3,-3]'
      write(*,*)'-r   run without pausing'
      write(*,*)'-i   toggle debugging information'
      write(*,*)'-w<name> toggle vtk file writing,'
     $     ,' optionally name the vector[',lvtk,' '
     $     ,fluxfilename(1:lentrim(fluxfilename))
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -h -?  Print usage and current switch values.'
      call exit(0)
 202  continue
      if(lentrim(diagfilename).lt.5)goto 203
      end
!********************************************************************
! End of diagexamargs
!********************************************************************
      subroutine minmax3(ifull,iuds,u,umin,umax)
! Find the minimum and maximum of a 3d array. Allocated dimensions ifull,
! used dimensions iuds.
      parameter (ndims=3)
      integer ifull(ndims),iuds(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      umin=1.e30
      umax=-1.e30
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               if(u(i,j,k).gt.umax)umax=u(i,j,k)
               if(u(i,j,k).lt.umin)umin=u(i,j,k)
            enddo
         enddo
      enddo
      end
!*****************************************************************
      subroutine minmaxN(ndims,ifull,iuds,u,umin,umax)
! Find the minimum and maximum of a ndims array. Allocated dimensions ifull,
! used dimensions iuds. Using the general mditerator. 
      implicit none
      integer ndims
      integer ifull(ndims),iuds(ndims)
      real u(*),umin,umax
! Need local storage.
      integer ndimsmax,indexcontract,mditerator
      parameter (ndimsmax=5)
      integer indi(ndimsmax),iview(3,ndimsmax),i
      external mditerator,indexcontract
      umin=1.e30
      umax=-1.e30
      i=mditerator(ndims,iview,indi,4,iuds)
 1    i=1+indexcontract(ndims,ifull,indi)
         if(u(i).gt.umax)umax=u(i)
         if(u(i).lt.umin)umin=u(i)
      if(mditerator(ndims,iview,indi,0,iuds).eq.0)goto 1
      end
!*************************************************************
! This demonstrates the minimalistic iterator.
      subroutine minmaxM(ndims,ifull,iuds,u,umin,umax)
      integer ndims,ifull(ndims),iuds(ndims)
      real u(*),umin,umax
      integer mditer
      external mditer
      integer ii,i
      umin=1.e30
      umax=-1.e30
      ii=mditer(ndims,ifull,iuds,i)
 1       if(u(i).gt.umax)umax=u(i)
         if(u(i).lt.umin)umin=u(i)
      if(mditer(ndims,ifull,iuds,i).eq.0)goto 1
      end
!**************************************************************
      subroutine examineall(diagsum,zminmax,mname,i1,isingle,istd
     $     ,ndiags,icontour,label)
! Default examination of all diagnostics.
      implicit none
      include 'examdecl.f'
! diagmax here must be 7+1 to accommodate potential possibly.
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
! Volumes are stored in ndiagmax+1
      real zminmax(2)
      integer i1,isingle,istd,ndiags,icontour
      character*20 mname(ndiagmax+1)
      character*70 label
      integer lentrim
      external lentrim

! Local variables
      real dum
      integer k,ifix
      do k=i1,ndiags
         zp(1,1,1)=99
         ifix=3
!            write(fluxfilename,'(''diagnorm('',i1,'')'')')k
         fluxfilename=mname(k)(1:lentrim(mname(k)))
     $        //'('//label(1:lentrim(label))//')'
!            write(*,*)k,isingle
         if(k.eq.ndiags)
     $        fluxfilename='!Af!@('//label(1:lentrim(label))//')'
         if(istd.gt.0.and.(k.eq.isingle.or.isingle.eq.0))write(*,*)k,
     $        fluxfilename(1:lentrim(fluxfilename))
         if(isingle.eq.0.or.isingle.eq.k)then
            if(zminmax(1).lt.zminmax(2))then
               write(*,*)'zminmax',zminmax
               call scale3(0.,1.,0.,1.,zminmax(1),zminmax(2))
               call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $              ixnp,xn,ifix+256+icontour*16+512
     $              ,fluxfilename(1:lentrim(fluxfilename)+2) ,dum
     $              ,dum)   
            elseif(zminmax(1).gt.zminmax(2))then
! If limits are crossed on entry. Find minimum and maximum.
               call minmaxM(ndims,ifull,iuds,diagsum(1,1,1,k),
     $              zminmax(1),zminmax(2))
               write(*,*)zminmax,'=zmin/max'
               call exit(0)
            else
               call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $              ixnp,xn,ifix+icontour*16+512
     $              ,fluxfilename(1:lentrim(fluxfilename)+2) ,dum
     $              ,dum)
            endif
         endif
      enddo
      end
!***********************************************************************
      subroutine lineout(i1d,ipp,iwr,xtitle,ytitle,xleg,diagsum)
      implicit none
      include 'examdecl.f'
      integer i1d,ipp,iwr
      character*70 xtitle,ytitle
      real xleg
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      

! 1-d plotting arrays.
      real z1d(na_m),u1d(na_m),dene1d(na_m),deni1d(na_m)
      real znorm

      integer i,j,k,ii
      integer lentrim
      external lentrim

         do i=1,iuds(i1d)
            z1d(i)=0.
            u1d(i)=0.
            deni1d(i)=0.
         enddo
         write(*,*)'i1d,iuds',i1d,iuds
         do k=2,iuds(3)-1
            do j=2,iuds(2)-1
               do i=2,iuds(1)-1
                  ii=1
                  if(i1d.eq.2)ii=j
                  if(i1d.eq.3)ii=k
                  z1d(ii)=z1d(ii)+diagsum(i,j,k,1+i1d)
                  u1d(ii)=u1d(ii)+u(i,j,k)
                  deni1d(ii)=deni1d(ii)+q(i,j,k)
               enddo
            enddo
         enddo
         znorm=(iuds(1)-2)*(iuds(2)-2)*(iuds(3)-2)/(iuds(i1d)-2)
!         write(*,*)phifilename
         if(phifilename(1:1).ne.' ')write(*,*)'   z           v       '
     $        ,'    phi         n_i        z_E     4pi grad(phi)'
         if(iwr.ne.0)then
            if(phifilename(1:1).ne.' ')then
               write(*,*)iuds(i1d)-2,5
            else
               write(*,*)iuds(i1d)-2
            endif
         endif
! Fix up first and last value for gradient calculation
         u1d(1)=(2.*u1d(2)-u1d(3))/znorm
         u1d(iuds(i1d))=(2.*u1d(iuds(i1d)-1)-u1d(iuds(i1d)-2))/znorm
         do i=2,iuds(i1d)-1
            z1d(i)=z1d(i)/znorm
            u1d(i)=u1d(i)/znorm
            deni1d(i)=deni1d(i)/znorm
            dene1d(i)=exp(u1d(i))
!            write(*,*)i,u1d(i),dene1d(i),deni1d(i)
            if(iwr.ne.0)then
               if(phifilename(1:1).ne.' ')then
                  write(*,'(6f12.5)')xn(ixnp(i1d)+i),z1d(i),u1d(i)
     $                 ,deni1d(i)
     $                 ,(xn(ixnp(i1d)+i)+xn(ixnp(i1d)+i-1))/2.
     $                 ,(u1d(i)-u1d(i-1))*4.*3.14159
     $                 /(xn(ixnp(i1d)+i)-xn(ixnp(i1d)+i-1))
               else
                  write(*,*)xn(ixnp(i1d)+i),z1d(i)
               endif
            endif
         enddo
! Line-out plot.
         call fitinit(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),z1d(2)
     $        ,max(z1d(iuds(i1d)-1),0.))
         call polyline(xn(ixnp(i1d)+2),z1d(2),ixnp(i1d+1)-ixnp(i1d)-2)
         call axptset(1.,0.)
! xticoff reverses the tics.
         call xaxis(0.,0.)
         call ticset(-.015,.015,-.03,.007,0,0,0,0)
         call yaxis(0.,0.)
         call axlabels(xtitle(1:lentrim(xtitle))
     $        ,ytitle(1:lentrim(ytitle)))
         call axlabels('','v!di!d')
         call ticset(-.015,.015,-.03,.007,-1,0,0,0)
         call axptset(0.,1.)
         call altxaxis(1.,1.)
         call axptset(0.,1.)
         call ticset(0.,0.,0.,0.,0,0,0,0)


         if(phifilename(1:1).ne.' ')then
            call legendline(xleg,.15,0,'v!di!d')
            if(ipp.eq.1)then
               call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),u1d(2)
     $              ,u1d(iuds(i1d)-1),.false.,.false.)
               call axptset(-.2,0.)
               call color(4)
               call dashset(1)
               call altyaxis(1.,1.)
               call polyline(xn(ixnp(i1d)+2),u1d(2),ixnp(i1d+1)
     $              -ixnp(i1d)-2)
               call ticset(.015,.015,-.03,-.007,0,0,0,0)
               call axlabels('','!Af!@')
               call legendline(xleg,.2,0,'!Af!@')
            endif
            call axptset(0.,0.)
            call color(5)
            call dashset(2)
!            write(*,*)i1d,iuds(i1d)-1,u1d(2),dene1d(2),u1d(iuds(i1d)-1)
!     $           ,dene1d(iuds(i1d)-1)
            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),dene1d(2)
     $           ,dene1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),dene1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
!            call ticrev()
            call altyaxis(1.,1.)
            call ticset(.015,.015,-.03,-.007,0,0,0,0)
!            call ticset(0.,0.,0.,0.,0,0,0,0)
            call axlabels('','Density')
            call legendline(xleg,.05,0,'n!de!d')
            call ticset(0.,0.,0.,0.,0,0,0,0)
         endif
         if(denfilename(1:1).ne.' ')then
            call color(6)
            call dashset(6)
!            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),deni1d(2)
!     $           ,deni1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),deni1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
            call legendline(xleg,.1,0,'n!di!d')
         endif
         call pltend()

         end
!*********************************************************************
      subroutine readextra(iuphi)
      implicit none
      include 'examdecl.f'
      integer iuphi(3)
      integer j,ied,ierr
      integer lentrim
      external lentrim
      if(lentrim(phifilename).gt.1)then
! Read in a separate potential as well.
         ied=1
         ierr=1
         call array3read(phifilename,ifull,iuphi,ied,u,ierr)
         do j=1,3
            if(iuphi(j).ne.iuds(j))then 
               write(*,*)'Potential array dimensions',iuphi
     $              ,'  incompatible with diagnostics',iuds
               stop
            endif
         enddo
      endif
      if(lentrim(denfilename).gt.1)then
! Read in a density as well.
         ied=1
         ierr=1
         call array3read(denfilename,ifull,iuphi,ied,q,ierr)
         do j=1,3
            if(iuphi(j).ne.iuds(j))then 
               write(*,*)'Density array dimensions',iuphi
     $              ,'  incompatible with diagnostics',iuds
               stop
            endif
         enddo
      endif
      end
!*****************************************************************
      subroutine unnormplot(diagsum,i1,ndiags)
! Unnormalized diagnostic plotting.
      implicit none
      include 'examdecl.f'
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
      integer i1,ndiags
      integer k,ifix
      real dum

      do k=i1,ndiags
! Suppress help.
         zp(1,1,1)=99
         ifix=2
         write(fluxfilename,'(''diagsum('',i1,'')'',a1)')k,char(0)
         call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $        ixnp,xn,ifix,fluxfilename,dum,dum)
      enddo
      end
!****************************************************************
      subroutine donormalize(ndiags,diagsum,mcell,isrs,iurs)
! Normalize by dividing by density, which is first diagnostic.
! But ndiags is phi, so don't normalize that.
      implicit none
      include 'examdecl.f'
      integer ndiags,mcell
      real diagsum(na_i,na_j,na_k,ndiagmax+1)
      integer iurs(3),isrs(3)

      integer id,i,j,k
      do id=2,ndiags-1
         do k=1,iuds(3)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  if(abs(diagsum(i,j,k,1)).gt.mcell+1.e-5)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)/diagsum(i,j,k
     $                    ,1)
                  else
                     diagsum(i,j,k,id)=0.
                  endif
! Subtract v^2 from the second moment to give temperature.             
                  if(id.gt.ndims+1)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)
     $                    -diagsum(i,j,k,id-3)**2
                  endif
               enddo
            enddo
         enddo
      enddo

! Normalize the density= diagsum(1)/volumes into the volumes 
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               if(diagsum(i,j,k,ndiags+1).ne.0.)then
                  diagsum(i,j,k,ndiagmax+1)
     $                 =diagsum(i,j,k,1)/diagsum(i,j,k,ndiags+1)
               else
                  diagsum(i,j,k,ndiagmax+1)=0.
!                  write(*,*)'Warning zero volumes',i,j,k
               endif
            enddo
         enddo
      enddo
!---------------------------------------------------------------
! Decide the actual ends and beginnings of the relevant data.
! If the middle of a face has zero density, that is a dummy face.
! e.g. because of periodicity.
      isrs(1)=1
      if(diagsum(1,iuds(2)/2,iuds(3)/2,1).eq.0)isrs(1)=2
      iurs(1)=iuds(1)+1-isrs(1)
      if(diagsum(iuds(1),iuds(2)/2,iuds(3)/2,1).eq.0)iurs(1)=iurs(1)-1
      isrs(2)=1
      if(diagsum(iuds(1)/2,1,iuds(3)/2,1).eq.0)isrs(2)=2
      iurs(2)=iuds(2)+1-isrs(2)
      if(diagsum(iuds(1)/2,iuds(2),iuds(3)/2,1).eq.0)iurs(2)=iurs(2)-1
      isrs(3)=1
      if(diagsum(iuds(1)/2,iuds(2)/2,1,1).eq.0)isrs(3)=2
      iurs(3)=iuds(3)+1-isrs(3)
      if(diagsum(iuds(1)/2,iuds(2)/2,iuds(3),1).eq.0)iurs(3)=iurs(3)-1
      end
!*****************************************************************
      subroutine diavtkwrite(iurs,isrs,diagsum,ndiags,istat)
! Write Visit-readable vtk file of density, velocity, and potential.
      implicit none
      include 'examdecl.f'
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
      integer iurs(3),isrs(3),ndiags,istat
      integer lentrim

      integer i,ibinary
      if(fluxfilename(1:1).eq.' ')then
         fluxfilename='Velocity'//char(0)
      else
         write(*,*)'Variable name: '
     $        ,fluxfilename(1:lentrim(fluxfilename))
! Spaces are not allowed in visit data names. Fix:
         do i=1,lentrim(fluxfilename)
            if(fluxfilename(i:i).eq.' ')fluxfilename(i:i)='_'
         enddo
         call termchar(fluxfilename)
      endif
      ibinary=0
      call vtkwritescalar(ifull,iurs
     $     ,diagsum(isrs(1),isrs(2),isrs(3),1)
     $     ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $     ,xn(isrs(3)+iuds(1)+iuds(2))
     $     ,ibinary
     $     ,'N'//diagfilename(1:lentrim(diagfilename))//char(0)
     $     ,'Cell_Particles'//char(0))
      if(istat.eq.1)call vtkwritescalar(ifull,iurs
     $     ,diagsum(isrs(1),isrs(2),isrs(3),ndiagmax+1)
     $     ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $     ,xn(isrs(3)+iuds(1)+iuds(2))
     $     ,ibinary
     $     ,'n'//diagfilename(1:lentrim(diagfilename))//char(0)
     $     ,'Density'//char(0))
!         write(*,*)'Finished density vtkwrite',ifull,iuds,u(1,1,1)
      if(ndiags.ge.4)then
         call vtkwritevector(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),2)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary,3
     $        ,'V'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,fluxfilename)
         write(*,*)'Finished velocity vtkwrite',ifull,iuds,u(1,1,1)
! This is supposed to be the potential. It is the last array ndiags.
         call vtkwritescalar(ifull,iurs
     $        ,diagsum(isrs(1),isrs(2),isrs(3),ndiags)
     $        ,xn(isrs(1)),xn(isrs(2)+iuds(1))
     $        ,xn(isrs(3)+iuds(1)+iuds(2))
     $        ,ibinary
     $        ,'u'//diagfilename(1:lentrim(diagfilename))//char(0)
     $        ,'Potential'//char(0))
!         write(*,*)'Finished potential vtkwrite',ifull,iuds,u(1,1,1)
      endif
      stop
      end
!*****************************************************************
      real function thecentroid(diagsum,idiag,id,idp1,idp2)
! Calculate the centroid, \sum_i diagsum()**2*x(idiag,i) of the idiag
! diagnostic, in direction id, at orthogonal index values idp1, idp2
! (weighted by its square value). 
      include 'examdecl.f'
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
      integer idims(ndims)
      idims(mod(id,ndims)+1)=idp1
      idims(mod(id+1,ndims)+1)=idp2
      thecentroid=0.
      thetotal=0.
      do i=1,iuds(id)
         idims(id)=i
         weight=diagsum(idims(1),idims(2),idims(3),idiag)**2
         thecentroid=thecentroid+xn(ixnp(id)+i)*weight
         thetotal=thetotal+weight
      enddo
!      write(*,*)'thecentroid,thetotal',thecentroid,thetotal,idiag
      thecentroid=thecentroid/thetotal
      end
!*****************************************************************
      subroutine hole2position(diagsum,xcentroid,ndiags)
! Detect the hole position in direction id defined as the centroid of the
! potential in the x-direction, for a 2-D x,y calculation.
! Assume that the z-position is 2 because it is ignorable. 
      include 'examdecl.f'
      real xcentroid(na_m)
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
      do j=1,iuds(2)
         xcentroid(j)=thecentroid(diagsum,ndiags,1,j,2)
      enddo
      end
!*****************************************************************
      subroutine getstats(xcentroid,ncen,xmean,xvar,xmax,xmin)
! Get the mean, variance, max and min of xcentroid(1:ncen)
      real xcentroid(ncen)
      xmean=0.
      xvar=0.
      xmax=-1.e20
      xmin=+1.e20
      do i=1,ncen
         x=xcentroid(i)
         xmean=xmean+x
         xvar=xvar+x**2
         if(x.gt.xmax)xmax=x
         if(x.lt.xmin)xmin=x
      enddo
      xmean=xmean/ncen
      xvar=xvar/ncen-xmean**2
      end
!*****************************************************************
      subroutine fitexp(time,xamp,ifile,work,it0,i,A,tau,ldebug,nopause
     $     ,ylabel)
! Fit an exponential line to the data of the growth phase of a trace.
! The line is xamp=A exp(time/tau), but is fitted to the range defined
! by the first maximum value of xamp after entry it0 at i through the
! two points at (i+2*it0)/3 and (2*i+it0)/3. it0 is adjusted
! if the fit is unsatisfactory, and it0 and i are set on return.
      real time(ifile),xamp(ifile),work(ifile)
      logical ldebug,nopause
      character*(*) ylabel

! Hack to smooth once first
      call boxcarover(ifile,1,work,xamp)         
      idone=0
! By default skip the first few files.
 10   it0=10
      if(it0.gt.ifile/2)it0=ifile/2+1
      if(ldebug)write(*,*)'fitexp',ifile,it0
 1    continue
      xamin=1.e20
      xamax=0.
      imin=0
! Find the first peak
      do i=it0,ifile-1
         if(xamp(i).lt.xamin)then
            xamin=xamp(i)
            imin=i
         endif
!         if(xamp(i).lt.xamax)goto 3  ! simple maximum
! Count as a peak only something that stays down for 2 steps:
         if(xamp(i).lt.xamax.and.xamp(i+1).lt.xamax)goto 3
         xamax=xamp(i)
      enddo
      i=i-1
 3    continue
      imax=i
      it1=(imax+2*it0)/3
      xa1=xamp(it1)
      it2=(2*imax+it0)/3
      xa2=xamp(it2)
      dt=time(it2)-time(it1)
      tau=dt/alog(xa2/xa1)
      A=xa1/exp(time(it1)/tau)
      if(ldebug)write(*,*)'it0,imax,A,tau',it0,imax,A,tau
! If it has given us a sensible fit accept it.
!      write(*,*)xamin,xamax
      if((imax-it0.gt.min(10,ifile/3)).and.(tau.gt.5)
     $     .and.(xamax.gt.1.6*xamin))goto 2
! If we are too far into the trace, quit
      if(it0.gt.ifile/2)goto 4
! Else try again, starting later.
!      it0=it0+max(2,it0/4)
      it0=imin+max(2,it0/4)
      goto 1
 2    continue
      write(*,'(2a)')'Fit: tau,    A    tmin,    tmax   ',ylabel
      write(*,'(f8.3,f7.4,2f8.2)')tau,A,time(it0),time(imax)
      do i=it0,imax
         work(i)=A*exp(time(i)/tau)
      enddo
      call lautoplot(time,xamp,ifile,.false.,.true.)
      call axlabels('time',ylabel)
      call dashset(1)
      call color(4)
      call polyline(time(it0),work(it0),imax-it0+1)
      call color(15)
      call dashset(0)
      if(nopause)then
         call prtend(' ')
      else
         call pltend()
      endif
      return
 4    continue
! Fit failure.
      if(idone.lt.5)then
! First attempt failed to fit. Try smoothing and do again.
         idone=idone+1
         nb=idone**2
         write(*,*)'Trying smoothing',nb
         call boxcarover(ifile,nb,work,xamp)         
         goto 10
      endif
      write(*,*)'fitexp failed to fit',it0,imax,ifile
      i=0
      end
!****************************************************************
      subroutine plotexcursion(ifile,time,xmean,xmax,xmin,dt)
      real xmean(ifile),xmax(ifile),xmin(ifile),time(ifile),dt
! Plot the hole position and excursion as a function of time.
         call pltinit(0.,time(ifile),-10.,10.)
         call axis()
         if(dt.eq.1)then
            call axlabels('File','x-positions')
         else
            call axlabels('time [!AXw!@!dp!d]','x-positions')
         endif
         call polyline(time,xmean,ifile)
         call color(1)
         call polyline(time,xmax,ifile)
         call color(2)
         call polyline(time,xmin,ifile)
         call pltend()
         end
!****************************************************************
      subroutine findargvalue(carg,varg,istat,ldebug)
! If there's a copticgeom.dat file, then look for the value of the 
! argument whose switch is carg, and return it in varg. istat=1
! denotes success, istat=0 failure
      character*(*) carg
      real varg
      logical ldebug
      
      character*80 line
      parameter (nlinemax=200)
      real vread

      istat=0
      open(14,file='copticgeom.dat',status='old',err=101)
      do i=1,nlinemax
         read(14,'(a)',end=102,err=102)line
!         write(*,*)line
         imatch=istrstr(line,carg)
         if(imatch.ne.0)then
            read(line(imatch+lentrim(carg):),*,err=102,end=102)vread
            varg=vread
            istat=1
            if(ldebug)write(*,*)'Found value of ',carg,varg
            goto 103
         endif
      enddo
 101  continue
      if(ldebug)write(*,*)'No copticgeom.dat file found'
      return
 102  if(ldebug)write(*,*)'Error reading argvalue in line',line
 103  close(14)
      end
!******************************************************************
      subroutine realtocomplexfft(n,rarray,cfft)
! Calculate the complex fft cfft of a real array rarray of length N.
! Library fftpack5.1 is used. Maximum length of array nmax parameter.
      integer n
      real rarray(n)
      complex cfft(n)

      integer nmax
      INTEGER INC, LENC, LENSAV, LENWRK, IER
      parameter (nmax=2000,inc=1,lensav=(2*nmax+30),lenwrk=2*nmax)
      REAL       WSAVE(LENSAV), WORK(LENWRK)

      if(n.gt.nmax)then
         write(*,*)'realtocomplexfft array length too great',n
         return
      endif
      lenc=n
      do j=1,n
         cfft(j)=rarray(j)
      enddo
      call CFFT1I (n, WSAVE, LENSAV, IER)
      call CFFT1F (n,INC,cfft,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
c B is a straight sum, F is sum divided by N.
c      call CFFT1B (n,INC,cfft,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
      if(ier.ne.0)then
         write(*,*)'CFFT1F error',ier
      endif
c$$$ fftpack5.1 Documentation:
c$$$ Output Arguments
c$$$ 
c$$$ C       For index J*INC+1 where J=0,...,N-1 (that is, for the Jth 
c$$$         element of the sequence),
c$$$ 
c$$$            C(J*INC+1) = 
c$$$            N-1
c$$$            SUM C(K*INC+1)*EXP(-I*J*K*2*PI/N)
c$$$            K=0
c Actually by experiment, forward transform is divided by N, backward not.
c$$$         where I=SQRT(-1).
c$$$ IER     =  0 successful exit
c$$$         =  1 input parameter LENC   not big enough
c$$$         =  2 input parameter LENSAV not big enough
c$$$         =  3 input parameter LENWRK not big enough
c$$$         = 20 input error returned by lower level routine
      end
!*******************************************************************
      subroutine getmodes(n,ifile,xcentroids,nmodes,xmodes,nml,xms
     $     ,xmcent)
! Fourier transform xcentroids and return nmode modes in xmodes
! Also sum the amplitudes of the lower order nml modes into xms.
! And the centroid of the lower modes into xmcent
! Getting the centroid using symmetric modes on either side of
! the peak mode.
      integer n,ifile,nmodes,nml
      real xcentroids(n,ifile),xmodes(nmodes,ifile),xms(ifile)
      real xmcent(ifile)
      include 'examdecl.f'
      complex cfft(na_m)
      
      do i=1,ifile
         call realtocomplexfft(n,xcentroids(1,i),cfft)
         xms(i)=0.
         xmcent(i)=0.
         do j=1,nmodes
            xmodes(j,i)=abs(cfft(j))
            if(j.le.nml+1)then 
! Fit xms to low order modes omitting uniform mode 1.
               if(j.gt.1)xms(i)=xms(i)+xmodes(j,i)
! Centroid of low order modes shifted to origin 0.
               xmcent(i)=xmcent(i)+xmodes(j,i)*(j-1)
            endif
         enddo
         xmcent(i)=xmcent(i)/(xms(i)+xmodes(1,i))
! Alternative limited to near the maximum
         if(.true.)then
! Find the maximum mode (mode of modes!)
            xcmax=0.
            kmax=0
            do k=1,nml
               if(xmodes(k,i).gt.xcmax)then
                  xcmax=xmodes(k,i)
                  kmax=k
               endif
            enddo
            if(kmax.le.2)then
! If the maximum is m=1, take that, not the xcentroid.
               xmcent(i)=kmax-1
               xms(i)=xmodes(kmax,i)
            else
! Use the centroid of the three modes centered on kmax.
               xmcent(i)=0
               xms(i)=0.
               do j=kmax-1,kmax+1
                  xmcent(i)=xmcent(i)+xmodes(j,i)*(j-1)
                  xms(i)=xms(i)+xmodes(j,i)
               enddo
               xmcent(i)=xmcent(i)/xms(i)
            endif
         endif
      enddo
      end
!*********************************************************************
      subroutine webinteract(y,time,xcentroids,na_m,ny,ifile,
     $     lx,ly,lz)
      character*(*) lx,ly,lz
      integer na_m,ifile
      real y(ny),time(ifile)
      real xcentroids(na_m,ifile)
 51   call accisinit()
      isw=1
      call hidweb(y,time,xcentroids,na_m,ny,ifile,isw)
      call ax3labels(lx,ly,lz)
      isw=irotatezoom()
      if(isw.ne.ichar('q').and.isw.ne.0)goto 51
      end
!*********************************************************************
      subroutine webinteract3(x,y,data,nxmax,nx,nymax,ny,nz,
     $     lx,ly,lz)
! Plot interactively data in the x,y plane for different z-levels
! controlled by u and d responses.
      character*(*) lx,ly,lz
      integer nxmax
      real x(nx),y(ny)
      real data(nxmax,nymax,nz)
      character*20 string
      im=1
      string=lz
      lenstr=lentrim(string)
      write(*,*)'Webinteract3. u: raise-z, d: lower-z, q: quit'
 51   call accisinit()
      isw=1
      call hidweb(x,y,data(1,1,im),nxmax,nx,ny,isw)
      call ax3labels(lx,ly,lz)
      call iwrite(im-1,iwidth,string(lenstr+1:))
      call jdrwstr(.1,.1,string,1.)
      isw=irotatezoom()
      if(isw.eq.ichar('u'))im=min(im+1,nz)
      if(isw.eq.ichar('d'))im=max(im-1,1)
      if(isw.ne.ichar('q').and.isw.ne.0)goto 51
      end
!***********************************************************************
      subroutine modefortime(time,xmcent,ifile,it0,imax,lplot)
! Find the dominant mode averaged over ilow to ihigh.
      real time(ifile),xmcent(ifile)
      logical lplot
      icount=0
      xmave=0.
      ilow=(2*it0+imax)/3
      ihigh=(it0+2*imax)/3
      do i=ilow,ihigh
         icount=icount+1
         xmave=xmave+xmcent(i)
      enddo
      xmave=xmave/icount
      write(*,'(a,f6.2,a,2f6.1)')'Mode centroid',xmave
     $     ,' averaged over time',time(ilow),time(ihigh)
      if(lplot)then
         call autoinit(time,xmcent,ifile)
         call axis()
         call polyline(time(it0),xmcent(it0),imax-it0)
         call pltend()
      endif
      end
!********************************************************************
      subroutine boxcarover(ifile,nb,work,xamp)
      real work(ifile),xamp(ifile)
! This leaves xamp smoothed, not in original form.
      do i=1,ifile
         work(i)=xamp(i)
      enddo
      call boxcarave(ifile,nb,work,xamp)
      end
!********************************************************************
      subroutine linewrite(diagsumN,linevec,dt)
      include 'examdecl.f'
      integer linevec(ndims)
      real diagsumN(na_i,na_j,na_k),dt
      integer index(3),idir,istepno
      character*15 lwstring,stepnostring

      lwstring=diagfilename(istrstr(diagfilename
     $     ,'dia'):lentrim(diagfilename))//'.dat'
!      write(*,*)'Diagfilename=',diagfilename,lwstring
      
      open(15,file=lwstring,status='unknown',err=101)
      stepnostring=lwstring(4:istrstr(lwstring,'.dat')-1)
      read(stepnostring,*)istepno
      write(15,*)'istepno,t=',istepno,dt*istepno
      idir=0
      do i=1,3
         index(i)=linevec(i)
         if(linevec(i).eq.0)idir=i
      enddo
      if(idir.eq.0)stop 'linewrite error no zero linevec component'
      
      write(15,*)'Lineout at',linevec,' total values:'
      write(15,*)iuds(idir)
      do i=1,iuds(idir)
         index(idir)=i
        write(15,*)xn(ixnp(idir)+i),diagsumN(index(1),index(2),index(3))
      enddo
      close(15)
      write(*,'(a,3i4,2a)')'Wrote lineout',linevec,' to:  ',lwstring
      return
 101  write(*,*)'Error opening file ',lwstring,'for output'
      end
!*********************************************************************
      subroutine getrundata(dt,istd,ndiags,ldebug)
      implicit none
      include 'examdecl.f'
      real diagsum(na_i,na_j,na_k,ndiagmax+1)      
      logical ldebug
      real dt
      integer istat,istd,ndiags
! Attempt to read dt from the copticgeom.dat file.
            call findargvalue('-dt',dt,istat,ldebug)
            if(dt.eq.0)then
               write(*,'(a)')'Timestep dt value undetermined; set to 1'
               dt=1.
            endif
!-------------------------------------
! Attempt to read the volumes data.
            istat=1
            call stored3geometry(diagsum(1,1,1,ndiags+1),iuds,ifull
     $           ,istat,.false.)
            if(istd.gt.0)then
! Only alert about volumes if we are doing more than getting min/max.
               if(istat.eq.1)then
                  if(ldebug)write(*,'(a,f8.4,$)'
     $                 )'Read volumes ',diagsum(2,2,2,ndiags+1)
               else
                  if(ldebug)then
                     write(*,*)'Storedgeom failed; returned',istat,iuds
     $                    ,ifull
                     write(*,*)
     $               '***** No volume-corrected density is available.'
                  endif
               endif
            endif
      end
!*********************************************************************
      subroutine fftdata(diagsumN,phimodes,nfiles,ifile,nmodes,linevec)
! For this ifile iterate over transverse positions to find fourier modes
! of the diagnostic passed.
! On entry: diagsumN is the data, ifile file index, linevec indicates
!           the line along which to fft, and zero for the variable coordinate.
! On exit: phimodes contains the fourier modes.
      include 'examdecl.f'
      integer nfiles,nmodes
      real diagsumN(na_i,na_j,na_k)
! na_m was too big in some cases, so using na_m. 
! The fix for that is to use -mcmodel=medium compile flag.
      complex phimodes(nfiles,na_m,nmodes),phit(na_m)
      integer ifile,linevec(ndims)
      integer index(ndims)
      real phireal(na_m)

! Set up index.
      idir=0
      do i=1,ndims
         index(i)=linevec(i)
         if(linevec(i).eq.0)idir=i ! The direction in which to FFT
      enddo
      if(idir.eq.0)stop 'fftdata error no zero linevec component'
! Do the transforming
      idx=mod(idir+1,ndims)+1       ! Transverse directions
      idz=mod(idir,ndims)+1
      index(idz)=2               ! Assume 2-d data z-variation absent.
      do k=1,iuds(idx)           ! over transverse
         index(idx)=k
         do i=1,iuds(idir)      ! over line
            index(idir)=i
            phireal(i)=diagsumN(index(1),index(2),index(3))
!            write(*,*)'index',index
         enddo
         call realtocomplexfft(iuds(idir),phireal,phit)
         do im=1,nmodes
            phimodes(ifile,k,im)=phit(im)
!            write(*,*)'Mode',im,iuds(idir)
!            write(*,'(i4,2f8.4,$)')k,phit(im)
         enddo
      enddo
! Now phimodes contains the complex amplitude of modes as a function
! of time (file) and the transverse position. 
      end
!***********************************************************************
      subroutine projectphi(nfiles,na_m,nmodes,ifile,ix1,ix2
     $     ,phimodes,rphimodes,theta)
! How to combine real and imaginary parts: choose a phase angle theta(t) such
! that for each time we maximize the mode amplitude squared averaged over
! all relevant x. theta must be the same for all x (but not t).
      complex phimodes(nfiles,na_m,nmodes)
      real rphimodes(nfiles,na_m,nmodes)
      real theta(nfiles,nmodes)

!      write(*,*)nfiles,na_m,nmodes,ifile,ix1,ix2
      do i=1,ifile
         do m=1,nmodes
            theta(i,m)=0.
            Srr=0.
            Sri=0.
            Sii=0.
            do ix=ix1,ix2
               pr=real(phimodes(i,ix,m))
               pi=imag(phimodes(i,ix,m))
               Srr=Srr+pr**2
               Sri=Sri+pr*pi
               Sii=Sii+pi**2
            enddo
            theta(i,m)=0.5*atan2(Sri,Srr-Sii)
            ct=cos(theta(i,m))
            st=sin(theta(i,m))
            do ix=ix1,ix2
               rphimodes(i,ix,m)=ct*real(phimodes(i,ix,m))
     $              +st*imag(phimodes(i,ix,m))
            enddo
         enddo
      enddo

      end
