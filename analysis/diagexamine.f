      program diagexamine

      include 'examdecl.f'

      parameter (ndiagmax=7)
      real diagsum(na_i,na_j,na_k,ndiagmax)      

c Extra work array for arrowplotting in sliceGweb.
      real vp(na_m,na_m,3,3)
c 1-d plotting arrays.
      real z1d(na_m),u1d(na_m),dene1d(na_m),deni1d(na_m)
      character*20 mname(7)

      character*70 xtitle,ytitle
      integer iuphi(3)
      integer iunp,i1d,isingle,i1,iwr
      data iunp/0/i1d/0/iwr/0/

      mname(1)='Density'
      mname(2)='v!d1!d'
      mname(3)='v!d2!d'
      mname(4)='v!d3!d'
      mname(5)='T!d1!d'
      mname(6)='T!d2!d'
      mname(7)='T!d3!d'
      xleg=.75
      ipp=0
c      pscale=3.
c 
      call diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle)
c      write(*,*)'ifull',ifull
      i1=1
      ied=ndiagmax
      call pfset(3)
      call array3read(diagfilename,ifull,iuds,ied,diagsum,ierr)
      if(ierr.eq.1)stop 'Error reading diag file'
      ndiags=ied
      if(isingle.ne.0)then
         isingle=min(isingle,ndiags)
         ndiags=isingle
         i1=isingle
      endif

c-------------------------------------
      if(lentrim(phifilename).gt.1)then
c Read in a potential as well.
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
c Read in a density as well.
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
c-------------------------------------

      if(iunp.ne.0)then
c Unnormalized diagnostic plotting.
         do k=i1,ndiags
c Suppress help.
            zp(1,1,1)=99
            ifix=2
            write(fluxfilename,'(''diagsum('',i1,'')'',a1)')k,char(0)
            call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $           ixnp,xn,ifix,fluxfilename,dum,dum)
         enddo
      endif

c Normalize by dividing by density, which is first diagnostic.
      do id=2,ndiags
         do k=1,iuds(3)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  if(diagsum(i,j,k,1).gt.5.)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)/diagsum(i,j,k
     $                    ,1)
                  else
                     diagsum(i,j,k,id)=0.
                  endif
c Subtract v^2 from the second moment to give temperature.             
                  if(id.gt.ndims+1)then
                     diagsum(i,j,k,id)=diagsum(i,j,k,id)
     $                    -diagsum(i,j,k,id-3)**2
                  endif
               enddo
            enddo
         enddo
      enddo
      write(*,*)'rs,debyelen,vd,Ti',rs,debyelen,vd,Ti

c      write(*,*)'Normalized diagnostics.'
      if(i1d.ne.0)then
c Calculate average profiles in direction i1d.
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
c         write(*,*)phifilename
         if(phifilename(1:1).ne.' ')write(*,*)'   z           v       '
     $        ,'    phi         n_i        z_E     4pi grad(phi)'
         if(iwr.ne.0)write(*,*)iuds(i1d)-2
c Fix up first and last value for gradient calculation
         u1d(1)=(2.*u1d(2)-u1d(3))/znorm
         u1d(iuds(i1d))=(2.*u1d(iuds(i1d)-1)-u1d(iuds(i1d)-2))/znorm
         do i=2,iuds(i1d)-1
            z1d(i)=z1d(i)/znorm
            u1d(i)=u1d(i)/znorm
            deni1d(i)=deni1d(i)/znorm
            dene1d(i)=exp(u1d(i))
c            write(*,*)i,u1d(i),dene1d(i),deni1d(i)
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
c Line-out plot.
         call fitinit(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),z1d(2)
     $        ,max(z1d(iuds(i1d)-1),0.))
         call pfset(3)
         call polyline(xn(ixnp(i1d)+2),z1d(2),ixnp(i1d+1)-ixnp(i1d)-2)
         call axptset(1.,0.)
c xticoff reverses the tics.
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
            write(*,*)i1d,iuds(i1d)-1,u1d(2),dene1d(2),u1d(iuds(i1d)-1)
     $           ,dene1d(iuds(i1d)-1)
            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),dene1d(2)
     $           ,dene1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),dene1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
c            call ticrev()
            call altyaxis(1.,1.)
            call ticset(.015,.015,-.03,-.007,0,0,0,0)
c            call ticset(0.,0.,0.,0.,0,0,0,0)
            call axlabels('','Density')
            call legendline(xleg,.05,0,'n!de!d')
            call ticset(0.,0.,0.,0.,0,0,0,0)
         endif
         if(denfilename(1:1).ne.' ')then
            call color(6)
            call dashset(6)
c            call fitscale(xn(ixnp(i1d)+1),xn(ixnp(i1d+1)),deni1d(2)
c     $           ,deni1d(iuds(i1d)-1),.false.,.false.)
            call polyline(xn(ixnp(i1d)+2),deni1d(2),ixnp(i1d+1)
     $           -ixnp(i1d)-2)
            call legendline(xleg,.1,0,'n!di!d')
         endif
         call pltend()
      elseif(lentrim(phifilename).gt.1)then
c Arrow plotting of velocity:
         fluxfilename=ytitle
         ifix=2+4
         call sliceGweb(ifull,iuds,u,na_m,zp,
     $        ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $        ,diagsum(1,1,1,2),vp)
      else
c Default examination of all diagnostics.
         do k=i1,ndiags
            zp(1,1,1)=99
            ifix=2
c         write(fluxfilename,'(''diagnorm('',i1,'')'')')k
            fluxfilename=mname(k)
c         write(*,*)mname(k),fluxfilename,lentrim(fluxfilename)
            call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $           ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $           ,dum,dum)
         enddo
      endif
c Arrow plotting of velocity:
c      fluxfilename=mname(1)
c      ifix=2+4
c      if(isingle.eq.0)call sliceGweb(ifull,iuds,diagsum(1,1,1,1),na_m,zp
c     $     ,ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
c     $     ,diagsum(1,1,1,2),vp)

      end


c*************************************************************
      subroutine diagexamargs(iunp,isingle,i1d,iwr,ipp,xtitle,ytitle)
      integer iunp,isingle,i1d
      character*70 xtitle,ytitle
      include 'examdecl.f'

         ifull(1)=na_i
         ifull(2)=na_j
         ifull(3)=na_k

c silence warnings:
      zp(1,1,1)=0.
      isingle=0
c Defaults
      diagfilename=' '
      phifilename=' '
      denfilename=' '
      ytitle=' '
      xtitle=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
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
            if(argument(1:3).eq.'-ly')
     $           read(argument(4:),'(a)',err=201)ytitle
            if(argument(1:3).eq.'-lx')
     $           read(argument(4:),'(a)',err=201)xtitle
            if(argument(1:2).eq.'-u')iunp=1
            if(argument(1:2).eq.'-d')
     $           read(argument(3:),*,err=201)isingle
            if(argument(1:2).eq.'-w')iwr=1
            if(argument(1:2).eq.'-f')ipp=1
            if(argument(1:2).eq.'-a')
     $           read(argument(3:),*,err=201)i1d
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)diagfilename
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
      write(*,301)'Usage: diagexamine [switches] <diagfile>'
      write(*,*)' Plot diagnostics from file'
      write(*,*)' If additional parameter file is set,'
     $     ,' do arrow plot of velocity on it.'
      write(*,301)' -p   set name of additional parameter file.'
     $     //' (Given <=twice: phi, den.)'
      write(*,301)' -ly   set label of parameter. -lx label of xaxis'
      write(*,301)' -d   set single diagnostic to be examined.',isingle
      write(*,301)' -a   set dimension number for average profile',i1d
      write(*,301)' -w   write out the profiles',iwr
      write(*,301)' -f   plot potential profile (if -p given)',ipp
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -u   plot un-normalized diagnostics.',iunp
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(diagfilename).lt.5)goto 203
      end
c*****************************************************************
