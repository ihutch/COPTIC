      program diagexamine

      include 'examdecl.f'

      parameter (ndiagmax=7)
      real diagsum(na_i,na_j,na_k,ndiagmax)      

c Extra work array for arrowplotting in sliceGweb.
      real vp(na_m,na_m,3,3)
      character*20 mname(7)
      integer iuphi(3)
      integer iunp
      data iunp/0/

      mname(1)='Density'
      mname(2)='v!d1!d'
      mname(3)='v!d2!d'
      mname(4)='v!d3!d'
      mname(5)='T!d1!d'
      mname(6)='T!d2!d'
      mname(7)='T!d3!d'

c 
      call diagexamargs(iunp)

c      write(*,*)ifull
      ied=ndiagmax
      call array3read(diagfilename,ifull,iuds,ied,diagsum,ierr)
      if(ierr.eq.1)stop 'Error reading diag file'
      ndiags=ied

      if(lentrim(phifilename).gt.1)then
c Read in a potential as well.
         ied=1
         call array3read(phifilename,ifull,iuphi,ied,u,ierr)
         do j=1,3
            if(iuphi(j).ne.iuds(j))then 
               write(*,*)'Potential array dimensions',iuphi
     $              ,'  incompatible with diagnostics',iuds
               stop
            endif
         enddo
      endif

      if(iunp.ne.0)then
         do k=1,ndiags
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
c      write(*,*)rs,debyelen,vd,Ti

c      write(*,*)'Normalized diagnostics.'
      do k=1,ndiags
         zp(1,1,1)=99
         ifix=2
c         write(fluxfilename,'(''diagnorm('',i1,'')'')')k
         fluxfilename=mname(k)
c         write(*,*)mname(k),fluxfilename,lentrim(fluxfilename)
         call sliceGweb(ifull,iuds,diagsum(1,1,1,k),na_m,zp,
     $        ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $        ,dum,dum)
      enddo
c Arrow plotting of velocity:
      fluxfilename=mname(1)
      ifix=2+4
      call sliceGweb(ifull,iuds,diagsum(1,1,1,1),na_m,zp,
     $     ixnp,xn,ifix,fluxfilename(1:lentrim(fluxfilename)+2)
     $     ,diagsum(1,1,1,2),vp)

      end


c*************************************************************
      subroutine diagexamargs(iunp)
      include 'examdecl.f'

         ifull(1)=na_i
         ifull(2)=na_j
         ifull(3)=na_k

c silence warnings:
      zp(1,1,1)=0.
c Defaults
      diagfilename=' '
      phifilename=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:13).eq.'--objfilename')
     $           read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-p')
     $           read(argument(3:),'(a)',err=201)phifilename
            if(argument(1:2).eq.'-u')iunp=1
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
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -p   set name of potential file.'
      write(*,301)' -u   plot un-normalized diagnostics.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(diagfilename).lt.5)goto 203
      end
c*****************************************************************
