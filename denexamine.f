      program denexamine

      include 'examdecl.f'

      parameter (ntheta=4,nr=30)
      real thetadist(nr,ntheta),thetaval(ntheta),rval(nr)
      integer ithetacount(nr,ntheta)
      
c silence warnings:
      fluxfilename=' '
c 
      call denexamargs

c      write(*,*)ifull
      call array3read(denfilename,ifull,iuds,q,ierr)
      if(ierr.eq.1)stop

      rs=sqrt(3.)*rs
c      write(*,*)na_i,na_j,na_k,na_m,ifull,iuds

      do j=1,nr
         do i=1,ntheta
            thetadist(j,i)=0.
            ithetacount(j,i)=0
            if(i.eq.1)thetaval(i)=-1.+ (i-0.5)/float(ntheta)
         enddo
         rval(j)=1.+(rs-1.)*(j-0.5)/float(nr)
c         write(*,*)j,rval(j)
      enddo

c      write(*,*)(q(16,16,k),k=1,ifull(3))
      ifix=2
      call sliceGweb(ifull,iuds,q,na_m,zp,
     $        ixnp,xn,ifix,'Density:'//'n')

c plot density versus radius.

      write(*,*)rs,debyelen,vd,Ti
c      write(*,*)rs
c      call pltinit(0.,rs,q(iuds(1)/2,iuds(2)/2,iuds(3)/2),0.)
      call pltinit(0.,rs,0.,2.)
      call axis()
      call axis2()
      call axlabels('radius','density')
      call charsize(.001,.001)
      denmin=0.
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               x=xn(ixnp(1)+i)
               y=xn(ixnp(2)+j)
               z=xn(ixnp(3)+k)
               r=sqrt(x**2+y**2+z**2)
               c=z/r
               itheta=nint(0.5+ntheta*(c+1.000001)/2.00001)
               ir=nint(0.5+nr*(r-.999999)/(rs-1.00001))
               if(ir.le.nr .and. ir.ge.1 .and. q(i,j,k).ne.0.)then
                  ithetacount(ir,itheta)=ithetacount(ir,itheta)+1
                  thetadist(ir,itheta)=thetadist(ir,itheta)+q(i,j,k)
               endif
               call polymark(r,q(i,j,k),1,10)
               if(r.gt.rs .and.
     $              .not.(q(i,j,k).eq.0. .or. q(i,j,k).eq.1.))then
c                  write(*,'(4f12.6,3i3)')x,y,z,q(i,j,k),i,j,k
               endif
               if(q(i,j,k).lt.denmin)denmin=q(i,j,k)
            enddo
         enddo
      enddo
      do j=1,nr
         do i=1,ntheta
            if(ithetacount(j,i).ne.0)then
               thetadist(j,i)=thetadist(j,i)/ithetacount(j,i)
            endif
         enddo
      enddo
      nmax=nr
      do j=nr,1,-1
         do i=1,ntheta
            if(ithetacount(j,i).eq.0)then
               ip=mod(i,ntheta)+1
               im=mod(i-2,ntheta)+1
               if(ithetacount(j,ip).ne.0)then
                  thetadist(j,i)=thetadist(j,ip)
               elseif(ithetacount(j,im).ne.0)then
                  thetadist(j,i)=thetadist(j,im)
               else
                  nmax=j-1
                  write(*,*)'Failed setting',j,i,im,ip
                  write(*,'(20i4)')(ithetacount(j,ii),ii=1,ntheta)
               endif
            endif
         enddo
      enddo

      do i=1,ntheta
         call color(mod(i,16))
         call polyline(rval,thetadist(1,i),nmax)
      enddo
      call pltend()

      write(*,*)'r, density distribution (ir,itheta)'
      write(*,*)nmax, ntheta
      do j=1,nmax
         write(*,'(10f8.4)')rval(j),(thetadist(j,i),i=1,ntheta)
      enddo

      end


c*************************************************************
      subroutine denexamargs()
      include 'examdecl.f'

         ifull(1)=na_i
         ifull(2)=na_j
         ifull(3)=na_k

c silence warnings:
      zp(1,1,1)=0.
c Defaults
      denfilename=' '

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
         if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)',err=201)denfilename
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)denfilename
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
      write(*,301)'Usage: denexamine [switches] <denfile>'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' -f   set name of denfile.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(denfilename).lt.5)goto 203
      end
c*****************************************************************
