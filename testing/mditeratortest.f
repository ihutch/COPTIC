c Testing:
c******************************************************************
      program mditeratortest
      integer ndims
      parameter (ndims=3)
      integer ifull(ndims),iLs(ndims+1)
      integer iused(ndims),iu2(ndims),indi(ndims)
      integer iview(3,ndims)
      parameter (istart=1,iend=2,istride=3)
c      external rbroutine,bdyroutine

      parameter (ifull1=100,ifull2=100,ifull3=100)
      parameter (if3=ifull1*ifull2*ifull3)
      real u(if3)


      ifull(1)=ifull1
      ifull(2)=ifull2
      ifull(3)=ifull3
      iused(1)=20
      iused(2)=10
      iused(3)=4
c      iused(1)=100
c      iused(2)=100
c      iused(3)=100
      iu2(1)=iused(1)-4
      iu2(2)=iused(2)-4
      iu2(3)=iused(3)

      iLs(1)=1
      do i=1,ndims
         iLs(i+1)=iLs(i)*ifull(i)
      enddo

      do j=1,if3
            u(j)=0.
      enddo
      mdims=3


      write(*,*)'******** first call ************'

      id=2
      do k=1,100
c Total Setup
      icomplete=mditerator(ndims,iview,indi,4,iused)
 1    u(1+indexcontract(ndims,ifull,indi))=mod(1+indi(1),10)
      if(mditerator(mdims,iview,indi,0,iused).eq.0)goto 1
      enddo

      iview(iend,id)=0
 101  ii=indexcontract(ndims,ifull,indi)
      write(*,*)indi,ii,id,iused(id)
     $    ,u(1+iLs(id)+ii),u(1+(iused(id)-2)*ils(id)+ii)
            u(1+iLs(id)+ii)=u(1+iLs(id)+ii)
     $           +u(1+(iused(id)-1)*iLs(id)+ii)
            u(1+(iused(id)-2)*ils(id)+ii)
     $           =u(1+(iused(id)-2)*iLs(id)+ii)+u(1+ii)
      write(*,*)indi,ii,id,iused(id)
     $    ,u(1+iLs(id)+ii),u(1+(iused(id)-2)*ils(id)+ii)
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 101

      write(*,*)'Complete'
      iform=1
      uscale=1
      if(iused(3).le.6)call udisplay(mdims,u,ifull,iused,iform ,uscale)




      end
c************************************************************************
c
c Test of 100x 100x100x100 takes 1.1s


