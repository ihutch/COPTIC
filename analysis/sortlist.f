      subroutine sortbottomlimit(v,vlist,nvlist)
      real vlist(nvlist)
! Insert the value v into its ordered place in vlist, retaining the bottom
! nvlist values.
!      write(*,*)'bot',v,nvlist,vlist
      if(v.ge.vlist(nvlist))return
      if(v.ge.vlist(1))then
         i1=1
         i2=nvlist
! Find my position by bisection.
 1       i=(i1+i2)/2
         if(v.lt.vlist(i))then
            i2=i
         else
            i1=i
         endif
         if(i2.gt.i1+1)goto 1
      else
         i2=1
      endif
! Here i2 is the position of list value just greater than v.
!      write(*,*)'botend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v

      end
!*********************************************************************
      subroutine sorttoplimit(v,vlist,nvlist)
      real vlist(nvlist)
! Insert the value v into its reverse-ordered place in vlist, retaining
! the top nvlist values.
!      write(*,*)'topstart',nvlist,v,vlist
      if(v.le.vlist(nvlist))return
      if(v.le.vlist(1))then
         i1=1
         i2=nvlist
! Find my position by bisection.
 1       i=(i1+i2)/2
         if(v.gt.vlist(i))then
            i2=i
         else
            i1=i
         endif
         if(i2.gt.i1+1)goto 1
      else
         i2=1
      endif
! Here i2 is the position of list value just less than v.
!      write(*,*)'topend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v

      end

!*************************************************************
      program testsorts
      real v
      parameter (nvlist=5,nnum=20)
      real vbotlist(nvlist)
      real vtoplist(nvlist)
      data vbotlist/nvlist*1.e20/vtoplist/nvlist*(-1.e20)/

      v=ran1(-1)
      do i=1,nnum
         v=ran1(1)
!         write(*,*)'Random',v
         call sortbottomlimit(v,vbotlist,nvlist)
!         write(*,*)i,v,nvlist,vtoplist
         call sorttoplimit(v,vtoplist,nvlist)
         write(*,101)vbotlist,vtoplist
      enddo
 101  format(10f7.4)
      end
