c New version for more than three fonts.
      program fontshow2
c Plot aligned sets of the three installed fonts.
      character*256 entered
      character*129 str1
      character*1 shift
      real js,px,py
      integer i,j
c
      call pfset(-1)
      call pltinit(0.,1.,0.,1.)
      do 8 k=1,3
         if(k.eq.3)then
            noff=1
         else
            noff=0
         endif
         do 9 j=1,3
            js=0.
            do 10 i=1,32
               py=0.6-(4*(k-1)+j)*0.043
               px=0.02+.03*i
               if(char(31+32*(j-1)+i).eq.'\\'
     $              .or.char(31+32*(j-1)+i).eq.'!')then 
                  str1='!'//char(63+k)//'!'//
     $                 char(31+32*(j-1)+i)//'!@'//char(0)
               else
                  str1='!'//char(63+k)//char(31+32*(j-1)+i)
     $                 //'!@'//char(0)
               endif
               call drwstr(px,py,str1)
   10       continue
    9    continue
    8 continue
      do 11 i=1,15
         call color(i)
         call iwrite(i,iwidth,str1)
         call termchar(str1)
         px=0.02+0.03*i
         call drwstr(px,0.7,str1)
         call polymark(xn2xw(px+.01),yn2yw(0.65),1,i)
c         call polymark(px,0.65,1,i)
 11   continue
c Display the 15 named colors.
      call color(idarkblue())
      call drwstr(.03,.08,'idarkblue')
      call color(idarkgreen())
      call drcstr(' idarkgreen')
      call color(iskyblue())
      call drcstr(' iskyblue')
      call color(ibrickred())
      call drcstr(' ibrickred')
      call color(ipurple())
      call drcstr(' ipurple')
      call color(igold())
      call drcstr(' igold')
      call color(igray())
      call drcstr(' igray')
      call color(ilightgray())
      call drcstr(' ilightgray')
      call color(iblue())
      call drwstr(.03,.04,' iblue')
      call color(igreen())
      call drcstr(' igreen')
      call color(icyan())
      call drcstr(' icyan')
      call color(ired())
      call drcstr(' ired')
      call color(imagenta())
      call drcstr(' imagenta')
      call color(iyellow())
      call drcstr(' iyellow')
      call color(iblack())
      call drcstr(' iblack')
      call pltend()
c Second plot.
      call pltinit(0.,1.,0.,1.)
      do k=4,6
         if (k.eq.4)then 
            shift='R'
         elseif(k.eq.5)then
            shift='E'
         elseif(k.eq.6)then
            shift='G'
         endif
         if(k.eq.3)then
            noff=1
         else
            noff=0
         endif
         do j=1,3
            js=0.
            do i=1,32
               py=0.6-(4*(k-4)+j)*0.043
               px=0.02+.03*i
               if(char(31+32*(j-1)+i).eq.'\\'
     $              .or.char(31+32*(j-1)+i).eq.'!')then 
                  str1='!'//shift//'!'//
     $                 char(31+32*(j-1)+i)//'!@'//char(0)
               else
                  str1='!'//shift//char(31+32*(j-1)+i)
     $                 //'!@'//char(0)
               endif
               call drwstr(px,py,str1)
            enddo
         enddo
      enddo
      call pltend()
      stop
      end

