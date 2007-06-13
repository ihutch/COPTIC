      subroutine alternate(is,ibeg,iodd,ieven,lodd,leven,id)
c Construct a new index array starting at ibeg in stack is,
c that consists of alternating odd and even indexes of the previous
c types which start at iodd, ieven (in is), 
c ibeg is altered.

c storage index space
      integer is(*)
c current address of beginning of unused storage
      integer ibeg
c addresses and lengths of the odd and even types for dimensions nd.
c The values for i=id-1 must exist on input. Values for id are returned.
      integer iodd(*),ieven(*),lodd(*),leven(*)
c dimension at present issue
      integer id
c block length of present dimension.
      integer idlen

      iodd(id)=ibeg
      do i=1,idlen/2
         call catenate(is,ibeg,iodd(id-1),lodd(id-1))
         call catenate(is,ibeg,ieven(id-1),leven(id-1))
      enddo
      lodd(id)=ibeg-iodd(id)

      ieven(id)=ibeg
      do i=1,idlen/2
         call catenate(is,ibeg,ieven(id-1),leven(id-1))
         call catenate(is,ibeg,iodd(id-1),lodd(id-1))
      enddo
      leven(id)=ibeg-ieven(id)

      end
c*****************************************************************
c Concatenate data from istart with length ilen into is at ibeg.
      subroutine catenate(is,ibeg,istart,ilen)
      integer is(*)
      do j=1,ilen
         is(ibeg)=is(istart-1+j)
         ibeg=ibeg+1
      enddo
      end
c****************************************************************
      subroutine facecreate(nn,nd,ibeg,is,iLs,myside,
     $     iodd,ieven,lodd,leven)
c nn: normal dimension, nd: total dimensions.
      integer nn,nd,id
c ibeg: stack counter, is: stack array, iLs: u-structure,
c myside: block length in dimension n.
      integer ibeg,is(*),iLs(nd+1),myside(nd)
c iodd,ieven: returns the vector representing the face indices of is,
c in the elements iodd(nd), ieven(nd) (subfaces earlier).
      integer iodd(nd),ieven(nd),lodd(nd),leven(nd)

c Create vector that gets elements from face normal to nn.
c Iterate over all the dimensions
        do  nc=nn,nn+nd-1
c The actual dimension: n
           n=mod(nc,nd)+1
c Count of the dimension id starting at one.
           id=nc-nn+1
           if(id.eq.1)then
c First create even and odd 1-d vectors: (1,3,5,..) and (2,4,6,...)
              do i=1,myside(n),2
                 iodd(id)=ibeg
                 is(ibeg)=i*iLs(n)
              enddo
              lodd(id)=ibeg-iodd(1)
              do i=2,myside(n),2
                 ieven(id)=ibeg
                 is(ibeg)=i*iLs(n)
              enddo
              leven(id)=ibeg-ieven(1)
           else
c Then bootstrap the lamination to all active dimensions
              call alternate(is,ibeg,iodd,ieven,lodd,leven,id)
           endif
        enddo
c Now iodd(id),lodd(id),ieven(id),leven(id) 
c are the even and odd vectors face normal to dimension nn.
        end
