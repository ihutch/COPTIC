C*******************************************************************
      subroutine minmax(array,idim,aamin,aamax)
      integer i,idim
      real array(idim),aamin,aamax,arrc
      aamin=array(1)
      aamax=aamin
      do 1 i=2,idim
         arrc=array(i)
         if(arrc.gt.aamax)aamax=arrc
         if(arrc.lt.aamin)aamin=arrc
    1 continue
      return
      end
C
C*******************************************************************
      subroutine minmax2(array,ldim1,idim1,idim2,aamin,aamax)
c Get the min & max of array with leading dimension ldim1.
      integer i,j,ldim1,idim1,idim2
      real array(ldim1,idim2),aamin,aamax,arrc
      aamin=array(1,1)
      aamax=aamin
      do 2 j=1,idim2
      do 1 i=1,idim1
         arrc=array(i,j)
         if(arrc.gt.aamax)aamax=arrc
         if(arrc.lt.aamin)aamin=arrc
    1 continue
    2 continue
      return
      end
C
C*******************************************************************
        function indmin(vector,ilen)
c Return the index of the minimum value of vector.
        integer indmin,ilen,i
        real vector(ilen)
        real aamin
        aamin=0
        do 1 i=1,ilen
                if(aamin.gt.vector(i))then
                        aamin=vector(i)
                        indmin=i
                endif
1       continue
        end
C*******************************************************************
        function indmax(vector,ilen)
c Return the index of the maximum value of vector.
        integer indmax,ilen,i
        real vector(ilen)
        real aamax
        aamax=0
        do 1 i=1,ilen
                if(aamax.lt.vector(i))then
                        aamax=vector(i)
                        indmax=i
                endif
1       continue
        end
