c**********************************************************************
      integer function ipof(ndims,ifull,ico)
c Return the address of the position given by coordinates ico(ndims)
c where ico are fortran (1-based indices) and so is the address ipof.
c Within an array whose full dimension are ifull(ndims)
      parameter (mdim=10)
      ipof=1
      iL=1
      do n=1,ndims
         ipof=(ico(n)-1)*iL
         iL=iL*ifull(n)
      enddo
      end
