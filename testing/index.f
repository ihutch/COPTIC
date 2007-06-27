      parameter (ndims=3)
      integer ifull(ndims),ix(ndims)

      ifull(1)=10
      ifull(2)=10
      ifull(3)=10

 1    write(*,*)'Enter an index:'
      read(*,*)index
      call indexexpand(ndims,ifull,index,ix)
      write(*,*)'index=',index,'  expanded=',ix
      if (index.ne.0)goto 1

      end
