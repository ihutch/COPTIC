c************************************************************************
c Print neatly a multidimensional array, u, of ndims up to 3-d.
c Full lengths ifull(ndims), used lengths iuds(ndims), 
c iform digits printed as integer. u multiplied by uscale.
      subroutine udisplay(ndims,u,ifull,iuds,iform,uscale)
      integer ndims
      real u(*)
      integer ifull(ndims),iuds(ndims)

      integer ildim(3),ilil(3)
      character*20 sform      

      if(uscale.eq.0)uscale=1
      iLs=1
      do id=1,3
         if(id.le.ndims)then
            ildim(id)=iuds(id)
            ilil(id)=iLs
         else
            ildim(id)=1
            ilil(id)=0
         endif
         if(id.lt.ndims)iLs=iLs*ifull(id)
      enddo
         
c      write(*,*)'u='
      write(sform,'(''(i4,1x,'',i3,''i'',i2'',i4)'')')
     $     ildim(1),iform
      do k=1,ildim(3)
         il3=(k-1)*ilil(3)
         write(*,*)'k=',k,', columns:i, rows:j'
         write(*,sform)1000000,(mod(i,10**iform),i=1,ildim(1))
         write(*,sform)(j,(int(uscale*u(i+(j-1)*ilil(2)+il3)),
     $        i=1,ildim(1)),j,
     $        j=1,ildim(2))
         write(*,sform)1000000,(mod(i,10**iform),i=1,ildim(1))
      enddo
      end
