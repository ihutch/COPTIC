c Test of versions of mditermults.
      parameter (ndims=3,n1=4,n2=5,n3=3)
      integer ifull(ndims),iused(ndims)
      real u(n1,n2,n3),s(n1,n2,n3)
      data ifull/n1,n2,n3/iused/n1,n2,n3/

      do k=1,n3
         do j=1,n2
            do i=1,n1
               u(i,j,k)=0.
               s(i,j,k)=j
            enddo
         enddo
      enddo
      


      uscale=1.
      iform=1
      write(*,*)'Initial'
      call udisplay(ndims,u,ifull,iused,iform,uscale)

      v=1.
      w=1.
      ipin=0
      call mditermults(u,ndims,ifull,iused,ipin,v,w)

      write(*,*)'After mditermults'
      call udisplay(ndims,u,ifull,iused,iform,uscale)

      call mditeradd(u,ndims,ifull,iused,ipin,s)

      write(*,*)'After mditeradd'
      call udisplay(ndims,u,ifull,iused,iform,uscale)

      end
