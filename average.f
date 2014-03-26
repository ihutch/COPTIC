      subroutine average3d(q,qave,ifull,iuds,istepave)
c 3d only version of averaging. Obsolete.
      integer ifull(3),iuds(3)
      real q(ifull(1),ifull(2),ifull(3))
      real qave(ifull(1),ifull(2),ifull(3))

      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               qave(i,j,k)=((istepave-1.)*qave(i,j,k)+q(i,j,k))
     $              /float(istepave)
            enddo
         enddo
      enddo

      end
c*********************************************************************
c General dimensional version of running average. See mditerate.f
      subroutine averagegd(q,qave,ifull,iuds,istepave)
      include 'ndimsdecl.f'
      ipin=0
      call mditerave(q,ndims,ifull,iuds,ipin,qave,istepave)
      end
