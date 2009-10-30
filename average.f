      subroutine average3d(q,qave,ifull,iuds,istepave)

      integer ifull(3),iuds(3)
      real q(ifull(1),ifull(3),ifull(3))
      real qave(ifull(1),ifull(3),ifull(3))

      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               qave(i,j,k)=((istepave-1)*qave(i,j,k)+q(i,j,k))
     $              /int(istepave)
            enddo
         enddo
      enddo

      end
