!***********************************************************************
! Contains dummy routine(s) that are substituted for the mpi versions
! when a serial code is desired.
!***********************************************************************
! Block boundary communication.
      subroutine bbdy(iLs,ifull,iuds,u,kc,
     $     ndims,idims,icoords,iLcoords,myside,myorig,
     $     icommcart,mycartid,myid,lperiod)
! Dimensional structure of u, for 2d should be (1,Li,Lj), 
! 3d (1,Li,Li*Lj,Li*Lj*Lk) etc (last element may not be used)
      integer iLs(ndims+1)
! ifull full dimensions of u
      integer ifull(ndims)
! iuds used dimensions of u
      integer iuds(ndims)
! Inside this routine, u and iorig are referenced linearly.
      real u(*)
! kc      is iteration count which determines odd or even, and also 
!      if kc=-1 this is the final call: gather only.
!      if kc=-2 only (re)initialize the topology, cartesian communicator
      integer kc
! The number of dimensions of the cartesian topology. (2 for 2d) (IN)
      integer ndims
! The length of each topology dimension (number of blocks) (IN)
      integer idims(ndims)
! For each topology dimension whether it is periodic or not (IN)
      logical lperiod(ndims)
! Cartesian topology coords of this process block (OUT)
      integer icoords(ndims)
! structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndims+1)
! The origin of my block (OUT)
      integer myorig 
! My side length data, accounting for my position in the topology (OUT).
      integer myside(ndims)
! icommcart is the id of the cartesian topology communicator (OUT).
      integer icommcart
! mycartid is the process id in cartesian topology communicator (OUT).
      integer mycartid
! myid returns the process id in the MPI_COMM_WORLD (OUT)
      integer myid
! End of arguments.

! Set some arguments for this dummy version.
      myorig=1
      icommcart=0 
      myid=0
      mycartid=0
      iLs(1)=1
      do id=1,ndims
         icoords(id)=0
         myside(id)=iuds(id)
         iLs(id+1)=iLs(id)*ifull(id)
      enddo

      end
!******************************************************************
!********************************************************************
! Abstraction to isolate mpi calls.
      subroutine mpifinalize(ierr)
      ierr=0
      end
!*******************************************************************
      subroutine mpicommsize(nprcsses,ierr)
      nprcsses=1
      ierr=0
      end
!********************************************************************
      subroutine mpigetmyid(myid,nprcsses,ierr)
      myid=0
      nprcsses=1
      ierr=0
      end
!********************************************************************
      subroutine mpiconvgreduce(convgd,icommcart,ierr)
      implicit none
      real convgd(3)
      integer ierr,icommcart
      ierr=0
      end
!********************************************************************
      subroutine mpibarrier(ierr)
      ierr=0
      end
