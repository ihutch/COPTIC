c***********************************************************************
c This file contains the extra routines that are needed for parallel
c processing the particles. (mpibbdy contains others for fields)
c 
c For a non-mpi installation, the subroutines should be dummies, simply
c returning. Thus the particle code can be immediately un-MPIed.
c
c***********************************************************************
      subroutine fluxreduce()
      include 'mpif.h'
      include '3dcom.f'
c ALLREDUCE communicates the result to all processes.      
      call MPI_ALLREDUCE(MPI_IN_PLACE,ff_data(nf_address(1,1,nf_step))
     $     ,nf_posno(1,1),MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      write(*,*)'Total flux number',sum
c      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
c     $     i=1,nf_posno(1,1))

      end
c********************************************************************
      subroutine psumreduce(psum,ndims,ifull,iuds,iLs)
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real psum(*)
      include 'mpif.h'
      include 'partcom.f'
c Operator
      external addsubarray_MPI

      logical lfirst
      integer iaddtype,iaddop,iporig
      data lfirst/.true./
      save lfirst,iaddtype,iaddop,iporig

      if(lfirst)then
c Create addtype and operator for reduce sum.
         call mpisubopcreate(ndims,ifull,iuds,addsubarray_MPI,
     $        iaddtype,iaddop)
         iporig=1
         do i=1,ndims
            iporig=iporig+iLs(i)
         enddo
      endif
c Pass the (2,2,...) element of psum.
      call MPI_ALLREDUCE(MPI_IN_PLACE,psum(iporig),1,iaddtype,
     $     iaddop,MPI_COMM_WORLD,ierr)
c Reduce also the nrein and phirein.
      call MPI_ALLREDUCE(MPI_IN_PLACE,nrein,1,MPI_INTEGER,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,phirein,1,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      phirein=phirein/numprocs

      end

c********************************************************************

c**************************************************************
      subroutine mpisubopcreate(nrd,ifull,iuds,iopfun,itype,iaddop)
c For an nrd (IN) dimensional structure 
c described by ifull (IN), iuds (IN)
      integer nrd
      integer iuds(nrd),ifull(nrd)
c Create a datatype itype (OUT) and an operator iaddop (OUT)         iporig

c for doing MPI_REDUCEs on a subarray described by iuds
c (Its boundaries are omitted in the reduce.)
c The operation function is iopfun (EXTERNAL).
c iopfun MUST be declared external in the calling routine.
c
c Local storage
      parameter (ims=10)
      integer iside(2,ims),ktype(2**(ims+1)),iLs(ims+1)
      
      include 'mpif.h'
c Common for passing the dimensionals structures:
      include 'mditcom.f'
      external iopfun

      if(nrd.gt.ims)stop 'mpisubopcreate dimension too large'
      iLs(1)=1
      do i=1,nrd
c These are needed for the blockcreate
         iside(1,i)=iuds(i)
         iside(2,i)=iuds(i)
         iLs(i+1)=iLs(i)*ifull(i)
c These are needed for the ioperator usage.
         iasfull(i)=ifull(i)
         iasum2(i)=iuds(i)-2
      enddo
      nasdims=nrd
c iasfull etc don't need to be set for this call.
      call MPI_OP_CREATE(iopfun,.false.,iaddop,ierr)
      call bbdyblockcreate(nrd,ktype,iLs,iside,iSIZEOFREAL)
      itype=ktype(2**nrd)

c MPI_IN_PLACE is ok for ALLREDUCE. Not for REDUCE.
c Subsequent calls will be of the form:
c      call MPI_ALLREDUCE(MPI_IN_PLACE,u(2,2,2),1,itype,
c     $     iaddop,icommcart,ierr)
      end

c**************************************************************
c These two routines provide facility for addition reduce 
c over subarray.
      subroutine addsubarray_MPI(invec,inoutvec)
c actually (invec,inoutvec,ilen,itype) but ilen,itype not used.
c MPI_sum type function over subarray. The input and inout arrays
c are from MPI's viewpoint single pointers to the start of the 
c arrays to be added.
      real invec(*),inoutvec(*)
c      integer ilen,itype
c Common for passing the dimensional structures. Needs to be
c set in the calling routine of the REDUCE that references this
c operator use or by using mpisubopcreate:
      include 'mditcom.f'
      external iteradd
      ipoint=0
      call mditerarg(iteradd,nasdims,iasfull,iasum2,ipoint,
     $        invec,inoutvec,dum3,dum4)
      end
c********************************************************************
      subroutine iteradd(inc,ipoint,indi,ndims,iused,
     $     a,b,c,d)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real a(*),b(*),c(*),d(*)
c This routine for use in mditerarg.
c Hardly anything is used. Only inc,ipoint,a,b.
c Add a(*) to b(*) and leave in b(*)
      ind=1+ipoint
      b(ind)=b(ind)+a(ind)
      inc=1
      end
c*******************************************************************
