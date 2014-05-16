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
      include 'ndimsdecl.f'
      include '3dcom.f'
c ALLREDUCE communicates the result to all processes.      

c Have to reduce data for all objects and all quantities.
c Because of the address definition:
c               if(i.gt.1)nf_address(i,j,k)=
c     $              nf_address(i-1,j,k)+nf_posno(i-1,j)
c The quantities are adjacent to one another in ff_data. Therefore
c we don't have to iterate over quantities, just pass more data, 
c encompassing them all at once.
      do io=1,mf_obj
         nquant=nf_posno(1,io)*mf_quant(io)
         call MPI_ALLREDUCE(MPI_IN_PLACE
     $        ,ff_data(nf_address(1,io,nf_step))
     $        ,nquant,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE
     $        ,partforce(1,io,nf_step)
     $        ,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE
     $        ,colnforce(1,io,nf_step)
     $        ,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      enddo
c      write(*,*)'Total flux number',sum
c      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
c     $     i=1,nf_posno(1,1))

      end
c********************************************************************
      subroutine psumreduce(psum,nrein,phirein,numprocs,ndims,ifull,iuds
     $     ,iLs)
      integer nrein,ndims,numprocs
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real psum(*),phirein
      include 'mpif.h'
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
c         write(*,*)'iaddtype,iaddop=',iaddtype,iaddop
         do i=1,ndims
            iporig=iporig+iLs(i)
         enddo
         lfirst=.false.
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
c********************************************************************
c Sum reduce the diagnostic arrays. 
      subroutine diagreduce(diagsum,ndims,ifull,iuds,iLs,ndiags)
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real diagsum(*)
c Number of diagnostics ndiags
      integer ndiags
      include 'mpif.h'
c      include 'ndimsdecl.f'
c      include 'partcom.f'
c Operator
      external addarray_MPI

      logical lfirst
      integer iaddtype,iaddop,iporig
      data lfirst/.true./
      save lfirst,iaddtype,iaddop,iporig

      if(lfirst)then
c         write(*,*)'Calling mpiopcreate',ndiags
c Create addtype and operator for reduce sum; full array.
         call mpiopcreate(ndims,ifull,iuds,addarray_MPI,
     $        iaddtype,iaddop)
c         write(*,*)'Returned from mpiopcreate',iaddtype
         lfirst=.false.
      endif
c We do each diagnostic reduce as a separate call so that we can
c use variable numbers of diagnostics.
      do k=1,ndiags
         iporig=1+(k-1)*iLs(ndims+1)
         call MPI_ALLREDUCE(MPI_IN_PLACE,diagsum(iporig),1,iaddtype,
     $     iaddop,MPI_COMM_WORLD,ierr)
      enddo
      end

c********************************************************************
c********************************************************************
      subroutine ptdiagreduce()
c Reduce the particle distribution diagnostics accumulations.
      include 'mpif.h'
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      integer nfv,nfsv
      parameter(nfv=2*nptdiag*mdims,
     $     nfsv=nsbins*mdims*(nsub_tot+1)+nsub_tot)
c Reduce the data for common /cartdiag/fv,px,... nfvaccum
c Only the fv,px are to be added together. And nfvaccum.
      call MPI_ALLREDUCE(MPI_IN_PLACE,fv,nfv,MPI_REAL
     $     ,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nfvaccum,1,MPI_INTEGER
     $     ,MPI_SUM,MPI_COMM_WORLD,ierr)
c Reduce the data for common /subdiag/ ... fsv,fvx,denfvx
      call MPI_ALLREDUCE(MPI_IN_PLACE,fsv,nfsv,MPI_REAL 
     $     ,MPI_SUM,MPI_COMM_WORLD,ierr)
      end
c********************************************************************
      subroutine minmaxreduce(mdims,xnewlim)
c Inefficient max and min of xnewlims structure.
      integer mdims
      real xnewlim(2,mdims)
      include 'mpif.h'
      do id=1,mdims
         call MPI_ALLREDUCE(MPI_IN_PLACE,xnewlim(1,id),1,MPI_REAL
     $        ,MPI_MIN,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xnewlim(2,id),1,MPI_REAL
     $        ,MPI_MAX,MPI_COMM_WORLD,ierr)
      enddo
      end
c**************************************************************
      subroutine mpisubopcreate(nrd,ifull,iuds,iopfun,itype,iaddop)
c For an nrd (IN) dimensional structure 
c described by ifull (IN), iuds (IN)
      integer nrd
      integer iuds(nrd),ifull(nrd)
c Create a datatype itype (OUT) and an operator iaddop (OUT)

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
         iasuds(i)=iuds(i)
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
c This is just like mpisubopcreate, except (1) for full blocks,
c not omitting the boundaries, (2) not setting iascommons.
      subroutine mpiopcreate(nrd,ifull,iuds,iopfun,itype,iaddop)
c For an nrd (IN) dimensional structure 
c described by ifull (IN), iuds (IN)
      integer nrd
      integer iuds(nrd),ifull(nrd)
c Create a datatype itype (OUT) and an operator iaddop (OUT)

c for doing MPI_REDUCEs on a subarray described by iuds
c (Its boundaries are NOT omitted in the reduce.)
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
c These are needed for the blockcreate. We add 2 to isides here to
c make the created structure appropriate for full arrays. Inside
c the bbdyblockcreate code is the presumption that boundaries are 
c omitted. We overrule that by this expedient.
         iside(1,i)=iuds(i)+2
         iside(2,i)=iuds(i)+2
         iLs(i+1)=iLs(i)*ifull(i)
c These are needed for the ioperator usage. They ought already to have
c been set by the mpisubopcreate call. However, it does no harm and
c ensures that this call could be used alone if desired to set them
c here.
         iasfull(i)=ifull(i)
         iasuds(i)=iuds(i)
         iasum2(i)=iuds(i)-2
      enddo
      nasdims=nrd
c iasfull etc don't need to be set for this call.
      call MPI_OP_CREATE(iopfun,.false.,iaddop,ierr)
      call bbdyblockcreate(nrd,ktype,iLs,iside,iSIZEOFREAL)
      itype=ktype(2**nrd)

c MPI_IN_PLACE is ok for ALLREDUCE. Not for REDUCE.
c Subsequent calls will be of the form:
c      call MPI_ALLREDUCE(MPI_IN_PLACE,u(1,1,1),1,itype,
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
      subroutine iteradd(inc,ipoint,indi,ndims,iLs,iused,
     $     a,b)
c ,c,d)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real a(*),b(*)
c ,c(*),d(*)
c This routine for use in mditerarg.
c Hardly anything is used. Only inc,ipoint,a,b.
c Add a(*) to b(*) and leave in b(*)
      ind=1+ipoint
      b(ind)=b(ind)+a(ind)
      inc=1
      end
c*******************************************************************
c This is just like addsubarray: addition reduce over array,
c except that it uses iasuds for whole array, instead of iasum2, which
c excludes the boundaries.
      subroutine addarray_MPI(invec,inoutvec)
c actually (invec,inoutvec,ilen,itype) but ilen,itype not used.
c MPI_sum type function over subarray. The input and inout arrays
c are from MPI's viewpoint single pointers to the start of the 
c arrays to be added.
      real invec(*),inoutvec(*)
      include 'mditcom.f'
      external iteradd
      ipoint=0
      call mditerarg(iteradd,nasdims,iasfull,iasuds,ipoint,
     $        invec,inoutvec,dum3,dum4)
      end
c**************************************************************
