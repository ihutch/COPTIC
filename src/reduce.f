!***********************************************************************
! This file contains the extra routines that are needed for parallel
! processing the particles. (mpibbdy contains others for fields)
! 
! For a non-mpi installation, the subroutines should be dummies, simply
! returning. Thus the particle code can be immediately un-MPIed.
!
!***********************************************************************
      subroutine fluxreduce()
      include 'mpif.h'
      include 'ndimsdecl.f'
      include '3dcom.f'
! ALLREDUCE communicates the result to all processes.      

! Have to reduce data for all objects and all quantities.
! Because of the address definition:
!               if(i.gt.1)nf_address(i,j,k)=
!     $              nf_address(i-1,j,k)+nf_posno(i-1,j)
! The quantities are adjacent to one another in ff_data. Therefore
! we don't have to iterate over quantities, just pass more data, 
! encompassing them all at once.
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
!      write(*,*)'Total flux number',sum
!      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
!     $     i=1,nf_posno(1,1))

      end
!********************************************************************
      subroutine psumreduce(psum,nrein,phirein,ndims,ifull,iuds
     $     ,iLs)
      integer nrein,ndims
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real psum(*),phirein
      include 'mpif.h'
      include 'myidcom.f'
! Operator
      external addsubarray_MPI

      logical lfirst
      integer iaddtype,iaddop,iporig
      data lfirst/.true./
      save lfirst,iaddtype,iaddop,iporig

      if(lfirst)then
! Create addtype and operator for reduce sum.
         call mpisubopcreate(ndims,ifull,iuds,addsubarray_MPI,
     $        iaddtype,iaddop)
         iporig=1
!         write(*,*)'iaddtype,iaddop=',iaddtype,iaddop
         do i=1,ndims
            iporig=iporig+iLs(i)
         enddo
         lfirst=.false.
      endif
! Pass the (2,2,...) element of psum.
      call MPI_ALLREDUCE(MPI_IN_PLACE,psum(iporig),1,iaddtype,
     $     iaddop,MPI_COMM_WORLD,ierr)
! Reduce also the nrein and phirein.
      call MPI_ALLREDUCE(MPI_IN_PLACE,nrein,1,MPI_INTEGER,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,phirein,1,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      phirein=phirein/nprocs
!      write(*,*)'numprocs,phirein',numprocs,phirein

      end

!********************************************************************
!********************************************************************
! Sum reduce the diagnostic arrays. 
      subroutine diagreduce(diagsum,ndims,ifull,iuds,iLs,ndiags)
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real diagsum(*)
! Number of diagnostics ndiags
      integer ndiags
      include 'mpif.h'
!      include 'ndimsdecl.f'
!      include 'partcom.f'
! Operator
      external addarray_MPI

      logical lfirst
      integer iaddtype,iaddop,iporig
      data lfirst/.true./
      save lfirst,iaddtype,iaddop,iporig

      if(lfirst)then
!         write(*,*)'Calling mpiopcreate',ndiags
! Create addtype and operator for reduce sum; full array.
         call mpiopcreate(ndims,ifull,iuds,addarray_MPI,
     $        iaddtype,iaddop)
!         write(*,*)'Returned from mpiopcreate',iaddtype
         lfirst=.false.
      endif
! We do each diagnostic reduce as a separate call so that we can
! use variable numbers of diagnostics.
      do k=1,ndiags
         iporig=1+(k-1)*iLs(ndims+1)
         call MPI_ALLREDUCE(MPI_IN_PLACE,diagsum(iporig),1,iaddtype,
     $     iaddop,MPI_COMM_WORLD,ierr)
      enddo
      end

!********************************************************************
!********************************************************************
      subroutine ptdiagreduce()
! Reduce the particle distribution diagnostics accumulations.
      include 'mpif.h'
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      integer nfv,nfsv
      parameter(nfv=2*nptdiagmax*ndimsmax)
      parameter(nfsv=nsbins*ndimsmax*(nsub_tot+1)+nsub_tot)
      parameter(nf2sv=nsbins*ndimsmax*(nsub_tot*(1+nsbins)+1)+nsub_tot)
! Reduce the data for common /cartdiag/fv,px,... nfvaccum
! Only the fv,px are to be added together. And nfvaccum.
      call MPI_ALLREDUCE(MPI_IN_PLACE,fv,nfv,MPI_REAL
     $     ,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nfvaccum,1,MPI_INTEGER
     $     ,MPI_SUM,MPI_COMM_WORLD,ierr)
! Reduce the data from common /subdiag/  fsv,fvx,denfvx [,f2vx] 
! Assuming they are contiguous.
      if(nptdiag.eq.nsbins)then
! We need all the data.
         call MPI_ALLREDUCE(MPI_IN_PLACE,fsv,nf2sv,MPI_REAL 
     $        ,MPI_SUM,MPI_COMM_WORLD,ierr)
      else
! Only the 1-D data needs to be reduced.
         call MPI_ALLREDUCE(MPI_IN_PLACE,fsv,nfsv,MPI_REAL 
     $        ,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
      end
!********************************************************************
      subroutine minmaxreduce(mdims,xnewlim)
! Inefficient max and min of xnewlims structure.
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
!**************************************************************
      subroutine mpisubopcreate(nrd,ifull,iuds,iopfun,itype,iaddop)
! For an nrd (IN) dimensional structure 
! described by ifull (IN), iuds (IN)
      integer nrd
      integer iuds(nrd),ifull(nrd)
! Create a datatype itype (OUT) and an operator iaddop (OUT)

! for doing MPI_REDUCEs on a subarray described by iuds
! (Its boundaries are omitted in the reduce.)
! The operation function is iopfun (EXTERNAL).
! iopfun MUST be declared external in the calling routine.
!
! Local storage
      parameter (ims=10)
      integer iside(2,ims),ktype(2**(ims+1)),iLs(ims+1)
      
      include 'mpif.h'
! Common for passing the dimensionals structures:
      include 'mditcom.f'
      external iopfun

      if(nrd.gt.ims)stop 'mpisubopcreate dimension too large'
      iLs(1)=1
      do i=1,nrd
! These are needed for the blockcreate
         iside(1,i)=iuds(i)
         iside(2,i)=iuds(i)
         iLs(i+1)=iLs(i)*ifull(i)
! These are needed for the ioperator usage.
         iasfull(i)=ifull(i)
         iasuds(i)=iuds(i)
         iasum2(i)=iuds(i)-2
      enddo
      nasdims=nrd
! iasfull etc don't need to be set for this call.
      call MPI_OP_CREATE(iopfun,.false.,iaddop,ierr)
      call bbdyblockcreate(nrd,ktype,iLs,iside,iSIZEOFREAL)
      itype=ktype(2**nrd)

! MPI_IN_PLACE is ok for ALLREDUCE. Not for REDUCE.
! Subsequent calls will be of the form:
!      call MPI_ALLREDUCE(MPI_IN_PLACE,u(2,2,2),1,itype,
!     $     iaddop,icommcart,ierr)
      end

!**************************************************************
! This is just like mpisubopcreate, except (1) for full blocks,
! not omitting the boundaries, (2) not setting iascommons.
      subroutine mpiopcreate(nrd,ifull,iuds,iopfun,itype,iaddop)
! For an nrd (IN) dimensional structure 
! described by ifull (IN), iuds (IN)
      integer nrd
      integer iuds(nrd),ifull(nrd)
! Create a datatype itype (OUT) and an operator iaddop (OUT)

! for doing MPI_REDUCEs on a subarray described by iuds
! (Its boundaries are NOT omitted in the reduce.)
! The operation function is iopfun (EXTERNAL).
! iopfun MUST be declared external in the calling routine.
!
! Local storage
      parameter (ims=10)
      integer iside(2,ims),ktype(2**(ims+1)),iLs(ims+1)
      
      include 'mpif.h'
! Common for passing the dimensionals structures:
      include 'mditcom.f'
      external iopfun

      if(nrd.gt.ims)stop 'mpisubopcreate dimension too large'
      iLs(1)=1
      do i=1,nrd
! These are needed for the blockcreate. We add 2 to isides here to
! make the created structure appropriate for full arrays. Inside
! the bbdyblockcreate code is the presumption that boundaries are 
! omitted. We overrule that by this expedient.
         iside(1,i)=iuds(i)+2
         iside(2,i)=iuds(i)+2
         iLs(i+1)=iLs(i)*ifull(i)
! These are needed for the ioperator usage. They ought already to have
! been set by the mpisubopcreate call. However, it does no harm and
! ensures that this call could be used alone if desired to set them
! here.
         iasfull(i)=ifull(i)
         iasuds(i)=iuds(i)
         iasum2(i)=iuds(i)-2
      enddo
      nasdims=nrd
! iasfull etc don't need to be set for this call.
      call MPI_OP_CREATE(iopfun,.false.,iaddop,ierr)
      call bbdyblockcreate(nrd,ktype,iLs,iside,iSIZEOFREAL)
      itype=ktype(2**nrd)

! MPI_IN_PLACE is ok for ALLREDUCE. Not for REDUCE.
! Subsequent calls will be of the form:
!      call MPI_ALLREDUCE(MPI_IN_PLACE,u(1,1,1),1,itype,
!     $     iaddop,icommcart,ierr)
      end

!**************************************************************
! These two routines provide facility for addition reduce 
! over subarray.
      subroutine addsubarray_MPI(invec,inoutvec)
! actually (invec,inoutvec,ilen,itype) but ilen,itype not used.
! MPI_sum type function over subarray. The input and inout arrays
! are from MPI's viewpoint single pointers to the start of the 
! arrays to be added.
      real invec(*),inoutvec(*)
!      integer ilen,itype
! Common for passing the dimensional structures. Needs to be
! set in the calling routine of the REDUCE that references this
! operator use or by using mpisubopcreate:
      include 'mditcom.f'
      external iteradd
      ipoint=0
      call mditerarg(iteradd,nasdims,iasfull,iasum2,ipoint,
     $        invec,inoutvec,dum3,dum4)
      end
!********************************************************************
      subroutine iteradd(inc,ipoint,indi,ndims,iLs,iused,
     $     a,b)
! ,c,d)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real a(*),b(*)
! ,c(*),d(*)
! This routine for use in mditerarg.
! Hardly anything is used. Only inc,ipoint,a,b.
! Add a(*) to b(*) and leave in b(*)
      ind=1+ipoint
      b(ind)=b(ind)+a(ind)
      inc=1
      end
!*******************************************************************
! This is just like addsubarray: addition reduce over array,
! except that it uses iasuds for whole array, instead of iasum2, which
! excludes the boundaries.
      subroutine addarray_MPI(invec,inoutvec)
! actually (invec,inoutvec,ilen,itype) but ilen,itype not used.
! MPI_sum type function over subarray. The input and inout arrays
! are from MPI's viewpoint single pointers to the start of the 
! arrays to be added.
      real invec(*),inoutvec(*)
      include 'mditcom.f'
      external iteradd
      ipoint=0
      call mditerarg(iteradd,nasdims,iasfull,iasuds,ipoint,
     $        invec,inoutvec,dum3,dum4)
      end
!**************************************************************
! Generic All Reduce Sum Inplace Real. 
! For contiguous real array (of any actual shape).
      subroutine mpiallreducesum(array,nsize,ierr)
      include 'mpif.h'
      include 'myidcom.f'
      integer nsize,ierr
      real array(*)
      call MPI_ALLREDUCE(MPI_IN_PLACE,array,nsize,MPI_REAL,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)
      end
!**************************************************************
! Generic All Reduce Max Inplace Real. 
! For contiguous real array (of any actual shape).
      subroutine mpiallreducemax(array,nsize,ierr)
      include 'mpif.h'
      include 'myidcom.f'
      integer nsize,ierr
      real array(*)
      call MPI_ALLREDUCE(MPI_IN_PLACE,array,nsize,MPI_REAL,MPI_MAX,
     $     MPI_COMM_WORLD,ierr)
      end
!**************************************************************
! Generic Reduce Sum Inplace Real. Single destination: iddest 
! For contiguous real array (of any actual shape).
      subroutine mpireducesum(array,nsize,iddest,ierr)
      include 'mpif.h'
      include 'myidcom.f'
      integer nsize,ierr,iddest
      real array(nsize)
      if(myid.eq.iddest)then
         call MPI_REDUCE(MPI_IN_PLACE,array,nsize,MPI_REAL,MPI_SUM,
     $        iddest,MPI_COMM_WORLD,ierr)
      else
! When I am not the destination, the destination array is irrelevant.
         call MPI_REDUCE(array,array,nsize,MPI_REAL,MPI_SUM,
     $        iddest,MPI_COMM_WORLD,ierr)
      endif
      end
!**************************************************************
