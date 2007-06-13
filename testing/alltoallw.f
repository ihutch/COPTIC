c Test alltoallw
      include 'mpif.h'
      parameter (nproc=3,iblen=10)
      real u(iblen),v(iblen)
      integer isc(nproc),isd(nproc),ist(nproc)
      integer irc(nproc),ird(nproc),irt(nproc)
      
      call MPI_INIT(ierr)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if(nproc.ne.numprocs) stop 'wrong number of processes'
c initialize
      do i=1,nproc
         isc(i)=0
         isd(i)=i*4
         ist(i)=MPI_REAL
         irc(i)=0
c        ird(i)=(nproc+i)*4
        ird(i)=i*4
         irt(i)=MPI_REAL
      enddo

      do i=1,iblen
         u(i)=myid
         v(i)=myid
      enddo

c copy from 0 to 1, and from 1 to 0 
      if(myid.eq.0)then
c         isc(2)=1
         irc(2)=1
      endif
      if(myid.eq.1)then
c         irc(1)=1
         isc(1)=1
      endif


      write(*,101)(i,isc(i),isd(i),mod(ist(i),10000),
     $     irc(i),ird(i),mod(irt(i),10000),i=1,nproc)

      call MPI_ALLTOALLW(u,isc,isd,ist,u,irc,ird,irt,
     $     MPI_COMM_WORLD,ierr)

 101  format(7i8)
c      write(*,101)(isc(i),isd(i),mod(ist(i),10000),
c     $     irc(i),ird(i),mod(irt(i),1000),i=1,nproc)
      write(*,*)'u',u
      write(*,*)'v',v

      call MPI_FINALIZE()
      end
