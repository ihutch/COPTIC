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
c         v(i)=myid
      enddo

c copy from 0 to 1, and from 1 to 0 
      if(myid.eq.0)then
         isc(2)=1
         irc(2)=1
         irc(3)=1
      endif

      if(myid.eq.1)then
         irc(1)=1
         isc(1)=1
      endif

      if(myid.eq.2)then
         isc(1)=1
      endif

      write(*,201)(i,isc(i),isd(i),mod(ist(i),10000),
     $     irc(i),ird(i),mod(irt(i),10000),i=1,nproc)

 201  format(7i8)
      call MPI_ALLTOALLW(u,isc,isd,ist,u,irc,ird,irt,
     $     MPI_COMM_WORLD,ierr)
