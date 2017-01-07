c Do parallel writing as a test of file systems access etc.
      include 'mpif.h'
c Number of integers to write:
      parameter (iblen=100000)
      character*80 filename

      if(iargc().eq.1)then
         call getarg(1,filename)
      else
         filename=' '
      endif
      
      call MPI_INIT(ierr)

      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      ilen=nbcat(filename,'T1Written.')
      write(filename(lentrim(filename)+1:),'(i4.4)')myid
      write(*,*)'Writing to file: ',filename
      
      open(1,file=filename,status='unknown',err=101)
      close(1,status='delete')
      open(1,file=filename,status='new',form='formatted',err=101)
      do i=1,iblen
         write(1,*)i
      enddo
      close(1)
      
      call MPI_FINALIZE()
      call exit
      
 101  continue
      write(*,*)'Error opening file:',filename
      close(1,status='delete')
      end
