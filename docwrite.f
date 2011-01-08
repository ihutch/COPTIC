      program docwrite
c Reads a plain text file from disk and writes a fortran routine that
c will print the same text out.
      character*100 infile,outfile,routinename
      parameter (ncmax=200)
      character*(ncmax) aline
      parameter (mlen=54)
      if(iargc().lt.2)then
         write(*,*)'Usage: docwrite inputfile.txt outputfile.f [name]'
         call exit(1)
      else
         call getarg(1,infile)
         call getarg(2,outfile)
         if(iargc().eq.3)call getarg(3,routinename)
      endif

      open(23,file=infile,status='old',form='formatted',err=101)
      open(22,file=outfile,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=outfile,status='new',form='formatted',err=102)

      if(iargc().eq.3)write(22,'(3a)')'      subroutine '
     $     ,routinename(1:lentrim(routinename)),'()'
      write(22,'(a)')'8     format(9a)'
 12   read(23,'(a)',err=103,end=104)aline
c Double single quote characters.
      ic=ncmax
      do i=lentrim(aline),1,-1
         if(aline(i:i).eq.'''')then
            aline(ic:ic)=''''
            ic=ic-1
         endif
         aline(ic:ic)=aline(i:i)
         ic=ic-1
         if(ic.lt.1)then
            write(*,*)'Line too long:',aline
            stop
         endif
      enddo
      do i=1,ncmax-ic
         aline(i:i)=aline(ic+i:ic+i)
         aline(ic+i:ic+i)=' '
      enddo
c End of quote fixing.
      index=1
      len=lentrim(aline)
      if(len.eq.0)write(22,'(''      write(*,8)'',a,a,a)')
 11   continue
      if(len.gt.0)then
         if(index.eq.1)then
            write(22,'(''      write(*,8)'',a,a,a)')
     $           '''',aline(index:index+min(len,mlen)-1),''''
         else
            write(22,'(a,a,a)')'     &         ,'''
     $           ,aline(index:index+min(len,mlen)-1),''''
         endif
         len=len-mlen
         index=index+mlen
         goto 11
      endif
      goto 12
 104  continue
      write(22,'(a)')'      end'

      stop

 101  write(*,*)'Could not open infile: ',infile(1:lentrim(infile))
      stop
 102  write(*,*)'Error opening outfile:',outfile(1:lentrim(outfile))
      stop
 103  write(*,*)'Error reading infile:',infile(1:lentrim(infile))
      end
