c Boxcar average a trace.
c Make the traceave(i) equal to the average from i-nb to i+nb of trace.
c If nb=0, just transfer
      subroutine boxcarave(nt,nb,trace,traceave)
      integer nt,nb
      real trace(nt),traceave(nt)

      if(nb.lt.0)return
      do i=1,nt
         accum=0.
         nac=0
         do j=max(1,i-nb),min(nt,i+nb)
            accum=accum+trace(i+j)
            nac=nac+1
         enddo
         traceave(i)=accum/nac
      enddo

      end
