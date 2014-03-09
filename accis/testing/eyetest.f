c Test the geteye response
      i=0
      call pltinit(0.,1.,0.,1.)
 1    continue
      i=i+1
      write(*,*)'Call number',i
      call eye3d(isw)
      write(*,*)'Returned',isw
      goto 1
      end
