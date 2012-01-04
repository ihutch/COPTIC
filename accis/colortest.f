c Test of colors.
      integer r1,r2,g1,g2,b1,b2
      real x(2),y(2)
      data r1,g1,b1,r2,g2,b2/10,10,10,64000,20000,6000/
      
      write(*,*)'Starting accisgradinit'
      write(*,*)'Completed accisgradinit'
      call pltinit(0.,1.,0.,1.)
c      call accisgradinit(r1,g1,b1,r2,g2,b2)
      call accisgradinit(-32000,-64000,0,128000,64000,192000)
      do i=0,255
         x(1)=1.*i/255.
         x(2)=1.*i/255.
         y(1)=0.
         y(2)=1.
         call gradcolor(i)
         call polyline(x,y,2)
      enddo
      call pltend()
      end
