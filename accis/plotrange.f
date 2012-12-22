      parameter (ndth=100)
      real vfv(ndth),fv(ndth)

         vrange=5.
         dvf=2.*vrange/ndth
         do i=1,ndth
            vfv(i)=(-vrange+((i-1)+.5)*dvf)
            fv(i)=vfv(i)
         enddo

         call automark(vfv,fv,ndth,1)

         call pltend()

         end
