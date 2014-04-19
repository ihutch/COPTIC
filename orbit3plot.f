c Packaged version of orbit plotting.
      subroutine orbitplot(ifull,iuds,u,phi,rc,rs)
      include 'ndimsdecl.f'
      integer ifull(ndimsmax),iuds(ndimsmax),itemp(ndimsmax)
      real u(ifull(1),ifull(2),ifull(3))
      integer ifmax
      parameter (ifmax=100,Li=ifmax)
      real uplot(ifmax,ifmax)
      character cwork(ifmax,ifmax)

      include 'meshcom.f'
      include 'partcom.f'
c Silence warnings with spurious assigment.
      zclv=phi
c Calculate some stuff for contour plot.
      idf=3
      id1=mod(idf,3)+1
      id2=mod(idf+1,3)+1
      ifixed=iuds(idf)/2
      if(.true.)then
         do i=1,iuds(id1)
            do j=1,iuds(id2)
               itemp(idf)=ifixed
               itemp(id1)=i
               itemp(id2)=j
               uplot(i,j)=u(itemp(1),itemp(2),itemp(3))
            enddo
c               write(*,'(10f8.4)')(uplot(i,j),j=1,iuds(id2))
         enddo
      endif
      call dashset(0)
      nf1=iuds(id1)
      nf2=iuds(id2)
c      call pltinit(-rs,rs,-rs,rs)
      call pltinit(xmeshstart(1),xmeshend(1),xmeshstart(2),xmeshend(2))
c Contour without labels, with coloring, using vector coordinates.
      zclv=20.
      icl=0
      call contourl(uplot,cwork,Li,nf1,nf2,zclv,icl,
     $        xn(ixnp(id1)+1),xn(ixnp(id2)+1),17)
      call color(15)
      call axis()
      call color(13)
      call circleplot(0.,0.,rs)
      call circleplot(0.,0.,rc)
      do kk=1,norbits
         call color(kk)
         call polyline(xorbit(1,kk),yorbit(1,kk),iorbitlen(kk))
         call polymark(xorbit(1,kk),yorbit(1,kk),iorbitlen(kk),3)
      enddo
      call pltend()
c      write(*,*)'Returning from orbit3plot.'
      end
c******************************************************************
c Called in the middle of a 3-d plot by cijplot.
      subroutine orbit3plot()
      include 'ndimsdecl.f'
      include 'partcom.f'
c      write(*,*)'norbits,length=',norbits,iorbitlen(1)
      do kk=1,norbits
         call color(kk)
         call poly3line(xorbit(1,kk),yorbit(1,kk),zorbit(1,kk),
     $        iorbitlen(kk))
         call poly3mark(xorbit(1,kk),yorbit(1,kk),zorbit(1,kk),
     $        iorbitlen(kk),1)
         if(iorbitlen(kk).gt.0)call drcstr('End')         
      enddo

      end
