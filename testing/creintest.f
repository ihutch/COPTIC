      program creintest
c Test the reinjection scheme by forming cartesian distributions.
      include 'plascom.f'
      include 'meshcom.f'
      include 'creincom.f'

      parameter (ndiag=100,mdims=3)
      real fv(ndiag,mdims)
      real px(ndiag,mdims)
      real diagv(ndiag)
      real diagx(ndiag,mdims)
      real xr(3*mdims)
      character*100 string
      logical normal,lwork

c defaults
      vd=1.2
      Ti=.1
      vrange=5.
      nin=1000000
c      normal=.false.
      normal=.true.

c Initialize random number generator
      v=ran1(-6)

c Set up mesh data.
      do id=1,mdims
         xmeshstart(id)=-1.
         xmeshend(id)=1.
      enddo

      do i=1,ndiag
         diagv(i)=vrange*(-1.+2.*(i-0.5)/ndiag)
         do id=1,mdims
            fv(i,id)=0.
            px(i,id)=0.
            diagx(i,id)=xmeshstart(id)+(i-0.5)*
     $           (xmeshend(id)-xmeshstart(id))/(ndiag)
         enddo
      enddo

      write(*,*)'Starting reinjects'
      nrein=0
      do k=1,nin
         call reinject(xr,nrein)
c         write(*,*)'Returned from reinject',xr
c Assign velocities to bins.
         do id=1,mdims
            lwork=id.eq.abs(idrein)
            if(.not.normal)lwork=.not.lwork
            if(lwork)then
               v=xr(mdims+id)
               v=sign(min(vrange,abs(v)),v)
               ibin=nint(0.5*(ndiag)*(1.+0.99999*v/vrange)+0.5)
               if(ibin.lt.1.or.ibin.gt.ndiag)
     $              write(*,*)k,nin,id,' ibin',ibin,v
               fv(ibin,id)=fv(ibin,id)+1
c Assign positions to bins
               x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
               ibin=nint(0.50000+x*(ndiag-.00000))
               if(ibin.lt.1 .or. ibin.gt.ndiag)then
                  write(*,*)'ibin=',ibin,id,x,xr(id)
               else
                  px(ibin,id)=px(ibin,id)+1.
               endif
               if(normal.and.(ibin.eq.1 .and. v.lt.0.
     $              .or. ibin.eq.ndiag .and. v.gt.0.))then
                  write(*,*)'Reinject velocity sign wrong',
     $                 (xr(kk),kk=1,2*mdims)
               endif
            endif
         enddo
      enddo

      write(*,*)'Finished',nin,' injects.'
      call ticnumset(10)
      do id=1,mdims
         call autoplot(diagv,fv(1,id),ndiag)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('velocity',string(1:30))
         call pltend()
      enddo
      do id=1,mdims
         call autoplot(diagx(1,id),px(1,id),ndiag)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('position',string(1:30))
         call pltend()
      enddo

      end
