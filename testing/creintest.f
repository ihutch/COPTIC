      program creintest
c Test the reinjection scheme by forming cartesian distributions.
      include '../plascom.f'
      include '../meshcom.f'
      include '../creincom.f'
      include '../colncom.f'
      include '../partcom.f'

      parameter (ndiag=100,mdims=3)
      real fv(ndiag,mdims)
      real fvc(ndiag,mdims)
      real px(ndiag,mdims)
      real diagv(ndiag)
      real diagx(ndiag,mdims)
      integer ipinj(mdims)
      real xr(3*mdims)
      character*100 string,string2
      logical normal,lwork

c defaults (generally vary for testing)
      vd=.4
      vneutral=0.
      colntime=100.
c      colntime=0.
      Eneutral=0.
c      if(colntime.ne.0.)Eneutral=(vd-vneutral)/colntime
      Ti=.1
      Tneutral=.1
      vrange=3.
      nin=1000000
c Whether we plot the normal or tangential velocity distributions.
c      normal=.false.
      normal=.true.

c Initialize random number generator change (negative) value for different
c random number selections.
      v=ran1(-100)

c Set up mesh data.
      do id=1,mdims
         ipinj(id)=0.
         xmeshstart(id)=-1.
         xmeshend(id)=1.
         ipartperiod(id)=0
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

      if(Eneutral.ne.0.)then
         call colninit(0)
         call colreinit()
      endif

      write(*,*)'Starting reinjects'
      caverein=0.
      nrein=0
      do k=1,nin
         call reinject(xr,nrein,caverein)
c         write(*,*)'Returned from reinject',xr
c Assign velocities to bins.
         do id=1,mdims
c idrein is set by reinject as the dimension on which reinjected.
c            lwork=id.eq.abs(idrein)
            lwork=id.eq.abs(idrein)
            if(.not.normal)lwork=.not.lwork
c Accumulate reinjection data only in dimensional directions id that
c either coincide with the face dimension idrein (normal) or do not
c coincide (.not.normal). 
c For .not.normal, each particle contributes to two distributions.
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
               ipinj(id)=ipinj(id)+1
            endif
         enddo
      enddo

      write(*,*)'Finished',nin,' injects.'
c      write(*,*)'grein=',grein
c normalize the same as the sample.
         dv=diagv(2)-diagv(1)
         write(*,*)ipinj,dv
         do i=1,ndiag
            u=(diagv(i)-vneutral)/sqrt(2.*Tneutral)
            ud=(vd-vneutral)/sqrt(2.*Tneutral)
            do idrein=1,3
               if(normal)then
                  fvc(i,idrein)=ffdrein(diagv(i))*ipinj(idrein)*dv
     $                 /(grein(2*idrein-1)+grein(2*idrein))
               else
                  fvc(i,idrein)=fvdrein(diagv(i))*ipinj(idrein)*dv
               endif
            enddo
c These are the direct normalizations for the .not.normal case:
c            fvc(i,3)=fvcx(u,ud)*ipinj(3)*dv/sqrt(2.*Tneutral)
c            fvc(i,2)=exp(-diagv(i)**2/(2.*Ti))/sqrt(2.*3.1415926*Ti)
c     $           *ipinj(2)*dv
         enddo
c      write(*,*)(fvc(i,1),i=1,ndiag)
      call ticnumset(10)
      call pfset(3)
      do id=1,mdims
         call iwrite(id,iwidth,string2)
         string='distrib.'//string2
         write(*,*)string
         open(22,file=string,status='unknown')
         close(22,status='delete')
         open(12,file=string,status='new')
         write(12,*)ndiag
         write(12,'(2g16.6)')(diagv(j),fv(j,id),j=1,ndiag)
         close(12)
         call autoplot(diagv,fv(1,id),ndiag)
         call polymark(diagv,fv(1,id),ndiag,3)
         call color(6)
         call dashset(2)
         call polyline(diagv,fvc(1,id),ndiag)
         call dashset(0)
         call color(15)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('velocity',string(1:30))
         call pltend()
      enddo
      call pfset(0)

      do id=1,mdims
         call autoplot(diagx(1,id),px(1,id),ndiag)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('position',string(1:30))
         call pltend()
      enddo
     
      if(Eneutral.eq.0.)then
      do id=1,mdims
         call yautoplot(prein(0,id),ncrein+1)
         write(string,'(a,i3)')'prein, dimension',id
         call axlabels('',string(1:lentrim(string)))
         call pltend()
      enddo
      endif
      end
