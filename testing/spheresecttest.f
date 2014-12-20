c Test the spheresect code to understand what can cause particle leakage.
c My guess is that there are problems associated with steps that are
c very close to the sphere/circle. 
      program spheresecttest
      parameter (ndims=3)
      real xp(ndims),xd(ndims),xp1(ndims),xp2(ndims),xc(ndims),rc(ndims)

      ntheta=1000
      nphi=1000
      tiny=1.e-6
c      tiny=0.
      sumtol=1.e-8
c      halfmin=.4999999
      halfmin=.5
c Hope the problem shows up with zero center and 1 radius.
      do id=1,ndims
         xc(id)=0.
         rc(id)=1.+tiny
      enddo
c Length of step.
      dd=1.e-1
c Cylinder really:
c      ida=3
      ida=0
      do j=1,nphi
         phi=2.*3.1415926*float(i-1)/nphi
         do i=1,ntheta
            theta=3.1415926*float(i-1)/ntheta
c Find a place on the unit sphere 
            xp(3)=cos(theta)
            s=sin(theta)
            xp(1)=s*cos(phi)
            xp(2)=s*sin(phi)
c Now examine steps from/to this.
c Get a random vector.
            do ipm=1,2
c off or on the sphere
            do id=1,ndims
               call ranlux(xd(id),1)
            enddo
            do id=1,ndims
               xp1(id)=xp(id)+(2-ipm)*dd*xd(id)
               xp2(id)=xp(id)+(ipm-1)*dd*xd(id)
            enddo
            call spheresect(ndims,ida,xp1,xp2,xc,rc,f11,f21,sd1,C1,D1)
            call spheresect(ndims,ida,xp2,xp1,xc,rc,f12,f22,sd2,C2,D2)
c Test that swapping xp1 and xp2 does not produce a different result.
c The test for this was not correct.
            if(.false.)then
               write(*,'(a,2i3,6f7.2,2f8.4)')'sd error j,i,sd1,sd2=',j,i
     $              ,sd1,sd2
               write(*,*)'fs',f11,f21,f12,f22
               write(*,*)'CDs',C1,D1,C2,D2
            endif
c When we swap xp1 and xp2, f -> 1-f. 
c The fs are ordered to start closest to 0. 
c Asymmetries arise because C and D are not treated on an equal footing.
c So rounding is able to caused differences.
            if((abs(f22-halfmin).lt..5.or.abs(f21-halfmin).lt..5))then
c Two valid intersection case
               if(abs(f11+f22-1).gt.sumtol)then
                  write(*,'(a,2i3,4f7.2,2f10.6)')'benign2 j,i,fs,sums',j
     $                 ,i,f11,f21,f12,f22,f11+f22-1.,f21+f12-1.
                  write(*,*)f11,f21,f12,f22
                  write(*,*)C1,D1,C2,D2
               endif
               if(sd1.eq.0.or.sd2.eq.0.)write(*,*)'sd error',sd1,sd2
            elseif((abs(f11-halfmin).lt..5.or.abs(f12-halfmin).lt..5)
     $              .and.abs(f21 -halfmin).gt..5.and.abs(f22
     $              -halfmin).gt.5)then
c One valid intersection case
               if(abs(f11+f12-1).gt.sumtol)then
                  write(*,'(a,2i3,4f7.2,2f10.6)')'benign j,i,fs,sums',j
     $                 ,i,f11,f21,f12,f22,f11+f12-1.
C,f21+f12-1.
                  write(*,*)f11,f21,f12,f22
                  write(*,*)C1,D1,C2,D2
               endif
               if(sd1.eq.0.or.sd2.eq.0.)write(*,*)'sd error',sd1,sd2
            else
c No intersections in interval. But there might be outside; so sd is not
c necessarily zero.
            endif
            if(f21.eq.0.and.f11.ne.0 .or. f22.eq.0.and.f12.ne.0.)then
c            if(.false.)then
               write(*,*)'21 fs=',f11,f21,f12,f22
               write(*,*)C1,D1,C2,D2
            endif
            enddo
         enddo

c Put cylindrical test here.
      enddo
      write(*,*)'Finished. ntheta,nphi,dd=',ntheta,nphi,dd

c 100  format(a,6f10.4)
      end
