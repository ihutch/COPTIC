c****************************************************************
c  Obsolete routines beyond this point.
c****************************************************************
c****************************************************************
c****************************************************************
      subroutine spheresect(id,ipm,ndims,indi,rc,xc,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c the ndims-dimensional spheroid of semi-radii rc(ndims), center xc.
c
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real xc(ndims),rc(ndims)
      real fraction,dp
      include 'meshcom.f'
      A=0.
      B=0.
      C=-1.
      D=-1.
c x1 and x2 are the coordinates in system in which sphere has radius 1.
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=(xn(ix1)-xc(i))/rc(i)
         x2=x1
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=(xn(ix2)-xc(i))/rc(i)
            A=A+(x2-x1)**2
            B=B+x1*(x2-x1)
            dp=abs(x2-x1)*rc(i)
         endif
         C=C+x1**2
         D=D+x2**2
      enddo
      fraction=1.
c This condition tests for a sphere crossing.
      if(D.ne.0. .and. D*C.le.0.)then
         if(C.lt.0.)then
            fraction=(-B+sqrt(B*B-A*C))/A
         elseif(C.gt.0.)then
            fraction=(-B-sqrt(B*B-A*C))/A
         else
            fraction=0.
         endif
      endif
      end
c****************************************************************
      subroutine cubesect1(id,ipm,ndims,indi,tc,bc,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c the ndims-dimensional coordinate aligned cuboid 
c of top coordinates tc(ndims), bottom coordinates bc.
c
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real tc(ndims),bc(ndims)
      real fraction,dp
      include 'meshcom.f'
c      write(*,*)'bc=',bc,' tc=',tc
      fraction=1.
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=(xn(ix1)-bc(i))
         side=tc(i)-bc(i)
c x1 and x2 are the coordinates relative to the bottom corner as origin.
c         x2=x1
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=(xn(ix2)-bc(i))
            dp=abs(x2-x1)
            if(x1.le.0.)then
               if(x2.gt.0.)then
                  fraction=x1/(x1-x2)
               endif
            elseif(x1.ge.side)then
               if(x2.lt.side)then
                  fraction=(x1-side)/(x1-x2)
               endif
            else
               if(x2.gt.side)then
                  fraction=(x1-side)/(x1-x2)
               elseif(x2.lt.0.)then
                  fraction=x1/(x1-x2)
               else
                  fraction=1.
                  return
               endif
            endif
         else
            if(x1.lt.0. .or. x1.gt.side)then
               fraction=1.
               return
            endif
         endif
      enddo
      end
c****************************************************************
c****************************************************************
      subroutine cubesect(id,ipm,ndims,indi,rc,xc,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c the ndims-dimensional coordinate aligned cuboid 
c of center coordinates xc(ndims), radius (face position) rc.
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real rc(ndims),xc(ndims)
      real fraction,dp
      include 'meshcom.f'
      fraction=1.
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=(xn(ix1)-xc(i))
         ri=abs(rc(i))
c x1 and x2 are the coordinates relative to the center.
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=(xn(ix2)-xc(i))
            dp=abs(x2-x1)
            if(x1.ge.ri.neqv.x2.ge.ri)then
c Crossed ri.
c This prefers the positive-ri face if both are crossed.               
               fraction=(x1-ri)/(x1-x2)
c               fraction=.999999*(x1-ri)/(x1-x2)
            elseif(x1.le.-ri.neqv.x2.le.-ri)then
c Crossed -ri
               fraction=(x1+ri)/(x1-x2)
c Not fixing cases of mesh clash. Better to get warning:
c               fraction=.999999*(x1+ri)/(x1-x2)
            else
c No parallel intersection.
c               fraction=1.
               return
            endif
         else
            if(abs(x1).ge.ri)then
c Outside box in orthogonal direction. No intersection.
               fraction=1.
               return
            endif
         endif
      enddo
      end
c****************************************************************
      subroutine cylsect(id,ipm,ndims,indi,rc,xc,ida,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c the ndims-dimensional coordinate-aligned cylinder
c with  radii rc(ndims), center xc, axial coordinate ida. 
c
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
      integer id,ipm,ndims
      integer indi(ndims),ida
      real xc(ndims),rc(ndims)
      real fraction,dp
      include 'meshcom.f'

      fraction=1.
      ix1=indi(ida)+ixnp(ida)+1
      z1=(xn(ix1)-xc(ida))
c z2 is only correct if id.eq.ida but it is only used then
      ix2=ix1+ipm
      z2=(xn(ix2)-xc(ida))
c (On the boundary face counts as outside.)
      if(abs(z1).ge.rc(ida))then
c z1 outside ends. if we are seeking radially no
         if(id.ne.ida)return
c Seeking axially. If z2 also outside, no.
         if(abs(z2).ge.rc(ida))return
c Here we are seeking axially and have found z2 inside z1 outside. 
      else
c z1 inside ends.
         if(id.ne.ida)then
c Searching radially and inside ends. Do circle intersection
            A=0.
            B=0.
            C=-1.
            D=-1.
            do k=1,ndims-1
               i=mod(ida+k-1,ndims)+1
               ix1=indi(i)+ixnp(i)+1
               x1=(xn(ix1)-xc(i))/rc(i)
               x2=x1
               if(i.eq.id)then
                  ix2=ix1+ipm
                  x2=(xn(ix2)-xc(i))/rc(i)
                  A=A+(x2-x1)**2
                  B=B+x1*(x2-x1)
                  dp=abs(x2-x1)*rc(i)
               endif
               C=C+x1**2
               D=D+x2**2
            enddo
            disc=B*B-A*C
            if(disc.gt.0)then
               if(C.lt.0.)then
                  fraction=(-B+sqrt(disc))/A
               elseif(C.gt.0.)then
                  fraction=(-B-sqrt(disc))/A
               else
                  fraction=0.
               endif
               if(fraction.lt.0. .or. fraction.gt.1.)fraction=1.
               return
            else
c No intersection with circle.
               return
            endif
         else
c Searching axially. If z2 also inside, no.
            if(abs(z2).lt.0.)return
c Found z1 inside, z2 outside.
         endif
      endif
c z1,z2 inside/outside, axial search. Check for inside orthogonal circle.
      r1=0
      dp=abs(z2-z1)
      do k=1,ndims-1
         i=mod(ida+k-1,ndims)+1
         ix1=indi(i)+ixnp(i)+1
         x1=(xn(ix1)-xc(i))/rc(i)
         r1=r1+x1**2
      enddo
      if(r1.lt.1.)then
c Intersection with ends
         if(z1.lt.rc(ida).neqv.z2.lt.rc(ida))then
            fraction=(z1-rc(ida))/(z1-z2)
         elseif(z1.gt.-rc(ida).neqv.z2.gt.-rc(ida))then
            fraction=(z1+rc(ida))/(z1-z2)
         endif
      endif

      end
c****************************************************************
      subroutine pllelosect(id,ipm,ndims,indi,pp,fraction,dp)
c For mesh point indi() find the nearest intersection of the leg from
c this point to the adjacent point in dimension id, direction ipm, with
c a parallelopiped defined by the data structure pp.
c Return the fractional distance in fraction (1 for no intersection),
c and total mesh spacing in dp, so that bdy distance is fraction*dp.
c
c The parallelopiped data structure in ppcom.f consists of
c 1 pp_orig : origin xc (3=ndims) 
c 4 pp_vec : 3 (covariant) vectors v_p equal half the edges.(3x3)
c 13 pp_contra : 3 contravariant vectors v^q such that v_p.v^q =
c \delta_{pq} (3x3)
c A pp_total of 21 reals (in 3-D), of which the last 9 can be calculated
c from the first 12, but must have been set prior to the call. 
c A point is inside the pp if |Sum_i(x_i-xc_i).v^p_i|<1 for all p.
c A point is on the surface if, in addition, equality holds in at least
c one of the (6) relations. 
c [i-k refers to cartesian components, p-r to pp basis.] 
c
      include '3dcom.f'
      real objg(pp_total)
      integer id,ipm,ndims
      integer indi(ndims)
      real fraction,dp
      include 'meshcom.f'
      real s1(ndims),s2(ndims)
      real small
      parameter (small=1.e-6)

c      write(*,*)(objg(k),k=1,pp_total)
 1    fraction=1.
      do j=1,ndims
         s1(j)=0.
         s2(j)=0.
      enddo
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=xn(ix1)-objg(pp_orig+i-1)
         x2=x1
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=xn(ix2)-objg(pp_orig+i-1)
            dp=abs(x2-x1)
         endif
c x1, x2 are the coordinates with respect to the pp origin.
         do j=1,ndims
            s1(j)=s1(j)+x1*objg(pp_contra+ndims*(j-1)+i-1)
            s2(j)=s2(j)+x2*objg(pp_contra+ndims*(j-1)+i-1)
         enddo
      enddo
c Now we have the contravariant coefficients in s1,s2. 
c Determine whether outside
      is1=0
      is2=0
      do j=1,ndims
         if(abs(s1(j)).eq.1. .or. abs(s2(j)).eq.1.)then
            write(*,*)j,' s1,2=',s1(j),s2(j)
     $           ,' mesh clash. Scaling contra'
            do i=1,ndims
               objg(pp_contra+ndims*(j-1)+i-1)=
     $              objg(pp_contra+ndims*(j-1)+i-1)*(1.+small)
            enddo
            goto 1
         endif
c (On the boundary counts as outside.)
         if(abs(s1(j)).ge.1.)is1=1
         if(abs(s2(j)).ge.1.)is2=1
      enddo
c      write(*,*)is1,is2,s1,s2
      if(is1.ne.is2)then
c One end inside and one outside. Find crossing. 
c We want the crossing closest to the point that's inside.
c         write(*,*)s1,s2
         do j=1,ndims
            if(s1(j).le.-1. .and. s2(j).gt.-1.)then
               f1=abs(s1(j)+1.)/(abs(s1(j)+1.)+abs(s2(j)+1.))
            elseif(s1(j).gt.-1. .and. s2(j).le.-1.)then
               f1=abs(s1(j)+1.)/(abs(s1(j)+1.)+abs(s2(j)+1.))
            elseif((s1(j)-1.).ge.0. .and. (s2(j)-1.).lt.0.)then
               f1=abs(s1(j)-1.)/(abs(s1(j)-1.)+abs(s2(j)-1.))
            elseif((s1(j)-1.).lt.0. .and. (s2(j)-1.).ge.0.)then
               f1=abs(s1(j)-1.)/(abs(s1(j)-1.)+abs(s2(j)-1.))
            else
               f1=1.
            endif
c
            if(f1.ne.1.)then
               if(is1.eq.0)then
c x1 is inside, use minimum crossing fraction.
                  fraction=min(fraction,f1)
               else
c x1 is outside, use max crossing fraction ne 1.
                  if(fraction.eq.1.)then
                     fraction=f1
                  else
                     fraction=max(fraction,f1)
                  endif
               endif
            endif
         enddo
      endif
c Now fraction= closest face-crossing fraction (or 1.)
      end
c*******************************************************************
c8888888888888888888888888888888888888888888888888888888888888888
c The following was replaced by cijdirect.
c****************************************************************
c Routine for UPDATING cij, which is used by the general mditerate.
c The geometry: fractions etc, are assumed not to have changed.
c It might be faster just to cover the intersections that are 
c known to be present.
      subroutine cijupdate(inc,ipoint,indi,ndims,iused,cij,debyelen
     $     ,error)
c If there is an object crossing next to the point at indi(ndims) whose
c pointer is ipoint, in the dimension id (plus or minus), this routine
c adjusts the cij values for situations where the cij are variables
c changing from step to step.  At present is is assumed only C changes
c and it is updated in a smoothed manner toward the potential that
c corresponds to floating.
      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
c      integer indm1(mdims)
c Shifted to correct for face omission.
      integer indi(mdims)
      integer ndims
c Error indicator
      real error
c Not used in this routine
      integer iused(ndims)
      real cij(*)
c--------------------------------------------------
c Object-data storage.
      include 'objcom.f'
c Object information
      include '3dcom.f'
c-----------------------------------------------------
      real tiny
      parameter (tiny=1.e-15)
      integer oisor
c Probably this ought to be set up in a common. But for now:
      integer iavemax
      data iavemax/50/

c----------------------------------------
c     Return increment of 1 to mditerate.
      inc=1

c ipoint here is the offset (zero-based pointer)
c The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
      icij=icb*(ipoint+1)
c oi_cij is replaced by the pointer for this extended data.
      oisor=cij(icij)
c Don't do this update unless this object's a,b,c depends on flux
c as indicated by there being extra data 
      if(oisor.eq.0 .or. idob_cij(iextra_cij,oisor).eq.0)then
         return


      else
c         if(oisor.lt.20)         write(*,*)'cijupdate',ipoint,oisor
c     $     ,idob_cij(iextra_cij,oisor),(indi(kw),kw=1,ndims)
c Prevent divide by zero issues with debyelen
         debyehere=debyelen
         if(debyehere.lt.tiny)debyehere=tiny
         ichain=1
         dibdy=dob_cij(ibdy_cij,oisor)
         dob_cij(ibdy_cij,oisor)=0.
c Iterate over dimensions.
         do id=1,ndims
c For each direction in this dimension,
            do i=1,2
               ipm=1-2*(i-1)
               im=mod(i,2)+1
               ioad=ndata_cij*(2*(id-1)+(i-1))+1
               if(dob_cij(ioad,oisor).lt.1
     $              .and.dob_cij(ioad,oisor).ge.0.)then
c We intersected an object in this direction. Adjust Cij and B_y
c Get back coef
                  iobj=idob_cij(ioad,oisor+ichain)
                  coefoa=dob_cij(ioad+2,oisor+ichain)
c These must get the information from somewhere.
c Assume that a and b are unchanged by the variability:
                  a=obj_geom(oabc,iobj)
                  b=obj_geom(oabc+1,iobj)/debyehere
c c is the thing that depends on flux, so it's going to be different.
c-------------
c Address the flux data
                  ifobj=nf_map(iobj)
                  ijbin=idob_cij(ioad+1,oisor+ichain)
c Pull the area of this facet into here
                  iaddress=ijbin+nf_address(nf_flux,ifobj,-2)
                  area=ff_data(iaddress)
c Pull the flux for this timestep and this facet.
                  flux=ff_data(ijbin
     $                 +nf_address(nf_flux,ifobj,nf_step))
c Calculate new potential
                  cnew=-a*phiofcount(flux,area)
c Smooth over nave steps. 
                  nave=min(nf_step,iavemax)
                  cold=a*dob_cij(ioad+2,oisor)
                  c=(cold*(nave-1)+cnew)/nave
c         if(oisor.lt.20)write(*,*)'cijupdate',nf_step,iobj,ijbin
c    $                 ,flux,cold,cnew,c
c or test that it gives the same answer from the old value.
c                  c=cold
c-------------

                  if(a.eq.0.) a=tiny
c Diagonal term (denominator) difference Cd-Cij stored
c This does not change
c                     dob_cij(idgs_cij,oisor)=dob_cij(idgs_cij
c     $                    ,oisor)+ coef
c Adjust potential sum (B_y or tau in new reference)
                  dob_cij(ibdy_cij,oisor)=dob_cij(ibdy_cij
     $                    ,oisor)- coefoa*c
c Now we need to update coa. Not boa or fraction 
c since they haven't changed.                  
                  if(c .ne. a*dob_cij(ioad+2,oisor))then
c                     write(*,*)'Updating coa from',dob_cij(ioad+2,oisor)
c     $                    ,' to',c/a
                     dob_cij(ioad+2,oisor)=c/a
                  endif
               endif
            enddo
c End of dimension iteration. cij coefficients now set.
         enddo
c         if(dibdy.ne.dob_cij(ibdy_cij,oisor))then
c            write(*,*)'Boundary updated from',dibdy,' to'
c     $           ,dob_cij(ibdy_cij,oisor)
c         endif
      endif

      end
