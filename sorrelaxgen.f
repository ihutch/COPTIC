c********************************************************************
c Do a single relaxation.
      subroutine sorrelaxgen(k_sor,ndims,iLs,iuds,
     $     intcij,
     $     cij,u,q,myorig,
     $     laddu,faddu,oaddu,
     $     relax,rdelta,umin,umax)
c General number of dimensions ndims call.
c iLs is dimensions structure, cij is really (2*ndims+1,iuds(1),...)
c intcij is a copy of cij, but referred to as an integer. 
c Not yet used, because much of cijroutine would have to be changed.
c But it illustrates how to access a passed variable both as real and
c as integer if necessary.
c q charge. All IN
      integer ndims
      integer iLs(ndims+1),iuds(ndims)
      real cij(*),q(*)
      integer intcij(*)
c Potential IN/OUT
      real u(*)
c Origin pointer of this block in the whole arrays, IN.
      integer myorig 
c laddu (IN) true if there's additional (real) function faddu
      logical laddu
      external faddu
c oaddu is maximum relative weight of faddu term OUT
      real oaddu
c relax is the relaxation fraction IN
c rdelta is the maximum change in this step OUT
c umin, umax are the minimum and maximum values of u OUT
      real relax,rdelta,umin,umax

      parameter (imds=3,imds2=2*imds)
c Offset to adjacent points in stencil.
      integer iind(imds2)
      integer indi(imds),iused(imds)
      logical lfirst
      data lfirst/.true./
      save

c      write(*,*)'sorrelaxgen',k_sor,iuds,iLs

      if(lfirst)then
         if(ndims.gt.imds)then
            write(*,*)'sorrelax error: too many dimensions',ndims
            return
         endif
         do i=1,ndims
            do j=1,2
               iind(2*(i-1)+j)=(1-2*mod(j+1,2))*iLs(i)
            enddo
         enddo
c         write(*,*)'iind=',iind
         lfirst=.false.
      endif

      addu=0.
      daddu=0.
c Cell step.
      inc=2
c We exclude the outer boundary by decreasing used length and offsetting.
      ipoint=0
      do id=1,ndims
         iused(id)=iuds(id)-2
c If this degenerates a dimension to zero, there's nothing to do.
c Assuming that only the last relevant dimension degenerates.
         if(iused(id).le.0) return
         ipoint=ipoint+iLs(id)
c indi is the dimensioned address relative to the non-boundary origin.
         indi(id)=0
      enddo
c Hard-wired red-black relaxation order. k_sor odd => block is odd;
c odd meaning that (1,1,...) is active.
c Because of the exclusion of the boundary, the starting parity is
c adjusted by ndims.
      km=mod(k_sor+1+ndims,2)
      if(km.eq.0)then
c Odd start
      else
c Even start
         ipoint=ipoint+1
         indi(1)=indi(1)+1
      endif
c      write(*,*)'starting',ipoint,(indi(j),iused(j),j=1,ndims)

 103  continue
c Track parity changes.
      ica=0
c Starting dimension
      n=1
c Iteration over the multidimensional array. 
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm normally 1) of next dimension.
c Every carry of an even length: change parity.
         if(mod(iused(n),2).eq.0)ica=ica+1
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n).gt.iused(n)-1)goto 102
c Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
c Remember the highest level incremented-1. Obsolete.
c         ica=n-1
         goto 101
      elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
         n=n-1
         goto 101
      else
c We're at the base level and have succeeded in incrementing.
c Do whatever we need to and increment indi(1) and ipoint

c We build in correction of the increment here for red-black At each
c carry for an even length we switch (level-1) parity.
         if(mod(ica,2).ne.0) then
            iaj=(1-2*mod(indi(1),2))
            indi(1)=indi(1)+iaj
            ipoint=ipoint+iaj
c            write(*,'(a,8i8)')'ipoint=',ipoint,(indi(i),i=1,ndims),iaj
c We must also do any necessary carries first.
            goto 103
         endif

         if(.true.)then
c Start of treatment
c ipoint here is the offset (zero-based pointer)            
c         write(*,*)k_sor,ipoint,indi(1),indi(2),u(ipoint+1)
         dnum=q(ipoint+1)
         if(.not.dnum.lt.1.e30)then
            write(*,*)'dnum q-error',dnum,ipoint,q(ipoint+1)
         endif

         csum=0.
c         write(*,*)'cij',(cij((2*ndims+1)*ipoint+ic),ic=1,2*ndims)
c         write(*,*)'address,u',(ipoint+1+iind(ic),
c     $        u(ipoint+1+iind(ic)),ic=1,2*ndims)
         icind0=(2*ndims+1)*ipoint
         do ic=1,2*ndims
            icind=icind0+ic
            csum=csum+cij(icind)
            dnum=dnum+cij(icind)*u(ipoint+1+iind(ic))
         enddo
c ------------------------------------------------------------------
c Here is where the extra boundary data is needed/used. For a simple mesh
c with no embedded objects, it would not be necessary.
c Pointer to object data (converted to integer)
         io=int(cij(icind0+2*ndims+1))
         if(io.ne.0) then
c Adjust the denominator and numerator using external call.
            call ddn_cij(io,csum,dnum)
         endif
c ------------------------------------------------------------------
         if(laddu) then
            addu=faddu(u(ipoint+1),daddu,ipoint+myorig)
c The following gave problems for fnodensity which puts addu=0.
c            dscl=abs(addu)/max(abs(daddu),1.e-6)
            dscl=1.
c     relative weight of f term versus L term. Use max for next iteration.
c This seemed to be an error 2 July 09. Also dden was being calculated after.
            raddu=abs(daddu/csum)
            if(raddu.gt.oaddu) oaddu=raddu
            dnum=dnum-addu
            dden=csum+daddu
         else
            dden=csum
            dscl=1.
         endif
         if(dden.eq.0)then
            write(*,*)'sorelax error: ipoint,indi,dden,cij',ipoint,indi
     $           ,dden,(cij(icind0+ic),ic=1,2*ndims)
            stop
         endif
         delta=relax*(dnum-csum*u(ipoint+1))/dden
c Aug 09 This prevents excessive nonlinear steps.
         if(abs(delta).gt.3.*dscl)then
            delta=sign(3.*dscl,delta)
c            write(*,*)addu,daddu,u(ipoint+1),csum,dnum,dden,delta
         endif
         if(abs(delta).gt.abs(rdelta))rdelta=delta
         uij=u(ipoint+1)+delta
         u(ipoint+1)=uij
c Inverted tests to catch nans.
         if(.not.uij.ge.umin)umin=uij
         if(.not.uij.le.umax)umax=uij
c Corresponds to iftrue
         else
c test 
            write(*,*)'ipoint,k_sor',ipoint,k_sor
            u(ipoint+1)=k_sor
         endif
c End of treatment.

         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 103
      endif

 201  continue


      end
c************************************************************************
