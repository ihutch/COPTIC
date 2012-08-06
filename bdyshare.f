c***********************************************************************
c Obsolete argument order.
c      subroutine bdyshare(idone,ndimsdecl,ifull,iuds,cij,u,q
c     $           ,iLs,idims,lperiod,icoords,iLcoords,myside,myorig
c     $           ,icommcart,mycartid,myid)
      subroutine bdyshare(iLs,ifull,iuds,u,idone,ndimsdecl,idims,
     $     icoords,iLcoords,myside,myorig,
     $     icommcart,mycartid,myid,lperiod)
c      parameter (ndimsdecl=3)
c Parallelized boundary setting routine.
c Set the boundary conditions for my faces which are true boundaries.
c A face is a boundary if my block's icoords position is 
c     0 (lower) or idims-1 (upper)
c On entry, the used information is
c    idone(2)        if ==1, set periodic boundaries, else not.
c    ndimsdecl       the number of dimensions
c    myside(ndims)   the length of this block's sides
c    ifull(ndims)    the full array lengths in each direction.
c    iLs(ndims+1)    the structure vector of u.
c    myorig          the starting position in the u-array of this block.
c    lperiod(ndims)  whether we are periodic in this dimension.
c    icoords(ndims)  the block coordinates of this block.
c    iLcoords(ndims+1) the structure vector of icoords.
c    idims(ndims)    the number of blocks  in each dimension.
c Unused:  iLcoords,icommcart,mycartid,myid
c On exit, return the value idone(1) [OUT]>0 if successful.
c Pass to the setting routine the dimension,u,idone.
c
c Passing the dimensions into this routine.  This ought to work because
c all the variables in bbdydecl.f are passed as arguments. That's why we
c have some redundant arguments. However, there is a major bug in that
c the logical array lperiod is then somehow incorrectly interpreted and
c the following arguments are misaligned or otherwise corrupted. This
c was found to happen even when it was a parameter. I could not fix this
c until lperiod was moved to the end of the argument list. Seems like a
c compiler bug but gfortran was also broken.
      include 'bbdydecl.f'
c Only the first element of idone is actually used to communicate to
c send information back to the calling routine. But we need more
c communication to bdyshrroutine.
      integer idone(3)
      external bdyshrroutine
c      real u(*)

      ndims=ndimsdecl
c      write(*,*)'In bdyshare:'
c             write(*,*)iLs,ifull,iuds,idone,ndimsdecl,idims,lperiod,
c     $           icoords,iLcoords,myside,myorig,icommcart,mycartid,myid
c      write(*,*)'Entered bdyshare   ndims,myorig,ifull,iuds'
c     $     ,',idims,icoords,lperiod,iLcoords'
c      write(*,*)ndims,myorig,ifull,myside,idims,icoords,lperiod,iLcoords
c     $     ,iLs

      idone(1)=1
      do id=1,ndims
         if(lperiod(id))then
            ioff=iLs(id)*(myside(id)-2)
         else
            ioff=iLs(id)
         endif
c Tell mditerarg to iterate over steps 1. I.e. every point in the block.
c But tell it the block has only length 1 in the id dimension.
         mysave=myside(id)
         myside(id)=1
         if(icoords(id).eq.0)then
c set lower face
            ipin=myorig-1
            ioffset=-ioff
            idn=id
c            write(*,*)'Entering mditerarg lower',myside,icoords(id),ipin
            call mditerarg(bdyshrroutine,ndims,ifull,myside,ipin
     $           ,idn,u,idone,ioffset)
         endif
         if(icoords(id).eq.idims(id)-1)then
c set upper face
            ipin=myorig-1+(mysave-1)*iLs(id)
            ioffset=ioff
            idn=id+ndims
c            write(*,*)'Entering mditerarg upper',myside,icoords(id),ipin
            call mditerarg(bdyshrroutine,ndims,ifull,myside,ipin
     $           ,idn,u,idone,ioffset)
         endif
         myside(id)=mysave
      enddo
      end
c**********************************************************************
      subroutine bdyshrroutine(inc,ipoint,indi,ndims,iused
     $     ,idn,u,idone,ioffset)
c Set the boundary value of u according to the boundary conditions.
c The position at which to set it is in indi(ndims) which is u(1+ipoint).
c idn is the face index. ioffset is the offset to the adjacent point.
c On entry idone(2)=1 indicates set periodic boundaries (default not).
      integer inc,ipoint,ndims,indi(ndims),iused(ndims)
      integer idn
      integer idone(2)
      real u(*)

      include 'meshcom.f'
      include 'facebcom.f'
      idone(1)=1
      if(LF)then
c Currently this only works for face rectangular boundary conditions.
         if(.not.LPF(mod(idn-1,ndims)+1).or.idone(2).eq.1)then
c Only if we are not on a periodic face:
            if(LCF(idn))then
c Variable C. Calculate:
               C=C0F(idn)
               do i=1,ndims
                  C=C+CxyzF(i,idn)*xn(ixnp(i)+indi(i)+1)
               enddo
            else
c Simple short cut.
               C=C0F(idn)
            endif
            if(BF(idn).eq.0)then
c Assume A=1.
               u(ipoint+1)=-C
            else
               u(ipoint+1)=-(C+AmBF(idn)*u(ipoint+1-ioffset))
     $              /ApBF(idn)
            endif
         endif
      else
c Non-face BCs. Return failure for fallback.
         idone(1)=0
      endif
c      write(*,*)'bdyshrroutine',ipoint,(indi(k),k=1,ndims),u(1+ipoint)
      end
