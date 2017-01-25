***********************************************************************
      subroutine bdyshare(ifull,iuds,u,idone,ndimsbbdy,idims,
     $     icoords,iLcoords,myside,myorig,
     $     icommcart,mycartid,myid,lperiod)
      implicit none
      integer ndimsbbdy
      real u(*)
c      parameter (ndimsbbdy=3)
c Parallelized boundary setting routine.
c Set the boundary conditions for my faces which are true boundaries.
c A face is a boundary if my block's icoords position is 
c     0 (lower) or idims-1 (upper)
c On entry, the used information is
c    idone(2)        if ==1, set periodic boundaries, else not.
c    ndimsbbdy       the number of dimensions
c    myside(ndims)   the length of this block's sides
c    ifull(ndims)    the full array lengths in each direction.
c    myorig          the starting position in the u-array of this block.
c    lperiod(ndims)  whether we are periodic in this dimension.
c    icoords(ndims)  the block coordinates of this block.
c    iLcoords(ndims+1) the structure vector of icoords.
c    idims(ndims)    the number of blocks  in each dimension.
c Unused:  iLcoords,icommcart,mycartid,myid
c On exit, return the value idone(1) [OUT]>0 if successful.
c Pass to the setting routine the dimension,u,idone.
c
c Passing the dimensions into this routine.  This works because
c all the variables in bbdydecl.f are passed as arguments. That's why we
c have some redundant arguments.
      include 'bbdydecl.f'
c Only the first element of idone is actually used to communicate to
c send information back to the calling routine. But we need more
c communication to bdyshrroutine.
      integer idone(3)
      external bdyshrroutine
c Local variables
      integer ndims,id,ioff,ipin,mysave,ioffset,idn,ilsid

      ndims=ndimsbbdy

      if(.false.)then
         write(*,*)'myid=',myid,'icoords',icoords
      endif
      
      idone(1)=1
      ilsid=1
      do id=1,ndims
         if(lperiod(id))then
            ioff=ilsid*(myside(id)-2)
         else
            ioff=ilsid
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
c            write(*,'(a,10i6)')'Entering mditerarg lower',id,myside
c     $           ,icoords(id),ipin
            call mditerarg(bdyshrroutine,ndims,ifull,myside,ipin
     $           ,idn,u,idone,ioffset)
         endif
         if(icoords(id).eq.idims(id)-1)then
c set upper face
            ipin=myorig-1+(mysave-1)*ilsid
            ioffset=ioff
            idn=id+ndims
c            write(*,'(a,10i6)')'Entering mditerarg upper',id,myside
c     $           ,icoords(id),ipin
            call mditerarg(bdyshrroutine,ndims,ifull,myside,ipin
     $           ,idn,u,idone,ioffset)
         endif
         myside(id)=mysave
         ilsid=ilsid*ifull(id)
      enddo
      end
c**********************************************************************
      subroutine bdyshrroutine(inc,ipoint,indi,mdims,ilsunused,iused
     $     ,idn,u,idone,ioffset)
c Set the boundary value of u according to the boundary conditions.
c The position at which to set it is in indi(ndims) which is u(1+ipoint).
c idn is the face index. ioffset is the offset to the adjacent point.
c On entry idone(2)=1 indicates set periodic boundaries (default not).
      integer inc,ipoint,mdims,indi(mdims),iused(mdims)
      integer idn
      integer idone(2)
      real u(*)

      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'facebcom.f'
c Silence warnings:
      idone(1)=ilsunused
      idone(1)=inc
      idone(1)=iused(1)
      idone(1)=1
      if(LF)then
c Currently this only works for face rectangular boundary conditions.
         if(.not.LPF(mod(idn-1,ndims)+1).or.idone(2).eq.1)then
c I don't really understand that test. idone(2) seems to contradict LPF.
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
c This is where periodic conditions get set using ioffset:
               u(ipoint+1)=-(C+AmBF(idn)*u(ipoint+1-ioffset))
     $              /ApBF(idn)
            endif
         endif
      else
c Non-face BCs. Return failure for fallback.
         idone(1)=0
      endif
      if(u(ipoint+1).eq.0.)then
c         write(*,'(a,i8,4i4,2f8.3)')'bdyshrroutine',ipoint,(indi(k),k=1
c     $        ,ndims),idn,AmBF(idn),u(ipoint+1-ioffset)
      endif
      end
