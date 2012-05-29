c***********************************************************************
      subroutine bdyshare(idone,ndimsdecl,ifull,iuds,cij,u,q
     $           ,iLs,idims,lperiod,icoords,iLcoords,myside,myorig
     $           ,icommcart,mycartid,myid)
c Parallelized boundary setting routine.
c Set the boundary conditions for my faces which are true boundaries.
c A face is a boundary if my block's icoords position is 
c     0 (lower) or idims-1 (upper)
c On entry, the used information is
c    ndimsdecl       the number of dimensions
c    myside(ndims)   the length of this block's sides
c    ifull(ndims)    the full array lengths in each direction.
c    iLs(ndims+1)    the structure vector of u.
c    myorig          the starting position in the u-array of this block.
c    icoords(ndims)  the block coordinates of this block.
c    idims(ndims)    the number of blocks in each dimension.
c Unused:  iLcoords,icommcart,mycartid,myid
c On exit, return the value idone [OUT]>0 if successful.
c Pass to the mditerated setting routine the dimension,u,idone.
c
c We pass the dimensions into this routine. This works because all the 
c variables in bbdydecl.f are passed as arguments. That's why we have
c some redundant arguments.
      include 'bbdydecl.f'
      integer ifull(ndimsdecl)
      external bdyshrroutine

      ndims=ndimsdecl
c      write(*,*)'Entered bdyshare   myorig,ifull,iuds'
c     $     ,',idims,icoords,lperiod,iLcoords'
c      write(*,*)myorig,ifull,myside,idims,icoords,lperiod,iLcoords,iLs

      idone=1
      do id=1,ndims
c Tell mditerarg to iterate over steps 1. I.e. every point in the block.
c But tell it the block has only length 1 in the id dimension.
         mysave=myside(id)
         myside(id)=1
         if(icoords(id).eq.0)then
c set lower face
            ipin=myorig-1
            ioffset=-iLs(id)
            idn=id
c            write(*,*)'Entering mditerarg lower',myside,icoords(id),ipin
            call mditerarg(bdyshrroutine,ndims,ifull,myside,ipin
     $           ,idn,u,idone,ioffset)
         endif
         if(icoords(id).eq.idims(id)-1)then
c set upper face
            ipin=myorig-1+(mysave-1)*iLs(id)
            ioffset=iLs(id)
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
c On exit idone=1 indicates success.
      integer inc,ipoint,ndims,indi(ndims),iused(ndims)
      integer idn
      real u(*)

      include 'meshcom.f'
      include 'facebcom.f'
      if(LF)then
c Currently this only works for face rectangular boundary conditions.
         if(.not.LPF(mod(idn-1,ndims)+1))then
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
c Non-face BCs, just say we failed, so the fallback bdyroutine is called.
         idone=0
      endif
c      write(*,*)ipoint,(indi(k),k=1,ndims),u(1+ipoint)
      end
