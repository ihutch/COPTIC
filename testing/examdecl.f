      include '../src/ndimsdecl.f'
c Storage array spatial count size
      include '../src/griddecl.f'
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
      real zp(na_m,na_m,ndims)
      integer ifull(ndims),iuds(ndims)
c Object data
      include '../src/objcom.f'
c Mesh spacing description structure
      include '../src/meshcom.f'
c Plasma common data
      include '../src/plascom.f'
c Filenames
      character*100 partfilename
      character*100 phifilename
      character*100 denfilename
      character*100 fluxfilename
      character*100 objfilename
      character*100 diagfilename
      character*100 argument
      
      common/examcom/ifull,iuds,u,q,zp
     $     ,partfilename,phifilename,denfilename,objfilename,argument
     $     ,fluxfilename,diagfilename
