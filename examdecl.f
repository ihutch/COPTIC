c Storage array spatial count size
      include 'griddecl.f'
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
      real zp(na_m,na_m,ndims)
      integer ifull(ndims),iuds(ndims)
c Object data
      include 'objcom.f'
c Mesh spacing description structure
      include 'meshcom.f'
c Plasma common data
      include 'plascom.f'
c Filenames
      character*100 partfilename
      character*100 phifilename
      character*100 denfilename
      character*100 fluxfilename
      character*100 objfilename
      character*100 argument
      
      common/examcom/ifull,iuds,u,q
     $     ,partfilename,phifilename,denfilename,objfilename,argument
     $     ,fluxfilename
