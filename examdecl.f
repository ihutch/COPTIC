
c Storage array spatial count size
      integer Li
      parameter (Li=100,ndims=3)
      real u(Li,Li,Li),q(Li,Li,Li)
      integer ifull(ndims),iuds(ndims)
      real zp(Li,Li,ndims)
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
c Logicals for control?

      
      common/examcom/ifull,iuds,u,q
     $     ,partfilename,phifilename,denfilename,objfilename,argument
