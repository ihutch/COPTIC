      integer ndims_face
      parameter (ndims_face=3)
      logical LF,LCF(2*ndims_face),LPF(ndims_face)
      real AF(2*ndims_face),BF(2*ndims_face)
      real C0F(2*ndims_face),CxyzF(ndims_face,2*ndims_face)
      real ApBF(2*ndims_face),AmBF(2*ndims_face)
      common /FaceBC/LF,LCF,LPF,AF,BF,C0F,CxyzF,ApBF,AmBF
