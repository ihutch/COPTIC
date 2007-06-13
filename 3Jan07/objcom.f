c object-data storage. Guess at needed size. 
c current usage of first index:
c The object data consists of data enumerated as
c ndims*(2=forward/backward from point)*(2=fraction,potential)
c + diagonal + potential terms.

      integer nobj
c 3D choice needed for sharing common, I think. 
      parameter (nobj=4*3+2,Lobjmax=20000)
      real obj(nobj,lobjmax)
      integer objindex
      common /objcom/objindex,obj
