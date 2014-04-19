      integer ndims,ndimsmax
c The actual used number of dimensions. This eventually might be adjustable:
      parameter (ndims=3)
c The maximum number of dimensions, which defines the allocated length of
c such variables as ifull. This can NEVER be changed because many parts
c of the code access up to three values of (e.g.) ifull.
      parameter (ndimsmax=3)
