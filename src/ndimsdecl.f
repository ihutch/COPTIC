      integer ndims,ndimsmax
! The actual used number of dimensions. This eventually might be adjustable:
      parameter (ndims=3)
! The maximum number of dimensions, which defines the allocated length of
! such variables as ifull. This can NEVER be changed because many parts
! of the code access up to three values of (e.g.) ifull.
      parameter (ndimsmax=3)
! The maximum number of species
      integer nspeciesmax
      parameter (nspeciesmax=2)
