c mesh position data needed for non-uniform meshes to calculate cij
c inline x(iused(1)),y(iused(2)), etc. Specify sufficient length,
c at least the sum of the dimension lengths.
      real xn(1000)
c Pointer to the start of the dimension vector, which is xn(ixnp(id)+1)
c ixnp(id)+1 is the start of each dimension vector.
c ixnp(id+1)-ixnp(id) is the length of dimension id.
c Last element points to the last element of the last dimension, which
c is equal to the total length used of xn.
      parameter (ndims_mesh=3)
      integer ixnp(ndims_mesh+1)
      common /sormesh/ixnp,xn
