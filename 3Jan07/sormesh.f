c mesh position data needed for non-uniform meshes to calculate cij
c inline x(iused(1)),y(iused(2)), etc. Specify sufficient length,
c at least the sum of the dimension lengths.
      real xn(300)
c Pointer to the start of the dimension vector xn(ixp(id)+1)
      integer ixnp(3)
      common /sormesh/ixnp,xn
