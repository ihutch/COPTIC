c Particle data common.
      integer n_partmax,n_part,npdim,nstepmax
      parameter (n_partmax=1000000,npdim=3)

c Particle position and velocity (3D cartesian) in the order:
c (x,y,z) (vx,vy,vz) (xm,ym,zm) where xm... is the mesh position.
      real x_part(3*npdim,n_partmax)
c Particle flag(s).
      integer if_part(n_partmax)
c Timestep (unperturbed).
      real dt
c Rho at infinity
      real rhoinf
c Iregion where particles belong.
      integer iregion_part
c Control of diagnostics
      logical ldiags
      common/particles/n_part,x_part,if_part,iregion_part,
     $     dt,ldiags,rhoinf,numprocs

c Orbit plotting storage for tracking the first norbits orbits.
      integer nobsmax,norbits
      parameter (nobsmax=100)
      parameter (nstepmax=10000)
      real xorbit(nstepmax,nobsmax),yorbit(nstepmax,nobsmax),
     $     zorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      common /orbits/norbits,iorbitlen,xorbit,yorbit,zorbit
