c Particle data common.
      integer npartmax,npart,npdim,nstepmax
      parameter (npartmax=1000000,npdim=3)

c Particle position and velocity (3D cartesian) in the order:
c (x,y,z) (vx,vy,vz) (xm,ym,zm) where xm... is the mesh position.
      real x_part(3*npdim,npartmax)
c Particle flag(s).
      integer ipf(npartmax)
c Timestep (unperturbed).
      real dt
      common/particles/npart,x_part,ipf,dt

c Orbit plotting storage for tracking the first norbits orbits.
      integer nobsmax,norbits
      parameter (nobsmax=100)
      parameter (nstepmax=10000)
      real xorbit(nstepmax,nobsmax),yorbit(nstepmax,nobsmax),
     $     zorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      common /orbits/norbits,iorbitlen,xorbit,yorbit,zorbit
