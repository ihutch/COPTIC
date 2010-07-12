c Particle data common.
      integer n_partmax,npdim
      parameter (n_partmax=40000000,npdim=3)

c Number of actual active particles < n_partmax
      integer n_part
c Maximum particle slot that we must examine
      integer ioc_part
c Particle position and velocity (3D cartesian) in the order:
c (x,y,z) (vx,vy,vz) (xm,ym,zm) where xm... is the mesh position.
      real x_part(3*npdim,n_partmax)
c Particle flag(s).
      integer if_part(n_partmax)
c Particle previous cycle time step
      real dtprec(n_partmax)
c Timestep (unperturbed).
      real dt
c Iregion where particles belong. Obsolete.
      integer iregion_part
c Control of diagnostics
      logical ldiags
c Rho at infinity totalled over processes
      real rhoinf
c Number of reinjections this step
      integer nrein
c Average potential of reinjections
      real phirein
c Number of independent processors
      integer numprocs
c Number of reinjections at each step (if non-zero)
      integer ninjcomp
c Rho at infinity per processor, relevant only in setup.
      real ripernode
c Factor by which we relax the rhoinf calculation. 1 immediate, 0 never.
      real crelax
      common/particles/n_part,x_part,if_part,iregion_part,ioc_part
     $     ,dtprec,dt,ldiags,rhoinf,nrein,phirein,numprocs,ninjcomp
     $     ,ripernode,crelax

c Orbit plotting storage for tracking the first norbits orbits.
      integer nobsmax,norbits,nstepmax
      parameter (nobsmax=100)
      parameter (nstepmax=10000)
      real xorbit(nstepmax,nobsmax),yorbit(nstepmax,nobsmax),
     $     zorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      common /orbits/norbits,iorbitlen,xorbit,yorbit,zorbit
