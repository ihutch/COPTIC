c Particle data common. 
c Requires ndimsdecl.f to define ndims, nspeciesmax
      integer n_partmax
      parameter (n_partmax=4000000)
      integer iflag,idtp
      parameter(iflag=3*ndims+1,idtp=3*ndims+2)
c Actual number of species (default=1, max nspeciesmax)
      integer nspecies
c Number of actual active particles < n_partmax
      integer nparta(nspeciesmax)
c Maximum particle slot that we must examine
      integer iocparta(nspeciesmax)
c Starting particle slot
      integer iicparta(nspeciesmax+1)
c Particle position and velocity (3D cartesian) in the order:
c (x,y,z) (vx,vy,vz) (xm,ym,zm) where xm... is the mesh position.
      real x_part(idtp,n_partmax)
c Ratio of the number of steps and inverse number of particles:
      integer numratioa(nspeciesmax)
c Whether or not the external distribution is separable
      integer notseparable(nspeciesmax)
c Timestep (unperturbed).
      real dt
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
      integer ninjcompa(nspeciesmax)
c Rho at infinity per processor, relevant only in setup.
      real ripernode
c Factor by which we relax the rhoinf calculation. 1 immediate, 0 never.
      real crelax,caverein,chi
c Flags for which dimensions are periodic or absorbing for particles.
c 0 open, 1 lower absorbs, 2 upper absorbs, 3 both absorb, 4 periodic
      integer ipartperiod(ndims)
c Effective face area for purposes of reinjection. Small if periodic.
      real fcarea(ndims)
c Whether not all directions of particles are periodic
      logical lnotallp
      integer n_part,iic_part,ioc_part,ninjcomp
      equivalence (n_part,nparta(1)),(iic_part,iicparta(1))
     $     ,(ioc_part,iocparta(1)),(ninjcomp,ninjcompa(1))
      common/particles/x_part,nspecies
     $     ,nparta,iicparta,iocparta,ninjcompa,numratioa
     $     ,dt,ldiags,rhoinf,nrein,phirein,numprocs
     $     ,ripernode,crelax,ipartperiod,fcarea,lnotallp
     $     ,caverein,chi,notseparable

c Orbit plotting storage for tracking the first norbits orbits.
c This nstepmax does NOT control the maximum number of steps.
c That is controlled by nf_maxsteps not nstepmax
      integer nobsmax,norbits,nstepmax
      parameter (nobsmax=100)
      parameter (nstepmax=5000)
      real xorbits(nstepmax,nobsmax,nspeciesmax)
      real yorbits(nstepmax,nobsmax,nspeciesmax)
      real zorbits(nstepmax,nobsmax,nspeciesmax)
      integer iorbitlens(nobsmax,nspeciesmax)
      real xorbit(nstepmax,nobsmax)
      real yorbit(nstepmax,nobsmax)
      real zorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      equivalence (xorbit,xorbits),(yorbit,yorbits),(zorbit,zorbits)
      equivalence (iorbitlen,iorbitlens)
      common /orbits/norbits,iorbitlens,xorbits,yorbits,zorbits
