! Particle data common. 
! Requires ndimsdecl.f to define ndims, nspeciesmax
      integer n_partmax
      parameter (n_partmax=5000000)
      integer iflag,idtp
      parameter(iflag=3*ndims+1,idtp=3*ndims+2)
! Actual number of species (default=1, max nspeciesmax)
      integer nspecies
! Number of actual active particles < n_partmax
      integer nparta(nspeciesmax)
! Maximum particle slot that we must examine
      integer iocparta(nspeciesmax)
! Starting particle slot
      integer iicparta(nspeciesmax+1)
! Particle position and velocity (3D cartesian) in the order:
! (x,y,z) (vx,vy,vz) (xm,ym,zm) where xm... is the mesh position.
      real x_part(idtp,n_partmax)
! Ratio of the number of steps and inverse number of particles:
      integer numratioa(nspeciesmax)
! Whether or not the external distribution is separable
      integer notseparable(nspeciesmax)
! Timestep (unperturbed).
      real dt
! Control of diagnostics
      logical ldiags
! Rho at infinity totalled over processes
      real rhoinf
! Number of reinjections this step
      integer nrein
! Average potential of reinjections
      real phirein
! Number of independent processors
      integer numprocs
! Number of reinjections at each step (if non-zero)
      integer ninjcompa(nspeciesmax)
      real pinjcompa(nspeciesmax)
! Number of object passthroughs since last reset
      integer npassthrough
! Rho at infinity per processor, relevant only in setup.
      real ripernode
! Quiet start particle initialization max block number [1=>Not quiet.]
      integer nqblkmax
! Factor by which we relax the rhoinf calculation. 1 immediate, 0 never.
      real crelax,caverein,chi
! Flags for which dimensions are periodic or absorbing for particles.
! 0 open, 1 lower absorbs, 2 upper absorbs, 3 both absorb, 4 periodic
      integer ipartperiod(ndims)
! Effective face area for purposes of reinjection. Small if periodic.
      real fcarea(ndims)
! Whether not all directions of particles are periodic
      logical lnotallp
! Ibool defining region of particles.
      integer ibtotal_part
      parameter (ibtotal_part=100)
      integer ibool_part(ibtotal_part)
! Parameters of multidimensional hole initialization:
      integer hspecies
      real holepsi,holelen,holeum,holespeed,holeeta,holepow,holerad
      real holetoplen
! Coefficients etc of spatially-varying background
      integer bgnmax,bgn(ndims)
      parameter (bgnmax=3)
      real bga(bgnmax,ndims),bgmax(ndims)      
! These are equivalences for single-species code.
      integer n_part,iic_part,ioc_part,ninjcomp
      equivalence (n_part,nparta(1)),(iic_part,iicparta(1))
     $     ,(ioc_part,iocparta(1)),(ninjcomp,ninjcompa(1))

      common/particles/x_part,nspecies
     $     ,nparta,iicparta,iocparta,ninjcompa,pinjcompa,numratioa
     $     ,dt,ldiags,rhoinf,nrein,phirein,numprocs,npassthrough
     $     ,ripernode,crelax,ipartperiod,fcarea,lnotallp,ibool_part
     $     ,caverein,chi,notseparable,nqblkmax
     $     ,holepsi,holelen,holeum,holespeed,holeeta,holepow,holerad
     $     ,holetoplen,hspecies,bgn,bga,bgmax

! Orbit plotting storage for tracking the first norbits orbits.
! This nstepmax does NOT control the maximum number of steps.
! That is controlled by nf_maxsteps not nstepmax
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
