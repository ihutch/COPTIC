c Common data containing the object geometry information. 
c Each object, i < 31 has: type, data(odata).
      integer ngeomobjmax,odata,ngeomobj
      parameter (ngeomobjmax=31)
c Reference to the index of certain object parameters:
      integer otype,ocenter,oradius,oabc,ocylaxis,ovec,ocontra
      integer ocylrad
      parameter (otype=1,oabc=2,ocenter=5,oradius=8,ocylaxis=11)
c parallelopiped vectors start at oradius, contravariant +9
      parameter (ovec=oradius,ocontra=oradius+9,ocylrad=oradius+6)
c Flux indices.
      integer ofluxtype,ofn1,ofn2,ofn3,offc,omag
c FluxType is really the number of quantities.
      parameter (ofluxtype=ocontra+9)
c Dimensions in up to 3 directions. 
      parameter (ofn1=ofluxtype+1,ofn2=ofluxtype+2,ofn3=ofluxtype+3)
c Other references. offc->number of facets (calculated from type).
      parameter (omag=ofluxtype+1,offc=ofluxtype+4)
c Gradients of the abc coefficients.
      integer oagrad,obgrad,ocgrad
      parameter (ocgrad=offc+1,obgrad=ocgrad+3,oagrad=obgrad+3)
c Total
      parameter (odata=oagrad+2)
c The parallelopiped data structure in ppcom.f consists of
c 1 pp_orig : origin x_0 (3) which points to ocenter
c 4 pp_vec : 3 (covariant) vectors v_p equal half the edges.(3x3)
c 13 pp_contra : 3 contravariant vectors v^q such that v_p.v^q =
c \delta_{pq} (3x3)
c A pp_total of 21 reals (in 3-D), of which the last 9 can be calculated
c from the first 12, but must have been set prior to the call. 
c A point is inside the pp if |Sum_i(x_i-xc_i).v^p_i|<1 for all p.
c A point is on the surface if, in addition, equality holds in at least
c one of the (6) relations. 
c [i-k refers to cartesian components, p-r to pp basis.] 
      integer pp_orig,pp_vec,pp_contra,pp_total
      parameter (pp_orig=ocenter)
      parameter (pp_vec=pp_orig+ndims)
      parameter (pp_contra=pp_vec+ndims*ndims)
      parameter (pp_total=pp_contra+ndims*ndims-1)
      real obj_geom(odata,ngeomobjmax)
c
c Mapping from obj_geom object number to nf_flux object (many->fewer)
c Zero indicates no flux tracking for this object.
      integer nf_map(ngeomobjmax)
c Ibool defining region of particles.
      integer ibtotal_part
      parameter (ibtotal_part=100)
      integer ibool_part(ibtotal_part)
c Mask defining the objects (bits) relevant to field regions.
      integer ifield_mask
c Mask defining objects that are of special point-charge type.
      integer iptch_mask
c Has the particle region got an enclosed region
      logical lboundp
c What is the reinjection scheme?
      character*50 rjscheme
      common /objgeomcom/ngeomobj,obj_geom,nf_map
     $     ,ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
c-------------------------------------------------------------------
c Data that describes the flux to positions on the objects:
      integer nf_quant,nf_obj,nf_maxsteps,nf_datasize,nf_posdim
c Number of slots needed for position descriptors. Dimensions.
      parameter (nf_posdim=4)
c Maximum (i.e. storage size) of array 
      parameter (nf_quant=5*nspeciesmax,nf_obj=20,nf_maxsteps=6000)
      parameter (nf_datasize=10000000)
c Mnemonics for quantities:
      integer nf_flux,nf_gx,nf_gy,nf_gz,nf_heat
      parameter (nf_flux=1,nf_gx=2,nf_gy=3,nf_gz=4,nf_heat=5)
c Mnemonics for positional variables.
      integer nf_p1,nf_p2,nf_p3,nf_p4,nf_pr,nf_pt,nf_pz,nf_pa
      parameter (nf_p1=0,nf_p2=-1,nf_p3=-2,nf_p4=-3)
      parameter (nf_pr=nf_p1,nf_pt=nf_p2,nf_pz=nf_p3,nf_pa=nf_p4)
c Actual numbers of quantities, objects and steps <= maxes.
      integer nf_step,mf_quant(nf_obj),mf_obj
c Array of number of quantities by species. 
c mf_quant=sum_nspecies(if_quant) for each object.
c We don't yet have a mechanism for setting it different for different
c species
      integer kf_quant(nf_obj,nspeciesmax),if_quant(nf_obj,nspeciesmax)
c The number of positions at which this quantity is measured:
      integer nf_posno(nf_quant,nf_obj)
c The dimensional structure of these: nf_posno = prod nf_dimlens
      integer nf_dimlens(nf_quant,nf_obj,ndimsmax)
c The offset index to the start of cube faces
      integer nf_faceind(nf_quant,nf_obj,2*ndimsmax)
c Reverse mapping to the geomobj number from nf_obj number
      integer nf_geommap(nf_obj)
c The address of the data-start for the quantity, obj, step.
c Positions 1-nf_posdim:0 are used for position information.
c Position maxsteps+1 is used for averaging. maxsteps+2 ave/area.
      integer nf_address(nf_quant,nf_obj,1-nf_posdim:nf_maxsteps+2)
c The heap where the data actually lies.
      real ff_data(nf_datasize)
c The rhoinfin for each step 
      real ff_rho(nf_maxsteps)
c The total number of particles for each step:
      integer nf_npart(nf_maxsteps)
c The dt for each step
      real ff_dt(nf_maxsteps)
c The number of species. Should equal nspecies
      integer nf_species

      common /fluxdata/nf_species,nf_step,ff_rho,ff_dt,mf_quant,kf_quant
     $     ,if_quant,mf_obj,nf_posno,nf_npart,nf_dimlens,nf_faceind
     $     ,nf_geommap,nf_address,ff_data

c Flux explanation:
c For
c   nf_posno positions on the object for each of 
c   mf_quant quantities to be recorded for each of
c   mf_obj objects, for each of
c   nf_step steps
c   nf_address points to the start of data for quant, obj, step.
c     the value of nf_address(1,1,1-nf_posdim) is 1
c     i.e. the address is 1-based, not 0-based.
c   For each step, there are 
c     Stride of positions is 1, quantities nf_posno,
c     objects nf_posno*mf_obj
c ff_data is the heap of data.
c
c Steps 1-nf_posdim:0 store the quantitative position information. 
c where nf_posdim is the number of coefficients in position info.
c------------------------------------------------------------------
c Data for storing integrated field quantities such as forces.
      integer ns_nt,ns_np
c the size of the stress-calculating mesh in theta and psi directions
      parameter (ns_nt=20,ns_np=20)
      integer ns_flags(nf_obj)
      real fieldforce(ndims,nf_obj,nf_maxsteps)
      real pressforce(ndims,nf_obj,nf_maxsteps)
      real partforce(ndims,nf_obj,nf_maxsteps)
      real colnforce(ndims,nf_obj,nf_maxsteps)
      real charge_ns(nf_obj,nf_maxsteps)
      real surfobj(2*ndims,ns_nt,ns_np,nf_obj)
      common /stress/ns_flags,surfobj,fieldforce,pressforce
     $     ,partforce,colnforce,charge_ns
c------------------------------------------------------------------
c External field data (when used)
      logical lextfield
      real extfield(ndims)
      common /extfieldcom/lextfield,extfield
