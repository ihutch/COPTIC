c Common data containing the object geometry information. 
c Each object, i < 31 has: type, data(odata).
      integer ngeomobjmax,odata,ngeomobj
      parameter (ngeomobjmax=31,odata=16)
      real obj_geom(odata,ngeomobjmax)
c
c Mapping from obj_geom object number to nf_flux object (many->1)
c Zero indicates no flux tracking for this object.
      integer nf_map(ngeomobjmax)

c Ibool defining region of particles.
      integer ibtotal_part
      parameter (ibtotal_part=100)
      integer ibool_part(ibtotal_part)
      logical linregion
c Mask defining the bits relevant to field regions.
      integer ifield_mask

      common /objgeomcom/ngeomobj,obj_geom,nf_map
     $     ,ibool_part,ifield_mask

c Reference to the offset of certain object parameters:
      integer otype,ocenter,oradius,oabc
      parameter (otype=1,ocenter=2,oradius=5,oabc=8)
c
c Data that describes the flux to positions on the objects:
      integer nf_quant,nf_obj,nf_maxsteps,nf_datasize
c Number of dimensions needed for position descriptors
      parameter (nf_posdim=1)
c Maximum (i.e. storage size) of arrays. 
      parameter (nf_quant=5,nf_obj=5,nf_maxsteps=1000)
      parameter (nf_datasize=10000000)
c Mnemonics for quantities:
      parameter (nf_flux=1,nf_gx=2,nf_gy=3,nf_gz=4,nf_heat=5)
c Actual numbers of quantities, objects and steps <= maxes.
      integer nf_step,mf_quant,mf_obj
c The number of positions at which this quantity is measured:
      integer nf_posno(nf_quant,nf_obj)
c The address of the data-start for the quantity, obj, step.
      integer nf_address(nf_quant,nf_obj,1-nf_posdim:nf_maxsteps)
c The heap where the data actually lies.
      real ff_data(nf_datasize)
c The rhoinfin for each step 
      real ff_rho(nf_maxsteps)
c The dt for each step
      real ff_dt(nf_maxsteps)

      common /fluxdata/nf_step,ff_rho,ff_dt,mf_quant,mf_obj,nf_posno
     $     ,nf_address,ff_data


c Flux explanation.
c There are 
c   mf_quant quantities to be recorded for each of
c   mf_obj objects, for each of
c   nf_step steps
c   nf_address points to the start of data for quant, obj, step.
c     the value of nf_address(1,1,1-nf_posdim) is 1
c     i.e. the address is 1-based, not 0-based.
c   For each step, there are 
c     nf_posno positions on the object where quantities are recorded.
c ff_data is the heap of data.
c
c Steps 1-nf_posdim:0 stores the quantitative position information. 
c where nf_posdim is the number of dimensions in position info.
