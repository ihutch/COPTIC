c Common data containing the object geometric information. 
c Each object, i < 31 has: type, data(odata).
      integer ngeomobjmax,odata,ngeomobj
      parameter (ngeomobjmax=31,odata=16)
      real obj_geom(odata,ngeomobjmax)
c
c Mapping from obj_geom object number to nf_flux object (many->1)
c Zero indicates no flux tracking for this object.
      integer nf_map(ngeomobjmax)
      common /objgeomcom/ngeomobj,obj_geom,nf_map
c Reference to the offset of certain object parameters:
      integer ocenter,oradius,oabc
      parameter (ocenter=2,oradius=5,oabc=8)
c
c Data that describes the flux to positions on the objects:
      integer nf_quant,nf_obj,nf_maxsteps,nf_datasize
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
      integer nf_address(nf_quant,nf_obj,0:nf_maxsteps)
c The heap where the data actually lies.
      real ff_data(nf_datasize)

      common /fluxdata/nf_step,mf_quant,mf_obj,nf_posno
     $     ,nf_address,ff_data
