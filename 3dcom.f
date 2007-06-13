c Common data containing the object geometric information. 
c Each object, i < 63 has: type, data(odata).
      integer ngeomobjmax,odata,ngeomobj
      parameter (ngeomobjmax=63,odata=16)
      real obj_geom(odata,ngeomobjmax)
c Index Positions of circle data:
      common /objgeomcom/ngeomobj,obj_geom

      parameter (ocenter=2,oradius=5,oabc=8)
