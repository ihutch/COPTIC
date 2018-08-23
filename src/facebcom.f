! Common data for specification of potential boundary conditions on the 
! faces of the box.
!  LF whether we are using face boundary conditions.
!  LPF whether the condition is periodic in this dimension.
!  LCF whether for a non-periodic face, the coefficient C of the BC varies 
!      with position. If so, then 
!  LNPF whether there are any non-periodic faces
!  CxyzF are the coefficients such that 
!  C = C0F + Sum_id CxyzF(id,idn)*xn(ixnp(id)+indi(id)+1)
!  AF BF are the other coefficients AF phi + BF dphi/dn + C =0
!  If BF=0, then A is assumed =1. Otherwise precalculated coeffs are used
!  AmBF= (A/2-B/dn), ApBF= (A/2+B/dn) where dn is the outward mesh step. 
!  The precalculations are done in bdyfaceinit. 
!  AmBF ApBF are internal and ought not to be altered.

! Requires ndimsdecl.f
      logical LF,LNPF,LCF(2*ndims),LPF(ndims)
      real AF(2*ndims),BF(2*ndims)
      real C0F(2*ndims),CxyzF(ndims,2*ndims)
      real ApBF(2*ndims),AmBF(2*ndims)
      common /FaceBC/LF,LNPF,LCF,LPF,AF,BF,C0F,CxyzF,ApBF,AmBF
