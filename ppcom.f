c The parallelopiped data structure in ppcom.f consists of
c 1 pp_orig : origin x_0 (3) 
c 4 pp_vec : 3 (covariant) vectors v_p defining the edges from the
c origin.(3x3)
c 13 pp_contra : 3 contravariant vectors v^q such that v_p.v^q =
c \delta_{pq} (3x3)
c A pp_total of 21 reals (in 3-D), of which the last 9 can be calculated
c from the first 12, but must have been set prior to the call. 
c A point is inside the pp if 0<=Sum_i(x_i-x_{0i}).v^p_i<=1 for all p.
c A point is on the surface if, in addition, equality holds in at least
c one of the (6) relations. 
c [i-k refers to cartesian components, p-r to pp basis.] 
c
      integer pp_ndims,pp_orig,pp_vec,pp_contra,pp_total
      parameter (pp_ndims=3,pp_orig=1)
      parameter (pp_vec=pp_orig+pp_ndims)
      parameter (pp_contra=pp_vec+pp_ndims*pp_ndims)
      parameter (pp_total=pp_contra+pp_ndims*pp_ndims-1)
