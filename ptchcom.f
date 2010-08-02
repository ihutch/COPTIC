c Data representing point charges in the domains vicinity which are
c to be treated analytically. 
      integer n_ptch,n_ptchmax
      integer ndims_ptch
      parameter (n_ptchmax=20,ndims_ptch=3)
      real x_ptch(ndims_ptch,n_ptchmax)
      real q_ptch(n_ptchmax)
      common /ptch/n_ptch,x_ptch,q_ptch
