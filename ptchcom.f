c Data that is used for handling point charges.
c Need to have done      include 'griddecl.f'
c Before this inclusion.
c Array for compensating potential, charge
      real uc(na_i,na_j,na_k),rhoc(na_i,na_j,na_k)
c Single-index access version:
      real uci(na_i*na_j*na_k),rhoci(na_i*na_j*na_k)
      equivalence (uc,uci),(rhoc,rhoci)
c Copy of iptch_mask tells if we need to compensate.
      integer iptch_copy
      common /ptchcom/iptch_copy,uc,rhoc
