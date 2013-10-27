c Data that is used for handling point charges, 
c and Temperature gradients.
c Need to have done      include 'griddecl.f'
c Before this inclusion.
c Arrays for compensating potential, charge, Te
      real uc(na_i,na_j,na_k),rhoc(na_i,na_j,na_k)
      real Tec(na_i,na_j,na_k)
c Single-index access version:
      real uci(na_i*na_j*na_k),rhoci(na_i*na_j*na_k)
      real Teci(na_i*na_j*na_k)
      equivalence (uc,uci),(rhoc,rhoci),(Tec,Teci)
c Copy of iptch_mask tells if we need to compensate.
      integer iptch_copy
c Copy of gtt if we need variable Te. Get from plascom?
      real gtt_copy
      common /ptchcom/iptch_copy,gtt_copy,uc,rhoc,Tec
