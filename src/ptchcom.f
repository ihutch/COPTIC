c Data that is used for handling point charges, Temperature gradients,
c and partial Boltzmann electrons.
c Need to have done      include 'griddecl.f'
c Before this inclusion.
c Arrays for compensating potential, charge, Te
      real uc(na_i,na_j,na_k),rhoc(na_i,na_j,na_k)
      real Tec(na_i,na_j,na_k),boltzwt(na_i,na_j,na_k)
c Single-index access version:
      real uci(na_i*na_j*na_k),rhoci(na_i*na_j*na_k)
      real Teci(na_i*na_j*na_k)
      real boltzamp,boltzsign,boltzwti(na_i*na_j*na_k)
      equivalence (uc,uci),(rhoc,rhoci),(Tec,Teci),(boltzwt,boltzwti)
c Mask defining objects that are of special point-charge type.
      integer iptch_mask
c Copy of gtt if we need variable Te. Get from plascom?
      real gtt_copy,gnt_copy
      common /ptchcom/iptch_mask,gtt_copy,gnt_copy,uc,rhoc,Tec
     $     ,boltzamp,boltzsign,boltzwt
