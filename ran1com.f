c Internal state of the fortran random number generator ran1 is made
c visible through this common block. It can then be saved by just
c writing ranstate(1:103) and restored by reading it back.
c This requires integers to be no longer than reals.
      real ranstate(103)
c Internal state of ran1:
      real rrnd(97)
      integer irx1,irx2,irx3,jrc
      equivalence (rrnd,ranstate)
      equivalence (irx1,ranstate(98)),(irx2,ranstate(99)),
     $     (irx3,ranstate(100)),(jrc,ranstate(101))
c Internal state of the gaussian random generator gasdev:
      integer gd_iset
      real gd_gset
      equivalence (gd_iset,ranstate(102)),(gd_gset,ranstate(103))
c The whole thing:
      common /ran1com/ranstate
