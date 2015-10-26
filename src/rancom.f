c Internal state of the random number generators.
c This requires integers to be no longer than reals.
c Actually only the internal state of the gaussian random generator gasdev:
c is needed in common. But it is convenient to put it at the end of
c the ranluxstate.
      integer gd_iset
      real gd_gset
      integer ranluxstate(27)
      equivalence(gd_iset,ranluxstate(26)),(gd_gset,ranluxstate(27))
c The whole thing:
      common /ranluxcom/ranluxstate
