c Common for passing the dimensional structures. Needs to be
c set in the calling routine of the routine that references mditerarg
      include 'ndimsdecl.f'
      integer nasdims,iasfull(ndimsmax),iasuds(ndimsmax)
     $     ,iasum2(ndimsmax)
      common /addsubcom/nasdims,iasfull,iasuds,iasum2
