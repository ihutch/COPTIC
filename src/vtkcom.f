c Common data for vtkoutput subroutine
      integer nbmaxpoints, vtkindex, vtkflag
      parameter (nbmaxpoints=100000)
      real vtkpoints(nbmaxpoints)
      real vtkflx(nbmaxpoints)
      common /vtkcom/ vtkpoints,vtkflx,vtkindex,vtkflag
         
