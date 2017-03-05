c Common data for vtkoutput subroutine
      integer nvtkindmax, vtkindex, vtkflag
      parameter (nvtkindmax=100000)
      real vtkpoints(nvtkindmax)
      real vtkflx(nvtkindmax)
      common /vtkcom/ vtkpoints,vtkflx,vtkindex,vtkflag
         
