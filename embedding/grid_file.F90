!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

! FILE    : grid_file.F90
!

MODULE GRID_FILE
  use memory_allocator
  IMPLICIT NONE
  PRIVATE
  public read_grid  

  interface read_grid
     module procedure read_grid_onecol
     module procedure read_grid_manycol
  end interface read_grid

CONTAINS 
  
! description : this subroutine will read some operator (i.e. an embedding
!               potential) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the last the quadrature weight;
!                                 optionally, the value of the operator/property 
!                                 (the 5th column) is also read in
!

  SUBROUTINE READ_GRID_onecol(file,points,vc)
#include "priunit.h"
    REAL(kind=8),POINTER :: points(:,:)
    REAL(kind=8),pointer,optional :: vc(:)
    INTEGER                 :: file
    INTEGER                 :: i,npoints

    REWIND(file)
    READ(file,*) npoints
    
    call allocp(points,4,npoints,id="grid_points")
    
    if (present(vc)) then
       nullify(vc)
       call allocp(vc,npoints,id="potential")
       DO i=1,npoints
          READ(file,*) points(:,i),vc(i)
       END DO
    else
       DO i=1,npoints
          READ(file,*) points(:,i)
       END DO
    end if
  END SUBROUTINE READ_GRID_onecol

! description : this subroutine will read one or more operators/properties (i.e. an embedding
!               potential, a density) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point, number of properties on file
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the fourth the quadrature weight.
!                                 these are followed by whatever properties we
!                                 have in file
!
  SUBROUTINE READ_GRID_manycol(file,points,properties)
#include "priunit.h"
    REAL(kind=8),POINTER :: points(:,:)
    REAL(kind=8),pointer :: properties(:,:)
    INTEGER              :: file
    INTEGER              :: i,npoints, nprop
   

    REWIND(file)
    READ(file,*) npoints, nprop
   
    call allocp(points,4,npoints, id="grid_points")
    call allocp(properties,nprop,npoints, id="properties_on_grid")

    DO i=1,npoints
       READ(file,*) points(:,i),properties(:,i)
    END DO

  END SUBROUTINE READ_GRID_manycol
  
END MODULE GRID_FILE


