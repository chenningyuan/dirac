!
!
!      Copyright (c) 2008 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
!
!
!
! VERSION : $Revision$
! DATE    : $Date$
! FILE    : dirone.F
!

MODULE GRID_FILE
  IMPLICIT NONE
  PRIVATE
  public read_grid  

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
!                                 and the last is the value of the
!                                 embedding potential at that point
!                                 times the weight of the quadrature
!                                 [= v(i)*w(i)]
!               line          n : a non-negative integer, like  -42 =)
!
!

  SUBROUTINE READ_GRID(file,points,vc)
#include "priunit.h"
    REAL(kind=8),POINTER :: points(:,:)
    REAL(kind=8),pointer,optional :: vc(:)
    INTEGER                 :: file
    INTEGER                 :: i,npoints
    

    REWIND(file)
    READ(file,*) npoints
    
    ALLOCATE(points(4,1:npoints))
    
    if (present(vc)) then
       ALLOCATE(vc(1:npoints))
       DO i=1,npoints
          READ(file,*) points(:,i),vc(i)
       END DO
    else
       DO i=1,npoints
          READ(file,*) points(:,i)
       END DO
    end if
  END SUBROUTINE READ_GRID
  
  
END MODULE GRID_FILE





