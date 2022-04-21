!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

! Contains a module for elementary 3-vector operations.
! 
! JHS, 17-04-2008
!      

module vector_fun
  implicit none      

  private
  public cross,norm

contains
  function cross(vector1,vector2)
    !     Returns the crossproduct between two vectors
    real(kind=8) :: cross(3)
    real(kind=8) :: vector1(3),vector2(3)
    cross(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    cross(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    cross(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)
  end function

  real(kind=8) function norm(vector1)
    !     Returns the norm of a vector
    real(kind=8) :: vector1(:)
    norm=sqrt(dot_product(vector1,vector1))
  end function norm

end module
