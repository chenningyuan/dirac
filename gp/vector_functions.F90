!dirac_copyright_start
!dirac_copyright_end
! Contains a module for elementairy vector operations.
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
