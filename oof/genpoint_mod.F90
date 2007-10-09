!
! Module with various types of pointers 
! which serves for passing pointers
!

MODULE GeneralPointers

real(kind=8), dimension(:), pointer :: p_real8 
integer, dimension(:), pointer :: p_inte
logical, dimension(:), pointer :: p_logi

real(kind=8), dimension(:,:), pointer :: p2_real8
integer, dimension(:,:), pointer :: p2_inte 
logical, dimension(:,:), pointer :: p2_logi 

! for quaternion matrixes
real(kind=8), dimension(:,:,:), pointer :: p3_real8

END MODULE GeneralPointers

