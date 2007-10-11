!
!
!  Module with various types of pointers 
! serveing for passing addresses
!
!

MODULE GeneralPointers

real(kind=8), dimension(:), pointer :: p_real8
integer,      dimension(:), pointer :: p_inte
logical,      dimension(:), pointer :: p_logi

real(kind=8), dimension(:,:), pointer :: p2_real8
integer,      dimension(:,:), pointer :: p2_inte
logical,      dimension(:,:), pointer :: p2_logi

! for quaternion matrixes of 3D arrays
real(kind=8), dimension(:,:,:), pointer :: p3_real8
integer,      dimension(:,:,:), pointer :: p3_inte
logical,      dimension(:,:,:), pointer :: p3_logi

 contains

 subroutine nullify_1Dpointers
   p_real8 => NULL()
   p_inte  => NULL()
   p_logi  => NULL()
 end subroutine nullify_1Dpointers

END MODULE GeneralPointers

