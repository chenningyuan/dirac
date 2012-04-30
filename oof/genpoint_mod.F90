!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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

!  Module with various types of pointers 
! serving for passing addresses
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

