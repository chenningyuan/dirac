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
!
!
! module containing KR-CI property functionality/ies.
!
! written by sknecht september 2010
!
module krci_properties

  implicit none

  public get_inner_product_sigmaprop_dot_c

  private

  real(8), parameter :: val_d1      =  1.0d0
  real(8), parameter :: val_d0      =  0.0d0
  real(8), parameter :: val_dm1     = -1.0d0

contains

!**********************************************************************
  subroutine get_inner_product_sigmaprop_dot_c(mat1,               &
                                               mat2,               &
                                               propval_real,       &
                                               propval_imag,       &
                                               sym_sigma,          &
                                               sym_state,          &
                                               sigma_fh,           &
                                               off_state_file_glb, &
                                               mz,                 &
                                               state_file_extension)   
!**********************************************************************
!
!    purpose: calculate the inner product of sigma and state vectors.
!
!----------------------------------------------------------------------
     real(8), intent(inout)        :: mat1(*)
     real(8), intent(inout)        :: mat2(*)
     real(8), intent(out)          :: propval_real
     real(8), intent(out)          :: propval_imag
     integer, intent(in)           :: sym_sigma
     integer, intent(in)           :: sym_state
     integer, intent(in)           :: sigma_fh
     integer, intent(inout)        :: off_state_file_glb
     integer, intent(in)           :: mz
     character (len=3), intent(in) :: state_file_extension
!----------------------------------------------------------------------
     real(8)                :: propval_real_tmp
     real(8)                :: propval_imag_tmp
     integer                :: state_fh
     character (len=14)     :: state_file_name
!**********************************************************************
      
                 propval_real_tmp = val_d0
       if(mz==2) propval_imag_tmp = val_d0

       if(sym_sigma == sym_state)then

!        step 1: open state vector file
!        ------------------------------------------------------------------------------
         state_fh = 15
         write(state_file_name,'(a11,a3)') 'KRCI_CVECS.',state_file_extension
         open(unit=state_fh,file=state_file_name,status='old',action="readwrite",     &
              form='unformatted',position='rewind')

!        step 2: obtain pct matrix in DIRAC-sorted SA-AO basis
!        ------------------------------------------------------------------------------
         off_state_file_glb = off_state_file_glb + 1

         call inprdd_real_cplx(mat1,               &
                               mat2,               &
                               sigma_fh,           &
                               state_fh,           &
                               1,                  &
                               off_state_file_glb, &
                               1,                  &
                               -1,                 &
                               propval_real_tmp,   &
                               propval_imag_tmp,   &
                               mz)

         close(unit=state_fh,status='keep')
       end if

                 propval_real = propval_real_tmp
       if(mz==2) propval_imag = propval_imag_tmp

  end subroutine get_inner_product_sigmaprop_dot_c
!**********************************************************************
end module krci_properties
