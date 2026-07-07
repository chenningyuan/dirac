module foldy_wouthuysen

! Small numerical helpers for Foldy-Wouthuysen/X2C coefficient transforms.
! This module deliberately avoids dependencies on ROSE data types and file I/O.

  use linear_algebra, only: construct_fw_x_matrix, construct_fw_r_matrix, hermitian_matrix_sqrt_pair, matrix_mul

  implicit none

  private

  public :: construct_fw_coefficients

contains

   subroutine construct_fw_coefficients(CL,CS,SLL,SSS,CFW)

      ! Transform CL and CS to orthonormal large- and small-component bases
      ! before constructing X. The resulting FW coefficients are transformed
      ! back to the original large-component AO basis.

      implicit none

      real(8), intent(in)              :: CL(:,:,:), CS(:,:,:)
      real(8), intent(in)              :: SLL(:,:,:), SSS(:,:,:)
      real(8), allocatable, intent(out):: CFW(:,:,:)

      real(8), allocatable             :: X(:,:,:), R(:,:,:)
      real(8), allocatable             :: SLL_sqrt(:,:,:), SLL_invsqrt(:,:,:), SSS_sqrt(:,:,:), SSS_invsqrt(:,:,:)
      real(8), allocatable             :: CL_orth(:,:,:), CS_orth(:,:,:), CFW_orth(:,:,:)
      real(8), allocatable             :: I_L(:,:,:), I_S(:,:,:)
      integer                          :: nL, nS, npos, rcq, i

      nL   = size(CL,1)
      npos = size(CL,2)
      rcq  = size(CL,3)
      nS   = size(CS,1)

      if (size(CS,2) .ne. npos) stop 'construct_fw_coefficients: CL and CS have different number of spinors.'
      if (size(CS,3) .ne. rcq) stop 'construct_fw_coefficients: CL and CS work in different algebra.'
      if (size(SLL,1) .ne. nL .or. size(SLL,2) .ne. nL) &
         stop 'construct_fw_coefficients: SLL dimensions are inconsistent with CL.'
      if (size(SSS,1) .ne. nS .or. size(SSS,2) .ne. nS) &
         stop 'construct_fw_coefficients: SSS dimensions are inconsistent with CS.'
      if (size(SLL,3) .ne. rcq .or. size(SSS,3) .ne. rcq) &
         stop 'construct_fw_coefficients: overlaps and coefficients work in different algebra.'

      allocate(X(nS,nL,rcq))
      allocate(R(nL,nL,rcq))
      allocate(CFW(nL,npos,rcq))
      allocate(SLL_sqrt(nL,nL,rcq))
      allocate(SLL_invsqrt(nL,nL,rcq))
      allocate(SSS_sqrt(nS,nS,rcq))
      allocate(SSS_invsqrt(nS,nS,rcq))
      allocate(CL_orth(nL,npos,rcq))
      allocate(CS_orth(nS,npos,rcq))
      allocate(CFW_orth(nL,npos,rcq))
      allocate(I_L(nL,nL,rcq))
      allocate(I_S(nS,nS,rcq))

      call hermitian_matrix_sqrt_pair(SLL,SLL_sqrt,SLL_invsqrt,'construct_fw_coefficients')
      call hermitian_matrix_sqrt_pair(SSS,SSS_sqrt,SSS_invsqrt,'construct_fw_coefficients')

      call matrix_mul(SLL_sqrt,CL,CL_orth)
      call matrix_mul(SSS_sqrt,CS,CS_orth)

      I_L = 0.D0
      I_S = 0.D0
      forall(i=1:nL) I_L(i,i,1) = 1.D0
      forall(i=1:nS) I_S(i,i,1) = 1.D0

      call construct_fw_x_matrix(CL_orth,CS_orth,I_L,X)
      call construct_fw_r_matrix(X,I_L,I_S,R)
      call matrix_mul(R,CL_orth,CFW_orth)
      call matrix_mul(SLL_invsqrt,CFW_orth,CFW)

      deallocate(X,R,SLL_sqrt,SLL_invsqrt,SSS_sqrt,SSS_invsqrt,CL_orth,CS_orth,CFW_orth,I_L,I_S)

   end subroutine construct_fw_coefficients

end module foldy_wouthuysen
