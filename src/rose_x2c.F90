module rose_x2c

! ROSE-level handling of X2C/Foldy-Wouthuysen transformations.

  use rose_ao, only: basis_set_info_t
  use foldy_wouthuysen, only: construct_fw_coefficients
  use get_overlap, only: get_ao_overlap
  use rose_mo, only: mo_basis

  implicit none

  private

  public :: rose_fw_transform

contains

   subroutine rose_fw_transform(spinors,ao_basis)

      implicit none

      type(mo_basis), intent(inout)      :: spinors
      type(basis_set_info_t), intent(in) :: ao_basis
      type(basis_set_info_t)             :: ao_basis_small
      real(8), allocatable               :: SLL(:,:,:), SSS(:,:,:)
      real(8), allocatable               :: CFW(:,:,:)
      integer                            :: nao_large, nao_small, nmo, rcq

      if (.not. spinors%has_parent) &
         error stop 'rose_fw_transform: parent coefficients are not available.'
      if (.not. spinors%parent_read) &
         error stop 'rose_fw_transform: parent coefficients were not read from an external source.'
      if (.not. spinors%parent_active) &
         error stop 'rose_fw_transform: active coefficients are already transformed.'
      if (.not. ao_basis%has_small_component) &
         error stop 'rose_fw_transform: small-component AO basis is not available.'
      if (spinors%nao_large /= ao_basis%nao) &
         error stop 'rose_fw_transform: inconsistent large-component AO dimension.'
      if (spinors%nao_small /= ao_basis%nao_small) &
         error stop 'rose_fw_transform: inconsistent small-component AO dimension.'
      if (spinors%nao /= ao_basis%nao + ao_basis%nao_small) &
         error stop 'rose_fw_transform: inconsistent parent AO dimension.'
      if (spinors%type /= 4) &
         error stop 'rose_fw_transform: only quaternion spinors are supported.'
      if (.not. associated(spinors%coeff_parent_q)) &
         error stop 'rose_fw_transform: parent metadata without parent coefficients.'

      nao_large = spinors%nao_large
      nao_small = spinors%nao_small
      nmo       = spinors%nmo
      rcq       = 4

      ao_basis_small%basis_type = ao_basis%basis_type
      ao_basis_small%nao = ao_basis%nao_small
      ao_basis_small%basis_angular = ao_basis%basis_angular_small
      ao_basis_small%nshells = ao_basis%nshells_small
      ao_basis_small%gtos => ao_basis%gtos_small
      ao_basis_small%shell_indices => ao_basis%shell_indices_small

      allocate(SLL(nao_large,nao_large,rcq))
      allocate(SSS(nao_small,nao_small,rcq))
      SLL = 0.D0
      SSS = 0.D0
      call get_ao_overlap(ao_basis,ao_basis,SLL(:,:,1))
      call get_ao_overlap(ao_basis_small,ao_basis_small,SSS(:,:,1))

      call construct_fw_coefficients( &
         spinors%coeff_parent_q(1:nao_large,1:nmo,1:rcq), &
         spinors%coeff_parent_q(nao_large+1:nao_large+nao_small,1:nmo,1:rcq), &
         SLL,SSS,CFW)

      deallocate(spinors%coeff_q)
      allocate(spinors%coeff_q(nao_large,nmo,rcq))
      spinors%coeff_q = CFW
      spinors%nao = nao_large
      spinors%parent_active = .false.

      deallocate(SLL,SSS,CFW)

   end subroutine rose_fw_transform

end module rose_x2c
