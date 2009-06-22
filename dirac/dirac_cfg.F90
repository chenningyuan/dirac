module dirac_cfg

!radovan: this module can contain info about the chosen method and hamiltonian
!         the individual methods and hamiltonians can store their configuration
!         in individual modules (like for instance dft_cfg)
!         if slaves need to access the content make sure it is broadcast to them
!         (now it's not) this should be done by a method
!         which should be part of this module

  implicit none

  save

  logical, public :: dirac_cfg_dft_calculation = .false.
! logical, public :: dirac_cfg_hf_calculation  = .false.
!                    etc ...

  private

end module
