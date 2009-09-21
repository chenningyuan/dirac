module dirac_cfg

!radovan: this module can contain info about the chosen method and hamiltonian
!         the individual methods and hamiltonians can store their configuration
!         in individual modules (like for instance dft_cfg)
!         if slaves need to access the content make sure it is broadcast to them
!         (now it's not) this should be done by a method
!         which should be part of this module

  implicit none

  save

!                    kohn-sham dft
  logical, public :: dirac_cfg_dft_calculation = .false.
!                    hartree-fock
  logical, public :: dirac_cfg_hf_calculation  = .false.
!                    hf or ks
  logical, public :: dirac_cfg_scf_calculation = .false.
!                    ao response
  logical, public :: dirac_cfg_aoresponse      = .false.
!                    visualization module
  logical, public :: dirac_cfg_visual          = .false.
!                    relcc module
  logical, public :: dirac_cfg_relcc           = .false.

  private

end module
