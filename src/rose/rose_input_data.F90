module rose_input_data

! Module to hold all data that can be set by the user.

  use jagged_array
  implicit none

! Integer-type input data
  integer, public :: expnt         = 2 ! Exponent to be used in the localization procedure
  integer, public :: charge        = 0 ! Molecular charge (to be depreciated and obtained from molecule data file)
  integer, public :: n_fragments   = 0 ! Number of fragments that is used (in case of molecular fragments)
  integer, public :: nfrag_avas    = 0 ! Number of fragments that is used for the avas option
  integer, allocatable, public :: list_frag_avas(:) ! Fragments used in avas
  integer, allocatable, public :: frag_nval(:)      ! Number of valence orbitals (valence occupied + valence virtuals) per fragment
  integer, allocatable, public :: frag_ncor(:)      ! Number of core orbitals per fragment

! AVAS-related input data
  type(JaggedArray), public :: avas_orbitals
  real(8), public :: avas_threshold=0.1

! Real-type input data
  real(8), public                           :: additional_virtuals_cutoff=2.0 ! add virtual orbitals with energies below this treshold
  real(8), public                           :: fragment_cutoff=10.0 ! define reference virtuals as those with energies below this treshold
  real(8), allocatable, public              :: frag_bias(:) ! bias fragments when assigning the localized orbitals (for non-polar bonding orbitals)

! Control flags
  !logical, public :: localize_occupied             = .true.
  !logical, public :: localize_virtuals             = .true.
  !logical, public :: add_correlating               = .false.
  logical, public :: localize_occupied             = .false.
  logical, public :: localize_virtuals             = .false.
  logical, public :: add_correlating               = .true.
  logical, public :: restricted                    = .true.
  logical, public :: spatial_orbitals              = .true.
  logical, public :: atom_fragments                = .false.
  logical, public :: test                          = .false.
  logical, public :: include_core                  = .false.

! String-type input data
  character(11)   , public           :: version ! Version of the IAO construction (Standard_2013, Simple_2013, Simple_2014)
  character(len=6), public           :: file_for_MOs='dirac',file_for_molfrag='dirac'
  character(len=6), public           :: fragfile ! fragment filenames for AVAS

  save

end module
