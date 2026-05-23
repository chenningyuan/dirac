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
  logical, public :: localize_occupied             = .true.
  logical, public :: localize_virtuals             = .true.
  logical, public :: add_correlating               = .false.
  logical, public :: restricted                    = .true.
  logical, public :: spatial_orbitals              = .true.
  logical, public :: atom_fragments                = .false.
  logical, public :: test                          = .false.
  logical, public :: include_core                  = .false.

! String-type input data
  character(11)   , public           :: version ! Version of the IAO construction (Standard_2013, Simple_2013, Simple_2014)
#ifdef PRG_DIRAC
  character(len=6), public           :: file_for_MOs = "dirac", file_for_molfrag = "dirac" ! dirac, h5, fchk, ams, adf
#else
  character(len=6), public           :: file_for_MOs = "h5", file_for_molfrag = "h5" ! dirac, h5, fchk, ams, adf
#endif
  character(len=6), public           :: fragfile ! fragment filenames for AVAS

  save

contains

  subroutine rose_check_input

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !           CHECK IF FILE TYPE IS SUPPORTED      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (((trim(file_for_MOs) .eq. "adf") .and. (trim(file_for_molfrag) .ne. "adf")) .or. &
            ((trim(file_for_MOs) .ne. "adf") .and. (trim(file_for_molfrag) .eq. "adf"))) then
           error stop 'ADF has to be used for both molecule and fragments, because ADF works with Slater type orbitals'
        endif

        if (trim(file_for_MOs) .ne. "dirac" .and. &
            trim(file_for_MOs) .ne. "dfcoef" .and. &
            trim(file_for_MOs) .ne. "h5" .and. &
            trim(file_for_MOs) .ne. "fchk" .and. &
            trim(file_for_MOs) .ne. "psi4" .and. &
            trim(file_for_MOs) .ne. "pyscf" .and. &
            trim(file_for_MOs) .ne. "ams" .and. &
            trim(file_for_MOs) .ne. "adf") then
            print*,'Wrong input for the file containing the MOs:',trim(file_for_MOs)
            error stop
        elseif (trim(file_for_MOs) .eq. "psi4" .or. file_for_MOs .eq. "pyscf") then
           print*," Psi4 and pyscf are now interfaced via .h5 files, please update your input"
           stop 'old psi4 and pyscf interfaces are no longer supported'
        elseif (trim(file_for_MOs) .eq. "dfcoef") then
           print*," Legacy format dfcoef no longer supported, please update your DIRAC installation"
           stop 'dfcoef no longer supported'
        endif
        if (trim(file_for_molfrag) .ne. "dirac" .and. &
            trim(file_for_molfrag) .ne. "dfcoef" .and. &
            trim(file_for_molfrag) .ne. "h5" .and. &
            trim(file_for_molfrag) .ne. "fchk" .and. &
            trim(file_for_molfrag) .ne. "psi4" .and. &
            trim(file_for_molfrag) .ne. "pyscf" .and. &
            trim(file_for_molfrag) .ne. "ams" .and. &
            trim(file_for_molfrag) .ne. "adf") then
            print*,'Wrong input for the file containing the MOs of molecular fragments:',trim(file_for_molfrag)
            error stop
        elseif (trim(file_for_molfrag) .eq. "psi4" .or. file_for_MOs .eq. "pyscf") then
           print*," Psi4 and pyscf are now interfaced via .h5 files, please update your input"
           error stop 'old psi4 and pyscf interfaces are no longer supported'
        elseif (trim(file_for_molfrag) .eq. "dfcoef") then
           print*," Legacy format dfcoef no longer supported, please update your DIRAC installation"
           stop 'dfcoef no longer supported'
        endif

   end subroutine rose_check_input

end module
