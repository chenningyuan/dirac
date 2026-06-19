module rose_input_data

! Module to hold all data that can be set by the user.

  use jagged_array
  use linear_algebra, only: define_algebra_type
  implicit none
  private :: lowercase

! Integer-type input data
  integer, public :: expnt         = 2 ! Exponent to be used in the localization procedure
  integer, public :: charge        = 0 ! Molecular charge (to be depreciated and obtained from molecule data file)
  integer, public :: n_fragments   = 0 ! Number of fragments that is used (in case of molecular fragments)
  integer, public :: nfrag_avas    = 0 ! Number of fragments that is used for the avas option
  integer, public :: rcq           = 1 ! Number of real components: real, complex, or quaternion
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
  logical, public :: add_correlating               = .false.
  logical, public :: restricted                    = .true.
  logical, public :: spatial_orbitals              = .true.
  logical, public :: atom_fragments                = .false.
  logical, public :: test                          = .false.
  logical, public :: include_core                  = .false.

! String-type input data
  character(11)   , public           :: version ! Version of the IAO construction (Standard_2013, Simple_2013, Simple_2014)
  character(len=12), public           :: occupied_mode = "localized" ! original, localized
  character(len=12), public           :: virtual_mode  = "localized" ! original, localized, canonical
  character(len=12), public           :: virtual_source = "valhard" ! valence, correlating, valcorr, valhard
  character(len=1024), public         :: correlation_basis_file = ""
#ifdef PRG_DIRAC
  character(len=6), public           :: file_for_MOs = "dirac", file_for_molfrag = "dirac" ! dirac, h5, fchk, ams, adf
#else
  character(len=6), public           :: file_for_MOs = "h5", file_for_molfrag = "h5" ! dirac, h5, fchk, ams, adf
#endif
  character(len=6), public           :: fragfile ! fragment filenames for AVAS

  save

contains

  subroutine rose_check_input

      call define_algebra_type(spatial_orbitals,restricted,rcq)
      occupied_mode = lowercase(occupied_mode)
      virtual_mode = lowercase(virtual_mode)
      virtual_source = lowercase(virtual_source)
      add_correlating = trim(virtual_source) .eq. "correlating" .or. &
                        trim(virtual_source) .eq. "valcorr"

      if (trim(occupied_mode) .ne. "original" .and. trim(occupied_mode) .ne. "localized") then
         print*,'Wrong input for occupied orbital mode:',trim(occupied_mode)
         error stop
      end if

      if (trim(virtual_mode) .ne. "original" .and. trim(virtual_mode) .ne. "localized" .and. &
          trim(virtual_mode) .ne. "canonical") then
         print*,'Wrong input for virtual orbital mode:',trim(virtual_mode)
         error stop
      end if

      if (trim(virtual_source) .ne. "valence" .and. trim(virtual_source) .ne. "correlating" .and. &
          trim(virtual_source) .ne. "valcorr" .and. trim(virtual_source) .ne. "valhard") then
         print*,'Wrong input for virtual orbital source:',trim(virtual_source)
         error stop
      end if

      if (trim(virtual_mode) .eq. "original" .and. trim(virtual_source) .ne. "valhard") then
         error stop 'Original virtual orbitals cannot be combined with generated virtual sources'
      end if

      if (trim(occupied_mode) .eq. "original" .and. trim(virtual_source) .eq. "correlating" .and. &
          trim(virtual_mode) .ne. "canonical") then
         error stop 'Correlating-only runs with original occupied orbitals require canonical virtual mode'
      end if

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

   function lowercase(s)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: lowercase
      integer                      :: i, ia, iz, off

      ia = iachar('A')
      iz = iachar('Z')
      off = iachar('a') - ia
      lowercase = s
      do i = 1, len(s)
         if (iachar(s(i:i)) >= ia .and. iachar(s(i:i)) <= iz) then
            lowercase(i:i) = achar(iachar(s(i:i)) + off)
         end if
      end do
   end function lowercase

end module
