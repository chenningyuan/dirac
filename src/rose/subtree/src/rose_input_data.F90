module rose_input_data

! Module to hold all data that can be set by the user.

  use jagged_array
  use linear_algebra, only: define_algebra_type
  implicit none
  private :: lowercase
  public  :: set_atom_levels, get_corr_levels
  public  :: set_default_levels_from_line, set_atom_levels_from_line

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
  real(8), public                           :: reference_rank_tol=1.0D-6 ! drop reference virtuals whose residual after projection falls below this
  real(8), public                           :: reference_tie_tol=1.0D-4  ! relative width of the group of tied candidates in the pivoted selection
  logical, public                           :: preproject_reference=.false. ! project the reference fragment orbitals on B1 before building the IFOs
  real(8), allocatable, public              :: frag_bias(:) ! bias fragments when assigning the localized orbitals (for non-polar bonding orbitals)
  integer, public              :: n_groups = 0         ! number of recanonization groups; 0 = identity (resolved to n_fragments by rose_init_groups)
  integer, allocatable, public :: frag_to_group(:)     ! frag_to_group(i) = group index (1..n_groups) for IAO fragment i
  character(len=256), public   :: fragment_xyz_file = ''  ! labelled xyz file defining the recanonization groups ('' = not used)

! Control flags
  logical, public :: add_correlating               = .false.
  logical, public :: restricted                    = .true.
  logical, public :: spatial_orbitals              = .true.
  logical, public :: atom_fragments                = .false.
  logical, public :: test                          = .false.
  logical, public :: include_core                  = .false.

! Levels of correlating functions that are read from the correlating basis file (see read_corrfuns).
! Level 0 is the valence level, -1 the shell below valence, etc. The default applies to all atoms,
! atoms can be given their own list (keyed by the atom name used in the molecular geometry, so for
! instance the label in an xyz file) with the .ATOMLEVELS keyword.
  integer, parameter, public :: max_corr_levels = 8
  integer, parameter, public :: max_atom_levels = 64
  integer, public            :: n_corr_levels_default = 1
  integer, public            :: corr_levels_default(max_corr_levels) = 0
  type, public :: atom_levels_t
     character(len=12) :: name = ''
     integer           :: n = 0
     integer           :: levels(max_corr_levels) = 0
  end type atom_levels_t
  type(atom_levels_t), public :: atom_levels(max_atom_levels)
  integer, public             :: n_atom_levels = 0

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
        elseif (trim(file_for_MOs) .eq. "adf") then
           print*," Reading ADF rkf files directly is no longer supported."
           print*," Convert the calculation with python_scripts/rose_kf_bridge.py"
           print*," and set .FILE_FORMAT to h5."
           stop 'the file-based ADF interface has been removed'
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
        elseif (trim(file_for_molfrag) .eq. "adf") then
           print*," Reading ADF rkf files directly is no longer supported."
           print*," Convert the calculation with python_scripts/rose_kf_bridge.py"
           print*," and set .FRAG_FORMAT to h5."
           stop 'the file-based ADF interface has been removed'
        elseif (trim(file_for_molfrag) .eq. "psi4" .or. file_for_MOs .eq. "pyscf") then
           print*," Psi4 and pyscf are now interfaced via .h5 files, please update your input"
           error stop 'old psi4 and pyscf interfaces are no longer supported'
        elseif (trim(file_for_molfrag) .eq. "dfcoef") then
           print*," Legacy format dfcoef no longer supported, please update your DIRAC installation"
           stop 'dfcoef no longer supported'
        endif

   end subroutine rose_check_input

   subroutine set_atom_levels(name,levels)
   ! Define (or redefine) the list of correlating function levels of the atom with this name
      character(len=*), intent(in) :: name
      integer, intent(in)          :: levels(:)
      integer :: i

      if (size(levels) < 1 .or. size(levels) > max_corr_levels) then
         print*,' Number of levels for atom ',trim(name),' should be between 1 and ',max_corr_levels
         stop 'invalid input'
      end if
      do i = 1, n_atom_levels
         if (lowercase(trim(atom_levels(i)%name)) == lowercase(trim(adjustl(name)))) exit
      end do
      if (i > n_atom_levels) then
         if (n_atom_levels == max_atom_levels) stop 'too many atom specific level definitions'
         n_atom_levels = n_atom_levels + 1
         i = n_atom_levels
      end if
      atom_levels(i)%name   = adjustl(name)
      atom_levels(i)%n      = size(levels)
      atom_levels(i)%levels = 0
      atom_levels(i)%levels(1:size(levels)) = levels
   end subroutine set_atom_levels

   function get_corr_levels(name) result(levels)
   ! Levels of correlating functions to be used for the atom with this name (case insensitive match),
   ! the default list is returned if the atom has not been given a list of its own.
      character(len=*), intent(in) :: name
      integer, allocatable         :: levels(:)
      integer :: i

      do i = 1, n_atom_levels
         if (lowercase(trim(atom_levels(i)%name)) == lowercase(trim(adjustl(name)))) then
            levels = atom_levels(i)%levels(1:atom_levels(i)%n)
            return
         end if
      end do
      levels = corr_levels_default(1:n_corr_levels_default)
   end function get_corr_levels


   subroutine set_default_levels_from_line(line)
   ! Input line (.CORRLEVELS) with the default levels of correlating functions, e.g. "0 -1"
      character(len=*), intent(in) :: line
      character(len=12) :: dummy
      integer :: levels(max_corr_levels), n

      call parse_levels(line,.false.,dummy,levels,n)
      corr_levels_default = 0
      corr_levels_default(1:n) = levels(1:n)
      n_corr_levels_default = n
   end subroutine set_default_levels_from_line

   subroutine set_atom_levels_from_line(line)
   ! Input line (.ATOMLEVELS) with an atom name followed by its levels, e.g. "Se 0 -1"
      character(len=*), intent(in) :: line
      character(len=12) :: name
      integer :: levels(max_corr_levels), n

      call parse_levels(line,.true.,name,levels,n)
      call set_atom_levels(name,levels(1:n))
   end subroutine set_atom_levels_from_line

   subroutine parse_levels(line, with_name, name, levels, n)
!     Read a list of integer levels from a line, separated by blanks or commas. If with_name is set
!     the first word is the name of an atom that is returned in name.
      character(len=*), intent(in)  :: line
      logical, intent(in)           :: with_name
      character(len=*), intent(out) :: name
      integer, intent(out)          :: levels(:)
      integer, intent(out)          :: n

      character(len=len(line)) :: work
      integer :: i, j, ios, first, last

      work = line
      do i = 1, len(work)
         if (work(i:i) == ',') work(i:i) = ' '
      end do
      name = ''
      n = 0
      i = 1
      do
         first = verify(work(i:), ' ')
         if (first == 0) exit
         first = i + first - 1
         last = scan(work(first:), ' ')
         if (last == 0) then
            last = len(work)
         else
            last = first + last - 2
         end if
         if (with_name .and. n == 0 .and. len_trim(name) == 0) then
            name = work(first:last)
         else
            n = n + 1
            if (n > size(levels)) then
               print*,' Too many levels in input line: ',trim(line)
               stop 'invalid input'
            end if
            read(work(first:last),*,iostat=ios) levels(n)
            if (ios /= 0) then
               print*,' Cannot read an integer level from ',work(first:last),' in input line: ',trim(line)
               stop 'invalid input'
            end if
         end if
         i = last + 1
         if (i > len(work)) exit
      end do
      if (n == 0 .or. (with_name .and. len_trim(name) == 0)) then
         print*,' Missing levels (or atom name) in input line: ',trim(line)
         stop 'invalid input'
      end if
      j = n
   end subroutine parse_levels

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

  subroutine rose_init_groups
  ! Called after n_fragments is finalised (by rose_define_fragments).
  ! When no grouping was requested, installs the identity mapping so that
  ! n_groups == n_fragments and every IAO fragment is its own recanonisation group.
  ! This keeps all downstream code behaviour-neutral when no grouping is requested.
    integer :: i
    if (n_groups == 0) then
      n_groups = n_fragments
      if (.not. allocated(frag_to_group)) then
        allocate(frag_to_group(n_fragments))
        do i = 1, n_fragments
          frag_to_group(i) = i
        end do
      end if
    end if
    ! Basic sanity checks
    if (size(frag_to_group) /= n_fragments) &
      stop 'rose_init_groups: length of frag_to_group does not match n_fragments'
    if (minval(frag_to_group) < 1 .or. maxval(frag_to_group) > n_groups) &
      stop 'rose_init_groups: frag_to_group contains out-of-range group indices'
  end subroutine rose_init_groups

end module
