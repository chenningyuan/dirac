module dirac_rsp_utils

! general DIRAC modules
  use labeled_storage
  use checkpoint

  implicit none

  public gp_info_t
  public response_info_t
  !public response_vector_t

  private

  type gp_info_t
    private
    ! this DT is intended for data not available in CHECKPOINT, but needed in response;
    ! now, it has all data needed in response, so there is some duplication...(gosia fixme)
    integer :: nmo           = 0
    integer :: nao           = 0
    integer :: n2bbasxq      = 0
    integer :: n2orbxq       = 0
    integer :: nz            = 0
    integer :: nocc(2)       = 0
    integer :: norb(2)       = 0
    integer :: npsh(2)       = 0
    integer :: nish(2)       = 0
    integer :: nesh(2)       = 0
    integer :: nfsym         = 0
    integer :: nbsym         = 0
    integer :: ncmotq        = 0
    integer :: nfbas(2, 0:2) = 0
    integer :: i2basx(2,2)   = 0
    integer :: i2orbx(2,2)   = 0
    integer :: iorb(2)       = 0
    integer :: jbtof(0:7, 2) = 0
    integer :: ihrmop(10)    = 0
    integer :: icmoq(2)      = 0
    integer :: ihqmat(4,-1:1)= 0
    integer :: ipqtoq(4,0:7) = 0
    logical :: spinfr        = .false.
    integer :: intflg        = 0
    logical :: data_is_read  = .false.
  contains
    private
    procedure, pass, public :: erase => erase_gp_data
    procedure, pass, public :: display => display_gp_data
  end type

  type response_vector_t
    private
    ! --- variables for a single response vector --- !
    character(len=:), allocatable :: desc
    real(8), allocatable  :: values(:,:)
    integer, allocatable  :: orbital_rotation_indices(:,:)
    integer           :: id      = 0
    integer           :: lao     = 0
    integer           :: nzxope  = 0
    integer           :: nzxopp  = 0
    integer           :: nzxopeq = 0
    integer           :: nzxoppq = 0
    character(len=16) :: lab
    character(2)      :: typ
    integer           :: ih       = 0
    real(8)           :: freq     = 0.0d0
    real(8)           :: rnorm    = 0.0d0
    integer           :: fsym     = 0
    integer           :: irep     = 0
    integer           :: trsym    = 0
    integer           :: ifckop   = 0
    integer           :: jz_spfexc   = -1
    integer           :: jrep_spfexc = 0
    logical           :: rspvec_data_is_initiated = .false.
    logical           :: data_is_read = .false.
!   FIXME - data transferred from gp_info_t...
    integer :: ipqtoq(4)= 0
    integer :: ihqmat(4)= 0
    integer :: spinfr= 0
  contains
    private
    !procedure, pass, public :: init => init_rspvec_data
    procedure, pass, public :: erase => erase_rspvec_data
    procedure, pass, public :: get_size_of_values
    procedure, pass, public :: get_size_of_orbital_rotation_indices
    procedure, pass, public :: get_copy_of_values
    procedure, pass, public :: get_symmetry_info
    procedure, pass, public :: scatter_to_matrix_q
    procedure, pass, public :: get_ipqtoq
    procedure, pass, public :: get_ihqmat
    procedure, pass, public :: get_ih_sym
    procedure, pass, public :: get_tr_sym
    procedure, pass, public :: get_irep
    procedure, pass, public :: should_use_lao
    procedure, pass   :: assign_lao_order
    procedure, pass   :: transfer_gp_data_to_rsp
    procedure, pass   :: allocated_values
    procedure, pass   :: allocated_orbital_rotation_indices
    procedure, pass   :: alloc_values => alloc_rspvec_values
    procedure, pass   :: alloc_orbital_rotation_indices => alloc_rspvec_orbital_rotation_indices
    procedure, pass   :: read_values => read_rspvec_values_from_checkpoint
    procedure, pass   :: read_orbital_rotation_indices => &
                         read_rspvec_orbital_rotation_indices_from_checkpoint
    procedure, pass, public :: display => display_rspvec_data
    !procedure, pass, public :: display_orbital_rotation_indices => &
    !                           display_rspvec_orbital_rotation_indices
    !procedure, pass, public :: values_are_read
    !procedure, pass, public :: orbital_rotation_indices_are_read
  end type

  type response_info_t
    private
    ! --- general variables used in linear response --- !
    integer                     :: nr_rsp_vectors     = 0       ! total number of response vectors (as read from/written to PAMXVC file)
    integer                     :: nr_rsp_vectors_ee  = 0       ! total number of response vectors corresponding to e-e orbital rotations
    integer                     :: nr_rsp_vectors_ep  = 0       ! total number of response vectors corresponding to e-p orbital rotations
    integer                     :: ndim_ee_total      = 0       ! total space needed for all e-e orbital rotations for all response vectors
    integer                     :: ndim_ep_total      = 0       ! total space needed for all e-p orbital rotations for all response vectors
    integer                     :: nr_B_operators     = 0       ! number of B-operators
    integer                     :: ndamp              = 0       ! 1 - no damping, 2 - damped response
    integer                     :: nfreq              = 0       ! number of frequencies
    integer                     :: static             = 0       ! 0/1 - static calculations
    integer                     :: xlr_skipep         = 0       ! 0/1 - skip e-p orbital rotations
    integer                     :: xlr_skipee         = 0       ! 0/1 - skip e-e orbital rotations
    logical                     :: data_is_read= .false.
    type(response_vector_t), allocatable, public :: vectors(:)
  contains
    private
    procedure, pass         :: allocated_vectors
    procedure, pass         :: alloc => alloc_rsp_vectors
    procedure, pass, public :: erase => erase_rsp_data
    procedure, pass, public :: display => display_rsp_data
    !procedure, pass, public :: get_vectors
    procedure, pass, public :: read_from_DIRAC => read_rsp_data_from_DIRAC
    procedure, pass, public :: read_from_checkpoint => read_rsp_data_from_checkpoint
    procedure, pass, public :: write_to_checkpoint => write_rsp_data_to_checkpoint
  end type

  interface flatten
    module procedure flatten_real_2d_to_1d
    module procedure flatten_int_2d_to_1d
  end interface

contains

! ------------------------------------------------------------------------------
  pure function allocated_vectors(self)
    class(response_info_t), intent(in):: self
    logical:: allocated_vectors
    allocated_vectors = allocated(self%vectors)
  end function

! ------------------------------------------------------------------------------
  subroutine alloc_rsp_vectors(self)
    class(response_info_t), intent(inout):: self
    integer :: nr_vec
    if (self%allocated_vectors()) return
    nr_vec = self%nr_rsp_vectors
    allocate(self%vectors(nr_vec))
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_rsp_data(self)
    class(response_info_t), intent(out) :: self
  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspdim_data_from_response_file(f_unit, gp_info, rsp_info)
    integer, intent(in) :: f_unit
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    integer :: ivec, ivec_ee, ivec_ep
    logical :: vec_pp, vec_pn
    character(len=16) :: lab1

    rewind(f_unit)

    ivec    = 0
    ivec_ee = 0
    ivec_ep = 0
 1  continue
    read(f_unit) lab1
    if(lab1 .eq. 'END_OF_THIS_FILE') then
       go to 2
    else
       ivec = ivec + 1
       backspace(f_unit)
       call read_one_record_of_pamxvc(f_unit, &
                                      read_dim_only = .true., &
                                      ivec = ivec, &
                                      gp_info  = gp_info, &
                                      rsp_info = rsp_info, &
                                      rspvec_info_update = .false., &
                                      vec_pp = vec_pp, &
                                      vec_pn = vec_pn)

       if (vec_pp) then
          ivec_ee = ivec_ee + 1
       end if
       if (vec_pn) then
          ivec_ep = ivec_ep + 1
       end if

       go to 1
    end if

 2  continue

    rsp_info%nr_rsp_vectors = ivec
    rsp_info%nr_rsp_vectors_ee = ivec_ee
    rsp_info%nr_rsp_vectors_ep = ivec_ep

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspvec_from_response_file(f_unit, gp_info, rsp_info)
    integer, intent(in) :: f_unit
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    integer :: ivec, ioff, fullvec_off
    character(len=16) :: lab1

    rewind(f_unit)
    ivec = 0
    ioff = 1
    fullvec_off = 0

 1  continue
    read(f_unit) lab1
    if (lab1  .eq. 'END_OF_THIS_FILE') then
       go to 2
    else
       ivec = ivec + 1
       backspace(f_unit)
       call read_one_record_of_pamxvc(f_unit, &
                                      read_dim_only = .false., &
                                      ivec     = ivec, &
                                      gp_info  = gp_info, &
                                      rsp_info = rsp_info, &
                                      rspvec_info_update = .true., &
                                      fullvec_off = fullvec_off)
       ioff = ioff + fullvec_off
       go to 1
    end if

 2  continue

    close(f_unit, status = 'keep')
  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_one_record_of_pamxvc(f_unit,           &
                                       read_dim_only,    &
                                       ivec,             &
                                       gp_info,          &
                                       rsp_info,         &
                                       vec_pp,           &
                                       vec_pn,           &
                                       rspvec_info_update,  &
                                       fullvec_off)

    integer, intent(in) :: f_unit
    logical, intent(in) :: read_dim_only
    integer, intent(in) :: ivec
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    logical, intent(in), optional :: rspvec_info_update
    logical, intent(out), optional :: vec_pp
    logical, intent(out), optional :: vec_pn
    integer, intent(out), optional :: fullvec_off

    character(len=100) :: errmsg
    character(len=16) :: lab1, lab2, lbuf
    character(len=2) :: typ1, typ2, vector_type
    character(len=1) :: pbuf
    real(8) :: freq1, freq2, ergrsp, rnorm, ih_factor
    integer :: jsymop1, jsymop2, jtimop1, jtimop2, ivector
    integer :: nbas, norbs, nbdim, lenrec, intflg
    integer :: i, idum, length, ih, ioff
    logical :: rspvec_update
    logical :: debug_print
!
!   read one record of PAMXVC file;
!

    debug_print = .false.

    if (present(vec_pp)) vec_pp = .false.
    if (present(vec_pn)) vec_pn = .false.

    if (present(rspvec_info_update)) then
      rspvec_update = rspvec_info_update
    else
      rspvec_update = .false.
    end if

    read(f_unit) lab1, lab2, typ1, typ2, freq1, freq2,  &
                 jsymop1, jsymop2, jtimop1, jtimop2,    &
                 rnorm, lenrec, intflg, ergrsp, nbas, norbs

    length  = lenrec/(gp_info%nz)

    if (debug_print) then
      write(*, *) 'Solution vector nr. ', ivec !, &
                  !' read from ', rsp_info%rspvec_file
      write(*, '(a30,a16)')   '.... lab              = ', lab1
      write(*, '(a30,a8)')    '.... type             = ', typ1
      write(*, '(a30,i8)')    '.... pg sym (jsymop1) = ', jsymop1
      write(*, '(a30,i8)')    '.... tr sym (jtimop1) = ', jtimop1
      write(*, '(a30,f10.6)') '.... freq             = ', freq1
      write(*, '(a30,f10.6)') '.... rnorm            = ', rnorm
      write(*, '(a30,i8)')    '.... lenrec           = ', lenrec
      write(*, '(a30,i8)')    '.... length           = ', length
    end if

    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
        if (rspvec_update) then
          rsp_info%vectors(ivec)%nzxope  = length
          rsp_info%vectors(ivec)%nzxopeq = lenrec
          if (ivec.eq.1) then
            rsp_info%ndim_ee_total = 0
          else
            rsp_info%ndim_ee_total = rsp_info%ndim_ee_total + lenrec
          end if
        end if
        if (present(vec_pp)) vec_pp = .true.
      case ('P')
        vector_type = 'pn'
        if (rspvec_update) then
          rsp_info%vectors(ivec)%nzxopp  = length
          rsp_info%vectors(ivec)%nzxoppq = lenrec
          if (ivec.eq.1) then
            rsp_info%ndim_ep_total = 0
          else
            rsp_info%ndim_ep_total = rsp_info%ndim_ep_total + lenrec
          end if
        end if
        if (present(vec_pn)) vec_pn = .true.
      case default
        write (errmsg,*) 'reading PAMXVC: unexpected typ1(1:1)'
        call quit(errmsg)
    end select

    select case (typ1(2:2))
      case ('+')
        ih        =  1
        ih_factor =  1.0d0
      case ('-')
        ih        = -1
        ih_factor = -1.0d0
      case default
        write (errmsg,*) 'reading PAMXVC: unexpected typ1(2:2)'
        call quit(errmsg)
    end select

    if (rspvec_update) then
      rsp_info%vectors(ivec)%id     = ivec
      rsp_info%vectors(ivec)%lab    = lab1
      rsp_info%vectors(ivec)%typ    = typ1
      rsp_info%vectors(ivec)%freq   = freq1
      rsp_info%vectors(ivec)%rnorm  = rnorm
      rsp_info%vectors(ivec)%ih     = ih    
      rsp_info%vectors(ivec)%fsym   = jsymop1
      rsp_info%vectors(ivec)%irep   = jsymop1-1
      rsp_info%vectors(ivec)%trsym  = jtimop1
      rsp_info%vectors(ivec)%ifckop = gp_info%jbtof(jsymop1-1, 1)
      rsp_info%vectors(ivec)%data_is_read = .true.
      call rsp_info%vectors(ivec)%assign_lao_order()
      call rsp_info%vectors(ivec)%transfer_gp_data_to_rsp(&
                               gp_info%ipqtoq(1:gp_info%nz,jsymop1-1), &
                               gp_info%ihqmat(1:gp_info%nz,ih),   &
                               gp_info%spinfr)
    end if

    if (read_dim_only) then
!     we only want to count the number of vectors
      read(f_unit)
      read(f_unit)
      return
    else
!     allocate and read response vector
      call rsp_info%vectors(ivec)%alloc_values(length, gp_info%nz)
      read(f_unit) rsp_info%vectors(ivec)%values
      
!     allocate and read orbital rotation array
      call rsp_info%vectors(ivec)%alloc_orbital_rotation_indices(2, length)
      read(f_unit) rsp_info%vectors(ivec)%orbital_rotation_indices

    end if

   end subroutine


! ------------------------------------------------------------------------------
  function open_rsp_file(fname)
    character(*), intent(in) :: fname
    integer :: open_rsp_file
    integer :: funit
    logical :: exists, is_opened
    character(len=:),allocatable :: errmsg

    inquire(file=trim(fname), exist=exists, opened=is_opened)
    if (exists) then
      if (.not.is_opened) then
        open(file = trim(fname), &
             newunit = funit, &
             status = 'old',  &
             access = 'sequential',&
             form   = 'unformatted')
      else
        inquire(file=trim(fname), number=funit)
      end if
    else
      write(errmsg, '(a)') 'file ',trim(fname),' does not exist. Stopping now!'
      call quit(errmsg)
    end if
    open_rsp_file = funit
  end function


! ------------------------------------------------------------------------------
  subroutine flatten_real_2d_to_1d(arr_2d, arr_1d)
    real(8), intent(in) :: arr_2d(:,:)
    real(8), pointer, intent(out) :: arr_1d(:)
    integer :: n

    n = size(arr_2d)
    nullify(arr_1d)
    allocate(arr_1d(n))
    arr_1d(1:n) = 0.0d0
    arr_1d(1:n) = reshape(arr_2d, (/n/)) 

  end subroutine

! ------------------------------------------------------------------------------
  subroutine flatten_int_2d_to_1d(arr_2d, arr_1d)
    integer, intent(in) :: arr_2d(:,:)
    integer, pointer, intent(out) :: arr_1d(:)
    integer :: n

    n = size(arr_2d)
    nullify(arr_1d)
    allocate(arr_1d(n))
    arr_1d(1:n) = 0
    arr_1d(1:n) = reshape(arr_2d, (/n/)) 

  end subroutine


! ------------------------------------------------------------------------------
  subroutine assign_lao_order(self)
!   gosia fixme this is a temporary solution; needs to be generalized for higher-order LAO derivatives
    class(response_vector_t), intent(inout):: self
    character(len=100) :: errmsg
    if (lao_available(self)) then
      self%lao = 1
    else
      write (errmsg,*) 'WARNING: LAOs not available'
    end if
  end subroutine

! ------------------------------------------------------------------------------
  pure function lao_available(rspvec)
    class(response_vector_t), intent(in):: rspvec
    logical:: lao_available
    lao_available = (rspvec%lab(1:3) == "LAO")
  end function

! ------------------------------------------------------------------------------
  pure function should_use_lao(self)
    class(response_vector_t), intent(in):: self
    logical:: should_use_lao
    should_use_lao = (self%lao > 0 .and. lao_available(self))
  end function


! ------------------------------------------------------------------------------
  pure function allocated_values(self)
    class(response_vector_t), intent(in):: self
    logical:: allocated_values
    allocated_values = allocated(self%values)
  end function

! ------------------------------------------------------------------------------
  pure function allocated_orbital_rotation_indices(self)
    class(response_vector_t), intent(in):: self
    logical:: allocated_orbital_rotation_indices
    allocated_orbital_rotation_indices = allocated(self%orbital_rotation_indices)
  end function


! ------------------------------------------------------------------------------
  subroutine alloc_rspvec_values(self, ndim1, ndim2)
    class(response_vector_t), intent(inout):: self
    integer, intent(in) :: ndim1, ndim2
    if (self%allocated_values()) return
    allocate(self%values(ndim1, ndim2))
    self%values = 0.0d0
  end subroutine


! ------------------------------------------------------------------------------
  subroutine alloc_rspvec_orbital_rotation_indices(self, ndim1, ndim2)
    class(response_vector_t), intent(inout):: self
    integer, intent(in) :: ndim1, ndim2
    if (self%allocated_orbital_rotation_indices()) return
    allocate(self%orbital_rotation_indices(ndim1, ndim2))
    self%orbital_rotation_indices = 0
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_rspvec_data(self)
    class(response_vector_t), intent(out) :: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_gp_data(self)
    class(gp_info_t), intent(out) :: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine display_gp_data(self)
    class(gp_info_t), intent(in) :: self

    write(*,*) 'data on gp_info_t:'
    write(*,*) ' nmo        = ', self%nmo
    write(*,*) ' nao        = ', self%nao
    write(*,*) ' n2bbasxq   = ', self%n2bbasxq
    write(*,*) ' n2orbxq    = ', self%n2orbxq
    write(*,*) ' nz         = ', self%nz
    write(*,*) ' nocc       = ', self%nocc
    write(*,*) ' norb       = ', self%norb
    write(*,*) ' npsh       = ', self%npsh
    write(*,*) ' nish       = ', self%nish
    write(*,*) ' nesh       = ', self%nesh
    write(*,*) ' nfsym      = ', self%nfsym
    write(*,*) ' nbsym      = ', self%nbsym
    write(*,*) ' nfbas      = ', self%nfbas
    write(*,*) ' i2basx     = ', self%i2basx
    write(*,*) ' i2orbx     = ', self%i2orbx
    write(*,*) ' iorb       = ', self%iorb
    write(*,*) ' jbtof      = ', self%jbtof
    write(*,*) ' ihrmop     = ', self%ihrmop
    write(*,*) ' icmoq      = ', self%icmoq
    write(*,*) ' ihqmat     = ', self%ihqmat
    write(*,*) ' ncmotq     = ', self%ncmotq
    write(*,*) ' ipqtoq     = ', self%ipqtoq
    write(*,*) ' spinfr     = ', self%spinfr
    write(*,*) ' intflg     = ', self%intflg

  end subroutine


! ------------------------------------------------------------------------------
  subroutine display_rsp_data(self)
    class(response_info_t), intent(in) :: self
    integer :: i

    write(*,*) 'data on response_info_t:'
    write(*,*) ' nr_rsp_vectors    =', self%nr_rsp_vectors
    write(*,*) ' nr_rsp_vectors_ee =', self%nr_rsp_vectors_ee
    write(*,*) ' nr_rsp_vectors_ep =', self%nr_rsp_vectors_ep
    write(*,*) ' ndim_ee_total     =', self%ndim_ee_total
    write(*,*) ' ndim_ep_total     =', self%ndim_ep_total
    write(*,*) ' nr_B_operators    =', self%nr_B_operators
    write(*,*) ' ndamp             =', self%ndamp
    write(*,*) ' nfreq             =', self%nfreq
    write(*,*) ' static            =', self%static
    write(*,*) ' xlr_skipep        =', self%xlr_skipep
    write(*,*) ' xlr_skipee        =', self%xlr_skipee

    if (self%nr_rsp_vectors > 0) then
      !if (self%allocated_vectors()) then
        do i = 1, self%nr_rsp_vectors
          write(*,*) 'data on response_vector_t ', i
          call self%vectors(i)%display()
        end do
      !else
      !  write(*,*) 'response vectors NOT ALLOCATED '
      !end if
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine display_rspvec_data(self)
    class(response_vector_t), intent(in) :: self
    integer :: i, j, k, n1, n2

    if (self%data_is_read) then
      i = self%id
      write(*,*) 'data on response_vector_t ', i
      write(*,*) ' lao      =', self%lao
      write(*,*) ' nzxope   =', self%nzxope
      write(*,*) ' nzxopp   =', self%nzxopp
      write(*,*) ' nzxopeq  =', self%nzxopeq
      write(*,*) ' nzxoppq  =', self%nzxoppq
      write(*,*) ' lab      =', self%lab
      write(*,*) ' typ      =', self%typ
      write(*,*) ' ih       =', self%ih
      write(*,*) ' freq     =', self%freq
      write(*,*) ' rnorm    =', self%rnorm
      write(*,*) ' fsym     =', self%fsym
      write(*,*) ' irep     =', self%irep
      write(*,*) ' trsym    =', self%trsym
      write(*,*) ' ifckop   =', self%ifckop
      write(*,*) ' jz_spfexc   =', self%jz_spfexc
      write(*,*) ' jrep_spfexc =', self%jrep_spfexc
      if (self%allocated_values()) then
        call self%get_size_of_values(ndim1=n1, ndim2=n2)
        do j = 1, n2
          do k = 1, n1
            write(*,*) ' response_vector value =', k, j, self%values(k, j)
          end do
        end do
      else
        write(*,*) ' response_vector value NOT ALLOCATED'
      end if
      if (self%allocated_orbital_rotation_indices()) then
        call self%get_size_of_orbital_rotation_indices(ndim1=n1, ndim2=n2)
        do j = 1, n2
          do k = 1, n1
            write(*,*) ' response_vector orbital_rotation_indices =', &
                  k, j, self%orbital_rotation_indices(k, j)
          end do
        end do
      else
        write(*,*) ' response_vector orbital_rotation_indices NOT ALLOCATED'
      end if
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_gp_data_from_common_blocks(self)
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "maxaqn.h"
#include "onecom.h"
#include "ccom.h"
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbham.h"
    type(gp_info_t), intent(out) :: self

    call self%erase()

    self%nmo         = norbt
    self%nao         = ntbas(0)
    self%n2bbasxq    = n2bbasxq
    self%n2orbxq     = n2orbxq
    self%nz          = nz
    self%nocc        = nocc
    self%norb        = norb
    self%npsh        = npsh
    self%nish        = nish
    self%nesh        = nesh
    self%nfsym       = nfsym
    self%nbsym       = nbsym
    self%nfbas       = nfbas
    self%i2basx      = i2basx
    self%i2orbx      = i2orbx
    self%iorb        = iorb
    self%jbtof       = jbtof
    self%ihrmop      = ihrmop
    self%icmoq       = icmoq
    self%ihqmat      = ihqmat
    self%ncmotq      = ncmotq
    self%ipqtoq      = ipqtoq
    self%spinfr      = spinfr
    self%intflg      = intflg

    self%data_is_read = .true.

  end subroutine


! ----------------------------------------------------------------------------
  subroutine read_rsp_data_from_common_blocks(self)
#include "implicit.h"
#include "dcbxrs.h" !nfreq, ndamp, static
#include "dcbxlr.h" !nlrbpt, xlr_skipep, xlr_skipee
    type(response_info_t), intent(out) :: self
    integer :: n_ee, n_ep

    self%nfreq           = nfreq
    self%ndamp           = ndamp
    self%nr_B_operators  = nlrbpt

    if (static) then
      self%static = 1
    else
      self%static = 0
    end if

    if (xlr_skipep) then
      self%xlr_skipep = 1
    else
      self%xlr_skipep = 0
    end if

    if (xlr_skipee) then
      self%xlr_skipee = 1
    else
      self%xlr_skipee = 0
    end if
  end subroutine


! ----------------------------------------------------------------------------
  subroutine read_rsp_data_from_DIRAC(self, f_name)
!   data is read from DIRAC's common blocks and from 'PAMXVC' file
!   caution: it should be called when doing response, otherwise the
!   data on common blocks is not set!
!   this will be fixed, once all necessary data is transfered to files
    class(response_info_t), intent(out) :: self
    character(*), intent(in) :: f_name
    type(gp_info_t) :: gp_info
    integer :: f_unit

    call read_gp_data_from_common_blocks(gp_info)
    call gp_info%display()

!   read: nfreq, ndamp, static, nlrbpt, xlr_skipep, xlr_skipee
    call read_rsp_data_from_common_blocks(self)

    f_unit = open_rsp_file(f_name)
!   read: nr_rsp_vectors, nr_rsp_vectors_ee, nr_rsp_vectors_ep
    call read_rspdim_data_from_response_file(f_unit, gp_info, self)
    self%data_is_read = .true.
!   alloc+read: data related to response vectors
    call self%alloc()
    call read_rspvec_from_response_file(f_unit, gp_info, self)
    close(f_unit, status = 'keep')

    call self%display()

  end subroutine

! ------------------------------------------------------------------------------
  subroutine write_rsp_data_to_checkpoint(self, f_name)
    class(response_info_t), intent(in) :: self
    character(*), intent(in) :: f_name
    type(file_info_t) :: rsp_file

    call init_rsp_checkpoint(f_name, rsp_file)
    call write_rsp_nonallocatables_to_checkpoint(self, rsp_file)
    call write_rsp_allocatables_to_checkpoint(self, rsp_file)
    call rsp_checkpoint_close(rsp_file)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rsp_data_from_checkpoint(self, f_name)
    class(response_info_t), intent(out) :: self
    character(*), intent(in) :: f_name
    type(file_info_t) :: rsp_file

    call init_rsp_checkpoint(f_name, rsp_file, valid=.true.)
    call read_rsp_nonallocatables_from_checkpoint(self, rsp_file)
    call read_rsp_allocatables_from_checkpoint(self, rsp_file)
    call rsp_checkpoint_close(rsp_file)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rsp_nonallocatables_from_checkpoint(rsp_info, rsp_file)
    type(response_info_t), intent(inout) :: rsp_info
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80):: group

    group = '/result/response'
    call lab_read(rsp_file, trim(group)//'/nr_rsp_vectors', &
                  idata=rsp_info%nr_rsp_vectors)
    call lab_read(rsp_file, trim(group)//'/nr_rsp_vectors_ee',&
                  idata = rsp_info%nr_rsp_vectors_ee)
    call lab_read(rsp_file, trim(group)//'/nr_rsp_vectors_ep',&
                  idata = rsp_info%nr_rsp_vectors_ep)
    call lab_read(rsp_file, trim(group)//'/ndim_ee_total',&
                  idata = rsp_info%ndim_ee_total)
    call lab_read(rsp_file, trim(group)//'/ndim_ep_total',&
                  idata = rsp_info%ndim_ep_total)
    call lab_read(rsp_file, trim(group)//'/nr_B_operators',&
                  idata = rsp_info%nr_B_operators)
    call lab_read(rsp_file, trim(group)//'/ndamp',&
                  idata = rsp_info%ndamp)
    call lab_read(rsp_file, trim(group)//'/nfreq',&
                  idata = rsp_info%nfreq)
    call lab_read(rsp_file, trim(group)//'/static',&
                  idata = rsp_info%static)
    call lab_read(rsp_file, trim(group)//'/xlr_skipep',&
                  idata = rsp_info%xlr_skipep)
    call lab_read(rsp_file, trim(group)//'/xlr_skipee',&
                  idata = rsp_info%xlr_skipee)

    rsp_info%data_is_read = .true.

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rsp_allocatables_from_checkpoint(rsp_info, rsp_file)
    type(response_info_t), intent(inout) :: rsp_info
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80):: group
    character(len=4):: ich
    real(8), pointer:: tmp_rdata(:)
    integer, pointer:: tmp_idata(:)
    integer:: ivec, nr_vec

    if (.not. rsp_info%allocated_vectors()) then
      call rsp_info%alloc()
    end if

    nr_vec = rsp_info%nr_rsp_vectors
    do ivec = 1, nr_vec
      write(ich, '(i4)') ivec
      group = '/result/response/vectors/vec_'//trim(adjustl(ich))
      call lab_read(rsp_file, trim(group)//'/id', &
                   idata=rsp_info%vectors(ivec)%id)
      call lab_read(rsp_file, trim(group)//'/lao', &
                   idata=rsp_info%vectors(ivec)%lao)
      call lab_read(rsp_file, trim(group)//'/nzxope', &
                   idata=rsp_info%vectors(ivec)%nzxope)
      call lab_read(rsp_file, trim(group)//'/nzxopp', &
                   idata=rsp_info%vectors(ivec)%nzxopp)
      call lab_read(rsp_file, trim(group)//'/nzxopeq', &
                   idata=rsp_info%vectors(ivec)%nzxopeq)
      call lab_read(rsp_file, trim(group)//'/nzxoppq', &
                   idata=rsp_info%vectors(ivec)%nzxoppq)
      call lab_read(rsp_file, trim(group)//'/lab', &
                   sdata=rsp_info%vectors(ivec)%lab)
      call lab_read(rsp_file, trim(group)//'/typ', &
                   sdata=rsp_info%vectors(ivec)%typ)
      call lab_read(rsp_file, trim(group)//'/ih', &
                   idata=rsp_info%vectors(ivec)%ih)
      call lab_read(rsp_file, trim(group)//'/freq', &
                   rdata=rsp_info%vectors(ivec)%freq)
      call lab_read(rsp_file, trim(group)//'/rnorm', &
                   rdata=rsp_info%vectors(ivec)%rnorm)
      call lab_read(rsp_file, trim(group)//'/fsym', &
                   idata=rsp_info%vectors(ivec)%fsym)
      call lab_read(rsp_file, trim(group)//'/irep', &
                   idata=rsp_info%vectors(ivec)%irep)
      call lab_read(rsp_file, trim(group)//'/trsym', &
                   idata=rsp_info%vectors(ivec)%trsym)
      call lab_read(rsp_file, trim(group)//'/ifckop', &
                   idata=rsp_info%vectors(ivec)%ifckop)
      call lab_read(rsp_file, trim(group)//'/jz_spfexc', &
                   idata=rsp_info%vectors(ivec)%jz_spfexc)
      call lab_read(rsp_file, trim(group)//'/jrep_spfexc', &
                   idata=rsp_info%vectors(ivec)%jrep_spfexc)

      !FIXME
      call lab_read(rsp_file, trim(group)//'/ipqtoq', &
                    idata=rsp_info%vectors(ivec)%ipqtoq)
      call lab_read(rsp_file, trim(group)//'/ihqmat', &
                    idata=rsp_info%vectors(ivec)%ihqmat)
      call lab_read(rsp_file, trim(group)//'/spinfr', &
                    idata=rsp_info%vectors(ivec)%spinfr)

!     all non-allocatables of vector ivec are read, so alloc and read allocatables
      rsp_info%vectors(ivec)%data_is_read = .true.

      call rsp_info%vectors(ivec)%read_values(rsp_file, group)
      call rsp_info%vectors(ivec)%read_orbital_rotation_indices(rsp_file, group)
    end do

  end subroutine


! ------------------------------------------------------------------------------
  subroutine init_rsp_checkpoint(f_name, rsp_file, valid)
    character(*), intent(in) :: f_name
    type(file_info_t), intent(inout) :: rsp_file
    logical, optional, intent(in) :: valid
    character(len=80) :: top_groups(1)
    logical :: exists
    character(100) :: errmsg

    inquire(file=trim(f_name), exist=exists)
    if (.not. exists) then

      if (present(valid)) then
        if (valid) then
          write(errmsg, *) 'file ',trim(f_name),' does not exist. Stopping now!'
          call quit(errmsg)
        end if
      end if

      write(*,*) " Opening new HDF5 file for rsp data export: ", rsp_file%name
      top_groups(1) = '/result'
      rsp_file%name = trim(f_name)
      rsp_file%type = 2
      rsp_file%status = -1
      call lab_init(rsp_file, top_groups)
    else
      write(*,*) " Opening existing HDF5 file for rsp data export: ", rsp_file%name
      rsp_file%name = trim(f_name)
      rsp_file%type = 2
      rsp_file%status = 0
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine write_rsp_nonallocatables_to_checkpoint(rsp_info, rsp_file)
    type(response_info_t), intent(in) :: rsp_info
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80):: group

    group = '/result/response'
    call lab_write(rsp_file, trim(group)//'/nr_rsp_vectors', &
                   idata=rsp_info%nr_rsp_vectors)
    call lab_write(rsp_file, trim(group)//'/nr_rsp_vectors_ee',&
                   idata = rsp_info%nr_rsp_vectors_ee)
    call lab_write(rsp_file, trim(group)//'/nr_rsp_vectors_ep',&
                   idata = rsp_info%nr_rsp_vectors_ep)
    call lab_write(rsp_file, trim(group)//'/ndim_ee_total',&
                   idata = rsp_info%ndim_ee_total)
    call lab_write(rsp_file, trim(group)//'/ndim_ep_total',&
                   idata = rsp_info%ndim_ep_total)
    call lab_write(rsp_file, trim(group)//'/nr_B_operators',&
                   idata = rsp_info%nr_B_operators)
    call lab_write(rsp_file, trim(group)//'/ndamp',&
                   idata = rsp_info%ndamp)
    call lab_write(rsp_file, trim(group)//'/nfreq',&
                   idata = rsp_info%nfreq)
    call lab_write(rsp_file, trim(group)//'/static',&
                   idata = rsp_info%static)
    call lab_write(rsp_file, trim(group)//'/xlr_skipep',&
                   idata = rsp_info%xlr_skipep)
    call lab_write(rsp_file, trim(group)//'/xlr_skipee',&
                   idata = rsp_info%xlr_skipee)

  end subroutine



! ------------------------------------------------------------------------------
  subroutine write_rsp_allocatables_to_checkpoint(rsp_info, rsp_file)
    type(response_info_t), intent(in) :: rsp_info
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80):: group
    character(len=4):: ich
    real(8), pointer:: tmp_rdata(:)
    integer, pointer:: tmp_idata(:)
    integer:: ivec, nr_vec

    if (.not. rsp_info%allocated_vectors()) return

    nr_vec = rsp_info%nr_rsp_vectors
    do ivec = 1, nr_vec
      write(ich, '(i4)') ivec
      group = '/result/response/vectors/vec_'//trim(adjustl(ich))
      call lab_write(rsp_file, trim(group)//'/id', &
                   idata=rsp_info%vectors(ivec)%id)
      call lab_write(rsp_file, trim(group)//'/lao', &
                   idata=rsp_info%vectors(ivec)%lao)
      call lab_write(rsp_file, trim(group)//'/nzxope', &
                   idata=rsp_info%vectors(ivec)%nzxope)
      call lab_write(rsp_file, trim(group)//'/nzxopp', &
                   idata=rsp_info%vectors(ivec)%nzxopp)
      call lab_write(rsp_file, trim(group)//'/nzxopeq', &
                   idata=rsp_info%vectors(ivec)%nzxopeq)
      call lab_write(rsp_file, trim(group)//'/nzxoppq', &
                   idata=rsp_info%vectors(ivec)%nzxoppq)
      call lab_write(rsp_file, trim(group)//'/lab', &
                   sdata=rsp_info%vectors(ivec)%lab)
      call lab_write(rsp_file, trim(group)//'/typ', &
                   sdata=rsp_info%vectors(ivec)%typ)
      call lab_write(rsp_file, trim(group)//'/ih', &
                   idata=rsp_info%vectors(ivec)%ih)
      call lab_write(rsp_file, trim(group)//'/freq', &
                   rdata=rsp_info%vectors(ivec)%freq)
      call lab_write(rsp_file, trim(group)//'/rnorm', &
                   rdata=rsp_info%vectors(ivec)%rnorm)
      call lab_write(rsp_file, trim(group)//'/fsym', &
                   idata=rsp_info%vectors(ivec)%fsym)
      call lab_write(rsp_file, trim(group)//'/irep', &
                   idata=rsp_info%vectors(ivec)%irep)
      call lab_write(rsp_file, trim(group)//'/trsym', &
                   idata=rsp_info%vectors(ivec)%trsym)
      call lab_write(rsp_file, trim(group)//'/ifckop', &
                   idata=rsp_info%vectors(ivec)%ifckop)
      call lab_write(rsp_file, trim(group)//'/jz_spfexc', &
                   idata=rsp_info%vectors(ivec)%jz_spfexc)
      call lab_write(rsp_file, trim(group)//'/jrep_spfexc', &
                   idata=rsp_info%vectors(ivec)%jrep_spfexc)
      if (rsp_info%vectors(ivec)%allocated_values()) then
        call flatten(rsp_info%vectors(ivec)%values, tmp_rdata)
        call lab_write(rsp_file, trim(group)//'/values', &
                   rdata=tmp_rdata)
        deallocate(tmp_rdata)
      end if
      if (rsp_info%vectors(ivec)%allocated_orbital_rotation_indices()) then
        call flatten(rsp_info%vectors(ivec)%orbital_rotation_indices, tmp_idata)
        call lab_write(rsp_file, trim(group)//'/orbital_rotation_indices', &
                     idata=tmp_idata)
        deallocate(tmp_idata)
      end if
      !FIXME
      call lab_write(rsp_file, trim(group)//'/ipqtoq', &
                     idata=rsp_info%vectors(ivec)%ipqtoq)
      call lab_write(rsp_file, trim(group)//'/ihqmat', &
                     idata=rsp_info%vectors(ivec)%ihqmat)
      call lab_write(rsp_file, trim(group)//'/spinfr', &
                     idata=rsp_info%vectors(ivec)%spinfr)
    end do

  end subroutine



! ------------------------------------------------------------------------------
  function get_ihqmat(self, iq)
    class(response_vector_t), intent(in) :: self
    integer, intent(in) :: iq
    integer :: get_ihqmat

    get_ihqmat = self%ihqmat(iq)
  end function

! ------------------------------------------------------------------------------
  function get_ipqtoq(self, iz)
    class(response_vector_t), intent(in) :: self
    integer, intent(in) :: iz
    integer :: get_ipqtoq

    get_ipqtoq = self%ipqtoq(iz)
  end function


! ------------------------------------------------------------------------------
  function get_ih_sym(self)
    class(response_vector_t), intent(in) :: self
    integer :: get_ih_sym
    if (self%data_is_read) then
      get_ih_sym = self%ih
    else
      print *, 'response_vector data is not read. get_ih_sym unknown. stop'
      stop
    end if
  end function


! ------------------------------------------------------------------------------
  function get_tr_sym(self)
    class(response_vector_t), intent(in) :: self
    integer :: get_tr_sym
    if (self%data_is_read) then
      get_tr_sym = self%trsym
    else
      print *, 'response_vector data is not read. get_tr_sym unknown. stop'
      stop
    end if
  end function


! ------------------------------------------------------------------------------
  function get_irep(self)
    class(response_vector_t), intent(in) :: self
    integer :: get_irep
    if (self%data_is_read) then
      get_irep = self%irep
    else
      print *, 'response_vector data is not read. get_irep unknown. stop'
      stop
    end if
  end function

! ------------------------------------------------------------------------------
  subroutine transfer_gp_data_to_rsp(self, ipqtoq, ihqmat, spinfr)
    class(response_vector_t), intent(inout) :: self
    integer, intent(in) :: ipqtoq(4)
    integer, intent(in) :: ihqmat(4)
    logical, intent(in) :: spinfr

    self%ipqtoq = ipqtoq
    self%ihqmat = ihqmat
    if (spinfr) then
      self%spinfr = 1
    else
      self%spinfr = 0
    end if
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_symmetry_info(self, irep,ifckop,trsym,typ)
    class(response_vector_t), intent(in) :: self
    integer, intent(out) :: irep,ifckop,trsym
    character(2), intent(out) :: typ

    irep = self%irep
    ifckop = self%ifckop
    trsym = self%trsym
    typ = self%typ
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_copy_of_values(self, mat, mat_ndim)
    class(response_vector_t), intent(in) :: self
    integer, intent(in):: mat_ndim
    real(8), intent(out):: mat(mat_ndim)
    integer :: n, n1, n2, i, i1, i2

    call self%get_size_of_values(ndim=n, ndim1=n1, ndim2=n2)
    if (mat_ndim < n) stop
    mat = 0.0d0
    i = 1
    do i2 = 1, n2
      do i1 = 1, n1
        mat(i) = self%values(i1,i2)
        i = i + 1
      end do
    end do
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_size_of_values(self, ndim, ndim1, ndim2)
    class(response_vector_t), intent(in) :: self
    integer, optional, intent(out) :: ndim, ndim1, ndim2
    integer :: n, n1, n2

    if (self%typ(1:1) .eq. 'E') then
      n  = self%nzxopeq
      n1 = self%nzxope
      if (n1 /= 0) n2 = n/n1
    else if (self%typ(1:1) .eq. 'P') then
      n  = self%nzxoppq
      n1 = self%nzxopp
      if (n1 /= 0) n2 = n/n1
    end if

    if (present(ndim))  ndim  = n
    if (present(ndim1)) ndim1 = n1
    if (present(ndim2)) ndim2 = n2
  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_size_of_orbital_rotation_indices(self, ndim, ndim1, ndim2)
    class(response_vector_t), intent(in) :: self
    integer, optional, intent(out) :: ndim, ndim1, ndim2
    integer :: n, n1, n2

    if (self%typ(1:1) .eq. 'E') then
      n1 = 2
      n2 = self%nzxope
      n = n1*n2
    else if (self%typ(1:1) .eq. 'P') then
      n1 = 2
      n2 = self%nzxopp
      n = n1*n2
    end if

    if (present(ndim))  ndim  = n
    if (present(ndim1)) ndim1 = n1
    if (present(ndim2)) ndim2 = n2
  end subroutine

! ------------------------------------------------------------------------------
  subroutine read_rspvec_values_from_checkpoint(self, rsp_file, group)
    class(response_vector_t), intent(inout) :: self
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80), intent(in):: group
    integer :: ndim, ndim1, ndim2
    real(8), pointer:: tmp_rdata(:)

    call self%get_size_of_values(ndim, ndim1, ndim2)

    nullify(tmp_rdata)
    allocate(tmp_rdata(ndim))
    tmp_rdata = 0.0d0
    call lab_read(rsp_file, trim(group)//'/values', &
                   rdata=tmp_rdata)

    if (.not.self%allocated_values()) then
      call self%alloc_values(ndim1, ndim2)
    end if

    self%values = reshape(tmp_rdata, (/ndim1,ndim2/))

    deallocate(tmp_rdata)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspvec_orbital_rotation_indices_from_checkpoint(self, &
                                                       rsp_file, group)
    class(response_vector_t), intent(inout) :: self
    type(file_info_t), intent(inout) :: rsp_file
    character(len=80), intent(in):: group
    integer :: ndim, ndim1, ndim2
    integer, pointer:: tmp_idata(:)

    call self%get_size_of_orbital_rotation_indices(ndim, ndim1, ndim2)

    nullify(tmp_idata)
    allocate(tmp_idata(ndim))
    tmp_idata = 0
    call lab_read(rsp_file, trim(group)//'/orbital_rotation_indices', &
                   idata=tmp_idata)

    if (.not.self%allocated_orbital_rotation_indices()) then
      call self%alloc_orbital_rotation_indices(ndim1, ndim2)
    end if

    self%orbital_rotation_indices = reshape(tmp_idata, (/ndim1,ndim2/))

    deallocate(tmp_idata)

  end subroutine

 
! ------------------------------------------------------------------------------
  subroutine rsp_checkpoint_close(rsp_file)
    type(file_info_t), intent(inout) :: rsp_file
    rsp_file%status = 0
    print*, " RSP Checkpoint file closed"

  end subroutine


! ------------------------------------------------------------------------------
  subroutine scatter_to_matrix_q(self, mat, mat_ndim, mat_nz)
    class(response_vector_t), intent(in) :: self
    integer, intent(in) :: mat_ndim, mat_nz
    real(8), intent(out) :: mat(mat_ndim, mat_ndim, mat_nz)

    real(8)              :: ih_factor
    character(100) :: errmsg
    integer :: length, rsp_nz,   ndim, nz, irep, ih
    integer :: i, s, is, iz, iq, im
    real(8) :: f
!
!   based on scatter_vector (from old visual) and dmoli* subroutines (prp);
!   the resulting matrix is already symmetric/antisymmetric
!
   
    call self%get_size_of_values(ndim1=length,ndim2=rsp_nz)

    if (rsp_nz /= mat_nz) then
      write(errmsg, *) 'Inconsistent nz in scatter_to_matrix_q. Stop.'
      call quit(errmsg)
    end if

    do is = 1, length

      i = self%orbital_rotation_indices(1, is)
      s = self%orbital_rotation_indices(2, is)

      do iz = 1, mat_nz

        iq = self%get_ipqtoq(iz)
        im = self%get_ihqmat(iq)

        if (im.eq.1) then
          f = -1.0d0
        else if (im.eq.2) then
          f =  1.0d0
        end if

        mat(s, i, iz) =  mat(s, i, iz) &
                      + self%values(is, iz)
        mat(i, s, iz) = mat(i, s, iz) &
                         + f*self%values(is, iz)
      end do
    end do

!   this matrix is scaled by -1.0 in dmolitx1 subroutine, so we do the same here:
    mat = -1.0d0*mat

  end subroutine



end module
