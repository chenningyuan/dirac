module dirac_rsp_utils

! general DIRAC modules
  use checkpoint

  implicit none

  public gp_info_t
  public response_info_t
!  public construct_perturbed_dmat
! temp:
  public read_gp_data
  public read_rsp_data
  public read_response_file
  public write_nonallocatables_to_helper
  public write_allocatables_to_helper
  public print_gp_info
  public print_xlr_info

  private

  type gp_info_t
    private
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
    integer :: nfbas(2, 0:2) = 0
    integer :: i2basx(2,2)   = 0
    integer :: i2orbx(2,2)   = 0
    integer :: iorb(2)       = 0
    integer :: jbtof(0:7, 2) = 0
    integer :: ihrmop(10)    = 0
    integer :: icmoq(2)      = 0
    integer :: ihqmat(4,-1:1)= 0
    integer :: ncmotq        = 0
    integer :: ipqtoq(4,0:7) = 0
    logical :: spinfr        = .false.
    integer :: intflg        = 0
    logical :: data_is_read  = .false.
  contains
    procedure, pass, public :: erase => erase_gp_data
!#ifdef VAR_MPI
!    procedure, pass, public :: bcast => bcast_gp_data
!#endif
  end type

  type response_info_t
    private
    ! --- these are general variables used in linear response --- !
    integer                     :: nr_rsp_vectors     = 0       ! total number of response vectors (as read from/written to PAMXVC file)
    integer                     :: nr_rsp_vectors_ee  = 0       ! total number of response vectors corresponding to e-e orbital rotations
    integer                     :: nr_rsp_vectors_ep  = 0       ! total number of response vectors corresponding to e-p orbital rotations
    integer                     :: ndim_ee_total      = 0       ! total space needed for all e-e orbital rotations for all response vectors
    integer                     :: ndim_ep_total      = 0       ! total space needed for all e-p orbital rotations for all response vectors
    integer                     :: nr_B_operators     = 0       ! number of B-operators
    integer                     :: ndamp              = 0       ! 1 - no damping, 2 - damped response
    integer                     :: nfreq              = 0       ! number of frequencies
    logical                     :: static             = .false. ! .true. if static calculations
    logical                     :: xlr_skipep         = .false. ! skip e-p orbital rotations
    logical                     :: xlr_skipee         = .false. ! skip e-e orbital rotations
    logical                     :: rsp_diminfo_is_read= .false.
    ! --- these are variables for each response vector --- !
    ! number of orbital rotations
    integer, pointer            :: nzxope(:)          => null()
    integer, pointer            :: nzxopp(:)          => null()
    integer, pointer            :: nzxopeq(:)         => null()
    integer, pointer            :: nzxoppq(:)         => null()
    ! variables read/written from/to PAMXVC file(s)
    ! (if used in FDE, there can be more than one PAMXVC)
    character(len=:),allocatable, public:: rspvec_file
    character(len=16), pointer  :: rspvec_lab(:)      => null()
    character(2), pointer       :: rspvec_type(:)     => null()
    integer, pointer            :: rspvec_ih(:)       => null()
    real(8), pointer            :: rspvec_freq(:)     => null()
    real(8), pointer            :: rspvec_rnorm(:)    => null()
    integer, pointer            :: rspvec_fsym(:)     => null()
    integer, pointer            :: rspvec_irep(:)     => null()
    integer, pointer            :: rspvec_trsym(:)    => null()
    integer, pointer            :: rspvec_ifckop(:)   => null()
    integer, pointer            :: jz_spfexc(:)       => null() ! -1
    integer, pointer            :: jrep_spfexc(:)     => null()
    integer, pointer            :: from_full_to_ee(:) => null()
    integer, pointer            :: from_full_to_ep(:) => null()
    logical, pointer            :: rspvec_data_is_initiated(:) => null()
    logical, pointer            :: rspvec_data_is_read(:) => null()
  contains
    procedure, pass, public :: init => init_rsp_data
    procedure, pass, public :: erase => erase_rsp_data
  end type

  interface read_gp_data
    module procedure read_gp_data_from_cbs
    module procedure read_gp_data_from_helper_file
  end interface

  interface read_response_file
    ! read response information from 'PAMXVC'
    module procedure read_rspdim_data_from_response_file
    module procedure read_rspvec_from_response_file
  end interface

  interface read_rsp_data
    ! read response information from 'helper' or common blocks
    module procedure read_rspdim_data_from_cbs
    module procedure read_rspdim_data_from_helper_file
    module procedure read_rspvec_data_from_helper_file
  end interface

contains

! ------------------------------------------------------------------------------
  subroutine init_rsp_data(self)
    class(response_info_t), intent(inout):: self
    integer :: ivec
    integer :: nr_vec

    if (.not. self%rsp_diminfo_is_read) then
      ! read from DIRAC's common blocks (non-allocatables, dimensions, ...)
      call read_rsp_data(self)
    end if

    nr_vec = self%nr_rsp_vectors

    do ivec = 1, nr_vec
      if (associated(self%rspvec_data_is_initiated)) then
        if (self%rspvec_data_is_initiated(ivec) .eqv. .true.) return
      end if
    end do

    nullify(self%nzxope)
    allocate(self%nzxope(nr_vec))
    self%nzxope(1:nr_vec) = 0

    nullify(self%nzxopeq)
    allocate(self%nzxopeq(nr_vec))
    self%nzxopeq(1:nr_vec) = 0

    nullify(self%nzxopp)
    allocate(self%nzxopp(nr_vec))
    self%nzxopp(1:nr_vec) = 0

    nullify(self%nzxoppq)
    allocate(self%nzxoppq(nr_vec))
    self%nzxoppq(1:nr_vec) = 0

    nullify(self%rspvec_lab)
    allocate(self%rspvec_lab(nr_vec))
    self%rspvec_lab(1:nr_vec) = '                '

    nullify(self%rspvec_type)
    allocate(self%rspvec_type(nr_vec))
    self%rspvec_type(1:nr_vec) = '  '

    nullify(self%rspvec_ih)
    allocate(self%rspvec_ih(nr_vec))
    self%rspvec_ih(1:nr_vec) = 0

    nullify(self%rspvec_freq)
    allocate(self%rspvec_freq(nr_vec))
    self%rspvec_freq(1:nr_vec) = 0.0d0

    nullify(self%rspvec_rnorm)
    allocate(self%rspvec_rnorm(nr_vec))
    self%rspvec_rnorm(1:nr_vec) = 0.0d0

    nullify(self%rspvec_fsym)
    allocate(self%rspvec_fsym(nr_vec))
    self%rspvec_fsym(1:nr_vec) = 0

    nullify(self%rspvec_irep)
    allocate(self%rspvec_irep(nr_vec))
    self%rspvec_irep(1:nr_vec) = 0

    nullify(self%rspvec_trsym)
    allocate(self%rspvec_trsym(nr_vec))
    self%rspvec_trsym(1:nr_vec) = 0

    nullify(self%rspvec_ifckop)
    allocate(self%rspvec_ifckop(nr_vec))
    self%rspvec_ifckop(1:nr_vec) = 0

    nullify(self%jz_spfexc)
    allocate(self%jz_spfexc(nr_vec))
    self%jz_spfexc(1:nr_vec) = 0

    nullify(self%jrep_spfexc)
    allocate(self%jrep_spfexc(nr_vec))
    self%jrep_spfexc(1:nr_vec) = 0

    nullify(self%from_full_to_ee)
    allocate(self%from_full_to_ee(nr_vec))
    self%from_full_to_ee(1:nr_vec) = 0

    nullify(self%from_full_to_ep)
    allocate(self%from_full_to_ep(nr_vec))
    self%from_full_to_ep(1:nr_vec) = 0

    nullify(self%rspvec_data_is_initiated)
    allocate(self%rspvec_data_is_initiated(nr_vec))
    self%rspvec_data_is_initiated(1:nr_vec) = .true.

    nullify(self%rspvec_data_is_read)
    allocate(self%rspvec_data_is_read(nr_vec))
    self%rspvec_data_is_read(1:nr_vec) = .false.


  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_gp_data_from_helper_file(tmp, f_name)
    type(gp_info_t), intent(inout) :: tmp
    character(*), intent(in) :: f_name
    integer :: f_unit
!
!   read from a "rsp helper" file
!
    call tmp%erase()

    open(newunit = f_unit,&
         file = trim(f_name),&
         status = 'OLD',&
         form   = 'formatted',&
         access = 'sequential',&
         position = 'rewind')

    read(f_unit, '(i8)')  tmp%nao
    read(f_unit, '(i8)')  tmp%nmo
    read(f_unit, '(i8)')  tmp%nz
    read(f_unit, '(l)')   tmp%spinfr

    close(f_unit, status = 'keep')

    tmp%data_is_read = .true.

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspdim_data_from_helper_file(tmp, f_name)
    type(response_info_t), intent(inout) :: tmp
    character(*), intent(in) :: f_name
    integer :: f_unit, idummy
    logical :: ldummy
!
!   read from a "rsp helper" file
!
    open(newunit = f_unit,&
         file = trim(f_name),&
         status = 'OLD',&
         form   = 'formatted',&
         access = 'sequential',&
         position = 'rewind')

    read(f_unit, '(i8)')  idummy
    read(f_unit, '(i8)')  idummy
    read(f_unit, '(i8)')  idummy
    read(f_unit, '(l)')   ldummy

    read(f_unit, '(i8)')  tmp%ndamp
    read(f_unit, '(i8)')  tmp%nfreq
    read(f_unit, '(i8)')  tmp%nr_B_operators
    read(f_unit, '(l)')   tmp%static
    read(f_unit, '(l)')   tmp%xlr_skipep
    read(f_unit, '(l)')   tmp%xlr_skipee
    read(f_unit, '(i8)')  tmp%nr_rsp_vectors
    read(f_unit, '(i8)')  tmp%nr_rsp_vectors_ee
    read(f_unit, '(i8)')  tmp%nr_rsp_vectors_ep

    close(f_unit, status = 'keep')

    tmp%rsp_diminfo_is_read = .true.

  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspvec_data_from_helper_file(self, f_unit, ivec)
    type(response_info_t), intent(inout) :: self
    integer, intent(in) :: f_unit
    integer, intent(inout) :: ivec
!
!   read from a "rsp helper" file
!
    read(f_unit, '(5i8,a20,a5,i8,2es28.18,8i8)')  &
                                      ivec,                      &
                                      self%nzxope(ivec),          &
                                      self%nzxopeq(ivec),         &
                                      self%nzxopp(ivec),          &
                                      self%nzxoppq(ivec),         &
                                      self%rspvec_lab(ivec),      &
                                      self%rspvec_type(ivec),     &
                                      self%rspvec_ih(ivec),       &
                                      self%rspvec_freq(ivec),     &
                                      self%rspvec_rnorm(ivec),    &
                                      self%rspvec_fsym(ivec),     &
                                      self%rspvec_irep(ivec),     &
                                      self%rspvec_trsym(ivec),    &
                                      self%rspvec_ifckop(ivec),   &
                                      self%jz_spfexc(ivec),       &
                                      self%jrep_spfexc(ivec),     &
                                      self%from_full_to_ee(ivec), &
                                      self%from_full_to_ep(ivec)

  end subroutine


! ------------------------------------------------------------------------------
  subroutine write_nonallocatables_to_helper(f_unit,gp_info,rsp_info) 
    integer, intent(in) :: f_unit
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(in) :: rsp_info

    write(f_unit, '(i8)')  gp_info%nao !on checkpoint!
    write(f_unit, '(i8)')  gp_info%nmo !on checkpoint!
    write(f_unit, '(i8)')  gp_info%nz  !on checkpoint!
    write(f_unit, '(l)')   gp_info%spinfr
  
    write(f_unit, '(i8)')  rsp_info%ndamp
    write(f_unit, '(i8)')  rsp_info%nfreq
    write(f_unit, '(i8)')  rsp_info%nr_B_operators
    write(f_unit, '(l)')   rsp_info%static
    write(f_unit, '(l)')   rsp_info%xlr_skipep
    write(f_unit, '(l)')   rsp_info%xlr_skipee
    write(f_unit, '(i8)')  rsp_info%nr_rsp_vectors
    write(f_unit, '(i8)')  rsp_info%nr_rsp_vectors_ee
    write(f_unit, '(i8)')  rsp_info%nr_rsp_vectors_ep
    !write(f_unit, '(I8)')  rsp_info%ndim_ee_total
    !write(f_unit, '(I8)')  rsp_info%ndim_ep_total

  end subroutine

! ------------------------------------------------------------------------------
  subroutine write_allocatables_to_helper(f_unit, rsp_info, ivec_selected)
    integer, intent(in) :: f_unit
    type(response_info_t), intent(in) :: rsp_info
    integer, intent(in), optional :: ivec_selected
    integer :: ivec, nr_vec
  
    nr_vec = rsp_info%nr_rsp_vectors

    do ivec = 1, nr_vec
      if (present(ivec_selected)) then
        if (ivec .ne. ivec_selected) cycle
      end if
      write(f_unit, '(5i8,a20,a5,i8,2es28.18,8i8)') &
                    ivec,                        &
                    rsp_info%nzxope(ivec),       &
                    rsp_info%nzxopeq(ivec),      &
                    rsp_info%nzxopp(ivec),       &
                    rsp_info%nzxoppq(ivec),      &
                    rsp_info%rspvec_lab(ivec),   &
                    rsp_info%rspvec_type(ivec),  &
                    rsp_info%rspvec_ih(ivec),    &
                    rsp_info%rspvec_freq(ivec),  &
                    rsp_info%rspvec_rnorm(ivec), &
                    rsp_info%rspvec_fsym(ivec),  &
                    rsp_info%rspvec_irep(ivec),  &
                    rsp_info%rspvec_trsym(ivec), &
                    rsp_info%rspvec_ifckop(ivec),&
                    rsp_info%jz_spfexc(ivec),    &
                    rsp_info%jrep_spfexc(ivec),  &
                    rsp_info%from_full_to_ee(ivec), &
                    rsp_info%from_full_to_ep(ivec)
    end do
  
  end subroutine

!! ------------------------------------------------------------------------------
!  subroutine construct_perturbed_dmat(level,dmat,bmat,  &
!                                      gp_info,rsp_info,&
!                                      cmo)
!
!    character(len=4),intent(in) :: level
!    type(gp_info_t), intent(in)  :: gp_info
!    type(response_info_t), intent(in) :: rsp_info
!    real(kind=8), intent(in) :: bmat(gp_info%nmo, gp_info%nmo, gp_info%nz, rsp_info%nr_rsp_vectors)
!
!    real(kind=8), intent(out) :: dmat(gp_info%nao*gp_info%nao,  gp_info%nz, rsp_info%nr_rsp_vectors)
!    real(kind=8), optional, intent(in) :: cmo(*)
!
!    integer :: ivec
!    integer :: irep,iopsy,itim
!    integer :: iprint
!    integer :: i1, i2, i, j, iz, ns, is, js, ni, ii, ji
!    logical :: debug_print
!!
!!   on input:
!!   level - fixme
!!   bmat - a matrix with scattered response parameters
!!   cmo  - MOs
!!
!!   on output:
!!   dmat - perturbed density matrix:
!!          symmetrized or antisymmetrized depending on rsp_info%ih
!!
!!   code adapted after DMOLI3 subroutine from prp/pamrvc.F90
!!
!
!    iprint = 0  ! needed for qtrans
!    debug_print = .false.
!
!    if (debug_print) then
!      write(*, *) 'rsp_info in construct_perturbed_dmat'
!      call fde_print_rsp_info(rsp_info)
!      write(*, *) 'cmo entering construct_perturbed_dmat'
!      !call dirac_print_arr(cmo, gp_info%nao, gp_info%nmo, gp_info%nz, &
!      !                     1, 6)
!    end if
!
!
!!   construct half dmat = Ci*bmat*dag(Ca)
!
!    call dzero(dmat, gp_info%n2bbasxq*rsp_info%nr_rsp_vectors)
!
!    do ivec = 1, rsp_info%nr_rsp_vectors
!
!       if (debug_print) then
!         write(*, *) 'bmat entering construct_perturbed_dmat for ivec=', ivec
!         write(*, *) 'ih = ', rsp_info%rspvec_ih(ivec)
!         write(*, *) 'tr_sym = ', rsp_info%rspvec_trsym(ivec)
!         write(*, *) 'irep   = ', rsp_info%rspvec_irep(ivec)
!         !call dirac_print_arr(bmat(1,1,1,ivec), &
!         !                     gp_info%nmo, gp_info%nmo, gp_info%nz, &
!         !                     1, 6)
!       end if
!
!       irep   = rsp_info%rspvec_irep(ivec)
!       iopsy  = rsp_info%rspvec_ifckop(ivec) ! 1-gerade, 2-ungerade perturbation
!       itim   = rsp_info%rspvec_trsym(ivec)
!
!       do i1 = 1, gp_info%nfsym
!          i2 = mod(i1+iopsy,2) + 1
!
!!         second index: inactive index:
!          ni = gp_info%nish(i2)
!          ii = gp_info%iorb(i2) + gp_info%npsh(i2) + 1
!          ji = gp_info%icmoq(i2) + gp_info%npsh(i2)*gp_info%nfbas(i2,0)+1
!
!!         first index: secondary index:
!          if (rsp_info%rspvec_type(ivec)(1:1) == 'E') then
!!            .. secondary index is electronic:
!             ns = gp_info%nesh(i1)
!             is = gp_info%iorb(i1) + gp_info%npsh(i1) + 1
!             js = gp_info%icmoq(i1) + gp_info%npsh(i1)*gp_info%nfbas(i1, 0) + 1
!          else if (rsp_info%rspvec_type(ivec)(1:1) == 'P') then
!!            .. secondary index is positronic:
!             ns = gp_info%npsh(i1)
!             is = gp_info%iorb(i1) + 1
!             js = gp_info%icmoq(i1) + 1
!          end if
!
!          if (ni > 0 .and. ns > 0) then
!
!             if (debug_print) then
!                write(*, *) 'check dimensions before qtrans in construct_perturbed_dmat:'
!                write(*, *) '.. ivec, i1, i2 = ', ivec, i1, i2
!                write(*, *) '.. ni, ii, ji = ', ni, ii, ji
!                write(*, *) '.. ns, is, js = ', ns, is, js
!                write(*, *) '.. gp_info%nfbas(i1,0)',gp_info%nfbas(i1,0)
!                write(*, *) '.. gp_info%nfbas(i2,0)',gp_info%nfbas(i2,0)
!                write(*, *) '.. gp_info%norb(i1)',gp_info%norb(i1)
!                write(*, *) '.. gp_info%norb(i2)',gp_info%norb(i2)
!                write(*, *) '.. gp_info%iorb(i1)',gp_info%iorb(i1)
!                write(*, *) '.. gp_info%i2basx(i1,i2)',gp_info%i2basx(i1,i2)
!                write(*, *) '.. gp_info%nao',gp_info%nao
!                write(*, *) '.. gp_info%nmo',gp_info%nmo
!                write(*, *) '.. gp_info%nz',gp_info%nz
!                write(*, *) '.. gp_info%ipqtoq(1, irep)',gp_info%ipqtoq(1, irep)
!                write(*, *) '.. gp_info%ipqtoq(1, 0)',gp_info%ipqtoq(1, 0)
!                write(*, *) '.. iprint',iprint
!             end if
!
!            call qtrans90('MOAO',                                              &
!                          'S',                                                 &
!                          0.0d0,                                               &
!                          gp_info%nfbas(i1,0), gp_info%nfbas(i2,0),          &
!                          ns, ni,                                              &
!!                         ------------------------------------------------------
!                          dmat(gp_info%i2basx(i1,i2)+1, 1, ivec),             &
!                          gp_info%nao, gp_info%nao, gp_info%nz,             &
!                          gp_info%ipqtoq(1, irep),                            &
!!                         ------------------------------------------------------
!                          bmat(is, ii, 1, ivec),                               &
!                          gp_info%nmo, gp_info%nmo, gp_info%nz,             &
!                          gp_info%ipqtoq(1, irep),                            &
!!                         ------------------------------------------------------
!                          cmo(js),                                             &
!                          gp_info%nfbas(i1,0), gp_info%norb(i1), gp_info%nz,&
!                          gp_info%ipqtoq(1, 0),                               &
!!                         ------------------------------------------------------
!                          cmo(ji),                                             &
!                          gp_info%nfbas(i2,0), gp_info%norb(i2), gp_info%nz,&
!                          gp_info%ipqtoq(1, 0),                               &
!!                         ------------------------------------------------------
!                          iprint)
!          end if
!
!       end do
!
!       if (debug_print) then
!          write(*, *) 'partial perturbed dmat (in construct_perturbed_dmat) for ivec = ', ivec
!          !call dirac_print_arr(dmat(1,1,ivec), gp_info%nao, gp_info%nao, gp_info%nz, &
!          !                     1, 6, &
!          !                     gp_info%ipqtoq(1,irep))
!       end if
!
!
!!      construct full dmat = (half dmat) +/- dag(half dmat)
!!      take symmetric/antisymmetric combinations depending on the hermiticity of response vectors
!       call fde_pertdenmat_sym(mat_dim = gp_info%nao,   &
!                               mat_nz  = gp_info%nz,    &
!                               mat     = dmat(1,1,ivec), &
!                               ih      = rsp_info%rspvec_ih(ivec),    &
!                               irep    = rsp_info%rspvec_irep(ivec),  &
!                               ipqtoq  = gp_info%ipqtoq,  &
!                               ihqmat  = gp_info%ihqmat,  &
!                               imat    = ivec)
!
!       if (debug_print) then
!          write(*, *) 'Symmetrized perturbed dmat for ivec = ', ivec
!!         fixme
!!         in FDE, the elements of this matrix will be 2 times smaller than in VISUAL;
!!         because we multiply it by 2 in emb_dirac_integration.F90
!          !call dirac_print_arr(dmat(1,1,ivec), gp_info%nao, gp_info%nao, gp_info%nz, &
!          !                     1, 6, &
!          !                     gp_info%ipqtoq(1,irep))
!       end if
!
!    end do
!
!  end subroutine
!
!
!! ------------------------------------------------------------------------------
!  subroutine fde_pertdenmat_sym(mat_dim,mat_nz,mat, &
!                                ih,irep,ipqtoq,ihqmat,imat)
!
!    integer, intent(in)    :: mat_dim
!    integer, intent(in)    :: mat_nz
!    real(8), intent(inout) :: mat(mat_dim, mat_dim, mat_nz)
!    integer, intent(in)    :: ih
!    integer, intent(in)    :: irep
!    integer, intent(in)    :: ipqtoq(4,0:7)
!    integer, intent(in)    :: ihqmat(4,-1:1)
!    integer, intent(in)    :: imat
!    integer :: i, j, iz, iq, im
!    real(8) :: fullmat
!    real(8) :: f
!!
!!   after dmoli3 in prp/pamrvc.F and fulmat in gp/gpsaue.F
!!
!    do iz = 1, mat_nz
!      iq = ipqtoq(iz, irep)
!      im = ihqmat(iq, ih)
!      if (im .eq. 1) then
!        f = 1.0d0
!      else if (im .eq. 2) then
!        f =-1.0d0
!      end if
!
!      do j = 1, mat_dim
!        do i = j, mat_dim
!          fullmat = mat(i, j, iz) + f*mat(j, i, iz)
!          mat(i, j, iz) =  fullmat
!          mat(j, i, iz) =  f*fullmat
!        end do
!      end do
!    end do
!
!   end subroutine fde_pertdenmat_sym

! ------------------------------------------------------------------------------
  subroutine read_rspdim_data_from_response_file(f_name, gp_info, rsp_info)
    character(*), intent(in) :: f_name
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    integer :: f_unit
    integer :: nr_vectors
    integer :: ivec
    integer :: ivec_ee, ivec_ep
    logical :: vec_pp, vec_pn, rspvec_update
    character(len=16) :: lab1

    rspvec_update = .false.
    f_unit = open_rsp_file(f_name)

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
                                      gp_info = gp_info, &
                                      rsp_info = rsp_info, &
                                      rspvec_info_update = rspvec_update, &
                                      vec_pp = vec_pp, &
                                      vec_pn = vec_pn)

       if (vec_pp) then
          ivec_ee = ivec_ee + 1
          if (rspvec_update) then
             rsp_info%from_full_to_ee(ivec) = ivec_ee
          end if
       end if
       if (vec_pn) then
          ivec_ep = ivec_ep + 1
          if (rspvec_update) then
             rsp_info%from_full_to_ep(ivec) = ivec_ep
          end if
       end if

       go to 1
    end if

 2  continue

    nr_vectors    = ivec
    rsp_info%nr_rsp_vectors = ivec
    rsp_info%nr_rsp_vectors_ee = ivec_ee
    rsp_info%nr_rsp_vectors_ep = ivec_ep
    rsp_info%rsp_diminfo_is_read = .true.

    close(f_unit, status = 'keep')
  end subroutine


! ------------------------------------------------------------------------------
  subroutine read_rspvec_from_response_file(f_name, gp_info, rsp_info, &
                  rspvec_info_full_update, &
                  rspvec_id, &
                                          bmat,             &
                                          rspvec_to_dxfile, &
                                          dxfile_unit,      &
                                          dxfile_name,      &
                                          fullvec_to_memory,&
                                          fullvec,          &
                                          ibeig,            &
                                          ind_ee,           &
                                          ind_ep,           &
                                          get_jz_spfexc)
    character(*), intent(in) :: f_name
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    logical, intent(in) :: rspvec_info_full_update
    integer, intent(in), optional :: rspvec_id

    real(8), optional, intent(out) :: bmat(*)
      logical,     intent(in),   optional :: rspvec_to_dxfile
      integer,     intent(in),   optional :: dxfile_unit
      character(*),intent(in),   optional :: dxfile_name
      logical,     intent(in),   optional :: fullvec_to_memory
      real(8),     intent(inout),optional :: fullvec(*)
      integer,     intent(in),   optional :: ibeig(*)
      integer,     intent(out),  optional :: ind_ee(*)
      integer,     intent(out),  optional :: ind_ep(*)
      logical,     intent(in),   optional :: get_jz_spfexc

      integer :: f_unit
      integer           :: ivec, ivec_selected, ioff, fullvec_off
      character(len=16) :: lab1
!  ----------------------------------------------------------------------------

    f_unit = open_rsp_file(f_name)
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
       if (present(rspvec_id)) then
         ivec_selected = rspvec_id
         if (ivec_selected .ne. ivec) then
           read(f_unit)
           read(f_unit)
           ioff = ioff + fullvec_off
           go to 1
         end if
       end if
       backspace(f_unit)
       write(*, *) 'READING RECORD ivec = ', ivec
       call read_one_record_of_pamxvc(f_unit, &
                                      read_dim_only = .false., &
                                      ivec     = ivec, &
                                      gp_info  = gp_info, &
                                      rsp_info = rsp_info, &
                                      bmat     = bmat, &
                                      rspvec_info_update = rspvec_info_full_update, &
                                      fullvec_to_memory = fullvec_to_memory,&
                                      fullvec = fullvec(ioff), &
                                      fullvec_off = fullvec_off, &
                                      ibeig = ibeig, &
                                      ind_ee = ind_ee, &
                                      ind_ep = ind_ep, &
                                      get_jz_spfexc = get_jz_spfexc)
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
                                       bmat,             &
                                       rspvec_info_update,  &
                                       fullvec_to_memory,&
                                       fullvec,          &
                                       fullvec_off,      &
                                       ibeig,            &
                                       ind_ee,           &
                                       ind_ep,           &
                                       get_jz_spfexc)

    integer, intent(in) :: f_unit
    logical, intent(in) :: read_dim_only
    integer, intent(in) :: ivec
    type(gp_info_t), intent(in) :: gp_info
    type(response_info_t), intent(inout) :: rsp_info
    logical, intent(in),  optional       :: rspvec_info_update
    logical, intent(in),  optional       :: fullvec_to_memory
    logical, intent(in),  optional       :: get_jz_spfexc
    logical, intent(out), optional :: vec_pp
    logical, intent(out), optional :: vec_pn
    real(8), intent(inout), optional     :: fullvec(*)
    real(8), intent(out), optional       :: bmat(*)
    integer, intent(in),  optional       :: ibeig(*)
    integer, intent(out), optional       :: fullvec_off
    integer, intent(out), optional       :: ind_ee(*)
    integer, intent(out), optional       :: ind_ep(*)

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
    integer, allocatable :: orbital_rotation_indices(:, :)
    real(8), allocatable :: vector(:, :)
    integer :: ii,jj
!
!   read ONE record of PAMXVC file;
!   save information about ONE response vector to rsp_info
!

    debug_print = .true.

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
      write(*, *) 'Solution vector nr. ', ivec, ' read from ', rsp_info%rspvec_file
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
          rsp_info%nzxope(ivec)  = length
          rsp_info%nzxopeq(ivec) = lenrec
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
          rsp_info%nzxopp(ivec)  = length
          rsp_info%nzxoppq(ivec) = lenrec
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
      rsp_info%rspvec_lab(ivec)    = lab1
      rsp_info%rspvec_type(ivec)   = typ1
      rsp_info%rspvec_freq(ivec)   = freq1
      rsp_info%rspvec_rnorm(ivec)  = rnorm
      rsp_info%rspvec_ih(ivec)     = ih    
      rsp_info%rspvec_fsym(ivec)   = jsymop1
      rsp_info%rspvec_irep(ivec)   = jsymop1-1
      rsp_info%rspvec_trsym(ivec)  = jtimop1
      rsp_info%rspvec_ifckop(ivec) = gp_info%jbtof(rsp_info%rspvec_irep(ivec), 1)

      if (debug_print) then
        write(*, *) 'rsp_info%rspvec_lab(ivec)   ', rsp_info%rspvec_lab(ivec)
        write(*, *) 'rsp_info%rspvec_type(ivec)  ', rsp_info%rspvec_type(ivec)
        write(*, *) 'rsp_info%rspvec_ih(ivec)    ', rsp_info%rspvec_ih(ivec)
        write(*, *) 'rsp_info%rspvec_fsym(ivec)  ', rsp_info%rspvec_fsym(ivec)
        write(*, *) 'rsp_info%rspvec_irep(ivec)  ', rsp_info%rspvec_irep(ivec)
        write(*, *) 'rsp_info%rspvec_trsym(ivec) ', rsp_info%rspvec_trsym(ivec)
        write(*, *) 'rsp_info%rspvec_ifckop(ivec)', rsp_info%rspvec_ifckop(ivec)
      end if
    end if

    if (read_dim_only) then
!     we only want to count the number of vectors
      read(f_unit)
      read(f_unit)
      return
    else
!     allocate and read response vector
      allocate(vector(length, gp_info%nz))
      vector=0.0d0
      read(f_unit) vector
      
!     allocate and read orbital rotation array
      allocate(orbital_rotation_indices(2, length))
      orbital_rotation_indices=0
      read(f_unit) orbital_rotation_indices
      
!     fixme: test this:
         !if (freq1 .ne. 0.0d0) then
         !   call dscal(lenrec, 1/freq1, vector, 1)
         !end if
      
      if (present(get_jz_spfexc)) then
!       jz is needed for spinfree case
!       here the code follows prp/pamrvc.F: jrep_spfexc is set to 0
        if (get_jz_spfexc) then
          if (gp_info%spinfr .and. (rsp_info%rspvec_trsym(ivec)==1)) then
            if (rsp_info%rspvec_type(ivec)(1:1) == 'E') then
              call spfexc(vector,                    &
                          orbital_rotation_indices,  &
                          rsp_info%nzxope(ivec),     &
                          ibeig,                     &
                          rsp_info%jz_spfexc(ivec),  &
                          rsp_info%jrep_spfexc(ivec))
      
            else if (rsp_info%rspvec_type(ivec)(1:1) == 'P') then
              call spfexc(vector,                    &
                          orbital_rotation_indices,  &
                          rsp_info%nzxopp(ivec),     &
                          ibeig,                     &
                          rsp_info%jz_spfexc(ivec),  &
                          rsp_info%jrep_spfexc(ivec))
            end if
          end if
        end if
      end if
    
      if (rspvec_update) rsp_info%rspvec_data_is_read(ivec) = .true.

!     fixme test offsets!
      if (present(ind_ee) .and. rsp_info%rspvec_type(ivec)(1:1) == 'E') then
!       save information on 'ee' orbital rotation indices
        ioff=length*2*(ivec-1)
        call icopy(length*2,                        &
                   orbital_rotation_indices,        &
                   1,                               &
                   ind_ee(1+ioff:ioff+length*2),    &
                   1)
      end if
      if (present(ind_ep) .and. rsp_info%rspvec_type(ivec)(1:1) == 'P') then
!       save information on 'ep' orbital rotation indices
        ioff=length*2*(ivec-1)
        call icopy(length*2,                        &
                   orbital_rotation_indices,        &
                   1,                               &
                   ind_ep(1+ioff:ioff+length*2),    &
                   1)
      end if
!
!      save the currently read vector
!
!      ioff=length*gp_info%nz*(ivec-1)
!      call dcopy(length*gp_info%nz,              &
!                 vector,                 &
!                 1,                               &
!                 rspvec(1+ioff:ioff+length*gp_info%nz),  &
!                 1)
!

      if (present(fullvec_to_memory)) then
        if (fullvec_to_memory) then
          call dcopy(length*gp_info%nz, &
                     vector,             &
                     1,                  &
                     fullvec,            &
                     1)
          fullvec_off = length*gp_info%nz
        end if
      end if

      if (present(bmat)) then
!        scatter vector -> matrix
         if (debug_print) then
            write(*, *) '.... response_vector as a vector=', ivec
            !call dirac_print_arr(vector, length, 1, gp_info%nz, &
            !                     1, 6)
         end if
         ioff=gp_info%nmo*gp_info%nmo*gp_info%nz*(ivec-1)+1
         call scatter_vector_onto_matrix(length, &
                                         gp_info%nmo, gp_info%nz,&
                                         gp_info%ipqtoq,             &
                                         gp_info%ihqmat,             &
                                         rsp_info%rspvec_irep(ivec),  &
                                         rsp_info%rspvec_ih(ivec),    &
                                         orbital_rotation_indices,    &
                                         vector,                      &
                                         bmat(ioff))
         if (debug_print) then
            write(*, *) '.... response_vector as a matrix: ',ivec
            !call dirac_print_arr(bmat(ioff), gp_info%nmo, gp_info%nmo, gp_info%nz, &
            !                     1, 6, &
            !                     gp_info%ipqtoq(1,jsymop1-1))
         end if
      end if

      deallocate(vector)
      deallocate(orbital_rotation_indices)

    end if

   end subroutine


! ------------------------------------------------------------------------------
  subroutine scatter_vector_onto_matrix(length, &
                                        ndim, nz, &
                                        ipqtoq, ihqmat, &
                                        irep, ih, &
                                        orbital_rotation_indices, &
                                        response_vector, &
                                        matrix)

    integer, intent(in)  :: length, ndim, nz, irep, ih
    integer, intent(in)  :: ipqtoq(4,0:7)
    integer, intent(in)  :: ihqmat(4,-1:1)
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: response_vector(length, *)
    real(8), intent(out) :: matrix(ndim, ndim, nz)
    integer :: i, s, is, iz, iq, im
    real(8) :: f
!
!   reused scatter_vector subroutine from visual/visual.F90
!   the resulting matrix is already symmetric/antisymmetric
!

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        iq = ipqtoq(iz, irep)
        im = ihqmat(iq, ih)

        if (im.eq.1) then
          f = -1.0d0
        else if (im.eq.2) then
          f =  1.0d0
        end if

        matrix(s, i, iz) = response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         + f*response_vector(is, iz)
      end do
    end do

!   this matrix is scaled by -1.0 in dmolitx1 subroutine, so we do the same here:
    call dscal(ndim*ndim*nz, -1.0d0, matrix, 1)

  end subroutine


! ------------------------------------------------------------------------------
  function open_rsp_file(fname)
!   mv to IO
    character(*), intent(in) :: fname
    integer :: open_rsp_file
    integer :: funit
    logical :: exists, is_opened, old_dx, debug_print
    character(len=100) :: errmsg

    inquire(file=fname, exist=exists, opened=is_opened)
    if (exists) then
      if (.not.is_opened) then
        open(file = fname, &
             newunit = funit,   &
             status = 'old',    &
             access = 'sequential',&
             form   = 'unformatted', &
             position = 'rewind')
      else
        inquire(file=fname, number=funit)
        rewind(funit)
      end if
    else
      write(errmsg, *) 'file ', fname, ' does not exist',&
                       ' Stopping now!'
      call quit(errmsg)
    end if
    open_rsp_file = funit
  end function


! ------------------------------------------------------------------------------
  subroutine read_gp_data_from_cbs(self)
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
!
!   read information from common blocks to a derived type
!   fixme - this information can be read from CHECKPOINT, but not
!   everything is avalable there, so I use the old code for now
!
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


! ------------------------------------------------------------------------------
  subroutine erase_rsp_data(self)
    class(response_info_t), intent(out) :: self
  end subroutine

! ------------------------------------------------------------------------------
  subroutine erase_gp_data(self)
    class(gp_info_t), intent(out) :: self
  end subroutine

! ----------------------------------------------------------------------------
  subroutine read_rspdim_data_from_cbs(self)
#include "implicit.h"
#include "dcbxlr.h"
#include "dcbxrs.h"
#include "dcbxpr.h"
    type(response_info_t), intent(out) :: self
    integer :: n_ee, n_ep
    logical :: debug_print
!
!   read information from common blocks to a derived type:
!     * transferred variables are set in SETXLR subroutine,
!       and passed through dcbxlr.h,dcbxrs.h,dcbxpr.h
!     * only nonallocatable variables of response_info_t are set
!
    debug_print = .false.

    call self%erase()

    self%nfreq           = nfreq
    self%ndamp           = ndamp
    self%static          = static
    self%nr_B_operators  = nlrbpt
    self%xlr_skipep      = xlr_skipep
    self%xlr_skipee      = xlr_skipee

    n_ee = 1
    n_ep = 1
    if (xlr_skipep) n_ep = 0
    if (xlr_skipee) n_ee = 0

    if (static) then
       !self%nr_rsp_vectors_ee = n_ee*nfreq*nlrbpt*ndamp
       !self%nr_rsp_vectors_ep = n_ep*nfreq*nlrbpt*ndamp
       self%nr_rsp_vectors_ee = n_ee*nfreq*nlrbpt
       self%nr_rsp_vectors_ep = n_ep*nfreq*nlrbpt
    else
       ! we will have vectors of + and - time-reversal symmetry
       !self%nr_rsp_vectors_ee = 2*n_ee*nfreq*nlrbpt*ndamp
       !self%nr_rsp_vectors_ep = 2*n_ep*nfreq*nlrbpt*ndamp
       self%nr_rsp_vectors_ee = 2*n_ee*nfreq*nlrbpt
       self%nr_rsp_vectors_ep = 2*n_ep*nfreq*nlrbpt
    end if
    self%nr_rsp_vectors    = self%nr_rsp_vectors_ee + self%nr_rsp_vectors_ep

    if (debug_print) then
       write(*,*) '--- output from fde_dirac_xlr_data_interface: ---'
       write(*,*) 'self%nfreq             = ', self%nfreq
       write(*,*) 'self%ndamp             = ', self%ndamp
       write(*,*) 'self%static            = ', self%static
       write(*,*) 'self%nr_B_operators    = ', self%nr_B_operators
       write(*,*) 'self%xlr_skipep        = ', self%xlr_skipep
       write(*,*) 'self%xlr_skipee        = ', self%xlr_skipee
       write(*,*) 'self%nr_rsp_vectors_ee = ', self%nr_rsp_vectors_ee
       write(*,*) 'self%nr_rsp_vectors_ep = ', self%nr_rsp_vectors_ep
       write(*,*) 'self%nr_rsp_vectors    = ', self%nr_rsp_vectors
    end if

    self%rsp_diminfo_is_read = .true.
  end subroutine

!  ----------------------------------------------------------------------------
   subroutine print_nonallocatables_from_helper(gp_info, rsp_info)
!  ----------------------------------------------------------------------------
      type(gp_info_t), intent(in)       :: gp_info
      type(response_info_t), intent(in) :: rsp_info
!  ----------------------------------------------------------------------------

      write(*,*) '   gp_info%nao                ',  gp_info%nao
      write(*,*) '   gp_info%nmo                ',  gp_info%nmo
      write(*,*) '   gp_info%nz                 ',  gp_info%nz
      write(*,*) '   gp_info%spinfr             ',  gp_info%spinfr
      write(*,*) '   rsp_info%ndamp              ',  rsp_info%ndamp
      write(*,*) '   rsp_info%nfreq              ',  rsp_info%nfreq
      write(*,*) '   rsp_info%nr_B_operators     ',  rsp_info%nr_B_operators
      write(*,*) '   rsp_info%static             ',  rsp_info%static
      write(*,*) '   rsp_info%xlr_skipep         ',  rsp_info%xlr_skipep
      write(*,*) '   rsp_info%xlr_skipee         ',  rsp_info%xlr_skipee
      write(*,*) '   rsp_info%nr_rsp_vectors     ',  rsp_info%nr_rsp_vectors
      write(*,*) '   rsp_info%nr_rsp_vectors_ee  ',  rsp_info%nr_rsp_vectors_ee
      write(*,*) '   rsp_info%nr_rsp_vectors_ep  ',  rsp_info%nr_rsp_vectors_ep
      !read(f_unit, '(I8)')  rsp_info%ndim_ee_total
      !read(f_unit, '(I8)')  rsp_info%ndim_ep_total

   end subroutine print_nonallocatables_from_helper


! ------------------------------------------------------------------------------
  subroutine print_xlr_info(tmp)
    type(response_info_t), intent(in) :: tmp
    integer :: i

    write(*,*) 'tmp%nr_rsp_vectors             =', tmp%nr_rsp_vectors
    write(*,*) 'tmp%nr_rsp_vectors_ee          =', tmp%nr_rsp_vectors_ee
    write(*,*) 'tmp%nr_rsp_vectors_ep          =', tmp%nr_rsp_vectors_ep
    write(*,*) 'tmp%ndim_ee_total              =', tmp%ndim_ee_total
    write(*,*) 'tmp%ndim_ep_total              =', tmp%ndim_ep_total
    write(*,*) 'tmp%nr_B_operators             =', tmp%nr_B_operators
    write(*,*) 'tmp%ndamp                      =', tmp%ndamp
    write(*,*) 'tmp%nfreq                      =', tmp%nfreq
    write(*,*) 'tmp%static                     =', tmp%static
    write(*,*) 'tmp%xlr_skipep                 =', tmp%xlr_skipep
    write(*,*) 'tmp%xlr_skipee                 =', tmp%xlr_skipee

    if (tmp%nr_rsp_vectors > 0) then
      do i = 1, tmp%nr_rsp_vectors
        write(*,*) 'data for response vector ', i
        if (associated(tmp%rspvec_data_is_read) .and. &
            tmp%rspvec_data_is_read(i)) then
          write(*,*) 'tmp%nzxope          =', tmp%nzxope(i)
          write(*,*) 'tmp%nzxopp          =', tmp%nzxopp(i)
          write(*,*) 'tmp%nzxopeq         =', tmp%nzxopeq(i)
          write(*,*) 'tmp%nzxoppq         =', tmp%nzxoppq(i)
          write(*,*) 'tmp%rspvec_lab      =', tmp%rspvec_lab(i)
          write(*,*) 'tmp%rspvec_type     =', tmp%rspvec_type(i)
          write(*,*) 'tmp%rspvec_ih       =', tmp%rspvec_ih(i)
          write(*,*) 'tmp%rspvec_freq     =', tmp%rspvec_freq(i)
          write(*,*) 'tmp%rspvec_rnorm    =', tmp%rspvec_rnorm(i)
          write(*,*) 'tmp%rspvec_fsym     =', tmp%rspvec_fsym(i)
          write(*,*) 'tmp%rspvec_irep     =', tmp%rspvec_irep(i)
          write(*,*) 'tmp%rspvec_trsym    =', tmp%rspvec_trsym(i)
          write(*,*) 'tmp%rspvec_ifckop   =', tmp%rspvec_ifckop(i)
          write(*,*) 'tmp%jz_spfexc       =', tmp%jz_spfexc(i)
          write(*,*) 'tmp%jrep_spfexc     =', tmp%jrep_spfexc(i)
          write(*,*) 'tmp%from_full_to_ee =', tmp%from_full_to_ee(i)
          write(*,*) 'tmp%from_full_to_ep =', tmp%from_full_to_ep(i)
        end if
      end do
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine print_gp_info(tmp)
    type(gp_info_t), intent(in) :: tmp

    write(*,*) 'tmp%nmo        = ', tmp%nmo
    write(*,*) 'tmp%nao        = ', tmp%nao
    write(*,*) 'tmp%n2bbasxq   = ', tmp%n2bbasxq
    write(*,*) 'tmp%n2orbxq    = ', tmp%n2orbxq
    write(*,*) 'tmp%nz         = ', tmp%nz
    write(*,*) 'tmp%nocc       = ', tmp%nocc
    write(*,*) 'tmp%norb       = ', tmp%norb
    write(*,*) 'tmp%npsh       = ', tmp%npsh
    write(*,*) 'tmp%nish       = ', tmp%nish
    write(*,*) 'tmp%nesh       = ', tmp%nesh
    write(*,*) 'tmp%nfsym      = ', tmp%nfsym
    write(*,*) 'tmp%nbsym      = ', tmp%nbsym
    write(*,*) 'tmp%nfbas      = ', tmp%nfbas
    write(*,*) 'tmp%i2basx     = ', tmp%i2basx
    write(*,*) 'tmp%i2orbx     = ', tmp%i2orbx
    write(*,*) 'tmp%iorb       = ', tmp%iorb
    write(*,*) 'tmp%jbtof      = ', tmp%jbtof
    write(*,*) 'tmp%ihrmop     = ', tmp%ihrmop
    write(*,*) 'tmp%icmoq      = ', tmp%icmoq
    write(*,*) 'tmp%ihqmat     = ', tmp%ihqmat
    write(*,*) 'tmp%ncmotq     = ', tmp%ncmotq
    write(*,*) 'tmp%ipqtoq     = ', tmp%ipqtoq
    write(*,*) 'tmp%spinfr     = ', tmp%spinfr
    write(*,*) 'tmp%intflg     = ', tmp%intflg

  end subroutine



end module
