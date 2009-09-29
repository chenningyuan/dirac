module dirac_interface

  use memory_allocator
  use matrix_defop
  use matrix_genop
  use aor_cfg
  use character_processing

  implicit none

  public print_q_vector
  public di_get_gbds
  public di_read_operator_int
  public do_dft
  public get_s
  public get_f
  public get_p
  public get_g
  public get_C
  public grcont
  public ask_for_integrals
  public response_solver
  public get_natoms
  public dipnuc_ifc
  public qdrnuc_ifc
  public gradnn_ifc
  public di_set_f77_work
  public di_select_wrk
  public di_deselect_wrk
  public daxpy_iz_jz
  public print_prp_common_blocks
  public init_dirac_interface
  public scatter_vector
  public set_orbital_rotation_indices

  save

  integer,      public              :: unit_input       = 0
  integer,      public              :: unit_output      = 0

  integer,      public              :: nr_atoms         = 0
  integer,      public              :: nr_sym_dep_nuc   = 0
  integer,      public              :: nr_sym_indep_nuc = 0
  character(6), public, allocatable :: sym_dep_nuc_name(:)

  logical,      public              :: nonzero_small_metric = .false.

  integer,      public              :: integral_flag        = 0
  logical,      public              :: include_ll_integrals = .false.
  logical,      public              :: include_ls_integrals = .false.
  logical,      public              :: include_ss_integrals = .false.
  logical,      public              :: include_gt_integrals = .false.

  logical,      public              :: include_pp_rotations = .false.
  logical,      public              :: include_pn_rotations = .false.
  logical,      public              :: diamagnetic_via_pn   = .false.

  private

  logical :: dirac_interface_data_allocated = .false.
  logical :: dirac_interface_data_set       = .false.
#ifdef VAR_MPI
  logical :: dirac_interface_data_synced    = .false.
#endif

  interface di_get_gbds
    module procedure get_G
  end interface

  interface di_read_operator_int
    module procedure get_P
  end interface

#include "priunit.h"
#include "mxcent.h"
#include "../abacus/cbinuc.h"
#include "cbihr2.h"
#include "dcbbas.h"
#include "dcbcls.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbprl.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"
#include "dipole.h"
#include "dorps.h"
#include "exeinf.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"

  integer, parameter :: max_aoproper_indices = 900
  character(16)      :: op_name_to_aoproper_index(max_aoproper_indices) = '----------------'

  character(2)       ::  i_iterator(3) = (/ 'X',  'Y',  'Z'/)
  character(2)       :: ij_iterator(6) = (/'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'/)

  integer            :: coulomb_mem_req   ! Memory requirements for the Coulomb build
  integer            :: dft_mem_req       ! Memory requirements for the DFT build
  real(8), pointer   :: f77_work(:)       ! Fortran 77 work array pointer
  integer, public    :: len_f77_work  = 0 ! Length of the work array
  integer            :: memory_usage  = 0 ! Amount of work array used
  integer            :: mem_wrk_count = 0 ! Counts

contains

  subroutine init_dirac_interface()

!   ----------------------------------------------------------------------------
    integer :: ill_int
    integer :: ils_int
    integer :: iss_int
    integer :: igt_int
!   ----------------------------------------------------------------------------

    if (.not. dirac_interface_data_allocated) then

      allocate(sym_dep_nuc_name(mxcent))

      dirac_interface_data_allocated = .true.

    end if

!   if (.not. dirac_interface_data_set) then

      unit_input       = lucmd
      unit_output      = lupri

      nr_atoms         = natoms
      nr_sym_dep_nuc   = nucdep
      nr_sym_indep_nuc = nucind
      sym_dep_nuc_name = namdep

      nonzero_small_metric = .false.
      if (mc == 2) then
        nonzero_small_metric = .true.
      end if
      if (bss .or. twocomp) then
        nonzero_small_metric = .false.
      end if

      ill_int = iand(intgen,   1)
      ils_int = iand(intgen/2, 1)
      iss_int = iand(intgen/4, 1)
      igt_int = iand(intgen/8, 1)

      integral_flag = 0
      if (ill_int == 1) then
        integral_flag = integral_flag + 1
        include_ll_integrals = .true.
      end if
      if (ils_int == 1) then
        integral_flag = integral_flag + 2
        include_ls_integrals = .true.
      end if
      if (iss_int == 1) then
        integral_flag = integral_flag + 4
        include_ss_integrals = .true.
      end if
      if (igt_int == 1) then
        integral_flag = integral_flag + 8
        include_gt_integrals = .true.
      end if

      if (.not. aor_cfg_skip_pp) then
        include_pp_rotations = .true.
      end if

      if (nonzero_small_metric .and. .not. aor_cfg_skip_pn) then
        include_pn_rotations = .true.
        diamagnetic_via_pn   = .true.
      end if

      dirac_interface_data_set = .true.

!   end if

  end subroutine

  subroutine daxpy_iz_jz(f, A, iz, B, jz)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: f
    integer,      intent(in)    :: iz, jz
    type(matrix), intent(inout) :: A, B
!   ----------------------------------------------------------------------------
    integer                     :: n
!   ----------------------------------------------------------------------------

    n = A%nrow*A%ncol
    call daxpy(n, f, A%elms((iz-1)*n + 1), 1, B%elms((jz-1)*n + 1), 1)

  end subroutine

  function index_on_aoproper(op_name)

!   ----------------------------------------------------------------------------
    integer                   :: index_on_aoproper
    character(*), intent(in)  :: op_name
!   ----------------------------------------------------------------------------
    integer                   :: i
!   ----------------------------------------------------------------------------

    index_on_aoproper = 0
    do i = 1, max_aoproper_indices
      if (op_name_to_aoproper_index(i) == op_name(1:len_trim(op_name))) then
        index_on_aoproper = i
      end if
    end do

    if (index_on_aoproper == 0) then
      call quit('index_on_aoproper: index for ' &
              //op_name(1:len_trim(op_name))    &
              //' not found')
    end if

  end function

  SUBROUTINE di_set_f77_work(wrk,lwork)
     real(8),target :: wrk(:)
     integer            :: lwork, ltempwork
     logical            :: ldummy1,ldummy2
     len_f77_work = lwork
     f77_work     => wrk
!    call get_linsca_work(ldummy1,ldummy2,ltempwork,lwork)
!    call di_set_coulomb_mem_req(ltempwork)
!    ! Simen: Refine!
!    call di_set_dft_mem_req(ltempwork)
  END SUBROUTINE di_set_f77_work

  subroutine di_select_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------

    if (need_mem > len_f77_work) then

       call di_create_wrk(work, need_mem)

    else

      work => f77_work

      memory_usage  = need_mem
      mem_wrk_count = mem_wrk_count + 1

      if (mem_wrk_count > 1) call quit('counter > 1 not implemented')

    endif

  end subroutine

  subroutine di_deselect_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------

    if (need_mem > len_f77_work) then

       call di_delete_wrk(work)

    else

      nullify(work)

      memory_usage  = 0
      mem_wrk_count = mem_wrk_count - 1

    endif

  end subroutine

  subroutine di_create_wrk(work, need_mem)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
    integer, intent(in) :: need_mem
!   ----------------------------------------------------------------------------
    integer             :: error
!   ----------------------------------------------------------------------------

    allocate(work(need_mem), stat=error)

    if (error /= 0) then
      call quit('unable to allocate work space in di_create_wrk')
    end if

  end subroutine

  subroutine di_delete_wrk(work)

!   ----------------------------------------------------------------------------
    real(8), pointer    :: work(:)
!   ----------------------------------------------------------------------------

    deallocate(work)

  end subroutine

  subroutine print_prp_common_blocks()

!   ----------------------------------------------------------------------------
    integer :: i
!   ----------------------------------------------------------------------------

    do i = 1, nprps
      write(*, *) 'dcbxpr:  --------------'
      write(*, *) 'dcbxpr:  i             ', i
      write(*, *) 'dcbxpr:  prpnam(i)     ', prpnam(i)
      write(*, *) 'dcbxpr:  iprptyp(i)    ', iprptyp(i)
      write(*, *) 'dcbxpr:  iprpsym(i)    ', iprpsym(i)
      write(*, *) 'dcbxpr:  iprptim(i)    ', iprptim(i)
      write(*, *) 'dcbxpr:  facprp(1, i)  ', facprp(1, i)
      write(*, *) 'dcbxpr:  facprp(2, i)  ', facprp(2, i)
      write(*, *) 'dcbxpr:  facprp(3, i)  ', facprp(3, i)
      write(*, *) 'dcbxpr:  iprplbl(1, i) ', iprplbl(1, i)
      write(*, *) 'dcbxpr:  iprplbl(2, i) ', iprplbl(2, i)
      write(*, *) 'dcbxpr:  iprplbl(3, i) ', iprplbl(3, i)
      write(*, *) 'dcbxpr:  --------------'
    end do

    do i = 1, nprplbl
      write(*, *) 'dcbprl:  -----------'
      write(*, *) 'dcbprl:  i          ', i
      write(*, *) 'dcbprl:  prplbl(i)  ', prplbl(i)
      write(*, *) 'dcbprl:  pdoint(i)  ', pdoint(i)
      write(*, *) 'dcbprl:  iprlrep(i) ', iprlrep(i)
      write(*, *) 'dcbprl:  iprltyp(i) ', iprltyp(i)
      write(*, *) 'dcbprl:  -----------'
    end do

    do i = 1, nprpcls
      write(*, *) 'dcbcls:  ----------'
      write(*, *) 'dcbcls:  i         ', i
      write(*, *) 'dcbcls:  clscal(i) ', clscal(i)
      write(*, *) 'dcbcls:  clsint(i) ', clsint(i)
      write(*, *) 'dcbcls:  clscmb(i) ', clscmb(i)
      write(*, *) 'dcbcls:  ----------'
    end do

  end subroutine

  subroutine request_1el_integrals(label, blocks)

!   ----------------------------------------------------------------------------
    character(7), intent(in) :: label
    character(4), intent(in) :: blocks
!   ----------------------------------------------------------------------------
    integer                  :: i
!   ----------------------------------------------------------------------------

    i = nprpcls + 1

    if (i > maxcls) then
      call quit('request_1el_integrals: i > maxcls')
    end if

    clscal(i) = .true.
    clsint(i) = label
    clscmb(i) = blocks

!   important only for cartesian or spherical moment integrals
!   here just set it to something
    iordcl(i) = 0

    nprpcls = i

  end subroutine

  subroutine grcont(work,     &
                    lwork,    &
                    g,        &
                    g_dim,    &
                    geo,      &
                    london,   &
                    max_derv, &
                    iatom_in, &
                    energy,   &
                    fock,     &
                    dmat,     &
                    nr_dmat)

!   ----------------------------------------------------------------------------
    real(8)              :: work(*)
    integer, intent(in)  :: lwork

    real(8), intent(out) :: g(*)
    integer, intent(in)  :: g_dim

    logical, intent(in)  :: geo
    logical, intent(in)  :: london

    integer, intent(in)  :: max_derv
    integer, intent(in)  :: iatom_in

    logical, intent(in)  :: energy
    logical, intent(in)  :: fock

    real(8), intent(in)  :: dmat(*)
    integer, intent(in)  :: nr_dmat
!   ----------------------------------------------------------------------------
    real(8), allocatable :: tmat(:), fmat(:)

    logical              :: gencnt

    integer              :: i2typ, itype, ifctyp(8), iatom, fmat_dim

    real(8)              :: skip(3, 14, 3)
    real(8)              :: screen, f
!   ----------------------------------------------------------------------------

    dopert = .true.

    ftronv = .true.
    relgrd = .true.
    call rsetsym(.false., relgrd, .false., .false.)
    call rtroini(.false.)

    call troinv(work, lwork)

    screen = -1.0d0

!   set ifctyp
    ifctyp    = 0
    ifctyp(1) = 13 !coulomb and exchange
    ifctyp(2) = 22 !only exchange
    ifctyp(3) = 22
    ifctyp(4) = 22
    ifctyp(5) = 13 !coulomb and exchange
    ifctyp(6) = 22 !only exchange
    ifctyp(7) = 22
    ifctyp(8) = 22

    iatom = 0

    if (energy) then

      gencnt   = .true.
      fmat_dim = g_dim

      if (geo) then
        itype = 2
      else if (london) then
        itype = -5
      else
        call quit('grcont: differentiated integrals not recognized')
      end if

    else if (fock) then

      gencnt   = .false.
      fmat_dim = 3*n2bbasxq

      if (geo) then
        itype = 8
        iatom = iatom_in
      else if (london) then
        itype = -5
      else
        call quit('grcont: differentiated integrals not recognized')
      end if

    else
      call quit('programming error in call to grcont')
    end if

    f = 1.0d0
    if (fock) then
      if (geo)                        f = 0.25d0
      if (london .and. max_derv == 1) f = 0.50d0
      if (london .and. max_derv == 2) f = 1.00d0
    end if

    allocate(tmat(fmat_dim))
    allocate(fmat(fmat_dim))

    fmat = 0.0d0

    if (include_ll_integrals) then

      i2typ    = 1
      intclass = 0

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, fmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    if (include_ls_integrals) then

      i2typ    = 2
      intclass = 1

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, fmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    if (include_ss_integrals) then

      i2typ    = 3
      intclass = 2

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, fmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    call dcopy(g_dim, fmat, 1, g, 1)

    deallocate(fmat)

  end subroutine

  subroutine dipnuc_ifc(n, f, fg)

    integer       :: n
    real(8) :: f(3), fg(3, 3*n)

    call dipnuc((/0.0d0/), (/0.0d0/), 0, .true.)

    f  = dipmn
    fg = ddipn(:, 1:3*n)

  end subroutine

  subroutine get_natoms(return_nr_atoms)

!   ----------------------------------------------------------------------------
    integer, intent(out) :: return_nr_atoms
!   ----------------------------------------------------------------------------

    return_nr_atoms = nr_atoms

  end subroutine

  subroutine qdrnuc_ifc(q)

!   nuclear contribution to quadrupole moments
!   based on Kenneth's nucqdr in Dalton

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: q(6)
!   ----------------------------------------------------------------------------
    real(8)              :: qn(3, 3), x, y, z, r2
    integer              :: ix, iy, iz, i, isymop, k
!   ----------------------------------------------------------------------------

    qn = 0.0d0

    ix = iptax(1, 1)
    iy = iptax(2, 1)
    iz = iptax(3, 1)

    k = 0
    do i = 1, nucind
      do isymop = 0, maxopr
        if (iand(isymop, istbnu(i)) == 0) then
          k = k + 1

          x  = cord(1, k)
          y  = cord(2, k)
          z  = cord(3, k)

          r2 = x*x + y*y + z*z
          r2 = r2/3.0d0

          qn(ix, ix) = qn(ix, ix) + charge(i)*(x*x - r2)
          qn(iy, iy) = qn(iy, iy) + charge(i)*(y*y - r2)
          qn(iz, iz) = qn(iz, iz) + charge(i)*(z*z - r2)
          qn(ix, iy) = qn(ix, iy) + charge(i)* x*y
          qn(ix, iz) = qn(ix, iz) + charge(i)* x*z
          qn(iy, iz) = qn(iy, iz) + charge(i)* y*z

        end if
      end do
    end do

    qn(iy, ix) = qn(ix, iy)
    qn(iz, ix) = qn(ix, iz)
    qn(iz, iy) = qn(iy, iz)

    qn = 1.5d0*qn

!   ordering is: xx, yx, yy, zx, zy, zx
    q = (/qn(1, 1), qn(1:2, 2), qn(1:3, 3)/)

  end subroutine

  subroutine gradnn_ifc(n, g)

    integer       :: n
    real(8) :: g(3*n)

    iprint = 0
    maxdif = 1

    call nucrep((0.0d0), (0.0d0), (0.0d0))

    g(1:3*n) = gradnn(1:3*n)

  end subroutine

  function do_dft() result(l)

    use dirac_cfg

!   ----------------------------------------------------------------------------
    logical :: l
!   ----------------------------------------------------------------------------

    l = dirac_cfg_dft_calculation

  end function

  subroutine get_F(F)

    type(matrix), intent(inout) :: F
!   ----------------------------------------------------------------------------
    type(matrix)                :: F1, F2
    integer,      parameter     :: file_unit = 70

    call init_mat(F1, ntbas(0), ntbas(0))
    call init_mat(F2, ntbas(0), ntbas(0))

    F1%elms = 0.0d0
    F2%elms = 0.0d0

    call opnfil(file_unit, 'DFFCK1', 'OLD', 'get_F')
    call reafck(file_unit, F1%elms, .true., 1)
    close(file_unit, status = 'keep')

    call opnfil(file_unit, 'DFFCK2', 'OLD', 'get_F')
    call reafck(file_unit, F2%elms, .true., 1)
    close(file_unit, status = 'keep')

    if (.not. mat_isdef(F)) call init_mat(F, ntbas(0), ntbas(0))

    F = F1 + F2

    F1 = 0
    F2 = 0

  end subroutine

  subroutine get_G(D, G)

    type(matrix),  intent(in)    :: D
    type(matrix),  intent(inout) :: G
!   ----------------------------------------------------------------------------
    integer                      :: npos_nr
    integer,       allocatable   :: npos(:)

    real(8), pointer   :: work(:)
    integer                      :: lwork

    if (.not. mat_isdef(G)) call init_mat(G, ntbas(0), ntbas(0))
    G = 1.0d0*D !inherit symmetry information
    G%elms = 0.0d0

    call stex_get_npos(npos_nr)
    call alloc(npos, npos_nr)

!   fixme:
!   call with proper limit
!   at the moment i don't know the limit
!   simply give entire work and hope for the best
    lwork = len_f77_work
    call di_select_wrk(work, lwork)

    call twofck((/G%irep + 1/),  &
!               (/G%ih_sym/),    &
                (/0/),           &
                (/1/),           &
                  G%elms,        &
                  D%elms,        &
                  1,             &
                  npos,          &
                  integral_flag, &
                  0,             &
                  work,          &
                  lwork)

    call di_deselect_wrk(work, lwork)

    call dealloc(npos)

  end subroutine

  subroutine get_full_c(C)

!   ----------------------------------------------------------------------------
    type(matrix),  intent(inout) :: C
!   ----------------------------------------------------------------------------
    real(8), allocatable   :: eigenvalues(:)
    integer,       parameter     :: file_unit = 66
    real(8)                :: dummy_dp
    integer                      :: dummy_i
!   ----------------------------------------------------------------------------

    if (.not. mat_isdef(C)) call init_mat(C, ntbas(0), norbt)

    C%elms = 0.0d0

    call alloc(eigenvalues, norbt)

    call opnfil(file_unit, 'DFCOEF', 'OLD', 'get_full_c')

    call reacmo(file_unit,   &
                'DFCOEF',    &
                C%elms,      &
                eigenvalues, &
                dummy_i,     &
                dummy_dp,    &
                6,           &
                dummy_dp,    &
                dummy_i,     &
                dummy_i)

    close(file_unit, status = 'keep')

    call dealloc(eigenvalues)

  end subroutine

  subroutine get_C(C, i, s, g, u)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: C
    real(8),      intent(in)    :: i, s, g, u
!   ----------------------------------------------------------------------------
    type(matrix)                :: C_file

    integer                     :: nr_g_mo_ns
    integer                     :: nr_g_mo_pi
    integer                     :: nr_g_mo_ps
    integer                     :: nr_u_mo_ns
    integer                     :: nr_u_mo_pi
    integer                     :: nr_u_mo_ps
    integer                     :: nr_g_mo
    integer                     :: nr_u_mo
    integer                     :: nr_g_ao
    integer                     :: nr_u_ao
    integer                     :: k, l, m, n, ir, iw1, iw2, iz
!   ----------------------------------------------------------------------------

    nr_g_mo_ns = npsh(1)
    nr_g_mo_pi = nish(1)
    nr_g_mo_ps = nesh(1) - nr_g_mo_pi
    nr_g_mo    = nr_g_mo_ns + nr_g_mo_pi + nr_g_mo_ps
    nr_u_mo_ns = npsh(2)
    nr_u_mo_pi = nish(2)
    nr_u_mo_ps = nesh(2) - nr_u_mo_pi
    nr_u_mo    = nr_u_mo_ns + nr_u_mo_pi + nr_u_mo_ps
    nr_g_ao    = nfbas(1, 0)
    nr_u_ao    = nfbas(2, 0)

    if (.not. mat_isdef(C)) call init_mat(C, ntbas(0), norbt)
    C%elms = 0.0d0

    call get_full_c(C_file)

    ir = 0

    do iz = 1, nz

      iw1 = (iz - 1)*(nr_g_mo + nr_u_mo)*(nr_g_ao + nr_u_ao)
      k   = 0

      do m = 1, nr_g_mo_ns
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*g*s
        end do
      end do

      do m = 1, nr_g_mo_pi
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*g*i
        end do
      end do

      do m = 1, nr_g_mo_ps
        k = k + 1
        do n = 1, nr_g_ao
          ir  = ir + 1
          iw2 = iw1 + (k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*g*s
        end do
      end do

    end do

    if (nfsym == 1) return

    do iz = 1, nz

      iw1 = (iz - 1)*(nr_g_mo + nr_u_mo)*(nr_g_ao + nr_u_ao)
      k   = 0

      do m = 1, nr_u_mo_ns
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*u*s
        end do
      end do

      do m = 1, nr_u_mo_pi
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*u*i
        end do
      end do

      do m = 1, nr_u_mo_ps
        k = k + 1
        do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
          ir  = ir + 1
          iw2 = iw1 + (nr_g_mo + k - 1)*(nr_g_ao + nr_u_ao) + n
          C%elms(iw2) = C_file%elms(ir)*u*s
        end do
      end do

    end do

  end subroutine

  subroutine get_S(S)

    type(matrix),  intent(inout) :: S
!   ----------------------------------------------------------------------------

    if (.not. mat_isdef(S)) call init_mat(S, ntbas(0), ntbas(0))
    S%elms = 0.0d0

!   call gtovlt(S%elms, ssmtrc)
    call gtovlx(S%elms, ssmtrc)

  end subroutine


  subroutine get_P(op_name, P)

    character(*),  intent(in)    :: op_name
    type(matrix),  intent(inout) :: P
!   ----------------------------------------------------------------------------
    real(8), allocatable   :: ptri(:), op1int(:)
    logical,       allocatable   :: first(:)
    logical                      :: lexpst(2) = .false.
    integer                      :: op_index, file_unit, idimension
    real(8)                :: dummy_dp
    logical                      :: file_exists, debug_me = .false.

    real(8), pointer   :: work(:)
    integer            :: lwork

    file_unit = lu1int
    call opnfil(file_unit, 'AOPROPER', 'UNKNOWN', 'get_P')
    rewind file_unit

    op_index = index_on_aoproper(op_name)

    idimension = nnbbasx
    if (iprptim(op_index) == 0) then
      idimension = n2bbasx
    end if

    call alloc(ptri, nz*idimension)
    call alloc(op1int, idimension)
    allocate(first(nz))

    if (.not. mat_isdef(P)) call init_mat(P, ntbas(0), ntbas(0))
    P%elms = 0.0d0

    P%irep   = iprpsym(op_index) - 1
    P%tr_sym = iprptim(op_index)
    P%ih_sym = P%tr_sym

    lwork = n2bbasxq*2
    call di_select_wrk(work, lwork)

    call prpex2(op_index,   &
                dummy_dp,   &
                .false.,    &
                lexpst,     &
                idimension, &
                ptri,       &
                P%elms,     &
                op1int,     &
                first,      &
                work, 1, lwork, 0)

    call di_deselect_wrk(work, lwork)
    close(file_unit, status = 'keep')

    if (debug_me) then
      call print_mat(P, label = 'debug P in AO basis')
    end if

    call dealloc(ptri)
    call dealloc(op1int)
    deallocate(first)

  end subroutine

  subroutine set_orbital_rotation_indices(irep,       &
                                          length,     &
                                          from_block, &
                                          to_block,   &
                                          orbital_rotation_indices)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: irep
    integer, intent(in)  :: length
    integer, intent(in)  :: from_block
    integer, intent(in)  :: to_block
    integer, intent(out) :: orbital_rotation_indices(2, length)
!   ----------------------------------------------------------------------------
    integer, parameter   :: ns = 1
    integer, parameter   :: pi = 2
    integer, parameter   :: ps = 3

    integer, parameter   :: g  = 1
    integer, parameter   :: u  = 2

    logical              :: gu_combine(2, 2)

    integer              :: nr_mo(2, 3)       = 0
    integer              :: range_mo(2, 3, 2) = 0

    integer              :: i, s, k
    integer              :: m, m1, m2
    integer              :: gu, gu1, gu2
    integer              :: block, block1, block2

    logical              :: debug_me = .false.
!   ----------------------------------------------------------------------------

    nr_mo(g, ns) = npsh(g)
    nr_mo(g, pi) = nish(g)
    nr_mo(g, ps) = nesh(g) - nish(g)

    nr_mo(u, ns) = npsh(u)
    nr_mo(u, pi) = nish(u)
    nr_mo(u, ps) = nesh(u) - nish(u)

    k = 0
    do gu = 1, 2
      do block = 1, 3
        do m = 1, nr_mo(gu, block)
          k = k + 1
          if (k > range_mo(gu, block, 1)) range_mo(gu, block, 1) = k
          range_mo(gu, block, 2) = k
        end do
      end do
    end do

    gu_combine = .false.
    if (jbtof(irep, 1) == 2) then
!     ungerade perturbation
      gu_combine(1, 2) = .true.
      gu_combine(2, 1) = .true.
    else
!     gerade perturbation
      gu_combine(1, 1) = .true.
      gu_combine(2, 2) = .true.
    end if

    k = 0
    i = 0
    do gu1 = 1, 2
      do block1 = 1, 3
        do m1 = 1, nr_mo(gu1, block1)
          i = i + 1
          s = 0
          do gu2 = 1, 2
            do block2 = 1, 3
              if (gu_combine(gu1, gu2)         &
                  .and. (block1 == from_block) &
                  .and. (block2 == to_block)) then
                do m2 = 1, nr_mo(gu2, block2)
                  k = k + 1
                  s = s + 1
                  orbital_rotation_indices(1, k) = i
                  orbital_rotation_indices(2, k) = s
                end do
              else
                s = s + nr_mo(gu2, block2)
              end if
            end do
          end do
        end do
      end do
    end do

    if (debug_me) then
      call header('debug orbital_rotation_indices', -1)
      do k = 1, length
        write(*, *) k, orbital_rotation_indices(1, k), &
                 '->', orbital_rotation_indices(2, k)
      end do
    end if

  end subroutine

  subroutine response_solver(freq, RHS, Dp)

!   1. RHS (AO) -> RHS (MO)
!   2. RHS (MO) -> property gradient
!   3. solve response equation, in:  property gradient
!                               out: response vector
!   4. response vector -> Wp (MO)
!   5. Wp (MO)         -> Dp (AO)

!   ----------------------------------------------------------------------------
    real(8),      intent(in)    :: freq
    type(matrix), intent(inout) :: RHS
    type(matrix), intent(inout) :: Dp
!   ----------------------------------------------------------------------------
    type(matrix)                :: RHS_mo, Wp
    type(matrix)                :: C, Cig, Ciu, Csg, Csu

    integer,      allocatable   :: ibtyp(:)
    integer,      allocatable   :: ibtyp_pointer_pp(:)
    integer,      allocatable   :: ibtyp_pointer_pn(:)

    real(8),      allocatable   :: eigval(:)
    real(8),      allocatable   :: eigvec(:)
    real(8),      allocatable   :: convergence(:)
                  
    real(8),      allocatable   :: response_vector_pph(:, :, :)
    real(8),      allocatable   :: response_vector_ppa(:, :, :)
    real(8),      allocatable   :: response_vector_pnh(:, :, :)
    real(8),      allocatable   :: response_vector_pna(:, :, :)

    real(8),      pointer       :: work(:)
    integer                     :: lwork, kfree, lfree

    real(8),      allocatable   :: prop_gradient_pp(:, :)
    real(8),      allocatable   :: prop_gradient_pn(:, :)

    integer,      allocatable   :: orbital_rotation_indices_pp(:, :)
    integer,      allocatable   :: orbital_rotation_indices_pn(:, :)

    integer                     :: length_pp
    integer                     :: length_pn

    integer                     :: iz, k, i, s, off
    integer                     :: nr_freq

    logical                     :: debug_me = .false.
!   ----------------------------------------------------------------------------

!   fixme: nssh not always properly set, change to nesh minus nish

    if (jbtof(RHS%irep, 1) == 2) then
!     ungerade perturbation
      length_pp = nish(1)*nssh(2) + nish(2)*nssh(1)
      length_pn = nish(1)*npsh(2) + nish(2)*npsh(1)
    else
!     gerade perturbation
      length_pp = nish(1)*nssh(1) + nish(2)*nssh(2)
      length_pn = nish(1)*npsh(1) + nish(2)*npsh(2)
    end if

    nzconf  = 0
    nzxope  = 0
    nzxopp  = 0

    if (include_pp_rotations) nzxope = length_pp
    if (include_pn_rotations) nzxopp = length_pn

    nzxopt  = nzxope + nzxopp

    nzconfq = nzconf*nz
    nzxopeq = nzxope*nz
    nzxoppq = nzxopp*nz
    nzxoptq = nzxopt*nz

    nzvar   = nzxopt
    nzvarq  = nzxoptq

    call alloc(orbital_rotation_indices_pp, 2, length_pp)
    call set_orbital_rotation_indices(RHS%irep,  &
                                      length_pp, &
                                      2,         &
                                      3,         &
                                      orbital_rotation_indices_pp)

    call alloc(orbital_rotation_indices_pn, 2, length_pn)
    call set_orbital_rotation_indices(RHS%irep,  &
                                      length_pn, &
                                      2,         &
                                      1,         &
                                      orbital_rotation_indices_pn)

    call get_C(C, i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)
    RHS_mo = dag(C)*(RHS*C)
    C      = 0

    call copy_matrix_sym(RHS_mo, RHS)

    if (debug_me) then
      call print_mat(RHS,    label = 'debug RHS in AO basis')
      call print_mat(RHS_mo, label = 'debug RHS in MO basis')
    end if


!   get positive -> positive property gradient
!   ==========================================

    if (include_pp_rotations) then

      call alloc(prop_gradient_pp, length_pp, nz)
      do iz = 1, nz
        off = (iz - 1)*norbt*norbt
        do k = 1, length_pp
          i = orbital_rotation_indices_pp(1, k)
          s = orbital_rotation_indices_pp(2, k)
          prop_gradient_pp(k, iz) = prop_gradient_pp(k, iz) &
                                  - 2.0d0*RHS_mo%elms(off + (i - 1)*norbt + s)
        end do
      end do

      if (debug_me) then
        call header('debug prop_gradient_pp', -1)
        call print_q_vector(prop_gradient_pp, length_pp)
      end if

    end if


!   get positive -> negative property gradient
!   ==========================================

    if (include_pn_rotations) then

      call alloc(prop_gradient_pn, length_pn, nz)
      do iz = 1, nz
        off = (iz - 1)*norbt*norbt
        do k = 1, length_pn
          i = orbital_rotation_indices_pn(1, k)
          s = orbital_rotation_indices_pn(2, k)
          prop_gradient_pn(k, iz) = prop_gradient_pn(k, iz) &
                                  - 2.0d0*RHS_mo%elms(off + (i - 1)*norbt + s)
        end do
      end do

      if (debug_me) then
        call header('debug prop_gradient_pn', -1)
        call print_q_vector(prop_gradient_pn, length_pn)
      end if

    end if

!   inherit shape
    Wp     = RHS_mo
    RHS_mo = 0

!   fixme:
!   call with proper limit
!   at the moment i don't know the limit
!   simply give entire work and hope for the best
    lwork = len_f77_work
    call di_select_wrk(work, lwork)

    kfree = 1
    lfree = len_f77_work
    call memget('INTE', kzxope, 2*length_pp, work, kfree, lfree)
    call memget('INTE', kzxopp, 2*length_pn, work, kfree, lfree)

    call fill_work(include_pp_rotations,                      &
                   include_pn_rotations,                      &
                   length_pp,                  &
                   length_pn,                  &
                   orbital_rotation_indices_pp, &
                   orbital_rotation_indices_pn, &
                   work(kzxope),               &
                   work(kzxopp))


!   parameters
!   ==========

    call setrsp

    static    = .false.

    lineq     = .true.
!   in setrsp
!   fkrmc     = .false.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .false.
    diaghe    = .true.
    iprxrs    = 0
    thcxrs    = aor_cfg_threshold_response
    resfac    = 1.0d3
    maxitr    = 150
    nredm     = 400
    n2redm    = nredm*nredm
    loffty    = 0
    cnvint(1) = 1.0d20
    cnvint(2) = 1.0d20
    itrint(1) = 1
    itrint(2) = 1
!   in setrsp
!   uncoup    = .false.
    sternh    = .false.

    intdef    = integral_flag

!   nr_freq   = size(freq)
    nr_freq   = 1
    jsymop    = RHS%irep + 1
    jtimop    = RHS%tr_sym
    jopsy     = jbtof(RHS%irep, 1)
    nfreq     = nr_freq
    nexsim    = nr_freq
    nexstv    = nr_freq
    nexcnv    = nr_freq
    ncred     = 0
    nered     = 0
    npred     = 0
    nzred     = ncred + nered + npred


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call alloc(ibtyp,            2*nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)
    call alloc(convergence,      nr_freq)
    call alloc(eigval,           nr_freq)
    call alloc(eigvec,           nredm*nr_freq)

    call dcopy(nr_freq, freq, 1, eigval, 1)

    call xrsctl((/0.0d0/),        &
                prop_gradient_pp, &
                prop_gradient_pn, &
                ibtyp,            &
                (/0/),            &
                ibtyp_pointer_pp, &
                ibtyp_pointer_pn, &
                convergence,      &
                eigval,           &
                eigvec,           &
                work, kfree, lfree)

    call dealloc(convergence)
    call dealloc(eigval)
    call di_deselect_wrk(work, lwork)


!   construct response vector
!   =========================

    if (include_pp_rotations) then

      call alloc(response_vector_pph, length_pp, nz, nr_freq)
      call alloc(response_vector_ppa, length_pp, nz, nr_freq)

      response_vector_pph = 0.0d0
      response_vector_ppa = 0.0d0

      call construct_response_vector('pp',                &
                                     1,                   &
                                     response_vector_pph, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      call construct_response_vector('pp',                &
                                     -1,                  &
                                     response_vector_ppa, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pph', -1)
        call prbvec(6, response_vector_pph, nr_freq, length_pp)
        call header('debug response_vector_ppa', -1)
        call prbvec(6, response_vector_ppa, nr_freq, length_pp)
      end if

    end if

    if (include_pn_rotations) then

      call alloc(response_vector_pnh, length_pn, nz, nr_freq)
      call alloc(response_vector_pna, length_pn, nz, nr_freq)

      response_vector_pnh = 0.0d0
      response_vector_pna = 0.0d0

      call construct_response_vector('pn',                &
                                     1,                   &
                                     response_vector_pnh, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      call construct_response_vector('pn',                &
                                     -1,                  &
                                     response_vector_pna, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pnh', -1)
        call prbvec(6, response_vector_pnh, nr_freq, length_pn)
        call header('debug response_vector_pna', -1)
        call prbvec(6, response_vector_pna, nr_freq, length_pn)
      end if

    end if

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)
    call dealloc(eigvec)


!   scatter response vectors
!   ========================

    Wp%elms = 0.0d0

    if (include_pp_rotations) then
      call scatter_vector(length_pp,                   &
                          orbital_rotation_indices_pp, &
                          1.0d0,                       &
                          response_vector_pph,         &
                          Wp%elms,                     &
                          Wp%irep)
      call scatter_vector(length_pp,                   &
                          orbital_rotation_indices_pp, &
                         -1.0d0,                       &
                          response_vector_ppa,         &
                          Wp%elms,                     &
                          Wp%irep)

      call dealloc(orbital_rotation_indices_pp)
      call dealloc(response_vector_pph)
      call dealloc(response_vector_ppa)
    end if

    if (include_pn_rotations) then
      call scatter_vector(length_pn,                   &
                          orbital_rotation_indices_pn, &
                          1.0d0,                       &
                          response_vector_pnh,         &
                          Wp%elms,                     &
                          Wp%irep)
      call scatter_vector(length_pn,                   &
                          orbital_rotation_indices_pn, &
                         -1.0d0,                       &
                          response_vector_pna,         &
                          Wp%elms,                     &
                          Wp%irep)

      call dealloc(orbital_rotation_indices_pn)
      call dealloc(response_vector_pnh)
      call dealloc(response_vector_pna)
    end if


!   get coefficients
!   ================

    call get_C(Cig, i=1.0d0, s=0.0d0, g=1.0d0, u=0.0d0)
    call get_C(Csg, i=0.0d0, s=1.0d0, g=1.0d0, u=0.0d0)

    if (nfsym == 2) then
      call get_C(Ciu, i=1.0d0, s=0.0d0, g=0.0d0, u=1.0d0)
      call get_C(Csu, i=0.0d0, s=1.0d0, g=0.0d0, u=1.0d0)
    end if


!   construct perturbed AO density matrix
!   =====================================

    if (nfsym == 2) then
      if (jbtof(Wp%irep, 1) == 2) then
!       ungerade perturbation
        Dp = (Cig*(Wp*dag(Csu))) &
           + (Ciu*(Wp*dag(Csg))) &
           - (Csg*(Wp*dag(Ciu))) &
           - (Csu*(Wp*dag(Cig)))
      else
!       gerade perturbation
        Dp = (Cig*(Wp*dag(Csg))) &
           + (Ciu*(Wp*dag(Csu))) &
           - (Csg*(Wp*dag(Cig))) &
           - (Csu*(Wp*dag(Ciu)))
      end if
      Ciu = 0
      Csu = 0
    else
      Dp = (Cig*(Wp*dag(Csg))) &
         - (Csg*(Wp*dag(Cig)))
    end if
    Cig = 0
    Csg = 0

    call copy_matrix_sym(Dp, Wp)

    if (debug_me) then
      call print_mat(Wp, label = 'debug Wp')
      call print_mat(Dp, label = 'debug Dp')
    end if

    Wp = 0

  end subroutine

  subroutine construct_response_vector(vector_type,     &
                                       ih,              &
                                       response_vector, &
                                       length,          &
                                       nr_freq,         &
                                       ibtyp,           &
                                       ibtyp_pointer,   &
                                       eigvec)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: vector_type
    integer,      intent(in)  :: ih
    real(8),      intent(out) :: response_vector(*)
    integer,      intent(in)  :: length
    integer,      intent(in)  :: nr_freq
    integer,      intent(in)  :: ibtyp(*)
    integer,      intent(in)  :: ibtyp_pointer(*)
    real(8),      intent(in)  :: eigvec(*)
!   ----------------------------------------------------------------------------
    integer,      allocatable :: ivecs(:)
    real(8),      allocatable :: buffer(:)
    character(6)              :: file_name
    integer                   :: file_unit, nbtyp
!   ----------------------------------------------------------------------------

    select case (vector_type)
      case ('pp')
        nbtyp     = 1
        file_unit = luboe
        file_name = 'PAMBOE'
      case ('pn')
        nbtyp     = 2
        file_unit = lubop
        file_name = 'PAMBOP'
      case default
        call quit('error in construct_response_vector: unknown vector type')
    end select

    call alloc(ivecs,  nr_freq)
    call alloc(buffer, length*nz)

    open(file_unit,              &
         file   = file_name,     &
         form   = 'unformatted', &
         access = 'direct',      &
         recl   = 8*length*nz,   &
         status = 'old')

    call xrsxv1(ih,              &
                nbtyp,           &
                response_vector, &
                eigvec,          &
                nr_freq,         &
                ibtyp,           &
                ibtyp_pointer,   &
                ivecs,           &
                buffer)

    close(file_unit, status = 'keep')

    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine

  subroutine print_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(:)
    integer,       intent(in) :: vector_dimension

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                1,                &
                ipqtoq(1, 0),     &
                6)

  end subroutine

  subroutine print_q_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(*)
    integer,       intent(in) :: vector_dimension

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                nz,               &
                ipqtoq(1, 0),     &
                6)

  end subroutine

  subroutine add_op(op_name,    &
                    blocks,     &
                    typ,        &
                    tr_sym,     &
                    pg_sym,     &
                    prefactors, &
                    labels)

!   ----------------------------------------------------------------------------
    character(16), intent(in) :: op_name
    character(4),  intent(in) :: blocks
    integer,       intent(in) :: typ
    integer,       intent(in) :: tr_sym
    integer,       intent(in) :: pg_sym
    real(8),       intent(in) :: prefactors(:)
    character(8),  intent(in) :: labels(:)
!   ----------------------------------------------------------------------------
    integer                   :: i, ixpr, iprl
!   ----------------------------------------------------------------------------

    ixpr = nprps + 1
    iprl = nprplbl

    if (ixpr > maxprps) then
      call quit('add_op: ixpr > maxprps')
    end if

    prpnam(ixpr)  = op_name
    iprptyp(ixpr) = typ
    iprpsym(ixpr) = pg_sym
    iprptim(ixpr) = tr_sym

    do i = 1, size(prefactors)

      iprl = iprl + 1

      if (iprl > maxprplbl) then
        call quit('add_op: iprl > maxprplbl')
      end if

      prplbl(iprl)  = labels(i)
      pdoint(iprl)  = blocks
      iprlrep(iprl) = pg_sym - 1
      iprltyp(iprl) = tr_sym

      facprp(i, ixpr)  = prefactors(i)
      iprplbl(i, ixpr) = iprl

    end do

    nprps   = ixpr
    nprplbl = iprl

!   private
    op_name_to_aoproper_index(ixpr) = op_name

  end subroutine

  subroutine ask_for_integrals()

!   ----------------------------------------------------------------------------
    character(1) :: ss
    character(4) :: blocks
    character(8) :: s, s1, s2
    real(8)      :: f
    integer      :: icomponent, iatom
!   ----------------------------------------------------------------------------
!   todo: 1. adapt for symmetry
!            only diplen is right at the moment
!         2. generalize for other hamiltonians
!            not every combination hamiltonian--property is tested

    ss = '0'
    if (nonzero_small_metric) ss = '+'

!   add to dcbcls

!                               for the electric field
    call request_1el_integrals('DIPLEN ', '+00'//ss)
    call request_1el_integrals('DPLGRA ', '+00'//ss)
    call request_1el_integrals('QUAGRA ', '+00'//ss)
!                               for the magnetic field
    call request_1el_integrals('DIPLEN ', '0++0')
    call request_1el_integrals('DPLGRA ', '0++0')

    call request_1el_integrals('THETA  ', '+00'//ss)
    call request_1el_integrals('DSUSCGO', '+00'//ss)
    call request_1el_integrals('G1O    ', '+00'//ss)
    call request_1el_integrals('G1N    ', '+00'//ss)
    call request_1el_integrals('G1B    ', '000+')
    call request_1el_integrals('G1KX   ', '0++0')
    call request_1el_integrals('G1KY   ', '0++0')
    call request_1el_integrals('G1KZ   ', '0++0')

!   add to dcbxpr and dcbprl

    f = 1.0d0
    blocks = '+00'//ss
    call add_op('XDIPLEN         ', &
                blocks,             &
                1,                  &
                1,                  &
                (isymax(1, 1) + 1), &
                (/f/),              &
                (/'XDIPLEN '/))
    call add_op('YDIPLEN         ', &
                blocks,             &
                1,                  &
                1,                  &
                (isymax(2, 1) + 1), &
                (/f/),              &
                (/'YDIPLEN '/))
    call add_op('ZDIPLEN         ', &
                blocks,             &
                1,                  &
                1,                  &
                (isymax(3, 1) + 1), &
                (/f/),              &
                (/'ZDIPLEN '/))

    f = 1.0d0
    blocks = '+00'//ss
    do icomponent = 1, size(ij_iterator)
      s = ij_iterator(icomponent)//'THETA '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                  1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

    f = 0.5d0*cval
    blocks = '0++0'
    call add_op('XANGMOM         ', &
                blocks,             &
                5,                  &
               -1,                  &
                1,                  &
                (/f, f/),           &
                (/'ZDIPLEN ', 'YDIPLEN '/))
    call add_op('YANGMOM         ', &
                blocks,             &
                6,                  &
               -1,                  &
                1,                  &
                (/f, f/),           &
                (/'XDIPLEN ', 'ZDIPLEN '/))
    call add_op('ZANGMOM         ', &
                blocks,             &
                7,                  &
               -1,                  &
                1,                  &
                (/f, f/),           &
                (/'YDIPLEN ', 'XDIPLEN '/))

    f = -cval
    blocks = '0++0'
    do iatom = 1, nr_atoms*3
      s1 = prefix_zeros(iatom, 2)//' DPG Z'
      s2 = prefix_zeros(iatom, 2)//' DPG Y'
      call add_op(prefix_zeros(iatom, 3)//'AMDRX        ', &
                  blocks,                                  &
                  5,                                       &
                 -1,                                       &
                  1,                                       &
                  (/f,  f/),                               &
                  (/s1, s2/))
      s1 = prefix_zeros(iatom, 2)//' DPG X'
      s2 = prefix_zeros(iatom, 2)//' DPG Z'
      call add_op(prefix_zeros(iatom, 3)//'AMDRY        ', &
                  blocks,                                  &
                  6,                                       &
                 -1,                                       &
                  1,                                       &
                  (/f,  f/),                               &
                  (/s1, s2/))
      s1 = prefix_zeros(iatom, 2)//' DPG Y'
      s2 = prefix_zeros(iatom, 2)//' DPG X'
      call add_op(prefix_zeros(iatom, 3)//'AMDRZ        ', &
                  blocks,                                  &
                  7,                                       &
                 -1,                                       &
                  1,                                       &
                  (/f,  f/),                               &
                  (/s1, s2/))
    end do

    f = 1.0d0
    blocks = '+00'//ss
    do icomponent = 1, size(ij_iterator)
      s = ij_iterator(icomponent)//'SUSCGO'
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                  1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

    f = 1.0d0
    blocks = '+00'//ss
    do iatom = 1, nr_atoms*3
      do icomponent = 1, size(i_iterator)
        s = prefix_zeros(iatom, 2)//' DPG '//i_iterator(icomponent)
        call add_op(s//'        ',      &
                    blocks,             &
                    1,                  &
                    1,                  &
                    1,                  &
                    (/f/),              &
                    (/s/))
      end do
    end do

    f = 1.0d0
    blocks = '+00'//ss
    do iatom = 1, nr_atoms*3
      do icomponent = 1, size(ij_iterator)
        s = prefix_zeros(iatom, 2)//'QDG '//ij_iterator(icomponent)
        call add_op(s//'        ',      &
                    blocks,             &
                    1,                  &
                    1,                  &
                    1,                  &
                    (/f/),              &
                    (/s/))
      end do
    end do

    f = 1.0d0
    blocks = '+00'//ss
    do iatom = 1, nr_atoms*3
      s = 'G1O'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                  1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

    f = 1.0d0
    blocks = '+00'//ss
    do iatom = 1, nr_atoms*3
      s = 'G1N'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                  1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

    f = 2.0d0*(cval**2.0d0)
    blocks = '000+'
    do iatom = 1, nr_atoms*3
      s = 'G1B'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                  1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

    f = cval
    blocks = '0++0'
    do iatom = 1, nr_atoms*3
      s = 'G1KX'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                 -1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do
    do iatom = 1, nr_atoms*3
      s = 'G1KY'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                 -1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do
    do iatom = 1, nr_atoms*3
      s = 'G1KZ'//prefix_zeros(iatom, 3)//'  '
      call add_op(s//'        ',      &
                  blocks,             &
                  1,                  &
                 -1,                  &
                  1,                  &
                  (/f/),              &
                  (/s/))
    end do

  end subroutine

  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out) :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine

end module
