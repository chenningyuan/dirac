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
  public get_i_or_s_c
  public debug_1grad
  public debug_2grad
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
  public init_cb_info

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

  private

  logical :: cb_info_allocated = .false.
  logical :: cb_info_set       = .false.
#ifdef VAR_MPI
  logical :: cb_info_synced    = .false.
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
  integer            :: len_f77_work  = 0 ! Length of the work array
  integer            :: memory_usage  = 0 ! Amount of work array used
  integer            :: mem_wrk_count = 0 ! Counts

contains

  subroutine init_cb_info()

!   ----------------------------------------------------------------------------
    integer :: ill_int
    integer :: ils_int
    integer :: iss_int
    integer :: igt_int
!   ----------------------------------------------------------------------------

    if (.not. cb_info_allocated) then

      allocate(sym_dep_nuc_name(mxcent))

      cb_info_allocated = .true.

    end if

!   if (.not. cb_info_set) then

      unit_input       = lucmd
      unit_output      = lupri

      nr_atoms         = natoms
      nr_sym_dep_nuc   = nucdep
      nr_sym_indep_nuc = nucind
      sym_dep_nuc_name = namdep

      if (mc == 2) then
        nonzero_small_metric = .true.
      else
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

      cb_info_set = .true.

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

  SUBROUTINE di_select_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req
!     if (mem_req .GT. len_f77_work - memory_usage) then
!    if (mem_req .GT. len_f77_work) then
!       call quit('code di_create_wrk first')
!       call di_create_wrk(wrk,mem_req)
!    else
!        wrk => f77_work((memory_usage+1):len_f77_work)
!        memory_usage = memory_usage + mem_req

        wrk => f77_work

        memory_usage  = mem_req
        mem_wrk_count = mem_wrk_count + 1
        if (mem_wrk_count > 1) call quit('counter > 1 not implemented')
!    endif
  end subroutine

  SUBROUTINE di_deselect_wrk(wrk, mem_req)
     real(8),pointer :: wrk(:)
     integer, intent(in) :: mem_req
!    if (mem_req .GT. len_f77_work) then
!       call di_delete_wrk(wrk)
!       call quit('code di_delete_wrk first')
!    else
        nullify(wrk)
!        memory_usage = memory_usage - mem_req
        memory_usage  = 0
        mem_wrk_count = mem_wrk_count - 1
!    endif
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
      write(*, *) 'dcbcls:  CLSCAL(i) ', CLSCAL(i)
      write(*, *) 'dcbcls:  CLSINT(i) ', CLSINT(i)
      write(*, *) 'dcbcls:  CLSCMB(i) ', CLSCMB(i)
      write(*, *) 'dcbcls:  IORDCL(i) ', IORDCL(i)
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






    call dcopy(g_dim, fmat, 1, g, 1)

    deallocate(fmat)

  end subroutine

  subroutine debug_1grad(D, F, nr_atoms)

!   ----------------------------------------------------------------------------
    type(matrix), intent(in) :: D
    type(matrix), intent(in) :: F
    integer                  :: nr_atoms
!   ----------------------------------------------------------------------------
    type(matrix)             :: C, D_folded, F_mo, F_temp, F_folded

    integer                  :: i, j, k, ifctyp(8), iz, ifsym, ioff

    real(8), allocatable     :: work(:)
    integer                  :: lwork

    real(8), allocatable     :: w01(:)
    real(8), allocatable     :: w02(:)
    real(8), allocatable     :: w03(:)
    real(8), allocatable     :: w04(:)
    integer, allocatable     :: w05(:)
    integer, allocatable     :: w06(:)
    integer, allocatable     :: w07(:)
    real(8), allocatable     :: w08(:)

    real(8)                  :: dummy
    integer                  :: idummy

    logical                  :: doint(2, 2) = .true.
    logical                  :: debug_me    = .false.
!   ----------------------------------------------------------------------------

    lwork = 20000000
    call alloc(work, lwork)

    D_folded = 4.0d0*D
    F_temp   = 0.5d0*F

    call bstobu(D_folded%elms, nz, work, lwork)
    call bstobu(F_temp%elms,   nz, work, lwork)

    call get_C(C)

    F_mo   = dag(C)*(F_temp*C)
    F_temp = 1.0d-50*F_mo

    do ifsym = 1, nfsym
      if (nish(ifsym) > 0) then
        ioff = i2orbt(ifsym) + npsh(ifsym)*norb(ifsym)
        call dcopy(nish(ifsym)*norb(ifsym), F_mo%elms(1 + ioff), 1, F_temp%elms(1 + ioff), 1)
      end if
    end do

    F_temp   = 2.0d0*F_temp
    F_folded = C*(F_temp*dag(C))

    C      = 0
    F_mo   = 0
    F_temp = 0

    do iz = 1, nz
      call rmolgrd_foldmat(D_folded%elms((iz - 1)*n2bbasx + 1), ntbas(0), (iz == 1))
      call rmolgrd_foldmat(F_folded%elms((iz - 1)*n2bbasx + 1), ntbas(0), (iz == 1))
    end do

!   fixme: half-phases ...
!   fixme: QO to AO

    if (debug_me) then
      call print_mat(F_mo,     label = 'debug F_mo')
      call print_mat(F_temp,   label = 'debug F_temp')
      call print_mat(F_folded, label = 'debug F_folded')
      call print_mat(D_folded, label = 'debug D_folded')
    end if

    grd_nosmlv = .false.
    grd_nosmls = .false.

    dopert = .true.
    ftronv = .true.
    relgrd = .true.

    call rsetsym(.false., relgrd, .false., .false.)
    call rtroini(.false.)

    call troinv(work, lwork)

!   set ifctyp
    ifctyp    = 0
    ifctyp(1) = 13 !coulomb and exchange
    ifctyp(2) = 22 !only exchange
    ifctyp(3) = 22
    ifctyp(4) = 22

    call alloc(w01, nnbbasx*3)
    call alloc(w02, nucdep*3)
    call alloc(w03, nucdep*3)
    call alloc(w04, nucdep*3)
    call alloc(w05, nucdep)
    call alloc(w06, nucdep)
    call alloc(w07, nucdep)
    call alloc(w08, nucdep)

    call onedr1(w01,                                                       &
                D_folded%elms, F_folded%elms,                              &
                w02, w03, w04, w05, w06, w07, w08,                         &
                dummy, work, lwork, 0, .true., 1, .true., .true., .false., &
                ntbas(0), nnbbasx, nnbbast, doint,                         &
                (/'FFFF', 'FFFF'/), .false.)

    D_folded = 0
    F_folded = 0

    call dealloc(work)

    call dealloc(w01)
    call dealloc(w02)
    call dealloc(w03)
    call dealloc(w04)
    call dealloc(w05)
    call dealloc(w06)
    call dealloc(w07)
    call dealloc(w08)

    do i = 1, nr_atoms
      j = (i - 1)*3 + 1
      write(*, '(a14, 3f21.10)') '1grad kin mo', (gradkn(j + k, 1), k = 0, 2)
      write(*, '(a14, 3f21.10)') '1grad kin be', (gradkn(j + k, 2), k = 0, 2)
      write(*, '(a14, 3f21.10)') '1grad nuc LL', (gradnu(j + k, 1), k = 0, 2)
      write(*, '(a14, 3f21.10)') '1grad nuc SS', (gradnu(j + k, 2), k = 0, 2)
      write(*, '(a14, 3f21.10)') '1grad reo LL', (gradro(j + k, 1), k = 0, 2)
      write(*, '(a14, 3f21.10)') '1grad reo SS', (gradro(j + k, 2), k = 0, 2)
      write(*, '(a14         )') '---------'
    end do

  end subroutine

  subroutine debug_2grad(D, nr_atoms)

!   ----------------------------------------------------------------------------
    type(matrix), intent(inout) :: D
    integer                     :: nr_atoms
!   ----------------------------------------------------------------------------
    logical                     :: debug_me = .false.

    real(8)               :: skip(3, 14, 3)
    real(8)               :: screen

    integer                     :: i, j, k, i2typ, itype, ifctyp(8)

    real(8), allocatable  :: r(:, :)

    real(8), allocatable  :: work(:)
    integer                     :: lwork

    real(8)               :: dummy
    integer                     :: idummy
!   ----------------------------------------------------------------------------

    D = 4.0d0*D
    if (debug_me) then
      call print_mat(D, label = 'debug D in AO basis')
    end if

    do i = 1, 3*mxcent
      dopert(i, 1) = .true.
      dopert(i, 2) = .true.
    end do

    ftronv = .true.
    relgrd = .true.
    call rsetsym(.false., relgrd, .false., .false.)
    call rtroini(.false.)

    lwork = 20000000
    call alloc(work, lwork)

    call troinv(work, lwork)

    itype  = 2
    screen = -1.0d0

!   set ifctyp
    ifctyp    = 0
    ifctyp(1) = 13 !coulomb and exchange
    ifctyp(2) = 22 !only exchange
    ifctyp(3) = 22
    ifctyp(4) = 22

    call alloc(r, nr_atoms*3, 5)
    r = 0.0d0


!   LL integrals
!   ============

    i2typ    = 1
    intclass = 0

    call dzero(gradee, mxcoor)
    call dzero(skip, 3*14*3)

    call twoint(work, lwork, dummy, D%elms, 4, idummy,     &
                ifctyp, dummy, idummy, idummy, itype,      &
                1, 0, .true., .true., .false., tktime,     &
                1, 0, 0, 0, 0, .false., idummy, i2typ,     &
                icedif, screen, dummy, dummy, dummy, skip, &
                relcal, .false., dummy, dummy)

    call dcopy(nr_atoms*3, gradee, 1, r(:, 1), 1)
    call daxpy(nr_atoms*3, 1.0d0, gradee, 1, r(:, 5), 1)


!   LS integrals
!   ============

    i2typ    = 2
    intclass = 1

    call dzero(gradee, mxcoor)
    call dzero(skip, 3*14*3)

    call twoint(work, lwork, dummy, D%elms, 4, idummy,     &
                ifctyp, dummy, idummy, idummy, itype,      &
                1, 0, .true., .true., .false., tktime,     &
                1, 0, 0, 0, 0, .false., idummy, i2typ,     &
                icedif, screen, dummy, dummy, dummy, skip, &
                relcal, .false., dummy, dummy)

    call dcopy(nr_atoms*3, gradee, 1, r(:, 2), 1)
    call daxpy(nr_atoms*3, 1.0d0, gradee, 1, r(:, 5), 1)


!   SS integrals
!   ============

    i2typ    = 3
    intclass = 2

    call dzero(gradee, mxcoor)
    call dzero(skip, 3*14*3)

    call twoint(work, lwork, dummy, D%elms, 4, idummy,     &
                ifctyp, dummy, idummy, idummy, itype,      &
                1, 0, .true., .true., .false., tktime,     &
                1, 0, 0, 0, 0, .false., idummy, i2typ,     &
                icedif, screen, dummy, dummy, dummy, skip, &
                relcal, .false., dummy, dummy)

    call dcopy(nr_atoms*3, gradee, 1, r(:, 3), 1)
    call daxpy(nr_atoms*3, 1.0d0, gradee, 1, r(:, 5), 1)

    call dealloc(work)

    do i = 1, 4
      j = (i - 1)*3 + 1

      write(*, '(a14, 3f27.16)') '2grad tot', r(j, 5), r(j + 1, 5), r(j + 2, 5)
      write(*, '(a14, 3f27.16)') '2grad  LL', r(j, 1), r(j + 1, 1), r(j + 2, 1)
      write(*, '(a14, 3f27.16)') '2grad  LS', r(j, 2), r(j + 1, 2), r(j + 2, 2)
      write(*, '(a14, 3f27.16)') '2grad  SS', r(j, 3), r(j + 1, 3), r(j + 2, 3)
      write(*, '(a14         )') '---------'

    end do

    call dealloc(r)

    D = 0.25d0*D

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

  subroutine get_C(C)

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

    call opnfil(file_unit, 'DFCOEF', 'OLD', 'get_C')

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

  subroutine get_i_or_s_C(C, kind)

    type(matrix), intent(inout) :: C
    character,    intent(in)    :: kind
!   ----------------------------------------------------------------------------
    type(matrix)                :: T
    integer                     :: pgs_nr,    egi_nr,    egs_nr
    integer                     :: pus_nr,    eui_nr,    eus_nr
    integer                     :: pgs_start, egi_start, egs_start
    integer                     :: pus_start, eui_start, eus_start
    integer                     :: i, s, iz
    logical                     :: debug_me = .false.

    pgs_nr = npsh(1)
    pus_nr = npsh(2)
    egi_nr = nish(1)
    eui_nr = nish(2)
    egs_nr = nssh(1)
    eus_nr = nssh(2)

    pgs_start = 1
    egi_start = pgs_start + pgs_nr
    egs_start = egi_start + egi_nr

    pus_start = pgs_nr    + pus_nr + egi_nr + 1
    eui_start = pus_start + pus_nr
    eus_start = eui_start + eui_nr

    call get_C(C)

    if (nfsym == 2) then

      select case (kind)

        case ('i')
          do iz = 1, nz
            do s = pgs_start, pgs_start + pgs_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
            do s = egs_start, egs_start + egs_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
            do s = pus_start, pus_start + pus_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
            do s = eus_start, eus_start + eus_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
          end do

        case ('s')
          do iz = 1, nz
            do i = egi_start, egi_start + egi_nr - 1
              call dzero(C%elms(1 + (i-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
            do i = eui_start, eui_start + eui_nr - 1
              call dzero(C%elms(1 + (i-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
          end do

      end select

    else

      select case (kind)

        case ('i')
          do iz = 1, nz
            do s = pgs_start, pgs_start + pgs_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
            do s = egs_start, egs_start + egs_nr - 1
              call dzero(C%elms(1 + (s-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
          end do

        case ('s')
          do iz = 1, nz
            do i = egi_start, egi_start + egi_nr - 1
              call dzero(C%elms(1 + (i-1)*ntbas(0) + (iz-1)*ntbas(0)*norbt), ntbas(0))
            end do
          end do

      end select

    end if

  end subroutine

  subroutine get_S(S)

    type(matrix),  intent(inout) :: S
!   ----------------------------------------------------------------------------
    real(8), pointer   :: work(:)
    integer            :: lwork

    if (.not. mat_isdef(S)) call init_mat(S, ntbas(0), ntbas(0))
    S%elms = 0.0d0

    lwork = 2*((nnbbasx + n2bbasx)*2)/2
    call di_select_wrk(work, lwork)

    call gtovlt(S%elms, ssmtrc, work, lwork)

    call di_deselect_wrk(work, lwork)

  end subroutine


  subroutine get_P(op_name, P)

    character(*),  intent(in)    :: op_name
    type(matrix),  intent(inout) :: P
!   ----------------------------------------------------------------------------
    real(8), allocatable   :: ptri(:), op1int(:)
    logical,       allocatable   :: first(:), lexpst(:)
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
    allocate(lexpst(nfsym))

    lexpst = .false.

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
    deallocate(lexpst)

  end subroutine

  subroutine set_orbital_rotation_vector(set,                        &
                                         length_pp,                  &
                                         length_pn,                  &
                                         orbital_rotation_vector_pp, &
                                         orbital_rotation_vector_pn)

    logical, intent(in)    :: set
    integer, intent(inout) :: length_pp
    integer, intent(inout) :: length_pn
    integer, intent(out)   :: orbital_rotation_vector_pp(2, length_pp)
    integer, intent(out)   :: orbital_rotation_vector_pn(2, length_pn)
!   ---------------------------------------------------------------------------
    integer                :: pgs_nr,    egi_nr,    egs_nr
    integer                :: pus_nr,    eui_nr,    eus_nr
    integer                :: pgs_start, egi_start, egs_start
    integer                :: pus_start, eui_start, eus_start
    integer                :: i, s
    logical                :: debug_me = .false.

    pgs_nr = npsh(1)
    pus_nr = npsh(2)
    egi_nr = nish(1)
    eui_nr = nish(2)
    egs_nr = nssh(1)
    eus_nr = nssh(2)

    pgs_start = 1
    egi_start = pgs_start + pgs_nr
    egs_start = egi_start + egi_nr

    pus_start = pgs_nr    + pus_nr + egi_nr + 1
    eui_start = pus_start + pus_nr
    eus_start = eui_start + eui_nr

    length_pp = 0
    length_pn = 0

    if (nfsym == 2) then

      do i = egi_start, egi_start + egi_nr - 1
        do s = eus_start, eus_start + eus_nr - 1
          length_pp = length_pp + 1
          if (set) then
            orbital_rotation_vector_pp(1, length_pp) = i
            orbital_rotation_vector_pp(2, length_pp) = s
          end if
        end do
        do s = pus_start, pus_start + pus_nr - 1
          length_pn = length_pn + 1
          if (set) then
            orbital_rotation_vector_pn(1, length_pn) = i
            orbital_rotation_vector_pn(2, length_pn) = s
          end if
        end do
      end do

      do i = eui_start, eui_start + eui_nr - 1
        do s = egs_start, egs_start + egs_nr - 1
          length_pp = length_pp + 1
          if (set) then
            orbital_rotation_vector_pp(1, length_pp) = i
            orbital_rotation_vector_pp(2, length_pp) = s
          end if
        end do
        do s = pgs_start, pgs_start + pgs_nr - 1
          length_pn = length_pn + 1
          if (set) then
            orbital_rotation_vector_pn(1, length_pn) = i
            orbital_rotation_vector_pn(2, length_pn) = s
          end if
        end do
      end do

    else

      do i = egi_start, egi_start + egi_nr - 1
        do s = egs_start, egs_start + egs_nr - 1
          length_pp = length_pp + 1
          if (set) then
            orbital_rotation_vector_pp(1, length_pp) = i
            orbital_rotation_vector_pp(2, length_pp) = s
          end if
        end do
        do s = pgs_start, pgs_start + pgs_nr - 1
          length_pn = length_pn + 1
          if (set) then
            orbital_rotation_vector_pn(1, length_pn) = i
            orbital_rotation_vector_pn(2, length_pn) = s
          end if
        end do
      end do

    end if

    if (debug_me) then
      call header('debug orbital_rotation_vector_pp', -1)
      do i = 1, size(orbital_rotation_vector_pp)/2
        write(*, *) i, orbital_rotation_vector_pp(1, i), '->', orbital_rotation_vector_pp(2, i)
      end do
      call header('debug orbital_rotation_vector_pn', -1)
      do i = 1, size(orbital_rotation_vector_pn)/2
        write(*, *) i, orbital_rotation_vector_pn(1, i), '->', orbital_rotation_vector_pn(2, i)
      end do
    end if

  end subroutine

  subroutine response_solver(freq, RHS, Dp)

    real(8), intent(in)    :: freq
    type(matrix),  intent(inout) :: RHS
    type(matrix),  intent(inout) :: Dp
!   -------------------------------------------------------------------------
    type(matrix)                 :: RHS_mo, C, Ci, Cs, Wp

    integer,       allocatable   :: ibtyp(:)
    integer,       allocatable   :: ibtyp_pointer_pp(:)
    integer,       allocatable   :: ibtyp_pointer_pn(:)

    real(8), allocatable   :: eigval(:)
    real(8), allocatable   :: eigvec(:)
    real(8), allocatable   :: convergence(:)

    real(8), allocatable   :: response_vector_pph(:, :, :)
    real(8), allocatable   :: response_vector_ppa(:, :, :)
    real(8), allocatable   :: response_vector_pnh(:, :, :)
    real(8), allocatable   :: response_vector_pna(:, :, :)

    real(8), pointer   :: work(:)
    integer            :: lwork
    integer            :: kfree, lfree

    real(8), allocatable   :: prop_gradient_pp(:)
    real(8), allocatable   :: prop_gradient_pn(:)

    integer,       allocatable   :: orbital_rotation_vector_pp(:)
    integer,       allocatable   :: orbital_rotation_vector_pn(:)

    integer                      :: length_pp
    integer                      :: length_pn

    integer                      :: nr_freq

    logical                      :: debug_me = .false.
    logical                      :: do_pp    = .false.
    logical                      :: do_pn    = .false.

    real(8)                :: dummy_dp
    integer                      :: dummy_i

    if (.not. aor_cfg_skip_pp) do_pp = .true.
    if (.not. aor_cfg_skip_pn) do_pn = .true.
    if (levyle .or. bss)       do_pn = .false.

    call set_orbital_rotation_vector(.false.,                    &
                                     length_pp,                  &
                                     length_pn,                  &
                                     orbital_rotation_vector_pp, &
                                     orbital_rotation_vector_pn)

    nzconf  = 0
    nzxope  = 0
    nzxopp  = 0

    if (do_pp) nzxope = length_pp
    if (do_pn) nzxopp = length_pn

    nzxopt  = nzxope + nzxopp

    nzconfq = nzconf*nz
    nzxopeq = nzxope*nz
    nzxoppq = nzxopp*nz
    nzxoptq = nzxopt*nz

    nzvar   = nzxopt
    nzvarq  = nzxoptq

    call alloc(orbital_rotation_vector_pp, 2*length_pp)
    call alloc(orbital_rotation_vector_pn, 2*length_pn)

    call set_orbital_rotation_vector(.true.,                     &
                                     length_pp,                  &
                                     length_pn,                  &
                                     orbital_rotation_vector_pp, &
                                     orbital_rotation_vector_pn)

    call get_C(C)
    RHS_mo = dag(C)*(RHS*C)
    C      = 0

!fixme
    RHS_mo%ih_sym = RHS%ih_sym
    RHS_mo%tr_sym = RHS%tr_sym
    RHS_mo%irep   = RHS%irep

    if (debug_me) then
      call print_mat(RHS,    label = 'debug RHS in AO basis')
      call print_mat(RHS_mo, label = 'debug RHS in MO basis')
    end if

    if (do_pp) then

      call alloc(prop_gradient_pp, length_pp*nz)
      prop_gradient_pp = 0.0d0

      call prporb(prop_gradient_pp,           &
                  orbital_rotation_vector_pp, &
                  length_pp,                  &
                  RHS_mo%elms,                &
                  norbt,                      &
                  dummy_dp,                   &
                  nz,                         &
                  RHS_mo%ih_sym,              &
                  0)

      call dscal(length_pp*nz, -1.0d0, prop_gradient_pp, 1)

      if (debug_me) then
        call header('debug prop_gradient_pp', -1)
        call print_q_vector(prop_gradient_pp, length_pp)
      end if

    end if

    if (do_pn) then

      call alloc(prop_gradient_pn, length_pn*nz)
      prop_gradient_pn = 0.0d0

      call prporb(prop_gradient_pn,           &
                  orbital_rotation_vector_pn, &
                  length_pn,                  &
                  RHS_mo%elms,                &
                  norbt,                      &
                  dummy_dp,                   &
                  nz,                         &
                  RHS_mo%ih_sym,              &
                  0)

      call dscal(length_pn*nz, -1.0d0, prop_gradient_pn, 1)

      if (debug_me) then
        call header('debug prop_gradient_pn', -1)
        call print_q_vector(prop_gradient_pn, length_pn)
      end if

    end if

    Wp     = 1.0d-100*RHS_mo
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

    call fill_work(do_pp,                      &
                   do_pn,                      &
                   length_pp,                  &
                   length_pn,                  &
                   orbital_rotation_vector_pp, &
                   orbital_rotation_vector_pn, &
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
!   jopsy     = RHS%fermi_pointer
    jopsy     = 1
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

    call xrsctl(dummy_dp,         &
                prop_gradient_pp, &
                prop_gradient_pn, &
                ibtyp,            &
                dummy_i,          &
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

    if (do_pp) then

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

    if (do_pn) then

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


!   construct perturbed AO density matrix Dp
!   ========================================

    Wp%elms = 0.0d0

    if (do_pp) then

      call scatter_vector(length_pp, orbital_rotation_vector_pp, 1.0d0, response_vector_pph, Wp%elms, Wp%irep + 1)
      call scatter_vector(length_pp, orbital_rotation_vector_pp,-1.0d0, response_vector_ppa, Wp%elms, Wp%irep + 1)

      call dealloc(orbital_rotation_vector_pp)

      call dealloc(response_vector_pph)
      call dealloc(response_vector_ppa)

    end if

    if (do_pn) then

      call scatter_vector(length_pn, orbital_rotation_vector_pn, 1.0d0, response_vector_pnh, Wp%elms, Wp%irep + 1)
      call scatter_vector(length_pn, orbital_rotation_vector_pn,-1.0d0, response_vector_pna, Wp%elms, Wp%irep + 1)

      call dealloc(orbital_rotation_vector_pn)

      call dealloc(response_vector_pnh)
      call dealloc(response_vector_pna)

    end if

    Dp      = RHS
    Dp%elms = 0.0d0
    if (.not. static) Dp%ih_sym = 0

    call get_i_or_s_C(Ci, 'i')
    call get_i_or_s_C(Cs, 's')

    Dp = (Ci*(Wp*dag(Cs))) - (Cs*(Wp*dag(Ci)))

!fixme
    Dp%ih_sym = Wp%ih_sym
    Dp%tr_sym = Wp%tr_sym
    Dp%irep   = Wp%irep

    if (debug_me) then
      call print_mat(Wp, label = 'debug Wp')
      call print_mat(Dp, label = 'debug Dp')
      call print_mat(Ci, label = 'debug Ci')
      call print_mat(Cs, label = 'debug Cs')
    end if

    Ci = 0
    Cs = 0
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

    prpnam(ixpr)  = op_name
    iprptyp(ixpr) = typ
    iprpsym(ixpr) = pg_sym
    iprptim(ixpr) = tr_sym

    do i = 1, size(prefactors)

      iprl = iprl + 1

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

end module
