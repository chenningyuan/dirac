module dirac_interface

  use kinds_m
  use memory_allocator
  use matrix_defop
  use matrix_genop
  use aor_cfg

  implicit none

  interface di_get_gbds
    module procedure get_G
  end interface

  interface di_read_operator_int
    module procedure get_P
  end interface

  integer, parameter            :: max_nr_op = 137
  character(16),           save :: op_name_to_aoproper(max_nr_op)

contains

  subroutine print_prp_common_blocks()

!   ----------------------------------------------------------------------------
    integer :: i 
!   ----------------------------------------------------------------------------

#include "dcbxpr.h"
#include "dcbprl.h"

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

  end subroutine

  subroutine request_1el_integrals(label, blocks)

!   can be used like this:
!   call request_1el_integrals('BETAMAT', '000+')
!   ----------------------------------------------------------------------------
    character(7), intent(in) :: label
    character(4), intent(in) :: blocks
!   ----------------------------------------------------------------------------
    integer                  :: last_position
!   ----------------------------------------------------------------------------

#include "dcbcls.h"

    last_position = nprpcls + 1
    if (last_position > maxcls) then
      call quit('request_1el_integrals: last_position > maxcls')
    end if

    clscal(last_position) = .true.
    clsint(last_position) = label
    clscmb(last_position) = blocks

!   important only for c/s moment integrals
!   here set it to something
    iordcl(last_position) = 0

    nprpcls = last_position

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

    real(kind_dp)              :: work(*)
    integer                    :: lwork

    real(kind_dp)              :: g(*)
    integer                    :: g_dim

    logical                    :: geo, london

    integer                    :: max_derv, iatom_in

    logical                    :: energy, fock

    real(kind_dp)              :: dmat(*)
    integer                    :: nr_dmat
!   ----------------------------------------------------------------------------
    real(kind_dp), allocatable :: fmat(:)

    real(kind_dp)              :: dummy
    integer                    :: idummy

    logical                    :: gencnt

    integer                    :: i, i2typ, itype, ifctyp(8), iatom, fmat_dim

    real(kind_dp)              :: skip(3, 14, 3)
    real(kind_dp)              :: screen, f

#include "mxcent.h"
#include "cbihr2.h"
#include "dcbgrd.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "exeinf.h"
#include "dorps.h"

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

    allocate(fmat(fmat_dim))
    fmat = 0.0d0

!   LL integrals
!   ============

    i2typ    = 1
    intclass = 0

!   call dzero(gradee, mxcoor)
    call dzero(skip, 3*14*3)

    call twoint(work, lwork, fmat, dmat, 4*nr_dmat, idummy,       &
                ifctyp, dummy, idummy, idummy, itype,             &
                max_derv, iatom, .true., .true., .false., tktime, &
                1, 0, 0, 0, 0, .false., idummy, i2typ,            &
                icedif, screen, dummy, dummy, dummy, skip,        &
                relcal, gencnt, dummy, dummy)

    if (fock) then
       if (geo)                          f = 0.25d0
       if (london .and. max_derv .eq. 1) f = 0.50d0
       if (london .and. max_derv .eq. 2) f = 1.00d0

       call dscal(fmat_dim, f, fmat, 1)
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

!fixme: make proper interfaces
#include "mxcent.h"
#include "cbihr2.h"
#include "dcbgrd.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "exeinf.h"
#include "dorps.h"
#include "dcbbas.h"
#include "nuclei.h"
#include "dgroup.h"

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

    real(kind_dp)               :: skip(3, 14, 3)
    real(kind_dp)               :: screen

    integer                     :: i, j, k, i2typ, itype, ifctyp(8)

    real(kind_dp), allocatable  :: r(:, :)

    real(kind_dp), allocatable  :: work(:)
    integer                     :: lwork

    real(kind_dp)               :: dummy
    integer                     :: idummy
!   ----------------------------------------------------------------------------

#include "mxcent.h"
#include "cbihr2.h"
#include "dcbgrd.h"
#include "dcbgen.h"
#include "exeinf.h"
#include "dorps.h"

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

      write(*, '(a14, 3f21.10)') '2grad tot', r(j, 5), r(j + 1, 5), r(j + 2, 5)
      write(*, '(a14, 3f21.10)') '2grad  LL', r(j, 1), r(j + 1, 1), r(j + 2, 1)
      write(*, '(a14, 3f21.10)') '2grad  LS', r(j, 2), r(j + 1, 2), r(j + 2, 2)
      write(*, '(a14, 3f21.10)') '2grad  SS', r(j, 3), r(j + 1, 3), r(j + 2, 3)
      write(*, '(a14         )') '---------'

    end do

    call dealloc(r)

    D = 0.25d0*D

  end subroutine

  subroutine dipnuc_ifc(n, f, fg)

#include "mxcent.h"
#include "dipole.h"

    integer       :: n
    real(kind_dp) :: f(3), fg(3, 3*n)

    call dipnuc((/0.0d0/), (/0.0d0/), 0, .true.)

    f  = dipmn
    fg = ddipn(:, 1:3*n)

  end subroutine

  subroutine qdrnuc_ifc(q)

#include "mxcent.h"
#include "nuclei.h"
#include "priunit.h"

    real(kind_dp) :: q(6)

    call quit('implement qdrnuc_ifc first')
!   call nucqdr(cord(:, 1:nucdep), (/0.0d0, 0.0d0, 0.0d0/), lupri, 0)
!   q = (/qdrnuc(1,1), qdrnuc(1:2, 2), qdrnuc(1:3, 3)/)

  end subroutine

  subroutine gradnn_ifc(n, g)

#include "mxcent.h"
#include "../abacus/cbinuc.h"
#include "energy.h"

    integer       :: n
    real(kind_dp) :: g(3*n)

    iprint = 0
    maxdif = 1

    call nucrep((0.0d0), (0.0d0), (0.0d0))

    g(1:3*n) = gradnn(1:3*n)

  end subroutine

  function do_dft() result(l)

    use dft_cfg

!   ----------------------------------------------------------------------------
    logical :: l
!   ----------------------------------------------------------------------------

    l = dft_cfg_is_dft_calculation

  end function

  subroutine get_F(F)

    type(matrix), intent(inout) :: F
!   ----------------------------------------------------------------------------
    type(matrix)                :: F1, F2
    integer,      parameter     :: file_unit = 70

#include "dcbbas.h"

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
    F = 2.0d0*F

    F1 = 0
    F2 = 0

  end subroutine

  subroutine get_G(D, G)

    type(matrix),  intent(in)    :: D
    type(matrix),  intent(inout) :: G
!   ----------------------------------------------------------------------------
    integer                      :: npos_nr
    integer,       allocatable   :: npos(:)

    real(kind_dp), allocatable   :: work(:)
    integer                      :: kfree, lfree, freemem, verbose = 0

#include "dcbbas.h"
#include "dcbdhf.h"

    if (.not. mat_isdef(G)) call init_mat(G, ntbas(0), ntbas(0))
    G = 1.0d0*D !inherit symmetry information
    G%elms = 0.0d0

    call stex_get_npos(npos_nr)
    call alloc(npos, npos_nr)

!   call alloc_maxsize(freemem, 1, verbose)
!   lfree = (freemem/2)*2
    lfree = n2bbasxq*20
    call alloc(work, lfree)
    kfree = 1

    call twofck(G%irep + 1, &
!               G%ih_sym, &
                0,        &
                1,        &
                G%elms,   &
                D%elms,   &
                1,        &
                npos,     &
                intflg,   &
                0,        &
                work(kfree), lfree)

    call dealloc(work)
    call dealloc(npos)

    G = 4.0d0*G

  end subroutine

  subroutine get_C(C)

!   ----------------------------------------------------------------------------
    type(matrix),  intent(inout) :: C
!   ----------------------------------------------------------------------------
    real(kind_dp), allocatable   :: eigenvalues(:)
    integer,       parameter     :: file_unit = 66
    real(kind_dp)                :: dummy_dp
    integer                      :: dummy_i
!   ----------------------------------------------------------------------------

#include "dcborb.h"
#include "dcbbas.h"

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

#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"

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
    real(kind_dp), allocatable   :: work(:)
    integer                      :: kfree, lfree, kwork, lwork

#include "dcbbas.h"
#include "dcbham.h"

    kfree = 1

#include "memint.h"

    if (.not. mat_isdef(S)) call init_mat(S, ntbas(0), ntbas(0))
    S%elms = 0.0d0

    lfree = 2*((nnbbasx + n2bbasx)*2)/2
    call alloc(work, lfree)

    call gtovlt(S%elms, ssmtrc, work, lfree)

    call dealloc(work)

    S = 2.0d0*S

  end subroutine

  subroutine get_P(op_name, P)

    character(*),  intent(in)    :: op_name
    type(matrix),  intent(inout) :: P
!   ----------------------------------------------------------------------------
    real(kind_dp), allocatable   :: ptri(:), op1int(:)
    logical,       allocatable   :: first(:), lexpst(:)
    integer                      :: op_index, file_unit, idimension
    real(kind_dp)                :: dummy_dp
    logical                      :: file_exists, debug_me = .false.

    real(kind_dp), allocatable   :: work(:)
    integer                      :: kfree, lfree, freemem, verbose = 0

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

!   write(*, *) 'debug', op_name
    op_index = find_op_index(op_name)

    idimension = nnbbasx
    if (iprptim(op_index) == 0) then
      idimension = n2bbasx
    end if

    call alloc(ptri, nz*idimension)
    call alloc(op1int, idimension)
    allocate(first(nz))
    allocate(lexpst(nfsym))

    lexpst = .false.

    file_unit = lu1int

    call opnfil(file_unit, 'AOPROPER', 'UNKNOWN', 'get_P')

    if (.not. mat_isdef(P)) call init_mat(P, ntbas(0), ntbas(0))
    P%elms = 0.0d0

    P%irep   = iprpsym(op_index) - 1
    P%tr_sym = iprptim(op_index)
    P%ih_sym = P%tr_sym

!   call alloc_maxsize(freemem, 1, verbose)
!   lfree = (freemem/2)*2
    lfree = n2bbasxq*2
    call alloc(work, lfree)
    kfree = 1

    call prpex2(op_index,   &
                dummy_dp,   &
                .false.,    &
                lexpst,     &
                idimension, &
                ptri,       &
                P%elms,     &
                op1int,     &
                first,      &
                work, kfree, lfree, 0)

    call dealloc(work)
    close(file_unit, status = 'keep')

    P = 2.0d0*P

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

#include "dcborb.h"
#include "dgroup.h"

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

    real(kind_dp), intent(in)    :: freq
    type(matrix),  intent(inout) :: RHS
    type(matrix),  intent(inout) :: Dp
!   -------------------------------------------------------------------------
    type(matrix)                 :: RHS_mo, C, Ci, Cs, Wp

    integer,       allocatable   :: ibtyp(:)
    integer,       allocatable   :: ibtyp_pointer_pp(:)
    integer,       allocatable   :: ibtyp_pointer_pn(:)

    real(kind_dp), allocatable   :: eigval(:)
    real(kind_dp), allocatable   :: eigvec(:)
    real(kind_dp), allocatable   :: convergence(:)

    real(kind_dp), allocatable   :: response_vector_pph(:, :, :)
    real(kind_dp), allocatable   :: response_vector_ppa(:, :, :)
    real(kind_dp), allocatable   :: response_vector_pnh(:, :, :)
    real(kind_dp), allocatable   :: response_vector_pna(:, :, :)

    real(kind_dp), allocatable   :: work(:)
    integer                      :: kfree, lfree, freemem, verbose = 0

    real(kind_dp), allocatable   :: prop_gradient_pp(:)
    real(kind_dp), allocatable   :: prop_gradient_pn(:)

    integer,       allocatable   :: orbital_rotation_vector_pp(:)
    integer,       allocatable   :: orbital_rotation_vector_pn(:)

    integer                      :: length_pp
    integer                      :: length_pn

    integer                      :: nr_freq
    integer                      :: ill_int, isl_int, iss_int, igt_int

    logical                      :: debug_me = .false.
    logical                      :: do_pp    = .false.
    logical                      :: do_pn    = .false.

    real(kind_dp)                :: dummy_dp
    integer                      :: dummy_i

#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "dgroup.h"

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

!   call alloc_maxsize(freemem, 1, verbose)
!   lfree = (freemem/2)*2
    lfree = n2bbasxq*20
    call alloc(work, lfree)
    kfree = 1

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

!fixme:check call setrsp


    static    = .false.

    lineq     = .true.
    fkrmc     = .false.
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
    uncoup    = .false.
    sternh    = .false.

    ill_int   = iand(intgen,   1)
    isl_int   = iand(intgen/2, 1)
    iss_int   = iand(intgen/4, 1)
    igt_int   = iand(intgen/8, 1)

    intdef    = ill_int + 2*isl_int + 4*iss_int + 8*igt_int

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
    call dealloc(work)


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
!   Dp = (Ci*(dag(Wp)*dag(Cs))) + (Cs*(Wp*dag(Ci)))
    Dp = 0.5d0*Dp

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

    character(*),  intent(in)    :: vector_type
    integer,       intent(in)    :: ih
    real(kind_dp), intent(inout) :: response_vector(:, :, :)
    integer,       intent(in)    :: length
    integer,       intent(in)    :: nr_freq
    integer,       intent(in)    :: ibtyp(2, *)
    integer,       intent(in)    :: ibtyp_pointer(*)
    real(kind_dp), intent(in)    :: eigvec(*)
!   ----------------------------------------------------------------------------
    integer,       allocatable   :: ivecs(:)
    real(kind_dp), allocatable   :: buffer(:)
    character(6)                 :: file_name
    integer                      :: file_unit, nbtyp

#include "dcbxrs.h"
#include "dgroup.h"

    call alloc(ivecs, nr_freq)

    select case (vector_type)
      case ('pp')
        nbtyp     = 1
        file_unit = luboe
        file_name = 'PAMBOE'
      case ('pn')
        nbtyp     = -1
        file_unit = lubop
        file_name = 'PAMBOP'
      case default
        call quit('error in construct_response_vector: unknown vector type')
    end select

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

    call dscal(length*nz,  0.5d0, response_vector, 1)

    close(file_unit, status = 'keep')

    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine

  subroutine print_vector(vector, vector_dimension)

    real(kind_dp), intent(in) :: vector(:)
    integer,       intent(in) :: vector_dimension

#include "dgroup.h"

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

    real(kind_dp), intent(in) :: vector(:)
    integer,       intent(in) :: vector_dimension

#include "dgroup.h"

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                nz,               &
                ipqtoq(1, 0),     &
                6)

  end subroutine

  function find_op_index(op_name) result(op_index)

    character(*), intent(in) :: op_name
    integer                  :: op_index, i

    op_index = -1

    do i = 1, max_nr_op
      if (op_name_to_aoproper(i) == op_name) then
        op_index = i
      end if
    end do

    if (op_index == -1) then
      call quit('find_op_index: index for '//op_name//' not found')
    end if

  end function

  subroutine add_op_by_hand(op_name,    &
                            typ,        &
                            pg_sym,     &
                            tr_sym,     &
                            prefactors, &
                            pointers,   &
                            blocks)

!   ----------------------------------------------------------------------------
    character(16), intent(in)  :: op_name
    integer,       intent(in)  :: typ
    integer,       intent(in)  :: pg_sym
    integer,       intent(in)  :: tr_sym
    real(8),       intent(in)  :: prefactors(3)
    integer,       intent(in)  :: pointers(3)
    character(4),  intent(in)  :: blocks
!   ----------------------------------------------------------------------------
    integer                    :: i
!   ----------------------------------------------------------------------------

#include "dcbxpr.h"
#include "dcbprl.h"

    i = nprps + 1

    prpnam(i)  = op_name
    iprptyp(i) = typ
    iprpsym(i) = pg_sym
    iprptim(i) = tr_sym

    facprp(1:3,  i) = prefactors(1:3)
    iprplbl(1:3, i) = pointers(1:3)

    nprps = i

    op_name_to_aoproper(i) = op_name

    i = nprplbl + 1

    prplbl(i)  = op_name
    pdoint(i)  = blocks
    iprlrep(i) = pg_sym - 1
    iprltyp(i) = tr_sym

    nprplbl = i

  end subroutine

  subroutine add_op(op_name,       &
                    typ,           &
                    nr_components, &
                    prefactors,    &
                    labels)

!   ----------------------------------------------------------------------------
    character(7), intent(in)  :: op_name
    integer,      intent(in)  :: typ
    integer,      intent(in)  :: nr_components
    real(8),      intent(in)  :: prefactors(nr_components)
    character(8), intent(in)  :: labels(nr_components)
!   ----------------------------------------------------------------------------
    integer                   :: i
    integer                   :: pg_sym
    integer                   :: tr_sym
    real(8),      allocatable :: work(:)
    integer                   :: lwork
!   ----------------------------------------------------------------------------

    lwork = 200000
    call alloc(work, lwork)

    call xprind(op_name//'         ', &                                
                typ,                  &              
                nr_components,        &                           
                prefactors,           &                                
                labels,               &                                
                i,                    &            
                pg_sym,               &                 
                tr_sym,               &                     
                work,                 &               
                lwork,                &                
                0)

    call dealloc(work)

    op_name_to_aoproper(i) = op_name

  end subroutine

  subroutine ask_for_integrals()
      integer :: i
      real(8) :: v(3)

      do i = 1, max_nr_op
        op_name_to_aoproper(i) = ''
      end do

      call add_op('XDIPLEN', 1, 1, (/1.0d0/), (/'XDIPLEN '/))
      call add_op('YDIPLEN', 1, 1, (/1.0d0/), (/'YDIPLEN '/))
      call add_op('ZDIPLEN', 1, 1, (/1.0d0/), (/'ZDIPLEN '/))

!     call request_1el_integrals('BETAMAT', '000+')
      call request_1el_integrals('GEO1KIN', '+00+')

      v = (/1.0d0, 0.0d0, 0.0d0/)

      call add_op_by_hand('GEO1KIN001      ',1,1,1,v, &
                          (/ 9, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN002      ',1,1,1,v,  &
                          (/10, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN003      ',1,1,1,v,  &
                          (/11, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN004      ',1,1,1,v,  &
                          (/12, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN005      ',1,1,1,v,  &
                          (/13, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN006      ',1,1,1,v,  &
                          (/14, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN007      ',1,1,1,v,  &
                          (/15, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN008      ',1,1,1,v,  &
                          (/16, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN009      ',1,1,1,v,  &
                          (/17, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN010      ',1,1,1,v,  &
                          (/18, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN011      ',1,1,1,v,  &
                          (/19, 0, 0/), '+00+')
      call add_op_by_hand('GEO1KIN012      ',1,1,1,v,  &
                          (/20, 0, 0/), '+00+')

      call print_prp_common_blocks()
  end subroutine

end module
