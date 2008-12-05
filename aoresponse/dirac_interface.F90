module dirac_interface_m

  use kinds_m
  use memory_allocator
  use matrix_oper
  use matrix_autop
  use aoresponse_input_m
  use aoresponse_xc_m

  implicit none

  interface di_get_gbds
    module procedure get_G
  end interface

  interface di_read_operator_int
    module procedure get_P
  end interface

contains

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

    

!     DO I=1,3
!        DO J=1,MXCENT
!           DOCOOR(I,J) = .TRUE.
!        END DO
!     END DO

    do i = 1, 3*mxcent
      dopert(i, 1) = .true.
      dopert(i, 2) = .true.
    end do

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

  subroutine mat_set_from_full(full, times, A, l)

    real(kind_dp), intent(in)           :: full(*)
    real(kind_dp), intent(in)           :: times
    type(matrix)                        :: A
    character(*),  intent(in), optional :: l
!   ----------------------------------------------------------------------------

    call dcopy(size(A%elms), full, 1, A%elms, 1)

    if (times .ne. 1.0d0) call dscal(size(A%elms), times, A%elms, 1)

  end subroutine

  subroutine mat_to_full(A, times, full, l)

    type(matrix),  intent(in)           :: A
    real(kind_dp), intent(in)           :: times
    real(kind_dp), intent(out)          :: full(*)
    character(*),  intent(in), optional :: l
!   ----------------------------------------------------------------------------

    call dcopy(size(A%elms), A%elms, 1, full, 1)

    if (times .ne. 1.0d0) call dscal(size(A%elms), times, full, 1)

  end subroutine

  function do_dft() result(l)

#include "dftcom.h"

    logical :: l

    l = dftadd

  end function

  subroutine get_F(F)

    type(matrix), intent(inout) :: F
!   ----------------------------------------------------------------------------
    type(matrix)                :: F1, F2
    integer,      parameter     :: file_unit = 70

#include "dcbbas.h"

    if (.not. F%is_defined) call mat_init(F, ntbas(0), ntbas(0))

    F1 = 1.0d-100*F
    F2 = 1.0d-100*F

    call opnfil(file_unit, 'DFFCK1', 'OLD', 'get_F')
    call reafck(file_unit, F1%elms, .true., 1)
    close(file_unit, status = 'keep')

    call opnfil(file_unit, 'DFFCK2', 'OLD', 'get_F')
    call reafck(file_unit, F2%elms, .true., 1)
    close(file_unit, status = 'keep')

    F = 2.0d0*F1 + 2.0d0*F2

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

    if (.not. G%is_defined) call mat_init(G, ntbas(0), ntbas(0))
    G = 1.0d0*D !inherit symmetry information
    G%elms = 0.0d0

    call stex_get_npos(npos_nr)
    call alloc(npos, npos_nr)

!   call alloc_maxsize(freemem, 1, verbose)
!   lfree = (freemem/2)*2
    lfree = n2bbasxq*20
    call alloc(work, lfree)
    kfree = 1

    call twofck(G%pg_sym, &
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

    type(matrix),  intent(inout) :: C
!   ----------------------------------------------------------------------------
    real(kind_dp), allocatable   :: eigenvalues(:)
    integer,       parameter     :: file_unit = 66
    real(kind_dp)                :: dummy_dp
    integer                      :: dummy_i

#include "dcborb.h"
#include "dcbbas.h"

    if (.not. C%is_defined) call mat_init(C, ntbas(0), norbt)
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

    if (.not. S%is_defined) call mat_init(S, ntbas(0), ntbas(0))
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

    if (.not. P%is_defined) call mat_init(P, ntbas(0), ntbas(0))
    P%elms = 0.0d0

    P%pg_sym        = iprpsym(op_index)
    P%tr_sym        = iprptim(op_index)
    P%fermi_pointer = jbtof(P%pg_sym - 1, 1)
    P%ih_sym        = P%tr_sym

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
      call print(P, 'P in AO basis')
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
      call header('orbital_rotation_vector_pp', -1)
      do i = 1, size(orbital_rotation_vector_pp)/2
        write(*, *) i, orbital_rotation_vector_pp(1, i), '->', orbital_rotation_vector_pp(2, i)
      end do
      call header('orbital_rotation_vector_pn', -1)
      do i = 1, size(orbital_rotation_vector_pn)/2
        write(*, *) i, orbital_rotation_vector_pn(1, i), '->', orbital_rotation_vector_pn(2, i)
      end do
    end if

  end subroutine

  subroutine response_solver(freq, RHS, Dp)

    real(kind_dp), intent(in)    :: freq(:)
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

    if (.not. aoresponse_skip_pp) do_pp = .true.
    if (.not. aoresponse_skip_pn) do_pn = .true.
    if (levyle .or. bss)          do_pn = .false.

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
    RHS_mo = dag(C)*RHS*C
    C      = 0

    if (debug_me) then
      call print(RHS,    'RHS in AO basis')
      call print(RHS_mo, 'RHS in MO basis')
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
        call header('prop_gradient_pp', -1)
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
        call header('prop_gradient_pn', -1)
        call print_q_vector(prop_gradient_pn, length_pn)
      end if

    end if

    Wp     = 1.0d-100*RHS_mo
    RHS_mo = 0

!   call alloc_maxsize(freemem, 1, verbose)
!   lfree = (freemem/2)*2
    lfree = n2bbasxq*100
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

    static    = .false.

    lineq     = .true.
    fkrmc     = .false.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .false.
    diaghe    = .true.
    iprxrs    = 0
    thcxrs    = aoresponse_threshold_response
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

    nr_freq   = size(freq)
    jsymop    = RHS%pg_sym
    jtimop    = RHS%tr_sym
    jopsy     = RHS%fermi_pointer
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
        call header('response_vector_pph', -1)
        call prbvec(6, response_vector_pph, nr_freq, length_pp)
        call header('response_vector_ppa', -1)
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
        call header('response_vector_pnh', -1)
        call prbvec(6, response_vector_pnh, nr_freq, length_pn)
        call header('response_vector_pna', -1)
        call prbvec(6, response_vector_pna, nr_freq, length_pn)
      end if

    end if

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)
    call dealloc(eigvec)


!   construct perturbed AO density matrix Dp
!   ========================================

    call get_i_or_s_C(Ci, 'i')
    call get_i_or_s_C(Cs, 's')

    Wp%elms = 0.0d0

    if (do_pp) then

      call scatter_vector(length_pp, orbital_rotation_vector_pp, 1.0d0, response_vector_pph, Wp%elms, Wp%pg_sym)
      call scatter_vector(length_pp, orbital_rotation_vector_pp,-1.0d0, response_vector_ppa, Wp%elms, Wp%pg_sym)

      call dealloc(orbital_rotation_vector_pp)

      call dealloc(response_vector_pph)
      call dealloc(response_vector_ppa)

    end if

    if (do_pn) then

      call scatter_vector(length_pn, orbital_rotation_vector_pn, 1.0d0, response_vector_pnh, Wp%elms, Wp%pg_sym)
      call scatter_vector(length_pn, orbital_rotation_vector_pn,-1.0d0, response_vector_pna, Wp%elms, Wp%pg_sym)

      call dealloc(orbital_rotation_vector_pn)

      call dealloc(response_vector_pnh)
      call dealloc(response_vector_pna)

    end if

    Dp      = RHS
    Dp%elms = 0.0d0
    if (.not. static) Dp%ih_sym = 0

    Dp = 0.5d0*(Ci*Wp*dag(Cs) - Cs*Wp*dag(Ci))

    if (debug_me) then
      call print(Wp, 'Wp')
      call print(Dp, 'Dp')
      call print(Ci, 'Ci')
      call print(Cs, 'Cs')
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

end module
