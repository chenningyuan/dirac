module dirac_interface_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use aoresponse_input_m

  implicit none

! private

contains

  subroutine get_D_mo(D_mo)

    type(matrix), intent(inout) :: D_mo

#include "dcborb.h"
#include "dgroup.h"

    integer(kind_i4)            :: i, j, offset

    if (.not. isdef(D_mo)) call mat_init(D_mo, norbt, norbt)

    do i = 1, nfsym
      if(nish(i) > 0) then
        do j = 1, nish(i)

          offset = (norb(1)*(i-1) + npsh(i))*D_mo%nrow
          offset = offset + norb(1)*(i-1) + npsh(i)
          offset = offset + (j-1)*D_mo%nrow + j - 1

          D_mo%elements(1 + offset) = 1.0d0

        end do
      end if
    end do

  end subroutine get_D_mo

  subroutine get_F(F)

    type(matrix), intent(inout) :: F

#include "dcbbas.h"

    type(matrix)                :: F1, F2
    integer, parameter          :: file_unit = 70

    if (.not. isdef(F)) call mat_init(F, ntbas(0), ntbas(0))

    F1 = 1.0d-100*F
    F2 = 1.0d-100*F

    call opnfil(file_unit, 'DFFCK1', 'OLD', 'get_F')
    call reafck(file_unit, F1%elements, .true., 1)
    close(file_unit, status = 'keep')

    call opnfil(file_unit, 'DFFCK2', 'OLD', 'get_F')
    call reafck(file_unit, F2%elements, .true., 1)
    close(file_unit, status = 'keep')

    F = 2.0d0*F1 + F2

    F1 = 0
    F2 = 0

  end subroutine get_F

  subroutine get_G(D, G)

    type(matrix), intent(in)      :: D
    type(matrix), intent(inout)   :: G

#include "dcbbas.h"
#include "dcbdhf.h"

    real(kind_dp),    allocatable :: work(:)
    integer(kind_i4)              :: kfree, lfree

    integer(kind_i4)              :: npos_nr
    integer(kind_i4), allocatable :: npos(:)

    call stex_get_npos(npos_nr)

    call alloc(npos, npos_nr)

!   if (.not. isdef(G)) call mat_init(G, ntbas(0), ntbas(0))
    G = 1.0d-100*D

    kfree = 1
    lfree = n2bbasxq*20
    call alloc(work, lfree)

    call twofck(G%pg_sym,   &
                G%ih_sym,   &
                1,          &
                G%elements, &
                D%elements, &
                1,          &
                npos,       &
                intflg,     &
                0,          &
                work(kfree), lfree)

    call dealloc(npos)
    call dealloc(work)

  end subroutine get_G

  subroutine get_C(C)

    type(matrix), intent(inout) :: C

#include "dcborb.h"
#include "dcbbas.h"

    real(kind_dp), allocatable :: eigenvalues(:)
    integer, parameter         :: file_unit = 66
    real(kind_dp)              :: dummy_dp
    integer(kind_i4)           :: dummy_i4

    if (.not. isdef(C)) call mat_init(C, ntbas(0), norbt)

    call alloc(eigenvalues, norbt)

    call opnfil(file_unit, 'DFCOEF', 'OLD', 'get_C')

    call reacmo(file_unit,   &
                'DFCOEF',    &
                C%elements,  &
                eigenvalues, &
                dummy_i4,    &
                dummy_dp,    &
                6,           &
                dummy_dp,    &
                dummy_i4,    &
                dummy_i4)

    close(file_unit, status = 'keep')

    call dealloc(eigenvalues)

  end subroutine get_C

  subroutine get_S(S)

    type(matrix),  intent(inout) :: S

#include "dcbbas.h"

    real(kind_dp), allocatable   :: work(:)
    integer(kind_i4)             :: kfree, lfree, kwork, lwork

    kfree = 1

#include "memint.h"

    if (.not. isdef(S)) call mat_init(S, ntbas(0), ntbas(0))

    lfree = nnbbasx + n2bbasx
!fixme: factor 2 shouldn't be needed
    lfree = lfree*2
    call alloc(work, lfree)

    call gtovlt(S%elements, 1.0d0, work, lfree)

    call dealloc(work)

    S = 2.0d0*S

  end subroutine get_S

  subroutine get_P(op_name, P)

    character(*), intent(in)    :: op_name
    type(matrix), intent(inout) :: P

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    real(kind_dp), allocatable :: ptri(:)
    real(kind_dp), allocatable :: op1int(:)
    logical,       allocatable :: first(:)

    real(kind_dp), allocatable :: work(:)
    integer(kind_i4)           :: kfree, lfree

    integer(kind_i4)           :: op_index
    logical                    :: file_exists, is_square
    real(kind_dp)              :: dummy_dp
    integer(kind_i4)           :: file_unit, dummy_i4, idimension

    is_square = .false.
    if (iprptim(find_op_index(op_name)) == 0) then
      is_square = .true.
    end if

    if (is_square) then
      idimension = n2bbasx
    else
      idimension = nnbbasx
    end if

    call alloc(ptri,   nz*idimension)
    call alloc(op1int, idimension)
    allocate(first(nz))

    kfree = 1
    lfree = 2*nz*n2bbasx
    call alloc(work, lfree)

    file_unit = lu1int

    call opnfil(file_unit, 'AOPROPER', 'UNKNOWN', 'get_P')

    if (.not. isdef(P)) call mat_init(P, ntbas(0), ntbas(0))

    op_index = find_op_index(op_name)

    call prpex2(op_index,   &
                dummy_dp,   &
                .false.,    &
                idimension, &
                ptri,       &
                P%elements, &
                op1int,     &
                first,      &
                work, kfree, lfree, 0)

    P%pg_sym        = iprpsym(op_index)
    P%tr_sym        = iprptim(op_index)
    P%fermi_pointer = jbtof(P%pg_sym - 1, 1)
    P%ih_sym        = P%tr_sym

    close(file_unit, status = 'keep')

    call dealloc(ptri)
    call dealloc(op1int)
    deallocate(first)
    call dealloc(work)

  end subroutine get_P

  subroutine set_orbital_rotation_vector(orbital_rotation_vector_pp, orbital_rotation_vector_pn, &
                                         length_pp, length_pn)

    integer         , intent(out) :: orbital_rotation_vector_pp(2, length_pp)
    integer         , intent(out) :: orbital_rotation_vector_pn(2, length_pn)
    integer(kind_i4), intent(in)  :: length_pp
    integer(kind_i4), intent(in)  :: length_pn

#include "dcborb.h"
#include "dgroup.h"

    integer(kind_i4)              :: counter_pp
    integer(kind_i4)              :: counter_pn

    integer(kind_i4)              :: pg_nr,    egi_nr,    egs_nr
    integer(kind_i4)              :: pu_nr,    eui_nr,    eus_nr
    integer(kind_i4)              :: pg_start, egi_start, egs_start
    integer(kind_i4)              :: pu_start, eui_start, eus_start
    integer(kind_i4)              :: i, a

    logical                       :: debug_me = .false.

    pg_nr  = npsh(1)
    pu_nr  = npsh(2)
    egi_nr = nish(1)
    eui_nr = nish(2)
    egs_nr = nesh(1) - nish(1)
    eus_nr = nesh(2) - nish(2)

    pg_start  = 1
    egi_start = pg_start  + pg_nr
    egs_start = egi_start + egi_nr

    pu_start  = pg_nr     + pu_nr  + egi_nr + 1
    eui_start = pu_start  + pu_nr
    eus_start = eui_start + eui_nr

    counter_pp = 0
    counter_pn = 0

    if (nfsym == 2) then

      do i = egi_start, egi_start + egi_nr - 1
        do a = eus_start, eus_start + eus_nr - 1
          counter_pp = counter_pp + 1
          orbital_rotation_vector_pp(1, counter_pp) = i
          orbital_rotation_vector_pp(2, counter_pp) = a
        end do
      end do
      do i = eui_start, eui_start + eui_nr - 1
        do a = egs_start, egs_start + egs_nr - 1
          counter_pp = counter_pp + 1
          orbital_rotation_vector_pp(1, counter_pp) = i
          orbital_rotation_vector_pp(2, counter_pp) = a
        end do
      end do
      do i = egi_start, egi_start + egi_nr - 1
        do a = pu_start, pu_start + pu_nr - 1
          counter_pn = counter_pn + 1
          orbital_rotation_vector_pn(1, counter_pn) = i
          orbital_rotation_vector_pn(2, counter_pn) = a
        end do
      end do
      do i = eui_start, eui_start + eui_nr - 1
        do a = pg_start, pg_start + pg_nr - 1
          counter_pn = counter_pn + 1
          orbital_rotation_vector_pn(1, counter_pn) = i
          orbital_rotation_vector_pn(2, counter_pn) = a
        end do
      end do

    else

      do i = egi_start, egi_start + egi_nr - 1
        do a = egs_start, egs_start + egs_nr - 1
          counter_pp = counter_pp + 1
          orbital_rotation_vector_pp(1, counter_pp) = i
          orbital_rotation_vector_pp(2, counter_pp) = a
        end do
      end do
      do i = egi_start, egi_start + egi_nr - 1
        do a = pg_start, pg_start + pg_nr - 1
          counter_pn = counter_pn + 1
          orbital_rotation_vector_pn(1, counter_pn) = i
          orbital_rotation_vector_pn(2, counter_pn) = a
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

  end subroutine set_orbital_rotation_vector

  subroutine get_length(length_pp, length_pn)

    integer(kind_i4), intent(out) :: length_pp, length_pn

#include "dcborb.h"
#include "dgroup.h"

    integer(kind_i4)              :: pg_nr,    egi_nr,    egs_nr
    integer(kind_i4)              :: pu_nr,    eui_nr,    eus_nr
    integer(kind_i4)              :: pg_start, egi_start, egs_start
    integer(kind_i4)              :: pu_start, eui_start, eus_start
    integer(kind_i4)              :: i, a

    pg_nr  = npsh(1)
    pu_nr  = npsh(2)
    egi_nr = nish(1)
    eui_nr = nish(2)
    egs_nr = nesh(1) - nish(1)
    eus_nr = nesh(2) - nish(2)

    pg_start  = 1
    egi_start = pg_start  + pg_nr
    egs_start = egi_start + egi_nr

    pu_start  = pg_nr     + pu_nr  + egi_nr + 1
    eui_start = pu_start  + pu_nr
    eus_start = eui_start + eui_nr

    length_pp = 0
    length_pn = 0

    if (nfsym == 2) then

      do i = egi_start, egi_start + egi_nr - 1
        do a = eus_start, eus_start + eus_nr - 1
          length_pp = length_pp + 1
        end do
      end do
      do i = eui_start, eui_start + eui_nr - 1
        do a = egs_start, egs_start + egs_nr - 1
          length_pp = length_pp + 1
        end do
      end do
      do i = egi_start, egi_start + egi_nr - 1
        do a = pu_start, pu_start + pu_nr - 1
          length_pn = length_pn + 1
        end do
      end do
      do i = eui_start, eui_start + eui_nr - 1
        do a = pg_start, pg_start + pg_nr - 1
          length_pn = length_pn + 1
        end do
      end do

    else

      do i = egi_start, egi_start + egi_nr - 1
        do a = egs_start, egs_start + egs_nr - 1
          length_pp = length_pp + 1
        end do
      end do
      do i = egi_start, egi_start + egi_nr - 1
        do a = pg_start, pg_start + pg_nr - 1
          length_pn = length_pn + 1
        end do
      end do

    end if

  end subroutine get_length

  subroutine response_solver(freq, RHS, D1)

    real(kind_dp),    intent(in)    :: freq(:)
    type(matrix),     intent(in)    :: RHS
    type(matrix),     intent(inout) :: D1

#include "dcbbas.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "dgroup.h"

    integer(kind_i4), allocatable   :: ibtyp(:)
    integer(kind_i4), allocatable   :: ibtyp_pointer_c(:)
    integer(kind_i4), allocatable   :: ibtyp_pointer_pp(:)
    integer(kind_i4), allocatable   :: ibtyp_pointer_pn(:)
                                    
    integer(kind_i4), allocatable   :: isymop(:)
    integer(kind_i4), allocatable   :: ihrmop(:)
    integer(kind_i4), allocatable   :: ifckop(:)
                                    
    real(kind_dp),    allocatable   :: eigval(:)
    real(kind_dp),    allocatable   :: eigvec(:)
    real(kind_dp),    allocatable   :: convergence(:)
                                    
    real(kind_dp),    allocatable   :: response_vector_pph(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_ppa(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_pnh(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_pna(:, :, :)
                                    
    real(kind_dp),    allocatable   :: temp_ao_dmat(:, :, :, :)
                                    
    real(kind_dp),    allocatable   :: work(:)
    integer(kind_i4)                :: kfree, lfree
                                    
    real(kind_dp),    allocatable   :: prop_gradient_pp(:)
    real(kind_dp),    allocatable   :: prop_gradient_pn(:)
                                    
    integer(kind_i4), allocatable   :: orbital_rotation_vector_pp(:, :)
    integer(kind_i4), allocatable   :: orbital_rotation_vector_pn(:, :)
                                    
    integer(kind_i4)                :: length_pp, length_pn
    integer(kind_i4)                :: i, nr_freq, nr_dmat
    integer(kind_i4)                :: ill_int, isl_int, iss_int, igt_int
    logical                         :: debug_me = .false.
    real(kind_dp)                   :: dummy_dp
    type(matrix)                    :: M, C

!!!!!!!!!!
    kfree = 1
    lfree = 100*ntbas(0)*ntbas(0) ! no idea how much is needed
    call alloc(work, lfree)
!!!!!!!!!!

    call get_length(length_pp, length_pn)

    nzconf  = 0
    nzconfq = 0
    nzxope  = length_pp
    nzxopeq = length_pp*nz
    nzxopp  = length_pn
    nzxoppq = length_pn*nz
    nzxopt  = length_pp + length_pn
    nzxoptq = nzxopt*nz
    nzvar   = nzxopt
    nzvarq  = nzxoptq

    call alloc(prop_gradient_pp, length_pp*nz)
    call alloc(prop_gradient_pn, length_pn*nz)

    allocate(orbital_rotation_vector_pp(2, length_pp))
    allocate(orbital_rotation_vector_pn(2, length_pn))

    call set_orbital_rotation_vector(orbital_rotation_vector_pp, orbital_rotation_vector_pn, &
                                     length_pp, length_pn)

    call get_C(C)

    M = 1.0d-100*RHS

    call qrfaomo(M%elements,        &
                 RHS%elements,      &
                 RHS%pg_sym,        &
                 RHS%tr_sym,        &
                 RHS%fermi_pointer, &
                 0,                 &
                 C%elements,        &
                 work(kfree),       &
                 lfree)

    call prporb(prop_gradient_pp,           &
                orbital_rotation_vector_pp, &
                length_pp,                  &
                M%elements,                 &
                norbt,                      &
                dummy_dp,                   &
                nz,                         &
                M%ih_sym,                   &
                0)

    call prporb(prop_gradient_pn,           &
                orbital_rotation_vector_pn, &
                length_pn,                  &
                M%elements,                 &
                norbt,                      &
                dummy_dp,                   &
                nz,                         &
                M%ih_sym,                   &
                0)

    M = 0

    call dscal(length_pp, -1.0d0, prop_gradient_pp, 1)
    call dscal(length_pn, -1.0d0, prop_gradient_pn, 1)

    if (debug_me) then
      call header('prop_gradient_pp', -1)
      call print_q_vector(prop_gradient_pp, length_pp)
      call header('prop_gradient_pn', -1)
      call print_q_vector(prop_gradient_pn, length_pn)
    end if


!   parameters
!   ==========

    static = .false.
    if (nr_freq == 1 .and. freq(1) == 0.0d0) then
      static = .true.
    end if

    lineq     = .true.
    fkrmc     = .false.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .false.
    diaghe    = .true.
    iprxrs    = 0
    thcxrs    = 1.0d-9
    resfac    = 1.0d3
    maxitr    = 30
    nredm     = 400
    n2redm    = 400*400
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

    call alloc(ibtyp,            2*nredm)
    call alloc(ibtyp_pointer_c,  nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)
    call alloc(convergence,      nr_freq)
    call alloc(eigval,           nr_freq)
    call alloc(eigvec,           nredm*nr_freq)

    call dcopy(nr_freq, freq, 1, eigval, 1)

    call memget('INTE', kzxope, 2*length_pp, work, kfree, lfree)
    call memget('INTE', kzxopp, 2*length_pn, work, kfree, lfree)

    call fill_work(length_pp,                  &
                   length_pn,                  &
                   orbital_rotation_vector_pp, &
                   orbital_rotation_vector_pn, &
                   work(kzxope),               &
                   work(kzxopp))
!   if you release kzxope and kzxopp now you enter a world of pain
!   it is unfortunately used in several places later

    deallocate(orbital_rotation_vector_pp)
    deallocate(orbital_rotation_vector_pn)


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call xrsctl(dummy_dp,         &
                prop_gradient_pp, &
                prop_gradient_pn, &
                ibtyp,            &
                ibtyp_pointer_c,  &
                ibtyp_pointer_pp, &
                ibtyp_pointer_pn, &
                convergence,      &
                eigval,           &
                eigvec,           &
                work, kfree, lfree)

    call dealloc(convergence)
    call dealloc(eigval)

    call dealloc(ibtyp_pointer_c)


!   construct response vector
!   =========================

    call alloc(response_vector_pph, length_pp, nz, nr_freq)
    response_vector_pph = 0.0

    call construct_response_vector('pp',                &
                                   1,                   &
                                   response_vector_pph, &
                                   length_pp,           &
                                   nr_freq,             &
                                   ibtyp,               &
                                   ibtyp_pointer_pp,    &
                                   eigvec)

    if (.not. static) then

      call alloc(response_vector_ppa, length_pp, nz, nr_freq)
      response_vector_ppa = 0.0

      call construct_response_vector('pp',                &
                                     -1,                  &
                                     response_vector_ppa, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

    end if

    call alloc(response_vector_pnh, length_pn, nz, nr_freq)
    response_vector_pnh = 0.0

    call construct_response_vector('pn',                &
                                   1,                   &
                                   response_vector_pnh, &
                                   length_pn,           &
                                   nr_freq,             &
                                   ibtyp,               &
                                   ibtyp_pointer_pn,    &
                                   eigvec)

    if (.not. static) then

      call alloc(response_vector_pna, length_pn, nz, nr_freq)
      response_vector_pna = 0.0

      call construct_response_vector('pn',                &
                                     -1,                  &
                                     response_vector_pna, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

    end if

    if (debug_me) then
      call header('response_vector_pph', -1)
      call prbvec(6, response_vector_pph, nr_freq, length_pp)
      call header('response_vector_ppa', -1)
      call prbvec(6, response_vector_ppa, nr_freq, length_pp)
      call header('response_vector_pnh', -1)
      call prbvec(6, response_vector_pnh, nr_freq, length_pn)
      call header('response_vector_pna', -1)
      call prbvec(6, response_vector_pna, nr_freq, length_pn)
    end if

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)
    call dealloc(eigvec)


!   construct perturbed ao density matrices
!   =======================================

    ncsim = 0
    nesim = nr_freq
    npsim = nr_freq

    ntsim = ncsim + nesim + npsim

    call alloc(isymop, ntsim)
    call alloc(ihrmop, ntsim)
    call alloc(ifckop, ntsim)

    do i = 1, ntsim
      isymop(i) = jsymop
      ifckop(i) = 1
    end do

    call alloc(temp_ao_dmat, ntbas(0), ntbas(0), nz, ntsim)

    if (.not. isdef(D1)) call mat_init(D1, ntbas(0), ntbas(0))

    do i = 1, ntsim
      ihrmop(i) = 1
    end do

    temp_ao_dmat = 0.0d0
    call dmolitx(temp_ao_dmat,        &
                 response_vector_pph, &
                 response_vector_pnh, &
                 C%elements,          &
                 dummy_dp,            &
                 dummy_dp,            &
                 ihrmop,              &
                 nr_dmat,             &
                 work, kfree, lfree)

    call dealloc(response_vector_pph)
    call dealloc(response_vector_pnh)

    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1),           1, D1%elements, 1)
    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, D1%elements, 1)

    D1%pg_sym        = RHS%pg_sym
    D1%tr_sym        = RHS%tr_sym
    D1%fermi_pointer = RHS%fermi_pointer
    D1%ih_sym        = 1
!!!!!!!!! do freq

    if (.not. static) then

      do i = 1, ntsim
        ihrmop(i) = -1
      end do

      temp_ao_dmat = 0.0d0
      call dmolitx(temp_ao_dmat,        &
                   response_vector_ppa, &
                   response_vector_pna, &
                   C%elements,          &
                   dummy_dp,            &
                   dummy_dp,            &
                   ihrmop,              &
                   nr_dmat,             &
                   work, kfree, lfree)

      call dealloc(response_vector_ppa)
      call dealloc(response_vector_pna)

      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1),           1, D1%elements, 1)
      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, D1%elements, 1)

    end if

!   don't release this earlier
    call memrel('aoresponse', work, 1, kzxopp, kfree, lfree)
    call memrel('aoresponse', work, 1, kzxope, kfree, lfree)

    call dealloc(isymop)
    call dealloc(ihrmop)
    call dealloc(ifckop)

    call dealloc(temp_ao_dmat)

    call dealloc(work)

  end subroutine response_solver

  subroutine construct_response_vector(vector_type,     &
                                       ih,              &
                                       response_vector, &
                                       length,          &
                                       nr_freq,         &
                                       ibtyp,           &
                                       ibtyp_pointer,   &
                                       eigvec)
    character(*),     intent(in)    :: vector_type
    integer(kind_i4), intent(in)    :: ih
    real(kind_dp),    intent(inout) :: response_vector(:, :, :)
    integer(kind_i4), intent(in)    :: length
    integer(kind_i4), intent(in)    :: nr_freq
    integer(kind_i4), intent(in)    :: ibtyp(2, *)
    integer(kind_i4), intent(in)    :: ibtyp_pointer(*)
    real(kind_dp),    intent(in)    :: eigvec(*)

#include "dcbxrs.h"
#include "dgroup.h"

    integer(kind_i4), allocatable   :: ivecs(:)
    real(kind_dp),    allocatable   :: buffer(:)

    character(6)                    :: file_name
    integer(kind_i4)                :: file_unit, nbtyp

    call alloc(ivecs, nr_freq)

    select case (vector_type)
      case ('pp')
        nbtyp     = 1
        file_unit = luboe
        file_name = 'PAMBOE'
        call alloc(buffer, nzxopeq)
      case ('pn')
        nbtyp     = 2
        file_unit = lubop
        file_name = 'PAMBOP'
        call alloc(buffer, nzxoppq)
      case default
        call quit('error in construct_response_vector, unknown vector type')
    end select

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

    close(file_unit, &
          status = 'keep')

    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine construct_response_vector

  subroutine print_vector(vector, vector_dimension)

    real(kind_dp),    intent(in) :: vector(:)
    integer(kind_i4), intent(in) :: vector_dimension

#include "dgroup.h"

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                1,                &
                ipqtoq(1, 0),     &
                6)

  end subroutine print_vector

  subroutine print_q_vector(vector, vector_dimension)

    real(kind_dp),    intent(in) :: vector(:)
    integer(kind_i4), intent(in) :: vector_dimension

#include "dgroup.h"

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                nz,               &
                ipqtoq(1, 0),     &
                6)

  end subroutine print_q_vector

end module dirac_interface_m
