module aoresponse_main_m

  use kinds_m
! use scf_density_class
  use memory_allocator
  use dirac_matrix_operations_m
  use dirac_interface_m
  use propint
  use pertdens
  use aoresponse_input_m

  implicit none

#include "priunit.h"

! public test_propint
! public test_matrix_operations
! public aoresponse_driver
! public get_ao_pmat
! public solve_perturbed_ao_dmat

! private

contains

  subroutine test_propint()

#include "dcbbas.h"
#include "dgroup.h"

    type(matrix):: C, D, F, S,ff(3)
    real*8:: ef(3)
    integer:: noi(0),i,j,k
    character*4:: nop(0)

!    write(*, *) 'blub', find_op_index('xdiplen')
!    write(*, *) 'blub', find_op_index('xdiplen')
!    write(*, *) 'blub', find_op_index('ydiplen')
!    write(*, *) 'blub', find_op_index('zdiplen')

    call mat_init(S, ntbas(0), ntbas(0))
    call mat_init(D, ntbas(0), ntbas(0))
    call mat_init(C, ntbas(0), ntbas(0))
    call mat_init(F, ntbas(0), ntbas(0))

    call get_mo_coef(C)
    call get_zo_ao_dmat(C, D)


    ff(1)=0*s; ff(2)=0*s; ff(3)=0*s
    call propfock((/'EL'/),(/D/),ff)
    ff(1)=0; ff(2)=0; ff(3)=0
    ef=0; call propenergy((/'EL'/),(/D/),ef)
    print '("dipole moment -ef"/3f18.13/)',-ef
    print*

    call mat_free(s)
    call mat_free(C)
    call mat_free(d)
    call mat_free(f)

  end subroutine test_propint

  subroutine test_matrix_operations()

    type(matrix)     :: d22, d23, d32, d33
    type(matrix)     :: a, b, c, d, e
    type(matrix)     :: m, n
    real(kind_dp)    :: accept = 1.0d-5
    real(kind_dp)    :: r, s

    call mat_init(d22, 2, 2)
    call mat_init(d23, 2, 3)
    call mat_init(d32, 3, 2)
    call mat_init(d33, 3, 3)

    d22%elements = 1.0d0
    d23%elements = 1.0d0
    d32%elements = 1.0d0
    d33%elements = 1.0d0

    a = 1.0d0*d22
    b = 0.5d0*d22
    c = 0.5d0*d23
    d = 0.5001d0*d23
    e = 0.500000001d0*d23

    call compare_matrices(a, b, accept, 'testing compare_matrices, should fail due to values')
    call compare_matrices(a, c, accept, 'testing compare_matrices, should fail due to shapes')
    call compare_matrices(b, c, accept, 'testing compare_matrices, should fail due to shapes')
    call compare_matrices(c, d, accept, 'testing compare_matrices, should fail due to values')
    call compare_matrices(c, e, accept, 'testing compare_matrices')

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    m = 1.0d0*d23
    n = 1.0d0*d32

    m%elements = 0.0d0
    n%elements = 0.0d0

    m%elements( 1) =  1.0d0
    m%elements( 2) =  2.0d0
    m%elements( 3) =  3.0d0
    m%elements( 4) =  4.0d0
    m%elements( 5) =  5.0d0
    m%elements( 6) =  6.0d0
    m%elements( 7) =  7.0d0
    m%elements( 8) =  8.0d0
    m%elements( 9) =  9.0d0
    m%elements(10) = 10.0d0
    m%elements(11) = 11.0d0
    m%elements(12) = 12.0d0
    m%elements(13) = 13.0d0
    m%elements(14) = 14.0d0
    m%elements(15) = 15.0d0
    m%elements(16) = 16.0d0
    m%elements(17) = 17.0d0
    m%elements(18) = 18.0d0
    m%elements(19) = 19.0d0
    m%elements(20) = 20.0d0
    m%elements(21) = 21.0d0
    m%elements(22) = 22.0d0
    m%elements(23) = 23.0d0
    m%elements(24) = 24.0d0

    n%elements( 1) =  1.0d0
    n%elements( 2) =  2.0d0
    n%elements( 3) =  3.0d0
    n%elements( 4) =  4.0d0
    n%elements( 5) =  5.0d0
    n%elements( 6) =  6.0d0
    n%elements( 7) =  7.0d0
    n%elements( 8) =  8.0d0
    n%elements( 9) =  9.0d0
    n%elements(10) = 10.0d0
    n%elements(11) = 11.0d0
    n%elements(12) = 12.0d0
    n%elements(13) = 13.0d0
    n%elements(14) = 14.0d0
    n%elements(15) = 15.0d0
    n%elements(16) = 16.0d0
    n%elements(17) = 17.0d0
    n%elements(18) = 18.0d0
    n%elements(19) = 19.0d0
    n%elements(20) = 20.0d0
    n%elements(21) = 21.0d0
    n%elements(22) = 22.0d0
    n%elements(23) = 23.0d0
    n%elements(24) = 24.0d0

    a = 1.0d0*m
    b = 1.0d0*n

    c = 2.0d0*a
    d = a + a + a + (-a)
    e = 4.0d0*a - a - a

    call compare_matrices(c, d, accept, 'testing matrix plus/minus')
    call compare_matrices(d, e, accept, 'testing matrix plus/minus')
    call compare_matrices(a, e, accept, 'testing matrix plus/minus, should fail due to values')

    r = a**a

    call compare_dp(r, 4900.0d0, accept, 'checking matrix ddot')

    r = b**b
    s = norm(a)*norm(b)

    call compare_dp(r, s, accept, 'checking matrix_norm')

    c = a/0.5d0
    d = a/2 + a/2 + a

    call compare_matrices(c, d, accept, 'testing matrix/r and matrix/i')












    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    m = 0
    n = 0

    call mat_free(d22)
    call mat_free(d23)
    call mat_free(d32)
    call mat_free(d33)

  end subroutine test_matrix_operations

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcbxrs.h"
#include "dgroup.h"

    integer(kind_i4), allocatable   :: op_index(:)
    real(kind_dp),    allocatable   :: freq(:)
    real(kind_dp),    allocatable   :: zo_ao_dmat(:, :, :)
    real(kind_dp),    allocatable   :: perturbed_ao_dmat(:, :, :, :, :, :)
    real(kind_dp),    allocatable   :: ao_pmat(:, :, :, :)
    real(kind_dp),    allocatable   :: mo_coef(:, :, :)

    integer(kind_i4)                :: nr_op, nr_freq, iop, ifreq
!!!!!!!!!!!!!!!!!!
    integer(kind_i4)                :: i, j
!!!!!!!!!!!!!!!!!!
    real(kind_dp)                   :: ddot, value

    nr_op   = 1
    nr_freq = 1

    call alloc(op_index,          nr_op)
    call alloc(freq,              nr_freq)
    call alloc(zo_ao_dmat,        ntbas(0), ntbas(0), nz)
    call alloc(perturbed_ao_dmat, ntbas(0), ntbas(0), nz, nr_freq, 2, nr_op)
    call alloc(ao_pmat,           ntbas(0), ntbas(0), nz, nr_op)
    call alloc(mo_coef,           ntbas(0), ntbas(0), nz)

    op_index          = 0
    freq              = 0.0d0
    zo_ao_dmat        = 0.0d0
    perturbed_ao_dmat = 0.0d0
    ao_pmat           = 0.0d0
    mo_coef           = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!   write(*, *) 'index ex: ', e_field%component(1)%op_index
!   write(*, *) 'index ey: ', e_field%component(2)%op_index
!   write(*, *) 'index ez: ', e_field%component(3)%op_index
    op_index(1) = 5
    freq(1) = 0.0d0
!   freq(2) = 0.2d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    static = .false.
    if (nr_freq == 1 .and. freq(1) == 0.0d0) then
      static = .true.
    end if

!   call get_mo_coef(mo_coef)
!   call get_zo_ao_dmat(zo_ao_dmat, mo_coef)

    do iop = 1, nr_op
!     call get_ao_pmat(op_index(iop), ao_pmat(:, :, :, iop))
      call solve_perturbed_ao_dmat(op_index(iop), freq, perturbed_ao_dmat(:, :, :, :, :, iop))
    end do

    do iop = 1, nr_op

      value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, zo_ao_dmat, 1)
      write(*, *) 'd(0)*p = ', value

      do ifreq = 1, nr_freq

        value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, perturbed_ao_dmat(:, :, :, ifreq, 1, iop), 1)
        write(*, *) 'freq = ', freq(ifreq), '    d(1)*p = ', value

      end do
    end do

    call dealloc(op_index)
    call dealloc(freq)
    call dealloc(zo_ao_dmat)
    call dealloc(perturbed_ao_dmat)
    call dealloc(ao_pmat)
    call dealloc(mo_coef)

  end subroutine aoresponse_driver



  subroutine get_ao_pmat_old(op_index, ao_pmat)
    integer(kind_i4), intent(in)  :: op_index
    real(kind_dp),    intent(out) :: ao_pmat(:, :, :)

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    real(kind_dp),    allocatable :: ptri(:)
    real(kind_dp),    allocatable :: op1int(:)
    logical,          allocatable :: first(:)

    real(kind_dp),    allocatable :: work(:)
    integer(kind_i4)              :: kfree, lfree

    logical                       :: file_exists, is_square
    real(kind_dp)                 :: dummy_dp
    integer(kind_i4)              :: file_unit, dummy_i4, idimension

    is_square = .false.
    if (iprptim(op_index) == 0) then
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

    file_exists = .false.
    inquire(file = 'AOPROPER', exist = file_exists)
    if (.not. file_exists) then
      call quit('file aoproper not found')
    end if

    open(file_unit,           &
         file   = 'AOPROPER', &
         status = 'unknown',  &
         form   = 'unformatted')

    rewind file_unit

    call prpex2(op_index,   &
                dummy_dp,   &
                .false.,    &
                idimension, &
                ptri,       &
                ao_pmat,    &
                op1int,     &
                first,      &
                work, kfree, lfree, 0)

    close(file_unit, &
          status = 'keep')

    call dealloc(ptri)
    call dealloc(op1int)
    deallocate(first)
    call dealloc(work)

  end subroutine get_ao_pmat_old








  subroutine solve_perturbed_ao_dmat(op_index, freq, perturbed_ao_dmat)
    integer(kind_i4), intent(in)  :: op_index
    real(kind_dp),    intent(in)  :: freq(:)
    real(kind_dp),    intent(out) :: perturbed_ao_dmat(:, :, :, :, :)

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    integer(kind_i4), allocatable :: ibtyp(:)
    integer(kind_i4), allocatable :: ibtyp_pointer_c(:)
    integer(kind_i4), allocatable :: ibtyp_pointer_pp(:)
    integer(kind_i4), allocatable :: ibtyp_pointer_pn(:)

    integer(kind_i4), allocatable :: isymop(:)
    integer(kind_i4), allocatable :: ihrmop(:)
    integer(kind_i4), allocatable :: ifckop(:)

    real(kind_dp),    allocatable :: eigval(:)
    real(kind_dp),    allocatable :: eigvec(:)
    real(kind_dp),    allocatable :: convergence(:)

    real(kind_dp),    allocatable :: prop_gradient_c(:)
    real(kind_dp),    allocatable :: prop_gradient_pp(:)
    real(kind_dp),    allocatable :: prop_gradient_pn(:)

    real(kind_dp),    allocatable :: response_vector_pph(:, :, :)
    real(kind_dp),    allocatable :: response_vector_ppa(:, :, :)
    real(kind_dp),    allocatable :: response_vector_pnh(:, :, :)
    real(kind_dp),    allocatable :: response_vector_pna(:, :, :)

    real(kind_dp),    allocatable :: temp_ao_dmat(:, :, :, :)
    real(kind_dp),    allocatable :: mo_coef(:, :, :)

    real(kind_dp),    allocatable :: work(:)
    integer(kind_i4)              :: kfree, lfree

    real(kind_dp)                 :: dummy_dp
    logical                       :: skip_pp, skip_pn, do_pp, do_pn
    logical                       :: debug_me = .false.
    character(72)                 :: orbital_string(3,2)

    integer(kind_i4)              :: i, ifsym, nr_freq, nr_dmat
    integer(kind_i4)              :: length_pp, length_pn
    integer(kind_i4)              :: ill_int, isl_int, iss_int, igt_int

    skip_pp = .false.
    skip_pn = .false.


!   set common block
!   ================

    call setdcbidx()


!   get orbital string
!   ==================

    do ifsym = 1, nfsym
      orbital_string(1, ifsym) = ' '
      orbital_string(2, ifsym) = ' '
      orbital_string(3, ifsym) = ' '
    end do

    kfree = 1
    lfree = 2*norbt
    call alloc(work, lfree)

    call ostring(orbital_string, &
                 skip_pp,        &
                 skip_pn,        &
                 work, kfree, lfree)

    call dealloc(work)


!   initialize symmetry independent variables
!   =========================================

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


!   initialize configurational and orbital parameters
!   =================================================

    nr_freq = size(freq)
    jsymop  = iprpsym(op_index)
    jtimop  = iprptim(op_index)
    jopsy   = jbtof(jsymop - 1, 1)
    nfreq   = nr_freq
    nexsim  = nr_freq
    nexstv  = nr_freq
    nexcnv  = nr_freq
    ncred   = 0
    nered   = 0
    npred   = 0
    nzred   = ncred + nered + npred

    kfree = 1
    lfree = 200*ntbas(0)*ntbas(0) ! no idea how much is needed
    call alloc(work, lfree)

    call xrspar(orbital_string, &
                skip_pp,        &
                skip_pn,        &
                work, kfree, lfree, 0)



!   get property gradient
!   =====================

    nzconfq = 0

    call alloc(prop_gradient_c,  nzconfq)
    call alloc(prop_gradient_pp, nzxopeq)
    call alloc(prop_gradient_pn, nzxoppq)

    prop_gradient_c  = 0.0d0
    prop_gradient_pp = 0.0d0
    prop_gradient_pn = 0.0d0

    call pamgrd(op_index,         &
                prop_gradient_c,  &
                prop_gradient_pp, &
                prop_gradient_pn, &
                jopsy,            &
                work(kzconf),     &
                work(kzxope),     &
                work(kzxopp),     &
                nzconf,           &
                nzxope,           &
                nzxopp,           &
                work, kfree, lfree, 0)

!   call dealloc(work)

    if (debug_me) then
      call header('prop_gradient_pp',-1)
!     call print_q_vector(prop_gradient_pp, nzxopeq)
      call header('prop_gradient_pn',-1)
!     call print_q_vector(prop_gradient_pn, nzxoppq)
    end if

!   todo: check the norm of the property gradient
!         and if below threshold then skip things


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call alloc(ibtyp,            2*nredm)
    call alloc(ibtyp_pointer_c,  nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)
    call alloc(convergence,      nr_freq)
    call alloc(eigval,           nr_freq)
    call alloc(eigvec,           nredm*nr_freq)

!   kfree = 1
!   lfree = 200*ntbas(0)*ntbas(0) ! no idea how much is needed
!   call alloc(work, lfree)

    call dcopy(nr_freq, freq, 1, eigval, 1)

    call xrsctl(prop_gradient_c,  &
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

    call dealloc(prop_gradient_c)
    call dealloc(prop_gradient_pp)
    call dealloc(prop_gradient_pn)
    call dealloc(convergence)
    call dealloc(eigval)
!   call dealloc(work)

!   construct the pp and pn response vectors
!   ========================================

    do_pp = .false.
    do_pn = .false.

    length_pp = nzxope
    length_pn = nzxopp

    if (length_pp > 0 .and. .not. skip_pp) then
      do_pp = .true.
    end if
    if (length_pn > 0 .and. .not. skip_pn) then
      do_pn = .true.
    end if

    if (do_pp) then

      call alloc(response_vector_pph, length_pp, nz, nr_freq)
      call alloc(response_vector_ppa, length_pp, nz, nr_freq)

      response_vector_pph = 0.0
      response_vector_ppa = 0.0

      call construct_response_vector('pp',                &
                                     1,                   &
                                     response_vector_pph, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      if (.not. static) then
        call construct_response_vector('pp',                &
                                       -1,                  &
                                       response_vector_ppa, &
                                       length_pp,           &
                                       nr_freq,             &
                                       ibtyp,               &
                                       ibtyp_pointer_pp,    &
                                       eigvec)
      end if

      if (debug_me) then
        call header('response_vector_pph',-1)
        call prbvec(lupri, response_vector_pph, nr_freq, length_pp)
        call header('response_vector_ppa',-1)
        call prbvec(lupri, response_vector_ppa, nr_freq, length_pp)
      end if

    end if

    if (do_pn) then

      call alloc(response_vector_pnh, length_pn, nz, nr_freq)
      call alloc(response_vector_pna, length_pn, nz, nr_freq)

      response_vector_pnh = 0.0
      response_vector_pna = 0.0

      call construct_response_vector('pn',                &
                                     1,                   &
                                     response_vector_pnh, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      if (.not. static) then
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
        call header('response_vector_pnh',-1)
        call prbvec(lupri, response_vector_pnh, nr_freq, length_pn)
        call header('response_vector_pna',-1)
        call prbvec(lupri, response_vector_pna, nr_freq, length_pn)
      end if

    end if

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_c)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)
    call dealloc(eigvec)


!   construct perturbed ao density matrices
!   =======================================

    nesim = 0
    npsim = 0

    if (do_pp) then
      nesim = nr_freq
    end if
    if (do_pn) then
      npsim = nr_freq
    end if

    ntsim = nesim + npsim

    call alloc(isymop, ntsim)
    call alloc(ihrmop, ntsim)
    call alloc(ifckop, ntsim)

    do i = 1, ntsim
      isymop(i) = jsymop
      ifckop(i) = 1
    end do

    call alloc(temp_ao_dmat, ntbas(0), ntbas(0), nz, ntsim)
    call alloc(mo_coef,      ntbas(0), ntbas(0), nz)

    perturbed_ao_dmat = 0.0d0
    mo_coef           = 0.0d0

!   call get_mo_coef(mo_coef)

    do i = 1, ntsim
      ihrmop(i) = 1
    end do

    temp_ao_dmat = 0.0d0

!   kfree = 1
!   lfree = 200*ntbas(0)*ntbas(0) ! no idea how much is needed
!   call alloc(work, lfree)

    call dmolitx(temp_ao_dmat,        &
                 response_vector_pph, &
                 response_vector_pnh, &
                 mo_coef,             &
                 dummy_dp,            &
                 dummy_dp,            &
                 ihrmop,              &
                 nr_dmat,             &
                 work, kfree, lfree)

    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 1.0d0, temp_ao_dmat(:, :, :, 1),           1, perturbed_ao_dmat(:, :, :, 1, 1), 1)
    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 1.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, perturbed_ao_dmat(:, :, :, 1, 1), 1)

    if (.not. static) then

      do i = 1, ntsim
        ihrmop(i) = -1
      end do

      temp_ao_dmat = 0.0d0

      call dmolitx(temp_ao_dmat,        &
                   response_vector_ppa, &
                   response_vector_pna, &
                   mo_coef,             &
                   dummy_dp,            &
                   dummy_dp,            &
                   ihrmop,              &
                   nr_dmat,             &
                   work, kfree, lfree)

      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 1.0d0, temp_ao_dmat(:, :, :, 1),           1, perturbed_ao_dmat(:, :, :, 1, 2), 1)
      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 1.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, perturbed_ao_dmat(:, :, :, 1, 2), 1)

    end if

    call dscal(ntbas(0)*ntbas(0)*nz*nr_freq*2, 2.0d0, perturbed_ao_dmat, 1)

    call dealloc(temp_ao_dmat)
    call dealloc(mo_coef)

    call dealloc(isymop)
    call dealloc(ihrmop)
    call dealloc(ifckop)

    if (do_pp) then
      call dealloc(response_vector_pph)
      call dealloc(response_vector_ppa)
    end if
    if (do_pn) then
      call dealloc(response_vector_pnh)
      call dealloc(response_vector_pna)
    end if

    call dealloc(work)

  end subroutine solve_perturbed_ao_dmat

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

end module aoresponse_main_m
