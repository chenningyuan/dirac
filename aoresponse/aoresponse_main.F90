module aoresponse_main_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use dirac_interface_m
  use propint
  use pertdens
  use aoresponse_input_m

  implicit none

! public something

! private

contains

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcbxrs.h"
#include "dgroup.h"

    type(matrix)               :: A, C, D, D1H, D1A

    real(kind_dp), allocatable :: freq(:)
    integer(kind_i4)           :: nr_freq = 1

    real(kind_dp), allocatable :: prop_gradient_pp(:)
    real(kind_dp), allocatable :: prop_gradient_pn(:)

    real(kind_dp), allocatable :: work(:)
    integer(kind_i4)           :: kfree, lfree

    call mat_init(C, ntbas(0), ntbas(0))
    call get_mo_coef(C)

    call alloc(freq, nr_freq)
    freq = 0.0d0

!!!!!!!!!!
    freq(1) = 0.02d0
!   freq(1) = 0.00d0
!!!!!!!!!!

    A   = 1.0d-100*C
    D   = 1.0d-100*C
    D1H = 1.0d-100*C
    D1A = 1.0d-100*C

    call get_zo_ao_dmat(C, D)
    call get_ao_pmat('ZDIPLEN', A)

!!!!!!!!!!
    kfree = 1
    lfree = 400*ntbas(0)*ntbas(0) ! no idea how much is needed
    call alloc(work, lfree)
!!!!!!!!!!

!!!!!!!!!!
    nzxopeq = 560
    nzxoppq = 660
    call alloc(prop_gradient_pp, nzxopeq)
    call alloc(prop_gradient_pn, nzxoppq)
    prop_gradient_pp = 0.0d0
    prop_gradient_pn = 0.0d0
!!!!!!!!!!

    call get_prop_gradient('ZDIPLEN', prop_gradient_pp, prop_gradient_pn, work, kfree, lfree)
    call solve_perturbed_ao_dmat('ZDIPLEN', prop_gradient_pp, prop_gradient_pn, freq, C, D1H, D1A, work, kfree, lfree)
!   call quit('over and out')

    call dealloc(work)
    call dealloc(prop_gradient_pp)
    call dealloc(prop_gradient_pn)

    write(*, *) 'A**D   = ', A**D
    write(*, *) 'A**D1H = ', A**D1H
    write(*, *) 'A**D1A = ', A**D1A

    write(*, *) 'compare: -4.828661928957 at 0.00'
    write(*, *) 'compare: -4.833040214626 at 0.02'

    A   = 0
    D   = 0
    D1H = 0
    D1A = 0

    call mat_free(C)
    call dealloc(freq)

  end subroutine aoresponse_driver

  subroutine get_prop_gradient(op_name, prop_gradient_pp, prop_gradient_pn, work, kfree, lfree)

    character(*),     intent(in)    :: op_name
    real(kind_dp),    intent(inout) :: prop_gradient_pp(:)
    real(kind_dp),    intent(inout) :: prop_gradient_pn(:)
    real(kind_dp),    intent(inout) :: work(:)
    integer(kind_i4), intent(inout) :: kfree, lfree

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"

    integer(kind_i4)                :: op_index, ifsym
    real(kind_dp)                   :: dummy_dp
    logical                         :: skip_pp, skip_pn
    logical                         :: debug_me = .false.
    character(72)                   :: orbital_string(3,2)

    op_index = find_op_index(op_name)
    jsymop   = iprpsym(op_index)
    jtimop   = iprptim(op_index)
    jopsy    = jbtof(jsymop - 1, 1)

    call setdcbidx()

    do ifsym = 1, nfsym
      orbital_string(1, ifsym) = ' '
      orbital_string(2, ifsym) = ' '
      orbital_string(3, ifsym) = ' '
    end do

    skip_pp  = .false.
    skip_pn  = .false.

    call ostring(orbital_string, &
                 skip_pp,        &
                 skip_pn,        &
                 work, kfree, lfree)

    call xrspar(orbital_string, &
                skip_pp,        &
                skip_pn,        &
                work, kfree, lfree, 0)

    nzconf = 0

    call pamgrd(op_index,         &
                dummy_dp,         &
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

    if (debug_me) then
      call header('prop_gradient_pp', -1)
      call print_q_vector(prop_gradient_pp, nzxope)
      call header('prop_gradient_pn', -1)
      call print_q_vector(prop_gradient_pn, nzxopp)
    end if

  end subroutine get_prop_gradient

  subroutine solve_perturbed_ao_dmat(op_name, prop_gradient_pp, prop_gradient_pn, freq, C, D1H, D1A, work, kfree, lfree)

    character(*),     intent(in)    :: op_name
    real(kind_dp),    intent(in)    :: prop_gradient_pp(:)
    real(kind_dp),    intent(in)    :: prop_gradient_pn(:)
    real(kind_dp),    intent(in)    :: freq(:)
    type(matrix),     intent(inout) :: C, D1H, D1A
    real(kind_dp),    intent(inout) :: work(:)
    integer(kind_i4), intent(inout) :: kfree, lfree

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    integer(kind_i4)                :: op_index

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

    real(kind_dp),    allocatable :: response_vector_pph(:, :, :)
    real(kind_dp),    allocatable :: response_vector_ppa(:, :, :)
    real(kind_dp),    allocatable :: response_vector_pnh(:, :, :)
    real(kind_dp),    allocatable :: response_vector_pna(:, :, :)

    real(kind_dp),    allocatable :: temp_ao_dmat(:, :, :, :)
!   real(kind_dp),    allocatable :: mo_coef(:, :, :)


    real(kind_dp)                 :: dummy_dp
    logical                       :: skip_pp, skip_pn, do_pp, do_pn
    logical                       :: debug_me = .true.
    character(72)                 :: orbital_string(3,2)

    integer(kind_i4)              :: i, ifsym, nr_freq, nr_dmat
    integer(kind_i4)              :: length_pp, length_pn
    integer(kind_i4)              :: ill_int, isl_int, iss_int, igt_int

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

    op_index = find_op_index(op_name)
    nr_freq   = size(freq)
    jsymop    = iprpsym(op_index)
    jtimop    = iprptim(op_index)
    jopsy     = jbtof(jsymop - 1, 1)
    nfreq     = nr_freq
    nexsim    = nr_freq
    nexstv    = nr_freq
    nexcnv    = nr_freq
    ncred     = 0
    nered     = 0
    npred     = 0
    nzred     = ncred + nered + npred

!   call xrspar(orbital_string, &
!               skip_pp,        &
!               skip_pn,        &
!               work, kfree, lfree, 0)


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call alloc(ibtyp,            2*nredm)
    call alloc(ibtyp_pointer_c,  nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)
    call alloc(convergence,      nr_freq)
    call alloc(eigval,           nr_freq)
    call alloc(eigvec,           nredm*nr_freq)



    call dcopy(nr_freq, freq, 1, eigval, 1)

    nzconfq = 0

    call alloc(prop_gradient_c,  nzconfq)

    prop_gradient_c  = 0.0d0


    call xrsctl(prop_gradient_c,  &
!   call xrsctl(0.0d0,  &
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
!   call dealloc(prop_gradient_pp)
!   call dealloc(prop_gradient_pn)
    call dealloc(convergence)
    call dealloc(eigval)

    call dealloc(ibtyp_pointer_c)





    length_pp = size(prop_gradient_pp)/nz
    length_pn = size(prop_gradient_pn)/nz

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

    nesim = 0
    npsim = 0

!   if (do_pp) then
      nesim = nr_freq
!   end if
!   if (do_pn) then
      npsim = nr_freq
!   end if

    ntsim = nesim + npsim

    call alloc(isymop, ntsim)
    call alloc(ihrmop, ntsim)
    call alloc(ifckop, ntsim)

    do i = 1, ntsim
      isymop(i) = jsymop
      ifckop(i) = 1
    end do

!   call alloc(mo_coef,      ntbas(0), ntbas(0), nz)

!   perturbed_ao_dmat = 0.0d0
!   mo_coef           = 0.0d0

!   call get_mo_coef(mo_coef)





    call alloc(temp_ao_dmat, ntbas(0), ntbas(0), nz, ntsim)

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

    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1),           1, D1H%elements, 1)
    call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, D1H%elements, 1)

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

      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1),           1, D1A%elements, 1)
      call daxpy(ntbas(0)*ntbas(0)*nz*nr_freq, 2.0d0, temp_ao_dmat(:, :, :, 1 + nr_freq), 1, D1A%elements, 1)

    end if

    call dealloc(isymop)
    call dealloc(ihrmop)
    call dealloc(ifckop)

    call dealloc(temp_ao_dmat)

    call dealloc(response_vector_ppa)
    call dealloc(response_vector_pna)

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

end module aoresponse_main_m
