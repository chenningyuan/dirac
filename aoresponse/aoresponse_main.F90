module aoresponse_main_m

  use kinds_m
  use scf_density_class
  use memory_allocator
  use dirac_matrix_operations_m
  use propint
  use pertdens
  use aoresponse_input_m

  implicit none

#include "priunit.h"

  public test_propint
  public test_matrix_operations
  public aoresponse_driver
  public get_ao_pmat
  public solve_perturbed_ao_dmat

  private

contains

  subroutine test_propint()

    type(Matrix):: S,D,F,Ff(3)
    real*8:: Ef(3)
    integer:: NOI(0),i,j,k
    character*4:: NOP(0)
    S=S; D=D; F=F
    Ff(1)=0*S; Ff(2)=0*S; Ff(3)=0*S
    call PropFock((/'EL'/),(/D/),Ff)
    Ff(1)=0; Ff(2)=0; Ff(3)=0
    Ef=0; call PropEnergy((/'EL'/),(/D/),Ef)
    print '("Dipole moment -Ef"/3f18.13/)',-Ef
    print*


  end subroutine test_propint

  subroutine test_matrix_operations()

    type(matrix)     :: D22, D23, D32, D33
    type(matrix)     :: A, B, C, D, E
    type(matrix)     :: M, N
    real(kind_dp)    :: accept = 1.0d-5
    real(kind_dp)    :: r, s

    call mat_init(D22, 2, 2)
    call mat_init(D23, 2, 3)
    call mat_init(D32, 3, 2)
    call mat_init(D33, 3, 3)

    D22%elements = 1.0d0
    D23%elements = 1.0d0
    D32%elements = 1.0d0
    D33%elements = 1.0d0

    A = 1.0d0*D22
    B = 0.5d0*D22
    C = 0.5d0*D23
    D = 0.5001d0*D23
    E = 0.500000001d0*D23

    call compare_matrices(A, B, accept, 'testing compare_matrices, should fail due to values')
    call compare_matrices(A, C, accept, 'testing compare_matrices, should fail due to shapes')
    call compare_matrices(B, C, accept, 'testing compare_matrices, should fail due to shapes')
    call compare_matrices(C, D, accept, 'testing compare_matrices, should fail due to values')
    call compare_matrices(C, E, accept, 'testing compare_matrices')

    A = 0
    B = 0
    C = 0
    D = 0
    E = 0

    M = 1.0d0*D23
    N = 1.0d0*D32

    M%elements = 0.0d0
    N%elements = 0.0d0

    M%elements( 1) =  1.0d0
    M%elements( 2) =  2.0d0
    M%elements( 3) =  3.0d0
    M%elements( 4) =  4.0d0
    M%elements( 5) =  5.0d0
    M%elements( 6) =  6.0d0
    M%elements( 7) =  7.0d0
    M%elements( 8) =  8.0d0
    M%elements( 9) =  9.0d0
    M%elements(10) = 10.0d0
    M%elements(11) = 11.0d0
    M%elements(12) = 12.0d0
    M%elements(13) = 13.0d0
    M%elements(14) = 14.0d0
    M%elements(15) = 15.0d0
    M%elements(16) = 16.0d0
    M%elements(17) = 17.0d0
    M%elements(18) = 18.0d0
    M%elements(19) = 19.0d0
    M%elements(20) = 20.0d0
    M%elements(21) = 21.0d0
    M%elements(22) = 22.0d0
    M%elements(23) = 23.0d0
    M%elements(24) = 24.0d0

    N%elements( 1) =  1.0d0
    N%elements( 2) =  2.0d0
    N%elements( 3) =  3.0d0
    N%elements( 4) =  4.0d0
    N%elements( 5) =  5.0d0
    N%elements( 6) =  6.0d0
    N%elements( 7) =  7.0d0
    N%elements( 8) =  8.0d0
    N%elements( 9) =  9.0d0
    N%elements(10) = 10.0d0
    N%elements(11) = 11.0d0
    N%elements(12) = 12.0d0
    N%elements(13) = 13.0d0
    N%elements(14) = 14.0d0
    N%elements(15) = 15.0d0
    N%elements(16) = 16.0d0
    N%elements(17) = 17.0d0
    N%elements(18) = 18.0d0
    N%elements(19) = 19.0d0
    N%elements(20) = 20.0d0
    N%elements(21) = 21.0d0
    N%elements(22) = 22.0d0
    N%elements(23) = 23.0d0
    N%elements(24) = 24.0d0

    A = 1.0d0*M
    B = 1.0d0*N

    C = 2.0d0*A
    D = A + A + A + (-A)
    E = 4.0d0*A - A - A

    call compare_matrices(C, D, accept, 'testing matrix plus/minus')
    call compare_matrices(D, E, accept, 'testing matrix plus/minus')
    call compare_matrices(A, E, accept, 'testing matrix plus/minus, should fail due to values')

    r = A**A

    call compare_dp(r, 4900.0d0, accept, 'checking matrix ddot')

    r = B**B
    s = norm(A)*norm(B)

    call compare_dp(r, s, accept, 'checking matrix_norm')

    C = A/0.5d0
    D = A/2 + A/2 + A

    call compare_matrices(C, D, accept, 'testing matrix/r and matrix/i')












    A = 0
    B = 0
    C = 0
    D = 0
    E = 0

    M = 0
    N = 0

    call mat_free(D22)
    call mat_free(D23)
    call mat_free(D32)
    call mat_free(D33)

  end subroutine test_matrix_operations

  subroutine aoresponse_driver(work, kfree, lfree)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer(kind_i4), intent(inout) :: kfree
    integer(kind_i4), intent(inout) :: lfree

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
    do i = 1, 2
      do j = 1, 3

        write(*, *) 'blub', i, j, perturbation_field(i)%component(j)%op_index

      end do
    end do



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

    call get_mo_coef(mo_coef)
    call get_zo_ao_dmat(zo_ao_dmat, mo_coef)

    do iop = 1, nr_op

      call get_ao_pmat(op_index(iop), ao_pmat(:, :, :, iop))

      call solve_perturbed_ao_dmat(op_index(iop),                         &
                                   freq,                                  &
                                   perturbed_ao_dmat(:, :, :, :, :, iop), &
                                   work, kfree, lfree)
    end do


    do iop = 1, nr_op

      value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, zo_ao_dmat, 1)
      write(*, *) 'D(0)*P = ', value

      do ifreq = 1, nr_freq

        value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, perturbed_ao_dmat(:, :, :, ifreq, 1, iop), 1)
        write(*, *) 'freq = ', freq(ifreq), '    D(1)*P = ', value

      end do
    end do

    call dealloc(op_index)
    call dealloc(freq)
    call dealloc(zo_ao_dmat)
    call dealloc(perturbed_ao_dmat)
    call dealloc(ao_pmat)
    call dealloc(mo_coef)

  end subroutine aoresponse_driver

  subroutine get_ao_pmat(op_index, ao_pmat)
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

    logical                       :: file_exists, is_square
    real(kind_dp)                 :: dummy_dp
    integer(kind_i4)              :: file_unit, dummy_i4, idimension, kfree, lfree

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
    call alloc(work, 2*nz*n2bbasx)
    kfree = 1
    lfree = 2*nz*n2bbasx

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

  end subroutine get_ao_pmat

  subroutine solve_perturbed_ao_dmat(op_index,          &
                                     freq,              &
                                     perturbed_ao_dmat, &
                                     work, kfree, lfree)
    integer(kind_i4), intent(in)    :: op_index
    real(kind_dp),    intent(in)    :: freq(:)
    real(kind_dp),    intent(out)   :: perturbed_ao_dmat(:, :, :, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer(kind_i4), intent(inout) :: kfree
    integer(kind_i4), intent(inout) :: lfree

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dcbbas.h"
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

    real(kind_dp),    allocatable   :: prop_gradient_c(:)
    real(kind_dp),    allocatable   :: prop_gradient_pp(:)
    real(kind_dp),    allocatable   :: prop_gradient_pn(:)

    real(kind_dp),    allocatable   :: response_vector_pph(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_ppa(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_pnh(:, :, :)
    real(kind_dp),    allocatable   :: response_vector_pna(:, :, :)

    real(kind_dp),    allocatable   :: temp_ao_dmat(:, :, :, :)
    real(kind_dp),    allocatable   :: mo_coef(:, :, :)

    real(kind_dp)                   :: dummy_dp
    logical                         :: skip_pp, skip_pn, do_pp, do_pn
    logical                         :: debug_me = .false.
    character(72)                   :: orbital_string(3,2)

    integer(kind_i4)                :: i, ifsym, nr_freq, nr_dmat
    integer(kind_i4)                :: length_pp, length_pn
    integer(kind_i4)                :: ill_int, isl_int, iss_int, igt_int

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

    call ostring(orbital_string, &
                 skip_pp,        &
                 skip_pn,        &
                 work, kfree, lfree)


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

    call xrspar(orbital_string, &
                skip_pp,        &
                skip_pn,        &
                work,kfree,lfree, 0)


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

    if (debug_me) then
      call header('prop_gradient_pp',-1)
      call print_q_vector(prop_gradient_pp, nzxopeq)
      call header('prop_gradient_pn',-1)
      call print_q_vector(prop_gradient_pn, nzxoppq)
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
                work,kfree,lfree)

    call dealloc(prop_gradient_c)
    call dealloc(prop_gradient_pp)
    call dealloc(prop_gradient_pn)
    call dealloc(convergence)
    call dealloc(eigval)


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

    call get_mo_coef(mo_coef)

    do i = 1, ntsim
      ihrmop(i) = 1
    end do

    temp_ao_dmat = 0.0d0

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
