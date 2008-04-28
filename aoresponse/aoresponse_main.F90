module aoresponse_main_m

  use kinds_m
  use scf_density_class
  use memory_allocator

  implicit none

#include "priunit.h"

  public aoresponse_driver
  public get_ao_pmat
  public solve_perturbed_ao_dmat

  private

contains

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
!   real(kind_dp),    allocatable   :: perturbed_ao_dmat(:, :, :, :, :)
    real(kind_dp),    allocatable   :: perturbed_ao_dmat(:, :, :, :)
    real(kind_dp),    allocatable   :: ao_pmat(:, :, :, :)
    real(kind_dp),    allocatable   :: mo_coef(:, :, :)

    integer(kind_i4)                :: nr_op, nr_freq, iop
    real(kind_dp)                   :: ddot, value

    nr_op   = 1
    nr_freq = 1

    allocate(op_index(nr_op))
    allocate(freq(nr_freq))
    allocate(zo_ao_dmat(ntbas(0), ntbas(0), nz))
!   allocate(perturbed_ao_dmat(ntbas(0), ntbas(0), nz, nesim + npsim))
!   allocate(perturbed_ao_dmat(ntbas(0), ntbas(0), nz, nr_op, nr_freq))
    allocate(perturbed_ao_dmat(ntbas(0), ntbas(0), nz, 2))
    allocate(ao_pmat(ntbas(0), ntbas(0), nz, nr_op))
    allocate(mo_coef(ntbas(0), ntbas(0), nz))

    op_index          = 0
    freq              = 0.0d0
    zo_ao_dmat        = 0.0d0
    perturbed_ao_dmat = 0.0d0
    ao_pmat           = 0.0d0
    mo_coef           = 0.0d0

    static = .false.
    if (nr_freq == 1 .and. freq(1) == 0.0d0) then
      static = .true.
    end if

    op_index(1) = 5

    call get_mo_coef(mo_coef)
    call get_zo_ao_dmat(zo_ao_dmat, mo_coef)

    do iop = 1, 1

      call get_ao_pmat(op_index(iop),         &
                       ao_pmat(:, :, :, iop), &
                       work, kfree, lfree)

      value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, zo_ao_dmat, 1)
      write(*, *) 'D(0)*P = ', value

      call solve_perturbed_ao_dmat(op_index(iop),     &
                                   freq,              &
                                   perturbed_ao_dmat, &
                                   work, kfree, lfree)

      value = ddot(ntbas(0)*ntbas(0)*nz, ao_pmat(:, :, :, iop), 1, perturbed_ao_dmat(:, :, :, 1), 1)
      write(*, *) 'D(1)*P = ', value*2.0d0

    end do

    deallocate(op_index)
    deallocate(freq)
    deallocate(zo_ao_dmat)
    deallocate(perturbed_ao_dmat)
    deallocate(ao_pmat)
    deallocate(mo_coef)

  end subroutine aoresponse_driver

  subroutine get_ao_pmat(op_index, &
                         ao_pmat,  &
                         work, kfree, lfree)
    integer(kind_i4), intent(in)    :: op_index
    real(kind_dp),    intent(out)   :: ao_pmat(:, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer,          intent(inout) :: kfree
    integer,          intent(inout) :: lfree

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    real(kind_dp),    allocatable   :: ptri(:)
    real(kind_dp),    allocatable   :: op1int(:)
    logical,          allocatable   :: first(:)

    logical                         :: file_exists, is_square
    real(kind_dp)                   :: dummy_dp
    integer(kind_i4)                :: file_unit, dummy_i4, idimension

    is_square = .false.
    if (iprptim(op_index) == 0) then
      is_square = .true.
    end if

    if (is_square) then
      idimension = n2bbasx
    else
      idimension = nnbbasx
    end if

    call alloc(ptri, nz*idimension)
    call alloc(op1int, idimension)
    allocate(first(nz))

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

  end subroutine get_ao_pmat

  subroutine solve_perturbed_ao_dmat(op_index,          &
                                     freq,              &
                                     perturbed_ao_dmat, &
                                     work, kfree, lfree)
    integer(kind_i4), intent(in)    :: op_index
    real(kind_dp),    intent(in)    :: freq(:)
    real(kind_dp),    intent(out)   :: perturbed_ao_dmat(:, :, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer(kind_i4), intent(inout) :: kfree
    integer(kind_i4), intent(inout) :: lfree

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"

    integer(kind_i4)                :: i, ifsym, half_length, max_full_length, nr_dmat, nr_freq
    integer(kind_i4)                :: ill_int, isl_int, iss_int, igt_int
    integer(kind_i4), allocatable   :: ibtyp(:), ibtyp_pointer_c(:), ibtyp_pointer_pp(:), ibtyp_pointer_pn(:)
    integer(kind_i4), allocatable   :: isymop(:), ihrmop(:), ifckop(:)
    real(kind_dp),    allocatable   :: eigval(:), eigvec(:)
    real(kind_dp),    allocatable   :: convergence(:)
    real(kind_dp),    allocatable   :: prop_gradient_c(:), prop_gradient_pp(:), prop_gradient_pn(:)
    real(kind_dp),    allocatable   :: response_vector_pp(:), response_vector_pn(:)
    real(kind_dp),    allocatable   :: mo_coef(:, :, :)
    real(kind_dp)                   :: dummy_dp
    logical                         :: skip_pp, skip_pn
    logical                         :: debug_me = .false.
    character(72)                   :: orbital_string(3,2)


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

    skip_pp = .false.
    skip_pn = .false.

    call ostring(orbital_string, &
                 skip_pp,        &
                 skip_pn,        &
                 work, kfree, lfree)

    if (debug_me) then
      write(*, *) 'debug solve_perturbed_ao_dmat'

      do ifsym = 1, nfsym
        write(*, *) 'orbital_string(1, ', ifsym, ') = ', orbital_string(1, ifsym)
        write(*, *) 'orbital_string(2, ', ifsym, ') = ', orbital_string(2, ifsym)
        write(*, *) 'orbital_string(3, ', ifsym, ') = ', orbital_string(3, ifsym)
      end do
    end if


!   initialize symmetry independent variables
!   =========================================

    lineq     = .true.
    fkrmc     = .false.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .false.
    diaghe    = .true.
    iprxrs    = 0
    thcxrs    = 1.0d-5
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

    if (debug_me) then
      write(*, *) 'debug solve_perturbed_ao_dmat'

      write(*, *) 'jsymop = ', jsymop
      write(*, *) 'jtimop = ', jtimop
      write(*, *) 'jopsy  = ', jopsy
    end if

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
      write(*, *) 'debug solve_perturbed_ao_dmat'

      write(*, *)
      write(*, *) 'prop_gradient_pp:'
      call print_q_vector(prop_gradient_pp, nzxopeq)

      write(*, *)
      write(*, *) 'prop_gradient_pn:'
      call print_q_vector(prop_gradient_pn, nzxoppq)
    end if

!   todo: check the norm of the property gradient
!         and if below threshold then skip things


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    call alloc(ibtyp, 2*nredm)
    call alloc(ibtyp_pointer_c,  nredm)
    call alloc(ibtyp_pointer_pp, nredm)
    call alloc(ibtyp_pointer_pn, nredm)

    call alloc(eigval, nexsim)
    call alloc(eigvec, nredm*nexsim)
    call alloc(convergence, nr_freq)

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


!   construct the pp and pn response vectors
!   ========================================

    nesim = 0
    npsim = 0
    ncsim = 0

    if (nzxope > 0) then

      if (static) then
        nesim = nr_freq
      else
        nesim = nr_freq*2
      end if

      half_length     = nzxope
      max_full_length = nzxope*2

      call alloc(response_vector_pp, max_full_length*nz)

      call construct_response_vector('pp',               &
                                     response_vector_pp, &
                                     half_length,        &
                                     nr_freq,            &
                                     ibtyp,              &
                                     ibtyp_pointer_pp,   &
                                     eigvec)

      if (debug_me) then
        call header('response_vector_pp',-1)
        call prbvec(lupri, response_vector_pp, nr_freq, max_full_length)
      end if

    end if

    if (nzxopp > 0) then

      if (static) then
        npsim = nr_freq
      else
        npsim = nr_freq*2
      end if

      half_length     = nzxopp
      max_full_length = nzxopp*2

      call alloc(response_vector_pn, max_full_length*nz)

      call construct_response_vector('pn',               &
                                     response_vector_pn, &
                                     half_length,        &
                                     nr_freq,            &
                                     ibtyp,              &
                                     ibtyp_pointer_pn,   &
                                     eigvec)

      if (debug_me) then
        call header('response_vector_pn',-1)
        call prbvec(lupri, response_vector_pn, nr_freq, max_full_length)
      end if

    end if

    call dealloc(eigval)
    call dealloc(eigvec)


!   construct perturbed ao density matrices
!   =======================================

    call alloc(isymop, ntsim)
    call alloc(ihrmop, ntsim)
    call alloc(ifckop, ntsim)

    call setgmt(isymop,           &
                ihrmop,           &
                ifckop,           &
                ibtyp_pointer_c,  &
                ibtyp_pointer_pp, &
                ibtyp_pointer_pn, &
                ibtyp, 0)

    call dealloc(ibtyp)
    call dealloc(ibtyp_pointer_c)
    call dealloc(ibtyp_pointer_pp)
    call dealloc(ibtyp_pointer_pn)

!   allocate(perturbed_ao_dmat(ntbas(0), ntbas(0), nz, nesim + npsim))
    allocate(mo_coef(ntbas(0), ntbas(0), nz))

    mo_coef = 0.0d0

    call get_mo_coef(mo_coef)

    call dmolitx(perturbed_ao_dmat,  &
                 response_vector_pp, &
                 response_vector_pn, &
                 mo_coef,            &
                 dummy_dp,           &
                 dummy_dp,           &
                 ihrmop,             &
                 nr_dmat,            &
                 work, kfree, lfree)

    deallocate(mo_coef)

    call dealloc(isymop)
    call dealloc(ihrmop)
    call dealloc(ifckop)

    if (nzxope > 0) then
      call dealloc(response_vector_pp)
    end if
    if (nzxopp > 0) then
      call dealloc(response_vector_pn)
    end if

    if (debug_me) then
      write(*, *) 'debug solve_perturbed_ao_dmat'

      do i = 1, nr_dmat
        write(*, *)
        write(*, *) 'perturbed_ao_dmat ', i
        call print_square_q_matrix(perturbed_ao_dmat(:, :, :, i), ntbas(0))
      end do
    end if

  end subroutine solve_perturbed_ao_dmat

  subroutine construct_response_vector(vector_type,     &
                                       response_vector, &
                                       half_length,     &
                                       nr_freq,         &
                                       ibtyp,           &
                                       ibtyp_pointer,   &
                                       eigvec)
    character(*),     intent(in)  :: vector_type
    real(kind_dp),    intent(out) :: response_vector(:)
    integer(kind_i4), intent(in)  :: half_length
    integer(kind_i4), intent(in)  :: nr_freq
    integer(kind_i4), intent(in)  :: ibtyp(2, *)
    integer(kind_i4), intent(in)  :: ibtyp_pointer(*)
    real(kind_dp),    intent(in)  :: eigvec(*)

#include "dcbxrs.h"
#include "dgroup.h"

    real(kind_dp),    allocatable :: response_vector_z(:, :, :)
    real(kind_dp),    allocatable :: response_vector_y(:, :, :)
    integer(kind_i4), allocatable :: ivecs(:)
    real(kind_dp),    allocatable :: buffer(:)

    character(6)                  :: file_name
    integer(kind_i4)              :: length, file_unit, ifreq, iz, nbtyp, ioff

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
        call quit('error in construct_response_vector, unknown vector type')
    end select

    allocate(response_vector_z(half_length, nz, nr_freq))
    allocate(response_vector_y(half_length, nz, nr_freq))
    call alloc(ivecs, nr_freq)
    call alloc(buffer, half_length*nz)

    response_vector_z = 0.0d0
    response_vector_y = 0.0d0

    open(file_unit,                 &
         file   = file_name,        &
         form   = 'unformatted',    &
         access = 'direct',         &
         recl   = 8*half_length*nz, &
         status = 'old')

    call xrsxv1( 1, nbtyp, response_vector_z, &
                    eigvec, nr_freq, ibtyp, ibtyp_pointer, ivecs, buffer)
    call xrsxv1(-1, nbtyp, response_vector_y, &
                    eigvec, nr_freq, ibtyp, ibtyp_pointer, ivecs, buffer)

    close(file_unit, &
          status = 'delete')

    response_vector = 0.0d0

    length = 2*half_length

    do ifreq = 1, nr_freq
      do iz = 1, nz

        ioff = (iz - 1)*length
        call daxpy(half_length, 1.0d0, response_vector_z(1, iz, ifreq), 1, response_vector(1 + ioff), 1)
        call daxpy(half_length, 1.0d0, response_vector_y(1, iz, ifreq), 1, response_vector(1 + ioff), 1)

        ioff = (iz - 1)*length + half_length
        call daxpy(half_length, 1.0d0, response_vector_z(1, iz, ifreq), 1, response_vector(1 + ioff), 1)
        call daxpy(half_length,-1.0d0, response_vector_y(1, iz, ifreq), 1, response_vector(1 + ioff), 1)

        if (ipqtoq(iz, jsymop - 1) > 1) then
          call dscal(half_length, -1.0d0, response_vector(1 + (iz - 1)*length + half_length), 1)
        end if

      end do
    end do

    deallocate(response_vector_z)
    deallocate(response_vector_y)
    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine construct_response_vector

end module aoresponse_main_m
