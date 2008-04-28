module aoresponse_main_m

  use kinds_m
  use scf_density_class
  use memory_allocator

  implicit none

#include "priunit.h"

! type gao_data_t
!   integer(kind_i4)              :: max_derivative
!   integer(kind_i4), allocatable :: gao_center(:)
!   real(kind_dp),    allocatable :: gao(:)
!   real(kind_dp),    allocatable :: gao1(:)
!   real(kind_dp),    allocatable :: gao2(:)
!   real(kind_dp),    allocatable :: buffer(:)
! end type gao_data_t

  public open_response_driver, get_ao_pmat, solve_perturbed_ao_dmat

  private ! everything below

contains














! subroutine open_response_driver(work, kfree, lfree, ao_dmat_fo)
  subroutine open_response_driver(work, kfree, lfree)
#include "dcbbas.h"
#include "dgroup.h"
    real(kind_dp), intent(inout)         :: work(lfree)
    integer,       intent(inout)         :: kfree
    integer,       intent(inout)         :: lfree
!   real(kind_dp), intent(in)            :: ao_dmat_fo(ntbas(0)*ntbas(0)*nz)

    real(kind_dp),           allocatable :: pmat(:, :, :)
    real(kind_dp),    allocatable        :: ao_dmat_zo(:, :, :)
    real(kind_dp),    allocatable        :: mo_coef(:, :, :)
    real(kind_dp)                        :: ddot
    real(kind_dp)                        :: dip
    integer                              :: i



      call init_mo_coef(mo_coef)
      call get_mo_coef(mo_coef)

      call init_ao_dmat(ao_dmat_zo)
      call get_ao_dmat_zo(ao_dmat_zo, mo_coef)

      call del_mo_coef(mo_coef)



!   write(*, *) 'calling PRPXL1_clone'
!   call PRPXL1_clone(WORK,KFREE,LFREE)



    do i = 1, 10
       call init_ao_dmat(pmat)
       call get_ao_pmat(i, pmat, work, kfree, lfree)
       dip = ddot(ntbas(0)*ntbas(0)*nz, pmat, 1, ao_dmat_zo, 1)
!      dip = ddot(ntbas(0)*ntbas(0)*nz, pmat, 1, ao_dmat_fo, 1)
!      call print_square_q_matrix(pmat, ntbas(0))
       write(*, *) 'blub', i, dip

       call del_ao_dmat(pmat)
    end do

    call del_ao_dmat(ao_dmat_zo)




  end subroutine open_response_driver







  subroutine get_ao_pmat(property_index, &
                         ao_pmat,        &
                         work, kfree, lfree)
    integer(kind_i4), intent(in)    :: property_index
    real(kind_dp),    intent(out)   :: ao_pmat(:, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer,          intent(inout) :: kfree
    integer,          intent(inout) :: lfree

#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbxpr.h"

    logical                         :: file_exists, is_square
    real(kind_dp)                   :: dummy_dp
    integer(kind_i4)                :: file_unit, dummy_i4, idimension
    real(kind_dp),    allocatable   :: ptri(:)
    real(kind_dp),    allocatable   :: op1int(:)
    logical,          allocatable   :: first(:)

    is_square = .false.
    if (iprptim(property_index) == 0) then
      is_square = .true.
    end if

    if (is_square) then
      idimension = n2bbasx
    else
      idimension = nnbbasx
    end if

    allocate(ptri(nz*idimension))
    allocate(op1int(idimension))
    allocate(first(nz))

    file_unit = lu1int

    file_exists = .false.
    inquire(file = 'AOPROPER', exist = file_exists)
    if (.not. file_exists) then
      call quit('file AOPROPER not found')
    end if

    open(file_unit,          &
         file = 'AOPROPER',  &
         status = 'unknown', &
         form = 'unformatted')

    rewind file_unit

    call prpex2(property_index, &
                dummy_dp,       &
                .false.,        &
                idimension,     &
                ptri,           &
                ao_pmat,        &
                op1int,         &
                first,          &
                work, kfree, lfree, 0)

    close(file_unit, &
          status = 'keep')

    deallocate(ptri)
    deallocate(op1int)
    deallocate(first)

  end subroutine get_ao_pmat

  subroutine solve_perturbed_ao_dmat(property_index, &
!                                    frequencies,    &
!                                    ao_dmat,        &
                                     work, kfree, lfree)
    integer(kind_i4), intent(in)    :: property_index
!   real(kind_dp),    intent(in)    :: frequencies(:)
    real(kind_dp)                   :: frequencies(1)
    real(kind_dp),    allocatable   :: ao_dmat(:, :, :, :)
    real(kind_dp),    intent(inout) :: work(lfree)
    integer(kind_i4), intent(inout) :: kfree
    integer(kind_i4), intent(inout) :: lfree

#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"

    integer(kind_i4)                :: ifsym, half_length, max_full_length, nr_dmat
    integer(kind_i4)                :: ill_int, isl_int, iss_int, igt_int
    integer(kind_i4), allocatable   :: ibtyp(:), ibtyp_pointer_c(:), ibtyp_pointer_pp(:), ibtyp_pointer_pn(:)
    integer(kind_i4), allocatable   :: isymop(:), ihrmop(:), ifckop(:)
    real(kind_dp),    allocatable   :: eigenvalues(:), eigenvectors(:)
    real(kind_dp),    allocatable   :: convergence(:)
    real(kind_dp),    allocatable   :: prop_gradient_c(:), prop_gradient_pp(:), prop_gradient_pn(:)
    real(kind_dp),    allocatable   :: response_vector_pp(:), response_vector_pn(:)
    real(kind_dp),    allocatable   :: mo_coef(:, :, :)
    real(kind_dp)                   :: dummy_dp
    logical                         :: skip_pp, skip_pn
    character(72)                   :: orbital_string(3,2)



    logical                         :: debug_me
!   character(16)                   :: operator_name
    integer(kind_i4)                :: i
!   integer(kind_i4)                :: ITRXPR, IND1OP, LAQROP,NAQROP, ISYXPR


    debug_me = .true.

    frequencies(1) = 0.0d0
!   property_index = property_index - 2




!   set common block dcbidx
!   =======================

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

    jsymop = iprpsym(property_index)
    jtimop = iprptim(property_index)
    jopsy  = jbtof(jsymop - 1, 1)
    nfreq  = size(frequencies)
    nexsim = nfreq
    nexstv = nfreq
    nexcnv = nfreq
    ncred  = 0
    nered  = 0
    npred  = 0
    nzred  = ncred + nered + npred

    if (debug_me) then
      write(*, *) 'debug solve_perturbed_ao_dmat'

      write(*, *) 'jsymop = ', jsymop
      write(*, *) 'jtimop = ', jtimop
      write(*, *) 'jopsy  = ', jopsy
    end if

    static = .false.
    if (nfreq == 1 .and. frequencies(1) == 0.0d0) then
      static = .true.
    end if

    call xrspar(orbital_string, &
                skip_pp,        &
                skip_pn,        &
                work,kfree,lfree, 0)


!   get property gradient
!   =====================

    nzconfq = 0

    allocate(prop_gradient_c(nzconfq))
    allocate(prop_gradient_pp(nzxopeq))
    allocate(prop_gradient_pn(nzxoppq))

    prop_gradient_c  = 0.0d0
    prop_gradient_pp = 0.0d0
    prop_gradient_pn = 0.0d0

    call pamgrd(property_index,   &
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

    allocate(ibtyp(2*nredm))
    allocate(ibtyp_pointer_c(nredm))
    allocate(ibtyp_pointer_pp(nredm))
    allocate(ibtyp_pointer_pn(nredm))

    allocate(eigenvalues(nexsim))
    allocate(eigenvectors(nredm*nexsim))
    allocate(convergence(nfreq))

    call dcopy(nfreq, frequencies, 1, eigenvalues, 1)

    call xrsctl(prop_gradient_c,      &
                prop_gradient_pp,     &
                prop_gradient_pn,     &
                ibtyp,                &
                ibtyp_pointer_c,      &
                ibtyp_pointer_pp,     &
                ibtyp_pointer_pn,     &
                convergence,          &
                eigenvalues,          &
                eigenvectors,         &
                work,kfree,lfree)

    deallocate(prop_gradient_c)
    deallocate(prop_gradient_pp)
    deallocate(prop_gradient_pn)

    deallocate(convergence)


!   construct the pp and pn response vectors
!   ========================================

    nesim = 0
    npsim = 0
    ncsim = 0

    if (nzxope > 0) then

      if (static) then
        nesim = nfreq
      else
        nesim = nfreq*2
      end if

      half_length     = nzxope
      max_full_length = nzxope*2

      allocate(response_vector_pp(max_full_length*nz))

      call construct_response_vector('pp',                 &
                                     property_index,       &
                                     response_vector_pp,   &
                                     half_length,          &
                                     ibtyp,                &
                                     ibtyp_pointer_pp,     &
                                     eigenvalues,          &
                                     eigenvectors)

      if (debug_me) then
        call header('response_vector_pp',-1)
        call prbvec(lupri, response_vector_pp, nfreq, max_full_length)
      end if

    end if

    if (nzxopp > 0) then

      if (static) then
        npsim = nfreq
      else
        npsim = nfreq*2
      end if

      half_length     = nzxopp
      max_full_length = nzxopp*2

      allocate(response_vector_pn(max_full_length*nz))

      call construct_response_vector('pn',                 &
                                     property_index,       &
                                     response_vector_pn,   &
                                     half_length,          &
                                     ibtyp,                &
                                     ibtyp_pointer_pn,     &
                                     eigenvalues,          &
                                     eigenvectors)

      if (debug_me) then
        call header('response_vector_pn',-1)
        call prbvec(lupri, response_vector_pn, nfreq, max_full_length)
      end if

    end if

    deallocate(eigenvalues)
    deallocate(eigenvectors)


!   construct perturbed ao density matrices
!   =======================================

    allocate(isymop(ntsim))
    allocate(ihrmop(ntsim))
    allocate(ifckop(ntsim))

    call setgmt(isymop,           &
                ihrmop,           &
                ifckop,           &
                ibtyp_pointer_c,  &
                ibtyp_pointer_pp, &
                ibtyp_pointer_pn, &
                ibtyp, 0)

    deallocate(ibtyp)
    deallocate(ibtyp_pointer_c)
    deallocate(ibtyp_pointer_pp)
    deallocate(ibtyp_pointer_pn)

    allocate(ao_dmat(ntbas(0), ntbas(0), nz, nesim + npsim))

    call init_mo_coef(mo_coef)
    call get_mo_coef(mo_coef)

!   ao_dmat order: 1...nesim, 1...npsim
    call dmolitx(ao_dmat,            &
                 response_vector_pp, &
                 response_vector_pn, &
                 mo_coef,            &
                 dummy_dp,           &
                 dummy_dp,           &
                 ihrmop,             &
                 nr_dmat,            &
                 work, kfree, lfree)

    if (debug_me) then
      write(*, *) 'debug solve_perturbed_ao_dmat'

      do i = 1, nesim + npsim
        write(*, *)
        write(*, *) 'ao_dmat ', i
        call print_square_q_matrix(ao_dmat(:, :, :, i), ntbas(0))
      end do
    end if

    call del_mo_coef(mo_coef)

    if (nzxope > 0) then
      deallocate(response_vector_pp)
    end if
    if (nzxopp > 0) then
      deallocate(response_vector_pn)
    end if

    deallocate(ao_dmat)

    deallocate(isymop)
    deallocate(ihrmop)
    deallocate(ifckop)

  end subroutine solve_perturbed_ao_dmat

  subroutine construct_response_vector(vector_type,          &
                                       property_index,       &
                                       response_vector,      &
                                       half_length,          &
                                       ibtyp,                &
                                       ibtyp_pointer,        &
                                       eigenvalues,          &
                                       eigenvectors)
#include "dcbxrs.h"
#include "dgroup.h"

    character(*),     intent(in)  :: vector_type
    integer(kind_i4), intent(in)  :: property_index
    real(kind_dp),    intent(out) :: response_vector(:)
    integer(kind_i4), intent(in)  :: half_length
    integer(kind_i4), intent(in)  :: ibtyp(2, *)
    integer(kind_i4), intent(in)  :: ibtyp_pointer(*)
    real(kind_dp),    intent(in)  :: eigenvalues(*)
    real(kind_dp),    intent(in)  :: eigenvectors(*)

    real(kind_dp),    allocatable :: response_vector_z(:, :, :)
    real(kind_dp),    allocatable :: response_vector_y(:, :, :)
    integer(kind_i4), allocatable :: ivecs(:)
    real(kind_dp),    allocatable :: buffer(:)

    character(6)                  :: file_name
    integer(kind_i4)              :: length, file_unit, ifreq, iz, nbtyp, nbred, ioff

    select case (vector_type)
      case ('pp')
        nbtyp = 1
        nbred  = nered
        file_unit = luboe
        file_name = 'PAMBOE'
      case ('pn')
        nbtyp = 2
        nbred  = npred
        file_unit = lubop
        file_name = 'PAMBOP'
      case default
        call quit('error in construct_response_vector, unknown vector type')
    end select

    call alloc((/response_vector_z/), half_length*nz*nfreq)
    call alloc((/response_vector_y/), half_length*nz*nfreq)
    call alloc(ivecs, nfreq)
    call alloc(buffer, half_length*nz)

    response_vector_z = 0.0d0
    response_vector_y = 0.0d0

    open(file_unit,               &
         file = file_name,        &
         form = 'unformatted',    &
         access = 'direct',       &
         recl = 8*half_length*nz, &
         status = 'old')

    call xrsxv1( 1, nbtyp, response_vector_z, &
                    eigenvectors, nfreq, ibtyp, ibtyp_pointer, ivecs, buffer)
    call xrsxv1(-1, nbtyp, response_vector_y, &
                    eigenvectors, nfreq, ibtyp, ibtyp_pointer, ivecs, buffer)

    close(file_unit, &
          status = 'delete')

    response_vector = 0.0d0

    length = 2*half_length

    do ifreq = 1, nfreq
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

    call dealloc((/response_vector_z/))
    call dealloc((/response_vector_y/))
    call dealloc(ivecs)
    call dealloc(buffer)

  end subroutine construct_response_vector

end module aoresponse_main_m
