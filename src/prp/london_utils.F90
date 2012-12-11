module london_utils

  use memory_allocator 
  use dirac_cfg
  use dft_cfg
  use num_grid_gen
  use dirac_interface

  implicit none

  public test_new_routines
  public check_symmetry 
  public combine_kappa_T
  public get_T_ai 
  public get_m_ai 
  public get_m_ij 
  public get_T_ji 
  public get_T_ji2 
  save

  private

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbfir.h"
#include "dcbham.h"
#include "iratdef.h"
#include "cbihr2.h"
#include "dcbdhf.h"
#include "shells.h"

#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
  dimension     istat(mpi_status_size)
#endif /* ifdef VAR_MPI */



  integer, parameter :: max_nr_xyz_comp = 3
  integer            ::     nr_xyz_comp

  integer, parameter :: tbmo_unit = 71
  logical            :: debug_me

contains

  subroutine get_m_ai(m, sum_ep_rn)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     logical, intent(in)    :: sum_ep_rn

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j

!    ---------------------------------------------------------------------------

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ai: entering T:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if

     
     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)

        debug_me = .false.
        if (debug_me) then
           write(*, *) 'ifrp, ip = ', ifrp, ip(ifrp)
           write(*, *) 'ifrp, ap = ', ifrp, ap(ifrp)
           write(*, *) 'ifrp, ie = ', ifrp, ie(ifrp)
           write(*, *) 'ifrp, ae = ', ifrp, ae(ifrp)
           write(*, *) 'ifrp, nish = ', ifrp, nish(ifrp)
           write(*, *) 'ifrp, npsh = ', ifrp, npsh(ifrp)
           write(*, *) 'ifrp, norb = ', ifrp, norb(ifrp)
        end if
     end do

!    'ai' blocks (zero-out 'ij' blocks)
     do iz = 1, nz
        do ifrp = 1, nfsym

!          zero-out 'ij' blocks
           do j = ie(ifrp), ae(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
           end do

!          zero-out 'aa' blocks
           do j = ip(ifrp), ie(ifrp) - 1
              do i = ip(ifrp), ie(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 m(i, j, iz) = 0.0d0
              end do
           end do
           do j = ae(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 m(i, j, iz) = 0.0d0
              end do
           end do

           if (.not. sum_ep_rn) then
              do j = ie(ifrp), ae(ifrp) - 1
                 do i = ip(ifrp), ie(ifrp) -1
                    m(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ip(ifrp), ie(ifrp) -1
                 do i = ie(ifrp), ae(ifrp) - 1
                    m(i, j, iz) = 0.0d0
                 end do
              end do
           end if

        end do
     end do

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ai: T_ai on output, sum_ep_rn = :', sum_ep_rn
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if


  end subroutine


  subroutine get_T_ai(T, T_ai, sum_ep_rn)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: T(norbt, norbt, nz)
     real(8), intent(inout) :: T_ai(norbt, norbt, nz)
     logical, intent(in)    :: sum_ep_rn

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j

!    on input:  T    - 'full' connection matrix in MO basis
!                      (note: for 1 component of magnetic field)
!    on output: T_ai - take only inactive-secondary blocks of T 
!                      and zero-out the rest
!    ---------------------------------------------------------------------------
     
     call dcopy(norbt*norbt*nz, T, 1, T_ai, 1)

     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)

        debug_me = .false.
        if (debug_me) then
           write(*, *) 'ifrp, ip = ', ifrp, ip(ifrp)
           write(*, *) 'ifrp, ap = ', ifrp, ap(ifrp)
           write(*, *) 'ifrp, ie = ', ifrp, ie(ifrp)
           write(*, *) 'ifrp, ae = ', ifrp, ae(ifrp)
           write(*, *) 'ifrp, nish = ', ifrp, nish(ifrp)
           write(*, *) 'ifrp, npsh = ', ifrp, npsh(ifrp)
           write(*, *) 'ifrp, norb = ', ifrp, norb(ifrp)
        end if
     end do

!    'ai' blocks (zero-out 'ij' blocks)
     do iz = 1, nz
        do ifrp = 1, nfsym

!          zero-out 'ij' blocks
           do j = ie(ifrp), ae(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do

!          zero-out 'aa' blocks
           do j = ip(ifrp), ie(ifrp) - 1
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do
           do j = ae(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ai(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ai(i, j, iz) = 0.0d0
              end do
           end do

           if (.not. sum_ep_rn) then
              do j = ie(ifrp), ae(ifrp) - 1
                 do i = ip(ifrp), ie(ifrp) -1
                    T_ai(i, j, iz) = 0.0d0
                 end do
              end do
              do j = ip(ifrp), ie(ifrp) -1
                 do i = ie(ifrp), ae(ifrp) - 1
                    T_ai(i, j, iz) = 0.0d0
                 end do
              end do
           end if

        end do
     end do

  end subroutine


  subroutine get_T_ji2(T, T_ji, sum_ep_rn)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: T(norbt, norbt, nz)
     real(8), intent(inout) :: T_ji(norbt, norbt, nz)
     logical, intent(in)    :: sum_ep_rn

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j

!    on input:  T    - 'full' connection matrix in MO basis
!                      (note: for 1 component of magnetic field)
!    on output: T_ji - take only inactive-inactive blocks of T
!                      and zero-out the rest
! gosia todo: adapt for symmetry!
!    ---------------------------------------------------------------------------
     
     call dzero(T_ji, norbt*norbt*nz)

     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)
     end do

     do iz = 1, nz
        do ifrp = 1, nfsym

           do j = ie(ifrp), ae(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 T_ji(i, j, iz) = T(i, j, iz)
              end do
           end do

        end do
     end do

  end subroutine


  subroutine get_m_ij(m)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j
!    ---------------------------------------------------------------------------
     
     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)

        debug_me = .false.
        if (debug_me) then
           write(*, *) 'ifrp, ip = ', ifrp, ip(ifrp)
           write(*, *) 'ifrp, ap = ', ifrp, ap(ifrp)
           write(*, *) 'ifrp, ie = ', ifrp, ie(ifrp)
           write(*, *) 'ifrp, ae = ', ifrp, ae(ifrp)
           write(*, *) 'ifrp, nish = ', ifrp, nish(ifrp)
           write(*, *) 'ifrp, npsh = ', ifrp, npsh(ifrp)
           write(*, *) 'ifrp, norb = ', ifrp, norb(ifrp)
        end if
     end do

     do iz = 1, nz
        do ifrp = 1, nfsym

           do j = ip(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 m(i, j, iz) = 0.0d0
              end do
           end do
           do j = ip(ifrp), ie(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
           end do
           do j = ae(ifrp), norb(ifrp)
              do i = ie(ifrp), ae(ifrp) - 1
                 m(i, j, iz) = 0.0d0
              end do
           end do

        end do
     end do

  end subroutine

  subroutine get_T_ji(T, T_ji, sum_ep_rn)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: T(norbt, norbt, nz)
     real(8), intent(inout) :: T_ji(norbt, norbt, nz)
     logical, intent(in)    :: sum_ep_rn

     integer :: ip(2), ap(2), ie(2), ae(2)
     integer :: ifrp, iz, i, j

!    on input:  T    - 'full' connection matrix in MO basis
!                      (note: for 1 component of magnetic field)
!    on output: T_ji - take only inactive-inactive blocks of T
!                      and zero-out the rest
!    ---------------------------------------------------------------------------
     
     call dcopy(norbt*norbt*nz, T, 1, T_ji, 1)

     do ifrp = 1, nfsym
        ip(ifrp) = iorb(ifrp) + 1
        ap(ifrp) = ip(ifrp) + nish(ifrp)
        ie(ifrp) = ip(ifrp) + npsh(ifrp)
        ae(ifrp) = ie(ifrp) + nish(ifrp)

        debug_me = .false.
        if (debug_me) then
           write(*, *) 'ifrp, ip = ', ifrp, ip(ifrp)
           write(*, *) 'ifrp, ap = ', ifrp, ap(ifrp)
           write(*, *) 'ifrp, ie = ', ifrp, ie(ifrp)
           write(*, *) 'ifrp, ae = ', ifrp, ae(ifrp)
           write(*, *) 'ifrp, nish = ', ifrp, nish(ifrp)
           write(*, *) 'ifrp, npsh = ', ifrp, npsh(ifrp)
           write(*, *) 'ifrp, norb = ', ifrp, norb(ifrp)
        end if
     end do

     do iz = 1, nz
        do ifrp = 1, nfsym

           do j = ip(ifrp), norb(ifrp)
              do i = ip(ifrp), ie(ifrp) - 1
                 T_ji(i, j, iz) = 0.0d0
              end do
              do i = ae(ifrp), norb(ifrp)
                 T_ji(i, j, iz) = 0.0d0
              end do
           end do
           do j = ip(ifrp), ie(ifrp) - 1
              do i = ie(ifrp), ae(ifrp) - 1
                 T_ji(i, j, iz) = 0.0d0
              end do
           end do
           do j = ae(ifrp), norb(ifrp)
              do i = ie(ifrp), ae(ifrp) - 1
                 T_ji(i, j, iz) = 0.0d0
              end do
           end do

        end do
     end do

  end subroutine

  subroutine test_new_routines(T, irep)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: T(norbt, norbt, nz)
     integer, intent(in)    :: irep

     real(8), allocatable :: T_ai(:, :, :)
     real(8), allocatable :: T_ji(:, :, :)

     integer :: i, j, iz
     real(8) :: a
!    ---------------------------------------------------------------------------

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'test_new_routines: entering T:'
        call prqmat(T, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, irep), lupri)
     end if

     call alloc(T_ai, norbt, norbt, nz)
     call get_T_ai(T, T_ai, doeprn)

     call alloc(T_ji, norbt, norbt, nz)
     call get_T_ji(T, T_ji, doeprn)

     call dzero(T, norbt*norbt*nz)

     do iz = 1, nz
        do j = 1, norbt
           do i = 1, norbt
              T(i, j, iz) = T_ai(i, j, iz) + T_ji(i, j, iz)
           end do
        end do
     end do

     call dealloc(T_ai)
     call dealloc(T_ji)

  end subroutine


  subroutine combine_kappa_T(kappa, length, nbtyp, lab1)
!    ---------------------------------------------------------------------------
     real(8), intent(inout)   :: kappa(length, nz)
     integer, intent(in)      :: length, nbtyp
     character*16, intent(in) :: lab1

     real(8), allocatable   :: T(:, :, :)
     real(8), allocatable   :: T_ai_vec(:, :)
     integer, allocatable   :: orbital_rotation_indices(:, :)
     integer :: which_icomp, chosen_icomp, irep, iz, i, j
     logical :: sum_ep_rn, ep, ee

!    this subroutine does:
!    kappa_{ai}^B := kappa_{ai}^B - T_{ai}^B

!    where: kappa_{ai}^B - elements of solution vector
!           T_{ai}^B     - elements of connection matrix (occ-virt blocks)
!                          the 'full' connection matrix is read from TBMO file

!    it is used to calculate the NMR shielding tensor
!    from Eq. 40-42 from JCP 131, 124119
!    ---------------------------------------------------------------------------

!    which component of magnetic field is that:
!    ------------------------------------------
     if (lab1 == 'LAO-XRM1H1') chosen_icomp = 1
     if (lab1 == 'LAO-YRM1H1') chosen_icomp = 2
     if (lab1 == 'LAO-ZRM1H1') chosen_icomp = 3

!    prepare info about solution vector:
!    -----------------------------------
     ep = .false.
     ee = .false.
     if (nbtyp == 1) then
       ee = .true.
       ep = .false.
     else if (nbtyp == 2) then
       ep = .true.
       ee = .false.
     else
        call quit('error in combine_kappa_T')
     end if

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_T: solution vector for operator = ', lab1
        write(*, *) 'which class of rotations: ee, ep = ', ee, ep
        write(*, *) 'length, nz = ', length, nz
        write(*, *) 'entering kappa:'
        call prbvec(lupri, kappa , 1, length)
     end if

!    get information about orbital rotation indices
!    ----------------------------------------------
     call alloc(orbital_rotation_indices, 2, length)
     call fill_orbital_rotation_indices_vector(ee,     &
                                               ep,     &
                                               length, & 
                                               orbital_rotation_indices)
     

!    read connection matrix from TBMO file:
!    --------------------------------------
     call alloc(T, norbt, norbt, nz)

     open(tbmo_unit,                 &
          file   = 'TBMO',           &
          status = 'old',            &
          form   = 'unformatted',    &
          access = 'sequential',     &
          position = 'rewind',       &
          action = 'read')
    
1    continue
    
     read(tbmo_unit) nr_xyz_comp
     read(tbmo_unit) which_icomp
     read(tbmo_unit) irep
     read(tbmo_unit) sum_ep_rn
     read(tbmo_unit) T

     if (nr_xyz_comp > 1) then
!       when calculating full shielding tensor
!       there will be three matrices on TBMO file: T^Bx, T^By and T^Bz
!       here we need to read the requested one:
        if (chosen_icomp == 2) then
           if (which_icomp == 2) then
              go to 2
           else
              go to 1
           end if
        else if (chosen_icomp == 3) then
           if (which_icomp == 3) then
              go to 2
           else
              go to 1
           end if
        end if
     end if
     
2    continue

     close(tbmo_unit, status = 'keep')


!    leave only occ-virt part of T:
!    ------------------------------
!gosia: todo fix sum_ep_rn: by default it should be fixed to \kappa's
!*) if e-p rotations are allowed, then sum_ep_rn=true,
!*) if only e-e rotations are allowed, then sum_ep_rn=false
!
     call get_m_ai(T, sum_ep_rn)

!    T_ai is a matrix, so now let's transform it to a vector
!    take care of the correct order -> should be the same as
!    on solution vector:
     call alloc(T_ai_vec, length, nz)
     T_ai_vec = 0.0d0

     irep = 0
     call collect_matrix_to_vector(length, &
                                   orbital_rotation_indices, &
!                                   h, &
                                   T_ai_vec, &
                                   T, &
                                   irep)

     call dealloc(T)
     call dealloc(orbital_rotation_indices)

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'T_ai_vec,ee, ep = ', ee, ep
        do iz = 1, nz
          do i = 1, length
             write(*, *) 'iz,i,T_ai_vec(i, iz) = ', iz, i, T_ai_vec(i, iz)
          end do
        end do
     end if

!    finally substract: kappa := kappa - T
!    -------------------------------------
     do iz = 1, nz
        do i = 1, length
           kappa(i, iz) = kappa(i, iz) - T_ai_vec(i, iz)
        end do
     end do

     call dealloc(T_ai_vec)

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'kappa on output'
        call prqmat(kappa(1,1), length, 1, length, 1, nz,  &
                    ipqtoq(1,irep), lupri)
     end if

  end subroutine


  subroutine fill_orbital_rotation_indices_vector(ee,     &
                                                  ep,     &
                                                  length, &
                                                  orbital_rotation_indices)
!    ---------------------------------------------------------------------------
     logical, intent(in)    :: ee, ep
     integer, intent(in)    :: length
     integer, intent(inout) :: orbital_rotation_indices(2, *)
     integer :: i1, i, j, isym
     integer :: ip(2), ap(2), ie(2), ae(2)
!    ---------------------------------------------------------------------------

     do isym = 1, nfsym
        ip(isym) = iorb(isym) + 1
        ap(isym) = ip(isym) + nish(isym)
        ie(isym) = ip(isym) + npsh(isym)
        ae(isym) = ie(isym) + nish(isym)
     end do

     orbital_rotation_indices(1:2, 1:length) = 0
!    fill orbital_rotation_indices vector:
     i1 = 0
     do isym = 1, nfsym

        if (ee) then
!          inactive:
           do j = ie(isym), ae(isym) - 1
!             secondary:
              do i = ae(isym), norb(isym)
                 i1 = i1 + 1
                 orbital_rotation_indices(1, i1) = j 
                 orbital_rotation_indices(2, i1) = i
              end do
           end do
        else if (ep) then
!          inactive:
           do j = ie(isym), ae(isym) - 1
!             secondary:
              do i = ip(isym), ie(isym) - 1
                 i1 = i1 + 1
                 orbital_rotation_indices(1, i1) = j 
                 orbital_rotation_indices(2, i1) = i
              end do
           end do
        else
           write(*, *) 'error: ee and ep are false!'
           call quit('error!')
        end if
     end do

     debug_me = .false.
     if (debug_me) then
        do i = 1, i1
           write(*, *) 'i, orbital_rotation_indices(1, i), orbital_rotation_indices(2, i) = ', &
                        i, orbital_rotation_indices(1, i), orbital_rotation_indices(2, i) 
        end do
     end if

  end subroutine


  subroutine collect_matrix_to_vector(length,         &
                            orbital_rotation_indices, &
!                            h,                        &
                            vector,                   &
                            matrix,                   &
                            irep)

!   subroutine "opposite" to scatter_vector in openrsp/openrsp_interface_response.F90
!   it makes a vector from a matrix
!   it takes indices from orbital_rotation_indices 
!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
!    real(8), intent(in)  :: h
    real(8), intent(out) :: vector(length, nz)
    real(8), intent(in)  :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    debug_me = .false.
    if (debug_me) then
      write(*, *) 'collect_matrix_to_vector: matrix'
      call prqmat(matrix, norbt, norbt, norbt, norbt, nz,  &
                  ipqtoq(1,irep), lupri)
    end if

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

!       gosia todo: it looks like this simple line is working,
!       but anyway i will have to check it for all symmetries
        vector(is, iz) = matrix(s, i, iz) 

      end do
    end do
 
    debug_me = .false.
    if (debug_me) then
      write(*, *) 'collect_matrix_to_vector: vector'
      call prqmat(vector, length, 1, length, 1, nz,  &
                  ipqtoq(1,irep), lupri)
    end if

  end subroutine


  subroutine check_symmetry(mat, ndim, nz)
    integer, intent(in) :: ndim, nz
    real(8), intent(in) :: mat(ndim, ndim, *)
    real(8) :: a
    integer :: i, j, iz


    a = 0.0d0
    do iz = 1, nz
      do i = 1, ndim - 1
        do j = i + 1, ndim
          a = a + (dabs(mat(i, j, iz)) - dabs(mat(j, i, iz)))
        end do
      end do
    end do
    a = a/dfloat(((ndim*ndim - ndim)/2)*nz)
    write(*, *) 'symmetry mat: ', a

  end subroutine

end module
