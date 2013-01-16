module london_utils

  use memory_allocator 
  use dirac_cfg
  use dft_cfg
  use num_grid_gen
  use dirac_interface

  implicit none

!  public test_new_routines
  public check_symmetry 
  public combine_kappa_T
  public get_m_ai 
  public get_m_ij 
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
#include "dcbxlr.h"

#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
  dimension     istat(mpi_status_size)
#endif /* ifdef VAR_MPI */



  integer, parameter :: max_nr_xyz_comp = 3
  integer            ::     nr_xyz_comp

  integer, parameter :: tbmo_unit = 71
  logical            :: debug_me

! info about virtual orbitals:
  logical :: virt_e, virt_p
  integer :: length_virt_p(2), length_virt_e(2), length_virt(2), length_virt_total

contains

  subroutine get_m_ai(m, ee, ep)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     logical, intent(in) :: ee, ep
     real(8), allocatable   :: m_cp(:,:,:)
     integer, allocatable :: indices_virt(:)
     integer, allocatable :: indices_occ(:)
     integer :: ifrp, iz, i, j
!    ---------------------------------------------------------------------------

!    get absolute indices of occuppied orbitals:
!    -------------------------------------------
     call alloc(indices_occ, nocct)
     indices_occ = 0

     call get_occ_orbitals_indices(indices_occ)

!    get absolute indices of virtual orbitals:
!    -----------------------------------------
     call virtual_orbital_indices_setup() 

     call alloc(indices_virt, length_virt_total)
     indices_virt = 0

     call get_virt_orbitals_indices(indices_virt, ee, ep)


!    copy only needed elements of matrix m to helper matrix m_cp
!    -----------------------------------------------------------

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do iz = 1, nz

        if (nfsym == 2) then
           do j = 1, length_virt_total
              do i = 1, nocct
                 if (indices_occ(i) .ne. 0 .and. indices_virt(j) .ne. 0) then
                    m_cp(indices_occ(i), indices_virt(j), iz) = &
                       m(indices_occ(i), indices_virt(j), iz)
                    m_cp(indices_virt(j), indices_occ(i), iz) = &
                       m(indices_virt(j), indices_occ(i), iz)
                 end if
              end do
           end do

        else  !nfsym = 1

           do j = 1, length_virt(1)
              do i = 1, nocc(1)
                 if (indices_occ(i) .ne. 0 .and. indices_virt(j) .ne. 0) then
                    m_cp(indices_occ(i), indices_virt(j), iz) = &
                       m(indices_occ(i), indices_virt(j), iz)
                    m_cp(indices_virt(j), indices_occ(i), iz) = &
                       m(indices_virt(j), indices_occ(i), iz)
                 end if
              end do
           end do

        end if

     end do

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ai, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m-ai, m_cp, ai only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if

!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)

     call dealloc(indices_occ)
     call dealloc(indices_virt)

  end subroutine


  subroutine get_m_ij(m)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     real(8), allocatable :: m_cp(:,:,:)

     integer :: ifrp, iz, i, j
     integer, allocatable :: indices_occ(:)
!    ---------------------------------------------------------------------------

!    get absolute indices of occuppied orbitals:
!    -------------------------------------------
     call alloc(indices_occ, nocct)
     indices_occ = 0

     call get_occ_orbitals_indices(indices_occ)

!    copy only needed elements of matrix m to helper matrix m_cp
!    -----------------------------------------------------------


     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do iz = 1, nz

        if (nfsym == 2) then
           do j = 1, nocct
              do i = 1, nocct
                 if (indices_occ(i) .ne. 0 .and. indices_occ(j) .ne. 0)then
                    m_cp(indices_occ(i), indices_occ(j), iz) = &
                       m(indices_occ(i), indices_occ(j), iz)
                 end if
              end do
           end do

        else ! nfsym = 1
 
           do j = 1, nocc(1)
              do i = 1, nocc(1)
                 if (indices_occ(i) .ne. 0 .and. indices_occ(j) .ne. 0)then
                    m_cp(indices_occ(i), indices_occ(j), iz) = &
                       m(indices_occ(i), indices_occ(j), iz)
                 end if
              end do
           end do

        end if

     end do

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ij, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m_ij, m_cp, ij only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if


!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)

     call dealloc(indices_occ)


  end subroutine


!  subroutine test_new_routines(T, irep)
!!    ---------------------------------------------------------------------------
!     real(8), intent(inout) :: T(norbt, norbt, nz)
!     integer, intent(in)    :: irep
!
!     real(8), allocatable :: T_ai(:, :, :)
!     real(8), allocatable :: T_ji(:, :, :)
!
!     integer :: i, j, iz
!     real(8) :: a
!!    ---------------------------------------------------------------------------
!
!     debug_me = .false.
!     if (debug_me) then
!        write(*, *) 'test_new_routines: entering T:'
!        call prqmat(T, &
!                    norbt, norbt, norbt, norbt, &
!                    nz, ipqtoq(1, irep), lupri)
!     end if
!
!     call alloc(T_ai, norbt, norbt, nz)
!     call get_T_ai(T, T_ai, doeprn)
!
!     call alloc(T_ji, norbt, norbt, nz)
!     call get_T_ji(T, T_ji, doeprn)
!
!     call dzero(T, norbt*norbt*nz)
!
!     do iz = 1, nz
!        do j = 1, norbt
!           do i = 1, norbt
!              T(i, j, iz) = T_ai(i, j, iz) + T_ji(i, j, iz)
!           end do
!        end do
!     end do
!
!     call dealloc(T_ai)
!     call dealloc(T_ji)
!
!  end subroutine


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

! gosia: todo: now it works only if we always add reorthonormalization contribution
! (run routines in london_reorth), as there the TBMO file is created...
!  i'll have to rewrite it to make it independent on that file
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
!
     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_t: full T, ee, ep = ', ee, ep
        call prqmat(T, norbt, norbt, norbt, norbt, nz,  &
                    ipqtoq(1,irep), lupri)
     end if

     call get_m_ai(T, ee, ep)

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_t: T_ai,ee, ep = ', ee, ep
        call prqmat(T, norbt, norbt, norbt, norbt, nz,  &
                    ipqtoq(1,irep), lupri)
     end if


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
        write(*, *) 'combine_kappa_t: T_ai_vec,ee, ep = ', ee, ep
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
        write(*, *) 'combine_kappa_t: kappa on output, ee, ep = ', ee, ep
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
     integer, allocatable   :: temp_occ(:), temp_virt(:) 
     integer :: ind, nr_occ_used, nr_virt_used
!    ---------------------------------------------------------------------------

     orbital_rotation_indices(1:2, 1:length) = 0

!    fill orbital_rotation_indices vector:
     call alloc(temp_occ, nocct) 
     temp_occ = 0
     call get_occ_orbitals_indices(temp_occ)

     call virtual_orbital_indices_setup() 
     call alloc(temp_virt, length_virt_total)
     temp_virt = 0
     call get_virt_orbitals_indices(temp_virt, ee, ep)


!    combine them:
!    -------------
!    orbital_rotation_indices(1, i) = occ1,      ...      occ1, occ2,       ...      occ2, etc.
!                                     |- 'length_virt' times -| |- 'length_virt' times -|                    
!
!    orbital_rotation_indices(1, i) = virt1, virt2, ..., virtn,     virt1, virt2, ..., virtn, etc.
!                                     |-      occ1           -|     |-          occ2       -| 

     nr_occ_used = 0
     do i = 1, nocct
        if (temp_occ(i) .ne. 0) nr_occ_used = nr_occ_used + 1
     end do
     nr_virt_used = 0
     do i = 1, length_virt_total
        if (temp_virt(i) .ne. 0) nr_virt_used = nr_virt_used + 1
     end do

     do i = 1, nr_occ_used
        do j = 1, nr_virt_used
           if (temp_occ(i) .ne. 0 .and. temp_virt(j) .ne. 0) then
              ind = (i-1)*nr_virt_used + j
              orbital_rotation_indices(1, ind) = temp_occ(i)
              orbital_rotation_indices(2, ind) = temp_virt(j)
           end if
        end do
     end do

     call dealloc(temp_occ)
     call dealloc(temp_virt)

     debug_me = .false.
     if (debug_me) then
        do i = 1, length
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

    do iz = 1, nz

       do is = 1, length

          i = orbital_rotation_indices(1, is)
          s = orbital_rotation_indices(2, is)

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


  subroutine get_occ_orbitals_indices(indices_occ)
!    ---------------------------------------------------------------------------
     integer, intent(inout) :: indices_occ(*)
     integer :: iz, i, j, ind, isym
!    ---------------------------------------------------------------------------

     do isym = 1, nfsym
        if (xlr_indstr(1, isym) .ne. ' ') then
           call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                                 indices_occ,                 &
                                                 .true., .false., .false.,    &
                                                 isym)
        end if
     end do

  end subroutine

  subroutine get_virt_orbitals_indices(indices_virt, e, p)
!    ---------------------------------------------------------------------------
     integer, intent(inout) :: indices_virt(*)
     logical, intent(in) :: e, p
     integer :: iz, i, j, ind, isym
     integer :: ip(2), ap(2), ie(2), ae(2)
!    ---------------------------------------------------------------------------

     do isym = 1, nfsym
        if (xlr_indstr(3, isym) .ne. ' ') then
          call get_orbital_indices_from_strings(xlr_indstr(3, isym),     &
                                                indices_virt,            &
                                                .false., e, p, &
                                                isym)
        end if  
     end do

  end subroutine


  subroutine get_orbital_indices_from_strings(indices_string, indices_vect, &
                                              occ, e, p, isym)
!    ---------------------------------------------------------------------------
     character*72, intent(inout) :: indices_string(3, 2)
     integer, intent(inout) :: indices_vect(*)
     logical, intent(in)    :: occ, e, p
     integer, intent(in)    :: isym
     integer :: ifrp, ioff, i, ip, nvec, nvirt, ind1, ind2, ishift
!    ---------------------------------------------------------------------------


     if (occ) then

!       occuppied orbitals indices
!       --------------------------

        nvec = 1
        call numlst(indices_string, indices_vect,     &
                    nocc(isym), 1, nocc(isym), isym, nvec)

!       we need absolute indices:
        i = 0
        ioff = npsh(isym) + iorb(isym)
        do while (i < nocc(isym))
           i = i + 1
           if (indices_vect(i) .ne. 0) then
             indices_vect(i) = ioff + indices_vect(i)
           end if
        end do

        debug_me = .false.
        if (debug_me) then
           do i = 1, nocc(isym)
              write(*, *) 'OCC: isym, i, indices_vect(i) = ', &
                                isym, i, indices_vect(i)
           end do
        end if

     else

!       virtual orbitals indices
!       ------------------------

        call virtual_orbital_indices_setup()

        nvirt = norb(isym) - nocc(isym)

        call numlst(indices_string, indices_vect,     &
                    nvirt, -npsh(isym), nesh(isym), isym, nvec)

        ioff = 0
        ip = 0

!       indices_vect contains ALL virtual indices, now, we need to decide what we keep
!       either only positronic indices (p = true, e = false)
!           or only electronic indices (p = false, e = true)
!           or both (p = true, e = true)

        if (p .and. e) then
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
                 ioff = npsh(isym) + iorb(isym) + 1
                 indices_vect(i) = ioff + indices_vect(i)
              else if (indices_vect(i) .gt. nocc(isym)) then
!                check for electronic virtual indices
                 ioff = npsh(isym) + iorb(isym)
                 indices_vect(i) = ioff + indices_vect(i)
              end if
           end do

        else if (p .and. .not. e) then
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
                 ioff = npsh(isym) + iorb(isym) + 1
                 indices_vect(i) = ioff + indices_vect(i)
              else if (indices_vect(i) .gt. nocc(isym)) then
!                if there are any electronic indices, we don't want them here, so put to 0
                 indices_vect(i) = 0
              end if
           end do

        else if (e .and. .not. p) then
           ishift = 0
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
!                we don't want them here, so just check how many there are (kept in "ishift")
                 ishift = ishift + 1
                 indices_vect(i) = 0
              end if
           end do
           do i = ishift+1, nvirt
              if (indices_vect(i) .ne. 0) then
              !if (indices_vect(i) .gt. nocc(isym)) then
!                the electronic indices we want to write from indices_vect(1:), so:
                 ioff = npsh(isym) + iorb(isym)
                 indices_vect(i-ishift) = ioff + indices_vect(i)
!                zero-out remainings:
                 if (ishift .gt. 0) indices_vect(i) = 0
              end if
           end do
           if (ishift .gt. 0) then
              indices_vect((nvirt - ishift + 1):nvirt) = 0
           end if 

        end if

        debug_me = .false.
        if (debug_me) then
           write(*, *) 'VIRT: length_virt(isym), e, p = ', length_virt(isym), e, p
           do i = 1, nvirt
              write(*, *) 'VIRT: isym, i, indices_vect(i) = ', &
                                 isym, i, indices_vect(i)
           end do
        end if

     end if

  end subroutine


  subroutine virtual_orbital_indices_setup()
     integer :: isym

     virt_e = .true.
     virt_p = .true.
     length_virt_total = 0

     do isym = 1, nfsym
        if (xlr_skipep) then
           virt_p = .false.
           length_virt_p(isym) = 0
           length_virt_e(isym) = nesh(isym) - nocc(isym)
        else if (xlr_skipee) then
           virt_e = .false.
           length_virt_p(isym) = npsh(isym)
           length_virt_e(isym) = 0
        else
           length_virt_p(isym) = npsh(isym)
           length_virt_e(isym) = nesh(isym) - nocc(isym)
        end if
        length_virt(isym) = length_virt_e(isym) + length_virt_p(isym) 
        length_virt_total = length_virt_total + length_virt(isym)
     end do


  end subroutine


end module
