module dft_density

  implicit none

  public dft_density_init
  public gga_sym
  public get_n
  public get_s
  public get_gn
  public get_gn_nonhermitian
  public get_gs
  public get_gs_nonhermitian
  public lda_distribute_r
  public gga_distribute_r
  public lda_distribute_i
  public gga_distribute_i
  public lda_london_distribute_r

  private

  integer :: jxyz_start = 1 !noncollinear
! integer :: jxyz_start = 3 !collinear

  save

  integer :: nr_density_blocks(0:3, 0:7)
  integer :: nr_ao(2, 16, 0:3, 0:7)
  integer :: st_ao(2, 16, 0:3, 0:7)
  integer :: iq(0:3)

  logical :: dft_density_set    = .false.
  logical :: dft_density_synced = .false.

contains

  subroutine gga_sym(mat_dim, &
                     nz,      &
                     nr_mat,  &
                     mat,     &
                     isym,    &
                     ih)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(inout) :: mat(mat_dim, mat_dim, nz, nr_mat)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
!   ----------------------------------------------------------------------------
    integer                :: irep, iblock, imat, jxyz
    integer                :: i, j, ii, jj
    integer                :: st1, st2, nr1, nr2
    logical                :: skip(0:3)
    real(8)                :: average
!   ----------------------------------------------------------------------------

    do imat = 1, nr_mat

      skip = .false.
      if (ih(imat) ==  1) skip(1:3) = .true.
      if (ih(imat) == -1) skip(0:0) = .true.

      irep = isym(imat) - 1

      do jxyz = 0, 3
        if (.not. skip(jxyz)) then

          do iblock = 1, nr_density_blocks(jxyz, irep)

            nr1 = nr_ao(1, iblock, jxyz, irep)
            st1 = st_ao(1, iblock, jxyz, irep)
            nr2 = nr_ao(2, iblock, jxyz, irep)
            st2 = st_ao(2, iblock, jxyz, irep)

            if (.not. (st2 > st1)) then

              do i = 1, nr1
                do j = 1, nr2

                  ii = st1 + i - 1
                  jj = st2 + j - 1

                  average = mat(ii, jj, iq(jxyz), imat) &
                          + mat(jj, ii, iq(jxyz), imat)
                  average = 0.5d0*average

                  mat(ii, jj, iq(jxyz), imat) = average
                  mat(jj, ii, iq(jxyz), imat) = average

                end do
              end do
            end if
          end do
        end if
      end do
    end do

  end subroutine

  subroutine get_n(n, irep, mat_dim, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(*)
    real(8), intent(in)  :: ao(*)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    n = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   mat_dim,                  &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        n = n + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_s(s, irep, mat_dim, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(*)
    real(8), intent(in)  :: ao(*)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    s = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   mat_dim,                  &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        s(jxyz) = s(jxyz) + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_gn(gn,      &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao,      &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    gn = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   mat_dim,                  &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        do ixyz = 1, 3
          gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
        end do

      end do
    end do

    gn = 2.0d0*gn

  end subroutine

  subroutine get_gn_nonhermitian(gn,      &
                                 irep,    &
                                 mat_dim, &
                                 dmat,    &
                                 buffer,  &
                                 ao,      &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gn = 0.0d0

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do ixyz = 1, 3

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     mat_dim,                  &
                     ao(st2),                  &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     mat_dim,                  &
                     ao1(st2, ixyz),           &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gn(ixyz) = gn(ixyz) + ddot(nr1, ao(st1), 1, buffer, 1)

        end do
      end do
    end do

  end subroutine

  subroutine get_gs(gs,      &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao,      &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    gs = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   mat_dim,                  &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        do ixyz = 1, 3
          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
        end do

      end do
    end do

    gs = 2.0d0*gs

  end subroutine

  subroutine get_gs_nonhermitian(gs,      &
                                 irep,    &
                                 mat_dim, &
                                 dmat,    &
                                 buffer,  &
                                 ao,      &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gs = 0.0d0

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do ixyz = 1, 3

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     mat_dim,                  &
                     ao(st2),                  &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)

          call dgemv('N',                      &
                     nr1,                      &
                     nr2,                      &
                     1.0d0,                    &
                     dmat(st1, st2, iq(jxyz)), &
                     mat_dim,                  &
                     ao1(st2, ixyz),           &
                     1,                        &
                     0.0d0,                    &
                     buffer,                   &
                     1)

          gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                         + ddot(nr1, ao(st1), 1, buffer, 1)

        end do
      end do
    end do

  end subroutine

  subroutine lda_distribute_r(s, irep, mat_dim, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f = s*ao(ii)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)
        end do
      end do
    end do

  end subroutine

  subroutine gga_distribute_r(s, v, irep, mat_dim, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    real(8), intent(in)    :: v(3)
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = 0, 0
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f =       s   *ao (ii)    &
            + 2.0d0*v(1)*ao1(ii, 1) &
            + 2.0d0*v(2)*ao1(ii, 2) &
            + 2.0d0*v(3)*ao1(ii, 3)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

  subroutine lda_distribute_i(s, irep, mat_dim, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f = s(jxyz)*ao(ii)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

  subroutine gga_distribute_i(s, v, irep, mat_dim, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    real(8), intent(in)    :: v(3, 3)
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do jxyz = jxyz_start, 3
      do iblock = 1, nr_density_blocks(jxyz, irep)

        nr1 = nr_ao(1, iblock, jxyz, irep)
        st1 = st_ao(1, iblock, jxyz, irep)
        nr2 = nr_ao(2, iblock, jxyz, irep)
        st2 = st_ao(2, iblock, jxyz, irep)

        do i = 1, nr2

          ii = st2 + i - 1

          f =       s(   jxyz)*ao (ii)    &
            + 2.0d0*v(1, jxyz)*ao1(ii, 1) &
            + 2.0d0*v(2, jxyz)*ao1(ii, 2) &
            + 2.0d0*v(3, jxyz)*ao1(ii, 3)

          call daxpy(nr1,                     &
                     f,                       &
                     ao(st1),                 &
                     1,                       &
                     fmat(st1, ii, iq(jxyz)), &
                     1)

        end do
      end do
    end do

  end subroutine

  subroutine dft_density_init(is_parallel_calculation)

#ifdef VAR_MPI
    use dirac_mpi
#include "infpar.h"
#include "mpif.h"
#endif

!   ----------------------------------------------------------------------------
    logical, intent(in) :: is_parallel_calculation
!   ----------------------------------------------------------------------------
    integer             :: i, ic, irep, irep1, irep2, nr1, nr2, jxyz
    integer             :: inter_irep(0:3)
!   ----------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"

#ifdef VAR_MPI
!   only master sets values
    if (mytid == mparid) then
#endif

    if (.not. dft_density_set) then

      iq(0) = 1
      do jxyz = 1, 3
        iq(jxyz) = jm4pos(4 - jxyz)
        iq(jxyz) = iqtopq(iq(jxyz), isymax(jxyz, 2))
      end do

      inter_irep(0) = 0
      do jxyz = 1, 3
        inter_irep(jxyz) = isymax(jxyz, 2)
      end do

      nr_density_blocks = 0
      do irep = 0, nbsym - 1
        do jxyz = 0, 3
          i = 0
          do ic = 1, mc
            do irep1 = 0, nbsym - 1
              irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))

              nr1 = nbbas(irep1, ic)
              nr2 = nbbas(irep2, ic)

              if ((nr1 > 0) .and. (nr2 > 0)) then
                i = i + 1
                nr_ao(1, i, jxyz, irep) = nr1
                st_ao(1, i, jxyz, irep) = ibbas(irep1, ic) + 1
                nr_ao(2, i, jxyz, irep) = nr2
                st_ao(2, i, jxyz, irep) = ibbas(irep2, ic) + 1
              end if

            end do
          end do
          nr_density_blocks(jxyz, irep) = i
        end do
      end do

      dft_density_set = .true.
    end if

#ifdef VAR_MPI
    end if !if (mytid == mparid) then
#endif

#ifdef VAR_MPI
    if (is_parallel_calculation .and. .not. dft_density_synced) then

      call dirac_mpi_bcast(nr_density_blocks)
      call dirac_mpi_bcast(nr_ao)
      call dirac_mpi_bcast(st_ao)
      call dirac_mpi_bcast(iq)

      dft_density_synced = .true.
    end if
#endif

  end subroutine



  subroutine lda_london_distribute_r(s, irep, mat_dim, fmat, ao, ao_b1, cx, cy, cz)

    use dft_cb_interface

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, 12)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao_b1(mat_dim,*)
    real(8), intent(in)    :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz, jj, j
    integer                :: a, b, c
    real(8)                :: f
    real(8)                :: f_left
    real(8)                :: f_right
    real(8)                :: r_co(3)
    real(8)                :: aob1(mat_dim, 3)
    integer                :: itest, jtest
!   ----------------------------------------------------------------------------

   r_co(     1) = cx - dft_origin(1) 
   r_co(     2) = cy - dft_origin(2)
   r_co(     3) = cz - dft_origin(3)

!    do jxyz = 0, 0
!      do iblock = 1, nr_density_blocks(jxyz, irep)
!
!        nr1 = nr_ao(1, iblock, jxyz, irep)
!        st1 = st_ao(1, iblock, jxyz, irep)
!        nr2 = nr_ao(2, iblock, jxyz, irep)
!        st2 = st_ao(2, iblock, jxyz, irep)

        do a = 1, 3
          select case (a)
            case (1)
              b = 2; c = 3
            case (2)
              b = 3; c = 1
            case (3)
              b = 1; c = 2
          end select

          do j = 1, mat_dim
            jj = j
            do i = 1, mat_dim
              ii = i

              f = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

              fmat(ii, jj, (a-1)*4 + 1) = fmat(ii, jj, (a-1)*4 + 1) + s*f
            end do
          end do


        end do

!      end do
!    end do


  end subroutine

end module
