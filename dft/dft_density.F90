module dft_density

! use blub

  implicit none

  public init_density_blocks
  public get_n, get_gn, get_gn_nonhermitian
  public get_s, get_gs, get_gs_nonhermitian
  public lda_distribute
  public gga_distribute
  public lda_distribute_qi
  public gga_distribute_qi

  private

  integer :: nr_density_blocks(0:3, 0:7)
  integer :: nr_ao(2, 16, 0:3, 0:7)
  integer :: st_ao(2, 16, 0:3, 0:7)
  integer :: iq(0:3)
  integer :: nr_ao_total

! noncollinear
  integer :: jxyz_start = 1
! collinear approximation
! integer :: jxyz_start = 3

contains

  subroutine init_density_blocks()

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    integer :: ic, irep, irep1, irep2, nr1, nr2, jxyz, icount
    integer :: inter_irep(0:3)
!   ----------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"

    nr_ao_total = ntbas(0)

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
        icount = 0
        do ic = 1, mc
          do irep1 = 0, nbsym - 1
            irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))

            nr1 = nbbas(irep1, ic)
            nr2 = nbbas(irep2, ic)

            if ((nr1 > 0) .and. (nr2 > 0)) then
              icount = icount + 1
              nr_ao(1, icount, jxyz, irep) = nr1
              st_ao(1, icount, jxyz, irep) = ibbas(irep1, ic) + 1
              nr_ao(2, icount, jxyz, irep) = nr2
              st_ao(2, icount, jxyz, irep) = ibbas(irep2, ic) + 1
            end if

          end do
        end do
        nr_density_blocks(jxyz, irep) = icount
      end do
    end do

  end subroutine

  subroutine get_n(n, irep, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8)              :: ddot
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
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        n = n + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_s(s, irep, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8)              :: ddot
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
                   nr_ao_total,              &
                   ao(st2),                  &
                   1,                        &
                   0.0d0,                    &
                   buffer,                   &
                   1)

        s(jxyz) = s(jxyz) + ddot(nr1, ao(st1), 1, buffer, 1)

      end do
    end do

  end subroutine

  subroutine get_gn(gn,     &
                    irep,   &
                    dmat,   &
                    buffer, &
                    ao,     &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
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

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
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

  subroutine get_gn_nonhermitian(gn,     &
                                 irep,   &
                                 dmat,   &
                                 buffer, &
                                 ao,     &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
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
                     nr_ao_total,              &
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
                     nr_ao_total,              &
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

  subroutine get_gs(gs,     &
                    irep,   &
                    dmat,   &
                    buffer, &
                    ao,     &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
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

        call dgemv('N',                      &
                   nr1,                      &
                   nr2,                      &
                   1.0d0,                    &
                   dmat(st1, st2, iq(jxyz)), &
                   nr_ao_total,              &
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

  subroutine get_gs_nonhermitian(gs,     &
                                 irep,   &
                                 dmat,   &
                                 buffer, &
                                 ao,     &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
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
                     nr_ao_total,              &
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
                     nr_ao_total,              &
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

  subroutine lda_distribute(s, irep, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
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

  subroutine gga_distribute(s, v, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    real(8), intent(in)    :: v(3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao1(nr_ao_total, 3)
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

  subroutine lda_distribute_qi(s, irep, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
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

  subroutine gga_distribute_qi(s, v, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    real(8), intent(in)    :: v(3, 3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao1(nr_ao_total, 3)
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

end module
