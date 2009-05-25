module dft_density

! use blub

  implicit none

  public init_density_blocks
  public get_n
  public get_gn
  public lda_distribute
  public gga_distribute

  private

  integer :: nr_density_blocks(0:3, 0:7)
  integer :: nr_ao(2, 16, 0:3, 0:7)
  integer :: st_ao(2, 16, 0:3, 0:7)
  integer :: iq(0:3)
  integer :: nr_ao_total

contains

  subroutine init_density_blocks()

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    integer :: ic, irep, irep1, irep2, nr1, nr2, ixyz, icount
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
    do ixyz = 1, 3
      iq(ixyz) = jm4pos(4 - ixyz)
      iq(ixyz) = iqtopq(iq(ixyz), isymax(ixyz, 2))
    end do

    inter_irep(0) = 0
    do ixyz = 1, 3
      inter_irep(ixyz) = isymax(ixyz, 2)
    end do

    nr_density_blocks = 0
    do irep = 0, nbsym - 1
      do ixyz = 0, 3
        icount = 0
        do ic = 1, mc
          do irep1 = 0, nbsym - 1
            irep2 = ieor(irep1, ieor(inter_irep(ixyz), irep))

            nr1 = nbbas(irep1, ic)
            nr2 = nbbas(irep2, ic)

            if ((nr1 > 0) .and. (nr2 > 0)) then
              icount = icount + 1
              nr_ao(1, icount, ixyz, irep) = nr1
              st_ao(1, icount, ixyz, irep) = ibbas(irep1, ic) + 1
              nr_ao(2, icount, ixyz, irep) = nr2
              st_ao(2, icount, ixyz, irep) = ibbas(irep2, ic) + 1
            end if

          end do
        end do
        nr_density_blocks(ixyz, irep) = icount
      end do
    end do

  end subroutine

  subroutine get_n(n, ixyz, irep, dmat, buffer, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    integer, intent(in)  :: ixyz
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    n = 0.0d0

    do iblock = 1, nr_density_blocks(ixyz, irep)

      nr1 = nr_ao(1, iblock, ixyz, irep)
      st1 = st_ao(1, iblock, ixyz, irep)
      nr2 = nr_ao(2, iblock, ixyz, irep)
      st2 = st_ao(2, iblock, ixyz, irep)

      call dgemv('N',            &
                 nr1,            &
                 nr2,            &
                 1.0d0,          &
                 dmat(st1, st2), &
                 nr_ao_total,    &
                 ao(st2),        &
                 1,              &
                 0.0d0,          &
                 buffer,         &
                 1)

      n = n + ddot(nr1, ao(st1), 1, buffer, 1)

    end do

  end subroutine

  subroutine get_gn(gn,     &
                    ixyz,   &
                    irep,   &
                    dmat,   &
                    buffer, &
                    ao,     &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: ixyz
    integer, intent(in)  :: irep
    real(8), intent(in)  :: dmat(nr_ao_total, nr_ao_total)
    real(8), intent(in)  :: buffer(nr_ao_total)
    real(8), intent(in)  :: ao(nr_ao_total)
    real(8), intent(in)  :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, i
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    gn = 0.0d0

    do iblock = 1, nr_density_blocks(ixyz, irep)

      nr1 = nr_ao(1, iblock, ixyz, irep)
      st1 = st_ao(1, iblock, ixyz, irep)
      nr2 = nr_ao(2, iblock, ixyz, irep)
      st2 = st_ao(2, iblock, ixyz, irep)

      call dgemv('N',            &
                 nr1,            &
                 nr2,            &
                 1.0d0,          &
                 dmat(st1, st2), &
                 nr_ao_total,    &
                 ao(st2),        &
                 1,              &
                 0.0d0,          &
                 buffer,         &
                 1)

      do i = 1, 3
        gn(i) = gn(i) &
              + ddot(nr1, ao1(st1, i), 1, buffer, 1)
      end do

    end do

    gn = 2.0d0*gn

  end subroutine

  subroutine lda_distribute(s, ixyz, irep, fmat, ao)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: ixyz
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)    :: ao(nr_ao_total)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i
    real(8)                :: f
!   ----------------------------------------------------------------------------

    do iblock = 1, nr_density_blocks(ixyz, irep)

      nr1 = nr_ao(1, iblock, ixyz, irep)
      st1 = st_ao(1, iblock, ixyz, irep)
      nr2 = nr_ao(2, iblock, ixyz, irep)
      st2 = st_ao(2, iblock, ixyz, irep)

      do i = 1, nr2

        f = s*ao(st2 + i - 1)

        call daxpy(nr1,                                   &
                   f,                                     &
                   ao(st1),                               &
                   1,                                     &
                   fmat((st2 + i - 2)*nr_ao_total + st1), &
                   1)

      end do
    end do

  end subroutine

  subroutine gga_distribute(s, v, ixyz, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    real(8)                :: v(3)
    integer, intent(in)    :: ixyz
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(*)
    real(8), intent(in)    :: ao(nr_ao_total)
    real(8), intent(in)    :: ao1(nr_ao_total, 3)
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i
    real(8)                :: f
!   ----------------------------------------------------------------------------

    v = 2.0d0*v

    do iblock = 1, nr_density_blocks(ixyz, irep)

      nr1 = nr_ao(1, iblock, ixyz, irep)
      st1 = st_ao(1, iblock, ixyz, irep)
      nr2 = nr_ao(2, iblock, ixyz, irep)
      st2 = st_ao(2, iblock, ixyz, irep)

      do i = 1, nr2

        f = s   *ao (st2 + i - 1)    &
          + v(1)*ao1(st2 + i - 1, 1) &
          + v(2)*ao1(st2 + i - 1, 2) &
          + v(3)*ao1(st2 + i - 1, 3)

        call daxpy(nr1,                                   &
                   f,                                     &
                   ao(st1),                               &
                   1,                                     &
                   fmat((st2 + i - 2)*nr_ao_total + st1), &
                   1)

      end do
    end do

  end subroutine

  subroutine h_anti_wrap_distribute(use_gga, s, v, irep, fmat, ao, ao1)

!   ----------------------------------------------------------------------------
    logical, intent(in)    :: use_gga
    real(8), intent(in)    :: s(3)
    real(8), intent(in)    :: v(3, 3)
    integer, intent(in)    :: irep
    real(8), intent(inout) :: fmat(nr_ao_total, nr_ao_total, *)
    real(8), intent(in)    :: ao(*)
    real(8), intent(in)    :: ao1(*)
!   ----------------------------------------------------------------------------
    integer                :: ixyz
!   ----------------------------------------------------------------------------

    do ixyz = 1, 3
      if (use_gga) then
        call gga_distribute(s(ixyz), v(1, ixyz), ixyz, irep, fmat(1, 1, iq(ixyz)), ao, ao1)
      else
        call lda_distribute(s(ixyz),             ixyz, irep, fmat(1, 1, iq(ixyz)), ao)
      end if
    end do

  end subroutine

end module
