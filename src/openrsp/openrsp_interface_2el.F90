module openrsp_interface_2el

   use dirac_interface

   implicit none

   public grcont

   private

#include "mxcent.h"
#include "cbihr2.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbgrd.h"
#include "dorps.h"
#include "exeinf.h"

contains

  subroutine grcont(work,     &
                    lwork,    &
                    g,        &
                    g_dim,    &
                    geo,      &
                    london,   &
                    max_derv, &
                    iatom_in, &
                    energy,   &
                    fock,     &
                    dmat,     &
                    nr_dmat)

!   ----------------------------------------------------------------------------
    real(8)              :: work(*)
    integer, intent(in)  :: lwork

    real(8), intent(out) :: g(*)
    integer, intent(in)  :: g_dim

    logical, intent(in)  :: geo
    logical, intent(in)  :: london

    integer, intent(in)  :: max_derv
    integer, intent(in)  :: iatom_in

    logical, intent(in)  :: energy
    logical, intent(in)  :: fock

    real(8), intent(in)  :: dmat(*)
    integer, intent(in)  :: nr_dmat
!   ----------------------------------------------------------------------------
    real(8), allocatable :: tmat(:), fmat(:)

    logical              :: gencnt

    integer              :: i2typ, itype, ifctyp(8), iatom, fmat_dim

    real(8)              :: skip(3, 14, 3)
    real(8)              :: screen, f
!   ----------------------------------------------------------------------------

    dopert = .true.

    ftronv = .true.
    relgrd = .true.
    call rsetsym(.false., relgrd, .false., .false.)
    call rtroini(.false.)

    call troinv(work, lwork)

    screen = -1.0d0

!   set ifctyp
    ifctyp    = 0
    ifctyp(1) = 13 !coulomb and exchange
    ifctyp(2) = 22 !only exchange
    ifctyp(3) = 22
    ifctyp(4) = 22
    ifctyp(5) = 13 !coulomb and exchange
    ifctyp(6) = 22 !only exchange
    ifctyp(7) = 22
    ifctyp(8) = 22

    iatom = 0

    if (energy) then

      gencnt   = .true.
      fmat_dim = g_dim

      if (geo) then
        itype = 2
      else if (london) then
        itype = -5
      else
        call quit('grcont: differentiated integrals not recognized')
      end if

    else if (fock) then

      gencnt   = .false.
      fmat_dim = 3*n2bbasxq*max_derv

      if (geo) then
        itype = 8
        iatom = iatom_in
      else if (london) then
        itype = -5
      else
        call quit('grcont: differentiated integrals not recognized')
      end if

    else
      call quit('programming error in call to grcont')
    end if

    f = 1.0d0
    if (fock) then
      if (geo)                        f = 0.25d0
      if (london .and. max_derv == 1) f = 2.00d0
      if (london .and. max_derv == 2) f = 4.00d0
    end if

    allocate(tmat(fmat_dim))
    allocate(fmat(fmat_dim))

    fmat = 0.0d0

    if (include_ll_integrals) then

      i2typ    = 1
      intclass = 0

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, tmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    if (include_ls_integrals) then

      i2typ    = 2
      intclass = 1

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, tmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    if (include_ss_integrals) then

      i2typ    = 3
      intclass = 2

      tmat = 0.0d0
      call dzero(skip, 3*14*3)

      call twoint(work, lwork, tmat, dmat, 4*nr_dmat, (/0/),             &
                  ifctyp, (/0.0d0/), (/0/), (/0/), itype,                &
                  max_derv, iatom, .true., .true., .false., tktime,      &
                  1, 0, 0, 0, 0, .false., (/0/), i2typ,                  &
                  icedif, screen, (/0.0d0/), (/0.0d0/), (/0.0d0/), skip, &
                  relcal, gencnt, (/0.0d0/), (/0.0d0/))

      call daxpy(fmat_dim, f, tmat, 1, fmat, 1)

    end if

    call dcopy(g_dim, fmat, 1, g, 1)

    deallocate(fmat)

  end subroutine

end module
