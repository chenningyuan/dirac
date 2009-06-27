module dft_cb_interface

! radovan: this module is the interface between the dft code and common blocks
!          ideally no common blocks should be included in the dft code
!          outside this module to keep the interface clean
!          hide as much as you can for public - it will save us lots of trouble in future

#ifdef VAR_MPI
  use dft_mpi
#endif

  implicit none

  public init_dft_cb
#ifdef VAR_MPI
  public sync_dft_cb
#endif

  save

  integer, parameter, public :: max_nr_shells = 750

  integer, public :: nr_qu
  integer, public :: nr_ao_total
  integer, public :: nr_ao_large
  integer, public :: nr_ao_small
  integer, public :: nr_boson_sym
  integer, public :: nr_ls_components
  integer, public :: nr_open_shells
  integer, public :: nr_electrons_orbital
  integer, public :: nr_electrons_orbital_active
  integer, public :: print_unit
  integer, public :: nr_density_blocks(0:3, 0:7)
  integer, public :: nr_ao(2, 16, 0:3, 0:7)
  integer, public :: st_ao(2, 16, 0:3, 0:7)
  integer, public :: iq(0:3)
  integer, public :: dft_cb_pq_to_uq(4, 0:7)
  integer, public :: dft_cb_uq_to_pq(4, 0:7)

  logical, public :: is_gga_functional
  logical, public :: is_parallel_calculation

  real(8), public :: dft_cb_cent(max_nr_shells, 3, 2)

  private

contains

  subroutine init_dft_cb()

#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "functionals.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "shells.h"
#include "symmet.h"
#include "nuclei.h"

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    integer :: i, j
    integer :: ic, irep, irep1, irep2, nr1, nr2, jxyz, icount
    integer :: inter_irep(0:3)
!   ----------------------------------------------------------------------------

    nr_qu                       = nz
    nr_ao_total                 = ntbas(0)
    nr_ao_large                 = ntbas(1)
    nr_ao_small                 = ntbas(2)
    nr_boson_sym                = nbsym
    nr_ls_components            = mc
    nr_open_shells              = nasht
    nr_electrons_orbital        = nelect
    nr_electrons_orbital_active = naelec
    print_unit                  = lupri

    is_gga_functional       = dftisgga()
    is_parallel_calculation = parcal

    dft_cb_cent(1:max_nr_shells, 1:3, 1:2) = cent(1:max_nr_shells, 1:3, 1:2)

    dft_cb_pq_to_uq = 0
    dft_cb_uq_to_pq = 0
    do i = 1, 4
      do j = 0, nbsym - 1
        dft_cb_pq_to_uq(i, j) = ipqtoq(i, j)
        dft_cb_uq_to_pq(i, j) = iqtopq(i, j)
      end do
    end do

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

#ifdef VAR_MPI
  subroutine sync_dft_cb()

    call dft_mpi_bcast(nr_qu                      )
    call dft_mpi_bcast(nr_ao_total                )
    call dft_mpi_bcast(nr_ao_large                )
    call dft_mpi_bcast(nr_ao_small                )
    call dft_mpi_bcast(nr_boson_sym               )
    call dft_mpi_bcast(nr_ls_components           )
    call dft_mpi_bcast(nr_open_shells             )
    call dft_mpi_bcast(nr_electrons_orbital       )
    call dft_mpi_bcast(nr_electrons_orbital_active)
    call dft_mpi_bcast(print_unit                 )
    call dft_mpi_bcast(nr_density_blocks          )
    call dft_mpi_bcast(nr_ao                      )
    call dft_mpi_bcast(st_ao                      )
    call dft_mpi_bcast(iq                         )
    call dft_mpi_bcast(dft_cb_pq_to_uq            )
    call dft_mpi_bcast(dft_cb_uq_to_pq            )
    call dft_mpi_bcast(is_gga_functional          )
    call dft_mpi_bcast(is_parallel_calculation    )
    call dft_mpi_bcast(dft_cb_cent                )

  end subroutine
#endif /* ifdef VAR_MPI */

end module
