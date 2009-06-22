module dft_cb_interface

! radovan: this module is the interface between the dft code and common blocks
!          ideally no common blocks should be included in the dft code
!          outside this module to keep the interface clean
!          hide as much as you can for public - it will save us lots of trouble in future

  implicit none

  public init_dft_cb
#ifdef VAR_MPI
  public sync_dft_cb
#endif

  save

  integer, parameter, public :: max_nr_shells = 750

  integer, public :: nr_qu
  integer, public :: nr_ao_total
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
  subroutine my_pack_integer(pack_integer, ipos, ndim, array, i_am_the_master)

!   ----------------------------------------------------------------------------
    integer, intent(inout) :: pack_integer(:)
    integer, intent(inout) :: ipos
    integer, intent(in)    :: ndim
    integer, intent(inout) :: array(ndim)
    logical, intent(in)    :: i_am_the_master
!   ----------------------------------------------------------------------------
    integer                :: i
!   ----------------------------------------------------------------------------

    do i = 1, ndim
      if (i_am_the_master) then
        pack_integer(ipos + i - 1) = array(i)
      else
        array(i) = pack_integer(ipos + i - 1)
      end if
    end do
    ipos = ipos + ndim

  end subroutine

  subroutine my_pack_double_precision(pack_double_precision, ipos, ndim, array, i_am_the_master)

!   ----------------------------------------------------------------------------
    real(8), intent(inout) :: pack_double_precision(:)
    integer, intent(inout) :: ipos
    integer, intent(in)    :: ndim
    real(8), intent(inout) :: array(ndim)
    logical, intent(in)    :: i_am_the_master
!   ----------------------------------------------------------------------------
    integer                :: i
!   ----------------------------------------------------------------------------

    do i = 1, ndim
      if (i_am_the_master) then
        pack_double_precision(ipos + i - 1) = array(i)
      else
        array(i) = pack_double_precision(ipos + i - 1)
      end if
    end do
    ipos = ipos + ndim

  end subroutine

  subroutine sync_dft_cb()

#include "mpif.h"
#include "infpar.h"

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    logical             :: i_am_the_master
    integer             :: idx, ierr, ipos
    logical             :: pack_logical(2)
    real(8)             :: pack_double_precision(max_nr_shells*3*2)
    integer             :: pack_integer(2156)
!                                   1000000 is a phantasy number
    integer, parameter  :: nbytes = 1000000
    character           :: scratch(nbytes)
!   ----------------------------------------------------------------------------

    if (mytid == mparid) then
      i_am_the_master = .true.
    else
      i_am_the_master = .false.
    end if

    if (i_am_the_master) then

      pack_integer(1) = nr_qu
      pack_integer(2) = nr_ao_total
      pack_integer(3) = nr_boson_sym
      pack_integer(4) = nr_ls_components
      pack_integer(5) = nr_open_shells
      pack_integer(6) = nr_electrons_orbital
      pack_integer(7) = nr_electrons_orbital_active
      pack_integer(8) = print_unit
      ipos = 9
      call my_pack_integer(pack_integer, ipos, size(nr_density_blocks), nr_density_blocks, i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(nr_ao            ), nr_ao            , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(st_ao            ), st_ao            , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(iq               ), iq               , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(dft_cb_pq_to_uq  ), dft_cb_pq_to_uq  , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(dft_cb_uq_to_pq  ), dft_cb_uq_to_pq  , i_am_the_master)

      pack_logical(1) = is_gga_functional
      pack_logical(2) = is_parallel_calculation

      ipos = 1
      call my_pack_double_precision(pack_double_precision, ipos, size(dft_cb_cent), dft_cb_cent, i_am_the_master)

      idx = 1
      call mpi_pack(pack_logical,                  &
                    size(pack_logical),            &
                    mpi_logical,                   &
                    scratch,                       &
                    nbytes,                        &
                    idx,                           &
                    mpi_comm_world,                &
                    ierr)
      call mpi_pack(pack_double_precision,         &
                    size(pack_double_precision),   &
                    mpi_double_precision,          &
                    scratch,                       &
                    nbytes,                        &
                    idx,                           &
                    mpi_comm_world,                &
                    ierr)
      call mpi_pack(pack_integer,                  &
                    size(pack_integer),            &
                    my_MPI_INTEGER,                &
                    scratch,                       &
                    nbytes,                        &
                    idx,                           &
                    mpi_comm_world,                &
                    ierr)

    end if

    call mpi_bcast(scratch,                        &
                   nbytes,                         &
                   mpi_packed,                     &
                   0,                              &
                   mpi_comm_world,                 &
                   ierr)

    if (.not. i_am_the_master) then

      idx = 1
      call mpi_unpack(scratch,                     &
                      nbytes,                      &
                      idx,                         &
                      pack_logical,                &
                      size(pack_logical),          &
                      mpi_logical,                 &
                      mpi_comm_world,              &
                      ierr)
      call mpi_unpack(scratch,                     &
                      nbytes,                      &
                      idx,                         &
                      pack_double_precision,       &
                      size(pack_double_precision), &
                      mpi_double_precision,        &
                      mpi_comm_world,              &
                      ierr)
      call mpi_unpack(scratch,                     &
                      nbytes,                      &
                      idx,                         &
                      pack_integer,                &
                      size(pack_integer),          &
                      my_MPI_INTEGER,              &
                      mpi_comm_world,              &
                      ierr)

      nr_qu                       = pack_integer(1)
      nr_ao_total                 = pack_integer(2)
      nr_boson_sym                = pack_integer(3)
      nr_ls_components            = pack_integer(4)
      nr_open_shells              = pack_integer(5)
      nr_electrons_orbital        = pack_integer(6)
      nr_electrons_orbital_active = pack_integer(7)
      print_unit                  = pack_integer(8)
      ipos = 9
      call my_pack_integer(pack_integer, ipos, size(nr_density_blocks), nr_density_blocks, i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(nr_ao            ), nr_ao            , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(st_ao            ), st_ao            , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(iq               ), iq               , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(dft_cb_pq_to_uq  ), dft_cb_pq_to_uq  , i_am_the_master)
      call my_pack_integer(pack_integer, ipos, size(dft_cb_uq_to_pq  ), dft_cb_uq_to_pq  , i_am_the_master)

      is_gga_functional       = pack_logical(1)
      is_parallel_calculation = pack_logical(2)

      ipos = 1
      call my_pack_double_precision(pack_double_precision, ipos, size(dft_cb_cent), dft_cb_cent, i_am_the_master)

    end if

  end subroutine
#endif /* ifdef VAR_MPI */

end module
