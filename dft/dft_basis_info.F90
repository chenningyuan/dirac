module dft_basis_info

  use memory_allocator
#ifdef VAR_MPI
  use dft_mpi
#endif

  implicit none

  public init_dft_basis_info

  save

  integer, public              :: nr_cartesian_ao     = 0
  integer, public              :: dft_nr_shells       = 0
  integer, public              :: dft_nr_large_shells = 0
  integer, public              :: dft_maxrep          = 0
  integer, public              :: dft_maxopr          = 0

  real(8), public, allocatable :: primitive_exp1(:)
  real(8), public, allocatable :: primitive_exp2(:, :)
  real(8), public, allocatable :: contraction_coef1(:)
  real(8), public, allocatable :: contraction_coef2(:, :)
  integer, public, allocatable :: nr_nonzero_contraction_coefs(:)
  integer, public, allocatable :: cartesian_deg(:)
  integer, public, allocatable :: cartesian_off(:)
  integer, public, allocatable :: spherical_deg(:)
  integer, public, allocatable :: spherical_off(:)
  integer, public, allocatable :: n_quantum_nr(:)
  logical, public, allocatable :: is_spherical(:)
  integer, public, allocatable :: dft_istbao(:)
  integer, public, allocatable :: dft_isymao(:, :)
  integer, public, allocatable :: dft_ipind(:)
  real(8), public, allocatable :: dft_pt(:)
  integer, public, allocatable :: dft_kstrt(:)
  real(8), public, allocatable :: dft_cent(:, :, :)
  integer, public              :: dft_isymax(3, 2)

  real(8), public, allocatable :: dft_csp(:)
  integer, public, allocatable :: dft_ispadr(:)

  integer, public, allocatable :: dft_icos(:, :)
  integer, public, allocatable :: dft_ibbas(:, :)
  integer, public, allocatable :: dft_nbbas(:, :)

  integer, public, allocatable :: l_exp(:)
  integer, public, allocatable :: m_exp(:)
  integer, public, allocatable :: n_exp(:)
  real(8), public, allocatable :: sbuffer(:)

  integer, parameter :: max_order_limit = 2

  integer, public    :: ndimfactor      = 0
  integer, public    :: max_ao_g_order  = 0
  integer, public    :: max_ao_m_order  = 0
  logical, public    :: need_ao_order(0:max_order_limit, 0:max_order_limit) = .false.
  integer, public    :: ao_start     (0:max_order_limit, 0:max_order_limit) = 0
  integer, public    :: ao_end       (0:max_order_limit, 0:max_order_limit) = 0

  private

  logical :: dft_basis_info_allocated = .false.
  logical :: dft_basis_info_set       = .false.
#ifdef VAR_MPI
  logical :: dft_basis_info_synced    = .false.
#endif

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "shells.h"
#include "symmet.h"
#include "aovec.h"
#include "primit.h"
#include "sphtrm.h"
#include "dcbbas.h"
#include "pincom.h"

#ifdef VAR_MPI
#include "infpar.h"
#include "mpif.h"
#endif

contains

  subroutine init_dft_basis_info(is_parallel_calculation)

!   ----------------------------------------------------------------------------
    logical, intent(in) :: is_parallel_calculation
!   ----------------------------------------------------------------------------
    integer             :: iprim, i, j, icount, i_c, i_s, contracted_shell_nr, ishell, isym
#ifdef VAR_MPI
    logical             :: i_am_the_master
#endif
!   ----------------------------------------------------------------------------

    if (.not. dft_basis_info_allocated) then

      call alloc(primitive_exp1, mxshel)
      call alloc(primitive_exp2, mxshel, mxcont)
      call alloc(contraction_coef1, mxshel)
      call alloc(contraction_coef2, mxshel, mxcont)
      call alloc(nr_nonzero_contraction_coefs, mxshel)
      call alloc(cartesian_deg, mxshel)
      call alloc(cartesian_off, mxshel)
      call alloc(spherical_deg, mxshel)
      call alloc(spherical_off, mxshel)
      call alloc(n_quantum_nr, mxshel)
      allocate(is_spherical(mxshel))
      call alloc(dft_istbao, mxshel)
      call alloc(dft_isymao, mxqn, mxaqn)
      call alloc(dft_ipind, mxcorb)
      allocate(dft_pt(0:7))
      call alloc(dft_kstrt, mxshel)
      call alloc(dft_cent, mxshel, 3, 2)

      nr_cartesian_ao = 0
      do ishell = 1, mxshel
        nr_cartesian_ao = nr_cartesian_ao + kckt(ishell)
      end do
      dft_nr_shells       = kmax
      dft_nr_large_shells = nlrgsh

      call alloc(dft_csp,    ncsp)
      call alloc(dft_ispadr, mxqn)

      call alloc(dft_icos, 8, 2)
      allocate(dft_ibbas(0:7, 0:2))
      allocate(dft_nbbas(0:7, 0:2))

#ifdef VAR_MPI
      if (is_parallel_calculation) then
        call dft_mpi_bcast(nr_cartesian_ao)
      end if
#endif

      call alloc(l_exp,   nr_cartesian_ao)
      call alloc(m_exp,   nr_cartesian_ao)
      call alloc(n_exp,   nr_cartesian_ao)
      call alloc(sbuffer, nr_cartesian_ao)

      dft_basis_info_allocated = .true.

    end if

#ifdef VAR_MPI
    i_am_the_master = .false.
    if (mytid == mparid) i_am_the_master = .true.

!   only master sets values
    if (i_am_the_master) then
#endif

    if (.not. dft_basis_info_set) then

      nr_nonzero_contraction_coefs = 0
      cartesian_deg                = 0
      cartesian_off                = 0
      spherical_deg                = 0
      spherical_off                = 0
      n_quantum_nr                 = 0
      is_spherical                 = .false.

      dft_ipind  = ipind
      dft_istbao = istbao
      dft_isymao = isymao
      dft_pt     = pt
      dft_kstrt  = kstrt
      dft_cent   = cent
      dft_isymax = isymax

      dft_csp    = csp
      dft_ispadr = ispadr

      dft_icos  = icos
      dft_ibbas = ibbas
      dft_nbbas = nbbas

      dft_maxrep = maxrep
      dft_maxopr = maxopr

      i_c = 1
      i_s = 1

      do ishell = 1, kmax

        n_quantum_nr(ishell)  = nhkt(ishell)

        cartesian_deg(ishell) = kckt(ishell)
        spherical_deg(ishell) = khkt(ishell)

        cartesian_off(ishell) = i_c - 1
        spherical_off(ishell) = i_s - 1

        call lmnval(n_quantum_nr(ishell), cartesian_deg(ishell), l_exp(i_c), m_exp(i_c), n_exp(i_c))

        i_c = i_c + cartesian_deg(ishell)
        i_s = i_s + spherical_deg(ishell)

        is_spherical(ishell)  = sphr(ishell)

        contracted_shell_nr = numcf(ishell)

        i = jstrt(ishell) + 1
        j = jstrt(ishell) + nuco(ishell)

        contraction_coef1(ishell) = priccf(i, contracted_shell_nr)
        primitive_exp1(ishell)    = priexp(i)

        icount = 0
        do iprim = i + 1, j
          if (dabs(priccf(iprim, contracted_shell_nr)) > 1.0d-30) then
            icount = icount + 1
            nr_nonzero_contraction_coefs(ishell) = icount
            contraction_coef2(ishell, icount)    = priccf(iprim, contracted_shell_nr)
            primitive_exp2(ishell, icount)       = priexp(iprim)
          end if
        end do

      end do

      dft_basis_info_set = .true.

    end if

#ifdef VAR_MPI
    end if !if (i_am_the_master) then
#endif

#ifdef VAR_MPI
    if (is_parallel_calculation .and. .not. dft_basis_info_synced) then

!     is done above
!     call dft_mpi_bcast(nr_cartesian_ao             )
      call dft_mpi_bcast(dft_nr_shells               )
      call dft_mpi_bcast(dft_nr_large_shells         )
      call dft_mpi_bcast(dft_maxrep                  )
      call dft_mpi_bcast(dft_maxopr                  )

      call dft_mpi_bcast(primitive_exp1              )
      call dft_mpi_bcast(primitive_exp2              )
      call dft_mpi_bcast(contraction_coef1           )
      call dft_mpi_bcast(contraction_coef2           )
      call dft_mpi_bcast(nr_nonzero_contraction_coefs)
      call dft_mpi_bcast(cartesian_deg               )
      call dft_mpi_bcast(cartesian_off               )
      call dft_mpi_bcast(spherical_deg               )
      call dft_mpi_bcast(spherical_off               )
      call dft_mpi_bcast(n_quantum_nr                )
      call dft_mpi_bcast(is_spherical                )
      call dft_mpi_bcast(dft_istbao                  )
      call dft_mpi_bcast(dft_isymao                  )
      call dft_mpi_bcast(dft_ipind                   )
      call dft_mpi_bcast(dft_pt                      )
      call dft_mpi_bcast(dft_kstrt                   )
      call dft_mpi_bcast(dft_cent                    )
      call dft_mpi_bcast(dft_isymax                  )

      call dft_mpi_bcast(dft_csp                     )
      call dft_mpi_bcast(dft_ispadr                  )

      call dft_mpi_bcast(dft_icos                    )
      call dft_mpi_bcast(dft_ibbas                   )
      call dft_mpi_bcast(dft_nbbas                   )

      call dft_mpi_bcast(l_exp                       )
      call dft_mpi_bcast(m_exp                       )
      call dft_mpi_bcast(n_exp                       )

      call dft_mpi_bcast(ndimfactor    )
      call dft_mpi_bcast(max_ao_g_order)
      call dft_mpi_bcast(max_ao_m_order)
      call dft_mpi_bcast(need_ao_order )
      call dft_mpi_bcast(ao_start      )
      call dft_mpi_bcast(ao_end        )

      dft_basis_info_synced = .true.

    end if
#endif

  end subroutine

end module
