module softness_m

! softness = \sum_{ia} [ (n_i * n_a) / (e_a - e_i) ]

! idea:                ulf
! first coded version: radovan

  use kinds_m
  use memory_allocator

  implicit none

contains

  subroutine get_softness(softness, orbital_density, gao_dim, gao, mo_coef, eigenvalue)

    real(kind_dp), intent(out) :: softness
    real(kind_dp), intent(out) :: orbital_density(:)
                                ! it is intent(out) and not just an internal
                                ! variable because i used it to debug it
                                ! by checking integrated orbital densities

    integer,       intent(in)  :: gao_dim
    real(kind_dp), intent(in)  :: gao(gao_dim)
    real(kind_dp), intent(in)  :: mo_coef(:)
    real(kind_dp), intent(in)  :: eigenvalue(:)
                                ! both orbital_density(:) and eigenvalue(:)
                                ! have the dimension norbt
!   ----------------------------------------------------------------------------
    integer                    :: a, i, k, kk, l, ll, j, p, pp, po, iz, izo, ifsym, istart, iend
    real(kind_dp)              :: factor
    logical                    :: is_there
    integer                    :: nr_i
    integer                    :: nr_a
    integer, allocatable       :: p_i(:)
    integer, allocatable       :: p_a(:)

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"

    nr_i = nish(1) + nish(2)
    nr_a = norbt - nr_i

    call alloc(p_i, nr_i)
    call alloc(p_a, nr_a)

!   the pointers p_i and p_a point to the orbital numbers
!   this can be done once for all times outside this routine
!   it does not have to be repeated like here in each point

    k = 0
    do i = 1, nfsym
      if (nish(i) > 0) then

        istart = norb(1)*(i-1) + npsh(i) + 1
        iend   = istart + nish(i) - 1

        do j = istart, iend

          k = k + 1
          p_i(k) = j

        end do
      end if
    end do

    k = 1
    do j = 1, norbt
      is_there = .false.
      do i = 1, nr_i
        if (j == p_i(i)) then
          is_there = .true.
        end if
      end do

      if (is_there) then
      else
        p_a(k) = j
        k = k + 1
      end if

    end do

    orbital_density = 0.0d0

!   i am sure there are better ways to calculate it than below
!   using optimized math routines
!   this is just a first working version

    do ifsym = 1, nfsym

      pp = 0
      do p = 1 + iorb(ifsym), iorb(ifsym) + norb(ifsym)
        po = pp*nfbas(ifsym, 0) + icmoq(ifsym)
        pp = pp + 1

        kk = 0
        do k = 1 + ibas(ifsym), ibas(ifsym) + nfbas(ifsym, 0)
          kk = kk + 1

          ll = 0
          do l = 1 + ibas(ifsym), ibas(ifsym) + nfbas(ifsym, 0)
            ll = ll + 1

            factor = 0.0d0
            do iz = 1, nz
              izo = (iz-1)*nfbas(ifsym, 0)*norb(ifsym)
              factor = factor + mo_coef(kk + po + izo)*mo_coef(ll + po + izo)
            end do

            orbital_density(p) = orbital_density(p) + gao(k)*gao(l)*factor*2.0d0
          end do
        end do
      end do
    end do

    softness = 0.0d0

    do i = 1, nr_i
      do a = 1, nr_a

!       this i haven't debugged - have checked only the orbital densities
        softness = softness                                          &
                 + (orbital_density(p_i(i))*orbital_density(p_a(a))) &
                  /(eigenvalue(p_a(a)) - eigenvalue(p_i(i)))

      end do
    end do

  end subroutine

end module
