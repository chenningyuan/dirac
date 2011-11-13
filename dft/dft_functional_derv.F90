!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

module dft_functional_derv

  use dft_cfg

  implicit none

  public get_functional_derv

  integer, parameter, public :: d0000000 =  0
  integer, parameter, public :: d1000000 =  1
  integer, parameter, public :: d0010000 =  2
  integer, parameter, public :: d0000100 =  3
  integer, parameter, public :: d2000000 =  4
  integer, parameter, public :: d0200000 =  5
  integer, parameter, public :: d0020000 =  6
  integer, parameter, public :: d0002000 =  7
  integer, parameter, public :: d1010000 =  8
  integer, parameter, public :: d0101000 =  9
  integer, parameter, public :: d1000100 = 10
  integer, parameter, public :: d0010100 = 11
  integer, parameter, public :: d3000000 = 12
  integer, parameter, public :: d0030000 = 13
  integer, parameter, public :: d1200000 = 14
  integer, parameter, public :: d1020000 = 15
  integer, parameter, public :: d1002000 = 16
  integer, parameter, public :: d2010000 = 17
  integer, parameter, public :: d0210000 = 18
  integer, parameter, public :: d1101000 = 19
  integer, parameter, public :: d0111000 = 20
  integer, parameter, public :: d0012000 = 21

  integer, parameter, public :: d0000010 = 22
  integer, parameter, public :: d0000020 = 23
  integer, parameter, public :: d1000010 = 24
  integer, parameter, public :: d0010010 = 25
  integer, parameter, public :: d0100001 = 26
  integer, parameter, public :: d0001001 = 27
  integer, parameter, public :: d0000002 = 28

  integer, parameter, public :: nr_nonzero_derv = 28

  private

#include "functionals.h"

  integer, parameter   :: max_dim = 99

contains


  subroutine get_functional_derv(f,            &
                                 order,        &
                                 block_length, &
                                 w,            &
                                 r_0,          &
                                 z_0,          &
                                 n,            &
                                 derv,         &
                                 alda_real, &
                                 alda_imag, &
                                 tau)

!   ----------------------------------------------------------------------------
    type(functional)              :: f
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(mxbllen, 0:n)
    logical, optional, intent(in) :: alda_real
    logical, optional, intent(in) :: alda_imag
    real(8), optional, intent(in) :: tau
!   ----------------------------------------------------------------------------
    real(8)                       :: t1(max_dim, mxbllen)
    real(8)                       :: t2(max_dim, mxbllen)
    integer                       :: k, l, id
    real(8), external             :: dftenergy
    logical                       :: do_alda_real
    logical                       :: do_alda_imag
!   ----------------------------------------------------------------------------

    do_alda_real = .false.
    do_alda_imag = .false.

    if (present(alda_real)) do_alda_real = alda_real
    if (present(alda_imag)) do_alda_imag = alda_imag

!   if less then linear response ignore alda
    if (order < 2) then
       do_alda_real = .false.
       do_alda_imag = .false.
    end if

    derv = 0.0d0
    t1   = 0.0d0
    t2   = 0.0d0


      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d0000000) = dftenergy(r_0(k), dsqrt(z_0(k)))
          end do

        case (1)
          do k = 1, block_length
            derv(k, d0000000) = dftenergy(r_0(k), dsqrt(z_0(k)))
            call dftpot0(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (2)
          do k = 1, block_length
            call  dftpot1(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
            call sdftpot1(t2(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (3)
          do k = 1, block_length
            call sdftpot2(t1(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

      do l = 1, n
        do k = 1, block_length
          derv(k, l) = derv(k, l) + t1(l, k)
        end do
      end do
      if (order == 2) then
         do l = 1, n
           do k = 1, block_length
             derv(k, l) = derv(k, l) + t2(l, k)
           end do
         end do
      end if


!   multiply by weight
    do l = 0, n
       do k = 1, block_length
          derv(k, l) = w(k)*derv(k, l)
       end do
    end do

  end subroutine




end module
