!dirac_copyright_start
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
!dirac_copyright_end

module xc_derv

  use dft_cfg
#ifdef MOD_XCFUN
  use xcfun
#endif
  use xc_max_block_length

  implicit none

  public get_functional_derv
#ifdef MOD_XCFUN
  public get_general_fun_derv
  public xcfun_set_functional
  public xcfun_is_lda
  public xcfun_is_gga
  public xcfun_is_tau_mgga
#endif

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

#include "../../dft/functionals.h"

  integer, parameter   :: max_dim = 99

contains

#ifdef MOD_XCFUN

   function xcfun_is_lda()
      logical :: xcfun_is_lda

      xcfun_is_lda = .false.
      if (xc_get_type(xc_fun%id) == 0) then
         xcfun_is_lda = .true.
      end if
   end function

   function xcfun_is_gga()
      logical :: xcfun_is_gga

      xcfun_is_gga = .false.
      if (xc_get_type(xc_fun%id) == 1) then
         xcfun_is_gga = .true.
      end if
   end function

   function xcfun_is_tau_mgga()
      logical :: xcfun_is_tau_mgga

      xcfun_is_tau_mgga = .false.
      if (xc_get_type(xc_fun%id) == 2) then
         xcfun_is_tau_mgga = .true.
      end if
   end function





   subroutine xcfun_set_functional(f,              &
                                   response_order, &
                                   parallel_xc)

!     --------------------------------------------------------------------------
      type(functional), intent(inout) :: f
      integer,          intent(in)    :: response_order
      logical,          intent(in)    :: parallel_xc
!     --------------------------------------------------------------------------
      integer                         :: i, l
!     --------------------------------------------------------------------------

      if (.not. f%xcfun_id_is_set) then
#ifdef VAR_MPI
         if (parallel_xc) then
            call sync_functional(f)
         end if
#endif
         f%id = xc_new_functional()
         call xc_set_mode(f%id, XC_VARS_NS)

         do i = 1, XC_NR_PARAMS
            if (dabs(f%w(i)) > tiny(0.0d0)) then
               call xc_set_param(f%id, i, f%w(i))
            end if
         end do
         f%xcfun_id_is_set = .true.
      end if

      if (allocated(f%input)) then
         l = size(f%input)
      else
         l = 0
      end if
      if (xc_input_length(f%id) > l) then
         l = xc_input_length(f%id)
         if (allocated(f%input)) then
            deallocate(f%input)
         end if
         allocate(f%input(l, 1))
      end if

      if (allocated(f%output)) then
         l = size(f%output)
      else
         l = 0
      end if
      if (xc_output_length(f%id, 1 + response_order) > l) then
         l = xc_output_length(f%id, 1 + response_order)
         if (allocated(f%output)) then
            deallocate(f%output)
         end if
         allocate(f%output(l, 1))
      end if

   end subroutine
#endif /* ifdef MOD_XCFUN */

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
    real(8),          intent(out) :: derv(max_block_length, 0:n)
    logical, optional, intent(in) :: alda_real
    logical, optional, intent(in) :: alda_imag
    real(8), optional, intent(in) :: tau
!   ----------------------------------------------------------------------------
    real(8)                       :: t1(max_dim, max_block_length)
    real(8)                       :: t2(max_dim, max_block_length)
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

#ifdef MOD_XCFUN
    if (dft_cfg_xcfun) then

      f%input = 0.0d0
      do k = 1, block_length
         f%input(1, k) = r_0(k)
      end do
      if (xcfun_is_gga()) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
      end if
      if (xcfun_is_tau_mgga()) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
         do k = 1, block_length
            f%input(6, k) = tau
         end do
      end if

      call xc_eval(f%id,         &
                   order,        &
                   block_length, &
                   f%input,      &
                   f%output)

      if (do_alda_real .or. do_alda_imag) then
         xc_fun_alda%input = 0.0d0
         do k = 1, block_length
            xc_fun_alda%input(1, k) = r_0(k)
         end do
         id = xc_fun_alda%id
         if (dft_cfg_xalda) then
            id = xc_fun_xalda%id
         end if
         call xc_eval(id,                &
                      order,             &
                      block_length,      &
                      xc_fun_alda%input, &
                      xc_fun_alda%output)
      end if

      do k = 1, block_length

         derv(k, d0000000) = f%output(1, k)

         select case (order)

            case (0)
               call quit('get_functional_derv: order 0; probably programming error')

            case (1)
               if (xcfun_is_lda()) then
                  derv(k, d1000000) = f%output(XC_D10, k)
               else if (xcfun_is_gga()) then
                  derv(k, d1000000) = f%output(XC_D10000, k)
                  derv(k, d0010000) = f%output(XC_D00100, k)
               else if (xcfun_is_tau_mgga()) then
                  derv(k, d1000000) = f%output(XC_D1000000, k)
                  derv(k, d0010000) = f%output(XC_D0010000, k)
                  derv(k, d0000010) = f%output(XC_D0000010, k)
               else
                  call quit('get_functional_derv: unknown fun type')
               end if

            case (2)
!              order 2 - spin unpolarized
               if (do_alda_real) then
                  derv(k, d2000000) = xc_fun_alda%output(XC_D20, k)
               else
                  if (xcfun_is_lda()) then
                     derv(k, d2000000) = f%output(XC_D20, k)
                  else if (xcfun_is_gga()) then
                     derv(k, d2000000) = f%output(XC_D20000, k)
                     derv(k, d1010000) = f%output(XC_D10100, k)
                     derv(k, d0020000) = f%output(XC_D00200, k)
                     derv(k, d0010000) = f%output(XC_D00100, k)
                  else if (xcfun_is_tau_mgga()) then
                     derv(k, d0010000) = f%output(XC_D0010000, k)
                     derv(k, d2000000) = f%output(XC_D2000000, k)
                     derv(k, d1010000) = f%output(XC_D1010000, k)
                     derv(k, d0020000) = f%output(XC_D0020000, k)

                     derv(k, d0000020) = f%output(XC_D0000020, k)
                     derv(k, d0010010) = f%output(XC_D0010010, k)
                     derv(k, d1000010) = f%output(XC_D1000010, k)
                  else
                     call quit('get_functional_derv: unknown fun type')
                  end if
               end if
!              order 2 - spin polarized
               if (do_alda_imag) then
                  derv(k, d0200000) = xc_fun_alda%output(XC_D02, k)
               else
                  if (xcfun_is_lda()) then
                     derv(k, d0200000) = f%output(XC_D02, k)
                  else if (xcfun_is_gga()) then
                     derv(k, d0200000) = f%output(XC_D02000, k)
                     derv(k, d0101000) = f%output(XC_D01010, k)
                     derv(k, d0002000) = f%output(XC_D00020, k)
                     derv(k, d0000100) = f%output(XC_D00001, k)
                  else if (xcfun_is_tau_mgga()) then
                     derv(k, d0000100) = f%output(XC_D0000100, k)
                     derv(k, d0200000) = f%output(XC_D0200000, k)
                     derv(k, d0101000) = f%output(XC_D0101000, k)
                     derv(k, d0002000) = f%output(XC_D0002000, k)

                     derv(k, d0000002) = f%output(XC_D0000002, k)
                     derv(k, d0001001) = f%output(XC_D0001001, k)
                     derv(k, d0100001) = f%output(XC_D0100001, k)
                  else
                     call quit('get_functional_derv: unknown fun type')
                  end if
               end if

            case (3)
               if (do_alda_real .or. do_alda_imag) then
                  call quit('get_functional_derv: xcfun alda qr not implemented')
               end if
               if (xcfun_is_lda()) then
                  derv(k, d2000000) = f%output(XC_D20, k)
                  derv(k, d0200000) = f%output(XC_D02, k)
                  derv(k, d1200000) = f%output(XC_D12, k)
                  derv(k, d3000000) = f%output(XC_D30, k)
               else if (xcfun_is_gga()) then
                  derv(k, d0010000) = f%output(XC_D00100, k)
                  derv(k, d0000100) = f%output(XC_D00001, k)
                  derv(k, d2000000) = f%output(XC_D20000, k)
                  derv(k, d1010000) = f%output(XC_D10100, k)
                  derv(k, d1000100) = f%output(XC_D10001, k)
                  derv(k, d0200000) = f%output(XC_D02000, k)
                  derv(k, d0101000) = f%output(XC_D01010, k)
                  derv(k, d0020000) = f%output(XC_D00200, k)
                  derv(k, d0010100) = f%output(XC_D00101, k)
                  derv(k, d0002000) = f%output(XC_D00020, k)
                  derv(k, d3000000) = f%output(XC_D30000, k)
                  derv(k, d2010000) = f%output(XC_D20100, k)
                  derv(k, d1200000) = f%output(XC_D12000, k)
                  derv(k, d1101000) = f%output(XC_D11010, k)
                  derv(k, d1020000) = f%output(XC_D10200, k)
                  derv(k, d1002000) = f%output(XC_D10020, k)
                  derv(k, d0210000) = f%output(XC_D02100, k)
                  derv(k, d0111000) = f%output(XC_D01110, k)
                  derv(k, d0030000) = f%output(XC_D00300, k)
                  derv(k, d0012000) = f%output(XC_D00120, k)
               else if (xcfun_is_tau_mgga()) then
                  call quit('get_functional_derv: tau mgga qr not implemented')
               else
                  call quit('get_functional_derv: unknown fun type')
               end if

            case default
               call quit('get_functional_derv: order too high')

         end select

      end do

    else
#endif /* ifdef MOD_XCFUN */

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

#ifdef MOD_XCFUN
    end if
#endif

!   multiply by weight
    do l = 0, n
       do k = 1, block_length
          derv(k, l) = w(k)*derv(k, l)
       end do
    end do

  end subroutine



#ifdef MOD_XCFUN
  subroutine get_general_fun_derv(f,            &
                                  order,        &
                                  block_length, &
                                  w,            &
                                  r_0,          &
                                  r_0_t,        &
                                  z_0,          &
                                  z_0_t,        &
                                  n,            &
                                  derv,         &
                                  parallel_xc,  &
                                  debug)
   use dirac_cfg
   use fde_mod

!   ----------------------------------------------------------------------------
    type(functional), intent(inout)  :: f
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    real(8),          intent(in)  :: r_0_t(*)
    real(8),          intent(in)  :: z_0_t(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(max_block_length, 0:n)
      logical,          intent(in)    :: parallel_xc
    logical,          intent(in)  :: debug
    real(8) :: dervtmp(max_block_length, 0:n)
    integer :: i, j
!   ----------------------------------------------------------------------------

    derv = 0.0d0

    if ( .not.have_frozen_density ) then
      call get_functional_derv(f,                  &
                               order,              &
                               block_length,       &
                               w,                  &
                               r_0,                &
                               z_0,                &
                               n,                  &
                               derv)

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): after nokern first xc'
    call debug_general_fun(f,               &
                           order,           &
                           block_length,    &
                           w,               &
                           r_0   ,          &
                           z_0   ,          &
                           n,               &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif

    else if ( dirac_cfg_fde_response.and.have_frozen_density ) then

! aspg, 2010-03-29
! in case we implement the use different xc for the total and active
! densities, we have to uncomment the line below, and suitably define
! tot_xc_fun in fde_mod (or some other place)
!      call xc_fun_set_functional(order, order, tot_xc_fun)
! but also replace fun below by tot_xc_fun

      if (.not.skip_xc_contribution) then
         call get_functional_derv(f,                &
                                  order,            &
                                  block_length,     &
                                  w,                &
                                  r_0_t,            &
                                  z_0_t,            &
                                  n,                &
                                  dervtmp)
      endif

      do i = 1, max_block_length
        do j = 0, n
           derv(i,j) = derv(i,j) + dervtmp(i,j)
         enddo
      enddo

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): after first xc'
    call debug_general_fun(f,              &
                           order,          &
                           block_length,   &
                           w,              &
                           r_0_t  ,        &
                           z_0_t  ,        &
                           n,              &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif

      if ( .not.dirac_cfg_dft_calculation ) then
! if we are not doing dft, we subtract what would be the xc contribution from
! the active subsystem from the response.

! aspg, 2010-03-29
! in case we implement the use of different xc for the total and active
! densities, we have to uncomment the line below
!      call xc_fun_set_functional(order, order, fun)

        if (.not.skip_xc_contribution) then
           call get_functional_derv(f,                 &
                                    order,             &
                                    block_length,      &
                                    w,                 &
                                    r_0,               &
                                    z_0,               &
                                    n,                 &
                                    dervtmp)
        endif

        do i = 1, max_block_length
          do j = 0, n
           derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
        enddo

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): after second xc'
    call debug_general_fun(f,            &
                           order,        &
                           block_length, &
                           w,            &
                           r_0  ,        &
                           z_0  ,        &
                           n,            &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif

      endif

! now we construct the contribution from the kinetic energy to the response
! of the embedding potential.
! first set functional to the of-kef, do the non-additive ke contribution
      call xcfun_set_functional(fde_kef, order, parallel_xc)

      if (.not.skip_kin_contribution) then
         call get_functional_derv(fde_kef,            &
                                  order,              &
                                  block_length,       &
                                  w,                  &
                                  r_0_t,              &
                                  z_0_t,              &
                                  n,                  &
                                  dervtmp)
      endif

      do i = 1, max_block_length
        do j = 0, n
           derv(i,j) = derv(i,j) + dervtmp(i,j)
         enddo
      enddo

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): after first kinetic'
    call debug_general_fun(fde_kef,         &
                           order,           &
                           block_length,    &
                           w,               &
                           r_0_t ,          &
                           z_0_t ,          &
                           n,               &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif

! now subtract the active subsystem's kinetic energy contribution

      if (.not.skip_kin_contribution) then
         call get_functional_derv(fde_kef,       &
                                  order,         &
                                  block_length,  &
                                  w,             &
                                  r_0,           &
                                  z_0,           &
                                  n,             &
                                  dervtmp)
      endif

      do i = 1, max_block_length
        do j = 0, n
           derv(i,j) = derv(i,j) - dervtmp(i,j)
         enddo
      enddo

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): after second kinetic'
    call debug_general_fun(fde_kef,         &
                           order,           &
                           block_length,    &
                           w,               &
                           r_0 ,            &
                           z_0 ,            &
                           n,               &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif

! and restore the default dft functional, if we are to use it
! down the line...
      call xcfun_set_functional(f, order, parallel_xc)

    endif

    if (debug) then
    write (*,*) ' '
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) 'DEBUG (functional derivative): at the end'
    call debug_general_fun(f,               &
                           order,           &
                           block_length,    &
                           w,               &
                           r_0 ,            &
                           z_0 ,            &
                           n,               &
                           derv, dervtmp)
    write (*,*) 'DEBUG (functional derivative)'
    write (*,*) ' '
    endif
!   call quit('debugging fde response. stop here.')

  end subroutine get_general_fun_derv





  subroutine debug_general_fun(f,            &
                               order,        &
                               block_length, &
                               w,            &
                               r_0,          &
                               z_0,          &
                               n,            &
                               derv, dervtmp)
   use dirac_cfg
   use fde_mod

!   ----------------------------------------------------------------------------
    type(functional), intent(inout)  :: f
    integer, intent(in)  :: order
    integer, intent(in)  :: block_length
    real(8), intent(in)  :: w(*)
    real(8), intent(in)  :: r_0(*)
    real(8), intent(in)  :: z_0(*)
    integer, intent(in)  :: n
    real(8), intent(out) :: derv(max_block_length, 0:n)
    real(8) :: dervtmp(max_block_length, 0:n)
    integer :: i, j
!   ----------------------------------------------------------------------------

      write (*,*) ''
      write (*,*) '   functional definition'
      do i = 1, XC_NR_PARAMS
         if (f%w(i) .gt. 0.0d0) then
            write (*,*) '      functional type :',i
            write (*,*) '            weight =',f%w(i)
         endif
      enddo
      write (*,*) ' '
      write (*,*) '   density and functional derivative values'
      do i = 1, 1
        write (*,*) '      r =',r_0(i),'  z =',z_0(i),'  w =',w(i)
        do j = 0, n
           write (*,*) '      derv(',i,',',j,') =',derv(i,j)
           write (*,*) '      dtmp(',i,',',j,') =',dervtmp(i,j)
           write (*,*) '      sum =',(derv(i,j)+dervtmp(i,j)),'  difference =',(derv(i,j)-dervtmp(i,j))
         enddo
      enddo
  end subroutine debug_general_fun
#endif /* ifdef MOD_XCFUN */

end module
