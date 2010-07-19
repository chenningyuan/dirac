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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!dirac_copyright_end

module dft_functional_derv

  use dft_cfg
#ifdef MOD_XCFUN
  use xcfun
#endif

  implicit none

  public get_functional_derv
#ifdef MOD_XCFUN
  public get_general_fun_derv
  public xcfun_set_functional
  public xcfun_is_gga
#endif

  integer, parameter, public :: d_0   =  0
  integer, parameter, public :: d_r   =  1
  integer, parameter, public :: d_z   =  2
  integer, parameter, public :: d_x   =  3
  integer, parameter, public :: d_rr  =  4
  integer, parameter, public :: d_ss  =  5
  integer, parameter, public :: d_zz  =  6
  integer, parameter, public :: d_yy  =  7
  integer, parameter, public :: d_rz  =  8
  integer, parameter, public :: d_sy  =  9
  integer, parameter, public :: d_rx  = 10
  integer, parameter, public :: d_zx  = 11
  integer, parameter, public :: d_rrr = 12
  integer, parameter, public :: d_zzz = 13
  integer, parameter, public :: d_rss = 14
  integer, parameter, public :: d_rzz = 15
  integer, parameter, public :: d_ryy = 16
  integer, parameter, public :: d_rrz = 17
  integer, parameter, public :: d_ssz = 18
  integer, parameter, public :: d_rsy = 19
  integer, parameter, public :: d_szy = 20
  integer, parameter, public :: d_zyy = 21

  integer, parameter, public :: nr_nonzero_derv(3) = (/2, 11, 21/)

  private

#include "functionals.h"

  integer, parameter   :: max_dim = 99

  real(8), allocatable :: xcfun_in(:, :)
  real(8), allocatable :: xcfun_out(:, :)
  integer              :: xcfun_in_length  = 0
  integer              :: xcfun_out_length = 0

  logical :: xcfun_id_is_set = .false.

contains

#ifdef MOD_XCFUN

  function xcfun_is_gga()

!   ----------------------------------------------------------------------------
    logical :: xcfun_is_gga
!   ----------------------------------------------------------------------------

    xcfun_is_gga = .false.
    if (xc_get_type(xc_fun%id) > 0) then
      xcfun_is_gga = .true.
    end if

  end function

   subroutine xcfun_set_functional(f,                 &
                                   response_order,    &
                                   response_order_ao, &
                                   parallel_xc)

!     --------------------------------------------------------------------------
      type(functional), intent(inout) :: f
      integer,          intent(in)    :: response_order
      integer,          intent(in)    :: response_order_ao
      logical,          intent(in)    :: parallel_xc
!     --------------------------------------------------------------------------
      integer                         :: i
      integer                         :: need_in_length
      integer                         :: need_out_length
!     --------------------------------------------------------------------------

      if (.not. xcfun_id_is_set) then
#ifdef VAR_MPI
         if (parallel_xc) then
            call sync_functional(f)
         end if
#endif
         f%id = xc_new_functional()
         call xc_set_mode(f%id, XC_VARS_RS)
       
         do i = 1, XC_NR_PARAMS
            if (dabs(f%w(i)) > tiny(0.0d0)) then
               call xc_set_param(f%id, i, f%w(i))
            end if
         end do
         xcfun_id_is_set = .true.
      end if

      need_in_length  = xc_input_length(f%id)
      need_out_length = xc_output_length(f%id, response_order)
  
      if (need_in_length > xcfun_in_length) then
         xcfun_in_length = need_in_length
         if (allocated(xcfun_in)) then
            deallocate(xcfun_in)        
         end if
         allocate(xcfun_in(xcfun_in_length, 1))
      end if
  
      if (need_out_length > xcfun_out_length) then
         xcfun_out_length = need_out_length
         if (allocated(xcfun_out)) then
            deallocate(xcfun_out)        
         end if
         allocate(xcfun_out(xcfun_out_length, mxbllen))
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
                                 derv)

!   ----------------------------------------------------------------------------
    type(functional), intent(in)  :: f
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(mxbllen, 0:n)
!   ----------------------------------------------------------------------------
    real(8)                       :: t(max_dim, mxbllen)
    integer                       :: k, l
    real(8), external             :: dftenergy
!   ----------------------------------------------------------------------------

    derv = 0.0d0

#ifdef MOD_XCFUN
    if (dft_cfg_xcfun) then

      do k = 1, block_length
        xcfun_in = 0.0d0
        xcfun_in(1, 1) = r_0(k)
!for mcscf ... 
!       here we will have to rotate to the spin density vector direction
!       xcfun_in(2) = z-component of spin density
        if (xcfun_in_length > 2) then
          xcfun_in(3, 1) = z_0(k)
        end if
      end do
      call xc_eval(f%id, order, block_length, xcfun_in, xcfun_out)

      do k = 1, block_length

        derv(k, d_0) = xcfun_out(1, k)
        if (order < 1) exit

        if (xcfun_is_gga()) then
           derv(k, d_r) = xcfun_out(XC_D10000, k)
           derv(k, d_z) = xcfun_out(XC_D00100, k)
        else
           derv(k, d_r) = xcfun_out(XC_D10, k)
        end if
        if (order < 2) exit

        if (xcfun_is_gga()) then
           derv(k, d_rr) = xcfun_out(XC_D20000, k)
           derv(k, d_ss) = xcfun_out(XC_D02000, k)
           derv(k, d_rz) = xcfun_out(XC_D10100, k)
           derv(k, d_zz) = xcfun_out(XC_D00200, k)
           derv(k, d_sy) = xcfun_out(XC_D01010, k)
           derv(k, d_yy) = xcfun_out(XC_D00020, k)
           derv(k, d_z ) = xcfun_out(XC_D00100, k)
           derv(k, d_x ) = xcfun_out(XC_D00001, k)
        else
           derv(k, d_rr) = xcfun_out(XC_D20, k)
           derv(k, d_ss) = xcfun_out(XC_D02, k)
        end if
        if (order < 3) exit

        if (xcfun_is_gga()) then
           derv(k, d_rrr) = xcfun_out(XC_D30000, k)
           derv(k, d_rss) = xcfun_out(XC_D12000, k)
           derv(k, d_rrz) = xcfun_out(XC_D20100, k)
           derv(k, d_rsy) = xcfun_out(XC_D11010, k)
           derv(k, d_rx ) = xcfun_out(XC_D10001, k)
           derv(k, d_ryy) = xcfun_out(XC_D10020, k)
           derv(k, d_rz ) = xcfun_out(XC_D10100, k)
           derv(k, d_rzz) = xcfun_out(XC_D10200, k)
           derv(k, d_ssz) = xcfun_out(XC_D02100, k)
           derv(k, d_sy ) = xcfun_out(XC_D01010, k)
           derv(k, d_szy) = xcfun_out(XC_D01110, k)
           derv(k, d_yy ) = xcfun_out(XC_D00020, k)
           derv(k, d_zx ) = xcfun_out(XC_D00101, k)
           derv(k, d_zyy) = xcfun_out(XC_D00120, k)
           derv(k, d_zz ) = xcfun_out(XC_D00200, k)
           derv(k, d_zzz) = xcfun_out(XC_D00300, k)
        else
           derv(k, d_rrr) = xcfun_out(XC_D30, k)
           derv(k, d_rss) = xcfun_out(XC_D12, k)
        end if
        if (order < 4) exit

        call quit('get_functional_derv: order too high')

      end do

    else
#endif /* ifdef MOD_XCFUN */

      select case (order)

        case (0)
          do k = 1, block_length
            derv(k, d_0) = dftenergy(r_0(k), dsqrt(z_0(k)))
          end do

        case (1)
          do k = 1, block_length
            derv(k, d_0) = dftenergy(r_0(k), dsqrt(z_0(k)))
            call dftpot0(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (2)
          do k = 1, block_length
            call sdftpot1(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case (3)
          do k = 1, block_length
            call sdftpot2(t(1, k), 1.0d0, r_0(k), dsqrt(z_0(k)))
          end do

        case default
          call quit('get_functional_derv: order too high')

      end select

      do l = 1, n
        do k = 1, block_length
          derv(k, l) = derv(k, l) + t(l, k)
        end do
      end do

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
                                  parallel_xc)
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
    real(8),          intent(out) :: derv(mxbllen, 0:n)
      logical,          intent(in)    :: parallel_xc
    real(8) :: dervtmp(mxbllen, 0:n)
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

!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): after nokern first xc' 
!   call debug_general_fun(fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0   ,                &
!                             z_0   ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '
                               
    else if ( dirac_cfg_fde_response.and.have_frozen_density ) then
! aspg, 2010-03-29      
! in case we implement the use different xc for the total and active
! densities, we have to uncomment the line below, and suitably define
! tot_xc_fun in fde_mod (or some other place)
!      call xc_fun_set_functional(order, order, tot_xc_fun)
! but also replace fun below by tot_xc_fun     
      call get_functional_derv(f,                &
                               order,              &
                               block_length,       &
                               w,                  &
                               r_0_t,              &
                               z_0_t,              &
                               n,                  &
                               dervtmp)

      do i = 1, mxbllen
        do j = 0, n
           derv(i,j) = derv(i,j) + dervtmp(i,j)
         enddo
      enddo
    
!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): after first xc' 
!   call debug_general_fun(fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0_t  ,                &
!                             z_0_t  ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '
    
      if ( .not.dirac_cfg_dft_calculation ) then
! if we are not doing dft, we subtract what would be the xc contribution from 
! the active subsystem from the response.

! aspg, 2010-03-29      
! in case we implement the use of different xc for the total and active
! densities, we have to uncomment the line below
!      call xc_fun_set_functional(order, order, fun)
        call get_functional_derv(f,                 &
                                 order,              &
                                 block_length,       &
                                 w,                  &
                                 r_0,                &
                                 z_0,                &
                                 n,                  &
                                 dervtmp)
                                 
        do i = 1, mxbllen
          do j = 0, n
           derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
        enddo

!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): after second xc' 
!   call debug_general_fun(fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0  ,                &
!                             z_0  ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '

      endif

! now we construct the contribution from the kinetic energy to the response
! of the embedding potential. 
! first set functional to the of-kef, do the non-additive ke contribution
      call xcfun_set_functional(kin_fun, order, order, parallel_xc)
      
      call get_functional_derv(kin_fun,            &
                               order,              &
                               block_length,       &
                               w,                  &
                               r_0_t,              &
                               z_0_t,              &
                               n,                  &
                               dervtmp)

      do i = 1, mxbllen
        do j = 0, n
           derv(i,j) = derv(i,j) + dervtmp(i,j)
         enddo
      enddo
      
!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): after first kinetic' 
!   call debug_general_fun(id_kin_fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0_t ,                &
!                             z_0_t ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '

! now subtract the active subsystem's kinetic energy contribution

      call get_functional_derv(kin_fun,            &
                               order,              &
                               block_length,       &
                               w,                  &
                               r_0,                &
                               z_0,                &
                               n,                  &
                               dervtmp)
                               
      do i = 1, mxbllen
        do j = 0, n
           derv(i,j) = derv(i,j) - dervtmp(i,j)
         enddo
      enddo

!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): after second kinetic' 
!   call debug_general_fun(id_kin_fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0 ,                &
!                             z_0 ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '
          
! and restore the default dft functional, if we are to use it 
! down the line...
      call xcfun_set_functional(f, order, order, parallel_xc)

    endif

!   write (*,*) ' '
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) 'DEBUG (functional derivative): at the end' 
!   call debug_general_fun(fun,               &
!                             order,              &
!                             block_length,       &
!                             w,                  &
!                             r_0 ,                &
!                             z_0 ,                &
!                             n,                  &
!                             derv, dervtmp)
!   write (*,*) 'DEBUG (functional derivative)'
!   write (*,*) ' '
!   call quit('debugging fde response. stop here.')

  end subroutine get_general_fun_derv





  subroutine debug_general_fun(order,        &
                               block_length, &
                               w,            &
                               r_0,          &
                               z_0,          &
                               n,            &
                               derv, dervtmp)
   use dirac_cfg
   use fde_mod

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: order
    integer, intent(in)  :: block_length
    real(8), intent(in)  :: w(*)
    real(8), intent(in)  :: r_0(*)
    real(8), intent(in)  :: z_0(*)
    integer, intent(in)  :: n
    real(8), intent(out) :: derv(mxbllen, 0:n)
    real(8) :: dervtmp(mxbllen, 0:n)
    integer :: i, j
!   ----------------------------------------------------------------------------

      write (*,*) ''
      write (*,*) '   functional definition'
!radovan: commented this out - we don't handle the weights outside xcfun
!         anymore
!     do i = 1, XC_PARAMS_LEN
!        if (fun%w(i) .gt. 0.0d0) then
!           write (*,*) '      functional type :',i
!           write (*,*) '            weight fulld =',fun%w(i)
!           write (*,*) '            weight alda  =',fun%w_alda(i)
!           write (*,*) '            weight xalda =',fun%w_xalda(i)
!        endif
!     enddo
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
#endif

end module
