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

module fde_xcfun_interface

   use fde_cfg
   use dft_cfg
   use dft_input
   use xcfun

   private

   public fde_calc_nadd_derv
   public fde_initializeFunctionals
   public fde_printFunctionals
   public set_fde_alda
   public is_fde_alda

!  public fde_kef,       fde_xcf
!  public fde_kef_alda,  fde_xcf_alda
!  public fde_kef_xalda, fde_xcf_xalda

   type(functional), save :: fde_kef,       fde_xcf
   type(functional), save :: fde_kef_alda,  fde_xcf_alda
   type(functional), save :: fde_kef_xalda, fde_xcf_xalda

! the default (orbital-free) kinetic energy functional for now is the thomas-fermi one
   character(len=80), save :: kefun_fde='kin_tf'
   character(len=80), save :: xcfun_fde='lda'

   logical, save :: fde_set_functionals(20)

   contains

      subroutine parse_xcf(input_unit)
         integer :: input_unit
         read(input_unit,'(a80)') xcfun_fde
         call parse_functional(xcfun_fde,fde_xcf)
      end subroutine parse_xcf

      subroutine parse_kef(input_unit)
         character(len=*) :: input_unit
         read(input_unit,'(a80)') kefun_fde
         call parse_functional(kefun_fde,fde_kef)
      end subroutine parse_kef

      subroutine set_nadd_alda()
         call set_kef_alda(fde_kef_alda)
         call set_xcf_alda(fde_xcf_alda)
      end subroutine set_nadd_alda


      subroutine fde_printFunctionals
         write(log_unit,'(A,A20)') ' * Non-additive KE functional : ',kefun_fde
         write(log_unit,'(A,A20)') ' * Non-additive XC functional : ',xcfun_fde
!        write(log_unit,'(A,A20)') '   if .ALDA is set, will use : ',kefun_alda_fde
     end subroutine fde_printFunctionals 


      subroutine fde_initializeFunctionals
         integer :: i

         reset_hfxfac = .false.

         call parse_functional(kefun_fde,fde_kef)
         call parse_functional(kefun_fde,fde_kef_alda)
!        call parse_functional(kefun_fde,fde_kef_xalda)

! aspg: in order to be able to set the xc functionals used in the
!       non-additive xc potential, we need to first let xcfun set more
!       than 5 functionals...
         call parse_functional(xcfun_fde,fde_xcf)
         call parse_functional(xcfun_fde,fde_xcf_alda)
         call parse_functional(xcfun_fde,fde_xcf_xalda)

         do i = 1, 20
            fde_set_functionals(i) = .true.
         enddo

!        write (*,*) 'debug bla'
!        call set_kin_fun_alda(fde_kef_alda)
!        call set_kin_fun_alda(fde_kef_xalda)
!        call set_xc_fun_alda(fde_xcf_alda)
!        call set_xc_fun_alda(fde_xcf_xalda)
!        write (*,*) 'debug bla'

         reset_hfxfac = .true.

      end subroutine fde_initializeFunctionals






      subroutine fde_calc_nadd_derv(f,           &
                                   f_alda,      &
                                   f_xalda, &
!
! subroutine get_fde_nadd_fun_derv(order,        &
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
   use xc_derv
   use xc_max_block_length

!   ----------------------------------------------------------------------------
    type(functional), intent(inout)  :: f, f_alda, f_xalda
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    real(8),          intent(in)  :: r_0_t(*)
    real(8),          intent(in)  :: z_0_t(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(max_block_length, 0:n)
    logical,          intent(in)  :: parallel_xc
    real(8)                       :: dervtmp(max_block_length, 0:n)
    integer :: i, j
!   ----------------------------------------------------------------------------

    derv = 0.0d0

    if (fde_set_functionals(order)) then
!      call xcfun_set_functional(fde_xcf     , order, parallel_xc) 
!      call xcfun_set_functional(fde_xcf_alda, order, parallel_xc) 
!      call xcfun_set_functional(fde_xcf_xalda,order, parallel_xc) 
       call xcfun_set_functional(fde_kef,      order, parallel_xc)
       call xcfun_set_functional(fde_kef_alda, order, parallel_xc)
!      call xcfun_set_functional(fde_kef_xalda,order, parallel_xc)
       fde_set_functionals(order) = .false.
   endif

    if (.not.fde_skip_xc_contribution) then

!      call get_functional_derv(fde_xcf,       &
!                               fde_xcf_alda,  &
!                               fde_xcf_xalda, &
       call get_functional_derv(f,       &
                                f_alda,  &
                                f_xalda, &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0_t,         &
                                z_0_t,         &
                                n,             &
                                dervtmp,       &
                                alda_real=fde_is_rsp_alda, &
                                alda_imag=fde_is_rsp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) + dervtmp(i,j)
          enddo
       enddo

!      call get_functional_derv(fde_xcf,       &
!                               fde_xcf_alda,  &
!                               fde_xcf_xalda, &
       call get_functional_derv(f,       &
                                f_alda,  &
                                f_xalda, &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                z_0,           &
                                n,             &
                                dervtmp,       &
                                alda_real=fde_is_rsp_alda, &
                                alda_imag=fde_is_rsp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
       enddo

    endif

! now we construct the contribution from the kinetic energy to the response
! of the embedding potential.
! first set functional to the of-kef, do the non-additive ke contribution

    if (.not.fde_skip_kin_contribution) then

       call get_functional_derv(fde_kef,       &
                                fde_kef_alda,  &
                                fde_kef_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0_t,         &
                                z_0_t,         &
                                n,             &
                                dervtmp,       &
                                alda_real=fde_is_rsp_alda, &
                                alda_imag=fde_is_rsp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) + dervtmp(i,j)
          enddo
       enddo

! now subtract the active subsystem's kinetic energy contribution

       call get_functional_derv(fde_kef,       &
                                fde_kef_alda,  &
                                fde_kef_alda,  &
                                order,         &
                                block_length,  &
                                w,             &
                                r_0,           &
                                z_0,           &
                                n,             &
                                dervtmp,       &
                                alda_real=fde_is_rsp_alda, &
                                alda_imag=fde_is_rsp_alda)

       do i = 1, max_block_length
          do j = 0, n
             derv(i,j) = derv(i,j) - dervtmp(i,j)
          enddo
       enddo

    endif

! ! and restore the default dft functional, if we are to use it
! ! down the line...
!     call xcfun_set_functional(f, order, parallel_xc)
!     call xcfun_set_functional(f_alda, order, parallel_xc)
!     call xcfun_set_functional(f_xalda, order, parallel_xc)

   end subroutine fde_calc_nadd_derv


   subroutine fde_calc_ngn_total(n_tot,gn_tot)
      real(kind=8) :: n_tot(*), gn_tot(*) 
      real(kind=8) :: gx, gy, gz
      integer :: i

      do i = 1, fde_grid_np 
        
         if (fde_frozen_n(i) > fde_cfg_screening ) then 
            n_tot(i) = fde_active_n(i) + fde_frozen_n(i) 
         
            gx = fde_active_gn(i,1) + fde_frozen_gn(i,1) 
            gy = fde_active_gn(i,2) + fde_frozen_gn(i,2) 
            gz = fde_active_gn(i,3) + fde_frozen_gn(i,3) 

            gn_tot(i) = gx*gx + gy*gy + gz*gz
         else
            n_tot(i)  = fde_active_n(i) 
            gn_tot(i) = fde_frozen_gn(i,1)*fde_frozen_gn(i,1) &
                      + fde_frozen_gn(i,2)*fde_frozen_gn(i,2) &
                      + fde_frozen_gn(i,3)*fde_frozen_gn(i,3) 
         endif

      enddo

   end subroutine fde_calc_ngn_total


   subroutine fde_calc_emb_derv()
      real(kind=8), pointer :: p
      real(kind=8), allocatable :: n_tot(:), gnn_tot(:)
      real(kind=8) :: w, n0, gn0, n0t, gn0t
      integer :: i

      call set_pointer_fde_grid(p)

      call alloc(n_tot,fde_grid_np)
      call alloc(gn_tot,fde_grid_np)

      call fde_calc_ngn_total(n_tot,gnn_tot)

      do i = 1, fde_grid_np 
         n0   = fde_active_n(i)
         gn0  = fde_active_gn(i,1)*fde_active_gn(i,1) &
              + fde_active_gn(i,2)*fde_active_gn(i,2) &
              + fde_active_gn(i,3)*fde_active_gn(i,3) 

         n0t  = n_tot(i) 
         gn0t = gnn_tot(i) 

         w = p(i,4)

         call set_pointer_fde_derv(fde_derv,i)
      
         call get_fde_nadd_fun_derv(max_fun_derv, &
                                    bllen,        &
                                    w,            &
                                    n0,           &
                                    n0t,          &
                                    gn0,          &
                                    gn0t,         &
                                    derv_length,  &
                                    derv_nadd,    &
                                    parallel_xc)

      enddo

      call dealloc(n_tot)
      call dealloc(gn_tot)

   end subroutine fde_calc_emb_derv

! we are assuming dirac
   subroutine fde_add_emb_derv(ipoint,derv)
      real(kind=8) :: derv(:,:)
      integer :: i, j, ipoint

      w = p(ipoint,4)

      do i = 1, max_block_length
         do j = 0, derv_length
            if (j >  1 ) &
               derv(i,j) = derv(i,j) + derv_nadd(ipoint,i,j)
! we add the electrostatic potential from the frozen subsystem(s) to the derivative
            if (j == 1 .and. update_vemb) &
               derv(i,j) = derv(i,j) + derv_nadd(ipoint,i,j) + w*fde_frozen_elpot(ipoint)
         enddo
      enddo

   end subroutine

end module fde_xcfun_interface
