
module fde_nadd_derv

      use dirac_cfg   ! why?
      use xcfun
      use interface_functional_read
      use xc_derv
      use xc_max_block_length


      public fde_initialize_nadd_functionals
      public fde_print_nadd_functionals
      public fde_set_skip_nadd_xc
      public fde_set_skip_nadd_ke
      public fde_set_nadd_kef
      public fde_set_nadd_xcf
      public fde_set_nadd_all_alda
      public get_fde_nadd_fun_derv
      public fde_xc_fun_is_lda
      public fde_xc_fun_is_gga
      public fde_ke_fun_is_lda
      public fde_ke_fun_is_gga

      logical,save :: skip_kin_contribution = .false.
      logical,save :: skip_xc_contribution  = .false.
      logical,save :: is_fdersp_alda        = .false.

      type(functional), save :: fde_kef,       fde_xcf
      type(functional), save :: fde_kef_alda,  fde_xcf_alda
      type(functional), save :: fde_kef_xalda, fde_xcf_xalda

! the default (orbital-free) kinetic energy functional for now is the thomas-fermi one
!      character(len=80), save :: kefun_fde='kin_tf'
!      character(len=80), save :: xcfun_fde='lda'
       character(len=80), save :: kefun_fde='kin_pw91'
       character(len=80), save :: xcfun_fde='pbe'

      character(len=80), save :: kefun_fde_alda='kin_tf'
      character(len=80), save :: xcfun_fde_alda='lda'

      logical, save :: fde_set_functionals(20) 

      real(8) :: hfx_out, mu_out, beta_out
      
      contains

!   ----------------------------------------------------------------------------
      function fde_xc_fun_is_lda()
!   ----------------------------------------------------------------------------
      logical :: fde_xc_fun_is_lda

         fde_xc_fun_is_lda = (xc_get_type(fde_xcf%id) == 0)
      end function


!   ----------------------------------------------------------------------------
      function fde_xc_fun_is_gga()
!   ----------------------------------------------------------------------------
      logical :: fde_xc_fun_is_gga

         fde_xc_fun_is_gga = (xc_get_type(fde_xcf%id) == 1)
      end function
   

!   ----------------------------------------------------------------------------
      function fde_ke_fun_is_lda()
!   ----------------------------------------------------------------------------
      logical :: fde_ke_fun_is_lda

         fde_ke_fun_is_lda = (xc_get_type(fde_kef%id) == 0)
      end function


!   ----------------------------------------------------------------------------
      function fde_ke_fun_is_gga()
!   ----------------------------------------------------------------------------
      logical :: fde_ke_fun_is_gga

         fde_ke_fun_is_gga = (xc_get_type(fde_kef%id) == 1)
      end function
      
         
!   ----------------------------------------------------------------------------
         subroutine fde_set_skip_nadd_xc(value)
!   ----------------------------------------------------------------------------
            logical :: value
            skip_xc_contribution  = value
         end subroutine 
    

!   ----------------------------------------------------------------------------
         subroutine fde_set_skip_nadd_ke(value)
!   ----------------------------------------------------------------------------
            logical :: value
            skip_kin_contribution  = value
         end subroutine 


!   ----------------------------------------------------------------------------
         subroutine fde_initialize_nadd_functionals
!   ----------------------------------------------------------------------------
            integer :: i, order
            logical :: parallel_fde
            
            parallel_fde = .false.
            order = 0

            fun_is_automatic = .true. 
            
            call parse_functional(kefun_fde,      fde_kef,      hfx_out, mu_out, beta_out, .false.)
            call parse_functional(kefun_fde_alda, fde_kef_alda, hfx_out, mu_out, beta_out, .false.)
            call parse_functional(xcfun_fde,      fde_xcf,      hfx_out, mu_out, beta_out, .false.)
            call parse_functional(xcfun_fde_alda, fde_xcf_alda, hfx_out, mu_out, beta_out, .false.)

!            call set_kin_fun_alda(fde_kef_alda)
!            call set_xc_fun_alda(fde_xcf_alda)             

            do i = 1, 20
               fde_set_functionals(i) = .true.
            enddo

            call xcfun_set_functional(fde_xcf     , order, parallel_fde) 
            call xcfun_set_functional(fde_xcf_alda, order, parallel_fde) 
            call xcfun_set_functional(fde_kef,      order, parallel_fde)
            call xcfun_set_functional(fde_kef_alda, order, parallel_fde)
            
         end subroutine fde_initialize_nadd_functionals


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_kef(id)
!   ----------------------------------------------------------------------------
            character(len=80) :: id

            kefun_fde = id
            call parse_functional(kefun_fde, fde_kef,      hfx_out, mu_out, beta_out, .false.)
         end subroutine fde_set_nadd_kef


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_xcf(id)
!   ----------------------------------------------------------------------------
            character(len=80) :: id

            xcfun_fde = id
            write (*,*) 'debug, nadd xc is ',id
            call parse_functional(xcfun_fde, fde_xcf,      hfx_out, mu_out, beta_out, .false.)
         end subroutine fde_set_nadd_xcf


!   ----------------------------------------------------------------------------
         subroutine fde_set_nadd_all_alda
!   ----------------------------------------------------------------------------
            is_fdersp_alda = .true.

!             call set_kin_fun_alda(fde_kef)
!             call set_xc_fun_alda(fde_xcf)
!            call parse_functional(kefun_fde_alda, fde_kef,      hfx_out, mu_out, beta_out, .false.)
!            call parse_functional(xcfun_fde_alda, fde_xcf,      hfx_out, mu_out, beta_out, .false.)
         end subroutine fde_set_nadd_all_alda 


!   ----------------------------------------------------------------------------
         subroutine fde_print_nadd_functionals(unit)
!   ----------------------------------------------------------------------------
            integer :: unit 

            write(unit,'(A,A20)') ' * Non-additive KE functional : ',kefun_fde
            write(unit,'(A,A20)') ' * Non-additive XC functional : ',xcfun_fde
!            write(unit,'(A,A20)') '   if .ALDA set, wil use : (XC) ',xcfun_fde_alda
!            write(unit,'(A,A20)') '                           (KE) ',kefun_fde_alda
         end subroutine fde_print_nadd_functionals


!   ----------------------------------------------------------------------------
         subroutine get_fde_nadd_fun_derv(order,        &
                                          block_length, &
                                          w,            &
                                          r_0,          &
                                          r_0_t,        &
                                          z_0,          &
                                          z_0_t,        &
                                          n,            &
                                          derv,         &
                                          parallel_fde)
!   ----------------------------------------------------------------------------
         integer,          intent(in)  :: order
         integer,          intent(in)  :: block_length
         real(8),          intent(in)  :: w(*)
         real(8),          intent(in)  :: r_0(*)
         real(8),          intent(in)  :: z_0(*)
         real(8),          intent(in)  :: r_0_t(*)
         real(8),          intent(in)  :: z_0_t(*)
         integer,          intent(in)  :: n
         real(8),          intent(out) :: derv(max_block_length, 0:n)
         logical,          intent(in)  :: parallel_fde
         real(8)                       :: dervtmp(max_block_length, 0:n)
         integer :: i, j

         derv = 0.0d0

         if (fde_set_functionals(order)) then
            call xcfun_set_functional(fde_xcf     , order, parallel_fde) 
            call xcfun_set_functional(fde_xcf_alda, order, parallel_fde) 
            call xcfun_set_functional(fde_kef,      order, parallel_fde)
            call xcfun_set_functional(fde_kef_alda, order, parallel_fde)
            fde_set_functionals(order) = .false.
         endif

         if (.not.skip_xc_contribution) then

         call get_functional_derv(fde_xcf,       &
                               fde_xcf_alda,  &
                               fde_xcf_alda,  &
                               order,         &
                               block_length,  &
                               w,             &
                               r_0_t,         &
                               z_0_t,         &
                               n,             &
                               dervtmp,       &
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

         do i = 1, max_block_length
            do j = 0, n
               derv(i,j) = derv(i,j) + dervtmp(i,j)
            enddo
         enddo

         call get_functional_derv(fde_xcf,       &
                               fde_xcf_alda,  &
                               fde_xcf_alda,  &
                               order,         &
                               block_length,  &
                               w,             &
                               r_0,           &
                               z_0,           &
                               n,             &
                               dervtmp,       &
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

         do i = 1, max_block_length
            do j = 0, n
               derv(i,j) = derv(i,j) - dervtmp(i,j)
            enddo
         enddo

      endif

! nw we construct the contribution from the kinetic energy to the response
! o the embedding potential.
! frst set functional to the of-kef, do the non-additive ke contribution

      if (.not.skip_kin_contribution) then

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
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

            do i = 1, max_block_length
               do j = 0, n
                  derv(i,j) = derv(i,j) + dervtmp(i,j)
               enddo
            enddo

! nw subtract the active subsystem's kinetic energy contribution

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
                               alda_real=is_fdersp_alda, &
                               alda_imag=is_fdersp_alda)

            do i = 1, max_block_length
               do j = 0, n
                  derv(i,j) = derv(i,j) - dervtmp(i,j)
               enddo
            enddo

         endif

      end subroutine get_fde_nadd_fun_derv

end module fde_nadd_derv

