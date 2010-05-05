module xc_driver_helpers

  use dft_cfg
  use ao_eval
  use dft_functional_derv
  use density_eval

contains

  subroutine lr(block_length,       &
                mat_dim,            &
                nz,                 &
                nr_mat,             &
                dmaterturbed,     &
                fmat,               &
                isym,               &
                ih,                 &
                allow_nonhermitian, &
                ao_p,               &
                r_0,                &
                gr_0,               &
                z_0,                &
                w,                  &
                is_gga_qr,          &
                is_gga_qi,          &
                buffer,             &
                n,                  &
                derv)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmaterturbed(*)
    real(8), intent(inout) :: fmat(*)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
    logical, intent(in)    :: allow_nonhermitian
    real(8), intent(in)    :: ao_p(*)
    real(8), intent(in)    :: r_0(*)
    real(8), intent(in)    :: gr_0(*)
    real(8), intent(in)    :: z_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(mxbllen, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi

    real(8)                :: u_qr,    u_qi(3)
    real(8)                :: v_qr(3), v_qi(3, 3)

    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    real(8)                :: ddot
    integer                :: k_start, k, im, im_off, j
!   ----------------------------------------------------------------------------

    k_start = 1
    if (dft_cfg_sdft_collinear) k_start = 3 !collinear approximation

    do im = 1, nr_mat

      calculate_qr = .false.
      calculate_qi = .false.

      if (ih(im) == +1) calculate_qr = .true.
      if (ih(im) == -1) calculate_qi = .true.

      if (ih(im) == 0 .and. allow_nonhermitian) then
        calculate_qr = .true.
        calculate_qi = .true.
      end if

      if (dft_cfg_no_sdft) then
        calculate_qi = .false.
      end if

      im_off = mat_dim*mat_dim*nz*(im - 1)

!     hermitian contribution (qr = quaternion real)
      if (calculate_qr) then

         if (is_gga_qr) then
            if (allow_nonhermitian) then
               call get_gn_nonhermitian(r_b,                      &
                                        gr_b,                     &
                                        isym(im) - 1,             &
                                        mat_dim,                  &
                                        dmaterturbed(1 + im_off), &
                                        buffer,                   &
                                        ao_p,                     &
                                        ao_p(ao_start(1, 0)))
            else
               call get_gn(r_b,                      &
                           gr_b,                     &
                           isym(im) - 1,             &
                           mat_dim,                  &
                           dmaterturbed(1 + im_off), &
                           buffer,                   &
                           ao_p,                     &
                           ao_p(ao_start(1, 0)))
            end if
         else
            call get_n(r_b,                        &
                       isym(im) - 1,               &
                       mat_dim,                    &
                       dmaterturbed(1 + im_off), &
                       buffer,                     &
                       ao_p)
         end if

        u_qr = derv(1, d_rr)*r_b
        
        if (is_gga_qr) then
          z_b  = 2.0d0*ddot(3, gr_0, 1, gr_b, 1)
          u_qr = u_qr + derv(1, d_rz)*z_b

          do j = 1, 3
            v_qr(j) = gr_0(j)*(                   & 
                                r_b*derv(1, d_rz) &
                              + z_b*derv(1, d_zz) &   
                              )                   &
                    + gr_b(j)*derv(1, d_z)
          end do
          v_qr = 2.0d0*v_qr
        end if

        if (is_gga_qr) then
          call gga_distribute_r(u_qr,             &
                                v_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p,             &
                                ao_p(ao_start(1, 0)))
        else
          call lda_distribute_r(u_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)
        end if

      end if !end of qr contribution

!     antihermitian contribution (qi = quaternion imaginary)
      if (calculate_qi) then

        
        if (is_gga_qi) then
          if (allow_nonhermitian) then
            call get_gs_nonhermitian(s_b,                       &
                                     gs_b,                       &
                                     isym(im) - 1,               &
                                     mat_dim,                    &
                                     dmaterturbed(1 + im_off), &
                                     buffer,                     &
                                     ao_p,                       &
                                     ao_p(ao_start(1, 0)))
          else
            call get_gs(s_b,                       &
                        gs_b,                       &
                        isym(im) - 1,               &
                        mat_dim,                    &
                        dmaterturbed(1 + im_off), &
                        buffer,                     &
                        ao_p,                       &
                        ao_p(ao_start(1, 0)))
          end if
        else
        call get_s(s_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        end if

!       initialize otherwise there can be numerical dirt
!       and then collinear numbers are wrong
        u_qi = 0.0d0
        v_qi = 0.0d0

        do k = k_start, 3

          u_qi(k) = derv(1, d_ss)*s_b(k)
         
          if (is_gga_qi) then
            y_b(k)  = ddot(3, gr_0, 1, gs_b(1, k), 1)
            u_qi(k) = u_qi(k) + derv(1, d_sy)*y_b(k)

            do j = 1, 3
              v_qi(j, k) = gr_0(j)*(                   & 
                                  s_b(k)*derv(1, d_sy) &
                                + y_b(k)*derv(1, d_yy) &   
                                )                      &
                      + 2.0d0*gs_b(j, k)*derv(1, d_x)
            end do
          end if
        end do

        if (is_gga_qi) then
          call gga_distribute_i(u_qi,             &
                                v_qi,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p,             &
                                ao_p(ao_start(1, 0)))
        else
          call lda_distribute_i(u_qi,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)
        end if

      end if !end of qi contribution

    end do

  end subroutine

  subroutine g_lao_contrib(block_length,       &
                mat_dim,            &
                nz,             &
                nr_mat,             &
                dmaterturbed,     &
                fmat,               &
                isym,               &
                ih,                 &
                allow_nonhermitian, &
                ao_p,               &
                r_0,                &
                gr_0,               &
                z_0,                &
                w,                  &
                is_gga_qr,          &
                is_gga_qi,          &
                buffer,             &
                n,                  &
                derv)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: block_length
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(in)    :: dmaterturbed(*)
    real(8), intent(inout) :: fmat(*)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
    logical, intent(in)    :: allow_nonhermitian
    real(8), intent(in)    :: ao_p(*)
    real(8), intent(in)    :: r_0(*)
    real(8), intent(in)    :: gr_0(*)
    real(8), intent(in)    :: z_0(*)
    real(8), intent(in)    :: w(*)
    logical, intent(in)    :: is_gga_qr
    logical, intent(in)    :: is_gga_qi
    real(8), intent(in)    :: buffer(*)
    integer, intent(in)    :: n
    real(8), intent(in)    :: derv(mxbllen, 0:n)
!   ----------------------------------------------------------------------------
    logical                :: calculate_qr, calculate_qi

    real(8)                :: u_qr,    u_qi(3)
    real(8)                :: v_qr(3), v_qi(3, 3)

    real(8)                ::  r_b,     s_b(3)
    real(8)                :: gr_b(3), gs_b(3, 3)

    real(8)                :: z_b
    real(8)                :: y_b(3)

    real(8)                :: ddot
    integer                :: k_start, k, im, im_off, j, off_mag
!   ----------------------------------------------------------------------------

    do im = 1, nr_mat

      im_off = mat_dim*mat_dim*nz*(im - 1)

        call get_n(r_b,                        &
                   isym(im) - 1,               &
                   mat_dim,                    &
                   dmaterturbed(1 + im_off), &
                   buffer,                     &
                   ao_p)
        
        u_qr = derv(1, d_rr)*r_b
        
          call lda_distribute_r(u_qr,             &
                                isym(im) - 1,     &
                                mat_dim,          &
                                fmat(1 + im_off), &
                                ao_p)

    end do

  end subroutine

end module
