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

module density_eval

   use interface_ao

   implicit none

   public mat_lao_sym
   public gga_sym
   public get_n
   public get_gamma5
   public get_s
   public get_j
   public get_t
   public get_tau
   public get_tau_imag
   public get_gn
   public get_gn_nonhermitian
   public get_gj
   public get_gs
   public get_gs_nonhermitian
   public omega_real
   public omega_imag
   public nabla_omega_real
   public nabla_omega_imag
   public lda_london_distribute_r
   public gga_london_distribute_r
   public get_s_lao
   public get_gs_lao

   private

  integer :: jxyz_start = 1 !noncollinear
! integer :: jxyz_start = 3 !collinear

  save

contains

  subroutine mat_lao_sym(mat_dim, &
                         nz,      &
                         nr_mat,  &
                         mat,     &
                         isym)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(inout) :: mat(mat_dim, mat_dim, nz, nr_mat)
    integer, intent(in)    :: isym(nr_mat)
!   ----------------------------------------------------------------------------
    integer                :: irep, iblock, imat, iz
    integer                :: i, j, ii, jj, h(4)
    integer                :: st1, st2, nr1, nr2
    real(8)                :: average
!   ----------------------------------------------------------------------------

    do imat = 1, nr_mat

       irep = isym(imat) - 1

       h(1) = - 1
       h(2:nz) = 1

       do iz = 1, nz
          do iblock = 1, llss_nr_density_blocks(0, irep)
             nr1 = llss_nr_ao(1, iblock, 0, irep)
             st1 = llss_st_ao(1, iblock, 0, irep)
             nr2 = llss_nr_ao(2, iblock, 0, irep)
             st2 = llss_st_ao(2, iblock, 0, irep)

             do i = 1, nr1
                do j = 1, nr2
                   ii = st1 + i - 1
                   jj = st2 + j - 1

                   average = mat(ii, jj, iz, imat) &
                           + h(iz)*mat(jj, ii, iz, imat)
                   average = 0.5d0*average

                   mat(ii, jj, iz, imat) = average
                   mat(jj, ii, iz, imat) = h(iz)*average
                end do
             end do
          end do
       end do
    end do

  end subroutine

  subroutine gga_sym(mat_dim, &
                     nz,      &
                     nr_mat,  &
                     mat,     &
                     isym,    &
                     ih)

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    integer, intent(in)    :: nr_mat
    real(8), intent(inout) :: mat(mat_dim, mat_dim, nz, nr_mat)
    integer, intent(in)    :: isym(nr_mat)
    integer, intent(in)    :: ih(nr_mat)
!   ----------------------------------------------------------------------------
    integer                :: irep, iblock, imat, jxyz
    integer                :: i, j, ii, jj
    integer                :: st1, st2, nr1, nr2
    logical                :: skip(0:3)
    real(8)                :: average
!   ----------------------------------------------------------------------------

    do imat = 1, nr_mat

       skip = .false.
       if (ih(imat) ==  1) skip(1:3) = .true.
       if (ih(imat) == -1) skip(0:0) = .true.

       irep = isym(imat) - 1

       do jxyz = 0, 3
          if (.not. skip(jxyz)) then
             do iblock = 1, llss_nr_density_blocks(jxyz, irep)
                nr1 = llss_nr_ao(1, iblock, jxyz, irep)
                st1 = llss_st_ao(1, iblock, jxyz, irep)
                nr2 = llss_nr_ao(2, iblock, jxyz, irep)
                st2 = llss_st_ao(2, iblock, jxyz, irep)

                if (.not. (st2 > st1)) then
                   do i = 1, nr1
                      do j = 1, nr2
                         ii = st1 + i - 1
                         jj = st2 + j - 1

                         average = mat(ii, jj, iq(jxyz), imat) &
                                 + mat(jj, ii, iq(jxyz), imat)
                         average = 0.5d0*average

                         mat(ii, jj, iq(jxyz), imat) = average
                         mat(jj, ii, iq(jxyz), imat) = average
                      end do
                   end do
                end if
             end do
          end if
       end do
    end do

  end subroutine

   subroutine get_n(n,       &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: n
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock
      real(8), external    :: ddot
!     --------------------------------------------------------------------------
   
      n = 0.0d0
      do iblock = 1, llss_nr_density_blocks(0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         nr2 = llss_nr_ao(2, iblock, 0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
    
         call dgemv('N',                   &
                    nr1,                   &
                    nr2,                   &
                    1.0d0,                 &
                    dmat(st1, st2, iq(0)), &
                    mat_dim,               &
                    ao(st2),               &
                    1,                     &
                    0.0d0,                 &
                    buffer,                &
                    1)
    
         n = n + ddot(nr1, ao(st1), 1, buffer, 1)
      end do
   
   end subroutine

   subroutine get_gamma5(n,       &
                         irep,    &
                         mat_dim, &
                         dmat,    &
                         buffer,  &
                         ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: n
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock
      real(8), external    :: ddot
!     --------------------------------------------------------------------------
   
      n = 0.0d0
      do iblock = 1, lssl_nr_density_blocks(0, irep)
         nr1 = lssl_nr_ao(1, iblock, 0, irep)
         st1 = lssl_st_ao(1, iblock, 0, irep)
         nr2 = lssl_nr_ao(2, iblock, 0, irep)
         st2 = lssl_st_ao(2, iblock, 0, irep)
    
         call dgemv('N',                   &
                    nr1,                   &
                    nr2,                   &
                    1.0d0,                 &
                    dmat(st1, st2, iq(0)), &
                    mat_dim,               &
                    ao(st2),               &
                    1,                     &
                    0.0d0,                 &
                    buffer,                &
                    1)
    
         n = n + ddot(nr1, ao(st1), 1, buffer, 1)
      end do
   
   end subroutine

   subroutine get_s(s,       &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: s(3)
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------
   
      s = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, llss_nr_density_blocks(jxyz, irep)
            nr1 = llss_nr_ao(1, iblock, jxyz, irep)
            st1 = llss_st_ao(1, iblock, jxyz, irep)
            nr2 = llss_nr_ao(2, iblock, jxyz, irep)
            st2 = llss_st_ao(2, iblock, jxyz, irep)
     
            call dgemv('N',                      &
                       nr1,                      &
                       nr2,                      &
                       1.0d0,                    &
                       dmat(st1, st2, iq(jxyz)), &
                       mat_dim,                  &
                       ao(st2),                  &
                       1,                        &
                       0.0d0,                    &
                       buffer,                   &
                       1)
     
            s(jxyz) = s(jxyz) &
                    + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao(st1), 1, buffer, 1)
         end do
      end do
   
   end subroutine

   subroutine get_j(j,       &
                    cval,    &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: j(3)
      real(8), intent(in)  :: cval
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------
   
      j = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, lssl_nr_density_blocks(jxyz, irep)
            nr1 = lssl_nr_ao(1, iblock, jxyz, irep)
            st1 = lssl_st_ao(1, iblock, jxyz, irep)
            nr2 = lssl_nr_ao(2, iblock, jxyz, irep)
            st2 = lssl_st_ao(2, iblock, jxyz, irep)
     
            call dgemv('N',                      &
                       nr1,                      &
                       nr2,                      &
                       1.0d0,                    &
                       dmat(st1, st2, iq(jxyz)), &
                       mat_dim,                  &
                       ao(st2),                  &
                       1,                        &
                       0.0d0,                    &
                       buffer,                   &
                       1)
     
            j(jxyz) = j(jxyz) + ddot(nr1, ao(st1), 1, buffer, 1)
         end do
      end do
      j = -cval*j
   
   end subroutine


  subroutine get_gn(n,       &
                    gn,      &
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao,      &
                    ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    n  = 0.0d0
    gn = 0.0d0

    do iblock = 1, llss_nr_density_blocks(0, irep)
       nr1 = llss_nr_ao(1, iblock, 0, irep)
       st1 = llss_st_ao(1, iblock, 0, irep)
       nr2 = llss_nr_ao(2, iblock, 0, irep)
       st2 = llss_st_ao(2, iblock, 0, irep)

       call dgemv('N',                      &
                  nr1,                      &
                  nr2,                      &
                  1.0d0,                    &
                  dmat(st1, st2, iq(0)), &
                  mat_dim,                  &
                  ao(st2),                  &
                  1,                        &
                  0.0d0,                    &
                  buffer,                   &
                  1)

       n = n + ddot(nr1, ao(st1), 1, buffer, 1)
       do ixyz = 1, 3
          gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
       end do
    end do

    gn = 2.0d0*gn

  end subroutine

   subroutine get_tau(t,       &
                      irep,    &
                      mat_dim, &
                      dmat,    &
                      buffer,  &
                      ao_g1)

!     --------------------------------------------------------------------------
!     radovan bast 2010-04-15
!                  this routine calculates the \tau = \nabla_i \phi_k \nabla_i \phi_k
!     
!                  needs ao_g1
!                  integral approaches the kinetic energy as c -> oo
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao_g1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      t = 0.0d0
      
      do iblock = 1, llss_nr_density_blocks(0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         nr2 = llss_nr_ao(2, iblock, 0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
      
         do ixyz = 1, 3
            call dgemv('N',                   &
                       nr1,                   &
                       nr2,                   &
                       1.0d0,                 &
                       dmat(st1, st2, iq(0)), &
                       mat_dim,               &
                       ao_g1(st2, ixyz),      &
                       1,                     &
                       0.0d0,                 &
                       buffer,                &
                       1)
      
            t = t + ddot(nr1, ao_g1(st1, ixyz), 1, buffer, 1)
         end do
      end do
      
      t = 0.5d0*t

   end subroutine

   subroutine get_tau_imag(t,       &
                           irep,    &
                           mat_dim, &
                           dmat,    &
                           buffer,  &
                           ao_g1)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: t(3)
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao_g1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, ixyz, jxyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      t = 0.0d0
     
      do jxyz = 1, 3 
         do iblock = 1, llss_nr_density_blocks(jxyz, irep)
            nr1 = llss_nr_ao(1, iblock, jxyz, irep)
            st1 = llss_st_ao(1, iblock, jxyz, irep)
            nr2 = llss_nr_ao(2, iblock, jxyz, irep)
            st2 = llss_st_ao(2, iblock, jxyz, irep)
         
            do ixyz = 1, 3
               call dgemv('N',                      &
                          nr1,                      &
                          nr2,                      &
                          1.0d0,                    &
                          dmat(st1, st2, iq(jxyz)), &
                          mat_dim,                  &
                          ao_g1(st2, ixyz),         &
                          1,                        &
                          0.0d0,                    &
                          buffer,                   &
                          1)
         
               t(jxyz) = t(jxyz) &
             + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao_g1(st1, ixyz), 1, buffer, 1)
            end do
         end do
      end do
      
      t = 0.5d0*t

   end subroutine

  subroutine get_t(t,       &
                   cval,    &
                   irep,    &
                   mat_dim, &
                   dmat,    &
                   buffer,  &
                   ao_g0,   &
                   ao_g1)

!   ----------------------------------------------------------------------------
!   radovan bast 2010-04-15
!                this routine calculates the kinetic energy density
!                c \psi+(r) \alpha \cdot p \psi(r)
!                you can call it as you would call for instance get_gn
!
!   caveat: (1)  it needs ao_g1 and they may not be
!                there in the naive lda case
!           (2)  it needs small components and by default
!                small components are not calculated for levy-leblond
!                this has to be reactivated in the ao evaluation
!                nonrel and x2c(4) fail for this reason
!                (will add nonrelativistic get_t soon)
!           (2)  with the default dft grid it can be very difficult to
!                integrate it up to the analytic kinetic energy
!                (with atoms i got it 99.8% right)
!   ----------------------------------------------------------------------------
    real(8), intent(out) :: t
    real(8), intent(in)  :: cval
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao_g0(mat_dim)
    real(8), intent(in)  :: ao_g1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    t = 0.0d0

    do jxyz = 1, 3
      do iblock = 1, lssl_nr_density_blocks(jxyz, irep)

        nr1 = lssl_nr_ao(1, iblock, jxyz, irep)
        st1 = lssl_st_ao(1, iblock, jxyz, irep)
        nr2 = lssl_nr_ao(2, iblock, jxyz, irep)
        st2 = lssl_st_ao(2, iblock, jxyz, irep)

        call dgemv('N',                           &
                   nr1,                           &
                   nr2,                           &
                   1.0d0,                         &
                   dmat(st1, st2, iq(jxyz)), &
                   mat_dim,                       &
                   ao_g0(st2),                    &
                   1,                             &
                   0.0d0,                         &
                   buffer,                        &
                   1)

        t = t + ddot(nr1, ao_g1(st1, jxyz), 1, buffer, 1)

      end do
    end do

    t = t*cval

  end subroutine

  subroutine get_gn_nonhermitian(n,       &
                                 gn,      &
                                 irep,    &
                                 mat_dim, &
                                 dmat,    &
                                 buffer,  &
                                 ao,      &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    n  = 0.0d0
    gn = 0.0d0

    do iblock = 1, llss_nr_density_blocks(0, irep)
       nr1 = llss_nr_ao(1, iblock, 0, irep)
       st1 = llss_st_ao(1, iblock, 0, irep)
       nr2 = llss_nr_ao(2, iblock, 0, irep)
       st2 = llss_st_ao(2, iblock, 0, irep)

       call dgemv('N',                      &
                  nr1,                      &
                  nr2,                      &
                  1.0d0,                    &
                  dmat(st1, st2, iq(0)), &
                  mat_dim,                  &
                  ao(st2),                  &
                  1,                        &
                  0.0d0,                    &
                  buffer,                   &
                  1)

       n = n + ddot(nr1, ao(st1), 1, buffer, 1)
       do ixyz = 1, 3
         gn(ixyz) = gn(ixyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
       end do

       do ixyz = 1, 3
         call dgemv('N',                      &
                    nr1,                      &
                    nr2,                      &
                    1.0d0,                    &
                    dmat(st1, st2, iq(0)), &
                    mat_dim,                  &
                    ao1(st2, ixyz),           &
                    1,                        &
                    0.0d0,                    &
                    buffer,                   &
                    1)

         gn(ixyz) = gn(ixyz) + ddot(nr1, ao(st1), 1, buffer, 1)
       end do
    end do

  end subroutine

   subroutine get_gj(v,       &
                     gv,      &
                     cval,    &
                     irep,    &
                     mat_dim, &
                     dmat,    &
                     buffer,  &
                     ao,      &
                     ao1)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: v(3)
      real(8), intent(out) :: gv(3, 3)
      real(8), intent(in)  :: cval
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: ao1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      v  = 0.0d0
      gv = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, lssl_nr_density_blocks(jxyz, irep)
            nr1 = lssl_nr_ao(1, iblock, jxyz, irep)
            st1 = lssl_st_ao(1, iblock, jxyz, irep)
            nr2 = lssl_nr_ao(2, iblock, jxyz, irep)
            st2 = lssl_st_ao(2, iblock, jxyz, irep)
  
            call dgemv('N',                      &
                       nr1,                      &
                       nr2,                      &
                       1.0d0,                    &
                       dmat(st1, st2, iq(jxyz)), &
                       mat_dim,                  &
                       ao(st2),                  &
                       1,                        &
                       0.0d0,                    &
                       buffer,                   &
                       1)
  
            v(jxyz) = v(jxyz) + ddot(nr1, ao(st1), 1, buffer, 1)
            do ixyz = 1, 3
               gv(ixyz, jxyz) = gv(ixyz, jxyz) + ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
            end do
         end do
      end do
      gv = -cval*2.0d0*gv

   end subroutine

   subroutine get_gs(v,       &
                     gv,      &
                     irep,    &
                     mat_dim, &
                     dmat,    &
                     buffer,  &
                     ao,      &
                     ao1)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: v(3)
      real(8), intent(out) :: gv(3, 3)
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(mat_dim)
      real(8), intent(in)  :: ao1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      v  = 0.0d0
      gv = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, llss_nr_density_blocks(jxyz, irep)
            nr1 = llss_nr_ao(1, iblock, jxyz, irep)
            st1 = llss_st_ao(1, iblock, jxyz, irep)
            nr2 = llss_nr_ao(2, iblock, jxyz, irep)
            st2 = llss_st_ao(2, iblock, jxyz, irep)
  
            call dgemv('N',                      &
                       nr1,                      &
                       nr2,                      &
                       1.0d0,                    &
                       dmat(st1, st2, iq(jxyz)), &
                       mat_dim,                  &
                       ao(st2),                  &
                       1,                        &
                       0.0d0,                    &
                       buffer,                   &
                       1)
  
            v(jxyz) = v(jxyz) &
                    + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao(st1), 1, buffer, 1)
            do ixyz = 1, 3
               gv(ixyz, jxyz) = gv(ixyz, jxyz) &
             + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
            end do
         end do
      end do
      gv = 2.0d0*gv

   end subroutine

  subroutine get_gs_nonhermitian(s,       &
                                 gs,      &
                                 irep,    &
                                 mat_dim, &
                                 dmat,    &
                                 buffer,  &
                                 ao,      &
                                 ao1)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3)
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao1(mat_dim, 3)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    s = 0.0d0
    gs = 0.0d0

    do jxyz = jxyz_start, 3

       do iblock = 1, llss_nr_density_blocks(jxyz, irep)
          nr1 = llss_nr_ao(1, iblock, jxyz, irep)
          st1 = llss_st_ao(1, iblock, jxyz, irep)
          nr2 = llss_nr_ao(2, iblock, jxyz, irep)
          st2 = llss_st_ao(2, iblock, jxyz, irep)

          call dgemv('N',                         &
                     nr1,                         &
                     nr2,                         &
                     1.0d0,                       &
                     dmat(st1, st2, iq(jxyz)), &
                     mat_dim,                     &
                     ao(st2),                     &
                     1,                           &
                     0.0d0,                       &
                     buffer,                      &
                     1)

          s(jxyz) = s(jxyz) + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao(st1), 1, buffer, 1)
          do ixyz = 1, 3
             gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                            + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao1(st1, ixyz), 1, buffer, 1)
          end do

          do ixyz = 1, 3
             call dgemv('N',                         &
                        nr1,                         &
                        nr2,                         &
                        1.0d0,                       &
                        dmat(st1, st2, iq(jxyz)), &
                        mat_dim,                     &
                        ao1(st2, ixyz),              &
                        1,                           &
                        0.0d0,                       &
                        buffer,                      &
                        1)

             gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                            + llss_prefactor(iblock, jxyz, irep)*ddot(nr1, ao(st1), 1, buffer, 1)
          end do
       end do

    end do

  end subroutine

   subroutine omega_real(fmat,    &
                         im,      &
                         irep,    &
                         mat_dim, &
                         nz,      &
                         u,       &
                         ao_g0)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u
      real(8), intent(in)    :: ao_g0(mat_dim)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do iblock = 1, llss_nr_density_blocks(0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         nr2 = llss_nr_ao(2, iblock, 0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
      
         do i = 1, nr2
            ii = st2 + i - 1
            f  = u*ao_g0(ii)
      
            call daxpy(nr1,                      &
                       f,                        &
                       ao_g0(st1),               &
                       1,                        &
                       fmat(st1, ii, iq(0), im), &
                       1)
         end do
      end do

   end subroutine

   subroutine nabla_omega_real(fmat,    &
                               im,      &
                               irep,    &
                               mat_dim, &
                               nz,      &
                               u,       &
                               v,       &
                               ao_g0,   &
                               ao_g1)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u
      real(8), intent(in)    :: v(3)
      real(8), intent(in)    :: ao_g0(mat_dim)
      real(8), intent(in)    :: ao_g1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do iblock = 1, llss_nr_density_blocks(0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         nr2 = llss_nr_ao(2, iblock, 0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
      
         do i = 1, nr2
            ii = st2 + i - 1
            f =       u   *ao_g0(ii)    &
              + 2.0d0*v(1)*ao_g1(ii, 1) &
              + 2.0d0*v(2)*ao_g1(ii, 2) &
              + 2.0d0*v(3)*ao_g1(ii, 3)
      
            call daxpy(nr1,                      &
                       f,                        &
                       ao_g0(st1),               &
                       1,                        &
                       fmat(st1, ii, iq(0), im), &
                       1)
         end do
      end do

   end subroutine

   subroutine omega_imag(fmat,    &
                         im,      &
                         irep,    &
                         mat_dim, &
                         nz,      &
                         u,       &
                         ao_g0)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u(3)
      real(8), intent(in)    :: ao_g0(mat_dim)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do jxyz = jxyz_start, 3

         do iblock = 1, llss_nr_density_blocks(jxyz, irep)
            nr1 = llss_nr_ao(1, iblock, jxyz, irep)
            st1 = llss_st_ao(1, iblock, jxyz, irep)
            nr2 = llss_nr_ao(2, iblock, jxyz, irep)
            st2 = llss_st_ao(2, iblock, jxyz, irep)
         
            do i = 1, nr2
               ii = st2 + i - 1
               f  = u(jxyz)*ao_g0(ii)*llss_prefactor(iblock, jxyz, irep)
         
               call daxpy(nr1,                         &
                          f,                           &
                          ao_g0(st1),                  &
                          1,                           &
                          fmat(st1, ii, iq(jxyz), im), &
                          1)
            end do
         end do
      end do

   end subroutine

   subroutine nabla_omega_imag(fmat,    &
                               im,      &
                               irep,    &
                               mat_dim, &
                               nz,      &
                               u,       &
                               v,       &
                               ao_g0,   &
                               ao_g1)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u(3)
      real(8), intent(in)    :: v(3, 3)
      real(8), intent(in)    :: ao_g0(mat_dim)
      real(8), intent(in)    :: ao_g1(mat_dim, 3)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do jxyz = jxyz_start, 3

         do iblock = 1, llss_nr_density_blocks(jxyz, irep)
            nr1 = llss_nr_ao(1, iblock, jxyz, irep)
            st1 = llss_st_ao(1, iblock, jxyz, irep)
            nr2 = llss_nr_ao(2, iblock, jxyz, irep)
            st2 = llss_st_ao(2, iblock, jxyz, irep)
         
            do i = 1, nr2
               ii = st2 + i - 1
               f =       u(   jxyz)*ao_g0(ii)    &
                 + 2.0d0*v(1, jxyz)*ao_g1(ii, 1) &
                 + 2.0d0*v(2, jxyz)*ao_g1(ii, 2) &
                 + 2.0d0*v(3, jxyz)*ao_g1(ii, 3)
              
               f = llss_prefactor(iblock, jxyz, irep)*f
         
               call daxpy(nr1,                         &
                          f,                           &
                          ao_g0(st1),                  &
                          1,                           &
                          fmat(st1, ii, iq(jxyz), im), &
                          1)
            end do
         end do
      end do

   end subroutine


   subroutine lda_london_distribute_r(s, isym, mat_dim, nz, fmat, ao, ao_b1, cx, cy, cz)
   
!    ----------------------------------------------------------------------------
     real(8), intent(in)    :: s
     integer, intent(in)    :: isym(*)
     integer, intent(in)    :: mat_dim
     integer, intent(in)    :: nz
     real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
     real(8), intent(in)    :: ao(mat_dim)
     real(8), intent(in)    :: ao_b1(mat_dim, *)
     real(8), intent(in)    :: cx, cy, cz
!    ----------------------------------------------------------------------------
     integer                :: nr1, nr2, st1, st2, iblock, ii, i, a, b, c
     real(8)                :: f
     real(8)                :: r_co(3)
!    ----------------------------------------------------------------------------
   
     r_co(1) = cx
     r_co(2) = cy
     r_co(3) = cz

!gosia: todo: "do a = 1, nr_xyz_comp" instead (now we add x/y/z on top)
     do a = 1, 3
       select case (a)
         case (1)
           b = 2; c = 3
         case (2)
           b = 3; c = 1
         case (3)
           b = 1; c = 2
       end select
     
         do iblock = 1, llss_nr_density_blocks(0, isym(a) - 1)
     
           nr1 = llss_nr_ao(1, iblock, 0, isym(a) - 1)
           st1 = llss_st_ao(1, iblock, 0, isym(a) - 1)
           nr2 = llss_nr_ao(2, iblock, 0, isym(a) - 1)
           st2 = llss_st_ao(2, iblock, 0, isym(a) - 1)
     
           do i = 1, nr2
             ii = st2 + i - 1
             f  = 2.0d0*(r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*s
     
             call daxpy(nr1,                     &
                        f,                       &
                        ao(st1),                 &
                        1,                       &
                        fmat(st1, ii, iq(0), a), &
                        1)
     
           end do
     
         end do
     end do

   end subroutine

   subroutine gga_london_distribute_r(d_r,      &
                                      d_z,      &
                                      n_g1,     &
                                      isym,     &
                                      mat_dim,  &
                                      nz,       &
                                      fmat,     &
                                      ao_g0_b0, &
                                      ao_g0_b1, &
                                      ao_g1_b0, &
                                      ao_g1_b1, &
                                      cx,       &
                                      cy,       &
                                      cz)
   
   
!    ----------------------------------------------------------------------------
     real(8), intent(in)    :: d_r
     real(8), intent(in)    :: d_z
     real(8), intent(in)    :: n_g1(3)
     integer, intent(in)    :: isym(*)
     integer, intent(in)    :: mat_dim
     integer, intent(in)    :: nz
     real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
     real(8), intent(in)    :: ao_g0_b0(mat_dim)
     real(8), intent(in)    :: ao_g0_b1(mat_dim, 3)
     real(8), intent(in)    :: ao_g1_b0(mat_dim, 3)
     real(8), intent(in)    :: ao_g1_b1(mat_dim, 3, 3)
     real(8), intent(in)    :: cx, cy, cz
!    ----------------------------------------------------------------------------
     integer                :: nr1, nr2, st1, st2, iblock, i, ii
     integer                :: a, b, c, ixyz
     real(8)                :: f, f_temp, f_vec(mat_dim)
     real(8)                :: r(3)
!    ----------------------------------------------------------------------------
   
     r(1) = cx
     r(2) = cy
     r(3) = cz
   
     do a = 1, 3
   
       select case (a)
         case (1)
           b = 2; c = 3
         case (2)
           b = 3; c = 1
         case (3)
           b = 1; c = 2
       end select
   
       do iblock = 1, llss_nr_density_blocks(0, isym(a) - 1)
   
         nr1 = llss_nr_ao(1, iblock, 0, isym(a) - 1)
         st1 = llss_st_ao(1, iblock, 0, isym(a) - 1)
         nr2 = llss_nr_ao(2, iblock, 0, isym(a) - 1)
         st2 = llss_st_ao(2, iblock, 0, isym(a) - 1)
   
         f_vec(:)  = r(b)*ao_g0_b1(:, c) - r(c)*ao_g0_b1(:, b)
   
         do i = 1, nr2
           ii = st2 + i - 1
   
           f_temp = 2.0d0*(r(c)*ao_g0_b1(ii, b) - r(b)*ao_g0_b1(ii, c))
           f = f_temp*d_r
   
!          so far it's lda type contribution, rest are gradient terms
   
           f_temp = 2.0d0*(n_g1(c)*ao_g0_b1(ii, b) - n_g1(b)*ao_g0_b1(ii, c))
           f = f + f_temp*d_z
   
           do ixyz = 1, 3
             f_temp = 2.0d0*(r(c)*ao_g1_b1(ii, b, ixyz) - r(b)*ao_g1_b1(ii, c, ixyz))
             f = f + f_temp*d_z*n_g1(ixyz)
           end do
   
           call daxpy(nr1,                          &
                      f,                            &
                      ao_g0_b0(st1),                &
                      1,                            &
                      fmat(st1, ii, iq(0), a),      &
                      1)
   
   
           f_temp = 2.0d0*(n_g1(1)*ao_g1_b0(ii, 1)  &
                  +        n_g1(2)*ao_g1_b0(ii, 2)  &
                  +        n_g1(3)*ao_g1_b0(ii, 3))
   
           f = f_temp*d_z
   
           call daxpy(nr1,                          &
                      f,                            &
                      f_vec(st1),                   &
                      1,                            &
                      fmat(st1, ii, iq(0), a),      &
                      1)
   
         end do ! i = 1, nr2
       end do   ! iblock = 1, llss_nr_density_blocks(0, isym(a) - 1)
     end do     ! a = 1, 3
   
   end subroutine


   subroutine get_s_lao(s, isym, mat_dim, dmat, ao, ao_b1, cx, cy, cz, buffer)
   
!    ----------------------------------------------------------------------------
     real(8), intent(out) :: s(3, 3)
     integer, intent(in)  :: isym(*)
     integer, intent(in)  :: mat_dim
     real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
     real(8), intent(in)  :: ao(mat_dim)
     real(8), intent(in)  :: ao_b1(mat_dim, *)
     real(8), intent(in)  :: cx, cy, cz
     real(8), intent(in)  :: buffer(*)
!    ----------------------------------------------------------------------------
     integer              :: nr1, nr2, st1, st2, iblock, jxyz
     integer              :: a, b, c, irp1, irp2, ii, jj
     real(8)              :: r_co(3), f(mat_dim), ff
     real(8), external    :: ddot
     integer :: i, j, k
!    ----------------------------------------------------------------------------
   
     r_co(1) = cx
     r_co(2) = cy
     r_co(3) = cz
   
     s = 0.0d0
     f = 0.0d0
     ff= 0.0d0

     do a = 1, 3   ! Bx, By, Bz
       select case (a)
         case (1)
           b = 2; c = 3
         case (2)
           b = 3; c = 1
         case (3)
           b = 1; c = 2
       end select
    
       do jxyz = jxyz_start, 3
    
         do iblock = 1, llss_nr_density_blocks(jxyz, isym(a) - 1)
    
           nr1 = llss_nr_ao(1, iblock, jxyz, isym(a) - 1)
           st1 = llss_st_ao(1, iblock, jxyz, isym(a) - 1)
           nr2 = llss_nr_ao(2, iblock, jxyz, isym(a) - 1)
           st2 = llss_st_ao(2, iblock, jxyz, isym(a) - 1)
    
           do ii = st1, st1 + nr1 - 1
             do jj = st2, st2 + nr2 - 1
    
               ff = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                  - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)


               s(jxyz, a) = s(jxyz, a) - dmat(jj, ii, iq(jxyz))*ff
    
             end do
           end do

         end do  ! do iblock = 1, llss_nr_density_blocks(jxyz, isym(a) - 1)

       end do    ! do jxyz = 1, jxyz_start, 3
     end do      ! do a = 1, 3

   end subroutine



   subroutine get_gs_lao(s,          &
                         gs,         &
                         isym,       &
                         mat_dim,    &
                         dmat,       &
                         buffer,     &
                         ao_g0_b0,   &
                         ao_g1_b0,   &
                         ao_g0_b1,   &
                         ao_g1_b1,   &
                         cx, cy, cz)
   
!    ----------------------------------------------------------------------------
     real(8), intent(out) :: s(3, 3)
     real(8), intent(out) :: gs(3, 3, 3)
     integer, intent(in)  :: isym(*)
     integer, intent(in)  :: mat_dim
     real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
     real(8), intent(in)  :: buffer(*)
     real(8), intent(in)  :: ao_g0_b0(mat_dim)
     real(8), intent(in)  :: ao_g1_b0(mat_dim, 3)
     real(8), intent(in)  :: ao_g0_b1(mat_dim, 3)
     real(8), intent(in)  :: ao_g1_b1(mat_dim, 3, 3)
     real(8), intent(in)  :: cx, cy, cz
!    ----------------------------------------------------------------------------
     integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
     integer              :: a, b, c, irp1, irp2, ii, jj
     real(8)              :: r_co(3), f(mat_dim), ff
!    ----------------------------------------------------------------------------

!gosia: fixme
     r_co(1) = cx
     r_co(2) = cy
     r_co(3) = cz

     gs = 0.0d0
     s  = 0.0d0
     f  = 0.0d0
     ff = 0.0d0

!gosia todo: do a = 1, nr_xyz_comp

     do a = 1, 3   ! Bx, By, Bz
       select case (a)
         case (1)
           b = 2; c = 3
         case (2)
           b = 3; c = 1
         case (3)
           b = 1; c = 2
       end select
    
       do jxyz = jxyz_start, 3
    
         do iblock = 1, llss_nr_density_blocks(jxyz, isym(a) - 1)
    
           nr1 = llss_nr_ao(1, iblock, jxyz, isym(a) - 1)
           st1 = llss_st_ao(1, iblock, jxyz, isym(a) - 1)
           nr2 = llss_nr_ao(2, iblock, jxyz, isym(a) - 1)
           st2 = llss_st_ao(2, iblock, jxyz, isym(a) - 1)


           do ii = st1, st1 + nr1 - 1
             do jj = st2, st2 + nr2 - 1
    
               ff = (r_co(c)*ao_g0_b1(jj, b) - r_co(b)*ao_g0_b1(jj, c))*ao_g0_b0(ii) &
                  - (r_co(c)*ao_g0_b1(ii, b) - r_co(b)*ao_g0_b1(ii, c))*ao_g0_b0(jj)


               s(jxyz, a) = s(jxyz, a) - dmat(jj, ii, iq(jxyz))*ff
    
             end do
           end do

         end do

       end do
     end do

   end subroutine

end module
