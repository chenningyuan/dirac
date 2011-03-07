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

   use basis_info

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
   public get_n_block
   public get_block_threshold

   private

  integer :: jxyz_start = 1 !noncollinear
! integer :: jxyz_start = 3 !collinear

  save

! cval, will not work in parallel!
#include "dcbgen.h"

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
                    irep,    &
                    mat_dim, &
                    dmat,    &
                    buffer,  &
                    ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: j(3)
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

   subroutine get_n_block(                 &
                          is_gga,          &
                          do_lr,           &
                          block_length,    &
                          block_threshold, &
                          mat_dim,         &
                          nz,              &
                          buffer,          &
                          ao_block,        &
                          nr_fmat,         &
                          fmat,            &
                          dmat_0,          &
                          nr_dmat,         &
                          dmat_p,          &
                          dmat_pg_sym,     &
                          dmat_ih_sym,     &
                          derv_length,     &
                          derv,            &
                          max_fun_derv,    &
                          w,               &
                          nr_electrons,    &
                          xc_energy        &
                         )

      use block_config
      use dft_cfg
      use dft_functional_derv

!     --------------------------------------------------------------------------
      logical, intent(in)    :: is_gga
      logical, intent(in)    :: do_lr
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: block_threshold
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8)                :: buffer(*)
      real(8), intent(in)    :: ao_block(mxbllen, *)
      integer, intent(in)    :: nr_fmat
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, nr_fmat)
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim, nz)
      integer, intent(in)    :: nr_dmat
      real(8), intent(in)    :: dmat_p(mat_dim, mat_dim, nz, nr_dmat)
      integer, intent(in)    :: dmat_pg_sym(nr_dmat)
      integer, intent(in)    :: dmat_ih_sym(nr_dmat)
      integer                :: derv_length
      real(8)                :: derv(mxbllen, 0:derv_length)
      integer, intent(in)    :: max_fun_derv
      real(8)                :: w(*)
      real(8)                :: nr_electrons
      real(8)                :: xc_energy
!     --------------------------------------------------------------------------
      real(8), allocatable   :: ao_sorted(:, :)
      integer, allocatable   :: ao_sort_index(:)
      integer, allocatable   :: ao_sort_index_reverse(:)
      integer                :: irep, imat
      integer                :: i, j, k
      integer                :: ixyz
      integer                :: ikoff, jkoff, ijoff
      integer, parameter     :: max_nr_blocks = 16
      integer                :: isoff(max_nr_blocks)
      integer                :: iksoff(max_nr_blocks)
      integer                :: nr_above_threshold(max_nr_blocks)
      real(8)                :: temp
      integer, parameter     :: blas_limit = 10
      integer                :: is, js
      integer                :: nr1, nr2, st1, st2
      integer                :: iblock, jblock
      real(8), allocatable   :: n(:, :)
      real(8), allocatable   :: n_g(:, :, :)
      real(8), allocatable   :: p_n(:, :)
      real(8), allocatable   :: p_n_g(:, :, :)
      real(8), allocatable   :: gnn(:, :, :)
      integer, allocatable   :: jblock_index(:)
!     --------------------------------------------------------------------------

      allocate(ao_sort_index(mat_dim))
      allocate(ao_sort_index_reverse(mat_dim))
      allocate(ao_sorted(mxbllen*mat_dim, 0:4))
      allocate(n(block_length, 0:nr_dmat))
      allocate(n_g(block_length, 3, 0:nr_dmat))
      allocate(p_n(block_length, maxval((/nr_dmat, 1/))))
      allocate(p_n_g(block_length, 3, maxval((/nr_dmat, 1/))))
      allocate(gnn(block_length, 0:1, 0:nr_dmat))
      allocate(jblock_index(mat_dim))

!     zero out density and gradient
      do imat = 0, nr_dmat
         do k = 1, block_length
            n(k, imat) = 0.0d0
         end do
         if (is_gga) then
            do ixyz = 1, 3
               do k = 1, block_length
                  n_g(k, ixyz, imat) = 0.0d0
               end do
            end do
         end if
      end do

      isoff              = 0
      iksoff             = 0
      nr_above_threshold = 0

      if (nr_dmat > 0) then
         irep = dmat_pg_sym(1) - 1
      else
         irep = 0
      end if

      do iblock = 1, llss_nr_density_blocks(0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)

!        get buffer for sorting
!        and get number of functions above block_threshold
         do i = 1, nr1
            buffer(i) = 0.0d0
            is = i + st1 - 1
            do k = 1, block_length
               temp = dabs(ao_block(k, is))
               if (temp > buffer(i)) then
                  buffer(i) = temp
               end if
            end do
            if (buffer(i) > block_threshold) then
               nr_above_threshold(iblock) = nr_above_threshold(iblock) + 1
            end if
         end do

!        check for early return
         if (nr_above_threshold(iblock) == 0) then
            cycle
         end if

!        shift offsets of all later blocks
         do jblock = iblock + 1, max_nr_blocks
            isoff(jblock) = isoff(jblock) + nr_above_threshold(iblock)
         end do
         do jblock = iblock + 1, max_nr_blocks
            iksoff(jblock) = isoff(jblock)*block_length
         end do

!        get sorting index
         call indexx(nr1, buffer, ao_sort_index_reverse)
         do i = 1, nr_above_threshold(iblock)
            ao_sort_index(isoff(iblock) + i) = ao_sort_index_reverse(nr1 - i + 1) + st1 - 1
         end do

!        sort ao matrix
         do i = 1, nr_above_threshold(iblock)
            is = ao_sort_index(isoff(iblock) + i)
            ikoff = (i - 1)*block_length
            do k = 1, block_length
               ao_sorted(iksoff(iblock) + ikoff + k, 0) = ao_block(k, is)
            end do
            if (is_gga) then
               do ixyz = 1, 3
                  do k = 1, block_length
                     ao_sorted(iksoff(iblock) + ikoff + k, ixyz) = ao_block(k, is + mat_dim*ixyz)
                  end do
               end do
            end if
         end do
      end do

      do iblock = 1, llss_nr_density_blocks(0, irep)
         st1 = llss_st_ao(1, iblock, 0, irep)
         jblock_index(st1) = iblock
      end do

      do iblock = 1, llss_nr_density_blocks(0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
         jblock = jblock_index(st2)
         if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
            call evaluate_densities(                            &
                                    do_lr,                      &
                                    is_gga,                     &
                                    block_length,               &
                                    mat_dim,                    &
                                    nz,                         &
                                    ao_sort_index,              &
                                    ao_sorted,                  &
                                    dmat_0,                     &
                                    nr_dmat,                    &
                                    dmat_p,                     &
                                    dmat_ih_sym,                &
                                    nr_above_threshold(iblock), &
                                    isoff(iblock),              &
                                    iksoff(iblock),             &
                                    nr_above_threshold(jblock), &
                                    isoff(jblock),              &
                                    iksoff(jblock),             &
                                    n,                          &
                                    n_g                         &
                                   )
         end if
      end do

      if (is_gga) then
         do k = 1, block_length
            gnn(k, 0, 0) = n_g(k, 1, 0)*n_g(k, 1, 0) &
                         + n_g(k, 2, 0)*n_g(k, 2, 0) &
                         + n_g(k, 3, 0)*n_g(k, 3, 0)
         
         end do
         do imat = 1, nr_dmat
            do k = 1, block_length
               gnn(k, 0, imat) = 2.0d0*n_g(k, 1, 0)*n_g(k, 1, imat) &
                               + 2.0d0*n_g(k, 2, 0)*n_g(k, 2, imat) &
                               + 2.0d0*n_g(k, 3, 0)*n_g(k, 3, imat)
            
            end do
         end do
      else
         do k = 1, block_length
            gnn(k, 0, 0) = 0.0d0
         end do
      end if

      call get_functional_derv(xc_fun,       &
                               max_fun_derv, &
                               block_length, &
                               w,            &
                               n(1, 0),      &
                               gnn(1, 0, 0), &
                               derv_length,  &
                               derv)

      do k = 1, block_length
         nr_electrons = nr_electrons + w(k)*n(k, 0)
         xc_energy    = xc_energy    + derv(k, d0000000)
      end do

!     form prefactors
      if (nr_dmat > 0) then
         do imat = 1, nr_dmat
            do k = 1, block_length
               p_n(k, imat) = derv(k, d2000000)*n(k, imat)
            end do
            if (is_gga) then
               do k = 1, block_length
                  p_n(k, imat) = p_n(k, imat) + derv(k, d1010000)*gnn(k, 0, imat)
               end do
               do ixyz = 1, 3
                  do k = 1, block_length
                     p_n_g(k, ixyz, imat) = 2.0d0*n_g(k, ixyz, 0)*n(k, imat)*derv(k, d1010000)      &
                                          + 2.0d0*n_g(k, ixyz, 0)*gnn(k, 0, imat)*derv(k, d0020000) &
                                          + 2.0d0*n_g(k, ixyz, imat)*derv(k, d0010000)
                  end do
               end do
            end if
         end do
      else
         do k = 1, block_length
            p_n(k, 1) = derv(k, d1000000)
         end do
         if (is_gga) then
            do ixyz = 1, 3
               do k = 1, block_length
                  p_n_g(k, ixyz, 1) = 2.0d0*n_g(k, ixyz, 0)*derv(k, d0010000)
               end do
            end do
         end if
      end if

      do iblock = 1, llss_nr_density_blocks(0, irep)
         st2 = llss_st_ao(2, iblock, 0, irep)
         jblock = jblock_index(st2)
         if ((nr_above_threshold(iblock) > 0) .and. (nr_above_threshold(jblock) > 0)) then
            call distribute_prefactors(                            &
                                       is_gga,                     &
                                       block_length,               &
                                       mat_dim,                    &
                                       nz,                         &
                                       ao_sort_index,              &
                                       nr_fmat,                    &
                                       fmat,                       &
                                       ao_sorted,                  &
                                       nr_above_threshold(iblock), &
                                       isoff(iblock),              &
                                       iksoff(iblock),             &
                                       nr_above_threshold(jblock), &
                                       isoff(jblock),              &
                                       iksoff(jblock),             &
                                       p_n,                        &
                                       p_n_g                       &
                                      )
         end if
      end do

      deallocate(ao_sort_index)
      deallocate(ao_sort_index_reverse)
      deallocate(ao_sorted)
      deallocate(p_n)
      deallocate(p_n_g)
      deallocate(n)
      deallocate(n_g)
      deallocate(gnn)
      deallocate(jblock_index)

   end subroutine

   subroutine distribute_prefactors(                      &
                                    is_gga,               &
                                    block_length,         &
                                    mat_dim,              &
                                    nz,                   &
                                    ao_sort_index,        &
                                    nr_fmat,              &
                                    fmat,                 &
                                    ao_sorted,            &
                                    nr_above_threshold_l, &
                                    isoff_l,              &
                                    iksoff_l,             &
                                    nr_above_threshold_r, &
                                    isoff_r,              &
                                    iksoff_r,             &
                                    p_n,                  &
                                    p_n_g                 &
                                   )

      use block_config
!     --------------------------------------------------------------------------
      logical                :: is_gga
      integer                :: block_length
      integer                :: mat_dim
      integer                :: nz
      integer                :: ao_sort_index(*)
      real(8)                :: ao_sorted(mxbllen*mat_dim, 0:4)
      integer                :: nr_fmat
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, nr_fmat)
      integer                :: nr_above_threshold_l
      integer                :: isoff_l
      integer                :: iksoff_l
      integer                :: nr_above_threshold_r
      integer                :: isoff_r
      integer                :: iksoff_r
      real(8)                :: p_n(block_length, nr_fmat)
      real(8)                :: p_n_g(block_length, 3, nr_fmat)
!     --------------------------------------------------------------------------
      integer :: i, k, ixyz, j, is, ijoff, js, ikoff, imat

      real(8), allocatable   :: xmat(:)
      real(8), allocatable   :: tmat(:)
!     --------------------------------------------------------------------------

      allocate(xmat(nr_above_threshold_l*nr_above_threshold_r))
      allocate(tmat(mxbllen*mat_dim))

      do imat = 1, nr_fmat

         do i = 1, nr_above_threshold_l
            ikoff = (i - 1)*block_length
            do k = 1, block_length
               tmat(ikoff + k) = ao_sorted(iksoff_l + ikoff + k, 0)*p_n(k, imat)
            end do
            if (is_gga) then
               do ixyz = 1, 3
                  do k = 1, block_length
                     tmat(ikoff + k) = tmat(ikoff + k) + 2.0d0*ao_sorted(iksoff_l + ikoff + k, ixyz)*p_n_g(k, ixyz, imat)
                  end do
               end do
            end if
         end do
         
         if (iksoff_l == iksoff_r) then !diagonal block
            call dsyr2k('u',                        & !uplo 
                        't',                        & !trans
                        nr_above_threshold_l,       & !n    
                        block_length,               & !k    
                        1.0d0,                      & !alpha
                        ao_sorted(iksoff_l + 1, 0), & !a    
                        block_length,               & !lda  
                        tmat,                       & !b    
                        block_length,               & !ldb  
                        0.0d0,                      & !beta 
                        xmat,                       & !c    
                        nr_above_threshold_l)         !ldc  
         else !off-diagonal block
            call dgemm('t',                         & !transa
                       'n',                         & !transb 
                       nr_above_threshold_r,        & !m      
                       nr_above_threshold_l,        & !n      
                       block_length,                & !k      
                       2.0d0,                       & !alpha  
                       ao_sorted(iksoff_r + 1, 0),  & !a    
                       block_length,                & !lda    
                       tmat,                        & !b      
                       block_length,                & !ldb    
                       0.0d0,                       & !beta   
                       xmat,                        & !c       
                       nr_above_threshold_r)          !ldc     
         end if
         
!        unsort onto fmat
         if (iksoff_l == iksoff_r) then !diagonal block
            do i = 1, nr_above_threshold_l
               is = ao_sort_index(isoff_l + i)
               ijoff = (i - 1)*nr_above_threshold_l
!              diagonal
               fmat(is, is, 1, imat) = fmat(is, is, 1, imat) + 0.5d0*xmat(ijoff + i)
               do j = 1, i - 1
                  js = ao_sort_index(isoff_l + j)
!                 upper triangle
                  fmat(js, is, 1, imat) = fmat(js, is, 1, imat) + xmat(ijoff + j)
               end do
            end do
         else !off-diagonal block
            do i = 1, nr_above_threshold_l
               is = ao_sort_index(isoff_l + i)
               ijoff = (i - 1)*nr_above_threshold_r
               do j = 1, nr_above_threshold_r
                  js = ao_sort_index(isoff_r + j)
                  fmat(js, is, 1, imat) = fmat(js, is, 1, imat) + 0.5d0*xmat(ijoff + j)
               end do
            end do
         end if

      end do

      deallocate(xmat)
      deallocate(tmat)

   end subroutine

   subroutine evaluate_densities(                      &
                                 do_lr,                &
                                 is_gga,               &
                                 block_length,         &
                                 mat_dim,              &
                                 nz,                   &
                                 ao_sort_index,        &
                                 ao_sorted,            &
                                 dmat_0,               &
                                 nr_dmat,              &
                                 dmat_p,               &
                                 dmat_ih_sym,          &
                                 nr_above_threshold_l, &
                                 isoff_l,              &
                                 iksoff_l,             &
                                 nr_above_threshold_r, &
                                 isoff_r,              &
                                 iksoff_r,             &
                                 n,                    &
                                 n_g                   &
                                )

      use block_config
!     --------------------------------------------------------------------------
      logical, intent(in)    :: do_lr
      logical, intent(in)    :: is_gga
      integer, intent(in)    :: block_length
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      integer, intent(in)    :: ao_sort_index(*)
      real(8), intent(in)    :: ao_sorted(mxbllen*mat_dim, 0:4)
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim, nz)
      integer, intent(in)    :: nr_dmat
      real(8), intent(in)    :: dmat_p(mat_dim, mat_dim, nz, nr_dmat)
      integer, intent(in)    :: dmat_ih_sym(nr_dmat)
      integer, intent(in)    :: nr_above_threshold_l
      integer, intent(in)    :: isoff_l
      integer, intent(in)    :: iksoff_l
      integer, intent(in)    :: nr_above_threshold_r
      integer, intent(in)    :: isoff_r
      integer, intent(in)    :: iksoff_r
      real(8), intent(inout) :: n(block_length, 0:nr_dmat)
      real(8), intent(inout) :: n_g(block_length, 3, 0:nr_dmat)
!     --------------------------------------------------------------------------
      integer :: i, ikoff, ixyz, k, jkoff, j, ijoff, is, js, imat
      real(8), allocatable   :: dmat_sorted(:)
      real(8), allocatable   :: xmat(:)
!     --------------------------------------------------------------------------

      allocate(xmat(mat_dim*mxbllen))

      allocate(dmat_sorted(nr_above_threshold_l*nr_above_threshold_l))

      do i = 1, nr_above_threshold_l
         is = ao_sort_index(isoff_l + i)
         ijoff = (i - 1)*nr_above_threshold_l
         do j = 1, i
            js = ao_sort_index(isoff_l + j)
            dmat_sorted(ijoff + j) = dmat_0(js, is, 1)
         end do
      end do

      call dsymm('r',                        & !side   
                 'u',                        & !uplo   
                 block_length,               & !m      
                 nr_above_threshold_l,       & !n      
                 1.0d0,                      & !alpha  
                 dmat_sorted,                & !a      
                 nr_above_threshold_l,       & !lda    
                 ao_sorted(iksoff_l + 1, 0), & !b      
                 block_length,               & !ldb    
                 0.0d0,                      & !beta   
                 xmat,                       & !c      
                 block_length)                 !ldc     

      deallocate(dmat_sorted)

      do i = 1, nr_above_threshold_l
         ikoff = (i - 1)*block_length
         do k = 1, block_length
            n(k, 0) = n(k, 0) + xmat(ikoff + k)*ao_sorted(iksoff_l + ikoff + k, 0)
         end do
         if (is_gga) then
            do ixyz = 1, 3
               do k = 1, block_length
                  n_g(k, ixyz, 0) = n_g(k, ixyz, 0) + 2.0d0*xmat(ikoff + k)*ao_sorted(iksoff_l + ikoff + k, ixyz)
               end do
            end do
         end if
      end do

      if (nr_dmat > 0) then

         allocate(dmat_sorted(nr_above_threshold_l*nr_above_threshold_r))

         do imat = 1, nr_dmat

            if (dmat_ih_sym(imat) == -1) cycle

            if (iksoff_l == iksoff_r) then !diagonal block
               do i = 1, nr_above_threshold_l
                  is = ao_sort_index(isoff_l + i)
                  ijoff = (i - 1)*nr_above_threshold_l
                  do j = 1, i
                     js = ao_sort_index(isoff_l + j)
                     dmat_sorted(ijoff + j) = dmat_p(js, is, 1, imat)
                  end do
               end do
            else !off-diagonal block
               do i = 1, nr_above_threshold_r
                  is = ao_sort_index(isoff_r + i)
                  ijoff = (i - 1)*nr_above_threshold_l
                  do j = 1, nr_above_threshold_l
                     js = ao_sort_index(isoff_l + j)
                     dmat_sorted(ijoff + j) = dmat_p(js, is, 1, imat)
                  end do
               end do
            end if

            if (iksoff_l == iksoff_r) then !diagonal block
               call dsymm('r',                            & !side   
                          'u',                            & !uplo   
                          block_length,                   & !m      
                          nr_above_threshold_l,           & !n      
                          1.0d0,                          & !alpha  
                          dmat_sorted,                    & !a      
                          nr_above_threshold_l,           & !lda    
                          ao_sorted(iksoff_l + 1, 0),     & !b      
                          block_length,                   & !ldb    
                          0.0d0,                          & !beta   
                          xmat,                           & !c      
                          block_length)                     !ldc     
            else !off-diagonal block
               call dgemm('n',                            & !transa
                          'n',                            & !transb 
                          block_length,                   & !m      
                          nr_above_threshold_r,           & !n      
                          nr_above_threshold_l,           & !k      
                          1.0d0,                          & !alpha  
                          ao_sorted(iksoff_l + 1, 0),     & !a      
                          block_length,                   & !lda    
                          dmat_sorted,                    & !b    
                          nr_above_threshold_l,           & !ldb
                          0.0d0,                          & !beta   
                          xmat,                           & !c       
                          block_length)                     !ldc     
            end if

            do i = 1, nr_above_threshold_r
               ikoff = (i - 1)*block_length
               do k = 1, block_length
                  n(k, imat) = n(k, imat) + xmat(ikoff + k)*ao_sorted(iksoff_r + ikoff + k, 0)
               end do
               if (is_gga) then
                  do ixyz = 1, 3
                     do k = 1, block_length
                        n_g(k, ixyz, imat) = n_g(k, ixyz, imat) + 2.0d0*xmat(ikoff + k)*ao_sorted(iksoff_r + ikoff + k, ixyz)
                     end do
                  end do
               end if
            end do

         end do

         deallocate(dmat_sorted)

      end if

      deallocate(xmat)

   end subroutine

   function get_block_threshold(mat_dim, &
                                dmat_0)

      use dft_cfg

!     --------------------------------------------------------------------------
      real(8)                :: get_block_threshold
!     --------------------------------------------------------------------------
      integer, intent(in)    :: mat_dim
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
!     --------------------------------------------------------------------------
      integer                :: i, j, iblock, st, nr
      real(8)                :: f, x
!     --------------------------------------------------------------------------

      x = 0.0d0
      iblock_loop: do iblock = 1, llss_nr_density_blocks(0, 0)
         st = llss_st_ao(1, iblock, 0, 0)
         nr = llss_nr_ao(1, iblock, 0, 0)
         do i = st, st + nr - 1
            do j = st, i
!              dmat_0 is positive definite but just to be sure
               f = dabs(dmat_0(j, i))
               if (f > x) then
                  x = f
               end if
            end do
         end do
      end do iblock_loop
!radovan: my own empirical estimate, no physics behind this
      get_block_threshold = dsqrt(10.0d0*dft_cfg_dfthr0/(mat_dim*mat_dim*x))

   end function

end module
