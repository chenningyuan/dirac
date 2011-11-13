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

module density_eval

   use basis_info

   implicit none

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
   public lda_london_distribute_i
   public get_s_lao
   public get_n_block
   public get_block_threshold

   private

  integer :: jxyz_start = 1 !noncollinear
! integer :: jxyz_start = 3 !collinear

  save

! cval, will not work in parallel!
#include "dcbgen.h"

contains

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


!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s
    integer, intent(in)    :: isym(*)
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao_b1(mat_dim, *)
    real(8), intent(in)    :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz, jj, j
    integer                :: a, b, c
    real(8)                :: f
    real(8)                :: r_co(3)
!   ----------------------------------------------------------------------------

!   r_co(1) = cx - dft_origin(1)
!   r_co(2) = cy - dft_origin(2)
!   r_co(3) = cz - dft_origin(3)
    r_co(1) = cx
    r_co(2) = cy
    r_co(3) = cz

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

          do ii = st1, st1 + nr1 - 1
            do jj = st2, st2 + nr2 - 1

              f = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

              fmat(ii, jj, (a-1)*nz + 1) = fmat(ii, jj, (a-1)*nz + 1) + s*f

            end do
          end do
        end do

    end do

  end subroutine

  subroutine gga_london_distribute_r(d_r,      &
                                     d_z,      &
                                     n_g1,     &
                                     irep,     &
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


!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: d_r
    real(8), intent(in)    :: d_z
    real(8), intent(in)    :: n_g1(3)
    integer, intent(in)    :: irep
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao_g0_b0(mat_dim)
    real(8), intent(in)    :: ao_g0_b1(mat_dim, 3)
    real(8), intent(in)    :: ao_g1_b0(mat_dim, 3)
    real(8), intent(in)    :: ao_g1_b1(mat_dim, 3, 3)
    real(8), intent(in)    :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, jxyz, i, j
    integer                :: a, b, c, ixyz
    real(8)                :: f, f_temp
    real(8)                :: r(3)
!   ----------------------------------------------------------------------------

!   r(1) = cx - dft_origin(1)
!   r(2) = cy - dft_origin(2)
!   r(3) = cz - dft_origin(3)
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

        do iblock = 1, llss_nr_density_blocks(0, irep)

          nr1 = llss_nr_ao(1, iblock, 0, irep)
          st1 = llss_st_ao(1, iblock, 0, irep)
          nr2 = llss_nr_ao(2, iblock, 0, irep)
          st2 = llss_st_ao(2, iblock, 0, irep)

          do i = st1, st1 + nr1 - 1
            do j = st2, st2 + nr2 - 1

              f_temp = (r(c)*ao_g0_b1(j, b) - r(b)*ao_g0_b1(j, c))*ao_g0_b0(i) &
                     - (r(c)*ao_g0_b1(i, b) - r(b)*ao_g0_b1(i, c))*ao_g0_b0(j)

              f = f_temp*d_r

!             so far it's lda type contribution, rest are gradient terms

              f_temp = (n_g1(c)*ao_g0_b1(j, b) - n_g1(b)*ao_g0_b1(j, c))*ao_g0_b0(i) &
                     - (n_g1(c)*ao_g0_b1(i, b) - n_g1(b)*ao_g0_b1(i, c))*ao_g0_b0(j)

              f_temp = f_temp                                                               &
                     + (n_g1(1)*ao_g1_b0(j, 1)                                              &
                     +  n_g1(2)*ao_g1_b0(j, 2)                                              &
                     +  n_g1(3)*ao_g1_b0(j, 3))*(r(b)*ao_g0_b1(i, c) - r(c)*ao_g0_b1(i, b)) &
                     - (n_g1(1)*ao_g1_b0(i, 1)                                              &
                     +  n_g1(2)*ao_g1_b0(i, 2)                                              &
                     +  n_g1(3)*ao_g1_b0(i, 3))*(r(b)*ao_g0_b1(j, c) - r(c)*ao_g0_b1(j, b))

              f = f + f_temp*d_z

              do ixyz = 1, 3
                f_temp = (r(c)*ao_g1_b1(j, b, ixyz) - r(b)*ao_g1_b1(j, c, ixyz))*ao_g0_b0(i) &
                       - (r(c)*ao_g1_b1(i, b, ixyz) - r(b)*ao_g1_b1(i, c, ixyz))*ao_g0_b0(j)

                f = f + f_temp*d_z*n_g1(ixyz)
              end do

              fmat(i, j, (a-1)*nz + 1) = fmat(i, j, (a-1)*nz + 1) + f

            end do
          end do
      end do

    end do

  end subroutine

  subroutine lda_london_distribute_i(s, isym, mat_dim, nz, fmat, ao, ao_b1, cx, cy, cz)

!   ----------------------------------------------------------------------------
    real(8), intent(in)    :: s(3)
    integer, intent(in)    :: isym(*)
    integer, intent(in)    :: mat_dim
    integer, intent(in)    :: nz
    real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
    real(8), intent(in)    :: ao(mat_dim)
    real(8), intent(in)    :: ao_b1(mat_dim, *)
    real(8), intent(in)    :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz, jj, j
    integer                :: a, b, c
    real(8)                :: f
    real(8)                :: r_co(3)
!   ----------------------------------------------------------------------------

    r_co(1) = cx
    r_co(2) = cy
    r_co(3) = cz

    do a = 1, 3
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

              f = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

              !fmat(ii, jj, (a-1)*nz + 1) = fmat(ii, jj, (a-1)*nz + 1) + s(jxyz)*f
              fmat(ii, jj, iq(jxyz)) = fmat(ii, jj, iq(jxyz)) + s(jxyz)*f

            end do
          end do
        end do
      end do
    end do

  end subroutine

  subroutine gga_london_distribute_i()
  end subroutine


  subroutine get_s_lao(s, isym, mat_dim, dmat, ao, ao_b1, cx, cy, cz)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3, 3)
    integer, intent(in)  :: isym(*)
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: ao(mat_dim)
    real(8), intent(in)  :: ao_b1(mat_dim, *)
    real(8), intent(in)  :: cx, cy, cz
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz
    integer              :: a, b, c, irp1, irp2, ii, jj
    real(8)              :: r_co(3), f
!   ----------------------------------------------------------------------------

    r_co(1) = cx
    r_co(2) = cy
    r_co(3) = cz

    s = 0.0d0
    f = 0.0d0

    do a = 1, 3
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

              f = (r_co(c)*ao_b1(jj, b) - r_co(b)*ao_b1(jj, c))*ao(ii) &
                - (r_co(c)*ao_b1(ii, b) - r_co(b)*ao_b1(ii, c))*ao(jj)

              s(jxyz, a) = s(jxyz, a) + dmat(ii, jj, iq(jxyz))*f

            end do
          end do
        end do
      end do

    end do


  end subroutine


   subroutine get_n_block( &
                          block_length,          &
                          block_threshold,       &
                          mat_dim,               &
                          n,                     &
                          buffer,                &
                          ao_block,              &
                          ao_sorted,             &
                          bmat,                  &
                          ao_sort_index,         &
                          ao_sort_index_reverse, &
                          fmat,                  &
                          dmat_0,                &
                          dmat_sorted,           &
                          xmat,                  &
                          irep,                  &
                          derv_length,           &
                          derv,                  &
                          max_fun_derv,          &
                          gnn,                   &
                          w)

      use block_config
      use dft_cfg
      use dft_functional_derv

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: block_threshold
      integer, intent(in)    :: mat_dim
      real(8), intent(out)   :: n(*)
      real(8)                :: buffer(*)
      real(8), intent(in)    :: ao_block(mxbllen, *)
      real(8)                :: ao_sorted(*)
      real(8)                :: bmat(*)
      integer                :: ao_sort_index(*)
      integer                :: ao_sort_index_reverse(*)
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, *)
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim)
      real(8)                :: dmat_sorted(*)
      real(8)                :: xmat(*)
      integer, intent(in)    :: irep
      integer                :: derv_length
      real(8)                :: derv(mxbllen, 0:derv_length)
      integer, intent(in)    :: max_fun_derv
      real(8)                :: gnn(*)
      real(8)                :: w(*)
!     --------------------------------------------------------------------------
      integer                :: i, j, k
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
!     --------------------------------------------------------------------------


      isoff = 0
      iksoff = 0

!     zero out density
      do k = 1, block_length
         n(k) = 0.0d0
      end do

      iblock_loop_1: do iblock = 1, llss_nr_density_blocks(0, irep)

         st1 = llss_st_ao(1, iblock, 0, irep)
         nr1 = llss_nr_ao(1, iblock, 0, irep)
!        st2 = llss_st_ao(2, iblock, 0, irep)
!        nr2 = llss_nr_ao(2, iblock, 0, irep)

!        if irep == 0 everything is along diagonal and square
!        for other ireps we will have to do some tricks
         if (irep /= 0) then
            call quit('implement irep /= 0')
         end if

!        get buffer for sorting
!        and get number of functions above block_threshold
         nr_above_threshold(iblock) = 0
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

!        save offsets
         do jblock = iblock + 1, max_nr_blocks
            isoff(jblock) = isoff(jblock) + nr_above_threshold(iblock)
         end do
         do jblock = iblock + 1, max_nr_blocks
            iksoff(jblock) = isoff(jblock)*block_length
         end do

!        get sorting index
         call indexx(nr1, buffer, ao_sort_index_reverse)
         do i = 1, nr_above_threshold(iblock)
            ao_sort_index(isoff(iblock) + i) = ao_sort_index_reverse(nr1 - i + 1)
         end do

!        sort ao matrix
         do i = 1, nr_above_threshold(iblock)
            is = ao_sort_index(isoff(iblock) + i) + st1 - 1
            ikoff = (i - 1)*block_length
            do k = 1, block_length
               ao_sorted(iksoff(iblock) + ikoff + k) = ao_block(k, is)
            end do
         end do

!        sort density matrix (upper triangle)
!        note that the lower triangle of the sorted density matrix
!        can contain numerical dirt from previous batch or previous
!        symmetry block but it is not used
         do i = 1, nr_above_threshold(iblock)
            is = ao_sort_index(isoff(iblock) + i) + st1 - 1
            ijoff = (i - 1)*nr_above_threshold(iblock)
            do j = 1, i
               js = ao_sort_index(isoff(iblock) + j) + st1 - 1
               dmat_sorted(ijoff + j) = dmat_0(js, is)
            end do
         end do

         if (nr_above_threshold(iblock) > blas_limit) then
            call dsymm('r',                            & !side   
                       'u',                            & !uplo   
                       block_length,                   & !m      
                       nr_above_threshold(iblock),     & !n      
                       1.0d0,                          & !alpha  
                       dmat_sorted,                    & !a      
                       nr_above_threshold(iblock),     & !lda    
                       ao_sorted(iksoff(iblock) + 1),  & !b      
                       block_length,                   & !ldb    
                       0.0d0,                          & !beta   
                       xmat,                           & !c      
                       block_length)                     !ldc     

            do i = 1, nr_above_threshold(iblock)
               ikoff = (i - 1)*block_length
               do k = 1, block_length
                  n(k) = n(k) + xmat(ikoff + k)*ao_sorted(iksoff(iblock) + ikoff + k)
               end do
            end do
         else
            do i = 1, nr_above_threshold(iblock)
               ijoff = (i - 1)*nr_above_threshold(iblock)
               ikoff = (i - 1)*block_length
!              triangle except diagonal
               do j = 1, i - 1
                  jkoff = (j - 1)*block_length
                  do k = 1, block_length
                     n(k) = n(k) + 2.0d0*ao_sorted(iksoff(iblock) + ikoff + k) &
                               *ao_sorted(iksoff(iblock) + jkoff + k)*dmat_sorted(ijoff + j)
                  end do
               end do
!              diagonal
               do k = 1, block_length
                  n(k) = n(k) + ao_sorted(iksoff(iblock) + ikoff + k)&
                                *ao_sorted(iksoff(iblock) + ikoff + k)*dmat_sorted(ijoff + i)
               end do
            end do
         end if

      end do iblock_loop_1

      call get_functional_derv(xc_fun,       &
                               max_fun_derv, &
                               block_length, &
                               w,            &
                               n,            &
                               gnn,          &
                               derv_length,  &
                               derv)

      iblock_loop_2: do iblock = 1, llss_nr_density_blocks(0, irep)

         st1 = llss_st_ao(1, iblock, 0, irep)

         do i = 1, nr_above_threshold(iblock)
            ikoff = (i - 1)*block_length
            do k = 1, block_length
               bmat(ikoff + k) = ao_sorted(iksoff(iblock) + ikoff + k)*derv(k, d1000000)
            end do
         end do
          
         if (nr_above_threshold(iblock) > blas_limit) then
            call dsyr2k('u',                           & !uplo 
                        't',                           & !trans
                        nr_above_threshold(iblock),    & !n    
                        block_length,                  & !k    
                        1.0d0,                         & !alpha
                        ao_sorted(iksoff(iblock) + 1), & !a    
                        block_length,                  & !lda  
                        bmat,                          & !b    
                        block_length,                  & !ldb  
                        0.0d0,                         & !beta 
                        xmat,                          & !c    
                        nr_above_threshold(iblock))      !ldc  
!           call dgemm('t',                            & !transa
!                      'n',                            & !transb 
!                      nr_above_threshold(iblock),     & !m      
!                      nr_above_threshold(iblock),     & !n      
!                      block_length,                   & !k      
!                      2.0d0,                          & !alpha  
!                      ao_sorted(iksoff(iblock) + 1),  & !a    
!                      block_length,                   & !lda    
!                      bmat,                           & !b      
!                      block_length,                   & !ldb    
!                      0.0d0,                          & !beta   
!                      xmat,                           & !c       
!                      nr_above_threshold(iblock))       !ldc     
         else
            do i = 1, nr_above_threshold(iblock)
               ijoff = (i - 1)*nr_above_threshold(iblock)
               ikoff = (i - 1)*block_length
!              triangle except diagonal
               do j = 1, i - 1
                  xmat(ijoff + j) = 0.0d0
                  jkoff = (j - 1)*block_length
                  do k = 1, block_length
                     xmat(ijoff + j) = xmat(ijoff + j) &
                                     + 2.0d0*ao_sorted(iksoff(iblock) + ikoff + k)*bmat(jkoff + k)
                  end do
               end do
!              diagonal
               xmat(ijoff + i) = 0.0d0
               do k = 1, block_length
                  xmat(ijoff + i) = xmat(ijoff + i) &
                                  + 2.0d0*ao_sorted(iksoff(iblock) + ikoff + k)*bmat(ikoff + k)
               end do
            end do
         end if

!        unsort onto fmat
         do i = 1, nr_above_threshold(iblock)
            is = ao_sort_index(isoff(iblock) + i) + st1 - 1
            ijoff = (i - 1)*nr_above_threshold(iblock)
!           diagonal
            fmat(is, is, 1) = fmat(is, is, 1) + 0.5d0*xmat(ijoff + i)
            do j = 1, i - 1
               js = ao_sort_index(isoff(iblock) + j) + st1 - 1
!              upper triangle
               fmat(js, is, 1) = fmat(js, is, 1) + xmat(ijoff + j)
            end do
         end do

      end do iblock_loop_2

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
