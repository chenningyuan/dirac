!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module density_eval

   use interface_ao
   use dirac_ao_eval

   implicit none

   public mat_lao_sym
   public gga_sym
   public get_n
   public get_gamma5
   public get_s
   public get_j
   public get_kin
   public get_kin_ls
   public get_kin_sl
   public get_kin_lap
   public get_kin_tau
   public get_kin_tau_imag
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
   public get_jblao
   public get_gjblao
   public get_jblao_ll

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
    integer                :: i, j, ii, jj
    integer, parameter     :: h(4) = (/-1, 1, 1, 1/)
    integer                :: st1, st2, nr1, nr2
    real(8)                :: average
!   ----------------------------------------------------------------------------

    do imat = 1, nr_mat

       irep = isym(imat) - 1

       do iz = 1, nz
          do iblock = 1, nr_ao_blocks
             nr1 = ao_block_nr(iblock)
             st1 = ao_block_start(iblock)
             nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
             st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
             if ((nr1 < 1) .or. (nr2 < 1)) cycle

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
       if (ih(imat) ==  1) skip = (/.false., .true.,  .true.,  .true./)
       if (ih(imat) == -1) skip = (/.true.,  .false., .false., .false./)

       irep = isym(imat) - 1

       do jxyz = 0, 3
          if (.not. skip(jxyz)) then
             do iblock = 1, nr_ao_blocks
                nr1 = ao_block_nr(iblock)
                st1 = ao_block_start(iblock)
                nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
                st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
                if ((nr1 < 1) .or. (nr2 < 1)) cycle

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
      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
         st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
    
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
      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(lssl_block_partner(iblock, 0, irep))
         st2 = ao_block_start(lssl_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
    
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
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
     
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
                    + llss_prefactor(iblock)*ddot(nr1, ao(st1), 1, buffer, 1)
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
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(lssl_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(lssl_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
     
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
                    ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(*)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8), external    :: ddot
!   ----------------------------------------------------------------------------

    n  = 0.0d0
    gn = 0.0d0

    do iblock = 1, nr_ao_blocks
       nr1 = ao_block_nr(iblock)
       st1 = ao_block_start(iblock)
       nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
       st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
       if ((nr1 < 1) .or. (nr2 < 1)) cycle

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

       gn(1) = gn(1) + ddot(nr1, ao(ao_off_g1_m0(1, 0) + st1), 1, buffer, 1)
       gn(2) = gn(2) + ddot(nr1, ao(ao_off_g1_m0(2, 0) + st1), 1, buffer, 1)
       gn(3) = gn(3) + ddot(nr1, ao(ao_off_g1_m0(3, 0) + st1), 1, buffer, 1)
    end do

    gn = 2.0d0*gn

  end subroutine

   subroutine get_kin_tau(t,       &
                      irep,    &
                      mat_dim, &
                      dmat,    &
                      buffer,  &
                      ao)

!     --------------------------------------------------------------------------
!     radovan bast 2010-04-15
!                  this routine calculates the \tau = \nabla_i \phi_k \nabla_i \phi_k
!     
!                  integral approaches the kinetic energy as c -> oo
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      t = 0.0d0
      
      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
         st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
      
         do ixyz = 1, 3
            call dgemv('N',                   &
                       nr1,                   &
                       nr2,                   &
                       1.0d0,                 &
                       dmat(st1, st2, iq(0)), &
                       mat_dim,               &
                       ao(ao_off_g1_m0(ixyz, 0) + st2),      &
                       1,                     &
                       0.0d0,                 &
                       buffer,                &
                       1)
      
            t = t + ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
         end do
      end do
      
      t = 0.5d0*t

   end subroutine

   subroutine get_kin_tau_imag(t,       &
                           irep,    &
                           mat_dim, &
                           dmat,    &
                           buffer,  &
                           ao)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: t(3)
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, ixyz, jxyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      t = 0.0d0
     
      do jxyz = 1, 3 
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
         
            do ixyz = 1, 3
               call dgemv('N',                      &
                          nr1,                      &
                          nr2,                      &
                          1.0d0,                    &
                          dmat(st1, st2, iq(jxyz)), &
                          mat_dim,                  &
                          ao(ao_off_g1_m0(ixyz, 0) + st2), &
                          1,                        &
                          0.0d0,                    &
                          buffer,                   &
                          1)
         
               t(jxyz) = t(jxyz) &
             + llss_prefactor(iblock)*ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
            end do
         end do
      end do
      
      t = 0.5d0*t

   end subroutine

   subroutine get_kin_lap(t,       &
                          irep,    &
                          mat_dim, &
                          dmat,    &
                          buffer,  &
                          ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz
      real(8), external    :: ddot
      integer, parameter   :: pos(9) = (/1, 2, 3, 2, 4, 5, 3, 5, 6/)
!     --------------------------------------------------------------------------
   
      t = 0.0d0
      
      do jxyz = 1, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
            call dgemv('N',               &
                       nr1,               &
                       nr2,               &
                       1.0d0,             &
                       dmat(st1, st2, 1), &
                       mat_dim,           &
                       ao(st2),           &
                       1,                 &
                       0.0d0,             &
                       buffer,            &
                       1)
            t = t - ddot(nr1, ao(ao_off_g2_m0(pos((jxyz-1)*3 + jxyz), 0) + st1), 1, buffer, 1)
         end do
      end do
      
      t = t*0.5d0
   
   end subroutine

   subroutine get_kin(t,       &
                      cval,    &
                      irep,    &
                      mat_dim, &
                      dmat,    &
                      buffer,  &
                      ao)
   
!     --------------------------------------------------------------------------
!     radovan bast 2010-04-15
!                  this routine calculates the kinetic energy density
!                  c \psi+(r) \alpha \cdot p \psi(r)
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
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
   
      t = 0.0d0
      
      do jxyz = 1, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(lssl_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(lssl_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
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
            t = t + ddot(nr1, ao(ao_off_g1_m0(jxyz, 0) + st1), 1, buffer, 1)
         end do
      end do
      
      t = t*cval*0.5d0
   
   end subroutine

   subroutine get_kin_ls(t,       &
                         cval,    &
                         irep,    &
                         mat_dim, &
                         dmat,    &
                         buffer,  &
                         ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
      real(8), intent(in)  :: cval
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: jxyz
      integer              :: k, l
!     --------------------------------------------------------------------------
   
      t = 0.0d0
      
      do k = 1, nr_ao_large
         do l = nr_ao_large+1, nr_ao
            do jxyz = 1, 3
               t = t - ao(k)*ao(ao_off_g1_m0(jxyz, 0) + l)*dmat(k, l, iq(jxyz))
            end do
         end do
      end do
      
      t = t*cval*0.5d0
   
   end subroutine

   subroutine get_kin_sl(t,       &
                         cval,    &
                         irep,    &
                         mat_dim, &
                         dmat,    &
                         buffer,  &
                         ao)
   
!     --------------------------------------------------------------------------
      real(8), intent(out) :: t
      real(8), intent(in)  :: cval
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(*)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: jxyz
      integer              :: k, l
!     --------------------------------------------------------------------------
   
      t = 0.0d0
      
      do l = 1, nr_ao_large
         do k = nr_ao_large+1, nr_ao
            do jxyz = 1, 3
               t = t - ao(k)*ao(ao_off_g1_m0(jxyz, 0) + l)*dmat(k, l, iq(jxyz))
            end do
         end do
      end do
      
      t = t*cval*0.5d0
   
   end subroutine

  subroutine get_gn_nonhermitian(n,       &
                                 gn,      &
                                 irep,    &
                                 mat_dim, &
                                 dmat,    &
                                 buffer,  &
                                 ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: n
    real(8), intent(out) :: gn(3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(*)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    n  = 0.0d0
    gn = 0.0d0

    do iblock = 1, nr_ao_blocks
       nr1 = ao_block_nr(iblock)
       st1 = ao_block_start(iblock)
       nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
       st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
       if ((nr1 < 1) .or. (nr2 < 1)) cycle

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
         gn(ixyz) = gn(ixyz) + ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
       end do

       do ixyz = 1, 3
         call dgemv('N',                      &
                    nr1,                      &
                    nr2,                      &
                    1.0d0,                    &
                    dmat(st1, st2, iq(0)), &
                    mat_dim,                  &
                    ao(ao_off_g1_m0(ixyz, 0) + st2), &
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
                     ao)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: v(3)
      real(8), intent(out) :: gv(3, 3)
      real(8), intent(in)  :: cval
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      v  = 0.0d0
      gv = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(lssl_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(lssl_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
  
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
               gv(ixyz, jxyz) = gv(ixyz, jxyz) + ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
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
                     ao)

!     --------------------------------------------------------------------------
      real(8), intent(out) :: v(3)
      real(8), intent(out) :: gv(3, 3)
      integer, intent(in)  :: irep
      integer, intent(in)  :: mat_dim
      real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
      real(8), intent(in)  :: buffer(mat_dim)
      real(8), intent(in)  :: ao(*)
!     --------------------------------------------------------------------------
      integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
      real(8), external    :: ddot
!     --------------------------------------------------------------------------

      v  = 0.0d0
      gv = 0.0d0
      do jxyz = jxyz_start, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
  
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
                    + llss_prefactor(iblock)*ddot(nr1, ao(st1), 1, buffer, 1)
            do ixyz = 1, 3
               gv(ixyz, jxyz) = gv(ixyz, jxyz) &
             + llss_prefactor(iblock)*ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
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
                                 ao)

!   ----------------------------------------------------------------------------
    real(8), intent(out) :: s(3)
    real(8), intent(out) :: gs(3, 3)
    integer, intent(in)  :: irep
    integer, intent(in)  :: mat_dim
    real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
    real(8), intent(in)  :: buffer(mat_dim)
    real(8), intent(in)  :: ao(*)
!   ----------------------------------------------------------------------------
    integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
    real(8)              :: ddot
!   ----------------------------------------------------------------------------

    s = 0.0d0
    gs = 0.0d0

    do jxyz = jxyz_start, 3

       do iblock = 1, nr_ao_blocks
          nr1 = ao_block_nr(iblock)
          st1 = ao_block_start(iblock)
          nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
          st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
          if ((nr1 < 1) .or. (nr2 < 1)) cycle

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

          s(jxyz) = s(jxyz) + llss_prefactor(iblock)*ddot(nr1, ao(st1), 1, buffer, 1)
          do ixyz = 1, 3
             gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                            + llss_prefactor(iblock)*ddot(nr1, ao(ao_off_g1_m0(ixyz, 0) + st1), 1, buffer, 1)
          end do

          do ixyz = 1, 3
             call dgemv('N',                         &
                        nr1,                         &
                        nr2,                         &
                        1.0d0,                       &
                        dmat(st1, st2, iq(jxyz)), &
                        mat_dim,                     &
                        ao(ao_off_g1_m0(ixyz, 0) + st2), &
                        1,                           &
                        0.0d0,                       &
                        buffer,                      &
                        1)

             gs(ixyz, jxyz) = gs(ixyz, jxyz) &
                            + llss_prefactor(iblock)*ddot(nr1, ao(st1), 1, buffer, 1)
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
      
      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
         st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
      
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
                               ao)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u
      real(8), intent(in)    :: v(3)
      real(8), intent(in)    :: ao(*)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, irep))
         st2 = ao_block_start(llss_block_partner(iblock, 0, irep))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
      
         do i = 1, nr2
            ii = st2 + i - 1
            f =       u   *ao(ii) &
              + 2.0d0*v(1)*ao(ao_off_g1_m0(1, 0) + ii) &
              + 2.0d0*v(2)*ao(ao_off_g1_m0(2, 0) + ii) &
              + 2.0d0*v(3)*ao(ao_off_g1_m0(3, 0) + ii)
      
            call daxpy(nr1,                      &
                       f,                        &
                       ao(st1),                  &
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

         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
         
            do i = 1, nr2
               ii = st2 + i - 1
               f  = u(jxyz)*ao_g0(ii)*llss_prefactor(iblock)
         
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
                               ao)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
      integer, intent(in)    :: im
      integer, intent(in)    :: irep
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8), intent(in)    :: u(3)
      real(8), intent(in)    :: v(3, 3)
      real(8), intent(in)    :: ao(*)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, i, ii, jxyz
      real(8)                :: f
!     --------------------------------------------------------------------------
      
      do jxyz = jxyz_start, 3

         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep))
            st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
         
            do i = 1, nr2
               ii = st2 + i - 1
               f =       u(   jxyz)*ao(ii)    &
                 + 2.0d0*v(1, jxyz)*ao(ao_off_g1_m0(1, 0) + ii) &
                 + 2.0d0*v(2, jxyz)*ao(ao_off_g1_m0(2, 0) + ii) &
                 + 2.0d0*v(3, jxyz)*ao(ao_off_g1_m0(3, 0) + ii)
              
               f = llss_prefactor(iblock)*f
         
               call daxpy(nr1,                         &
                          f,                           &
                          ao(st1),                  &
                          1,                           &
                          fmat(st1, ii, iq(jxyz), im), &
                          1)
            end do
         end do
      end do

   end subroutine


   subroutine lda_london_distribute_r(s, isym, mat_dim, nz, fmat, ao, cx, cy, cz)
   
!    ----------------------------------------------------------------------------
     real(8), intent(in)    :: s
     integer, intent(in)    :: isym(*)
     integer, intent(in)    :: mat_dim
     integer, intent(in)    :: nz
     real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, *)
     real(8), intent(in)    :: ao(*)
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
     
         do iblock = 1, nr_ao_blocks
     
           nr1 = ao_block_nr(iblock)
           st1 = ao_block_start(iblock)
           nr2 = ao_block_nr(llss_block_partner(iblock, 0, isym(a) - 1))
           st2 = ao_block_start(llss_block_partner(iblock, 0, isym(a) - 1))
           if ((nr1 < 1) .or. (nr2 < 1)) cycle
     
           do i = 1, nr2
             ii = st2 + i - 1
             f  = 2.0d0*(r_co(c)*ao(ao_off_g0_m1(0, b) + ii) - r_co(b)*ao(ao_off_g0_m1(0, c) + ii))*s
     
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
                                      ao, &
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
     real(8), intent(in)    :: ao(*)
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

       do ii = 1, mat_dim
          f_vec(ii)  = r(b)*ao(ao_off_g0_m1(0, c) + ii) &
                     - r(c)*ao(ao_off_g0_m1(0, b) + ii)
       end do
   
       do iblock = 1, nr_ao_blocks
   
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, isym(a) - 1))
         st2 = ao_block_start(llss_block_partner(iblock, 0, isym(a) - 1))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
   
         do i = 1, nr2
           ii = st2 + i - 1
   
           f_temp = 2.0d0*(r(c)*ao(ao_off_g0_m1(0, b) + ii) &
                         - r(b)*ao(ao_off_g0_m1(0, c) + ii))
           f = f_temp*d_r
   
!          so far it's lda type contribution, rest are gradient terms
   
           f_temp = 2.0d0*(n_g1(c)*ao(ao_off_g0_m1(0, b) + ii) &
                         - n_g1(b)*ao(ao_off_g0_m1(0, c) + ii))
           f = f + f_temp*d_z
   
           do ixyz = 1, 3
             f_temp = 2.0d0*(r(c)*ao(ao_off_g1_m1(ixyz, b) + ii) &
                           - r(b)*ao(ao_off_g1_m1(ixyz, c) + ii))
             f = f + f_temp*d_z*n_g1(ixyz)
           end do
   
           call daxpy(nr1,                     &
                      f,                       &
                      ao(st1),                 &
                      1,                       &
                      fmat(st1, ii, iq(0), a), &
                      1)
   
   
           f_temp = 2.0d0*(n_g1(1)*ao(ao_off_g1_m0(1, 0) + ii) &
                         + n_g1(2)*ao(ao_off_g1_m0(2, 0) + ii) &
                         + n_g1(3)*ao(ao_off_g1_m0(3, 0) + ii))
   
           f = f_temp*d_z
   
           call daxpy(nr1,                          &
                      f,                            &
                      f_vec(st1),                   &
                      1,                            &
                      fmat(st1, ii, iq(0), a),      &
                      1)
   
         end do ! i = 1, nr2
       end do   ! iblock = 1, nr_ao_blocks
     end do     ! a = 1, 3
   
   end subroutine


   subroutine get_s_lao(s, irep, mat_dim, dmat, ao, cx, cy, cz)
   
!    ----------------------------------------------------------------------------
     real(8), intent(out) :: s(3, 3)
     integer, intent(in)  :: irep(*)
     integer, intent(in)  :: mat_dim
     real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
     real(8), intent(in)  :: ao(*)
     real(8), intent(in)  :: cx, cy, cz
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
    
         do iblock = 1, nr_ao_blocks
    
           nr1 = ao_block_nr(iblock)
           st1 = ao_block_start(iblock)
           nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, irep(a)))
           st2 = ao_block_start(llss_block_partner(iblock, jxyz, irep(a)))
           if ((nr1 < 1) .or. (nr2 < 1)) cycle
    
           do ii = st1, st1 + nr1 - 1
             do jj = st2, st2 + nr2 - 1
    
               ff = (r_co(c)*ao(ao_off_g0_m1(0, b) + jj) - r_co(b)*ao(ao_off_g0_m1(0, c) + jj))*ao(ii) &
                  - (r_co(c)*ao(ao_off_g0_m1(0, b) + ii) - r_co(b)*ao(ao_off_g0_m1(0, c) + ii))*ao(jj)

               s(jxyz, a) = s(jxyz, a) &
                          + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff

    
             end do
           end do

         end do  ! do iblock = 1, nr_ao_blocks

       end do    ! do jxyz = 1, jxyz_start, 3
     end do      ! do a = 1, 3

   end subroutine

   subroutine get_gs_lao(s, gs, gsgn, gn_0, isym, mat_dim, dmat, ao, cx, cy, cz)

!    ----------------------------------------------------------------------------
     real(8), intent(out) :: s(3, 3)
     real(8), intent(out) :: gs(3, 3, 3)
     real(8), intent(out) :: gsgn(3, 3)
     real(8), intent(in)  :: gn_0(3)
     integer, intent(in)  :: isym(*)
     integer, intent(in)  :: mat_dim
     real(8), intent(in)  :: dmat(mat_dim, mat_dim, *)
     real(8), intent(in)  :: ao(*)
     real(8), intent(in)  :: cx, cy, cz
!    ----------------------------------------------------------------------------
     integer              :: nr1, nr2, st1, st2, iblock, jxyz, ixyz
     integer              :: a, b, c, ii, jj
     real(8)              :: r_co(3), ff, ff_gn
     real(8), external    :: ddot
     integer :: i, j, k
!    ----------------------------------------------------------------------------

     r_co(1) = cx
     r_co(2) = cy
     r_co(3) = cz

     s = 0.0d0
     gs = 0.0d0
     gsgn = 0.0d0
     ff= 0.0d0

     do a = 1, 3   ! B_{x,y,z}
       select case (a)
         case (1)
           b = 2; c = 3
         case (2)
           b = 3; c = 1
         case (3)
           b = 1; c = 2
       end select

       do jxyz = jxyz_start, 3  ! Sigma_{x,y,z}

         do iblock = 1, nr_ao_blocks

           nr1 = ao_block_nr(iblock)
           st1 = ao_block_start(iblock)
           nr2 = ao_block_nr(llss_block_partner(iblock, jxyz, isym(a) - 1))
           st2 = ao_block_start(llss_block_partner(iblock, jxyz, isym(a) - 1))
           if ((nr1 < 1) .or. (nr2 < 1)) cycle

           do ii = st1, st1 + nr1 - 1
             do jj = st2, st2 + nr2 - 1

!              lda term
               ff = (r_co(c)*ao(ao_off_g0_m1(0, b) + jj) - r_co(b)*ao(ao_off_g0_m1(0, c) + jj))*ao(ii) &
                  - (r_co(c)*ao(ao_off_g0_m1(0, b) + ii) - r_co(b)*ao(ao_off_g0_m1(0, c) + ii))*ao(jj)

               ff_gn = (gn_0(c)*ao(ao_off_g0_m1(0, b) + jj) - gn_0(b)*ao(ao_off_g0_m1(0, c) + jj))*ao(ii) &
                     - (gn_0(c)*ao(ao_off_g0_m1(0, b) + ii) - gn_0(b)*ao(ao_off_g0_m1(0, c) + ii))*ao(jj)

               s(jxyz, a) = s(jxyz, a) &
                          + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff

               gsgn(jxyz, a) = gsgn(jxyz, a) &
                             + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff_gn

!              start gradient terms
               ff_gn = 0.0d0
               do ixyz = 1, 3 ! dR_{x,y,z}

                  ff = (r_co(c)*ao(ao_off_g1_m1(ixyz, b) + jj) - r_co(b)*ao(ao_off_g1_m1(ixyz, c) + jj))*ao(ii)  &
                     - (r_co(c)*ao(ao_off_g1_m1(ixyz, b) + ii) - r_co(b)*ao(ao_off_g1_m1(ixyz, c) + ii))*ao(jj)  &
                     + (r_co(c)*ao(ao_off_g0_m1(0, b) + jj) - r_co(b)*ao(ao_off_g0_m1(0, c) + jj))*ao(ao_off_g1_m0(ixyz, 0) + ii)  &
                     - (r_co(c)*ao(ao_off_g0_m1(0, b) + ii) - r_co(b)*ao(ao_off_g0_m1(0, c) + ii))*ao(ao_off_g1_m0(ixyz, 0) + jj)

                  ff_gn = ff_gn + ff*gn_0(ixyz)

                  gs(ixyz, jxyz, a) = gs(ixyz, jxyz, a) &
                                    + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff

                  if (ixyz .eq. b) then

                     ff = -(ao(ao_off_g0_m1(0, c) + jj)*ao(ii) - ao(ao_off_g0_m1(0, c) + ii)*ao(jj))

                     gs(ixyz, jxyz, a) = gs(ixyz, jxyz, a) &
                                       + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff

                  else if (ixyz .eq. c) then

                     ff = ao(ao_off_g0_m1(0, b) + jj)*ao(ii) - ao(ao_off_g0_m1(0, b) + ii)*ao(jj)

                     gs(ixyz, jxyz, a) = gs(ixyz, jxyz, a) &
                                       + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff

                  end if

               end do ! do ixyz = 1, 3

               gsgn(jxyz, a) = gsgn(jxyz, a) &
                             + llss_prefactor(iblock)*dmat(jj, ii, iq(jxyz))*ff_gn

             end do
           end do

         end do  ! do iblock = 1, nr_ao_blocks

       end do    ! do jxyz = 1, jxyz_start, 3
     end do      ! do a = 1, 3

   end subroutine


   subroutine get_jblao(perturbation_component, cval, D_0, mat_dim, ao, r, j)

!     --------------------------------------------------------------------------
      character, intent(in)  :: perturbation_component
      real(8),   intent(in)  :: cval
      integer,   intent(in)  :: mat_dim
      real(8),   intent(in)  :: D_0(mat_dim, mat_dim, *)
      real(8),   intent(in)  :: ao(*)
      real(8),   intent(in)  :: r(3)
      real(8),   intent(out) :: j(3)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, k, l, b, c, ixyz
      real(8)                :: f
!     --------------------------------------------------------------------------

      select case (perturbation_component)
         case ('X')
            b = 2
            c = 3
         case ('Y')
            b = 3
            c = 1
         case ('Z')
            b = 1
            c = 2
      end select
    
      j = 0.0d0

      do ixyz = 1, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(lssl_block_partner(iblock, ixyz, 0))
            st2 = ao_block_start(lssl_block_partner(iblock, ixyz, 0))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
            do k = st1, st1 + nr1 - 1
               do l = st2, st2 + nr2 - 1
                  f = (r(c)*ao(ao_off_g0_m1(0, b) + l) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + l))*ao(k) &
                    - (r(c)*ao(ao_off_g0_m1(0, b) + k) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + k))*ao(l)
                  
                  j(ixyz) = j(ixyz) + f*D_0(l, k, iq(ixyz))
               end do
            end do
         end do
      end do

      j = -0.5d0*cval*j

   end subroutine

   subroutine get_gjblao(perturbation_component, cval, D_0, mat_dim, ao, r, gj)

!     --------------------------------------------------------------------------
      character, intent(in)  :: perturbation_component
      real(8),   intent(in)  :: cval
      integer,   intent(in)  :: mat_dim
      real(8),   intent(in)  :: D_0(mat_dim, mat_dim, *)
      real(8),   intent(in)  :: ao(*)
      real(8),   intent(in)  :: r(3)
      real(8),   intent(out) :: gj(3, 3)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, k, l, b, c, ixyz, jxyz
      real(8)                :: f
!     --------------------------------------------------------------------------

      select case (perturbation_component)
         case ('X')
            b = 2
            c = 3
         case ('Y')
            b = 3
            c = 1
         case ('Z')
            b = 1
            c = 2
      end select
    
      gj = 0.0d0

      do ixyz = 1, 3
         do iblock = 1, nr_ao_blocks
            nr1 = ao_block_nr(iblock)
            st1 = ao_block_start(iblock)
            nr2 = ao_block_nr(lssl_block_partner(iblock, ixyz, 0))
            st2 = ao_block_start(lssl_block_partner(iblock, ixyz, 0))
            if ((nr1 < 1) .or. (nr2 < 1)) cycle
            do jxyz = 1, 3
               do k = st1, st1 + nr1 - 1
                  do l = st2, st2 + nr2 - 1
                  f = (r(c)*ao(ao_off_g0_m1(0, b) + l) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + l))*ao(ao_off_g1_m0(jxyz, 0) + k) &
                    - (r(c)*ao(ao_off_g0_m1(0, b) + k) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + k))*ao(ao_off_g1_m0(jxyz, 0) + l) &
                    + (r(c)*ao(ao_off_g1_m1(jxyz, b) + l) &
                     - r(b)*ao(ao_off_g1_m1(jxyz, c) + l))*ao(k) &
                    - (r(c)*ao(ao_off_g1_m1(jxyz, b) + k) &
                     - r(b)*ao(ao_off_g1_m1(jxyz, c) + k))*ao(l)

                     gj(jxyz, ixyz) = gj(jxyz, ixyz) + f*D_0(l, k, iq(ixyz))

                     if (jxyz == b) then
                        f = - ao(ao_off_g0_m1(0, c) + l)*ao(k) &
                            + ao(ao_off_g0_m1(0, c) + k)*ao(l)
                        gj(jxyz, ixyz) = gj(jxyz, ixyz) + f*D_0(l, k, iq(ixyz))
                     end if
                     if (jxyz == c) then
                        f = + ao(ao_off_g0_m1(0, b) + l)*ao(k) &
                            - ao(ao_off_g0_m1(0, b) + k)*ao(l)
                        gj(jxyz, ixyz) = gj(jxyz, ixyz) + f*D_0(l, k, iq(ixyz))
                     end if
                  end do
               end do
            end do
         end do
      end do

      gj = -0.5d0*cval*gj

   end subroutine

   subroutine get_jblao_ll(perturbation_component, &
                           D_0,                    &
                           mat_dim,                &
                           ao,                     &
                           r,                      &
                           j)

!     --------------------------------------------------------------------------
      character, intent(in)  :: perturbation_component
      integer,   intent(in)  :: mat_dim
      real(8),   intent(in)  :: D_0(mat_dim, mat_dim, *)
      real(8),   intent(in)  :: ao(*)
      real(8),   intent(in)  :: r(3)
      real(8),   intent(out) :: j(3)
!     --------------------------------------------------------------------------
      integer                :: nr1, nr2, st1, st2, iblock, k, l, b, c, ixyz
      real(8)                :: f
!     --------------------------------------------------------------------------

      select case (perturbation_component)
         case ('X')
            b = 2
            c = 3
         case ('Y')
            b = 3
            c = 1
         case ('Z')
            b = 1
            c = 2
      end select
    
      j = 0.0d0

      do iblock = 1, nr_ao_blocks
         nr1 = ao_block_nr(iblock)
         st1 = ao_block_start(iblock)
         nr2 = ao_block_nr(llss_block_partner(iblock, 0, 0))
         st2 = ao_block_start(llss_block_partner(iblock, 0, 0))
         if ((nr1 < 1) .or. (nr2 < 1)) cycle
         do k = st1, st1 + nr1 - 1
            do l = st2, st2 + nr2 - 1
               do ixyz = 1, 3
                  f = (r(c)*ao(ao_off_g0_m1(0, b) + l) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + l))*ao(ao_off_g1_m0(ixyz, 0) + k) &
                    + (r(c)*ao(ao_off_g0_m1(0, b) + k) &
                     - r(b)*ao(ao_off_g0_m1(0, c) + k))*ao(ao_off_g1_m0(ixyz, 0) + l) &
                    - (r(c)*ao(ao_off_g1_m1(ixyz, b) + l) &
                     - r(b)*ao(ao_off_g1_m1(ixyz, c) + l))*ao(k) &
                    - (r(c)*ao(ao_off_g1_m1(ixyz, b) + k) &
                     - r(b)*ao(ao_off_g1_m1(ixyz, c) + k))*ao(l)

                  j(ixyz) = j(ixyz) + f*D_0(l, k, iq(0))
               end do
            end do
         end do
      end do

      j = -0.5d0*j

   end subroutine

end module
