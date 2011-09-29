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

module xc_blas3

   use interface_ao
   use xc_max_block_length
   use xc_derv
!fixme dft_cfg should not be used here, move xc derv business to separate module
   use dft_cfg
   use extra

   implicit none

   public xc_contribution_blas3
   public get_block_threshold

   private

contains

   subroutine xc_contribution_blas3(                         &
                                    is_gga,                  &
                                    do_lr,                   &
                                    block_length,            &
                                    block_threshold,         &
                                    mat_dim,                 &
                                    nz,                      &
                                    buffer,                  &
                                    ao_block,                &
                                    nr_fmat,                 &
                                    fmat,                    &
                                    dmat_0,                  &
                                    dmat_0_saop_outer,       &
                                    dmat_0_saop_resp,        &
                                    nr_dmat,                 &
                                    dmat_p,                  &
                                    dmat_pg_sym,             &
                                    dmat_ih_sym,             &
                                    derv_length,             &
                                    derv,                    &
                                    max_fun_derv,            &
                                    w,                       &
                                    nr_electrons_integrated, &
                                    xc_energy                &
                                   )

!     --------------------------------------------------------------------------
      logical, intent(in)    :: is_gga
      logical, intent(in)    :: do_lr
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: block_threshold
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      real(8)                :: buffer(*)
      real(8), intent(in)    :: ao_block(block_length, *)
      integer, intent(in)    :: nr_fmat
      real(8), intent(inout) :: fmat(mat_dim, mat_dim, nz, nr_fmat)
      real(8), intent(in)    :: dmat_0(mat_dim, mat_dim, nz)
      real(8), intent(in)    :: dmat_0_saop_outer(*)
      real(8), intent(in)    :: dmat_0_saop_resp(*)
      integer, intent(in)    :: nr_dmat
      real(8), intent(in)    :: dmat_p(mat_dim, mat_dim, nz, nr_dmat)
      integer, intent(in)    :: dmat_pg_sym(nr_dmat)
      integer, intent(in)    :: dmat_ih_sym(nr_dmat)
      integer                :: derv_length
      real(8)                :: derv(max_block_length, 0:derv_length)
      integer, intent(in)    :: max_fun_derv
      real(8)                :: w(*)
      real(8)                :: nr_electrons_integrated
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
      allocate(ao_sorted(max_block_length*mat_dim, 0:4))
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
                               xc_fun_alda,  &
                               xc_fun_xalda, &
                               max_fun_derv, &
                               block_length, &
                               w,            &
                               n(1, 0),      &
                               gnn(1, 0, 0), &
                               derv_length,  &
                               derv)

      do k = 1, block_length
         nr_electrons_integrated = nr_electrons_integrated + w(k)*n(k, 0)
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

!     --------------------------------------------------------------------------
      logical                :: is_gga
      integer                :: block_length
      integer                :: mat_dim
      integer                :: nz
      integer                :: ao_sort_index(*)
      real(8)                :: ao_sorted(max_block_length*mat_dim, 0:4)
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
      allocate(tmat(max_block_length*mat_dim))

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

!     --------------------------------------------------------------------------
      logical, intent(in)    :: do_lr
      logical, intent(in)    :: is_gga
      integer, intent(in)    :: block_length
      integer, intent(in)    :: mat_dim
      integer, intent(in)    :: nz
      integer, intent(in)    :: ao_sort_index(*)
      real(8), intent(in)    :: ao_sorted(max_block_length*mat_dim, 0:4)
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

      allocate(xmat(mat_dim*max_block_length))

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

   function get_block_threshold(         &
                                dfthr0,  &
                                mat_dim, &
                                dmat_0   &
                               )


!     --------------------------------------------------------------------------
      real(8)             :: get_block_threshold
!     --------------------------------------------------------------------------
      real(8), intent(in) :: dfthr0
      integer, intent(in) :: mat_dim
      real(8), intent(in) :: dmat_0(mat_dim, mat_dim)
!     --------------------------------------------------------------------------
      integer             :: i, j, iblock, st, nr
      real(8)             :: f, x
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
      get_block_threshold = dsqrt(10.0d0*dfthr0/(mat_dim*mat_dim*x))

   end function

end module
