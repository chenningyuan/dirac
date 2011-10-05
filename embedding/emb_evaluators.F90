

module fde_evaluators

   use fde_types
   use fde_data
   use fde_eval_dirac
   use fde_nadd_derv 

   public fde_get_density_matrix
   public fde_calculate_active_n_gn
   public fde_calculate_active_elpot
   public fde_cleanup_density_matrix
   public fde_calculate_embpot 

!   interface fde_get_density_matrix
!      module procedure fde_get_dmat_dirac
!   end interface

!   interface fde_calculate_active_n_gn
!      module procedure fde_get_ngn_dirac
!   end interface

!   interface fde_calculate_active_elpot
!      module procedure fde_get_vc_dirac
!   end interface

!   interface fde_cleanup_density_matrix
!      module procedure fde_free_dmat_dirac
!   end interface

   interface fde_calculate_embpot
      module procedure fde_calculate_vemb_dirac
   end interface

!  interface fde_calculate_embrsp
!     module procedure fde_calculate_lr_dirac
!  end interface

!  interface fde_calculate_embqrp
!     module procedure fde_calculate_qr_dirac
!  end interface


   contains
   
      subroutine calc_fde_interaction_energy(energy,gf_a,gf_b)
         type(grid_function), intent(in) :: gf_a
         type(grid_function), intent(in) :: gf_b
         real(kind=8),       intent(out) :: energy

         type(grid_function) :: gf_tot
         real(kind=8), allocatable :: derv_a(:,:)
         real(kind=8), allocatable :: derv_b(:,:)
         real(kind=8), allocatable :: derv_tot(:,:)
         real(kind=8) :: vca, vcb, w
         
         energy = 0.0d0

         call new_grid_function(gf_tot,gf_a%npoints)

         gf_tot%n       = gf_a%n  + gf_b%n
         gf_tot%gn      = gf_a%gn + gf_b%gn
         
         call fde_new_derv_result_partial(gf_a%npoints,derv_a)
         call fde_new_derv_result_partial(gf_b%npoints,derv_b)
         call fde_new_derv_result_partial(gf_tot%npoints,derv_tot)
                  
         call calc_xcke_energy_density(gf_a,  derv_a)
         call calc_xcke_energy_density(gf_b,  derv_b)
         call calc_xcke_energy_density(gf_tot,derv_tot)
         
         
         do i = 1, gf_a%npoints
            vca = gf_a%elpot(i)
            vcb = gf_b%elpot(i)
            w   = fde_grid(4,i)
            energy = energy                                          &
                   + w*(vcb + derv_tot(i,0) - derv_a(i,0))*gf_a%n(i) &
                   + w*(vca + derv_tot(i,0) - derv_b(i,0))*gf_b%n(i)
         enddo
         
         call fde_del_derv_result_partial(derv_a)
         call fde_del_derv_result_partial(derv_b)
         call fde_del_derv_result_partial(derv_tot)
         
      end subroutine calc_fde_interaction_energy

end module fde_evaluators


#ifdef FOOBAR

module fde_dens_eval

   use fde_data

   use interface_ao
   use interface_mo

   use dft_cfg
   use ao_eval
   use density_eval
   use xc_max_block_length

   implicit none

   public fde_get_act_density
   
   private
   
   real(8), allocatable :: rx(:)
   real(8), allocatable :: ry(:)
   real(8), allocatable :: rz(:)
   real(8), allocatable :: rw(:)
   
   integer :: npoints_used
   integer :: npoints_total
   integer :: nr_el_int
   integer :: max_ao_g_order
   integer :: max_ao_m_order
   integer :: nr_points_on_this_proc

   real(8), pointer             :: dmat_0(:, :)
   real(8), allocatable, target :: dmat_0_container(:, :)

   real(8), pointer             :: dmat(:, :, :, :)
   real(8), allocatable, target :: dmat_container(:, :, :, :)
   
   integer                      :: mat_0_dim
   
   contains
!   ----------------------------------------------------------------------------
      subroutine fde_get_act_density(dmat, mat_dim, get_gradient)
!   ----------------------------------------------------------------------------
         real(kind=8)  :: dmat(mat_dim,mat_dim)
         integer, intent(in)      :: mat_dim
         logical, intent(in) :: get_gradient
     
!         dmat_0 => dmat
!         mat_0_dim = mat_dim
         
         max_ao_m_order = 0
         max_ao_g_order = 0
         if (get_gradient) then
           max_ao_g_order = 1
         end if
     
         call ao_eval_init(max_ao_g_order, &
                           max_ao_m_order, &
                           .true.)

!        call insert_half_phases()
         call interface_mo_read()
         call loop_over_batches()

      end subroutine fde_get_act_density

!   ----------------------------------------------------------------------------
      subroutine loop_over_batches()
!   ----------------------------------------------------------------------------
         real(8) :: dummy
         integer :: idummy
         integer :: ibatch
         integer :: start_batch, i
         integer :: npoints_batch

         ibatch        = 0
         npoints_total = 0
         npoints_used  = 0
         
         start_batch = 1
         npoints_batch = fde_grid_np

!        loop over batches
         do
            if (npoints_batch < 0) exit

            allocate(rx(npoints_batch))
            allocate(ry(npoints_batch))
            allocate(rz(npoints_batch))
            allocate(rw(npoints_batch))

!              copy data over to rx, ry, rz, rw, fde_froz_n, fde_froz_ng
            call get_fde_grid_coord(start_batch,npoints_batch,rx,ry, rz)
            call get_fde_grid_weight(start_batch,npoints_batch,rw)
                 
            ibatch = ibatch + 1
            npoints_total = npoints_total + npoints_batch
            nr_points_on_this_proc = npoints_batch

            call integrate()

            deallocate(rx)
            deallocate(ry)
            deallocate(rz)
            deallocate(rw)

!          don't try to read another batch
           goto 999

!     end loop over batches
         end do
 999     continue
      end subroutine loop_over_batches

! ------------------------------------------------------------------------------
      subroutine integrate()
! ------------------------------------------------------------------------------
!        unperturbed density
         real(8) :: n_0(max_block_length)
!        unperturbed density gradient
         real(8) :: gn_0_pointwise(3)

         real(8) :: w(max_block_length)
         integer :: ipoint
         integer :: i
         real(8), allocatable :: ao(:, :)
         real(8), allocatable :: buffer(:, :)

         call alloc(ao,     1, nr_ao_slices*nr_ao_cartesian)
         call alloc(buffer, 1, nr_ao_slices*nr_ao_cartesian)
         ao     = 0.0d0
         buffer = 0.0d0
   
         do ipoint = 1, nr_points_on_this_proc
         
            gn_0_pointwise = 0.0d0

            call get_ao(1,          &
                        rx(ipoint), &
                        ry(ipoint), &
                        rz(ipoint), &
                        ao,         &
                        buffer)

            if (need_ao_order(1, 0)) then
!              density and density gradient
               call get_gn(n_0(1), gn_0_pointwise, 0, mat_0_dim, dmat_0, buffer, ao)
            else
!              density
               call get_n(n_0(1), 0, mat_0_dim, dmat_0, buffer, ao)
            endif

            if (n_0(1) > dft_cfg_small_density_threshold) then
               npoints_used = npoints_used + 1

               w(1) = rw(ipoint)

               fde_active_n(ipoint)    = n_0(1)
               fde_active_gn(ipoint,1) = gn_0_pointwise(1)
               fde_active_gn(ipoint,2) = gn_0_pointwise(2)
               fde_active_gn(ipoint,3) = gn_0_pointwise(3)

!       number of electrons
               nr_el_int = nr_el_int + w(1)*n_0(1)
             endif
          enddo
          deallocate(ao)
          deallocate(buffer)
   
      end subroutine

end module fde_dens_eval  

#endif
