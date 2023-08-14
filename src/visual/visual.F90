module visual

  use visual_cfg
  use visual_utils
  use visual_dirac
  use visual_in_point
  use visual_dirac
  use visual_grids
  use visual_grid_functions
  use visual_datatypes
  use visual_io
  use visual_mpi
#ifdef VAR_MPI
  use interface_to_mpi
  use integer_kind_mpilib
#endif
  use checkpoint

!   ---- legacy ----    
  use memory_allocator
  use dircmo
  use dirac_ao_eval
  use interface_ao
  use matrix_defop_old
  use matrix_genop_old, only: copy_matrix_sym
  use dirac_interface
  use visual_london
  use xc_london_c1
  use visual_integrate
  use electrostatic_potential
  use num_grid_cfg

  use labeled_storage
  use visual_txt_storage
  use visual_dirac_input
!   ---- end ----    

  implicit none

!   ---- legacy ----    
#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbxlr.h"
#include "codata.h"
#include "orgcom.h"

    integer, parameter:: file_unit = 137

!   ---- end ----    


  public run_visual


  private

contains
  subroutine run_visual(launch_only)
    logical, intent(in), optional :: launch_only

    type(mol_t):: mol
    type(ao_t):: ao_data
    type(mo_t):: mo_data
    real(8) :: walltime
    logical :: legacy_prep

!   ---- legacy ----    
    type(matrix)  :: D_node
    !type(matrix), pointer  :: D_node
    real(8), allocatable:: ao_node(:)
    real(8), allocatable:: buffer_node(:)
    real(8), allocatable, target:: D_container(:)


    type(matrix)         :: D_mo
    type(matrix)         :: D(max_para)
    !type(matrix),target         :: D(max_para)
    type(matrix)         :: D_0
    type(matrix)         :: D_con_mo, D_con_mo_ai, D_con, D_test
    type(matrix)         :: C, Cig, Ciu, Csg, Csu

    real(8), allocatable:: eigenvalues(:)
    real(8), allocatable:: occupation(:)
    real(8), allocatable:: ao(:)
    real(8), allocatable:: buffer(:)
    real(8), allocatable:: cmo_from_file(:)
    real(8), allocatable:: occ_from_file(:)

    real(8)              :: energy, time_start, second
    integer              :: max_ao_g_order, max_ao_m_order
    integer              :: ioff, imat, iq, iz, m, ivector
    integer, parameter   :: file_unit = 137

    logical              :: debug_me = .false.
    logical              :: krmcscf_found
    logical              :: natorb_found
    logical, external    :: fndlab

    integer              :: q_max_dim
    integer              :: irank

    integer:: i, j, off_ij, off_iz, off_ij_s, off_ij_e, jrdmo
    type(matrix):: D_mo_dagger, TP, TM, T
    logical, allocatable:: matrix_is_unperturbed(:)
!   ---- end ----    

    walltime = visual_walltime()


#ifdef VAR_MPI
    ! We should have initial information from user input stored in rs_input dt
    ! and other config information in visual_cfg* variables (legacy) available at this point.
    ! So first, we make sure that all workers have this information to start.
    ! Also, run_visual() is called from master (from subroutine dirac) 
    ! and from workers (from subroutine dirnod);
    ! in the latter case - we only want to wake up nodes, 
    ! pass initial information, and leave the rest to master - this is handled with "launch+only" flag
    call visual_initialize_and_sync_nodes(launch_only)
    print *, 'VISUAL parallel: VAR_MPI', is_parallel()
    print *, 'size of communicator', mpi_get_nr_proc()
    print *, 'rank of the process on master ', get_my_rank()
    print *, 'will call visual_initialize_nodes from visual'
#endif


    !legacy_prep = .true.
    !if (legacy_prep) then


    if (mpi_is_master()) then !temp
    max_ao_g_order = 0
    do imat = 1, visual_cfg_nr_dmat
       if (                                                        &
           (visual_cfg_property(imat) == iq_elf)              .or. &
           (visual_cfg_property(imat) == iq_kin)              .or. &
           (visual_cfg_property(imat) == iq_kin_ls)           .or. &
           (visual_cfg_property(imat) == iq_kin_sl)           .or. &
           (visual_cfg_property(imat) == iq_kin_tau)          .or. &
           (visual_cfg_property(imat) == iq_divs)             .or. &
           (visual_cfg_property(imat) == iq_rots)             .or. &
           (visual_cfg_property(imat) == iq_divj)             .or. &
           (visual_cfg_property(imat) == iq_rotj)             .or. &
           (visual_cfg_property(imat) == iq_density_gradient) .or. &
           (visual_cfg_property(imat) == iq_rdg)              .or. &
           (visual_cfg_property(imat) == iq_gradj)                 &
          ) then
          if (max_ao_g_order < 1) then
             max_ao_g_order = 1
          end if
       end if
       if (                                                       &
           (visual_cfg_property(imat) == iq_kin_lap)          .or. &
           (visual_cfg_property(imat) == iq_density_hessian)  .or. &
           (visual_cfg_property(imat) == iq_density_laplacian).or. &
           (visual_cfg_property(imat) == iq_density_sign_lambda2).or. &
           (visual_cfg_property(imat) == iq_kin_nr)               &
          ) then
          if (max_ao_g_order < 2) then
             max_ao_g_order = 2
          end if
       end if
       if (levyle) then
!         here we obtain the paramagnetic part of the current density
!         via the nr paramagnetic current density expression
!         for this we need the cartesian gradient
          if (                                             &
              (visual_cfg_property(imat) == iq_j)     .or. &
              (visual_cfg_property(imat) == iq_ndipx) .or. &
              (visual_cfg_property(imat) == iq_ndipy) .or. &
              (visual_cfg_property(imat) == iq_ndipz) .or. &
              (visual_cfg_property(imat) == iq_bdipx) .or. &
              (visual_cfg_property(imat) == iq_bdipy) .or. &
              (visual_cfg_property(imat) == iq_bdipz)      &
             ) then
             if (max_ao_g_order < 1) then
                max_ao_g_order = 1
             end if
          end if
       end if
    end do

    max_ao_m_order = 0
    if (visual_cfg_london) then
       max_ao_m_order = 1
    end if

    write(LUPRI, '(/,2x, a)') '*** Visualisation module ***'
    call verify_input_keywords()

!   get coefficients

    allocate(cmo_from_file(n2bbasxq))
    cmo_from_file = 0.0d0
    allocate(occ_from_file(norbt))
    occ_from_file = 0.0d0

    inquire(file='KRMCSCF', exist = krmcscf_found)

    if (krmcscf_found) then
!      read from KRMCSCF

       call opnfil(file_unit, 'KRMCSCF', 'OLD', 'visual.F90')
       rewind(file_unit)

       natorb_found = .FALSE.
       natorb_found = fndlab('MCCINATO', file_unit)
       rewind(file_unit)
       jrdmo = -1

       IF (natorb_found) THEN
          call rreadmo(cmo_from_file, jrdmo, 4, file_unit)
          call reakrmc(file_unit, 'MCNATOCC', occ_from_file, norbt)
          if (jrdmo == 0) then
             print *, 'orbitals read from label mccinato on file krmcscf'
          end if
       else
          call rreadmo(cmo_from_file, jrdmo, 1, file_unit)
          if (jrdmo == 0) then
             print *, 'orbitals read from label neworb on file krmcscf'
          end if
       end if

       close(file_unit)
    else
!      read from CHECKPOINT FILE

       call alloc(eigenvalues, norbt)
       eigenvalues = 0.0d0
       call reacmo_new(cmo = cmo_from_file, eig = eigenvalues, toterg = energy)

       write(LUPRI, '(2x, a)') 'Visualisation : CHECKPOINT file read'

!fixme verify whether this is really needed
!      when i introduced it i needed it so that certain arrays were properly set
!      when restarting from DFCOEF
       call fndocc(eigenvalues, .true.)
       call dealloc(eigenvalues)
    end if

!   get orbital string
    call alloc(occupation, norbt)
    call setup_orbital_string(occupation)

!warning: when reading from KRMCSCF many assumptions below are either incorrect
!         or become a hack (occupation 2, only inactive and secondary)
!         it's ok if you plot a single orbital but not ok for many things
!         need to add careful checks

    call get_c(C,   cmo_from_file, i = 1.0d0, s = 1.0d0, g = 1.0d0, u = 1.0d0)
    call get_c(Cig, cmo_from_file, i = 1.0d0, s = 0.0d0, g = 1.0d0, u = 0.0d0)
    call get_c(Csg, cmo_from_file, i = 0.0d0, s = 1.0d0, g = 1.0d0, u = 0.0d0)
    if (debug_me) then
       call print_mat(C,   label='debug C')
       call print_mat(Cig, label='debug Cig')
       call print_mat(Csg, label='debug Csg')
    end if

    if (nfsym == 2) then
       call get_c(Ciu, cmo_from_file, i = 1.0d0, s = 0.0d0, g = 0.0d0, u = 1.0d0)
       call get_c(Csu, cmo_from_file, i = 0.0d0, s = 1.0d0, g = 0.0d0, u = 1.0d0)
       if (debug_me) then
          call print_mat(Ciu, label='debug Ciu')
          call print_mat(Csu, label='debug Csu')
       end if
    end if
    deallocate(cmo_from_file)
    deallocate(occ_from_file)

    call init_mat(D_mo, norbt, norbt)

    if (visual_cfg_london) then

       D_mo%elms = 0.0d0
       do m = 1, norbt
          D_mo%elms((m-1)*norbt+m) = occupation(m)
       end do
       D(imat) = C*(D_mo*dag(C))

       D_0 = C*(D_mo*dag(C))
       D_0 = 2.0d0*D_0

       if (.not. visual_cfg_london_skip_ro) then
          ! get T^B matrix in MO basis

          call init_mat(D_con_mo, norbt, norbt)
          D_con_mo%elms = 0.0d0

          if (visual_cfg_london_none) then

             call init_mat(D_con_mo_ai, norbt, norbt)
             D_con_mo_ai%elms = 0.0d0

             call get_d_con_mo(D_con_mo, occupation, &
                               D_con_mo_ai = D_con_mo_ai)

             call copy_matrix_sym(D_con_mo_ai, D_con_mo)

          else

             call get_d_con_mo(D_con_mo, occupation)

          end if

!gosia:
! comparison with getmdm shows that this factor is not needed for 4cDC Hamiltonian; yet this works...
          if (.not. levyle) then
!            radovan: magic factor for levyle
!                     need to find the reason for it
             D_con_mo = 0.5d0*D_con_mo
             if (visual_cfg_london_none) D_con_mo_ai = 0.5d0*D_con_mo_ai
          end if

!         build 'modified' (by T^B) density matrices:

          call init_mat(D_con, norbt, norbt)
          D_con%elms = 0.0d0
          call copy_matrix_sym(D_con, D_con_mo)

          if (nfsym == 2) then
             call get_mdm_con(D_con,                     &
                              D_con_mo,                  &
                              C, Cig, Csg,               &
                              Ciu = Ciu,                 &
                              Csu = Csu)
          else
             call get_mdm_con(D_con,                     &
                              D_con_mo,                  &
                              C, Cig, Csg)
          end if
       end if

    end if

    allocate(matrix_is_unperturbed(visual_cfg_nr_dmat))
    matrix_is_unperturbed = .false.

    do imat = 1, visual_cfg_nr_dmat

       if (visual_cfg_dmat_file(imat) == 'DFCOEF') then
          matrix_is_unperturbed(imat) = .true.
!         get unperturbed density matrix

          D_mo%elms = 0.0d0
          do m = 1, norbt
             D_mo%elms((m-1)*norbt+m) = occupation(m)
          end do

          D_mo%ih_sym = 1
          D_mo%tr_sym = 1
          D_mo%irep   = 0

          D(imat) = C*(D_mo*dag(C))
          if (debug_me) then
             call print_mat(D_mo,    label = 'debug MO D matrix')
             call print_mat(D(imat), label = 'debug AO D matrix')
          end if

       else
!         get perturbed density matrix

          call setup_fo_dmat_ao(imat, occupation, D_mo)

          if (visual_cfg_london .and. visual_cfg_london_none) then
             TP = D_mo-dag(D_con_mo_ai)
             TM = D_mo-D_con_mo_ai
          else
             TP = 1.0d0*D_mo
             TM = 1.0d0*D_mo
          end if

          if (nfsym == 2) then
             if (jbtof(D_mo%irep, 1) == 2) then
!               ungerade perturbation
                D(imat) = (Cig*(TP*dag(Csu))) &
                        + (Ciu*(TP*dag(Csg))) &
                        - (Csg*(TM*dag(Ciu))) &
                        - (Csu*(TM*dag(Cig)))
             else
!               gerade perturbation
                D(imat) = (Cig*(TP*dag(Csg))) &
                        + (Ciu*(TP*dag(Csu))) &
                        - (Csg*(TM*dag(Cig))) &
                        - (Csu*(TM*dag(Ciu)))
             end if
          else
             D(imat) = (Cig*(TP*dag(Csg))) &
                     - (Csg*(TM*dag(Cig)))
          end if
          call copy_matrix_sym(D(imat), D_mo)
       end if
    end do

!   insert half phases
    if (nz < 4) then
       do iz = 1, nz
         ioff = ntbas(0)*ntbas(0)*(iz-1)
         if (visual_cfg_london) then
            iq   = ipqtoq(iz, D_0%irep)
            call q2bphase('D', iq, 1, D_0%elms(1+ioff))
            iq   = ipqtoq(iz, D_con%irep)
            call q2bphase('D', iq, 1, D_con%elms(1+ioff))
         end if
         do imat = 1, visual_cfg_nr_dmat
            iq   = ipqtoq(iz, D(imat)%irep)
            call q2bphase('D', iq, 1, D(imat)%elms(1+ioff))
         end do
       end do
    end if

    do imat = 1, visual_cfg_nr_dmat
       if (visual_cfg_london) then
          if (visual_cfg_london_skip_kappa) then
!            zero out kappa contribution
!            this is useful to isolate direct and/or reortho contributions
             D(imat) = tiny(0.0d0)*D(imat)
          end if

          if (.not. visual_cfg_london_skip_ro) then
!            dirty trick add D_con on top of D
!            on trunk we will have to check that this only
!            happes to the current, mag field, and so on
             if (.not. matrix_is_unperturbed(imat)) then
                ! don't do it with unp matrices, otherwise double counting
                D(imat) = D(imat) + D_con
             end if
          end if

       end if

!      scale density matrix
       D(imat) = 2.0d0*D(imat)

!      collect matrices for identical densities to save time
       if (visual_cfg_property(imat) /= iq_elf) then
       !if (visual_cfg_property(imat) /= iq_elf .and. visual_cfg_property(imat) /= iq_gradj) then
       ! gosia TODO: i need to figure out a better way to sum up contribs involving the same densities
       ! for a moment it's messing up the "GRADJ" calculations if done in the same input as "J" calculations

!         elf has nonlinear dependence
          j = visual_cfg_property_to_matrix(visual_cfg_property(imat))
          if (j == 0) then
             visual_cfg_property_to_matrix(visual_cfg_property(imat)) = imat
          else
             D(j) = D(j) + D(imat)
             if (visual_cfg_london) then
                T   = 1.0d0*D_0
                D_0 = D_0+T
                T   = 0
             end if
             visual_cfg_skip(imat) = .true.
          end if
       end if
    end do

!   deallocate various things
!   =========================

    deallocate(matrix_is_unperturbed)
    call dealloc(occupation)

    C    = 0
    Cig  = 0
    Csg  = 0

    if (nfsym == 2) then
      Ciu = 0
      Csu = 0
    end if


!   initialize ao module
!   ====================

    call dirac_ao_eval_init(max_ao_g_order, &
                      max_ao_m_order, &
                      .false.)

    call alloc(ao,     nr_ao_slices*nr_ao_cartesian*max_batch_length)
    call alloc(buffer, nr_ao_slices*nr_ao_cartesian*max_batch_length)

    if (visual_cfg_london) then
       if (.not. visual_cfg_london_skip_ro) then
          D_con_mo = 0
          D_con    = 0
          if (visual_cfg_london_none) D_con_mo_ai = 0
       end if
    end if


!   plot and integrate
!   ==================
    if (visual_cfg_tensor) then
      q_max_dim = 12
    else
      q_max_dim = 3
    end if

    if (visual_cfg_nics) then
      write(*,*)"VISUAL: nics value in ", visual_cfg_nics_origin(1), &
                                          visual_cfg_nics_origin(2), &
                                          visual_cfg_nics_origin(3)
    end if
    time_start = second()

    end if ! if (mpi_is_master()) then

    !end if ! if (legacy_prep) then

    if (visual_cfg_global_refactor) then

    ! gosia 19/5: now we do calculations on one molecule at a time, 
    ! so here we read the molecular data from DIRAC's main checkpoint;
    ! but this can be moved inside the loop over grids

!   todo - when developing a separate schema
!    write (lupri,'(/A,I2)')," * Opening checkpoint file for real-space data (VISUAL) "
!    call checkpoint_rs_open


    !if (is_parallel()) then
    !  if (mpi_is_master()) walltime = interface_MPI_WTIME()
    !else
    !  call cpu_time(walltime)
    !end if

    ! gosia: now this only reads mol data
    if (mpi_is_master()) call visual_initialize_data_from_dirac(mol, ao_data, mo_data)

    !temp:
! here: temp: pass to workers
#ifdef VAR_MPI
    if (is_parallel()) then
      call interface_mpi_bcast(nr_ao_slices,1,irank_master,global_communicator)
      call interface_mpi_bcast(nr_ao_cartesian,1,irank_master,global_communicator)
      !todo: track max_batch_length
      call interface_mpi_bcast(visual_cfg_nr_dmat,1,irank_master,global_communicator)
    end if
#endif

#ifdef VAR_MPI
    if (is_parallel()) then
      if (.not. mpi_is_master()) then
        !call init_mat(D_node, reset=.true.)
        allocate(D_container(n2bbasxq_node))
        D_node%elms => D_container
        write(*,*) 'D_node init on node: ', get_my_rank(), size(D_container), size(D_node%elms)
      else
        D_node = 2.0d0*D(1)
        D_node%irep = 0
        D_node%nrow = ntbas(0)
        D_node%ncol = ntbas(0)
        D_node%nz   = nz
        write(*,*) 'D_node init on master: ', get_my_rank(), size(D_node%elms)
      end if
    else
      D_node = 2.0d0*D(1)
      D_node%irep = 0
      D_node%nrow = ntbas(0)
      D_node%ncol = ntbas(0)
      D_node%nz   = nz
    end if
#else
    D_node = 2.0d0*D(1)
    D_node%irep = 0
    D_node%nrow = ntbas(0)
    D_node%ncol = ntbas(0)
    D_node%nz   = nz
#endif

    allocate(ao_node(nr_ao_slices*nr_ao_cartesian))
    ao_node = 0.0d0
    !write(*, *) 'allocated ao on node, size: ', get_my_rank(), nr_ao_slices*nr_ao_cartesian

    allocate(buffer_node(nr_ao_slices*nr_ao_cartesian))
    buffer_node = 0.0d0
    !write(*, *) 'allocated buffer on node, size: ', get_my_rank(), nr_ao_slices*nr_ao_cartesian

#ifdef VAR_MPI
    if (is_parallel()) then
      call interface_mpi_bcast(D_node%irep,1,irank_master,global_communicator)
      call interface_mpi_bcast(D_node%nrow,1,irank_master,global_communicator)
      call interface_mpi_bcast(D_node%ncol,1,irank_master,global_communicator)
      call interface_mpi_bcast(D_node%nz,1,irank_master,global_communicator)
      call interface_mpi_bcast(D_node%elms,size(D_node%elms),irank_master,global_communicator)
    end if
#endif
! here: end temp: pass to workers
    !if (get_my_rank() == 0) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_node%elms(i),i=1,n2bbasxq_node)
    !if (get_my_rank() == 0) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_node%elms(1:3))
    !write(*, *) 'TESTME D_node%irep ', get_my_rank(), D_node%irep
    !write(*, *) 'TESTME D_node%nrow ', get_my_rank(), D_node%nrow
    !write(*, *) 'TESTME D_node%ncol ', get_my_rank(), D_node%ncol
    !write(*, *) 'TESTME D_node%nz   ',   get_my_rank(), D_node%nz
    !if (get_my_rank() == 7) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_node%elms(i),i=1,n2bbasxq_node)
    !if (get_my_rank() == 7) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_node%elms(1:3))


#ifdef VAR_MPI
      if (is_parallel()) then
        call synchronize_data_from_dirac(mol, ao_data)
      end if
#endif
! 
    call visual_driver_refactored(visual_cfg_nr_dmat, &
                       D,                  &
                       D_node,                &
                       mol, ao_data,       &
                       ao_node,                 &
                       buffer_node)

!    call checkpoint_rs_close

    !if (mpi_is_master()) then
    !  walltime = interface_MPI_Wtime() - walltime
    !  write (*, '(a)') ''
    !  write (*, '(a, g14.6, a)') &
    !        'Elapsed wall clock time in VISUAL = ', walltime, ' seconds.'
    !  write (*, '(a)') ''
    !end if

    if (associated(D_node%elms)) then
      nullify(D_node%elms)
      !D_node = 0
    end if

    if (allocated(ao_node)) deallocate(ao_node)
    if (allocated(buffer_node)) deallocate(buffer_node)

    if (.not. mpi_is_master()) then
      if (allocated(D_container)) deallocate(D_container)
    end if


    else
    call visual_driver(visual_cfg_nr_dmat, &
                       D,                  &
                       D_0,                &
                       ao,                 &
                       buffer)
    end if

    if (allocated(rs_input%grids)) deallocate(rs_input%grids)
    if (allocated(rs_input%grid_functions)) deallocate(rs_input%grid_functions)
    call rs_input%erase()

    if (mol%in_use()) then
      call mol%erase()
    end if

! ---- legacy ----    
    if (visual_cfg_2d_integration .or. visual_cfg_3d_integration) then
       call visual_integration(visual_cfg_nr_dmat, &
                               D,                  &
                               D_0,                &
                               ao,                 &
                               buffer, &
                               q_max_dim)
    end if
    call timtxt('>>>> time used in integration:', (second() - time_start), lupri)

    if (mpi_is_master()) then !temp
    call dealloc(ao)
    call dealloc(buffer)


!   deallocate density matrices

    do imat = 1, visual_cfg_nr_dmat
       D(imat) = 0
    end do
    end if !if (mpi_is_master()) then !temp


! ---- legacy end ----    

    walltime = walltime - visual_walltime()
    if (mpi_is_master()) print *, 'Elapsed walltime in run_visual: ', walltime

  end subroutine


  subroutine visual_driver_refactored(nr_dmat, &
                            D,       &
                            D_0,     &
                            mol, ao_data, &
                            ao,      &
                            buffer)

!   --------------------------------------------------------------------------
    integer, intent(in)  :: nr_dmat
    type(matrix)         :: D(nr_dmat)
    type(matrix)         :: D_0
    real(8), intent(in)  :: ao(*)
    real(8)              :: buffer(*)
    type(mol_t), intent(inout):: mol
    type(ao_t), intent(inout):: ao_data
!   --------------------------------------------------------------------------
    type(grid_t):: grid
    type(grid_t), allocatable:: grids(:)
    type(grid_function_t), allocatable:: grid_functions(:)
    real(8), allocatable:: grid_function_values(:,:)
    integer:: ig, igf
    integer:: g_ndim, gf_ndim, npoints
    real(8):: extents(0:1, 3)
    integer:: my_rank
    real(8) :: walltime
    integer:: n_this_proc, idim_grid, idim_gf
    integer:: i,j
!   --------------------------------------------------------------------------

    if (mpi_is_master()) call rs_input%display

    do ig = 1, rs_input%nr_grids

!
!     get necessary information about the grid and share with workers
!

      if (mpi_is_master()) then
        call get_grid_info(ig, rs_input%grids(ig), mol)
      end if
#ifdef VAR_MPI
      if (is_parallel()) call synchronize_grid_data(rs_input%grids(ig))
#endif

      g_ndim  = rs_input%grids(ig)%ndim
      npoints = rs_input%grids(ig)%nr_points_total

      allocate(rs_input%grids(ig)%points(g_ndim, npoints))
      rs_input%grids(ig)%points = 0.0d0

!
!     get grid points
!
!     gosia:
!     * TODO: for optimization - try making all arrays 1-dimensional, then reshape?
!     grid_points = reshape(rs_input%grids(ig)%points(g_ndim, npoints))
!     * FIXME: grid can be read from file or generated;
!       * reading from file should be done on master
!       * grid generation should be done in parallel
!       now I do everything on master

      if (mpi_is_master()) then
        call get_grid_points(rs_input%grids(ig))
        if (rs_input%print_level > 0) then
          write(*, *) 'Grid points of grid ', ig
          do j = 1, npoints
            write(*, *) (rs_input%grids(ig)%points(i,j),i=1,g_ndim)
          end do
        end if
      end if

      if (rs_input%grids(ig)%export) then
!       also export grid to a separate file 
        if (mpi_is_master()) then
          call store_grid(rs_input%grids(ig))
          if (rs_input%print_level > 0) then
            write(*, *) 'Grid ', rs_input%grids(ig)%id, ' exported to file ', rs_input%grids(ig)%file_out
          end if
        end if
      end if

      do igf = 1, rs_input%nr_grid_functions

!
!       get necessary information about the grid function and share with workers
!
!       gosia:
!       * FIXME: now, this is used only to assign labels for hdf5 storage,
!         TODO: is there sth that could be shared now?

        if (mpi_is_master()) then 
          call get_grid_function_info(igf, rs_input%grid_functions(igf))
        end if
#ifdef VAR_MPI
        if (is_parallel()) call synchronize_grid_function_data(rs_input%grid_functions(igf))
#endif

        if (rs_input%grid_functions(igf)%id_assigned_grid == rs_input%grids(ig)%id) then

          ! gosia todo: it would be cleaner to use "associate", 
          ! but some compilers have troubles with nested associate
          ! https://stackoverflow.com/questions/66562058/associate-chains-in-fortran

          if (rs_input%print_level > 0) then
            if (mpi_is_master()) then 
              write(*, *) 'Start working on grid function ', igf, ' on grid ', ig
              call rs_input%grids(ig)%display
              call rs_input%grid_functions(igf)%display
              write(*, *) 'Grid points to work on: '
              do j = 1, npoints
                write(*, *) (rs_input%grids(ig)%points(i,j), i=1,g_ndim)
              end do
            end if
          end if
!
!         distribute grid points over nodes
!

          n_this_proc = rs_input%grids(ig)%nr_points_total
#ifdef VAR_MPI
          if (is_parallel()) then
            ! redefine n_this_proc
            !  call visual_interface_mpi_distribute_points(&
            !             rs_input%grids(ig)%points(1,:), &
            !             rs_input%grids(ig)%points(2,:), &
            !             rs_input%grids(ig)%points(3,:), &
            !             rs_input%grids(ig)%nr_points_total, n_this_proc)
            do idim_grid = 1, rs_input%grids(ig)%ndim
              call visual_interface_mpi_distribute_points(&
                         rs_input%grids(ig)%points(idim_grid,:), &
                         rs_input%grids(ig)%nr_points_total, n_this_proc)
            end do

            if (get_my_rank() == 7) then
              do j =1, n_this_proc
                write(*, '(a,i,3f8.3)') 'points on node 7: ', j, &
                     (rs_input%grids(ig)%points(i,j), i=1,g_ndim)
              end do
            end if
          end if
#endif

!
!         do work on grid
!

          gf_ndim  = rs_input%grid_functions(igf)%get_ndim()

          allocate(grid_function_values(gf_ndim, n_this_proc))
          grid_function_values = 0.0d0

          ! gosia todo:
          ! add "visualize"/"integrate" procedures in this call; now - visualization only
          ! add constructs for DIRAC data (D, D_0, ...)
          !if (get_my_rank() == 0) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_0%elms(1:3))
          !if (get_my_rank() == 7) write(*, *) 'TESTME D_node%elms   ',   get_my_rank(), (D_0%elms(1:3))
          call work_on_grid(rs_input%grids(ig), rs_input%grid_functions(igf), &
                            n_this_proc, &
                            grid_function_values, &
                            nr_dmat, D, D_0, ao, buffer)

          allocate(rs_input%grid_functions(igf)%values(gf_ndim, rs_input%grids(ig)%nr_points_total))
          rs_input%grid_functions(igf)%values = 0.0d0

            !write(*, *) 'AFTER WORK_ON_GRID'
            !write(*, *) 'values: ', ((grid_function_values(i,j),j=1,n_this_proc),i=1,gf_ndim), get_my_rank()
 
          if (is_parallel()) then
            do idim_gf = 1, gf_ndim
              call visual_interface_mpi_gatherv(grid_function_values(idim_gf,:),&
                                                rs_input%grid_functions(igf)%values(idim_gf,:),&
                                                n_this_proc)
            end do
          else
            do idim_gf = 1, gf_ndim
              rs_input%grid_functions(igf)%values(idim_gf,:) = grid_function_values(idim_gf,:)
            end do
          end if

          if (mpi_is_master()) then
            write(*, *) 'AFTER GATHERING'
            write(*, *) 'gathered values: ', rs_input%grid_functions(igf)%values
          end if
 
          deallocate(grid_function_values)

!
!         save data
!

          if (rs_input%grid_functions(igf)%export) then
             if (mpi_is_master()) then
               call store_data(rs_input%grids(ig), rs_input%grid_functions(igf))
             end if
          end if

          if (mpi_is_master()) then
            write(*, *) 'AFTER STORING'
            write(*, *) 'stored values?: ', rs_input%grid_functions(igf)%values
          end if
 
          if (allocated(rs_input%grid_functions(igf)%values))   deallocate(rs_input%grid_functions(igf)%values)
        end if ! if (rs_input%grid_functions(igf)%id_assigned_grid == ig) then

        if (.not. mpi_is_master()) then
        if (allocated(rs_input%grid_functions(igf)%name)) deallocate(rs_input%grid_functions(igf)%name)
        if (allocated(rs_input%grid_functions(igf)%file_inp)) deallocate(rs_input%grid_functions(igf)%file_inp)
        if (allocated(rs_input%grid_functions(igf)%file_out)) deallocate(rs_input%grid_functions(igf)%file_out)
        end if

      end do ! do igf = 1, rs_input%nr_grid_functions

      if (allocated(rs_input%grids(ig)%points)) deallocate(rs_input%grids(ig)%points)
      if (allocated(rs_input%grids(ig)%file_inp)) deallocate(rs_input%grids(ig)%file_inp)
      if (allocated(rs_input%grids(ig)%file_out)) deallocate(rs_input%grids(ig)%file_out)

    end do ! do ig = 1, rs_input%nr_grids


    do igf = 1, rs_input%nr_grid_functions
      call rs_input%grid_functions(igf)%erase()
    end do
    do ig = 1, rs_input%nr_grids
      call rs_input%grids(ig)%erase()
    end do


  end subroutine visual_driver_refactored


  subroutine work_on_grid(grid, grid_function, &
                          nr_points, &
                          grid_function_values, &
                          nr_dmat, D, D_0, ao, buffer)

    type(grid_t), intent(in):: grid
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out)  :: grid_function_values(:,:)
    integer, intent(in)  :: nr_points
    integer, intent(in)  :: nr_dmat
    type(matrix)         :: D(nr_dmat)
    type(matrix)         :: D_0
    real(8), intent(in)  :: ao(*)
    real(8)              :: buffer(*)
    integer:: ipoint, icomp
    integer:: i

!   here starts the expensive loop over points    
    do ipoint = 1, nr_points

      !call get_quantity_in_point(nr_dmat,                   &
      !                           D,                         &
      !                           D_0,                       &
      !                           grid%points(1, ipoint),    &
      !                           grid%points(2, ipoint),    &
      !                           grid%points(3, ipoint),    &
      !                           ao,                        &
      !                           buffer,                    &
      !                           grid_function%get_ndim(),  &
      !                           grid_function%values(1, ipoint),      &
      !                           grid_function, ipoint)


      call get_quantity_in_point_refactored(nr_dmat,                   &
                                 D,                         &
                                 D_0,                       &
                                 ao,                        &
                                 buffer,                    &
                                 grid, &
                                 grid_function, &
                                 grid_function_values, &
                                 ipoint)

      !write(*,*) 'grid_function%values on node', get_my_rank()
      !do i = 1, 147
      ! if (ao(i) .ne. 0.0d0) write(*,'(A,3I3,F8.3)') 'ao work: ', get_my_rank(), ipoint, i, ao(i)
      !end do
      !write(*,*) 'grid_function%values in ', ipoint, &
      !            grid%points(1, ipoint),& 
      !            grid%points(2, ipoint),&
      !            grid%points(3, ipoint),&
      !            grid_function_values(1:9,ipoint)
    end do
    !do icomp=1, grid_function%get_ndim()
    !  write(*,*) 'grid_function%values', (grid_function%values(icomp,ipoint), ipoint=1,grid%nr_points_total)
    !end do

  end subroutine work_on_grid

!   ---- legacy ----    

   subroutine visual_driver(nr_dmat, &
                            D,       &
                            D_0,     &
                            ao,      &
                            buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_dmat
      type(matrix)         :: D(nr_dmat)
      type(matrix)         :: D_0
      real(8), intent(in)  :: ao(*)
      real(8)              :: buffer(*)
!     --------------------------------------------------------------------------
      real(8), allocatable:: q(:, :)
      integer              :: q_max_dim
      real(8), allocatable:: atom_xyz(:, :)
      real(8)              :: step_3d(3, 0:3)
      integer              :: i, j, k, l, m, istep, n
      integer              :: icube, jcube, kcube
      integer              :: ix, iy
      real(8)              :: vx(3), vx_norm
      real(8)              :: vy(3), vy_norm
      real(8)              :: vz(3), vz_norm
      real(8)              :: r(3), v(3), t(9), w(0:12), s, rad, radtot, qbuf(0:12)
      integer              :: lunit_2d_scalar
      integer              :: lunit_2d_vector
      integer              :: lunit_2d_flow
      integer              :: lunit_2d_tensor
      integer              :: lunit_3d_cube
      integer              :: lunit_3d_scalar
      integer              :: lunit_3d_vector
      integer              :: lunit_3d_tensor
      integer              :: lunit_3d_vector_vtk
      integer              :: lunit_3d_jbfile
      integer              :: lunit_line_scalar
      integer              :: lunit_line_vector
      integer              :: lunit_line_tensor
      integer              :: lunit_3d_gridimp
      integer              :: lunit_radial_scalar
      integer              :: lunit_radial_vector
      integer              :: lunit_radial_tensor
      integer              :: lunit_outfil_custom_scalar
      integer              :: lunit_outfil_custom_vector
      integer              :: lunit_outfil_custom_2ndtensor
      real(8), allocatable:: temp(:, :)
      real(8), allocatable:: temp2(:, :, :, :)
      real(8), allocatable:: res(:)
      real(8), allocatable:: res_temp(:)
      real(8), allocatable:: xyz(:, :)
      real(8), allocatable:: xyzk(:, :)
      real(8)              :: acos, dsum
      real(8)              :: norm_grad
      integer              :: nleb
      integer              :: nr_points(3), nr_points_total
      integer              :: nr_points_gridimp
      logical              :: is_open

      integer              :: ifile
      logical              :: refactoring
      type(visual_plot_file), allocatable:: visual_inp_data_files(:)
      type(visual_plot_file), allocatable:: visual_plot_files(:)
      type(visual_grid_file), allocatable:: visual_grid_files(:)
      real(8), allocatable:: grid_point(:)
      real(8), allocatable:: grid_function(:)
      real(8), allocatable:: rx(:)
      real(8), allocatable:: ry(:)
      real(8), allocatable:: rz(:)
!     --------------------------------------------------------------------------

      i = 1
      if (visual_cfg_ncube(3) > 0) then
         i = visual_cfg_ncube(3)
      end if
      if (visual_cfg_tensor) then
        q_max_dim = 12
        allocate(q(0:12, i))
        q(0:12, 1:i) = 0.0d0
      else
        q_max_dim = 3
        allocate(q(0:3, i))
        q(0:3, 1:i) = 0.0d0
      end if

      allocate(atom_xyz(4, nucdep))


!     prepare input data (incoming grids and grid functions)
!     ======================================================

      if (visual_cfg_3d_gridimp .or. visual_cfg_3d_readjb) then
!        gosia: if we import grid on which we plot, then don't do*cube or*vtk files
!        as this grid is not necessarily uniform.
!        If it is uniform, then use the script lines_to_vtk.py or lines_to_cube.py in DIRAC's utils directory 
         visual_cfg_london_skip_cube = .true.
         visual_cfg_london_skip_vtk  = .true.

         call open_plot_files_on_import(visual_cfg_inp_data_files_use_legacy, &
                                        visual_inp_data_files)
      end if
         


!     get grid points
!     ===============

      !refactoring = .false.
      refactoring = .true.
      if (visual_cfg_3d .and. refactoring) then

         ! generate 3d cube grid

         allocate(visual_grid_files(visual_cfg_nr_grid_files))

         call open_grid_files(visual_grid_files)

         call print_grid_files_info(visual_grid_files)

         call generate_grid_file(visual_grid_files, &
               (/visual_cfg_ncube(1), visual_cfg_ncube(2), visual_cfg_ncube(3)/))
         call close_grid_files(visual_grid_files)

         ! read grid
         allocate(rx(visual_cfg_ncube(1)))
         allocate(ry(visual_cfg_ncube(2)))
         allocate(rz(visual_cfg_ncube(3)))
         rx = 0.0d0
         ry = 0.0d0
         rz = 0.0d0

         call open_old_grid_txt_file(visual_grid_files(1))
         call read_grid_files(rx, ry, rz, visual_grid_files, &
               (/visual_cfg_ncube(1), visual_cfg_ncube(2), visual_cfg_ncube(3)/))

         call close_grid_files(visual_grid_files)

         deallocate(visual_grid_files)

         do i = 1, visual_cfg_ncube(1)*visual_cfg_ncube(2)*visual_cfg_ncube(3)
            write(*,*) 'POINTS GRID: ', i, rx(i), ry(i), rz(i)
         end do
      end if



!     open files for outputting data on meshes
!     ========================================

      if (visual_cfg_2d .or. &
          visual_cfg_3d .or. visual_cfg_3d_fast .or. &
          visual_cfg_3d_gridimp .or. visual_cfg_3d_readjb .or. &
          visual_cfg_list .or. visual_cfg_line .or. visual_cfg_radial) then

         if (visual_cfg_print_data_to_dirac_output) then
            ! write to DIRAC output
            call write_data_to_dirac_output()
         else
            ! write to files
            write(*,*) 'BUBA ', visual_cfg_txt_plot_file_use_legacy, &
                                visual_cfg_cube_plot_file_use_legacy, &
                               (visual_cfg_txt_plot_file_use_legacy &
                           .or. visual_cfg_cube_plot_file_use_legacy)
            call open_plot_files_for_export((visual_cfg_txt_plot_file_use_legacy &
                                        .or. visual_cfg_cube_plot_file_use_legacy), &
                                            visual_plot_files)
            if (allocated(visual_plot_files)) then
               call print_plot_files_info(visual_plot_files)
            end if
         end if

      end if


      if (visual_cfg_list) then

         do i = 1, visual_cfg_nr_points_in_list
           call get_quantity_in_point(visual_cfg_nr_dmat,        &
                                      D,                         &
                                      D_0,                       &
                                      visual_cfg_xyz_list(i, 1), &
                                      visual_cfg_xyz_list(i, 2), &
                                      visual_cfg_xyz_list(i, 3), &
                                      ao,                        &
                                      buffer,                    &
                                      q_max_dim,                 &
                                      q)

           call write_plot_data(visual_plot_files, &
                                (/visual_cfg_xyz_list(i, 1), visual_cfg_xyz_list(i, 2), visual_cfg_xyz_list(i, 3)/), &
                                q(0:q_max_dim, 1))
         end do
           
      end if


      if (visual_cfg_line) then

         v = visual_cfg_line_to-visual_cfg_line_from
         v = v/visual_cfg_line_nr_steps
     
         do istep = 0, visual_cfg_line_nr_steps
            r = visual_cfg_line_from+istep*v
            call get_quantity_in_point(visual_cfg_nr_dmat, &
                                       D,                  &
                                       D_0,                &
                                       r(1),               &
                                       r(2),               &
                                       r(3),               &
                                       ao,                 &
                                       buffer,             &
                                       q_max_dim,          &
                                       q)
           
            call write_plot_data(visual_plot_files, &
                                (/r(1), r(2), r(3)/), &
                                q(0:q_max_dim, 1))
         end do
      end if

      if (visual_cfg_radial) then

!   Get angular points and weights
!   ==============================

         m = NLEB(num_grid_cfg_angint)
         write(LUPRI, '(/,4x, a)')  '** Radial distribution **'
         write(LUPRI, '(6x, A, I3)')  '   Lebedev angular grid L = ',num_grid_cfg_angint
         write(LUPRI, '(6x, A, I10)') '   Angular points:   ',m       
         allocate(temp(m, 4))
         CALL LEBEDEV(num_grid_cfg_angint, temp(1, 1), temp(1, 2), temp(1, 3), temp(1, 4), n)
         CALL DSCAL(m, 4*ACOS(-1.0D0), temp(1, 4), 1)

!   Loop over radial points
!   =======================

         radtot = 0.0D0
         s = visual_cfg_radial_length/visual_cfg_radial_nr_steps
         do istep = 0, visual_cfg_radial_nr_steps

            rad = istep*s

!           Angular integration
!           ===================
!gosia TODO: adapt to tensors
            q(0:3, 1) = 0.0d0
            do i = 1, m
              r   = visual_cfg_radial_from+rad*temp(i, 1:3)
              call get_quantity_in_point(visual_cfg_nr_dmat, &
                                         D,                  &
                                         D_0,                &
                                         r(1),               &
                                         r(2),               &
                                         r(3),               &
                                         ao,                 &
                                         buffer,             &
                                         3,                  &
                                         qbuf)
              q(0:3, 1) = q(0:3, 1) + qbuf(0:3)*temp(i, 4)
            end do
            q(0:3, 1) = q(0:3, 1)*rad*rad
            radtot = radtot+q(0, 1)
            
!           Write results
!           =============

            call write_plot_data(visual_plot_files, &
                                (/rad/), &
                                q(0:12, 1))


         end do
         write(lupri, *) '        The radial distribution integrates to: ',radtot*s
         deallocate(temp)
      end if

      if (visual_cfg_2d) then

         vx = (visual_cfg_2d_p_right-visual_cfg_2d_p_origin)/visual_cfg_2d_nr_right
         vy = (visual_cfg_2d_p_top   - visual_cfg_2d_p_origin)/visual_cfg_2d_nr_top

         vz(1) = vx(2)*vy(3) - vx(3)*vy(2)
         vz(2) = vx(3)*vy(1) - vx(1)*vy(3)
         vz(3) = vx(1)*vy(2) - vx(2)*vy(1)

         vx_norm = dsqrt(vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))
         vy_norm = dsqrt(vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))
         vz_norm = dsqrt(vz(1)*vz(1) + vz(2)*vz(2) + vz(3)*vz(3))

         allocate(grid_point(2))
         i = 0
         do ix = 0, visual_cfg_2d_nr_right
            do iy = 0, visual_cfg_2d_nr_top
               i = i+1
               r = visual_cfg_2d_p_origin+ix*vx+iy*vy
               call get_quantity_in_point(visual_cfg_nr_dmat, &
                                          D,                  &
                                          D_0,                &
                                          r(1),               &
                                          r(2),               &
                                          r(3),               &
                                          ao,                 &
                                          buffer,             &
                                          3,                  &
                                          q)

               s       = q(0,   1)
               v(1:3)  = q(1:3, 1)
!              gosia TODO add tensors
               !t(4:12) = q(4:12, 1)

               grid_point(1) = ix*vx_norm-vx_norm*visual_cfg_2d_nr_right/2
               grid_point(2) = iy*vy_norm-vy_norm*visual_cfg_2d_nr_top/2
               ! FIXME!
               !if flow then
               !   allocate(grid_function(1))
               !   grid_function(1) = (v(1)*vz(1) + v(2)*vz(2) + v(3)*vz(3))/vz_norm
               !else if scalar then
               !   allocate(grid_function(1))
               !   grid_function(1) = s
               !else if vector then
               !   allocate(grid_function(2))
               !   grid_function(1) = (v(1)*vx(1) + v(2)*vx(2) + v(3)*vx(3))/vx_norm
               !   grid_function(2) = (v(1)*vy(1) + v(2)*vy(2) + v(3)*vy(3))/vy_norm
               !end if
               !call write_plot_data(visual_plot_files, &
               !                     grid_point, &
               !                     grid_function)
               !deallocate(grid_function)
            end do
         end do
         deallocate(grid_point)
      end if

      if (visual_cfg_3d_gridimp) then

!        read grid points from visual_cfg_3d_gridfil file
         rewind(lunit_3d_gridimp)

         do 
            read(lunit_3d_gridimp, *, end = 98) r(1), r(2), r(3)
            if (visual_cfg_3d_grid_units .eq. 'angstrom') then
                r = r/xtang
            end if
            call get_quantity_in_point(visual_cfg_nr_dmat, &
                                       D,                  &
                                       D_0,                &
                                       r(1),               &
                                       r(2),               &
                                       r(3),               &
                                       ao,                 &
                                       buffer,             &
                                       q_max_dim,          &
                                       w)

            call write_plot_data(visual_plot_files, &
                                r(1:3), &
                                w(0:12))

         end do
 98      continue

      end if

      if (visual_cfg_3d_readjb) then

!        read grid points and magnetically-induced current density from visual_cfg_3d_jbfile
         rewind(lunit_3d_jbfile)

         do 
            read(lunit_3d_jbfile, *, end = 99) r(1), r(2), r(3), v(1), v(2), v(3)
            call get_quantity_in_point(visual_cfg_nr_dmat, &
                                       D,                  &
                                       D_0,                &
                                       r(1),               &
                                       r(2),               &
                                       r(3),               &
                                       ao,                 &
                                       buffer,             &
                                       q_max_dim,          &
                                       w, &
                                       jb_input = v)

            call write_plot_data(visual_plot_files, &
                                r(1:3), &
                                w(0:12))

         end do
 99      continue

      end if


      if (visual_cfg_3d .or. visual_cfg_3d_fast) then

      refactoring=.true.
      !refactoring=.false.
      if (refactoring) then
         write(*,*) 'refactoring-3d!'

         !do icube = 1, visual_cfg_ncube(1)
         !   do jcube = 1, visual_cfg_ncube(2)
         !      do kcube = 1, visual_cfg_ncube(3)
         !         r(1) = rx(icube)
         !         r(2) = ry(jcube)
         !         r(3) = rz(kcube)

         !         call get_quantity_in_point(visual_cfg_nr_dmat, &
         !                                    D,                  &
         !                                    D_0,                &
         !                                    r(1),               &
         !                                    r(2),               &
         !                                    r(3),               &
         !                                    ao,                 &
         !                                    buffer,             &
         !                                    q_max_dim,          &
         !                                    q(0, kcube))

         !         call write_plot_data(visual_plot_files, &
         !                             r(1:3), &
         !                             q(0:q_max_dim, kcube))

         !      end do
         !   end do
         !end do

      else

!        step_3d(*, 0)  coordinates of initial point
!        step_3d(*, 1)  step vector for the "x"-direction
!        step_3d(*, 2)  step vector for the "y"-direction
!        step_3d(*, 3)  step vector for the "z"-direction

!        atom_xyz(1, i) x-coordinate of atom i
!        atom_xyz(2, i) y-coordinate of atom i
!        atom_xyz(3, i) z-coordinate of atom i
!        atom_xyz(4, i) charge of atom i

         step_3d = 0.0d0

         do i = 1, 3
            step_3d(i, 0) = cord(i, 1)
            step_3d(i, i) = cord(i, 1)
         enddo

         i = 0
         do j = 1, nucind
            k = istbnu(j)
            do l = 0, maxopr
               if (iand(l, k) == 0) then
                  i = i+1
                  do m = 1, 3
                     atom_xyz(m, i) = pt(iand(isymax(m, 1), l))*cord(m, j)

                     step_3d(m, 0) = min(step_3d(m, 0), atom_xyz(m, i))
                     step_3d(m, m) = max(step_3d(m, m), atom_xyz(m, i))
                  end do
                  atom_xyz(4, i) = charge(j)
               end if
            end do
         end do

!        add some space, calculate step sizes and convert step vectors to angstrom

         do i = 1, 3
            step_3d(i, 0) = (step_3d(i, 0) - (visual_cfg_add_3d*2.0d0))*xtang
            step_3d(i, i) = (step_3d(i, i) + (visual_cfg_add_3d*2.0d0))*xtang
            step_3d(i, i) = (step_3d(i, i) -  step_3d(i, 0))/(visual_cfg_ncube(i) - 1)
         end do

!        write header

         if (.not. visual_cfg_london_skip_cube) then
            write(lunit_3d_cube, '(a)') 'cube file'
            write(lunit_3d_cube, '(a)') 'generated by DIRAC'
            
            write(lunit_3d_cube, '(i5, 4f12.6)') &
               nucdep,                          &
              (step_3d(j, 0), j = 1, 3)
            do i = 1, 3
              write(lunit_3d_cube, '(i5, 4f12.6)') &
                 visual_cfg_ncube(i),             &
                (step_3d(j, i), j = 1, 3)
            end do
            
            do i = 1, nucdep
              write(lunit_3d_cube, '(i5, 4f12.6)') &
                   nint(atom_xyz(4, i)),           &
                   atom_xyz(4, i),                 &
                  (atom_xyz(j, i), j = 1, 3)
            end do
         end if

         if (.not. visual_cfg_london_skip_vtk) then

            write(lunit_3d_vector_vtk, '(a)')      '# vtk DataFile Version 2.0'
            write(lunit_3d_vector_vtk, '(a)')      'Sample rectilinear grid'
            write(lunit_3d_vector_vtk, '(a)')      'ASCII'
            write(lunit_3d_vector_vtk, '(a)')      'DATASET RECTILINEAR_GRID'
            write(lunit_3d_vector_vtk, '(a, 3i8)') 'DIMENSIONS',        &
                                                   visual_cfg_ncube(1), &
                                                   visual_cfg_ncube(2), &
                                                   visual_cfg_ncube(3)
         end if

         n = maxval((/visual_cfg_ncube(1), visual_cfg_ncube(2), visual_cfg_ncube(3)/))
         allocate(temp(n, 3))

         do icube = 1, visual_cfg_ncube(1)
            do jcube = 1, visual_cfg_ncube(2)
               do kcube = 1, visual_cfg_ncube(3)
                  r(1) = step_3d(1, 0) + (icube-1)*step_3d(1, 1) &
                                       + (jcube-1)*step_3d(1, 2) &
                                       + (kcube-1)*step_3d(1, 3)
                  r(2) = step_3d(2, 0) + (icube-1)*step_3d(2, 1) &
                                       + (jcube-1)*step_3d(2, 2) &
                                       + (kcube-1)*step_3d(2, 3)
                  r(3) = step_3d(3, 0) + (icube-1)*step_3d(3, 1) &
                                       + (jcube-1)*step_3d(3, 2) &
                                       + (kcube-1)*step_3d(3, 3)
                  temp(icube, 1) = r(1)
                  temp(jcube, 2) = r(2)
                  temp(kcube, 3) = r(3)
               end do
            end do
         end do

         if (.not. visual_cfg_london_skip_vtk) then
            j = visual_cfg_ncube(1)
            write(lunit_3d_vector_vtk, '(a, i8, a)') 'X_COORDINATES', j, ' float'
            write(lunit_3d_vector_vtk, '(6f13.5)')   (temp(m, 1), m = 1, j)
            j = visual_cfg_ncube(2)
            write(lunit_3d_vector_vtk, '(a, i8, a)') 'Y_COORDINATES', j, ' float'
            write(lunit_3d_vector_vtk, '(6f13.5)')   (temp(m, 2), m = 1, j)
            j = visual_cfg_ncube(3)
            write(lunit_3d_vector_vtk, '(a, i8, a)') 'Z_COORDINATES', j, ' float'
            write(lunit_3d_vector_vtk, '(6f13.5)')   (temp(m, 3), m = 1, j)
         end if

         deallocate(temp)

         if (.not. visual_cfg_london_skip_vtk) then
            n = visual_cfg_ncube(1)*visual_cfg_ncube(2)*visual_cfg_ncube(3)
            write(lunit_3d_vector_vtk, '(a, i8)')    'POINT_DATA', n
            write(lunit_3d_vector_vtk, '(a)')        'SCALARS scalars float'
            write(lunit_3d_vector_vtk, '(a)')        'LOOKUP_TABLE default'
            
            do i = 1, n
               write(lunit_3d_vector_vtk, '(f13.5)') 0.0d0
            end do
            
            write(lunit_3d_vector_vtk, '(a)')
            write(lunit_3d_vector_vtk, '(a)')
            write(lunit_3d_vector_vtk, '(a)')       'VECTORS vectors float'
         end if

         if (visual_cfg_tensor) then
           allocate(temp2(visual_cfg_ncube(1), visual_cfg_ncube(2), visual_cfg_ncube(3), 12))
         else
           allocate(temp2(visual_cfg_ncube(1), visual_cfg_ncube(2), visual_cfg_ncube(3), 3))
         end if

         if (visual_cfg_3d_fast) then
            if (visual_cfg_ncube(3) > max_batch_length) then
               print *, 'error: visual_cfg_ncube(3) > max_batch_length'
               stop 1
            end if
            allocate(res(visual_cfg_ncube(3)))
            allocate(res_temp(visual_cfg_ncube(3)))
            allocate(xyz(3, visual_cfg_ncube(3)))
            allocate(xyzk(visual_cfg_ncube(3), 3))
         end if

         do icube = 1, visual_cfg_ncube(1)
            do jcube = 1, visual_cfg_ncube(2)

               if (visual_cfg_3d_fast) then
                  do kcube = 1, visual_cfg_ncube(3)
                     xyz(1, kcube) = step_3d(1, 0) + (icube-1)*step_3d(1, 1) &
                                                   + (jcube-1)*step_3d(1, 2) &
                                                   + (kcube-1)*step_3d(1, 3)
                     xyz(2, kcube) = step_3d(2, 0) + (icube-1)*step_3d(2, 1) &
                                                   + (jcube-1)*step_3d(2, 2) &
                                                   + (kcube-1)*step_3d(2, 3)
                     xyz(3, kcube) = step_3d(3, 0) + (icube-1)*step_3d(3, 1) &
                                                   + (jcube-1)*step_3d(3, 2) &
                                                   + (kcube-1)*step_3d(3, 3)
                  end do
                  do kcube = 1, visual_cfg_ncube(3)
                     xyzk(kcube, 1) = xyz(1, kcube)
                     xyzk(kcube, 2) = xyz(2, kcube)
                     xyzk(kcube, 3) = xyz(3, kcube)
                  end do
                  call get_ao(visual_cfg_ncube(3), &
                              xyzk(1, 1),           &
                              xyzk(1, 2),           &
                              xyzk(1, 3),           &
                              ao,                  &
                              buffer)
                  call get_esp(visual_cfg_ncube(3),      &
                               res_temp,                 &
                               D(1)%irep,                &
                               D(1)%nrow,                &
                               D(1)%elms,                &
                               xyz,                      &
                               include_nuc_part=.false., &
                               include_el_part=.true.)
                  res = res_temp*0.5d0
                  call get_esp(visual_cfg_ncube(3),     &
                               res_temp,                &
                               D(1)%irep,               &
                               D(1)%nrow,               &
                               D(1)%elms,               &
                               xyz,                     &
                               include_nuc_part=.true., &
                               include_el_part=.false.)
                  res = res+res_temp
                  do kcube = 1, visual_cfg_ncube(3)
                     q(0, kcube) = res(kcube)
                  end do
               else
                  do kcube = 1, visual_cfg_ncube(3)

                     r(1) = step_3d(1, 0) + (icube-1)*step_3d(1, 1) &
                                          + (jcube-1)*step_3d(1, 2) &
                                          + (kcube-1)*step_3d(1, 3)
                     r(2) = step_3d(2, 0) + (icube-1)*step_3d(2, 1) &
                                          + (jcube-1)*step_3d(2, 2) &
                                          + (kcube-1)*step_3d(2, 3)
                     r(3) = step_3d(3, 0) + (icube-1)*step_3d(3, 1) &
                                          + (jcube-1)*step_3d(3, 2) &
                                          + (kcube-1)*step_3d(3, 3)

                     call get_quantity_in_point(visual_cfg_nr_dmat, &
                                                D,                  &
                                                D_0,                &
                                                r(1),               &
                                                r(2),               &
                                                r(3),               &
                                                ao,                 &
                                                buffer,             &
                                                q_max_dim,          &
                                                q(0, kcube))

                     call write_plot_data(visual_plot_files, &
                                         r(1:3), &
                                         q(0:q_max_dim, kcube))


                     if (visual_cfg_tensor) then
                       temp2(icube, jcube, kcube, 1:12) = q(1:12, kcube)
                     else
                       temp2(icube, jcube, kcube, 1:3) = q(1:3, kcube)
                     end if

                     if ((visual_cfg_add_grid_points_grad_ed .or. &
                          visual_cfg_add_grid_points_grad_ed_scaled).and.      &
                        (visual_cfg_property(1) == iq_density_gradient)) then
                        if (visual_cfg_add_grid_points_grad_ed_scaled) then
                           norm_grad = q(1, kcube)*q(1, kcube) &
                                     + q(2, kcube)*q(2, kcube) &
                                     + q(3, kcube)*q(3, kcube)
                           norm_grad = sqrt(norm_grad)
                           r(1) = r(1) + q(1, kcube)/norm_grad
                           r(2) = r(2) + q(2, kcube)/norm_grad
                           r(3) = r(3) + q(3, kcube)/norm_grad
                        else
                           r(1) = r(1) + q(1, kcube)
                           r(2) = r(2) + q(2, kcube)
                           r(3) = r(3) + q(3, kcube)
                        end if
                        call get_quantity_in_point(visual_cfg_nr_dmat, &
                                                   D,                  &
                                                   D_0,                &
                                                   r(1),               &
                                                   r(2),               &
                                                   r(3),               &
                                                   ao,                 &
                                                   buffer,             &
                                                   q_max_dim,          &
                                                   q(0, kcube))

                        call write_plot_data(visual_plot_files, &
                                            r(1:3), &
                                            q(0:q_max_dim, kcube))

                     end if

                  end do
               end if
               if (.not. visual_cfg_london_skip_cube) then
                  write(lunit_3d_cube, '(6ES18.8E3)')  &
                     (q(0, m), m = 1, visual_cfg_ncube(3))
               end if
            end do
         end do

         if (visual_cfg_3d_fast) then
            deallocate(res)
            deallocate(res_temp)
            deallocate(xyz)
            deallocate(xyzk)
         end if

         if (.not. visual_cfg_london_skip_vtk) then
            do kcube = 1, visual_cfg_ncube(3)
               do jcube = 1, visual_cfg_ncube(2)
                  do icube = 1, visual_cfg_ncube(1)
                     write(lunit_3d_vector_vtk, '(3ES18.8E3)') &
                        temp2(icube, jcube, kcube, 1:3)
                  end do
               end do
            end do
         end if
         deallocate(temp2)
      end if  ! if (refactoring) then
      end if

      if (allocated(visual_plot_files)) then
         call close_plot_txt_files(visual_plot_files)
         deallocate(visual_plot_files)
      end if

      is_open=.false.
      inquire(lunit_3d_gridimp, opened = is_open)
      if (is_open) close(lunit_3d_gridimp, status = 'keep')

      is_open=.false.
      inquire(lunit_3d_jbfile, opened = is_open)
      if (is_open) close(lunit_3d_jbfile, status = 'keep')

      deallocate(q)
      deallocate(atom_xyz)
      if (allocated(rx)) deallocate(rx)
      if (allocated(ry)) deallocate(ry)
      if (allocated(rz)) deallocate(rz)

   end subroutine

   subroutine setup_orbital_string(occupation)

      ! this routine will set up occupation vector
      ! based on**VISUAL
      !          .OCCUPATION
      ! if this is not given, then it will occupy all
      ! "occupied" electrons with 1.0
      real(8):: occupation(norbt)

      integer:: i, ifsym, ioff

      occupation = 0.0d0

      ioff = 0
      do ifsym = 1, nfsym
         ioff = ioff+npsh(ifsym)
         if (visual_cfg_use_orbital_string) then
            do i = 1, nesh(ifsym)
               occupation(ioff+i) = visual_cfg_occupation(i, ifsym)
            end do
         else
            if (nash(ifsym) > 0) then
               print *, 'ERROR: with open shells you should use .OCCUPATION'
               print *, '       otherwise occupation will be wrong'
               print *, '       so we better stop here'
               stop 1
            end if
            do i = 1, nish(ifsym)
               occupation(ioff+i) = 1.0d0
            end do
         end if
         ioff = norb(ifsym)
      end do

   end subroutine

      SUBROUTINE SETUP_FO_DMAT_AO(IXVECTOR, &
     &                            STRING, D_mo)

      CHARACTER LAB1*16, &
     &          LAB2*16, &
     &          TYP1*2, &
     &          TYP2*2, &
     &          LBUF*16, &
     &          PBUF*1
      LOGICAL   FILE_EXISTS
      type(matrix):: D_mo
      real(8):: STRING(*)

      real(8), allocatable:: bvec(:)

      integer, allocatable:: orbital_rotation_indices(:, :)
      real(8), allocatable:: response_vector(:, :)
      integer              :: length
      character(2)         :: vector_type
      real(8)              :: ih_factor

      integer              :: nr_vectors, LUNIT_XVECTOR
      integer              :: ixvector, NBAS, NORBS, len
      integer              :: JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, IVECTOR
      real(8)              :: FREQ1, FREQ2, ERGRSP, unused

!gosia, i'll remove it later
      integer:: i, j
      character(10):: pattern
      character(16):: text

      INQUIRE(FILE = visual_cfg_dmat_file(IXVECTOR), EXIST = FILE_EXISTS)
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI, '(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI, '(A, A6, A)') 'solution vector file ', &
     &    visual_cfg_dmat_file(IXVECTOR), ' not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR, &
     &     FILE = visual_cfg_dmat_file(IXVECTOR), &
     &     STATUS = 'OLD', &
     &     ACCESS = 'SEQUENTIAL', &
     &     FORM = 'UNFORMATTED')

!     get number of vectors on file

      NR_VECTORS = 0
      REWIND LUNIT_XVECTOR
 1    CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 2
      ELSE
        NR_VECTORS = NR_VECTORS+1
!       skip solution vector
        READ(LUNIT_XVECTOR)
!       skip orbital rotations
        READ(LUNIT_XVECTOR)
        GOTO 1
      ENDIF
 2    CONTINUE

!     errors

      if ((visual_cfg_dmat_file_record(ixvector) > nr_vectors) &
          .or. (visual_cfg_dmat_file_record(ixvector) < 1)) then
         write(lupri, '(a)') 'error in visual:'
         write(lupri, '(a)') 'you try to read an invalid record number from file '//visual_cfg_dmat_file(ixvector)
         call quit('error in visual')
      end if

!     move to the right record

      IVECTOR = 0
      REWIND LUNIT_XVECTOR
 3    CONTINUE
      READ(LUNIT_XVECTOR)
      IVECTOR = IVECTOR+1
      IF(IVECTOR .EQ. visual_cfg_dmat_file_record(IXVECTOR)) THEN
        BACKSPACE(LUNIT_XVECTOR)
        GOTO 4
      ELSE
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        GOTO 3
      ENDIF
 4    CONTINUE


!     get type and dimensions of record IWHICH_XVECTOR(IXVECTOR)

      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2, &
     &                    JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, &
     &                    UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS



      visual_cfg_irep_xvector(IXVECTOR) = JSYMOP1-1
      length  = LEN/NZ
      D_mo%irep   = JSYMOP1-1


    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
      case ('P')
        vector_type = 'pn'
      case default
        call quit('unexpected typ1(1:1)')
    end select

    select case (typ1(2:2))
      case ('+')
        D_mo%ih_sym =  1
        ih_factor   =  1.0d0
      case ('-')
        D_mo%ih_sym = -1
        ih_factor   = -1.0d0
      case default
        call quit('unexpected typ1(2:2)')
    end select


!   allocate and read response vector
    call alloc(response_vector, length, nz)
    call readt(lunit_xvector, length*nz, response_vector)
!   allocate and read orbital rotation array
    call alloc(orbital_rotation_indices, 2, length)
    call readi(lunit_xvector, length*2, orbital_rotation_indices)
!   close file of solution vectors
    close(lunit_xvector, status = 'keep')

!fixme: this is a good place to zero elements according to occupation

    D_mo%elms = 0.0d0
    call scatter_vector(length,                   &
                        orbital_rotation_indices, &
                        ih_factor,                &
                        response_vector,          &
                        D_mo%elms,                &
                        D_mo%irep)

    call dealloc(orbital_rotation_indices)
    call dealloc(response_vector)

  end subroutine

  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out):: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              |
!              |  column
!              |  |
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine


  subroutine write_data_to_dirac_output()

!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------

    if (visual_cfg_2d .or. &
        visual_cfg_3d .or. visual_cfg_3d_fast .or. &
        visual_cfg_list .or. visual_cfg_line) then

       call around('Requested density at point(s)')
       
       write(lupri, *)
       if (visual_cfg_tensor) then
         write(lupri, '(a200)') &
              'in columns: r_x, r_y, r_z, '//&
              'scalar, x-component, y-component, z-component, '//&
              'xx-component, xy-component, xz-component, '//&
              'yx-component, yy-component, yz-component, '//&
              'zx-component, zy-component, zz-component'
       else
         write(lupri, '(a80)') &
              'in columns: r_x, r_y, r_z, '//&
              'scalar, x-component, y-component, z-component'
       end if
       
    else if (visual_cfg_radial) then
    end if

    write(lupri, *)

    if (visual_cfg_list) then
    end if
 
    !write(lupri, *)
    !call prsymb(lupri, '-', 100, 0)
    !write(lupri, *)
    !write(lupri, *)

  end subroutine


!   ---- end ----    
end module
