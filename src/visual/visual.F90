module visual

! This module contains routines that calculate 
! various densities on real-space grids,
! Written by Radovan Bast.
! Refactored by Gosia Olejniczak.
! Full description of features: LINK
!
! current limitations:
! * calculations are done on one molecular system,
! * data is read from DIRAC's CHECKPOINT.h5 (molecular information+wave function) and from RESPONSE.h5 (response parameters)
! * but AOs are generated in the "old" way (using `xcint/` modules) 

  use visual_cfg
  use visual_dirac
  use visual_in_point
  use visual_dirac
  use visual_grids
  use visual_utils
  use visual_grid_functions
  use visual_datatypes
  use visual_io
  use visual_txt_storage
  use visual_mpi
  use dirac_dmat_factory

! DIRAC modules
  use checkpoint
  use basis_set_datatypes
  use dirac_ao_eval
  use interface_ao

#ifdef VAR_MPI
  use interface_to_mpi
  use integer_kind_mpilib
#endif

  implicit none

  public run_visual


  private

contains

! -----------------------------------------------------------------------------
  subroutine run_visual(launch_only)
    logical, intent(in), optional :: launch_only
    type(mol_t):: mol
    type(cmo_t):: mo_data
    type(basis_set_info_t):: ao_basis
    real(8) :: walltime
    type(dmat_t):: dmat_0
! -----------------------------------------------------------------------------

    walltime = visual_walltime()

!
!   data available on input:
!     * information from user input stored in rs_input dt
!     * config information in visual_cfg* variables (legacy)
!   
    call visual_initialize_and_sync_nodes(launch_only)
#ifdef VAR_MPI
    if ( rs_input%print_level > 2 ) then
      print *, 'VISUAL parallel? ', is_parallel()
      print *, 'size of communicator: ', mpi_get_nr_proc()
      print *, 'rank of the process: ', get_my_rank()
    end if
#endif

    if (mpi_is_master()) then
      call visual_input_show_setup()
    end if

    if (mpi_is_master()) then
!
!     fixme:
!     ao_basis is now a placeholder:
!     now, AOs are still calculated in an external dirac_ao_eval_init,
!     as it also sets up all offsets used in xcint/density_eval.F90;
!     to fix this, one would need to think of decoupling visual from xcint,
!     which would require rewriting visual_in_point.F90
!
      call visual_initialize_data_from_dirac(mol, mo_data, ao_basis)
      call visual_update_data(mo_data)
      call dirac_dmat_init('unperturbed', mo_data, dmat_ao=dmat_0)
    end if

#ifdef VAR_MPI
    if (is_parallel()) then
      call synchronize_data_from_dirac(mol,mo_data)
      call synchronize_dmat(dmat_0)
    end if
#endif

!   todo
!   open rs-checkpoint (once a separate RS schema is ready)
!   write (lupri,'(/A,I2)'),&
!     "* Opening checkpoint file for real-space data (VISUAL)"
!   call checkpoint_rs_open
! 
    call visual_driver(mol,mo_data,dmat_0)

!   todo
!   call checkpoint_rs_close (once a separate RS schema is ready)
!
    call dmat_0%erase()

    if (allocated(rs_input%grids)) deallocate(rs_input%grids)
    if (allocated(rs_input%grid_functions)) deallocate(rs_input%grid_functions)
    call rs_input%erase()
    if (mol%in_use()) then
      call mol%erase()
    end if

    walltime = walltime - visual_walltime()
    if (mpi_is_master()) print *, 'Elapsed walltime in run_visual: ', walltime

  end subroutine


! ------------------------------------------------------------------------------
  subroutine visual_driver(mol,cmo,dmat_0)
    type(mol_t), intent(in):: mol
    type(cmo_t), intent(in):: cmo
    type(dmat_t), intent(in):: dmat_0

    type(dmat_t):: dmat, dmat_lao
    real(8), allocatable:: ao(:,:)
    real(8), allocatable:: buffer(:,:)
    type(grid_t):: grid
    type(grid_t), allocatable:: grids(:)
    type(grid_function_t), allocatable:: grid_functions(:)
    real(8), allocatable:: grid_function_values(:,:)
    real(8), allocatable:: grid_function_integrated_values(:)
    real(8), allocatable:: radial_values(:,:)
    integer:: ig, igf
    integer:: g_ndim, gf_ndim, npoints, npoints_radial
    integer:: my_rank
    real(8) :: walltime
    real(8):: rad, radial_spacing
    integer:: n_this_proc, idim_grid, idim_gf
    integer:: i,j,ipoint,istep,ii
    integer:: iz, off_ij
! ------------------------------------------------------------------------------

    do ig = 1, rs_input%nr_grids
!
!     get necessary information about the grid and share with workers
!

      if (mpi_is_master()) then
        call get_grid_info(ig, rs_input%grids(ig), mol, dmat_0)
      end if
#ifdef VAR_MPI
      if (is_parallel()) call synchronize_grid_data(rs_input%grids(ig))
#endif

      g_ndim  = rs_input%grids(ig)%ndim
      npoints = rs_input%grids(ig)%nr_points_total

      if (.not. rs_input%grids(ig)%allocated_points()) then
        allocate(rs_input%grids(ig)%points(g_ndim, npoints))
        rs_input%grids(ig)%points = 0.0d0
      end if

      if (rs_input%grids(ig)%requires_weights() .and. &
        .not. rs_input%grids(ig)%allocated_weights()) then
        allocate(rs_input%grids(ig)%weights(npoints))
        rs_input%grids(ig)%weights = 0.0d0
      end if
      if (rs_input%grids(ig)%typ == 2) then
        npoints_radial=rs_input%grids(ig)%radial_nr_steps
        if (.not. rs_input%grids(ig)%allocated_radial_coord()) then
          allocate(rs_input%grids(ig)%radial_coord(npoints_radial))
          rs_input%grids(ig)%radial_coord = 0.0d0
        end if
      end if

!
!     get grid points
!
!     todo:
!     * test whether reshaping to 1-dimensional arrays is better for optimization:
!     grid_points = reshape(rs_input%grids(ig)%points(g_ndim, npoints))
!
!     fixme:
!     * the grid can be read from file or generated;
!       now, grid generation is done sequentially (on master), it could be done in parallel
!
      if (mpi_is_master()) then
        call get_grid_points(rs_input%grids(ig))
        if (rs_input%print_level > 1) then
          write(*,*) 'GRID POINTS (AND WEIGHTS): '
          do j = 1, npoints
            if (rs_input%grids(ig)%allocated_weights()) then
              write(*, *) (rs_input%grids(ig)%points(i,j),i=1,g_ndim), rs_input%grids(ig)%weights(j)
            else
              write(*, *) (rs_input%grids(ig)%points(i,j),i=1,g_ndim)
            end if
          end do
          if (rs_input%grids(ig)%allocated_radial_coord()) then
            write(*,*) 'GRID RADIAL POINTS: '
            do j = 1, npoints_radial
              write(*, *) rs_input%grids(ig)%radial_coord(j)
            end do
          end if
        end if
      end if

      if (rs_input%grids(ig)%export) then
!       export this grid to a separate file 
        if (mpi_is_master()) then
          call store_grid(rs_input%grids(ig))
          if (rs_input%print_level > 2) then
            write(*, *) 'Grid ', rs_input%grids(ig)%id, ' exported to file ', &
                         rs_input%grids(ig)%file_out
          end if
        end if
      end if

      do igf = 1, rs_input%nr_grid_functions

!
!       get necessary information about the grid function and share with workers

        if (mpi_is_master()) then 
          call get_grid_function_info(igf, rs_input%grid_functions(igf))
        end if
#ifdef VAR_MPI
        if (is_parallel()) call synchronize_grid_function_data(rs_input%grid_functions(igf))
#endif

        if (rs_input%grid_functions(igf)%id_assigned_grid == rs_input%grids(ig)%id) then

          ! it would be cleaner to use "associate", 
          ! but some compilers have troubles with nested associate
          ! https://stackoverflow.com/questions/66562058/associate-chains-in-fortran

          if (rs_input%print_level > 0) then
            if (mpi_is_master()) then 
              write(*, *) 'Working on grid function ', igf, ', grid ', ig
              call rs_input%grids(ig)%display
              call rs_input%grid_functions(igf)%display
            end if
          end if

!
!         initialize ao module
!
          call dirac_ao_eval_init(rs_input%grid_functions(igf)%max_ao_g_order, &
                                  rs_input%grid_functions(igf)%max_ao_m_order, &
                                  .true.)

          allocate(ao(1,nr_ao_slices*nr_ao_cartesian))
          ao = 0.0d0
          allocate(buffer(1,nr_ao_slices*nr_ao_cartesian))
          buffer = 0.0d0

!
!         prepare density matrices
!
          if (rs_input%grid_functions(igf)%wf_order > 0) then
            if (mpi_is_master()) then
              call dirac_dmat_init('xlr', cmo, &
                                   grid_function=rs_input%grid_functions(igf), &
                                   dmat_ao=dmat)

              if (rs_input%grid_functions(igf)%lao > 0) then

                if (rs_input%grid_functions(igf)%lao_skip(1) == 1) then
!                 skip response contribution ('nokappa')
                  dmat%elements = 0.0d0
                end if

                if (rs_input%grid_functions(igf)%lao_skip(3) == 0) then
!                 do not skip reorthonormalization contribution ('noreorth')
                  call dirac_dmat_init('lao-connection-1', cmo, &
                                       grid_function=rs_input%grid_functions(igf), &
                                       dmat_ao=dmat_lao)
                  dmat%elements = dmat%elements + dmat_lao%elements
                  call dmat_lao%erase()
                else
                  if (rs_input%grid_functions(igf)%lao_skip(2) == 0) then
!                   do not skip direct contribution ('nodirect')
                    call get_xyzcomp_from_rspvec(rs_input%grid_functions(igf))
                  end if
                end if

              end if

            end if
#ifdef VAR_MPI
            if (is_parallel()) then
              call synchronize_dmat(dmat)
              if (rs_input%grid_functions(igf)%lao_skip(2)==0 &
                .or. rs_input%grid_functions(igf)%lao_skip(3)==0) then
                call interface_mpi_bcast(rs_input%grid_functions(igf)%lao_comp,&
                     1,irank_master,global_communicator)
              end if
            end if
#endif
          end if

!
!         distribute grid points over nodes
!

          n_this_proc = npoints

#ifdef VAR_MPI
          if (is_parallel()) then
            ! redefine n_this_proc
            do idim_grid = 1, g_ndim
              call visual_interface_mpi_distribute_points(&
                         rs_input%grids(ig)%points(idim_grid,:), &
                         npoints, n_this_proc)
            end do
            if (rs_input%grids(ig)%allocated_weights()) then
              call visual_interface_mpi_distribute_points(&
                         rs_input%grids(ig)%weights(:), &
                         npoints, n_this_proc)
            end if
!           no need to distribute rs_input%grids(ig)%radial_coord
!           we take care of it after sync
          end if
#endif

!
!         work on grid
!         note: here, we always work on 3d grid points
!

          gf_ndim  = rs_input%grid_functions(igf)%get_ndim()

          allocate(grid_function_values(gf_ndim, n_this_proc))
          grid_function_values = 0.0d0
          if (rs_input%grid_functions(igf)%purpose == 1) then
            allocate(grid_function_integrated_values(gf_ndim))
            grid_function_integrated_values = 0.0d0
          end if
 
          if (rs_input%grid_functions(igf)%purpose == 1.or. &
              rs_input%grids(ig)%requires_weights()) then
            call work_on_grid(g_ndim, n_this_proc,&
                              mol,&
                              rs_input%grids(ig),&
                              rs_input%grids(ig)%points(1:g_ndim,1:n_this_proc),&
                              rs_input%grid_functions(igf),&
                              grid_function_values,&
                              dmat_0,&
                              dmat,&
                              ao, buffer,&
                              rs_input%grids(ig)%weights(1:n_this_proc),&
                              grid_function_integrated_values)
          else
            call work_on_grid(g_ndim, n_this_proc,&
                              mol,&
                              rs_input%grids(ig),&
                              rs_input%grids(ig)%points(1:g_ndim,1:n_this_proc),&
                              rs_input%grid_functions(igf),&
                              grid_function_values,&
                              dmat_0,&
                              dmat,&
                              ao, buffer)
          end if

          deallocate(ao)
          deallocate(buffer)

          allocate(rs_input%grid_functions(igf)%values(gf_ndim, npoints))
          rs_input%grid_functions(igf)%values = 0.0d0
          if (rs_input%grid_functions(igf)%purpose == 1) then
            allocate(rs_input%grid_functions(igf)%integrated_values(gf_ndim))
            rs_input%grid_functions(igf)%integrated_values = 0.0d0
          end if
 
          if (rs_input%print_level > 1) then
            write(*, *) 'Grid_function values on this proc: '
            do j = 1, n_this_proc
              write(*, *) 'values: ', j, (grid_function_values(i,j),i=1,gf_ndim)
            end do
          end if
 
#ifdef VAR_MPI
          if (is_parallel()) then
            do idim_gf = 1, gf_ndim
              call visual_interface_mpi_gatherv(grid_function_values(idim_gf,:),&
                                                rs_input%grid_functions(igf)%values(idim_gf,:),&
                                                n_this_proc)
              if (rs_input%grid_functions(igf)%purpose == 1) then
                call visual_interface_mpi_reduce(&
                             grid_function_integrated_values(idim_gf), &
                             rs_input%grid_functions(igf)%integrated_values(idim_gf))
              end if
            end do
          else
            do idim_gf = 1, gf_ndim
              rs_input%grid_functions(igf)%values(idim_gf,:) = grid_function_values(idim_gf,:)
              if (rs_input%grid_functions(igf)%purpose == 1) then
                rs_input%grid_functions(igf)%integrated_values(idim_gf) = &
                         grid_function_integrated_values(idim_gf)
              end if
            end do
          end if
#else
          do idim_gf = 1, gf_ndim
            rs_input%grid_functions(igf)%values(idim_gf,:) = grid_function_values(idim_gf,:)
              if (rs_input%grid_functions(igf)%purpose == 1) then
                rs_input%grid_functions(igf)%integrated_values(idim_gf) = &
                         grid_function_integrated_values(idim_gf)
              end if
          end do
#endif

          if (mpi_is_master()) then

!           post-process non-3D grids

            if (rs_input%grids(ig)%typ == 2) then
!             radial grids
              if (rs_input%grid_functions(igf)%purpose == 0) then
!               visualization: plot data on the radial axis
              
                allocate(radial_values(gf_ndim, npoints_radial))
                radial_values = 0.0d0

                ipoint=1
                do istep = 0, npoints_radial-1
                  do ii = 1, rs_input%grids(ig)%nr_points_total/npoints_radial
                    radial_values(1:gf_ndim, istep+1) = &
                      radial_values(1:gf_ndim, istep+1) + &
                      rs_input%grid_functions(igf)%values(1:gf_ndim,ipoint)
                    ipoint = ipoint+1
                  end do
                end do

                deallocate(rs_input%grid_functions(igf)%values)
                allocate(rs_input%grid_functions(igf)%values(gf_ndim, npoints_radial))
                rs_input%grid_functions(igf)%values = 0.0d0
                rs_input%grid_functions(igf)%values(1:gf_ndim,1:npoints_radial) = radial_values(1:gf_ndim,1:npoints_radial)
                deallocate(radial_values)

              end if
            end if

!           printout

            if (rs_input%print_level > 1) then
              write(*, *) 'Grid_function values on master, after gather: '
              write(*, *) 'values: ', i,j, (rs_input%grid_functions(igf)%values(i,j),i=1,gf_ndim)
              if (rs_input%grid_functions(igf)%purpose == 1) then
                write(*, *) 'integrated values: ', i, (rs_input%grid_functions(igf)%integrated_values(i),i=1,gf_ndim)
              end if
            end if
          end if

          deallocate(grid_function_values)
          if (rs_input%grid_functions(igf)%purpose == 1) then
            deallocate(grid_function_integrated_values)
          end if
 
          call dmat%erase()

!
!         save and print data
!

          if (rs_input%grid_functions(igf)%export) then
             if (mpi_is_master()) then
               call store_data(rs_input%grids(ig), rs_input%grid_functions(igf))
             end if
          end if

          if (rs_input%grid_functions(igf)%outpri) then
            if (mpi_is_master()) then
              call write_plot_data_to_DIRAC_output(rs_input%grids(ig), &
                                          rs_input%grid_functions(igf))
            end if
          end if

!
!         clean up
!

          if (rs_input%grid_functions(igf)%allocated_values()) &
                        deallocate(rs_input%grid_functions(igf)%values)
          if (rs_input%grid_functions(igf)%allocated_integrated_values()) &
                        deallocate(rs_input%grid_functions(igf)%integrated_values)
        end if ! if (rs_input%grid_functions(igf)%id_assigned_grid == ig) then

        if (.not. mpi_is_master()) then
          if (allocated(rs_input%grid_functions(igf)%name)) deallocate(rs_input%grid_functions(igf)%name)
          if (allocated(rs_input%grid_functions(igf)%file_inp)) deallocate(rs_input%grid_functions(igf)%file_inp)
          if (allocated(rs_input%grid_functions(igf)%file_out)) deallocate(rs_input%grid_functions(igf)%file_out)
        end if

      end do ! do igf = 1, rs_input%nr_grid_functions

      if (rs_input%grids(ig)%allocated_points()) deallocate(rs_input%grids(ig)%points)
      if (rs_input%grids(ig)%allocated_weights()) deallocate(rs_input%grids(ig)%weights)
      if (rs_input%grids(ig)%allocated_radial_coord()) deallocate(rs_input%grids(ig)%radial_coord)
      if (allocated(rs_input%grids(ig)%file_inp)) deallocate(rs_input%grids(ig)%file_inp)
      if (allocated(rs_input%grids(ig)%file_out)) deallocate(rs_input%grids(ig)%file_out)

    end do ! do ig = 1, rs_input%nr_grids


    do igf = 1, rs_input%nr_grid_functions
      call rs_input%grid_functions(igf)%erase()
    end do
    do ig = 1, rs_input%nr_grids
      call rs_input%grids(ig)%erase()
    end do

  end subroutine

! ------------------------------------------------------------------------------
  subroutine work_on_grid(ndim, nr_points, &
                          mol,&
                          grid,   &
                          grid_points,   &
                          grid_function, &
                          grid_function_values, &
                          dmat_0,     &
                          dmat,     &
                          ao, buffer, &
                          grid_weights, &
                          grid_function_integrated_values)

    integer, intent(in):: ndim
    integer, intent(in):: nr_points
    real(8), intent(in):: grid_points(ndim, nr_points)
    type(mol_t), intent(in):: mol
    type(grid_t), intent(inout):: grid
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    real(8), intent(in), optional:: grid_weights(nr_points)
    real(8), intent(out), optional:: grid_function_integrated_values(:)
    type(dmat_t), intent(in):: dmat_0
    type(dmat_t), intent(in):: dmat
    real(8), intent(inout):: ao(*)
    real(8), intent(inout):: buffer(*)
    real(8):: radial_spacing
    integer:: ipoint

    do ipoint = 1, nr_points
    
      call get_ao(1,                     &
                  (/grid_points(1, ipoint)/),&
                  (/grid_points(2, ipoint)/),&
                  (/grid_points(3, ipoint)/),&
                  ao,                    &
                  buffer)
    
      select case (trim(grid_function%gf_type))
!       fixme:
!       this is a clumsy way to handle grid functions
!       that depend on a single density matrix (dmat_0 .or. dmat)
!       or on multiple density matrices (dmat_0 .and. dmat), as in LAO case;
!       this causes some code duplication in visual_in_point.F90 and should be refactored
        case ("composite")
          call get_quantity_in_point_refactored(dmat_0, dmat,&
                                                ao,buffer,&
                                                grid_function,grid_function_values,&
                                                grid_points(1, ipoint),grid_points(2, ipoint),grid_points(3, ipoint),&
                                                mol,&
                                                ipoint=ipoint)
        case ("response")
          call get_quantity_in_point_refactored(dmat, &
                                                ao,buffer,&
                                                grid_function,grid_function_values,&
                                                grid_points(1, ipoint),grid_points(2, ipoint),grid_points(3, ipoint),&
                                                mol,&
                                                ipoint=ipoint)
        case ("unperturbed")
          call get_quantity_in_point_refactored(dmat_0, &
                                                ao,buffer,&
                                                grid_function,grid_function_values,&
                                                grid_points(1, ipoint),grid_points(2, ipoint),grid_points(3, ipoint),&
                                                mol,&
                                                ipoint=ipoint)
      end select

      if (grid%typ == 2) then
!       radial grids
        grid_function_values(:,ipoint) = &
            grid_function_values(:,ipoint)*grid_weights(ipoint)
        if (grid_function%purpose==1) then
          if (grid%radial_nr_steps > 1) then
            radial_spacing = grid%radial_length/grid%radial_nr_steps
          else
            radial_spacing = grid%radial_length
          end if  
          grid_function_integrated_values(:)= grid_function_integrated_values(:) &
               + grid_function_values(:,ipoint)*radial_spacing
        end if
      else
        if (grid_function%purpose==1) then
!         integration on other grids
          grid_function_integrated_values(:) = grid_function_integrated_values(:) &
                + grid_function_values(:,ipoint)*grid_weights(ipoint)
        end if
      end if
    
    end do

  end subroutine

end module
