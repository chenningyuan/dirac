module visual

! This module contains routines that calculate 
! various quantities on the real-space grids,
! Written by Radovan Bast.
! Refactored by Gosia Olejniczak, spring-summer 2023.
! Full description of features: LINK

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
  use dirac_dmat_factory

! DIRAC modules
  use checkpoint
  use basis_set_datatypes

#ifdef VAR_MPI
  use interface_to_mpi
  use integer_kind_mpilib
#endif


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

! -----------------------------------------------------------------------------
  subroutine run_visual(launch_only)
    logical, intent(in), optional :: launch_only

    type(mol_t):: mol
    type(cmo_t):: mo_data
    type(basis_set_info_t):: ao_basis
    real(8) :: walltime
    logical :: legacy_prep
    type(dmat_t):: dmat_0



!   ---- legacy ----    
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

    logical              :: debug_print
    logical              :: krmcscf_found
    logical              :: natorb_found
    logical, external    :: fndlab

    integer              :: q_max_dim
    integer              :: irank

    integer:: i, j, off_ij, off_iz, off_ij_s, off_ij_e, jrdmo
    type(matrix):: D_mo_dagger, TP, TM, T
    logical, allocatable:: matrix_is_unperturbed(:)
!   ---- end ----    
! -----------------------------------------------------------------------------

    walltime = visual_walltime()

#ifdef VAR_MPI
!
!   data available on input:
!     * information from user input stored in rs_input dt
!     * config information in visual_cfg* variables (legacy)
!
!   purpose:
!     * make sure that all workers have all information to start
!     * if launch_only then wake up nodes, pass initial information,
!       and leave the rest to master 
!       (happens when called from dirnod; on workers),
!       otherwise - calculate (happens when called from dirac; on master)
!    
    call visual_initialize_and_sync_nodes(launch_only)
    if (debug_print) then
      print *, 'VISUAL parallel? ', is_parallel()
      print *, 'size of communicator: ', mpi_get_nr_proc()
      print *, 'rank of the process: ', get_my_rank()
    end if
#endif

    if (mpi_is_master()) then
      write(*, '(/,2x, a)') '*** Visualisation module ***'
      call verify_input_keywords()
    end if


    if (mpi_is_master()) then
!
!     limitations:
!     * calculations are done on one molecular system,
!     * data (molecular information+wave function) is read from DIRAC's CHECKPOINT file

!
!     fixme list: 
!     * ao_basis is kept as a placeholder:
!       now, AOs are still calculated via dirac_ao_eval_init,
!       as it also sets up all offsets used in xcint/density_eval.F90;
!       to fix this, we would need to think of decoupling visual from xcint
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

!   gosia todo
!   open rs-checkpoint (once a separate RS schema is ready)
!   write (lupri,'(/A,I2)')," * Opening checkpoint file for real-space data (VISUAL) "
!   call checkpoint_rs_open
! 
    call visual_driver_refactored(mol,mo_data,dmat_0)

!   call checkpoint_rs_close
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
  subroutine visual_driver_refactored(mol,cmo,dmat_0)
    type(mol_t), intent(in):: mol
    type(cmo_t), intent(in):: cmo
    type(dmat_t), intent(in):: dmat_0
    type(dmat_t):: dmat, dmat_lao

    real(8), allocatable              :: ao(:,:)
    real(8), allocatable              :: buffer(:,:)
    type(grid_t):: grid
    type(grid_t):: this_grid
    type(grid_t), allocatable:: grids(:)
    type(grid_function_t):: this_grid_function
    type(grid_function_t), allocatable:: grid_functions(:)
    real(8), allocatable:: grid_function_values(:,:)
    integer:: ig, igf
    integer:: g_ndim, gf_ndim, npoints
    real(8):: extents(0:1, 3)
    integer:: my_rank
    real(8) :: walltime
    integer:: n_this_proc, idim_grid, idim_gf
    integer:: i,j
    integer:: iz, off_ij
! ------------------------------------------------------------------------------

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
!     todo list:
!     * test whether reshaping to 1-dimensional arrays is better for optimization:
!     grid_points = reshape(rs_input%grids(ig)%points(g_ndim, npoints))
!
!     fixme list:
!     * the grid can be read from file or generated;
!       now, grid generation is done sequentially (on master), it should be done in parallel
!
      if (mpi_is_master()) then
        call get_grid_points(rs_input%grids(ig))
        if (rs_input%print_level > 0) then
          write(*, *) 'Grid points of grid ', ig
          do j = 1, npoints
            write(*, *) (rs_input%grids(ig)%points(i,j),i=1,g_ndim)
            !write(*, '(3es28.18e2)') (rs_input%grids(ig)%points(i,j),i=1,g_ndim)
          end do
        end if
      end if

      if (rs_input%grids(ig)%export) then
!       export this grid to a separate file 
        if (mpi_is_master()) then
          call store_grid(rs_input%grids(ig))
          !if (rs_input%print_level > 0) then
          !  write(*, *) 'Grid ', rs_input%grids(ig)%id, ' exported to file ', rs_input%grids(ig)%file_out
          !end if
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

          ! gosia: it would be cleaner to use "associate", 
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
!         gosia fixme(?)
!         VISUAL uses AO evaluation routines from xcint;
!         which also set up all offsets used in subroutines in xcint/density_eval.F90,
!         called from visual_in_point.F90
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
                                   wf_order=rs_input%grid_functions(igf)%wf_order,&
                                   dmat_ao=dmat, &
                                   rspvec_file=rs_input%grid_functions(igf)%rspvec_file,&
                                   rspvec_id=rs_input%grid_functions(igf)%rspvec_id)


              if (rs_input%grid_functions(igf)%lao > 0) then

                if (visual_cfg_london_skip_kappa) then
                  dmat%elements = 0.0d0
                end if

                if (.not. visual_cfg_london_skip_ro) then
                call dirac_dmat_init('lao-connection-1', cmo, &
                                     wf_order=rs_input%grid_functions(igf)%wf_order,&
                                     rspvec_file=rs_input%grid_functions(igf)%rspvec_file,&
                                     dmat_ao=dmat_lao, &
                                     rspvec_id=rs_input%grid_functions(igf)%rspvec_id)
                dmat%elements = dmat%elements + dmat_lao%elements
                end if

              end if

            end if
!#ifdef VAR_MPI
!            if (is_parallel()) call synchronize_dmat(dmat)
!#endif
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
          end if
#endif

!
!         work on grid
!

          gf_ndim  = rs_input%grid_functions(igf)%get_ndim()

          allocate(grid_function_values(gf_ndim, n_this_proc))
          grid_function_values = 0.0d0
 
          ! gosia todo:
          ! add "visualize"/"integrate" procedures in this call; now - visualization only
          !
!fixme:
          !if (rs_input%grid_functions(igf)%wf_order > 0) then
!         pass dmat:
          call work_on_grid(g_ndim, n_this_proc, &
                            rs_input%grids(ig)%points(1:g_ndim,1:n_this_proc), &
                            rs_input%grid_functions(igf), &
                            grid_function_values, &
                            dmat_0,           &
                            dmat,           &
                            ao, buffer)
          !else
!         !pass dmat_0:
          !call work_on_grid(g_ndim, n_this_proc, &
          !                  rs_input%grids(ig)%points(1:g_ndim,1:n_this_proc), &
          !                  rs_input%grid_functions(igf), &
          !                  grid_function_values, &
          !                  dmat_0,           &
          !                  ao, buffer)
          !end if

          deallocate(ao)
          deallocate(buffer)

          allocate(rs_input%grid_functions(igf)%values(gf_ndim, npoints))
          rs_input%grid_functions(igf)%values = 0.0d0

          !write(*, *) 'AFTER WORK_ON_GRID'
          !!write(*, *) 'values: ', ((grid_function_values(i,j),j=1,n_this_proc),i=1,gf_ndim), get_my_rank()
          !do j = 1, n_this_proc
          !  write(*, '(a,i4,3es28.18e2)') 'values: ', j, (grid_function_values(i,j),i=1,gf_ndim)
          !end do
 
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

          !if (mpi_is_master()) then
          !  write(*, *) 'AFTER GATHERING'
          !  write(*, *) 'gathered values: ', rs_input%grid_functions(igf)%values
          !end if
 
          deallocate(grid_function_values)
          call dmat%erase()

!
!         save data
!

          if (rs_input%grid_functions(igf)%export) then
             if (mpi_is_master()) then
               call store_data(rs_input%grids(ig), rs_input%grid_functions(igf))
             end if
          end if

          !if (mpi_is_master()) then
          !  write(*, *) 'AFTER STORING'
          !  write(*, *) 'stored values?: ', rs_input%grid_functions(igf)%values
          !end if
 
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

  end subroutine


! ------------------------------------------------------------------------------
  subroutine work_on_grid(ndim, nr_points, &
                          grid_points,   &
                          grid_function, &
                          grid_function_values, &
                          dmat_0,     &
                          dmat,     &
                          ao, buffer)

    integer, intent(in):: ndim
    integer, intent(in):: nr_points
    real(8), intent(in):: grid_points(ndim, nr_points)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    type(dmat_t), intent(in):: dmat_0
    type(dmat_t), intent(in):: dmat
    real(8), intent(inout):: ao(*)
    real(8), intent(inout):: buffer(*)
    integer:: ipoint

    if (grid_function%name .eq. 'espe') then

      !do ipoint = 1, nr_points
      !  call get_ao(1,                     &
      !              grid_points(1, ipoint),&
      !              grid_points(2, ipoint),&
      !              grid_points(3, ipoint),&
      !              ao,                    &
      !              buffer)
      !end do

      !call get_esp_functions_in_points(get_my_rank(),&
      !                                 dmat_0,&
      !                                 ao,buffer,&
      !                                 ndim,nr_points,&
      !                                 grid_points,&
      !                                 grid_function,grid_function_values)
      !call TEMP_test_values_in_points('grid_function_values in point', &
      !                       nr_points, grid_points, v1=grid_function_values)

    else
      do ipoint = 1, nr_points
      
        call get_ao(1,                     &
                    (/grid_points(1, ipoint)/),&
                    (/grid_points(2, ipoint)/),&
                    (/grid_points(3, ipoint)/),&
                    ao,                    &
                    buffer)
      
        call get_quantity_in_point_refactored(dmat_0, dmat,&
                                              ao,buffer,&
                                              grid_function,grid_function_values,&
                                              grid_points(1, ipoint),grid_points(2, ipoint),grid_points(3, ipoint),&
                                              ipoint)
      
        !call TEMP_test_values_in_point('grid_function_values in point', &
        !                       ipoint, grid_points(1:3, ipoint), v3=grid_function_values(1:3,ipoint))
      
      end do

    end if

  end subroutine work_on_grid



























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


  subroutine TEMP_test_values_in_point(text, ipoint, r3, v1, v3, v9)
    character(*), intent(in) :: text
    integer, intent(in) :: ipoint
    real(8), intent(in) :: r3(3)
    real(8), intent(in), optional :: v1
    real(8), intent(in), optional :: v3(3)
    real(8), intent(in), optional :: v9(9)

    !if (present(v1)) then
    !  write(*, '(A,2I6,3E16.8,ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                       r3(1:3), &
    !                                       v1
    !end if

    !if (present(v3)) then
    !  write(*, '(A,2I6,3E16.8,3ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                       r3(1:3), &
    !                                       v3(1:3)
    !end if

    !if (present(v9)) then
    !  write(*, '(A,2I6,3E16.8,9ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                       r3(1:3), &
    !                                       v9(1:9)
    !end if

  end subroutine

  subroutine TEMP_test_values_in_points(text, nr_points, r3, v1, v3, v9)
    character(*), intent(in) :: text
    integer, intent(in) :: nr_points
    real(8), intent(in) :: r3(3,nr_points)
    real(8), intent(in), optional :: v1(nr_points)
    real(8), intent(in), optional :: v3(3,nr_points)
    real(8), intent(in), optional :: v9(9,nr_points)
    integer :: ipoint

    !if (present(v1)) then
    !  do ipoint = 1, nr_points
    !    write(*, '(A,2I6,3E16.8,ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                         r3(1:3,ipoint), &
    !                                         v1(ipoint)
    !  end do
    !end if

    !if (present(v3)) then
    !  do ipoint = 1, nr_points
    !    write(*, '(A,2I6,3E16.8,3ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                          r3(1:3,ipoint), &
    !                                          v3(1:3,ipoint)
    !  end do
    !end if

    !if (present(v9)) then
    !  do ipoint = 1, nr_points
    !    write(*, '(A,2I6,3E16.8,9ES30.20E3)') text, get_my_rank(), ipoint, &
    !                                          r3(1:3,ipoint), &
    !                                          v9(1:9,ipoint)
    !  end do
    !end if

  end subroutine


end module
