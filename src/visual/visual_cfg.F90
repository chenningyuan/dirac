module visual_cfg

! change: visual_datatypes->rs_datatypes and visual_cfg -> rs_cfg
  use visual_datatypes, only : rs_input_t, grid_t, grid_function_t
  implicit none

  type(rs_input_t), public:: rs_input


  !public setup_grid
  !public setup_grid_function
  public grid_input_cfg
  public grid_function_input_cfg
  public verify_grid_input_cfg
  public verify_grid_function_input_cfg

  public import_dirac_cb_variables

!  public import_dirac_cb_variables
#ifdef VAR_MPI
  public sync_visual_cfg
#endif

   save


!  stays; remember to bcast!!!
   real(8), public, allocatable :: visual_cfg_occupation(:, :)
   logical, public              :: visual_cfg_use_orbital_string = .false.

   real(8), public:: visual_cfg_cval

!

   integer,      public, parameter:: max_batch_length = 500
   integer,      public, parameter:: max_para = 100
   integer,      public            :: visual_cfg_nr_dmat = 0 !FIXME
   character(6), public            :: visual_cfg_dmat_file(max_para)
   integer,      public            :: visual_cfg_dmat_file_record(max_para)
   integer,      public            :: visual_cfg_ref_nucleus(max_para)
   integer,      public            :: visual_cfg_property(max_para)
   logical,      public            :: visual_cfg_skip(max_para) = .false.
   integer,      public            :: visual_cfg_property_to_matrix(max_para) = 0
   integer,      public            :: visual_cfg_irep_xvector(max_para)
   logical,      public            :: visual_cfg_out_file_custom = .false.
   logical,      public            :: visual_cfg_out_scalar_file_custom = .false.
   logical,      public            :: visual_cfg_out_vector_file_custom = .false.
   logical,      public            :: visual_cfg_out_2ndtensor_file_custom = .false.
   logical,      public            :: visual_cfg_add_grid_points_grad_ed = .false.
   logical,      public            :: visual_cfg_add_grid_points_grad_ed_scaled = .false.
   character(60), public           :: visual_cfg_out_scalar_file_name
   character(60), public           :: visual_cfg_out_vector_file_name
   character(60), public           :: visual_cfg_out_2ndtensor_file_name
   logical,      public            :: visual_cfg_london_skip_cube = .false.
   logical,      public            :: visual_cfg_london_skip_vtk  = .false.

   logical,      public            :: visual_cfg_global_refactor = .false.


!
!  options for exporting real-space data to files
!
! WIP
!  character(len = :), public               :: visual_cfg_plot_file_formats(2) = (/'txt', 'cube', 'csv', 'hdf5'/)
   character(len = 6), public               :: visual_cfg_txt_plot_file_domains(4) = &
                                              [character(len = 6):: '2D', '3D', 'line', 'radial']

   logical,      public                     :: visual_cfg_plot_file_custom = .false.
   character(len = :), allocatable, public  :: visual_cfg_plot_file_format
   character(len = :), allocatable, public  :: visual_cfg_plot_domain
   character(len = :), allocatable, public  :: visual_cfg_plot_property

!  options for exporting to txt files
!  by default, export real-space data to txt files and use legacy file names ("plot.3d.scalar", etc.)
   logical,      public                     :: visual_cfg_txt_plot_file_use_legacy = .true.
   character(len = :), allocatable, public  :: visual_cfg_txt_plot_file_name
   character(len = :), allocatable, public  :: visual_cfg_txt_plot_file_float_format

!  options for exporting to cube files (applies to scalar fields on 3D cube grids)
!  by default, use old export to cube files
   logical,      public                     :: visual_cfg_cube_plot_file_use_legacy = .true.
   character(len = :), allocatable, public  :: visual_cfg_cube_plot_file_name
   character(len = :), allocatable, public  :: visual_cfg_cube_plot_file_float_format

!  options for grid manipulation
   integer, public:: visual_cfg_nr_grid_files = 1

   real(8), public  :: visual_cfg_txt_print_threshold = 1.0e-8

!
!  option for input data files
!
   logical,      public                     :: visual_cfg_inp_data_files_use_legacy= .false.
   integer,      public                     :: visual_cfg_nr_plot_files = 0

   integer, public:: iq_elf      =  3
   integer, public:: iq_j        =  4
   integer, public:: iq_rotj     =  5
   integer, public:: iq_divj     =  6
   integer, public:: iq_s        =  7
   integer, public:: iq_rots     =  8
   integer, public:: iq_divs     =  9
   integer, public:: iq_jdia     = 10
   integer, public:: iq_edipx    = 11
   integer, public:: iq_edipy    = 12
   integer, public:: iq_edipz    = 13
   integer, public:: iq_bdipx    = 14
   integer, public:: iq_bdipy    = 15
   integer, public:: iq_bdipz    = 16
   integer, public:: iq_ndipx    = 17
   integer, public:: iq_ndipy    = 18
   integer, public:: iq_ndipz    = 19
   integer, public:: iq_bdipxdia = 20
   integer, public:: iq_bdipydia = 21
   integer, public:: iq_bdipzdia = 22
   integer, public:: iq_ndipxdia = 23
   integer, public:: iq_ndipydia = 24
   integer, public:: iq_ndipzdia = 25
   integer, public:: iq_esp      = 26
   integer, public:: iq_espe     = 27
   integer, public:: iq_espn     = 28
   integer, public:: iq_esprho   = 29
   integer, public:: iq_esperho  = 30
   integer, public:: iq_espnrho  = 31
   integer, public:: iq_kin      = 32
   integer, public:: iq_kin_ls   = 33
   integer, public:: iq_kin_sl   = 34
   integer, public:: iq_kin_tau  = 35
   integer, public:: iq_kin_lap  = 36
   integer, public:: iq_kin_nr   = 37
   integer, public:: iq_jx       = 38
   integer, public:: iq_jy       = 39
   integer, public:: iq_jz       = 40
   integer, public:: iq_tcos     = 41
   integer, public:: iq_tsin     = 42   
   integer, public:: iq_gradj    = 48
   

   logical, public:: visual_cfg_print_data_to_dirac_output = .false.

!                     visualization in 2d
   logical, public:: visual_cfg_2d             = .false.
   real(8), public:: visual_cfg_2d_p_origin(3) = 0.0d0
   real(8), public:: visual_cfg_2d_p_right(3)  = 0.0d0
   real(8), public:: visual_cfg_2d_p_top(3)    = 0.0d0
   integer, public:: visual_cfg_2d_nr_right    = 0
   integer, public:: visual_cfg_2d_nr_top      = 0
   character(len = 6), public:: visual_cfg_2d_data_types(4) = &
                       [character(len = 6):: 'scalar', 'vector', 'tensor', 'flow']

!                     integration in 2d
   logical, public:: visual_cfg_2d_integration             = .false.
   real(8), public:: visual_cfg_2d_integration_p_origin(3) = 0.0d0
   real(8), public:: visual_cfg_2d_integration_p_right(3)  = 0.0d0
   real(8), public:: visual_cfg_2d_integration_p_top(3)    = 0.0d0
   integer, public:: visual_cfg_2d_integration_nr_right    = 0
   integer, public:: visual_cfg_2d_integration_nr_top      = 0
   integer, public:: visual_cfg_2d_integration_order       = 0

!                     visualization in 3d
   logical, public:: visual_cfg_3d         = .false.
   logical, public:: visual_cfg_3d_fast    = .false.
   integer, public:: visual_cfg_ncube(1:3) = 0
   real(8), public:: visual_cfg_add_3d     = 4.0d0
   logical, public:: visual_cfg_3d_gridimp = .false.
   logical, public:: visual_cfg_3d_readjb  = .false.
   character(80), public:: visual_cfg_3d_gridfil
   character(80), public:: visual_cfg_3d_grid_units
   character(80), public:: visual_cfg_3d_jbfile
   character(len = 10), public:: visual_cfg_3d_data_types(5) = &
                       [character(len = 10)::'scalar', 'vector', 'tensor', 'cube', 'vector_vtk']

!                     integration in 3d
   logical, public:: visual_cfg_3d_integration = .false.

!                     visualization along line
   logical, public:: visual_cfg_line          = .false.
   real(8), public:: visual_cfg_line_from(3)  = 0.0d0
   real(8), public:: visual_cfg_line_to(3)    = 0.0d0
   integer, public:: visual_cfg_line_nr_steps = 0
   character(len = 6), public:: visual_cfg_line_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']

!                                  densities for a list of points
   logical, public              :: visual_cfg_list              = .false.
   integer, public              :: visual_cfg_nr_points_in_list = 0
   real(8), public, allocatable:: visual_cfg_xyz_list(:, :)
   character(len = 6), public:: visual_cfg_list_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']


!                     radial distributions
   logical, public:: visual_cfg_radial          = .false.
   real(8), public:: visual_cfg_radial_from(3)  = 0.0d0
   integer, public:: visual_cfg_radial_nr_steps = 0
   real(8), public:: visual_cfg_radial_length   = 0
   character(len = 6), public:: visual_cfg_radial_data_types(3) = &
                       [character(len = 6)::'scalar', 'vector', 'tensor']
   
   real(8),       public:: visual_cfg_scale              = 1.0d0
   real(8),       public:: visual_cfg_cartesian_power(3) = 0.0   
   real(8),       public:: visual_cfg_gauge_origin(3)    = 0.0d0
   real(8),       public:: visual_cfg_wave_vector(3)     = 0.0d0
   real(8),       public:: visual_cfg_freq               = 0.0d0   
   real(8),       public:: visual_cfg_pol_vector(3)      = 0.0d0   
   logical,       public:: visual_cfg_force_small_ao     = .false.

   logical,       public:: visual_cfg_nics              = .false.
   real(8),       public:: visual_cfg_nics_origin(3)    = 0.0d0

   logical,      public:: visual_cfg_london               = .false.
   logical,      public:: visual_cfg_london_skip_kappa    = .false.
   logical,      public:: visual_cfg_london_skip_direct   = .false.
   logical,      public:: visual_cfg_london_skip_ro       = .false.
   logical,      public:: visual_cfg_london_none          = .false.
   character(1), public:: visual_cfg_london_component

   integer,      public:: visual_cfg_irep_conmat = 0
   integer,      public:: visual_cfg_itim_conmat = 0


   logical,      public:: visual_cfg_tensor          = .false.

   private

contains

  elemental subroutine str2int(str,int1,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int1
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int1
  end subroutine str2int

  elemental subroutine str2flt(str,flt,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    real(8),intent(out)         :: flt
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  flt
  end subroutine str2flt


   function word_count(s,sep)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len=1), intent(in) :: sep
      integer                  :: word_count
!     --------------------------------------------------------------------------
      integer                  :: i
      logical                  :: is_blank
!     --------------------------------------------------------------------------

      word_count = 0

      if (len(s) <= 0) return

      is_blank = .true.

      do i = 1, len(s)
         if (s(i:i) == sep) then
            is_blank = .true.
         else if (is_blank) then
            word_count = word_count + 1
            is_blank = .false.
         end if
      end do

   end function

   subroutine split_string(str, nr_words, word_array, separator)

!     --------------------------------------------------------------------------
      integer, intent(in)   :: nr_words
      character(len=*), intent(in) :: str
      character(len=*), intent(inout):: word_array(nr_words)
      character(len=1), optional, intent(in) :: separator
!     --------------------------------------------------------------------------
      integer                  :: i,j
      integer,allocatable::pos(:)
      integer                  :: len_str
      character(len=1):: sep
!     --------------------------------------------------------------------------

      len_str = len_trim(str)
      if (len_str <= 0) return

      if (present(separator)) then
        sep=separator
      else
        sep=' '
      end if

      allocate(pos(nr_words-1))
      pos = 0

      j = 0
      do i = 1, len_str
         if (str(i:i) == sep) then
            j = j + 1
            pos(j) = i
         end if
      end do

      word_array(1)(1:) = str(1:pos(1)-1)
      do i = 1, nr_words-2
        word_array(i+1)(1:) = str(pos(i)+1:pos(i+1)-1)
      end do
      word_array(nr_words)(1:) = str(pos(nr_words-1)+1:)

      deallocate(pos)

   end subroutine

  subroutine grid_input_cfg(key, val, grid, id_grid)
#include "codata.h"
    character(len=*), intent(in):: key, val
    type(grid_t),intent(inout)::grid
    integer, intent(in)::id_grid
    character(len=len_trim(val))::temp
    integer::stat, np
    real(8)::s
    integer::nr_words, i
    character(len=10),allocatable::chvalues(:)

!   gosia: possibly do sth tht mimics python parser here:
!   write(*, *) ' ".GRIDS" arguments available in this version:'

    select case (key)
      case ('id')
        call str2int(val, grid%id, stat)
      case ('ndim')
        call str2int(val, grid%ndim, stat)
        grid%ndim_known = .true.
      case ('npoints')
        temp = val
        call trim_ch(temp, val,  '[')
        call trim_ch(temp, temp, ']')
        if (index(val, ',') > 0) then
          nr_words = word_count(temp, ',')
          allocate(chvalues(nr_words))
          call split_string(temp, nr_words, chvalues, separator=',')
          do i = 1, nr_words
            call str2int(trim(chvalues(i)), np, stat)
            grid%nr_points(i) = np
            grid%nr_points_known(i) = .true.
          end do
          deallocate(chvalues)
        else
          call str2int(temp, grid%nr_points(1), stat)
          grid%nr_points_known(1) = .true.
        end if
      case ('spacing')
        temp = val
        call trim_ch(temp, val,  '[')
        call trim_ch(temp, temp, ']')
        if (index(val, ',') > 0) then
          nr_words = word_count(temp, ',')
          allocate(chvalues(nr_words))
          call split_string(temp, nr_words, chvalues, separator=',')
          do i = 1, nr_words
            call str2flt(trim(chvalues(i)), s, stat)
            ! scale to a.u.
            grid%spacing(i) = s/xtang
            grid%spacing_known(i) = .true.
          end do
          deallocate(chvalues)
        else
          call str2flt(temp, grid%spacing(1), stat)
          grid%spacing_known(1) = .true.
          ! scale to a.u.
          grid%spacing(1) = grid%spacing(1)/xtang
        end if
      case ('input')
        if (val == 'create') then
          grid%input = 0
        else if (val == 'import') then
          grid%input = 1
        else
          write(*, *) 'wrong argument; argument to "input" must be one of "create"/"import"'
          stop
        end if
      case ('typ')
        call str2int(val, grid%typ, stat)
      case ('weights')
        if (val == 'no') then
          grid%weights = 0
        else if (val == 'yes') then
          grid%weights = 1
        else
          write(*, *) 'wrong argument; argument to "weights" must be one of "yes"/"no"}'
          stop
        end if
      case ('margin')
        call str2flt(val, grid%margin, stat)
        ! scale to a.u.
        grid%margin = grid%margin/xtang
      case ('export')
        if (val == 'no') then
          grid%export = .false.
        else if (val == 'yes') then
          grid%export = .true.
        else
          write(*, *) 'wrong argument; argument to "weights" must be one of "yes"/"no"}'
          stop
        end if
      case ('file_inp')
        grid%file_inp = val
      case ('file_out')
        grid%file_out = val
      case default
        print *, 'unknown key'
        stop
    end select
    if (grid%id == 0) grid%id = id_grid
  end subroutine grid_input_cfg

  subroutine verify_grid_input_cfg(grid)
    type(grid_t),intent(inout)::grid
    integer::i
    logical:: err

    if (grid%input == 1) then
      ! when reading grid from file
      write(*, *) 'WARNING: grid information will be read from file'
      write(*, *) 'These input options will be ignored (if present): ndim, npoints, spacing, margin.'

    else if (grid%input == 0) then
      ! when creating grid
      if (grid%ndim == 0) then
        write(*, *) 'ERROR: grid "ndim" information is missing.'
        stop
      end if
      do i = 1, grid%ndim
        if (grid%nr_points(i) ==0 .and. grid%spacing(i) ==0) then
          write(*, *) 'ERROR: grid "npoints" or "spacing" information is missing for dimension ', i
          stop
        end if
      end do
    end if

    do i = 1, grid%ndim
      if (grid%nr_points(i) /=0 .and. grid%spacing(i) /=0) then
        write(*, *) 'ERROR: exclusive options: "npoints" and "spacing".'
        stop
      end if
    end do

  end subroutine verify_grid_input_cfg


  subroutine verify_grid_function_input_cfg(grid_function)
    type(grid_function_t),intent(in)::grid_function

    if (grid_function%id_assigned_grid == 0) then
      write(*, *) 'ERROR: missing information about the grid ("id_grid").'
      stop
    end if

    if (grid_function%input == 1) then
      if (len_trim(grid_function%file_inp) == 0) then
        write(*, *) 'ERROR: missing file for grid_function import.'
        stop
      end if
    !else if (grid_function%input == 0) then
    !  call verify_grid_function_predefined(grid_function)
    end if

  end subroutine verify_grid_function_input_cfg


  subroutine grid_function_input_cfg(key, val, grid_function, id_grid_function)
    character(len=*), intent(in):: key, val
    type(grid_function_t),intent(inout)::grid_function
    integer, intent(in)::id_grid_function
    character(len=len(key))::temp
    integer::stat

    select case (key)
      case ('name')
        grid_function%name = val
        ! gosia fixme: this call will be moved up to verify_grid_function_input_cfg subroutine
        call verify_grid_function_predefined(grid_function)
      case ('id_grid')
        call str2int(val, grid_function%id_assigned_grid, stat)
      case ('input')
        if (val == 'calculate') then
          grid_function%input = 0
        else if (val == 'import') then
          grid_function%input = 1
        else
          write(*, *) 'wrong argument; argument to "input" must be one of "create"/"import"'
          stop
        end if
      case ('rspvec_file')
!       fixme: verify, generalize (now: expected 'RESPONSE.h5' file; linear response only)
        grid_function%rspvec_file = val
        grid_function%wf_order = 1
      case ('rspvec_id')
!       fixme
          if (grid_function%wf_order /= 1) then
            write(*, *) 'rsp is not set up'
            stop
          end if
          call str2int(val, grid_function%rspvec_id, stat)
      case ('lao')
        if (grid_function%wf_order /= 1) then
          write(*, *) 'LAO activated for linear response only'
          stop
        end if
        call str2int(val, grid_function%lao, stat)
      case ('purpose')
        if (val == 'visualization') then
          grid_function%purpose = 0
        else if (val == 'integration') then
          grid_function%purpose = 1
        end if
      case ('export')
        if (val == 'no') then
          grid_function%export = .false.
        else if (val == 'yes') then
          grid_function%export = .true.
        else
          write(*, *) 'wrong argument; argument to "weights" must be one of "yes"/"no"}'
          stop
        end if
      case ('file_inp')
        grid_function%file_inp = val
      case ('file_out')
        grid_function%file_out = val
      case default
        print *, 'unknown key'
        stop
    end select

    call modify_grid_function_due_to_method(grid_function)

    if (grid_function%id == 0) grid_function%id = id_grid_function

  end subroutine


! ------------------------------------------------------------------------------
  subroutine verify_grid_function_predefined(grid_function)
    class(grid_function_t), intent(inout):: grid_function
!
!   gosia, summer 2023:
!   * initial setup of grid_function based on its name;
!   * grid_function%name also defines a label on HDF5
!

    select case (grid_function%name)
      case ('density', &
            'gamma5', &
            'esp', &
            'espe', &
            'espn', &
            'esprho', &
            'esperho', &
            'espnrho')
        grid_function%wf_order = 0
        grid_function%tensor_order = 0
        grid_function%max_ao_g_order = 0
        grid_function%max_ao_m_order = 0
      case ('rdg', &
            'elf', &
            'kin', &
            'kinls', &
            'kinsl', &
            'kintau', &
            'divj', &
            'rotj', &
            'gradj',&
            'divs') 
        grid_function%wf_order = 0
        grid_function%tensor_order = 0
        grid_function%max_ao_g_order = 1
        grid_function%max_ao_m_order = 0
      case ('s', &
            'j')
        grid_function%wf_order = 0
        grid_function%tensor_order = 1
        grid_function%max_ao_g_order = 0
        grid_function%max_ao_m_order = 0
      case ('density_gradient',&
            'rots')
        grid_function%wf_order = 0
        grid_function%tensor_order = 1
        grid_function%max_ao_g_order = 1
        grid_function%max_ao_m_order = 0
      case ('density_hessian')
        grid_function%wf_order = 0
        grid_function%tensor_order = 2
        grid_function%max_ao_g_order = 2
        grid_function%max_ao_m_order = 0
      case ('ed_laplacian',&
            'ed_sign_lambda2',&
            'kinlap',&
            'kinnr') 
        grid_function%wf_order = 0
        grid_function%tensor_order = 0
        grid_function%max_ao_g_order = 2
        grid_function%max_ao_m_order = 0
      case default
        write(*, *) 'unknown grid function name. stop'
        stop
    end select

  end subroutine


! ------------------------------------------------------------------------------
  subroutine modify_grid_function_due_to_method(grid_function)
    class(grid_function_t), intent(inout):: grid_function

    if (visual_cfg_london .and. grid_function%max_ao_m_order < 1) then
      grid_function%max_ao_m_order = 1
    end if

    !if (levyle) then
    !  if ((grid_function%name == 'j') .or. &
    !      (grid_function%name == 'ndipx') .or. &
    !      (grid_function%name == 'ndipy') .or. &
    !      (grid_function%name == 'ndipz') .or. &
    !      (grid_function%name == 'bdipx') .or. &
    !      (grid_function%name == 'bdipy') .or. &
    !      (grid_function%name == 'bdipz')) then
    !    if (grid_function%max_ao_g_order < 1) then
    !      grid_function%max_ao_g_order = 1
    !    end if
    !  end if
    !end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine should_depend_on_rsp_data(grid_function)
    class(grid_function_t), intent(inout):: grid_function

    if (grid_function%wf_order == 1) then
      !
    end if

  end subroutine


! strings

  subroutine trim_ch(s_out, s_inp, ch)

!     --------------------------------------------------------------------------
      character(*), intent(out) :: s_out
      character(*), intent(in) :: s_inp
      character(1), intent(in) :: ch
!     --------------------------------------------------------------------------
      integer                  :: i,j,k
!     --------------------------------------------------------------------------

      if (len(s_inp) <= 0) return

      k = 0
      do i = 1, len(s_inp)
        if (s_inp(i:i) .ne. ch) then
          k = k + 1
          s_out(k:k) = s_inp(i:i)
        end if
      end do
      s_out = s_out(1:k)

   end subroutine


  subroutine import_dirac_cb_variables()
#include "dcbgen.h"
   visual_cfg_cval = cval  ! speed of light used in Dirac                             
  end subroutine 

#ifdef VAR_MPI
  subroutine sync_visual_cfg()
    !
    ! this subroutine syncs data stored in rs_input
    ! if any new variable has been added to this dt, please also add it here
    !
    use visual_mpi
    use interface_to_mpi
    use integer_kind_mpilib

    call interface_mpi_bcast_l0(rs_input%outpri,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%print_level,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%nr_grids,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%nr_grid_functions,1,irank_master,global_communicator)
    if (.not.mpi_is_master()) then
      ! the allocation on master happens during DIRAC input reading (input/input_reader_sections.F90)
      allocate(rs_input%grids(rs_input%nr_grids))
      allocate(rs_input%grid_functions(rs_input%nr_grid_functions))
    end if

    ! legacy and variables from DIRAC's common blocks
    call interface_mpi_bcast_l0(visual_cfg_global_refactor,1,irank_master,global_communicator)
    call interface_mpi_bcast(visual_cfg_cval,1,irank_master,global_communicator)

    !if (rs_input%print_level > 0) call logprint
  end subroutine sync_visual_cfg


  subroutine logprint
    use visual_mpi
    write(*, *) 'Data on node ', get_my_rank(), 'print_level =       ', rs_input%print_level
    write(*, *) 'Data on node ', get_my_rank(), 'nr_grids            ', rs_input%nr_grids
    write(*, *) 'Data on node ', get_my_rank(), 'nr_grid_functions   ', rs_input%nr_grid_functions
    write(*, *) 'Data on node ', get_my_rank(), 'outpri              ', rs_input%outpri
    !write(*, *) 'Data on node ', get_my_rank(), 'allocated_grids:    ', rs_input%allocated_grids()
  end subroutine logprint

#endif

end module
