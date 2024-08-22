module visual_cfg

  use visual_datatypes, only : rs_input_t, rs_wf_t, grid_t, grid_function_t
  use visual_utils

  implicit none

  type(rs_input_t), public:: rs_input
  type(rs_wf_t), public:: rs_wf

  public grid_input_cfg
  public grid_user_input_cfg
  public grid_function_input_cfg
  public verify_grid_input_cfg
  public verify_grid_function_input_cfg

#ifdef VAR_MPI
  public sync_visual_cfg
#endif

   private

contains

! ------------------------------------------------------------------------------
  subroutine grid_input_cfg(key, val, grid, id_grid)
    character(len=*), intent(in):: key, val
    type(grid_t),intent(inout)::grid
    integer, intent(in)::id_grid
    character(len=len_trim(val))::temp
    integer::stat, np
    real(8)::s
    integer::nr_words, i
    character(len=10),allocatable::chvalues(:)

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
!           units: a.u.
            call str2flt(trim(chvalues(i)), s, stat)
            grid%spacing(i) = s
            grid%spacing_known(i) = .true.
          end do
          deallocate(chvalues)
        else
!         units: a.u.
          call str2flt(temp, grid%spacing(1), stat)
          grid%spacing_known(1) = .true.
        end if
      case ('input')
        if (lowercase(val) == 'create') then
          grid%input = 0
        else if (lowercase(val) == 'import') then
          grid%input = 1
        else if (lowercase(val) == 'list') then
          grid%input = 2
        else if (lowercase(val) == 'line') then
          grid%input = 3
        else if (lowercase(val) == '2d') then
          grid%input = 4
        else if (lowercase(val) == '3d') then
          grid%input = 5
        else if (lowercase(val) == 'radial') then
          grid%input = 6
        else
          write(*, *) 'wrong argument; argument to "input" must be one of "create"/"import"/"list"/"line"/"2d"/"3d"/"radial"'
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
!       units: a.u.
        call str2flt(val, grid%margin, stat)
      case ('export')
        if (lowercase(val) == 'no') then
          grid%export = .false.
        else if (lowercase(val) == 'yes') then
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

! ------------------------------------------------------------------------------
  subroutine grid_user_input_cfg(grid, domain, user_input_lines)
    type(grid_t),intent(inout)::grid
    character(*), intent(in):: domain
    character(*), dimension(:), intent(in):: user_input_lines
    integer::stat, np
    real(8)::s
    integer::nr_words, i,j
    character(len=len_trim(user_input_lines(1))),allocatable::chvalues(:)

    select case (domain)
      case ('list')
        grid%nr_points_total = size(user_input_lines)
        if (index(user_input_lines(1), ',') > 0) then
          nr_words = word_count(user_input_lines(1), ',')
          grid%ndim = nr_words
        end if
        allocate(grid%points(grid%ndim, grid%nr_points_total))
        grid%points = 0.0d0
        do i = 1, grid%nr_points_total
          allocate(chvalues(grid%ndim))
          call split_string(user_input_lines(i),grid%ndim,chvalues,separator=',')
          do j = 1, grid%ndim
            call str2flt(trim(chvalues(j)), s, stat)
            grid%points(j,i) = s
          end do
          deallocate(chvalues)
        end do
        grid%nr_points_total_known = .true.

      case ('line')
        grid%nr_points_total_known = .true.

      case ('2d')
        grid%ndim = 2
        grid%nr_points_total_known = .true.

      case ('3d')
        grid%ndim = 3
        grid%typ = 1
        allocate(chvalues(grid%ndim))
        call split_string(user_input_lines(1),grid%ndim,chvalues,separator=',')
        grid%nr_points_total = 1
        do i = 1, grid%ndim
          call str2int(trim(chvalues(i)), np, stat)
          grid%nr_points(i) = np
          grid%nr_points_known(i) = .true.
          grid%nr_points_total = grid%nr_points_total*grid%nr_points(i)
        end do
        deallocate(chvalues)
        grid%nr_points_total_known = .true.


      case ('radial')

        grid%ndim = 3
        !allocate(chvalues(grid%ndim))
        !call split_string(user_input_lines(1),grid%ndim,chvalues,separator=',')
        !do j = 1, grid%ndim
        !  call str2flt(trim(chvalues(j)), s, stat)
        !  tmp_origin(j) = s
        !end do
        !deallocate(chvalues)
        !call str2flt(trim(user_input_lines(2)), tmp_length, stat)
        !call str2int(trim(user_input_lines(2)), tmp_nr_steps, stat)
        !call generate_radial_grid(grid, tmp_origin, tmp_length, tmp_nr_steps)
        grid%nr_points_total_known = .true.

      case default
        print *, 'unknown key'
        stop
    end select
  end subroutine


! ------------------------------------------------------------------------------
  subroutine verify_grid_input_cfg(grid)
    type(grid_t),intent(inout)::grid
    integer::i
    logical:: err

    if (grid%input == 1) then
      write(*, *) 'VISUAL: Grid information will be read from file'
      write(*, *) '  WARNING: These input options will be ignored (if present):'//&
                  ' ndim, npoints, spacing, margin.'

    else if (grid%input == 0) then
      write(*, *) 'VISUAL: Grid will be created'
      if (grid%ndim == 0) then
        write(*, *) ' ERROR: grid "ndim" information is missing.'
        stop
      end if
      do i = 1, grid%ndim
        if (grid%nr_points(i) ==0 .and. grid%spacing(i) ==0) then
          write(*, *) ' ERROR: grid "npoints" or "spacing" information is missing for dimension ', i
          stop
        else if (grid%nr_points(i) /=0 .and. grid%spacing(i) /=0) then
          write(*, *) ' ERROR: grid exclusive options: "npoints" and "spacing".'
        end if
      end do

    else
      write(*, *) 'VISUAL: Grid will be created according to information passed in a separate line'
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine verify_grid_function_input_cfg(grid_function)
    type(grid_function_t),intent(in)::grid_function

    if (grid_function%id_assigned_grid == 0) then
      write(*, *) 'Input error. Missing information about the grid.'
      stop
    end if

    if (grid_function%input == 1) then
      if (len_trim(grid_function%file_inp) == 0) then
        write(*, *) 'Input error: Missing file for grid_function import.'
        stop
      end if
    end if

    if (grid_function%rspvec_id /= 0) then
      if (len_trim(grid_function%rspvec_file) == 0) then
        write(*, *) 'Input error: Missing response file for grid_function.'
        stop
      end if
    end if

    if (grid_function%lao /= 0) then
      if (grid_function%wf_order /= 1) then
        write(*, *) 'Input error. LAOs can be used for linear response properties only.'
        stop
      end if
      if (len_trim(grid_function%rspvec_file) == 0) then
        write(*, *) 'Input error: Requested LAOs, but no information '// &
                    'about response file for grid_function.'
        stop
      end if
    end if

    if (any(grid_function%lao_skip /= 0)) then
      if (grid_function%lao == 0) then
        write(*, *) 'Input error: "lao_skip" option presented, but LAOs are not requested.'
        stop
      end if
    end if

  end subroutine verify_grid_function_input_cfg


! ------------------------------------------------------------------------------
  subroutine grid_function_input_cfg(key, val, grid_function, id_grid_function)
    character(len=*), intent(in):: key, val
    type(grid_function_t),intent(inout)::grid_function
    integer, intent(in)::id_grid_function
    character(len=len_trim(val))::temp
    integer::stat, np
    real(8)::s,fac
    integer::nr_words, i
    character(len=10),allocatable::chvalues(:)

    select case (key)
      case ('name')
        grid_function%name = lowercase(val)
        call verify_grid_function_predefined(grid_function)
      case ('id_grid')
        call str2int(val, grid_function%id_assigned_grid, stat)
      case ('input')
        if (lowercase(val) == 'calculate') then
          grid_function%input = 0
        else if (lowercase(val) == 'import') then
          grid_function%input = 1
        else
          write(*, *) 'Input error. Wrong argument: '// &
                      'argument to "input" must be one of "calculate"/"import"'
          stop
        end if
      case ('rspvec_file')
!       fixme: expected 'RESPONSE.h5' file; linear response only
        grid_function%rspvec_file = val
        grid_function%wf_order = 1
      case ('rspvec_id')
          call str2int(val, grid_function%rspvec_id, stat)
      case ('nuc')
          call str2int(val, grid_function%nuc, stat)
      case ('lao')
        call str2int(val, grid_function%lao, stat)
      case ('lao_skip')
        temp = val
        call trim_ch(temp, val,  '[')
        call trim_ch(temp, temp, ']')
        if (index(val, ',') > 0) then
          nr_words = word_count(temp, ',')
          if (nr_words .ne. 3) then
            write(*,*) 'Input error. More than three values are in the "lao_skip".'
            stop
          end if
          allocate(chvalues(nr_words))
          call split_string(temp, nr_words, chvalues, separator=',')
          do i = 1, nr_words
            call str2int(trim(chvalues(i)), np, stat)
            grid_function%lao_skip(i) = np
          end do
          deallocate(chvalues)
        else
          write(*,*) 'Input error. "lao_skip" expects an array of integer values.'
          stop
        end if
      case ('purpose')
        if (lowercase(val) == 'visualization') then
          grid_function%purpose = 0
        else if (lowercase(val) == 'integration') then
          grid_function%purpose = 1
        else
          write(*, *) 'Input error. Wrong argument: '// &
                      'argument to "purpose" must be one of '//&
                      '"visualization"/"integration".'
          stop
        end if
      case ('export')
        if (lowercase(val) == 'no') then
          grid_function%export = .false.
        else if (lowercase(val) == 'yes') then
          grid_function%export = .true.
        else
          write(*, *) 'Input error. Wrong argument: '// &
                      'argument to "export" must be one of "yes"/"no".'
          stop
        end if
      case ('file_inp')
        grid_function%file_inp = val
      case ('file_out')
        grid_function%file_out = val
      case ('uscale')
        grid_function%modify = 1
        call str2flt(val, grid_function%up_scale, stat)
      case ('dscale')
        grid_function%modify = 1
        call str2flt(val, grid_function%down_scale, stat)
      case ('carpow')
        grid_function%modify = 1
        temp = val
        call trim_ch(temp, val,  '[')
        call trim_ch(temp, temp, ']')
        if (index(val, ',') > 0) then
          nr_words = word_count(temp, ',')
          if (nr_words .ne. 3) then
            write(*,*) 'Input error. More than three values are in the "lao_skip".'
            stop
          end if
          allocate(chvalues(nr_words))
          call split_string(temp, nr_words, chvalues, separator=',')
          do i = 1, nr_words
            call str2int(trim(chvalues(i)), np, stat)
            grid_function%cartesian_power(i) = np
          end do
          deallocate(chvalues)
        else
          write(*,*) 'Input error. "carpow" expects an array of floats.'
          stop
        end if
      case default
        print *, 'Input error. Unknown key.'
        stop
    end select

    !call modify_grid_function_due_to_method(grid_function)

    if (grid_function%id == 0) grid_function%id = id_grid_function

  end subroutine

! ------------------------------------------------------------------------------
  subroutine verify_grid_function_predefined(grid_function)
    class(grid_function_t), intent(inout):: grid_function
!
!   note: grid_function%name also defines a label on HDF5
!
    select case (grid_function%name)
      case ('ed',      &
            'gamma5',  &
            'esp',     &
            'espe',    &
            'espn',    &
            'esprho',  &
            'esperho', &
            'espnrho', & 
            'bedcos',  & 
            'bedsin')
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
            'divs') 
        grid_function%tensor_order = 0
        grid_function%max_ao_g_order = 1
        grid_function%max_ao_m_order = 0
      case ('s', &
            'jdiaBz', &
            'j', &
            'bdip',   &
            'ndip',   &
            'edip')
        grid_function%tensor_order = 1
        grid_function%max_ao_g_order = 0
        grid_function%max_ao_m_order = 0
      case ('ed_gradient',&
            'rots', &
            'rotj')
        grid_function%tensor_order = 1
        grid_function%max_ao_g_order = 1
        grid_function%max_ao_m_order = 0
      case ('gradj')
        grid_function%tensor_order = 2
        grid_function%max_ao_g_order = 1
        grid_function%max_ao_m_order = 0
      case ('ed_hessian')
        grid_function%tensor_order = 2
        grid_function%max_ao_g_order = 2
        grid_function%max_ao_m_order = 0
      case ('ed_laplacian',&
            'ed_sign_lambda2',&
            'kinlap',&
            'kinnr') 
        grid_function%tensor_order = 0
        grid_function%max_ao_g_order = 2
        grid_function%max_ao_m_order = 0
      case default
        write(*, *) 'Input error. Unknown grid function name.'
        stop
    end select

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
    call interface_mpi_bcast(rs_input%print_unit,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%cvalue,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%gauge_origin,3,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%nr_grids,1,irank_master,global_communicator)
    call interface_mpi_bcast(rs_input%nr_grid_functions,1,irank_master,global_communicator)
    if (.not.mpi_is_master()) then
      ! the allocation on master happens during DIRAC input reading (input/input_reader_sections.F90)
      allocate(rs_input%grids(rs_input%nr_grids))
      allocate(rs_input%grid_functions(rs_input%nr_grid_functions))
    end if

    ! legacy and variables from DIRAC's common blocks

    if (rs_input%print_level > 1) call logprint
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

! ------------------------------------------------------------------------------
! string mangling subroutines
! probably deserve a separate module and a better interface

   function lowercase(s)

!     copied from input/input_reader.F90
!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len(s))        :: lowercase
!     --------------------------------------------------------------------------
      integer                  :: off, i, ia
!     --------------------------------------------------------------------------

      lowercase = s

      off = iachar('a') - iachar('A')

      do i = 1, len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('A') .and. ia <= iachar('Z')) then
            lowercase(i:i) = achar(ia + off)
         end if
      enddo

   end function


! ------------------------------------------------------------------------------
  elemental subroutine str2int(str,int1,stat)
    implicit none
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int1
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int1
  end subroutine str2int

! ------------------------------------------------------------------------------
  elemental subroutine str2flt(str,flt,stat)
    implicit none
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

end module
