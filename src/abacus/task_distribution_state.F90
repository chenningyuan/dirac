module task_distribution_state
  implicit none

  integer, allocatable, save, public :: task_npos(:)
  integer, save :: task_npos_size = 0

contains

  integer function task_array_dimension(is_parallel) result(tsk_arr_dim)
#include "implicit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
    logical, intent(in) :: is_parallel

    tsk_arr_dim = 0
    if (is_parallel) then
      tsk_arr_dim = tsk_arr_dim + nlrgbl * (nlrgbl + 1) / 2
      tsk_arr_dim = tsk_arr_dim + nlrgbl * (nlrgbl + 1) / 2
      tsk_arr_dim = tsk_arr_dim + nsmlbl * (nsmlbl + 1) / 2
      tsk_arr_dim = tsk_arr_dim + nsmlbl * nlrgbl
    end if
  end function task_array_dimension

  subroutine ensure_task_order_storage(is_parallel)
    logical, intent(in) :: is_parallel
    integer :: required_size

    required_size = task_array_dimension(is_parallel)
    if (required_size <= 0) return

    if (.not. allocated(task_npos)) then
      allocate(task_npos(required_size))
      task_npos_size = required_size
    else if (task_npos_size /= required_size) then
      deallocate(task_npos)
      allocate(task_npos(required_size))
      task_npos_size = required_size
    end if
  end subroutine ensure_task_order_storage

  integer function task_class_offset(i2typ) result(offset)
#include "implicit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
    integer, intent(in) :: i2typ

    select case (i2typ)
    case (1)
      offset = 1
    case (2)
      offset = 1 + nlrgbl * (nlrgbl + 1) / 2
    case (3)
      offset = 1 + 2 * (nlrgbl * (nlrgbl + 1) / 2)
    case (4)
      offset = 1 + 2 * (nlrgbl * (nlrgbl + 1) / 2) + nsmlbl * (nsmlbl + 1) / 2
    case default
      offset = 1
    end select
  end function task_class_offset

end module task_distribution_state
