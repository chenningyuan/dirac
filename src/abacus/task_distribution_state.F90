module task_distribution_state
  implicit none

  private
  public :: ensure_task_order_storage, task_order, task_class_offset
  public :: reset_task_order, task_order_is_valid, mark_task_order_valid

  integer, allocatable, save :: task_order(:)
  integer, save :: task_order_size = 0
  logical, save :: task_order_valid(4) = .false.

contains

  integer function required_task_order_size(is_parallel) result(required_size)
      use aovec_mod, only: MXAOVC
      use maxorb_mod, only: MXCORB, MXSHEL
      use blocks_mod, only: NLRGBL, NSMLBL
#include "implicit.h"
!#include "dcbgen.h"
!#include "dcbdhf.h"
!#include "dcbbas.h"
    logical, intent(in) :: is_parallel

    required_size = 0
    if (is_parallel) then
      required_size = required_size + nlrgbl * (nlrgbl + 1) / 2
      required_size = required_size + nlrgbl * (nlrgbl + 1) / 2
      required_size = required_size + nsmlbl * (nsmlbl + 1) / 2
      required_size = required_size + nsmlbl * nlrgbl
    end if
  end function required_task_order_size

  subroutine ensure_task_order_storage(is_parallel)
    logical, intent(in) :: is_parallel
    integer :: required_size

    required_size = required_task_order_size(is_parallel)
    if (required_size <= 0) return

    if (.not. allocated(task_order)) then
      allocate(task_order(required_size))
      task_order_size = required_size
    else if (task_order_size /= required_size) then
      deallocate(task_order)
      allocate(task_order(required_size))
      task_order_size = required_size
    end if
  end subroutine ensure_task_order_storage

  subroutine reset_task_order()
    task_order_valid = .false.
  end subroutine reset_task_order

  logical function task_order_is_valid(i2typ) result(is_valid)
    integer, intent(in) :: i2typ

    is_valid = .false.
    if (i2typ >= 1 .and. i2typ <= size(task_order_valid)) then
      is_valid = task_order_valid(i2typ)
    end if
  end function task_order_is_valid

  subroutine mark_task_order_valid(i2typ)
    integer, intent(in) :: i2typ

    if (i2typ >= 1 .and. i2typ <= size(task_order_valid)) then
      task_order_valid(i2typ) = .true.
    end if
  end subroutine mark_task_order_valid

  integer function task_class_offset(i2typ) result(offset)
      use aovec_mod, only: MXAOVC
      use maxorb_mod, only: MXCORB, MXSHEL
      use blocks_mod, only: NLRGBL, NSMLBL
#include "implicit.h"
!#include "dcbgen.h"
!#include "dcbdhf.h"
!#include "dcbbas.h"
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
