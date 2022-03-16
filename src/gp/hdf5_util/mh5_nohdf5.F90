module mh5

! Stub in case no hdf5 library is available

use, intrinsic :: iso_fortran_env, only: int32, int64, real64, output_unit
use, intrinsic :: iso_c_binding, only: c_char, c_null_char
use, intrinsic :: iso_c_binding, only: c_double
use, intrinsic :: iso_c_binding, only: c_long

implicit none

integer(kind=int64), parameter :: iwp = int64, MOLCAS_C_INT = c_long
integer(kind=iwp), parameter :: wp = real64, MOLCAS_C_REAL = c_double

private

public :: mh5_create_file, mh5_open_file_a mh5_open_file_rw, mh5_open_file_r, mh5_close_file, &
          mh5_is_hdf5, mh5_open_group, mh5_exists_group, mh5_create_group, mh5_close_group, &
          mh5_exists_attr, mh5_open_attr, mh5_close_attr, mh5_exists_dset, mh5_open_dset, mh5_close_dset, mh5_create_attr_int, &
          mh5_create_attr_real, mh5_create_attr_str, mh5_put_attr, mh5_get_attr, mh5_init_attr, mh5_fetch_attr, &
          mh5_create_dset_int, mh5_create_dset_real, mh5_create_dset_str, mh5_put_dset, mh5_get_dset, mh5_init_dset, &
          mh5_fetch_dset, mh5_resize_dset, mh5_get_dset_dims

!======================
! Overloaded interfaces
!======================

! attr

interface mh5_create_attr_int
  module procedure :: mh5_create_attr_scalar_int, &
                      mh5_create_attr_array_int
end interface mh5_create_attr_int

interface mh5_create_attr_real
  module procedure :: mh5_create_attr_scalar_real, &
                      mh5_create_attr_array_real
end interface mh5_create_attr_real

interface mh5_create_attr_str
  module procedure :: mh5_create_attr_scalar_str, &
                      mh5_create_attr_array_str
end interface mh5_create_attr_str

interface mh5_put_attr
  module procedure :: mh5_put_attr_scalar_int, &
                      mh5_put_attr_scalar_real, &
                      mh5_put_attr_scalar_str, &
                      mh5_put_attr_array_int, &
                      mh5_put_attr_array_real, &
                      mh5_put_attr_array_str
end interface mh5_put_attr

interface mh5_get_attr
  module procedure :: mh5_get_attr_scalar_int, &
                      mh5_get_attr_scalar_real, &
                      mh5_get_attr_scalar_str, &
                      mh5_get_attr_array_int, &
                      mh5_get_attr_array_real, &
                      mh5_get_attr_array_str
end interface mh5_get_attr

! init = create + put + close
interface mh5_init_attr
  module procedure :: mh5_init_attr_scalar_int, &
                      mh5_init_attr_scalar_real, &
                      mh5_init_attr_scalar_str, &
                      mh5_init_attr_array_int, &
                      mh5_init_attr_array_real, &
                      mh5_init_attr_array_str
end interface mh5_init_attr

! fetch = open + get + close
interface mh5_fetch_attr
  module procedure :: mh5_fetch_attr_scalar_int, &
                      mh5_fetch_attr_scalar_real, &
                      mh5_fetch_attr_scalar_str, &
                      mh5_fetch_attr_array_int, &
                      mh5_fetch_attr_array_real, &
                      mh5_fetch_attr_array_str
end interface mh5_fetch_attr

! dset

interface mh5_create_dset_int
  module procedure :: mh5_create_dset_scalar_int, &
                      mh5_create_dset_array_int
end interface mh5_create_dset_int

interface mh5_create_dset_real
  module procedure :: mh5_create_dset_scalar_real, &
                      mh5_create_dset_array_real
end interface mh5_create_dset_real

interface mh5_create_dset_str
  module procedure :: mh5_create_dset_scalar_str, &
                      mh5_create_dset_array_str
end interface mh5_create_dset_str

interface mh5_put_dset
  module procedure :: mh5_put_dset_scalar_int, &
                      mh5_put_dset_scalar_real, &
                      mh5_put_dset_scalar_str, &
                      mh5_put_dset_array_int, &
                      mh5_put_dset_array_real, &
                      mh5_put_dset_array_str
end interface mh5_put_dset

interface mh5_get_dset
  module procedure :: mh5_get_dset_scalar_int, &
                      mh5_get_dset_scalar_real, &
                      mh5_get_dset_scalar_str, &
                      mh5_get_dset_array_int, &
                      mh5_get_dset_array_real, &
                      mh5_get_dset_array_str
end interface mh5_get_dset

! init = create + put + close
interface mh5_init_dset
  module procedure :: mh5_init_dset_scalar_int, &
                      mh5_init_dset_scalar_real, &
                      mh5_init_dset_scalar_str, &
                      mh5_init_dset_array_int, &
                      mh5_init_dset_array_real, &
                      mh5_init_dset_array_str
end interface mh5_init_dset

! fetch = open + get + close
interface mh5_fetch_dset
  module procedure :: mh5_fetch_dset_scalar_int, &
                      mh5_fetch_dset_scalar_real, &
                      mh5_fetch_dset_scalar_str, &
                      mh5_fetch_dset_array_int, &
                      mh5_fetch_dset_array_real, &
                      mh5_fetch_dset_array_str
end interface mh5_fetch_dset

interface mh5_resize_dset
  module procedure :: mh5_extend_dset_array
end interface mh5_resize_dset

interface mh5_get_dset_dims
  module procedure :: mh5_get_dset_array_dims
end interface mh5_get_dset_dims

contains

!==========================
! Specific Fortran wrappers
!==========================

#define MH5_MAX_LBL_LEN 256

! file

function mh5_create_file(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
end function mh5_create_file

function mh5_open_file_a(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
end function mh5_open_file_a

function mh5_open_file_rw(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
end function mh5_open_file_rw

function mh5_open_file_r(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
end function mh5_open_file_r

subroutine mh5_close_file(lu)
  integer(kind=iwp), intent(in) :: lu
end subroutine mh5_close_file

function mh5_is_hdf5(filename) result(ishdf5)
  logical(kind=iwp) :: ishdf5
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: rc, lrealname
  character(len=4096) :: realname, c_name
end function mh5_is_hdf5

! group

function mh5_open_group(lu,groupname) result(groupid)
  integer(kind=iwp) :: groupid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: groupname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_open_group

function mh5_exists_group(lu,groupname) result(groupid)
  integer(kind=iwp) :: groupid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: groupname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_exists_group

function mh5_create_group(lu,groupname) result(groupid)
  integer(kind=iwp) :: groupid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: groupname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_group

subroutine mh5_close_group(id)
  integer(kind=iwp), intent(in) :: id
end subroutine mh5_close_group

! attr

function mh5_exists_attr(id,attrname) result(exists)
  logical(kind=iwp) :: exists
  integer(kind=iwp), intent(in) :: id
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  integer(kind=iwp) :: rc
end function mh5_exists_attr

function mh5_open_attr(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_open_attr

subroutine mh5_close_attr(attrid)
  integer(kind=iwp), intent(in) :: attrid
end subroutine mh5_close_attr

! dset

function mh5_exists_dset(id,dsetname) result(exists)
  logical(kind=iwp) :: exists
  integer(kind=iwp), intent(in) :: id
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  integer(kind=iwp) :: rc
end function mh5_exists_dset

function mh5_open_dset(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_open_dset

subroutine mh5_close_dset(dsetid)
  integer(kind=iwp), intent(in) :: dsetid
end subroutine mh5_close_dset

! attr types

function mh5_create_attr_scalar_int(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_scalar_int

function mh5_create_attr_scalar_real(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_scalar_real

function mh5_create_attr_scalar_str(lu,attrname,length) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, length
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_scalar_str

subroutine mh5_put_attr_scalar_int(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(in) :: val
end subroutine mh5_put_attr_scalar_int

subroutine mh5_put_attr_scalar_real(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(in) :: val
end subroutine mh5_put_attr_scalar_real

subroutine mh5_put_attr_scalar_str(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  character(len=*), intent(in) :: val
end subroutine mh5_put_attr_scalar_str

subroutine mh5_get_attr_scalar_int(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(out) :: val
end subroutine mh5_get_attr_scalar_int

subroutine mh5_get_attr_scalar_real(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(out) :: val
end subroutine mh5_get_attr_scalar_real

subroutine mh5_get_attr_scalar_str(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  character(len=*), intent(out) :: val
end subroutine mh5_get_attr_scalar_str

function mh5_create_attr_array_int(lu,attrname,rank,dims) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_array_int

function mh5_create_attr_array_real(lu,attrname,rank,dims) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_array_real

function mh5_create_attr_array_str(lu,attrname,rank,dims,length) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_attr_array_str

subroutine mh5_put_attr_array_int(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(in) :: buffer(*)
end subroutine mh5_put_attr_array_int

subroutine mh5_put_attr_array_real(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(in) :: buffer(*)
end subroutine mh5_put_attr_array_real

subroutine mh5_put_attr_array_str(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  character, intent(in) :: buffer(*)
end subroutine mh5_put_attr_array_str

subroutine mh5_get_attr_array_int(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(inout) :: buffer(*)
end subroutine mh5_get_attr_array_int

subroutine mh5_get_attr_array_real(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(inout) :: buffer(*)
end subroutine mh5_get_attr_array_real

subroutine mh5_get_attr_array_str(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  character, intent(inout) :: buffer(*)
end subroutine mh5_get_attr_array_str

subroutine mh5_init_attr_scalar_int(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(in) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_scalar_int

subroutine mh5_init_attr_scalar_real(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(in) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_scalar_real

subroutine mh5_init_attr_scalar_str(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=*), intent(in) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_scalar_str

subroutine mh5_init_attr_array_int(lu,attrname,rank,dims,buffer)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_array_int

subroutine mh5_init_attr_array_real(lu,attrname,rank,dims,buffer)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_array_real

subroutine mh5_init_attr_array_str(lu,attrname,rank,dims,buffer,length)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: attrname
  character, intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_init_attr_array_str

subroutine mh5_fetch_attr_scalar_int(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(out) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_scalar_int

subroutine mh5_fetch_attr_scalar_real(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(out) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_scalar_real

subroutine mh5_fetch_attr_scalar_str(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=*), intent(out) :: val
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_scalar_str

subroutine mh5_fetch_attr_array_int(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(inout) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_array_int

subroutine mh5_fetch_attr_array_real(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(inout) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_array_real

subroutine mh5_fetch_attr_array_str(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character, intent(inout) :: buffer(*)
  integer(kind=iwp) :: attrid
end subroutine mh5_fetch_attr_array_str

! dset types

function mh5_create_dset_scalar_int(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_dset_scalar_int

function mh5_create_dset_scalar_real(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_dset_scalar_real

function mh5_create_dset_scalar_str(lu,dsetname,length) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, length
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
end function mh5_create_dset_scalar_str

subroutine mh5_put_dset_scalar_int(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: val
end subroutine mh5_put_dset_scalar_int

subroutine mh5_put_dset_scalar_real(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: val
end subroutine mh5_put_dset_scalar_real

subroutine mh5_put_dset_scalar_str(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  character(len=*), intent(in) :: val
end subroutine mh5_put_dset_scalar_str

subroutine mh5_get_dset_scalar_int(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(out) :: val
end subroutine mh5_get_dset_scalar_int

subroutine mh5_get_dset_scalar_real(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(out) :: val
end subroutine mh5_get_dset_scalar_real

subroutine mh5_get_dset_scalar_str(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  character(len=*), intent(out) :: val
end subroutine mh5_get_dset_scalar_str

function mh5_create_dset_array_int(lu,dsetname,rank,dims,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
end function mh5_create_dset_array_int

function mh5_create_dset_array_real(lu,dsetname,rank,dims,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
end function mh5_create_dset_array_real

function mh5_create_dset_array_str(lu,dsetname,rank,dims,length,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
end function mh5_create_dset_array_str

subroutine mh5_put_dset_array_int(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_put_dset_array_int

subroutine mh5_put_dset_array_real(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_put_dset_array_real

subroutine mh5_put_dset_array_str(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  character, intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_put_dset_array_str

subroutine mh5_get_dset_array_int(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(inout) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_get_dset_array_int

subroutine mh5_get_dset_array_real(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(inout) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_get_dset_array_real

subroutine mh5_get_dset_array_str(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  character, intent(inout) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
end subroutine mh5_get_dset_array_str

subroutine mh5_extend_dset_array(dsetid,dims)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: dims(*)
end subroutine mh5_extend_dset_array

subroutine mh5_get_dset_array_dims(dsetid,dims)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(inout) :: dims(*)
end subroutine mh5_get_dset_array_dims

subroutine mh5_init_dset_scalar_int(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(in) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_init_dset_scalar_int

subroutine mh5_init_dset_scalar_real(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(in) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_init_dset_scalar_real

subroutine mh5_init_dset_scalar_str(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=*), intent(in) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_init_dset_scalar_str

subroutine mh5_init_dset_array_int(lu,dsetname,rank,dims,buffer,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
end subroutine mh5_init_dset_array_int

subroutine mh5_init_dset_array_real(lu,dsetname,rank,dims,buffer,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
end subroutine mh5_init_dset_array_real

subroutine mh5_init_dset_array_str(lu,dsetname,rank,dims,buffer,length,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: dsetname
  character, intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
end subroutine mh5_init_dset_array_str

subroutine mh5_fetch_dset_scalar_int(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(out) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_fetch_dset_scalar_int

subroutine mh5_fetch_dset_scalar_real(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(out) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_fetch_dset_scalar_real

subroutine mh5_fetch_dset_scalar_str(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=*), intent(out) :: val
  integer(kind=iwp) :: dsetid
end subroutine mh5_fetch_dset_scalar_str

subroutine mh5_fetch_dset_array_int(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(inout) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
end subroutine mh5_fetch_dset_array_int

subroutine mh5_fetch_dset_array_real(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(inout) :: buffer(*)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
end subroutine mh5_fetch_dset_array_real

subroutine mh5_fetch_dset_array_str(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character, intent(inout) :: buffer(*)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_str

end module mh5
