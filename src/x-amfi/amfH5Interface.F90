!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end
!
!
! interaface to create an hdf5 file directory
!
! written by Kai N. Spauszus (and stknecht) March 2022
!

module amfH5Interface


        use amfDataStorage
        use labeled_storage, only: lab_read, file_info_t
        use mh5

        implicit none

        public amf_write_hdf5
        public amf_read_hdf5

        interface amf_write_hdf5
            module procedure amf_write_1d
            module procedure amf_write_3d
        end interface amf_write_hdf5

        interface amf_read_hdf5
            module procedure amf_read_1d
            module procedure amf_read_3d
        end interface amf_read_hdf5


contains

!> public interface routines

subroutine amf_write_1d(hdfData,data_1d,nrow,ncol,ndim,dataName)

        type(hdfDataStorage),intent(inout)      :: hdfData
        real*8              ,intent(in)         :: data_1d
        integer             ,intent(in)         :: nrow, ncol, ndim
        character(len=8)    ,intent(in)         :: dataName

        real*8, target                          :: data_1d_local(1,1,1)

        data_1d_local(1,1,1) = data_1d

        nullify(hdfData%dataStorage%mat)
        hdfData%dataStorage%mat => data_1d_local

        hdfData%dataStorage%nrow=nrow; hdfData%dataStorage%ncol     = ncol
        hdfData%dataStorage%ndim=ndim; hdfData%dataStorage%dataName = dataName

        call amf_write_general(hdfData)

        if(associated(hdfData%dataStorage%mat)) nullify(hdfData%dataStorage%mat)

end subroutine amf_write_1d

subroutine amf_write_3d(hdfData,data_3d,nrow,ncol,ndim,dataName)

        type(hdfDataStorage),intent(inout)      :: hdfData
        real*8, target      ,intent(in)         :: data_3d(nrow,ncol,ndim)
        integer             ,intent(in)         :: nrow, ncol, ndim
        character(len=8)    ,intent(in)         :: dataName

        nullify(hdfData%dataStorage%mat)
        hdfData%dataStorage%mat => data_3d

        hdfData%dataStorage%nrow=nrow; hdfData%dataStorage%ncol     = ncol
        hdfData%dataStorage%ndim=ndim; hdfData%dataStorage%dataName = dataName

        call amf_write_general(hdfData)

        if(associated(hdfData%dataStorage%mat)) nullify(hdfData%dataStorage%mat)

end subroutine amf_write_3d

subroutine amf_read_1d(hdfData,data_1d,dataName)

        type(hdfDataStorage),intent(inout)      :: hdfData
        real*8              ,intent(out)        :: data_1d
        character(len=8)    ,intent(in)         :: dataName

        real*8, target                          :: data_1d_local(1,1,1)

        nullify(hdfData%dataStorage%mat)
        hdfData%dataStorage%mat => data_1d_local
        hdfData%dataStorage%dataName = dataName

        call amf_read_general(hdfData)

        data_1d = data_1d_local(1,1,1)

        if(associated(hdfData%dataStorage%mat)) nullify(hdfData%dataStorage%mat)

end subroutine amf_read_1d

subroutine amf_read_3d(hdfData,data_3d,dataName)

        type(hdfDataStorage),intent(inout)      :: hdfData
        real*8, target      ,intent(out)        :: data_3d(:,:,:)
        character(len=8)    ,intent(in)         :: dataName

        nullify(hdfData%dataStorage%mat)
        hdfData%dataStorage%mat => data_3d
        hdfData%dataStorage%dataName = dataName

        call amf_read_general(hdfData)

        if(associated(hdfData%dataStorage%mat)) nullify(hdfData%dataStorage%mat)

end subroutine amf_read_3d

!> private routines

subroutine amf_write_general(hdfData)

        type(hdfDataStorage),intent(inout)      :: hdfData
        character(len = 80)                     :: description
        character(len =  8)                     :: rows, cols, dims
        integer                                 :: i
        character(len=250)                      :: label

        !> open or append file
        hdfData%fileId = mh5_open_file_a(hdfData%fname)

        !> initialize
        allocate(hdfData%groupId(0:amfpathDim)); hdfData%groupId(0) = hdfData%fileId
        hdfData%dataStorage%dsetId = -1

        !> set full data label
        label = trim(create_data_label(hdfData))

        !> create directory according to path
        do i=1,amfpathDim
             hdfData%groupId(i) = mh5_create_group(hdfData%groupId(i-1), trim(hdfData%path(i)))
        end do

        !> check if data set exists (will be replaced)
        if (mh5_exists_dset(hdfData%fileId, label)) then
            hdfData%dataStorage%dsetId = mh5_open_dset(hdfData%fileId, label)
        else

            if(size(hdfData%dataStorage%mat) == 1) then
                hdfData%dataStorage%dsetId = mh5_create_dset_real(hdfData%groupId(amfpathDim)        ,&
                                                                  hdfData%dataStorage%dataName       ,&
                                                                  int8(1)                            ,&
                                                                  [size(hdfData%dataStorage%mat      ,&
                                                                  kind=8)                             &
                                                                  ])
                call mh5_init_attr(hdfData%dataStorage%dsetId, "DESCRIPTION", "scalar:f*8")
            else
                hdfData%dataStorage%dsetId  = mh5_create_dset_real(hdfData%groupId(amfpathDim)       ,&
                                                                   hdfData%dataStorage%dataName      ,&
                                                                   int8(3)                           ,&
                                                                  [int8(hdfData%dataStorage%nrow)    ,&
                                                                   int8(hdfData%dataStorage%ncol)    ,&
                                                                   int8(hdfData%dataStorage%ndim)     &
                                                                  ])
                write(rows, "(i0)") hdfData%dataStorage%nrow
                write(cols, "(i0)") hdfData%dataStorage%ncol
                write(dims, "(i0)") hdfData%dataStorage%ndim
                description = "arranged matrix of size["//trim(rows)//","//trim(cols)//","//trim(dims)//"]"
                call mh5_init_attr(hdfData%dataStorage%dsetId, "DESCRIPTION" , trim(description))
            end if
        end if
        call mh5_put_dset(hdfData%dataStorage%dsetId,hdfData%dataStorage%mat)

        call mh5_close_file(hdfData%fileId)

        deallocate(hdfData%groupId)

end subroutine amf_write_general

subroutine amf_read_general(hdfData)

        type(hdfDataStorage),intent(inout)      :: hdfData

        character(len=250)                      :: label
        type(file_info_t)                       :: local_info

        local_info%type = 2; local_info%status = 0; local_info%name = hdfData%fname

        label = trim(create_data_label(hdfData))

        call lab_read(local_info,label,rdata=hdfData%dataStorage%mat)

        deallocate(local_info%name)

end subroutine amf_read_general

function create_data_label(hdfData) result(label)
        type(hdfDataStorage),intent(in)         :: hdfData
        character(len=250)                      :: label
        if(amfpathDim /= 5) call quit('ERROR in create_data_label - PATH dimension /= 5!')
        label = trim(hdfData%path(1))//"/"//trim(hdfData%path(2))//"/"//trim(hdfData%path(3))//"/"//&
                trim(hdfData%path(4))//"/"//trim(hdfData%path(5))//"/"//trim(hdfData%dataStorage%dataName)
        write(6,*) ' full label is ',label; call flush(6)

end function create_data_label

end module amfH5Interface
