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
        use mh5

        implicit none

        public write_hdf5
        public set_data_hdf5
        public release_data_hdf5

contains

subroutine set_data_hdf5(hdfData,xdata,nrow,ncol,ndim,dataName)
        type(hdfDataStorage),intent(inout)      :: hdfData
        real*8 ,             intent(in), target :: xdata(:,:,:)
        integer,             intent(in)         :: nrow, ncol, ndim
        character(len=8),    intent(in)         :: dataName

        allocate(hdfData%dataStorage(1))
        nullify(hdfData%dataStorage(1)%mat)
        hdfData%dataStorage(1)%mat => xdata

        hdfData%dataStorage(1)%nrow=nrow; hdfData%dataStorage(1)%ncol     = ncol
        hdfData%dataStorage(1)%ndim=ndim; hdfData%dataStorage(1)%dataName = dataName

end subroutine set_data_hdf5

subroutine release_data_hdf5(hdfData)
        type(hdfDataStorage),intent(inout)      :: hdfData

        if(associated(hdfData%dataStorage(1)%mat)) nullify(hdfData%dataStorage(1)%mat)
        deallocate(hdfData%dataStorage)

end subroutine release_data_hdf5

subroutine write_hdf5(hdfData)

        type(hdfDataStorage),intent(inout)      :: hdfData
        character(len = 80)                     :: description
        character(len =  8)                     :: rows, cols, dims
        integer                                 :: i

        hdfData%fileId = mh5_open_file_a(hdfData%fname)   ! append file
        allocate(hdfData%groupId(0:amfpathDim)); hdfData%groupId(0) = hdfData%fileId

        do i=1,amfpathDim   ! create directory according to path
                hdfData%groupId(i) = mh5_create_group(hdfData%groupId(i-1), trim(hdfData%path(i)))
        end do

        allocate(hdfData%dsetId(size(hdfData%dataStorage))); hdfData%dsetId = -1

        do i=1, size(hdfData%dataStorage)   !write data according to size of dset and add size as attribute
                if(size(hdfData%dataStorage(i)%mat) == 1) then
                        hdfData%dsetId(i) = mh5_create_dset_real(hdfData%groupId(amfpathDim)              ,&
                                                                 hdfData%dataStorage(i)%dataName          ,&
                                                                 int8(1)                                  ,&
                                                                 [size(hdfData%dataStorage(i)%mat         ,&
                                                                  kind=8)                                  &
                                                                 ])
                        call mh5_init_attr(hdfData%dsetId(i), "DESCRIPTION", "scalar:f*8") 
                else
                        hdfData%dsetId(i) = mh5_create_dset_real(hdfData%groupId(amfpathDim)              ,&
                                                                 hdfData%dataStorage(i)%dataName          ,&
                                                                 int8(3)                                  ,&
                                                                 [int8(hdfData%dataStorage(i)%nrow)       ,&
                                                                  int8(hdfData%dataStorage(i)%ncol)       ,&
                                                                  int8(hdfData%dataStorage(i)%ndim)]       &
                                                                 )
                        write(rows, "(i0)") hdfData%dataStorage(i)%nrow
                        write(cols, "(i0)") hdfData%dataStorage(i)%ncol
                        write(dims, "(i0)") hdfData%dataStorage(i)%ndim
                        description = "arranged matrix of size["//trim(rows)//","//trim(cols)//","//trim(dims)//"]"
                        call mh5_init_attr(hdfData%dsetId(i), "DESCRIPTION" , trim(description))
                end if
                call mh5_put_dset(hdfData%dsetId(i),hdfData%dataStorage(i)%mat)
        end do

        call mh5_close_file(hdfData%fileId)

        deallocate(hdfData%groupId,hdfData%dsetId)

end subroutine write_hdf5

end module amfH5Interface
