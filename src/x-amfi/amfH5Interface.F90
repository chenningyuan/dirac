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

contains

subroutine write_hdf5(hdfData)

        type(hdfDataStorage),intent(inout)      :: hdfData
        character(len = 50)                     :: description
        character(len =  3)                     :: rows, cols, dims
        integer                                 :: i 

        print *, 'hdfData%fname = ',hdfData%fname

        hdfData%fileId = mh5_open_file_a(hdfData%fname)   ! append file
        allocate(hdfData%groupId(0:amfpathDim)); hdfData%groupId(0) = hdfData%fileId
        
        do i=1,amfpathDim   ! create directory according to path                
                hdfData%groupId(i) = mh5_create_group(hdfData%groupId(i-1), trim(hdfData%path(i)))
        end do
                
        print *, 'hdfData%sizedatastorage = ',size(hdfData%dataStorage)

        do i=1, size(hdfData%dataStorage)   !write data according to size of dset and add size as attribute
                if(size(hdfData%dataStorage(i)%mat) == 1) then
                        hdfData%dsetId(i) = mh5_create_dset_real(hdfData%groupId(amfpathDim)              ,&
                                                                 hdfData%dataStorage(i)%dataName          ,&
                                                                 int8(1)                                  ,&
                                                                 [size(hdfData%dataStorage(i)%mat         ,&
                                                                  kind=8)                                  &
                                                                 ])
                        call mh5_init_attr(hdfData%dsetId(i), "DESCRIPTION", "scalar:f*8") 
                else!if(size(hdfData%dataStorage(i)%mat)==3) then
                        print *, 'hdfData%sizedatastorageMAT = ',size(hdfData%dataStorage(i)%mat)
                        print *, 'hdfData%groupID            = ',hdfData%groupId(0:amfpathDim), amfpathDim
                        print *, 'hdfData%dataNAME           = ',hdfData%dataStorage(i)%dataName
                        print *, 'hdfData%datanrow           = ',hdfData%dataStorage(i)%nrow
                        print *, 'hdfData%datancol           = ',hdfData%dataStorage(i)%ncol
                        print *, 'hdfData%datandim           = ',hdfData%dataStorage(i)%ndim
                        print *, 'hdfData%datamat(1:10)      = ',hdfData%dataStorage(i)%mat(1:10)
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
                        description = "arranged matrix of size["//rows//","//cols//","//dims//"]"
                        call mh5_init_attr(hdfData%dsetId(i), "DESCRIPTION" , description)
                end if
                call mh5_put_dset(hdfData%dsetId(i),hdfData%dataStorage(i)%dataName)
        end do

        call mh5_close_file(hdfData%fileId)

end subroutine write_hdf5

end module amfH5Interface
