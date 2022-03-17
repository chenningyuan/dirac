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
! global data storage type that is passed to amfH5Interface
!
! written by Kai N. Spauszus (and stknecht) March 2022
!
module amfDataStorage
        implicit none

        type, public                    :: dataStorage
                character(len= 8)       :: dataName
                integer                 :: nrow = -1
                integer                 :: ncol = -1
                integer                 :: ndim = -1
                real(8), pointer        :: mat(:,:,:)
        end type dataStorage


        type, public                            :: hdfDataStorage
                character(len=10)               :: fname = "amfPCEC.h5"
                character(len = 40),allocatable :: path(:)
               ! integer                         :: pathDim = 5
                integer(kind=8)                 :: fileId
                integer(kind=8),allocatable     :: groupId(:)
                integer(kind=8),allocatable     :: dsetId(:)
                type(dataStorage),allocatable   :: dataStorage(:)

        end type

        integer, parameter                      :: amfpathDim = 5

end module amfDataStorage
