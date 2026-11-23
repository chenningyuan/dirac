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
module fragment_x2c_cfg

! stefan: this module contains all required internal basis set, offset and 
!         pointer informations for the fragment-/or atomic-X2C approach.

  !use x2c_fio

  implicit none

  public fragment_x2c
  public fragment_x2c_init
  public fragment_x2c_free

! type
  type fragment_x2c

!   double precision block
    real(8), allocatable ::      &
    pctmat(:)

!   integer block
    integer              ::      &
    charge,                      &
    funit
    integer, allocatable ::      &
    naosh_all(:),                &
    naosh_L(:)

!   logical block
    logical              ::                &
    fragment_approach_enabled,             &
    fragment_approach_ismolecule = .false.,&
    initialized
    
  end type

! fragment_x2c type
  type(fragment_x2c), public, save :: fragment_x2c_info
! ----------------------------------------------------------------------------

! parameters

  integer, parameter, public :: max_nr_symm_ind_centers = 2000

contains 

!*******************************************************************************
  subroutine fragment_x2c_init(A,                     & 
                               nfsym,                 &
                               nz,                    &
                               charge_cent,           &
                               naosh_all_cent,        &
                               naosh_L_cent)
!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    use labeled_storage
    type(fragment_x2c)           :: A
    integer, intent(in)          :: nfsym
    integer, intent(in)          :: nz
    integer, intent(in)          :: charge_cent
    integer, intent(in)          :: naosh_all_cent(nfsym)
    integer, intent(in)          :: naosh_L_cent(nfsym)
!   ----------------------------------------------------------------------------
    integer                      :: pct_frag_dimension
    integer                      :: i
    character(len=3)             :: extension
    character(len=10)            :: filename
    character(len=12)            :: flabel
    type(file_info_t)            :: fragfile
!   ----------------------------------------------------------------------------

!   reset old type information
    call fragment_x2c_free(A)

    A%initialized              = .true.

    allocate(A%naosh_all(nfsym))
    allocate(A%naosh_L(nfsym))

    A%naosh_all(1:nfsym) = -1
    A%naosh_L(1:nfsym)   = -1
    A%charge             = charge_cent
    A%funit              = 89

    pct_frag_dimension   = 0
    do i = 1, nfsym
      A%naosh_all(i)     = naosh_all_cent(i)
      A%naosh_L(i)       = naosh_L_cent(i)
      pct_frag_dimension = pct_frag_dimension + naosh_all_cent(i) * naosh_L_cent(i) * nz
    end do

    allocate(A%pctmat(pct_frag_dimension))

    A%pctmat         =  0

!   open file with fragment U matrix
    extension        = '000' 
    write(extension,'(i3)') A%charge

    if(extension(1:1) == ' ') write(extension(1:1),'(i1)') 0
    if(extension(2:2) == ' ') write(extension(2:2),'(i1)') 0
    write(filename,'(a7,a3)') 'X2CMAT_',extension

    fragfile%type = 2
    fragfile%name = filename//'.h5'
    fragfile%status = 0
    !open(A%funit,file=filename,status='old',form='unformatted',    &
    !     access='sequential',action='readwrite',position='rewind')

!   read atomic pctmat from file
    write(flabel,'(a7,i4,i1)') 'pctmtAO',1,1
    call lab_read(fragfile,'/result/hamiltonian/x2c/ao_matrices/'//flabel,rdata=A%pctmat)

    !call x2c_read(flabel,                 &
    !              A%pctmat,               &
    !              A%naosh_all(1) *        &
    !              A%naosh_L(1)   *        &
    !              nz,                     &
    !              A%funit)

  end subroutine fragment_x2c_init
!*******************************************************************************

  subroutine fragment_x2c_free(A)

!   ----------------------------------------------------------------------------
    type(fragment_x2c) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%initialized) return

    close(A%funit, status='keep')

    A%initialized = .false.
    A%charge      = -1
    A%funit       = -1
    deallocate(A%naosh_all)
    deallocate(A%naosh_L)
    deallocate(A%pctmat)

  end subroutine fragment_x2c_free
!*******************************************************************************

end module fragment_x2c_cfg
