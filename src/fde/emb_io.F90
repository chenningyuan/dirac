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

module fde_io

   use fde_cfg
      
   implicit none

   public fde_open_file
   public fde_close_file
   public read_grid
   public read_grid_dim_only


   interface read_grid
      module procedure read_grid_onecol
      module procedure read_grid_manycol
      module procedure read_grid_proponly
   end interface read_grid

   real(kind=8) :: threshold = 1.0d-18

!  pertden info
   integer, public :: fde_pertden_ncol
   integer, public :: fde_pertden_nrow

   contains

     subroutine fde_open_file(name,unit)
        character(len=60), intent(in) :: name
        character(len=60)             :: errmsg
        logical                       :: file_found 
        integer, intent(in)           :: unit

        inquire(FILE=name,EXIST=file_found)
  
        if (file_found) then
            open(unit,                   &
                 FORM   = 'FORMATTED',   &
                 STATUS = 'UNKNOWN',     &
                 ACCESS = 'SEQUENTIAL',  &
                 POSITION = 'REWIND',    &
                 FILE   =  name)
         else
            write (errmsg,'(a60)') "FDE-related file: "//name(1:10)//" not found. Stopping now!"
            call fde_quit(errmsg)
         endif
     end subroutine fde_open_file

     subroutine fde_close_file(unit)
        integer, intent(in) :: unit
        close(unit)
     end subroutine fde_close_file



! description : this subroutine will read some operator (i.e. an embedding
!               potential) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the last the quadrature weight;
!                                 optionally, the value of the operator/property 
!                                 (the 5th column) is also read in
!

  SUBROUTINE READ_GRID_onecol(file,points,vc)
    REAL(kind=8),POINTER :: points(:,:)
    REAL(kind=8),pointer,optional :: vc(:)
    INTEGER                 :: file
    INTEGER                 :: i,npoints

    REWIND(file)
    READ(file,*) npoints
    
    if (associated(points)) nullify(points)
    allocate(points(4,npoints))
    
    if (present(vc)) then
       if (associated(vc)) nullify(vc)
       allocate(vc(npoints))
       DO i=1,npoints
          READ(file,*) points(:,i),vc(i)
       END DO

    else
       DO i=1,npoints
          READ(file,*) points(:,i)
       END DO
    end if
  END SUBROUTINE READ_GRID_onecol

! description : this subroutine will read one or more operators/properties (i.e. an embedding
!               potential, a density) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point, number of properties on file
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the fourth the quadrature weight.
!                                 these are followed by whatever properties we
!                                 have in file
!
  SUBROUTINE READ_GRID_manycol(file,points,properties)
    REAL(kind=8),POINTER :: points(:,:), ptmp_many(:,:)
    REAL(kind=8),pointer :: properties(:,:)
    INTEGER              :: file
    INTEGER              :: i,npoints, nprop, ierr
   
    REWIND(file)
    READ(file,*) npoints, nprop

    if (associated(points)) nullify(points)
    allocate(points(4,npoints))
    if (associated(properties)) nullify(properties)
    allocate(properties(nprop,npoints))

    DO i=1,npoints
       READ(file,*, IOSTAT=ierr) points(:,i),properties(:,i)
       if (ierr .ne. 0) write(*,*) 'error ierr in point', ierr, i
    END DO

  END SUBROUTINE READ_GRID_manycol
  
  SUBROUTINE READ_GRID_proponly(file,npoints,nprop,properties)
    REAL(kind=8),pointer :: properties(:,:)
    REAL(kind=8)         :: dummy(4)
    INTEGER              :: file
    INTEGER              :: i,npoints, nprop, ierr
   
    REWIND(file)
    if (associated(properties)) nullify(properties)
    allocate(properties(nprop,npoints))

    DO i=1,npoints
       READ(file,*, IOSTAT=ierr) dummy(:), properties(:,i)
       if (ierr .ne. 0) write(*,*) 'error ierr in point', ierr, i
    END DO

  END SUBROUTINE READ_GRID_proponly

   subroutine read_grid_dim_only(file,nr_points, nr_cols)
    integer    :: nr_points, nr_cols
    integer    :: file

    rewind(file)
    read(file,*) nr_points, nr_cols

  end subroutine


end module fde_io

