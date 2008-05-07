! 
!...   Copyright (c) 2004 by the authors of Dirac (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   Dirac, a relativistic ab initio electronic structure program,
!...   release DIRAC04.0 (2004),
!...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
!...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
!...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
!...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
!...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
!...   (http://dirac.chem.sdu.dk).
!...
!...   For a suitable BibTEX entry, see:
!...      http://dirac.chem.sdu.dk/doc/reference.shtml
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dirac,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For questions concerning this copyright write to:
!...      dirac-admin@dirac.chem.sdu.dk
!...
!...   For information on how to get a licence see:
!...      http://dirac.chem.sdu.dk
! 
MODULE GRID_FILE
  IMPLICIT NONE
  PRIVATE
  PUBLIC READ_GRID,WRITE_GRID, WRITE_XML
  PUBLIC GRID_DATA
  
  TYPE grid_data
     REAL(KIND=8) :: x(3)
     REAL(KIND=8) :: weight
     REAL(KIND=8) :: vc=0.0
     REAL(KIND=8) :: rho=0.0
     REAL(KIND=8) :: grad_rho(3)=0.0
     REAL(KIND=8) :: hes_rho(6)=0.0
  END TYPE grid_data
  
  
  !     XML output settings.
  !     Table that gives the tags that belong to existing variable names.

  CHARACTER(LEN=80),parameter :: tags(3,7)=reshape( (/ &
       &     'grid',      'grid',   '', &
       &     'grad_rho',  'dataset','name="gradient"', &
       &     'gridpoints','dataset','name="gridpoints"', &
       &     'hes_rho',   'dataset','name="hessian"', &
       &     'rho',       'dataset','name="density"', &
       &     'vc',        'dataset','name="vc"', &
       &     'weight',    'dataset','name="weight"'/),(/3,7/))
  INTEGER :: indent=0

CONTAINS
  
  
  SUBROUTINE WRITE_GRID(filename,points)
    USE FILE, ONLY: OPEN_FILE        
#include "priunit.h"
    CHARACTER(LEN=*)           :: filename
    TYPE(grid_data),POINTER    :: points(:)
    INTEGER                    :: i, fileid
    
    fileid=OPEN_FILE(filename)
    if (fileid<0) then
       WRITE(LUPRI,'(a)') "Error opening output file. file=",filename,"error", fileid
       RETURN
    end if
    
    WRITE(fileid,'(i8, a160)') size(points),'x y z weight rho V_Cintegrals+V_Cnuclei gradrho_x &
         &     gradrho_y gradrho_z hesrho_xx hesrho_xy hesrho_xz hesrho_yy  hesrho_yz hesrho_zz'
    DO i=1,size(points)
       WRITE(fileid,'(100e16.8)') points(i)%x,points(i)%weight,points(i)%rho,points(i)%vc, &
            &          points(i)%grad_rho, points(i)%hes_rho 
    END DO
    CLOSE(fileid)
    
  END SUBROUTINE WRITE_GRID
  
  CHARACTER(LEN=80) FUNCTION TAG(key,where,str)
    CHARACTER(LEN=*) :: key
    CHARACTER(LEN=1),OPTIONAL :: where
    CHARACTER(LEN=*),optional :: str
    INTEGER :: i
    TAG=''
    DO i=1,size(tags,2)
       if (trim(key)==trim(tags(1,i))) then 
          select case (where)
          case ('B'); if (present(str)) then; TAG='<'//trim(tags(2,i))//' '//trim(tags(3,i))//' '//trim(str)//'>'
                      else; TAG='<'//trim(tags(2,i))//' '//trim(tags(3,i))//'>'; endif
          case ('E'); TAG='</'//trim(tags(2,i))//'>'
          case default; TAG=tags(2,i)
          end select
       end if
    END DO
  END FUNCTION TAG
  
  CHARACTER(LEN=80) FUNCTION tab(var)
    USE conversion

    CHARACTER(LEN=*),OPTIONAL :: var
    if (present(var)) then
       tab='('//var//')'
       if (indent>0) tab='(t'//trim(inttochar(indent))//', '//var//')'                 
    else
       tab='(a)'
       if (indent>0) tab='(t'//trim(inttochar(indent))//', a)'                 
    end if
  END FUNCTION tab
  
  SUBROUTINE XML_SECTION(fileid,points,key)
    TYPE(grid_data) :: points(:)
    CHARACTER(LEN=*) :: key       
    character(len=80) :: str
    INTEGER :: fileid,i
#include "priunit.h"
    
   
    select case(key)
    case('gridpoints'); str='width="3"'
    case('grad_rho'); str='width="3"'
    case('hes_rho'); str='width="6"'
    case default; str=''
    end select

    WRITE(LUPRI,*),tab(),TAG(key,'B',str)
    WRITE(fileid,tab()) TAG(key,'B',str)
    indent=indent+3
    
    SELECT CASE(key)
    CASE('gridpoints')
       DO i=1,size(points)
          WRITE(fileid,tab('4e24.16')) points(i)%x
       END DO
    CASE('weight')
       DO i=1,size(points)
          WRITE(fileid,tab('4e24.16')) points(i)%weight
       END DO
    CASE('vc')
       DO i=1,size(points)
          WRITE(fileid,tab('1e24.16')) points(i)%vc
       END DO
    CASE('rho')
       DO i=1,size(points)
          WRITE(fileid,tab('1e24.16')) points(i)%rho
       END DO
    CASE('grad_rho')
       DO i=1,size(points)
          WRITE(fileid,tab('3e24.16')) points(i)%grad_rho
       END DO
    CASE('hes_rho')
       DO i=1,size(points)
          WRITE(fileid,tab('9e24.16')) points(i)%hes_rho
       END DO
    END SELECT
    
    indent=indent-3
    WRITE(fileid,tab()) TAG(key,'E')
    
  END SUBROUTINE XML_SECTION
  
  SUBROUTINE WRITE_XML(filename,points, grid_name)
    USE FILE, ONLY: OPEN_FILE 
    USE conversion       
#include "priunit.h"
    CHARACTER(LEN=*)           :: filename
    CHARACTER(LEN=*),optional  :: grid_name
    CHARACTER(LEN=80)          :: str=''
    TYPE(grid_data),POINTER    :: points(:)
    INTEGER                    :: i, fileid
    
    indent=0
    !     First write the gridpoints.
    
    fileid=OPEN_FILE(filename)
    if (fileid<0) then
       WRITE(LUPRI,'(a)') "Error opening XML file. file=",filename,"error", fileid
       RETURN
    end if       
    
!   Setup the information of the grid
    if (present(grid_name)) str=trim(str)//' name="'//trim(grid_name)//'"'
    str=trim(str)//' size="'//trim(inttochar(size(points)))//'"'

    WRITE(fileid,tab()) TAG('grid','B',str)
!    WRITE(fileid,tab()) "// XML file generated by Dirac"
!   Add the output of the date.

    indent=indent+3
    
    CALL XML_SECTION(fileid, points, 'gridpoints')
    CALL XML_SECTION(fileid, points, 'weight')
    CALL XML_SECTION(fileid, points, 'rho')
    CALL XML_SECTION(fileid, points, 'vc')
    CALL XML_SECTION(fileid, points, 'grad_rho')
    CALL XML_SECTION(fileid, points, 'hes_rho')
    
    indent=indent-3
    WRITE(fileid,tab()) TAG('grid',where='E')
    
    CLOSE(fileid)
  END SUBROUTINE WRITE_XML
  
  
  SUBROUTINE READ_GRID(file,points,npoints)
#include "priunit.h"
    TYPE(grid_data),POINTER :: points(:)
    INTEGER                 :: file
    INTEGER                 :: i,npoints
    
    REWIND(file)
    READ(file,*) npoints
    
    ALLOCATE(points(1:npoints))
    points(:)%vc=0.0
    points(:)%rho=0.0
    !        points(:)%grad_rho(:)=0.0
    !        points(:)%hes_rho=0.0
    
    DO i=1,npoints
       READ(file,*) points(i)%x(1),points(i)%x(2),points(i)%x(3),points(i)%weight
    END DO
    
  END SUBROUTINE READ_GRID
  
  
END MODULE GRID_FILE





