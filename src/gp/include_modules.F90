!-----------------------------------------------------------------------------------------
!
! Collection of Fortran90 modules encapsulating "old" Fortran77 include files
! (like dcbham.h) among Dirac source codes.
!
! Selected variables inside these h-include files are private and are
! accessible out of module scope only via dedicated subroutines (methods).
!
! Notation examples:  
! ------------------
!
!             use include_dcbgen_h 
!             use include_dcbham_h
!             use include_dgroup_h, only : get_IPQTOQ
!             use include_dcbdhf_h
!
!  function get_<variable>(parameter1,parameter2,...) - function to return include file (module) variable
!
!  subroutine set_<variable>(parameter1,parameter2,...) - method to change
!  include file (module) variable
!
!
!
! Written by Miro Ilias, August 2014, Banska Bystrica
!
!
!-----------------------------------------------------------------------------------------

module include_dcbgen_h
implicit none
private
#include "dcbgen.h"
!> specify public variables and public methods
public :: get_IPRGEN
contains
function get_IPRGEN() result (IPRGEN_out)
 implicit none
 integer :: IPRGEN_out
 IPRGEN_out = IPRGEN
end function get_IPRGEN
end module include_dcbgen_h

module include_dgroup_h
!> module for picking up/setting certain variables from dgroup.h
implicit none
private
#include"dgroup.h"
!> specify public variables and public methods
public :: get_IPQTOQ  
contains
function get_IPQTOQ(i,j) result (IPQTOQ_out)
 implicit none
 integer, intent(in) :: i,j
 integer :: IPQTOQ_out
 if (.not. (i>=1.and.i<=4.and.j>=0.and.j<=7)) then
   call quit('get_IPQTOQ(i,j): wrong values i,j !')
 endif
 IPQTOQ_out=IPQTOQ(i,j)
end function get_IPQTOQ
end module include_dgroup_h

module include_dcbham_h
!> module for picking up/setting certain variables from dcbham.h
implicit none
private
#include"dcbham.h"
!> specify public variables and public methods
public :: get_IPRHAM
contains
function get_IPRHAM() result (IPRHAM_out)
 implicit none
 integer :: IPRHAM_out
 IPRHAM_out = IPRHAM
end function get_IPRHAM
end module include_dcbham_h

module include_dcbdhf_h
implicit none
private
#include "dcbdhf.h"
!> specify public variables and public methods
public :: get_IPRSCF
contains
function get_IPRSCF() result (IPRSCF_out)
 implicit none
 integer :: IPRSCF_out
 IPRSCF_out = IPRSCF
end function get_IPRSCF
end module include_dcbdhf_h

module include_dcbprp_h
implicit none
private

!> needed for parameters used in dcbprp.h
#include "mxcent.h"

#include "dcbprp.h"
!> specify public variables and public methods
public :: get_IPRPRP
contains
function get_IPRPRP() result (IPRPRP_out)
 implicit none
 integer :: IPRPRP_out
 IPRPRP_out = IPRPRP
end function get_IPRPRP
end module include_dcbprp_h
