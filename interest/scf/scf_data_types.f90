! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Module:       Kramers unrestricted SCF module
!
! Description:  Public derived data types for the module 
!
! Contains:     
!
! Licensing:    dirac_copyright_start
!                     Copyright (c) 2010 by the authors of DIRAC.
!                     All Rights Reserved.
!                     
!                     This source code is part of the DIRAC program package.
!                     It is provided under a written license and may be used,
!                     copied, transmitted, or stored only in accordance to the
!                     conditions of that written license.
!                     
!                     In particular, no part of the source code or compiled modules may
!                     be distributed outside the research group of the license holder.
!                     This means also that persons (e.g. post-docs) leaving the research
!                     group of the license holder may not take any part of Dirac,
!                     including modified files, with him/her, unless that person has
!                     obtained his/her own license.
!                     
!                     For information on how to get a license, as well as the
!                     author list and the complete list of contributors to the
!                     DIRAC program, see: http://dirac.chem.vu.nl
!               dirac_copyright_end
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Logs:         
!
! Comments:     Note: do not change the order of components in derived data types
!                     unless you appropriately modify the intrinsic structure constructor 
!                     for the derived data type initialisation
!                            
! ------------------------------------------------------------------------------------

MODULE SCF_DATA_TYPES

  implicit none

  public module_interface


  !> definition of "scf" data type
  type type_scf 
       integer :: ndim
       integer :: nshl
       integer :: naos
       integer :: naoc
       integer :: max_scf_iterations
       integer :: max_diis_dimension
       real(8) :: damping_factor 
       real(8) :: diis_error_treshold
       real(8) :: energy_change_treshold
  end type type_scf
  type (type_scf) :: scf


  !> definition of "element" data type
  type type_element 
       character(len=7) :: name 
       real(8)          :: mass  
       real(8)          :: number
       character(len=4) :: symbol 
  end type type_element 
  type (type_element)   :: element(1:109)


  !> definition of "atom" data type 
  type type_atom
       real(8) :: charge 
       real(8) :: coordinate_x
       real(8) :: coordinate_y
       real(8) :: coordinate_z
       real(8) :: gnu_exponent 
  end  type
  type(type_atom), allocatable :: atom(:)


  !> definition of "gto" data type
  type type_gto
       integer :: index
       integer :: offset
       integer :: lvalue
       integer :: origin
       integer :: sdegen 
       integer :: cdegen 
       real(8) :: norm(1)
       real(8) :: exponent(1)
       real(8) :: coefficient(1)
  end  type
  type(type_gto), allocatable, target :: gto(:)


  !> definition of "real_vector" data type
  type real_vector 
       real(8), pointer :: elms(:) 
  end  type


  !> definition of "real_matrix" data type 
  type real_matrix
       integer          :: nrow 
       integer          :: ncol  
       real(8), pointer :: elms(:,:) 
  end  type


  !> definition of "complex_pointer" data type 
  type complex_pointer
       complex(8), pointer :: p => NULL()
  end  type


  !> definition of "complex_matrix" data type 
  type complex_matrix
       integer                 :: nrow 
       integer                 :: ncol  
       type(complex_pointer)   :: info(4,4)  
       complex(8), allocatable :: elms(:,:,:) 
  end  type


  !> naos(:) - number of atomic orbitals, spherical coordinates
  !> naoc(:) - number of atomic orbitals, cartesian coordinates 
  integer :: naos(1:4)
  integer :: naoc(1:4)


  !> occupation 
  integer :: nocc(2)


CONTAINS

! ------------------------------------------------------------------------------------
!> Interface to DIRAC
  SUBROUTINE module_interface()

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbdhf.h"
#include "nucdata.h"

    !> initialization of local variables
       integer :: i
       integer :: j
       integer :: ier
       type( type_gto ) :: tmpgto
    !> end
  
  
    !> test if the large component basis is uncontracted
       if( nplrg.ne.nlarge )then
         write(6,*) ' Error: contracted basis set found => stop!'
         write(6,*) ' Note:  uncontracted basis sets are only supported for kramers_unrestricted_scf_solver'
         stop
       endif
    !> end
  

    !> interface scf settings 
       !scf = type_scf( 1, 1, 1, 1, maxitr, mxdiis, scfcnv(1), scfcnv(1) )
       scf = type_scf( 1, 1, 1, 1, maxitr, mxdiis, dampfc, 1.d-9, 1.d-9 )
    !> end

  
    !> interface nuclear data
       do i=1,109 
         element(i) = type_element( 'not_set', 0.0d0, nucmass(i), nuclabel(i) )
       enddo
    !> end
  
  
    !> interface molecular data
       allocate( atom(nucind), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: atom(:)'
       do i=1,size(atom)
         atom(i) = type_atom( charge(i), cord(1,i), cord(2,i), cord(3,i), gnuexp(i) )
       enddo
    !> end
  
  
    !> interface basis set data
       allocate( gto(nlrgsh), stat=ier ); if( ier.ne.0 )stop ' Error in allocation: gto(:)'
       i = 1
       do while( i <= size(gto) )  
         do j=1,nrco(i)
           gto(i)  = type_gto( i, 0, nhkt(i), ncent(i), 0, 0, priccf(i,j), priexp(i), priccf(i,j) )
           i = i + 1
         enddo
       enddo
    !> end


    !> reorder basis functions to the angular momentum 
    !> ascending order to increase the integral performance
       do i=1,size(gto)-1
         do j=(i+1),size(gto)
           if( gto(j)%lvalue >= gto(i)%lvalue )cycle
           tmpgto = gto(i)
           gto(i) = gto(j)
           gto(j) = tmpgto
         enddo
       enddo
    !> end


    !> set some data + AO-vector offsets 
       naos(1) = 0
       naoc(1) = 0
       do i=1,size(gto)
         gto(i)%offset = naoc(1) 
         gto(i)%sdegen = 2*gto(i)%lvalue - 1 
         gto(i)%cdegen = ( gto(i)%lvalue*(gto(i)%lvalue+1) )/2 
         naos(1) = naos(1) + gto(i)%sdegen 
         naoc(1) = naoc(1) + gto(i)%cdegen
       enddo
       do i=2,4
         naos(i) = i*naos(1)
         naoc(i) = i*naoc(1)
       enddo
    !> end


    !> fixme: temporary set some data to scf
       scf%ndim = naoc(1)                           ! matrix dimension
       scf%nshl = 2                                 ! always open-shells
       nocc(1)  = nelect_dhf/2 + naelec_dhf 
       nocc(2)  = nelect_dhf/2
       print*,' num. alpha elec = ',nocc(1),' num. beta elec = ',nocc(2)
    !> end


    print*,'NELEC_DHF(1) = ',NELEC_DHF(1)
    print*,'NELEC_DHF(2) = ',NELEC_DHF(2)
    print*,'NELECT_DHF   = ',NELECT_DHF
    print*,'NAELEC_DHF   = ',NAELEC_DHF
!   stop


    do i=1,size(gto)
      write(6,'(1x,5(a,i3,3x),3(a,d15.5,3x))') ' i = ',i,                            &
                                               ' index  = ',gto(i)%index,            &
                                               ' offset = ',gto(i)%offset,           &
                                               ' lvalue = ',gto(i)%lvalue,           &
                                               ' origin = ',gto(i)%origin,           &
                                               ' norm = ',  gto(i)%norm(1),          &
                                               ' exponent = ',gto(i)%exponent(1),    &
                                               ' coefficient = ',gto(i)%coefficient(1)
    enddo

  
!!! write(6,*) 
!!! write(6,'(2x,a,i5)') 'Basis set informations:'
!!! write(6,'(2x,a,i5)') '-----------------------'
!!! write(6,*) 
!!! write(6,*) 
!!! write(6,'(2x,a)')    'label    atom   charge   prim   cont   L-basis   RKB     S-basis'
!!! write(6,'(2x,70a1)') ('-',i=1,70)
!!!
!!! write(6,'(2x,70a1)') ('-',i=1,70)
    write(6,*) 
    write(6,'(2x,a,i5  )') 'Total number of basis function shells:    ',size(gto)
    write(6,'(2x,a,i5,a)') 'Total number of spherical basis functions:',naos(1),' (used for calculation)' 
    write(6,'(2x,a,i5  )') 'Total number of cartesian basis functions:',naoc(1)
    write(6,*) 
  
  END SUBROUTINE

END MODULE
