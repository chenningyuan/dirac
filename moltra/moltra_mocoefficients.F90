!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

module dirac_MO

!
! Contains all the necessary information for a range of molecular orbitals
! In particular:  - Orbital energy
!                 - Boson characterization (for spinfree runs)
!                 - Coefficients of AO-functions

   implicit none
  
!  Private

  Public :: create, delete, set_coefficients, get_coefficients

  type MO
     
    integer ::  &
      nao,      &          ! number of aos
      nmo,      &          ! number of mos
      nz                   ! number of active quaternion units

    real(8), pointer ::  &
      energy    (:),     & ! energy
      coeff (:,:,:)        ! coefficients (AO,MO,IZ)
     
    integer, pointer ::  &
      boson_irrep(:),    & ! for spinfree calculations: boson irrep to which the spinor belongs
      fermion_irrep(:)     ! fermion irrep (1=gerade, 2=ungerade) to which the spinor belongs

  end type MO

  interface create
    module procedure AllocMO
  end interface  

  interface delete
    module procedure DeAllocMO
  end interface  
  
  interface set_coefficients
    module procedure SetMOcoefficients
  end interface  

  interface get_coefficients
    module procedure GetMOcoefficients
  end interface  

contains

  subroutine AllocMO (M, nao, nmo, nz)
   
!   Prepare to work with a block of Mos

    type (MO), intent (out) :: M
    integer, intent(in)           :: nao
    integer, intent(in)           :: nmo
    integer, intent(in), optional :: nz

    M%nao = nao  ! Size of AO basis. Is kept here as we may want to use non-global defs
    M%nmo = nmo  ! Size of MO basis. May vary (e.g. in the 4-index when we have different ranges)
    if (present(nz)) then
       M%nz = nz
    else
       M%nz = 1  ! We default to one, easy for spinfree applications
    endif

!   Reserve the memory needed to work with a block of Mos
    nullify  (M%energy)
    allocate (M%energy(M%nmo))
    M%energy = 0.0

    nullify  (M%coeff)
    allocate (M%coeff(M%nao,M%nmo,M%nz))
    M%coeff = 0.0

    nullify  (M%boson_irrep)
    allocate (M%boson_irrep(M%nmo))
    M%boson_irrep = 0

    nullify  (M%fermion_irrep)
    allocate (M%fermion_irrep(M%nmo))
    M%fermion_irrep = 0
    
  end subroutine AllocMO

  subroutine DeAllocMO (M)

!   Clean up workspace

    type (MO), intent (inout) :: M

!   Set the dimensions to zero (not really necessary, but may prevent errors due to abuse of the MO-type)
    M%nao = 0
    M%nmo = 0
    M%nz = 0

!   Free the memory associated with a block of Mos
    deallocate (M%energy)
    deallocate (M%coeff)
    deallocate (M%boson_irrep)
    deallocate (M%fermion_irrep)
    
  end subroutine DeAllocMO

  subroutine SetMOcoefficients (M, coefficients)
   
!   Fill the block of MOcoefficients with values
!   Presently this is a simple copy, but more advanced routines that
!   operate on subsets of coefficients are also conceivable
!
    type (MO), intent (inout) :: M
    real(8), intent(in)       :: coefficients(M%nao,M%nmo,M%nz)

    M%coeff = coefficients

  end subroutine SetMOcoefficients

  subroutine GetMOcoefficients (M, coefficients)
   
!   Get a block of MOcoefficients
 
    type (MO), intent (in)   :: M
    real(8), intent(inout)   :: coefficients(M%nao,M%nmo,M%nz)

    coefficients = M%coeff

  end subroutine GetMOcoefficients

end module dirac_MO
