!This module provides auxiliary routines required for
!integrating the ExaTensor library with the DIRAC code.
       module dirac_exatensor
        implicit none
        private
!PARAMETERS:
 !Generic:
        integer, parameter, private:: REALD=8 !default real for this module
!DERIVED TYPES:
 !Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       !orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        !atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        !atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        !number of primitives in the contraction
         real(REALD):: coord(1:3)        !coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type basis_func_info_t
 !Molecular spinor info:
        type, public:: mo_spinor_info_t
         integer:: symmetry=0            !molecular spinor irreducible representation (0: no symmetry)
        end type mo_spinor_info_t
 !Tensor descriptor:
        type, public:: tens_descr_t
         integer:: rank=-1                          !number of tensor dimensions
         integer, pointer, contiguous:: seg_id(:)   !segment id for each tensor dimension (determines the orbital subrange)
         integer, pointer, contiguous:: seg_size(:) !segment size for each tensor dimension (in case it is reduced)
        end type tens_descr_t
!GENERIC INTERFACES:
 !Compute 1e Overlap integrals:
        interface compute_1e_overlap_integrals
         module procedure compute_1e_overlap_integrals_r8
        end interface compute_1e_overlap_integrals
 !Compute 1e Hamiltonian integrals:
        interface compute_1e_hamiltonian_integrals
         module procedure compute_1e_hamiltonian_integrals_r8  ! real for scalar relativistic
         module procedure compute_1e_hamiltonian_integrals_c8  ! complex for soc relativistic
        end interface compute_1e_hamiltonian_integrals
 !Compute 2e Coulomb type (1/r12) integrals:
        interface compute_2e_coulomb_integrals
         module procedure compute_2e_coulomb_integrals_r8
        end interface compute_2e_coulomb_integrals
!VISIBILITY:
        public generate_ao_segments
        public generate_mo_segments
        public compute_1e_overlap_integrals
        public compute_1e_hamiltonian_integrals
        public compute_2e_coulomb_integrals
        public get_mo_fock_diagonal
        public compress_unc_shell

       contains
!-------------------------------------------------------------------------------------------
        subroutine generate_ao_segments(num_large_ao,num_small_ao,desired_seg_size,ao_info,&
                                       &num_large_segs,num_small_segs,segs_large,segs_small,ierr)
!This subroutine partitions the range of large/small basis functions into segments
!such that the size of each segment is close to the <desired_seg_size> as much as possible.
!Different segments may have different sizes in order to align the segments to atomic shell boundaries.
         implicit none
         integer, intent(in):: num_large_ao                !in: total number of large component basis functions
         integer, intent(in):: num_small_ao                !in: total number of small component basis functions
         integer, intent(in):: desired_seg_size(0:2)       !in: desired size of a segment (0: target value, 1: minimal value, 2: maximum value)
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: num_large_segs             !out: number of large component segments
         integer, intent(out):: num_small_segs             !out: number of small component segments
         integer, intent(out), allocatable:: segs_large(:) !out: size of each segment in the large component range
         integer, intent(out), allocatable:: segs_small(:) !out: size of each segment in the small component range
         integer, intent(out):: ierr                       !out: error code (0:success)

!        place holder initializations for testing compilation and cmake inclusion of of this module
         num_large_segs = 1 
         allocate (segs_large(num_large_segs))
         segs_large(1) = num_large_ao
         num_small_segs = 0 
         ierr=0
         return
        end subroutine generate_ao_segments
!----------------------------------------------------------------------------------
        subroutine generate_mo_segments(num_occ,num_virt,desired_seg_size,mo_info,&
                                       &num_occ_segs,num_virt_segs,segs_occ,segs_virt,ierr)
!This subroutine partitions the range of occupied/virtual molecular spinors into segments
!such that the size of each segement is close to the <desired_seg_size> as much as possible.
!Different segements may have different sizes in order to align the segements by symmetry IR.
         implicit none
         integer, intent(in):: num_occ                    !in: total number of occupied spinors
         integer, intent(in):: num_virt                   !in: total number of virtual spinors
         integer, intent(in):: desired_seg_size(0:2)      !in: desired size of a segment (0: target value, 1: minimal value, 2: maximum value)
         type(mo_spinor_info_t), intent(in):: mo_info(1:) !in: information on each molecular spinor
         integer, intent(out):: num_occ_segs              !out: number of segments in the occupied range
         integer, intent(out):: num_virt_segs             !out: number of segments in the virtual range
         integer, intent(out), allocatable:: segs_occ(:)  !out: size of each segment in the occupied range
         integer, intent(out), allocatable:: segs_virt(:) !out: size of each segment in the virtual range
         integer, intent(out):: ierr                      !out: error code (0:success)

!        place holder initializations for testing compilation and cmake inclusion of of this module
         num_occ_segs  = 1 
         allocate (segs_occ(num_occ_segs))
         segs_occ(1) = num_occ
         num_virt_segs = 1 
         allocate (segs_virt(num_virt_segs))
         segs_virt(1) = num_virt
         ierr=0
         return
        end subroutine generate_mo_segments
!-------------------------------------------------------------------------------
        subroutine compute_1e_hamiltonian_integrals_r8(tens_descr,tens_body,ao_info,ierr)
!This subroutine computes a block of 1e Hamiltonian integrals in the AO basis.
!Argument <tens_descr> specifies the tensor slice for which the integrals to be computed.
         implicit none
         type(tens_descr_t), intent(in):: tens_descr       !in: tensor descriptor (specifies which Hamiltonian tensor slice to compute)
         real(REALD), intent(inout):: tens_body(*)             !out: Hamiltonian integrals
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: ierr                       !out: error code (0:success)

         ierr=0
         return
        end subroutine compute_1e_hamiltonian_integrals_r8
!-------------------------------------------------------------------------------
        subroutine compute_1e_hamiltonian_integrals_c8(tens_descr,tens_body,ao_info,ierr)
!This subroutine computes a block of 1e Hamiltonian integrals in the AO basis.
!Argument <tens_descr> specifies the tensor slice for which the integrals to be computed.
         implicit none
         type(tens_descr_t), intent(in):: tens_descr       !in: tensor descriptor (specifies which Hamiltonian tensor slice to compute)
         complex(16), intent(inout):: tens_body(*)         !out: Hamiltonian integrals
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: ierr                       !out: error code (0:success)

         ierr=0
         return
        end subroutine compute_1e_hamiltonian_integrals_c8
!-------------------------------------------------------------------------------
        subroutine compute_1e_overlap_integrals_r8(tens_descr,tens_body,ao_info,ierr)
!This subroutine computes a block of 1e overlap integrals in the AO basis.
!Argument <tens_descr> specifies the tensor slice for which the integrals to be computed.
         implicit none
         type(tens_descr_t), intent(in):: tens_descr       !in: tensor descriptor (specifies which Hamiltonian tensor slice to compute)
         real(REALD), intent(inout):: tens_body(*)             !out: Hamiltonian integrals
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: ierr                       !out: error code (0:success)

         ierr=0
         return
        end subroutine compute_1e_overlap_integrals_r8
!-------------------------------------------------------------------------------
        subroutine compute_2e_coulomb_integrals_r8(tens_descr,tens_body,ao_info,ierr)
!This subroutine computes a block of 2e Hamiltonian integrals in the AO basis.
!Argument <tens_descr> specifies the tensor slice for which the integrals to be computed.
         implicit none
         type(tens_descr_t), intent(in):: tens_descr       !in: tensor descriptor (specifies which Hamiltonian tensor slice to compute)
         real(REALD), intent(inout):: tens_body(*)             !out: Hamiltonian integrals
         type(basis_func_info_t), intent(in):: ao_info(1:) !in: information on each basis function
         integer, intent(out):: ierr                       !out: error code (0:success)

         ierr=0
         return
        end subroutine compute_2e_coulomb_integrals_r8
!------------------------------------------------------
        subroutine get_mo_fock_diagonal(fock_diag,ierr)
!This subroutine returns an array of diagonal Fock elements
!for the entire molecular spinor range.
         implicit none
         complex(REALD), intent(out), allocatable:: fock_diag(:) !diagonal fock elements
         integer, intent(out):: ierr                         !out: error code (0:success)

         ierr=0
         return
        end subroutine get_mo_fock_diagonal

        integer function compress_unc_shell (shell)
!If basis functions were grouped together in the input, they appear as quasi-contracted functions
!with a lot of zero coefficients. Remove these zeroes and reduce the length to just one function.
!Routine will produce an error, if the function is a true contraced function (more than one non-zero coefficient)

          integer i, j, non_zero 
          type(basis_func_info_t), intent(inout) :: shell
          real(REALD), allocatable :: exponents(:), coefficients(:)   
          real(REALD), parameter   :: treshold=1.D-12   

          allocate(exponents(shell%n_primitives))
          allocate(coefficients(shell%n_primitives))
          exponents    = shell%exponent
          coefficients = shell%coefficient

          non_zero = 0
          do i = 1, shell%n_primitives
             if (abs(coefficients(i)) > treshold) then
                non_zero = non_zero + 1
                j = i
             end if
          end do
         
          if (non_zero > 1) then
             print*, " WARNING: compress_unc_shell called for contracted shell"
             compress_unc_shell = 1
             return   ! leave shell untouched and return with error condition
          else if (non_zero == 0) then
             print*, " WARNING: compress_unc_shell found empty shell"
             compress_unc_shell = 2
             j = 1    ! set j to the first function, to avoid complications with zero length arrays
          end if

          deallocate (shell%exponent)
          deallocate (shell%coefficient)
          allocate   (shell%exponent(1))
          allocate   (shell%coefficient(1))
          
          shell%n_primitives = 1
          shell%exponent     = exponents(j)
          shell%coefficient  = coefficients(j)
      
           deallocate (exponents)
           deallocate (coefficients)

           compress_unc_shell = 0

        end function compress_unc_shell

       end module dirac_exatensor
