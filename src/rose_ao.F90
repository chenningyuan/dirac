!This module defines the atomic orbital data types used in the ROSE program
       module rose_ao

        implicit none
        private

        public :: delete_adf_objects, dealloc_ao_basis, print_basis
        public :: collapse_ao_centers, match_center

        integer, parameter, private:: REALD=8 ! default real for this module

!       Basis function info:
        type, public:: gaussian_func_info_t
         integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        ! atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        ! atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        ! number of primitives in the contraction
         real(REALD):: coord(1:3)        ! coordinate of the basis function center (typically a nucleus)
         real(REALD), pointer :: exponent(:)
         real(REALD), pointer :: coefficient(:)
        end type gaussian_func_info_t

!       Basis set info:
        type, public:: basis_set_info_t

         ! general
         integer:: basis_type=1          ! 1=gaussian, 2=slater
         integer:: nao=0                 ! number of basis functions

         ! specific for Gaussians (DIRAC, pyscf, psi4, ...)
         integer:: basis_angular=1       ! 1=cartesian, 2=spherical
         integer:: nshells=0             ! number of shells in the basis set
         type(gaussian_func_info_t), pointer:: gtos(:)
         integer, pointer :: shell_indices(:)

         ! optional small-component Gaussian basis paired with the large-component basis
         logical :: has_small_component = .false.
         integer :: nao_small = 0
         integer :: nshells_small = 0
         integer :: basis_angular_small = 1
         type(gaussian_func_info_t), pointer :: gtos_small(:)
         integer, pointer :: shell_indices_small(:)

        end type basis_set_info_t

        contains

           integer function match_center (coord,center_coord,n_centers,tol)
!          Return the index of the center lying within tol of coord, or 0 if there
!          is none. The test is Euclidean distance; tol is in the same units as the
!          coordinates (bohr internally).
              real(REALD), intent(in) :: coord(3)         ! coordinate to look up
              real(REALD), intent(in) :: center_coord(:,:)! (3,>=n_centers) known centers
              integer,     intent(in) :: n_centers        ! number of centers currently defined
              real(REALD), intent(in) :: tol              ! matching tolerance
              integer     :: i
              real(REALD) :: tol2

              tol2 = tol * tol
              match_center = 0
              do i = 1, n_centers
                 if (sum((coord - center_coord(1:3,i))**2) <= tol2) then
                    match_center = i
                    return
                 end if
              end do
           end function match_center

           subroutine collapse_ao_centers (basis,tol,n_centers,center_coord, &
                                           shell_to_center,shell_to_center_small)
!          Collapse the AO shell coordinates onto the list of distinct centers that
!          carry a basis. This list, rather than the nucleus list, is the correct
!          enumeration for atomic IAO fragments: a ghost center carries a basis and
!          appears here, while a point nucleus carries no basis and does not.
!
!          Centers are numbered by first appearance in the shell list, so for a
!          conventional molecule the ordering is identical to the nucleus ordering.
!
!          When a small component is present every small-component shell must sit on
!          one of the large-component centers. Verifying that and building
!          shell_to_center_small are the same loop, so the check costs nothing beyond
!          the map that is needed anyway.
              type(basis_set_info_t),   intent(in)  :: basis
              real(REALD),              intent(in)  :: tol
              integer,                  intent(out) :: n_centers
              real(REALD), allocatable, intent(out) :: center_coord(:,:)
              integer,     allocatable, intent(out) :: shell_to_center(:)
              integer,     allocatable, intent(out) :: shell_to_center_small(:)

              integer                  :: ishell, icen
              real(REALD), allocatable :: scratch(:,:)

              allocate(scratch(3,basis%nshells))
              allocate(shell_to_center(basis%nshells))
              n_centers = 0

              do ishell = 1, basis%nshells
                 icen = match_center(basis%gtos(ishell)%coord,scratch,n_centers,tol)
                 if (icen == 0) then
                    n_centers = n_centers + 1
                    scratch(1:3,n_centers) = basis%gtos(ishell)%coord
                    icen = n_centers
                 end if
                 shell_to_center(ishell) = icen
              end do

              allocate(center_coord(3,n_centers))
              center_coord = scratch(1:3,1:n_centers)
              deallocate(scratch)

              if (basis%has_small_component) then
                 allocate(shell_to_center_small(basis%nshells_small))
                 do ishell = 1, basis%nshells_small
                    icen = match_center(basis%gtos_small(ishell)%coord,center_coord,n_centers,tol)
                    if (icen == 0) then
                       write(*,'(A,3F14.8)') ' Small-component shell at ', &
                          basis%gtos_small(ishell)%coord
                       write(*,'(A)') ' has no matching large-component center.'
                       stop 'collapse_ao_centers: small-component basis is inconsistent with the large component'
                    end if
                    shell_to_center_small(ishell) = icen
                 end do
              else
                 allocate(shell_to_center_small(0))
              end if

           end subroutine collapse_ao_centers

           subroutine delete_adf_objects (basis)
              type(basis_set_info_t), intent(inout) :: basis
           end subroutine delete_adf_objects

           subroutine dealloc_ao_basis (basis)
              type(basis_set_info_t), intent(inout) :: basis
              integer :: i
              if (basis%basis_type == 1) then
                 do i = 1, basis%nshells
                    deallocate(basis%gtos(i)%exponent)
                    deallocate(basis%gtos(i)%coefficient)
                    nullify(basis%gtos(i)%exponent)
                    nullify(basis%gtos(i)%coefficient)
                 end do
                 deallocate(basis%gtos,basis%shell_indices)
                 nullify(basis%gtos)
                 nullify(basis%shell_indices)
              end if
              if (basis%has_small_component) then
                 do i = 1, basis%nshells_small
                    deallocate(basis%gtos_small(i)%exponent)
                    deallocate(basis%gtos_small(i)%coefficient)
                    nullify(basis%gtos_small(i)%exponent)
                    nullify(basis%gtos_small(i)%coefficient)
                 end do
                 deallocate(basis%gtos_small,basis%shell_indices_small)
                 nullify(basis%gtos_small)
                 nullify(basis%shell_indices_small)
              end if
              basis%nao = 0
              basis%nshells = 0
              basis%has_small_component = .false.
              basis%nao_small = 0
              basis%nshells_small = 0
              basis%basis_angular_small = 1
           end subroutine dealloc_ao_basis

           subroutine print_basis (print_unit,basis)
              integer, intent(in) :: print_unit
              type(basis_set_info_t), intent(in):: basis
              integer:: i,j

              if (basis%basis_type==1) then
                 if (basis%basis_angular==1) then
                    write (print_unit, '(/A)') ' Cartesian Gaussian-type basis set'
                 elseif (basis%basis_angular==2) then
                    write (print_unit, '(/A)') ' Spherical Gaussian-type basis set'
                 else
                    write (print_unit, '(/A,I4)') ' WARNING: unknown angular type encountered: ',basis%basis_angular
                 end if
              elseif (basis%basis_type==2) then
                 write (print_unit, '(/A)') ' Slater-type basis set'
                 return
              else
                 write (print_unit, '(/A,I4)') ' WARNING: unknown basis type encountered: ',basis%basis_type
              end if

              write (print_unit, '(A,I6,A,I8)') ' Basis consists of ',basis%nshells,' shells and has dimension', &
              & basis%nao
              write (print_unit, '(/A,16X,A,17X,A,5X,A)') '  shell l-value','coordinates','exponent','coefficient'
              do i = 1, basis%nshells
                 write (print_unit,'(I6,5X,I1,3X,3F12.4,F16.4,F16.6)') i,basis%gtos(i)%orb_momentum-1,&
                 & basis%gtos(i)%coord,basis%gtos(i)%exponent(1),basis%gtos(i)%coefficient(1)
                 if (basis%gtos(i)%n_primitives>1) then
                    do j = 2, basis%gtos(i)%n_primitives
                       write (print_unit,'(51X,F16.4,F16.6)') basis%gtos(i)%exponent(j),basis%gtos(i)%coefficient(j)
                    end do
                 end if
              end do
              if (basis%has_small_component) then
                 write (print_unit, '(A,I6,A,I8)') ' Small component consists of ',basis%nshells_small, &
                 & ' shells and has dimension', basis%nao_small
                 write (print_unit, '(/A,16X,A,17X,A,5X,A)') '  shell l-value','coordinates','exponent','coefficient'
                 do i = 1, basis%nshells_small
                    write (print_unit,'(I6,5X,I1,3X,3F12.4,F16.4,F16.6)') i,basis%gtos_small(i)%orb_momentum-1,&
                    & basis%gtos_small(i)%coord,basis%gtos_small(i)%exponent(1),basis%gtos_small(i)%coefficient(1)
                    if (basis%gtos_small(i)%n_primitives>1) then
                       do j = 2, basis%gtos_small(i)%n_primitives
                          write (print_unit,'(51X,F16.4,F16.6)') basis%gtos_small(i)%exponent(j),basis%gtos_small(i)%coefficient(j)
                       end do
                    end if
                 end do
              end if


           end subroutine print_basis

       end module rose_ao
