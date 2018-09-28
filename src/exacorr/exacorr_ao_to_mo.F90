module exacorr_ao_to_mo

!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.

        use exatensor

        implicit none
        private
        public exacorr_compute_mo_integrals
        public test_integral_generation

       contains

#if defined (VAR_MPI)
        subroutine exacorr_compute_mo_integrals (nmo, mo_list)

!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor

!        Written by Lucas Visscher, march 2018

         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use exacorr_tensor_methods
         use module_interest_interface
         use module_interest_eri
         use exatensor
         use, intrinsic:: ISO_C_BINDING

#ifdef USE_MPI_MOD_F90
         use mpi
         implicit none
#else
         implicit none
#include "mpif.h"
#endif

         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: i, ierr
         integer             :: number_aos                 ! the size of the atomic orbital basis
         integer             :: nshells                    ! the number of shells (orbitals with the same center and l-value)
         type(basis_func_info_t), allocatable :: gto(:)    ! array with information about the shells

         type(subspace_basis_t):: basis_ao
         type(color_symmetry_t), allocatable :: color(:)
         integer(INTL):: nao_l,labs, ao_space_root
         class(h_space_t), pointer:: ao_space
         integer:: mpi_th_provided, my_rank,my_role, ao_space_id
         integer:: nao, lsh, l, nl

         type(tens_rcrsv_t):: aoint_tensor
!later   type(tens_rcrsv_t):: moint_tensor
         type(compute_2e_ao_tensor_t) :: aoint_calculator

!        Make some noise so that we know we are here
         call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
         print*, " Process", my_rank," arrived in compute_mo_integrals_exa routine"

!        Initialize libraries
         call interest_initialize
         call initialize_global_data

         print*, " Process", my_rank," leaves compute_mo_integrals_exa routine"
         call delete_global_data
         return ! quit here because the stuff below still requires debugging 

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)
         nao_l   = nao ! copy of nao to a long integer

!        Get the basis set information
         call fetch_gtos(1,nao,gto,nshells)

         !Create a basis for a vector space:
         call basis_ao%subspace_basis_ctor(nao_l,ierr)
         !This needs to be an array that is persistent on all nodes
         allocate(color(nshells))
         if(ierr.ne.0) call quit('subspace_basis_t.subspace_basis_ctor() failed!')

         labs = 0
         do lsh = 1, nshells !set basis functions
            nl   =  ncartesian(gto(lsh)%orb_momentum)
            call color(lsh)%color_symmetry_ctor(lsh) ! for the moment take one shell as a color
            do  l = 1, nl
                labs = labs + 1
                call basis_ao%set_basis_func(labs,BASIS_ABSTRACT,ierr,symm=color(lsh))
                if(ierr.ne.0) call quit('subspace_basis_t.set_basis_func() failed!')
            end do
         enddo

         call basis_ao%finalize(ierr)
         if(ierr.ne.0) call quit('subspace_basis_t.finalize() failed!')

         !Register a vector space:
         ierr=exatns_space_register('AO_space',basis_ao,ao_space_id,ao_space)
         if(ierr.ne.0) call quit('exatns_space_register() failed!')

         !Register the ao integral generation as a method
         ierr=exatns_method_register('Aoint_calculator',aoint_calculator)

!        Run ExaTENSOR within MPI_COMM_WORLD:
         ierr=exatns_start(MPI_COMM_WORLD)
         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role)
          if(my_role.eq.EXA_DRIVER) then

           !Create ao-tensor:
           ao_space_root=ao_space%get_root_id(ierr)
           ierr=exatns_tensor_create(aoint_tensor,'ao_tensor',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_C8)

           !Initialize ao-tensor:
           ierr=exatns_tensor_init(aoint_tensor,'Aoint_calculator')

           !Destroy ao-tensor:
           ierr=exatns_tensor_destroy(aoint_tensor)
          end if

          !Stop ExaTENSOR runtime:
          ierr=exatns_stop()

         else
          write(*,*) ' Process ',my_rank,' terminated with error ',ierr
         endif

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Process ",my_rank, " leaving compute_mo_integrals_exa routine"
         return
        end subroutine exacorr_compute_mo_integrals

#else
        subroutine exacorr_compute_mo_integrals (nmo, mo_list)
         implicit none
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)
         print*,  'called exatensor version of compute_mo_integrals in a serial version'
         print*,  'this functionality only works with MPI, we need to stop here'
         stop ' activate MPI'
        end subroutine exacorr_compute_mo_integrals

#endif

        subroutine test_integral_generation (print_level)

!This subroutine tests the parallel AO integral generation

!        Written by Lucas Visscher, march 2018

         use dirac_MO
         use exacorr_datatypes
         use exacorr_global
         use module_interest_interface
         use module_interest_eri

         implicit none
         integer, intent(in) :: print_level  ! 0: only show that we running this test
                                             ! 1: print only overall summary
                                             ! 2: print information for each batch

         integer             :: nao, number_aos            ! the size of the atomic orbital basis

         integer :: ierr
         integer :: i, j, k, l, ish, jsh, ksh, lsh, ijkl
         integer :: li,lj,lk,ll,ni,nj,nk,nl
         real(8) :: ei,ej,ek,el,ci,cj,ck,cl,xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
         real(8) :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
         real(8), parameter :: fijkl = 1.0
         real(8) :: batch_value, ddot
         real(8), allocatable :: batch_values(:,:,:,:)
         integer            :: nshells ! number of AO shells
         type(basis_func_info_t), allocatable :: gto(:)   ! array with basis function information
         real(8) :: relative_error

!        Make some noise so that we know we are here
         print*, " Testing AO integral generation by calculating all integrals once"

!        Initialize libraries
         call interest_initialize(print_level.gt.1)
         if (print_level > 0) print*, " Initialized interest library"
         call initialize_global_data
         if (print_level > 0) print*, " Initialized global data"

!        Retrieve basis set information from DIRAC common blocks
         nao     = number_aos(1)

!        Get AOs and allocate array to store batch values
         call fetch_gtos(1,nao,gto,nshells)
         allocate(batch_values(nshells,nshells,nshells,nshells))

         do lsh = 1, nshells
            ll   =  gto(lsh)%orb_momentum
            el   =  gto(lsh)%exponent(1)
            xl   =  gto(lsh)%coord(1)
            yl   =  gto(lsh)%coord(2)
            zl   =  gto(lsh)%coord(3)
            cl   =  gto(lsh)%coefficient(1)
            nl   =  ncartesian(ll)
            do ksh = 1, nshells
               lk   =  gto(ksh)%orb_momentum
               ek   =  gto(ksh)%exponent(1)
               xk   =  gto(ksh)%coord(1)
               yk   =  gto(ksh)%coord(2)
               zk   =  gto(ksh)%coord(3)
               ck   =  gto(ksh)%coefficient(1)
               nk   =  ncartesian(lk)
               do jsh = 1, nshells
                  lj   =  gto(jsh)%orb_momentum
                  ej   =  gto(jsh)%exponent(1)
                  xj   =  gto(jsh)%coord(1)
                  yj   =  gto(jsh)%coord(2)
                  zj   =  gto(jsh)%coord(3)
                  cj   =  gto(jsh)%coefficient(1)
                  nj   =  ncartesian(lj)
                  do ish = 1, nshells
                     li   =  gto(ish)%orb_momentum
                     ei   =  gto(ish)%exponent(1)
                     xi   =  gto(ish)%coord(1)
                     yi   =  gto(ish)%coord(2)
                     zi   =  gto(ish)%coord(3)
                     ci   =  gto(ish)%coefficient(1)
                     ni   =  ncartesian(li)
                     call interest_eri('llll',fijkl,gout,&
                                  lk,ek,xk,yk,zk,ck,&
                                  ll,el,xl,yl,zl,cl,&
                                  li,ei,xi,yi,zi,ci,&
                                  lj,ej,xj,yj,zj,cj )
                     batch_value = ddot(ni*nj*nk*nl,gout,1,gout,1)
                     batch_values(ish,jsh,ksh,lsh) = batch_value
                     if (print_level > 1) then
                        write (*,'(a,4i6,e20.8)') 'batch',ish,jsh,ksh,lsh,batch_value 
                     endif
                  end do
               end do
            end do
         end do

!        Test correctness by checking permutational symmetry
         ierr = 0
         do lsh = 1, nshells
            do ksh = 1, lsh
               do jsh = 1, lsh
                  do ish = 1, jsh
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(jsh,ish,ksh,lsh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #1',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(jsh,ish,ksh,lsh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(jsh,ish,ksh,lsh), &
                        ish,jsh,ksh,lsh
                     end if
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(ish,jsh,lsh,ksh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #2',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(ish,jsh,lsh,ksh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(ish,jsh,lsh,ksh), &
                        ish,jsh,ksh,lsh
                     end if
                     if (relative_error(batch_values(ish,jsh,ksh,lsh),batch_values(ksh,lsh,ish,jsh)).gt.1.e-10) then
                        ierr = ierr + 1
                        write (*,'(a,3e20.10,4i4)') 'error in perm. symmetry #3',  &
                        batch_values(ish,jsh,ksh,lsh),batch_values(ksh,lsh,ish,jsh), &
                        batch_values(ish,jsh,ksh,lsh)-batch_values(ksh,lsh,ish,jsh), &
                        ish,jsh,ksh,lsh
                     end if
                  end do
               end do
            end do
         end do

         if (ierr == 0) Then
            print*," passed permutational symmetry test"
         else
            print*,ierr," failures in permutational symmetry test"
         end if

         deallocate(gto)
         deallocate(batch_values)

!        Clean up global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*, " Leaving test_integral_generation routine"
         return

        end subroutine test_integral_generation

end module exacorr_ao_to_mo
