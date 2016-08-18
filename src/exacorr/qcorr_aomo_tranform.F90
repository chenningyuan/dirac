!This module drives the index transformation of two electron integrals
!integrating the ExaTensor library with the DIRAC code.
       module exacorr_ao_to_mo

        implicit none
        private
!PARAMETERS:
        integer, parameter, private:: REALD=8 !default real for this module
!DERIVED TYPES:
!VISIBILITY:
        public compute_mo_integrals

       contains
!-------------------------------------------------------------------------------------------
        subroutine compute_mo_integrals (nao, nmo, mo_list)
!This subroutine drives the AO to MO integral transformation using AO integrals from InteRest and the tensor contraction library from ExaTensor
!
         Use dirac_MO
!         Use talsh

         implicit none
         integer, intent(in) :: nao 
         integer, intent(in) :: nmo(4)       ! the size of the mo basis for each of the 4 indices (mulliken ordering)
         integer, intent(in) :: mo_list(:)   ! and their indices (4 lists consecutively)

         integer             :: ifrp = 1     ! we do not use symmetry, so ifrp is always 1
         integer             :: nz = 4       ! we do not use symmetry, so nz is always 4
         integer             :: i, imo, ierr
         type(MO)            :: qmo
         integer             :: ipqtoq(4)   ! only needed for the debug print, don't forget to remove later

!        Make some noise so that we know we are here
         print*, " Arrived in compute_mo_integrals routine"
!        Initialize library
!         ierr = talsh_init()
         
!        Create empty tensors for MO coefficients

!        Initialize tensors with MO coefficients

         do i = 1, 4
            ipqtoq(i) = i
         end do

         imo = 1 
         do i = 1, 4
            call create (qmo, nao, nmo(i), nz)
            call get_dirac_mo_coefficients (ifrp,qmo,mo_list(imo:imo+nmo(i)-1),nmo(i))

!           test that we have read the coefficients correctly
            print*,"Mo coefficients for index ",i
            CALL PRQMAT(qmo%coeff,nao,nmo(i),nao,nmo(i),nz,ipqtoq,6)

            call delete (qmo)
            imo = imo + nmo(i)
         end do

!        Initialize AO integral tensors

!        Shut down library
!         ierr = talsh_shutdown()

!        Make some noise so that we know we are leaving
         print*, " Leaving compute_mo_integrals routine"
         return
        end subroutine compute_mo_integrals

       end module exacorr_ao_to_mo
