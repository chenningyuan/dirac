module talsh_mp2lap
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 
!using Cholesky decompostion and Laplace transformation

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils

        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)

        private

        public talsh_mp2lap_driver

       contains

!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_mp2lap_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

        implicit none

        type(exacc_input), intent(in) :: exa_input

        real(8)                       :: threshold
        type(talsh_tens_t)            :: chol_tensor

        call print_date("Entering MP2_LAP Module")
 
        threshold = exa_input%t_cholesky

!       perform cholesk decompostion

        call get_cholesky_vec (exa_input%nocc,exa_input%nvir,exa_input%mo_occ,exa_input%mo_vir,uov,exa_input%t_cholesky,exa_input%print_level)


        call print_date('Leaving MP2_LAP Module')
      end subroutine talsh_mp2lap_driver


end module talsh_mp2lap
