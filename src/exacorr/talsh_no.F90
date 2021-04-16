module talsh_no
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use talsh_common_routines
        use exacorr_mo
        use exacorr_global,only: nfunctions
        use, intrinsic:: ISO_C_BINDING

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
        real(8), parameter    :: ONE_QUARTER=0.25D0

        type(talsh_tens_t) :: one_tensor

        private

        public talsh_no_driver

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_no_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.  
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan
        use exacorr_global
        type(exacc_input), intent(in) :: exa_input

        integer                       :: NORBT
      

        NORBT = exa_input%nkr_occ + exa_input%nkr_vir 
        call get_transform_Fock(NORBT)

        call print_date('Leave MP2_NO Module')
      end subroutine talsh_no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (NORBT)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO 
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix. 
         ! This routine doesn't truncate the new basis set.    
         ! Written by Xiang Yuan   

            use matrix_complex_to_quaternion 
            use write_orbital_to_file

                  integer, intent(in)                :: NORBT

                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  real(8), pointer                   :: density_quat(:,:,:)     ! quaternion MP2 density matrix 
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis
                  real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital 
                  integer                            :: sz,flb, fl              ! size of matrix and flag for block

!                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
!                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
!                  integer                            :: sz,flb                  ! size of matrix and flag for block
!                  character                          :: mat                     ! flag used for matrix basis transform
!                  integer                            :: i

                                                        
                  call quater_read(NORBT,density_quat)
                  call print_date('Read Quaternion Density matrix')

                  flb = 2                                      ! we read all orbital 
                  sz = exa_input%nkr_vir + exa_input%nkr_occ   
                  call quater_read_from_dirac_orbital(HF_AO,sz,flb)
                  call print_date('Read HF orbital coefficient from DFCOEF')

                  call diagoanl_matrix_quater(density_quat,N_OCCUPATION,NOs_MO)    
                  call print_date('Qdiag Density matrix')       

                  call sort_values_vectors (N_OCCUPATION,NOs_MO)
                  call print_date('Sort quaterbion natural orbital')

                  fl = 0 
                  call basis_transform(NOs_MO,HF_AO,NOs_AO,fl)
                  call print_date('Get new orbital in AO basis')

                  call write_to_newfile(NOs_AO,N_OCCUPATION,'NOs_MO_dirac_form')
                  call print_date('Write the all natural orbital and occupation number in file NOs_MO_dirac_form')

                    
!                  flb = 0                     ! we read vv block of Fock matrix
!                  sz = exa_input%nkr_vir      ! the size of the FM is the same as DM
!                  call quater_read_from_dirac_fock(Fock_HF,sz,flb)
!                  call print_date('Read Fock matrix from DFCOEF')

!                  mat = 'Matrix'
!                  call basis_transform(Fock_HF,NOs_MO,Fock_NO,mat)
!                  call print_date('Get full vv Fock matrix in natural orbital basis')

!                  call quater_write(Fock_NO,"FM_in_NatOrb")
!                  call print_date('Write the Full Fock matrix of natural orbital in file FM_in_NatOrb')
                        
      end subroutine get_transform_Fock

!-----------------------------------------------------------------------------------------------------------------------------

end module talsh_no

