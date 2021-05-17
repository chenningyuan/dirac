module get_natural_orbital
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

        public get_no_driver

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine get_no_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh.
!       The complex density matrix are saved on the file 'DM_complex'.  
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan
        use exacorr_global, only : get_nsolutions
        type(exacc_input), intent(in) :: exa_input

        integer                :: NORBT, NROOT
        real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
        real(kind=8), pointer  :: N_OCCUPATION(:)         ! occupation number from qdia of DM
        logical                :: alive
      

        inquire(file="CCDENS", exist=alive)
        if (alive) then
            call print_date('File CCDENS existed')
            NORBT = exa_input%nkr_occ + exa_input%nkr_vir 
            NROOT = get_nsolutions() / 2
            call get_transform_Fock(NROOT,NORBT,NOs_AO,N_OCCUPATION)
            call write_no_to_newfile (exa_input,NOs_AO,N_OCCUPATION,'DFCOEF','GENOS')
            call print_date('Generated Natural orbital')
        else
            stop "Error: density matrix File CCDENS not found" 
        end if 
      end subroutine get_no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (NROOT,NORBT,NOs_AO,N_OCCUPATION)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO 
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix. 
         ! This routine doesn't truncate the new basis set.    
         ! Written by Xiang Yuan   

            use matrix_complex_to_quaternion 
            use write_orbital_to_file

                  integer, intent(in)                :: NROOT,NORBT

                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  real(8), pointer                   :: density_quat_read(:,:,:)! quaternion MP2 density matrix read from CCDENS 
                  real(8), pointer                   :: density_quat(:,:,:)     ! quaternion MP2 density matrix 
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis
                  real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital 
                  integer                            :: sz,flb, fl              ! size of matrix and flag for block

!                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix 
!                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
!                  integer                            :: sz,flb                  ! size of matrix and flag for block
!                  character                          :: mat                     ! flag used for matrix basis transform
                  integer                            :: i,j,k
                  real                               :: trace,compact

                                                
                  allocate(density_quat(NORBT,NORBT,4))
                  allocate(density_quat_read(NROOT,NROOT,4))
                  call quater_read(NROOT,density_quat_read)
                  call print_date('Read Quaternion Density matrix')
                  trace = 0
      !            do i=1, size(density_quat,1)
      !                  trace = trace+density_quat(i,i,1)
      !            end do

      !            print *, "print DM_trace in talsh_no", trace
                  
                  do i=1, NORBT
                     do j =1, NORBT
                        density_quat(i,j,1) = density_quat_read(i,j,1)
                        density_quat(i,j,2) = density_quat_read(i,j,2)
                        density_quat(i,j,3) = density_quat_read(i,j,3)
                        density_quat(i,j,4) = density_quat_read(i,j,4)
                     end do
                  end do

                  flb = 2                                      ! we read all orbital 
!                  sz = exa_input%nkr_vir + exa_input%nkr_occ   
                  sz = NORBT
                  call quater_read_from_dirac_orbital(HF_AO,sz,flb)
                  call print_date('Read HF orbital coefficient from DFCOEF')

                  call diagoanl_matrix_quater(density_quat,N_OCCUPATION,NOs_MO)    
                  call print_date('Qdiag Density matrix')       

!                  call sort_values_vectors (N_OCCUPATION,NOs_MO)
!                  call print_date('Sort quaterbion natural orbital')

                  fl = 0 
                  call basis_transform(NOs_MO,HF_AO,NOs_AO,fl)
                  call print_date('Get new orbital in AO basis')

                  print *, " Print several orbitals' occupation number"
                  trace = 0.0
                  compact = 0.0 
                  do i=1, size(N_OCCUPATION)
                        if (N_OCCUPATION(i) < 0.0) then 
                              N_OCCUPATION(i) = N_OCCUPATION(i) + 2.0
                              compact = compact + N_OCCUPATION(i)
                        end if 
                        trace = trace + N_OCCUPATION(i)
                        print *, ' Orbital',i,'  :Occupation number: ', N_OCCUPATION(i)
                  end do

                  print *, "Sum of occupation number which should be the same as total electron number", trace
                  print *, "Sum of large component of occupation number", compact

                    
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
      subroutine write_no_to_newfile (exa_input,NOs_AO,N_OCCUPATION,file_old,file_new)
            !      write new orbital information into new file 'MP2NOs_AO', which is similar to DFCOEF 
            !      copy the basis information from DFCOEF and use new orbital energy and coefficient to replace the old orbital.   
            !      revise subroutine write_mos_to_dirac of exacorr_dirac.F
          
                   use exacorr_mo
                   use exacorr_datatypes
                   use exacorr_utils
                   use write_orbital_to_file
            
                   implicit none
              
                   type(exacc_input), intent(in)         :: exa_input
                   character(len=*),intent(in)           :: file_old,file_new
                   real(8),pointer,intent(in)            :: NOs_AO(:,:,:)          ! new canonical orbital in AO basis
                   real(kind=8),pointer,intent(in)       :: N_OCCUPATION(:)        ! occupation number
!                   real(kind=8),pointer,intent(in)       :: New_orbital_energy(:)  ! new orbital energy 
          
          
                   type(basis_set_info_t)   :: ao_basis
                   type(qumo)               :: qorbital        ! save the new canonical orbital 
                   type(qumo)               :: qorbital_hf     ! save the original HF canonical orbital 
          
                   integer                  :: nocc            ! occupied orbital 
                   integer                  :: nao, nmo, nno   ! number of atomic orbital; molecular orbital; natural orbital
                   integer,allocatable      :: index(:)
                   integer                  :: i,j
                   integer                  :: ierr
                   integer                  :: nmo_hf
                   real(kind=8)             :: diff_1, diff_2, diff_3, diff_4
          
          
          !        read information of from DFCOEF
                   call quater_read_from_dirac_all(ao_basis,qorbital_hf,file_old,ierr)
          
                   nocc = exa_input%nkr_occ
                   nao = qorbital_hf%nao
                   nno = size(NOs_AO,2)
                   nmo = nocc+nno
                   nmo_hf = qorbital_hf%nmo
                   allocate(index(nmo))
                   do i = 1, nmo
                      index(i) = i ! we set the index compatible with the fact that we store all mo's
                   end do
          
                   call alloc_mo(qorbital,nao,nmo_hf,4)
                   qorbital%restricted = qorbital_hf%restricted
                   qorbital%complex_algebra = qorbital_hf%complex_algebra
                   qorbital%total_energy = qorbital_hf%total_energy
                   qorbital%fermion_irrep = qorbital_hf%fermion_irrep
                   qorbital%boson_irrep = 0
                   qorbital%index = index
                   
                   
                   do i = 1, nocc
                      qorbital%energy(i)=N_OCCUPATION(i) !we use occupation number to replace the orbital energy
                   end do
 
                             
                   do  j = 1, nmo
                      do  i = 1, qorbital%nao
                            qorbital%coeff(i,j,1)=NOs_AO(i,j,1)
                            qorbital%coeff(i,j,2)=NOs_AO(i,j,2)
                            qorbital%coeff(i,j,3)=NOs_AO(i,j,3)
                            qorbital%coeff(i,j,4)=NOs_AO(i,j,4)
                      enddo
                   enddo
          
          
                   call quater_write_mos_to_dirac(ao_basis,qorbital,file_new)
                   call print_date('Write natural orbitals into DFNOSAO file like GENOS')
          
      end subroutine write_no_to_newfile
!----------------------------------------------------------------------------------------------------------

end module get_natural_orbital

