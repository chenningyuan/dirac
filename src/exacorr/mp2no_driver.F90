module mp2no_driver
!This module contains the routines needed to calculate the natural orbitals, provided quaternion density matrix.  

        use exacorr_datatypes
        use exacorr_global,only: get_nsolutions
        use mp2no_matrixtransforms 
        use mp2no_fileinterfaces
        use, intrinsic:: ISO_C_BINDING

        implicit none

        public get_no_driver

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine get_no_driver (exa_input)
!       This module is used to construct the natural orbital by reading and diagonalizing a quaternion density matrix.
!       Written by Xiang Yuan
        type(exacc_input), intent(in) :: exa_input

        integer                :: NORBT, NROOT, N_CORR_OCC
        real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
        real(kind=8), pointer  :: N_OCCUPATION(:)         ! occupation number from qdia of DM
        logical                :: alive
      

        call print_date('Entering get_no_driver')
        inquire(file="CCDENS", exist=alive)
        if (alive) then
            call print_date('File CCDENS exists')
            NORBT = exa_input%nkr_occ + exa_input%nkr_vir 
            NROOT = get_nsolutions() / 2
            N_CORR_OCC = exa_input%mokr_occ(1)
            call get_transform_NO(NROOT,NORBT,N_CORR_OCC,NOs_AO,N_OCCUPATION)
            call print_date('Generated natural orbitals in AO basis')
            call write_nos (exa_input%nkr_occ,NOs_AO,N_OCCUPATION,'scf','mbpt')
        else
            stop "Error: density matrix File CCDENS not found" 
        end if 
        call print_date('Leaving get_no_driver')
      end subroutine get_no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_NO (NROOT,NORBT,N_CORR_OCC,NOs_AO,N_OCCUPATION)
         ! This routine read the quternion density matrix and do diagonalization.
         ! Then do basis transform to go to atomic basis set.  
         ! NROOT used in dimension of ccden file and norbt is the total correlated orbital number. 
         ! NROOT is always lager than NORBT 
         ! Written by Xiang Yuan   

                  integer, intent(in)                :: NROOT,NORBT,N_CORR_OCC

                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  real(8), pointer                   :: density_quat_read(:,:,:)! quaternion MP2 density matrix read from CCDENS 
                  real(8), pointer                   :: density_quat(:,:,:)     ! quaternion MP2 density matrix 
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis
                  real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital 
                  integer                            :: sz,flb, fl,sz_occ       ! size of matrix and flag for block
                  integer                            :: nocc, nvirt

                  type(cmo)                          :: qorbital
                  type(basis_set_info_t)             :: ao_basis
                  integer                            :: i,j,k,q,sz_print,ierr
                  real                               :: trace,compact
                  real(kind=8)                       :: temp

                                                
                  allocate(density_quat(NORBT,NORBT,4))
                  allocate(density_quat_read(NROOT,NROOT,4))
                  call quater_read(NROOT,density_quat_read)
                  call print_date('Read density matrix (quaternion basis)')
                  trace = 0
                  
                  do i=1, NORBT
                     do j =1, NORBT
                        density_quat(i,j,1) = density_quat_read(i,j,1)
                        density_quat(i,j,2) = density_quat_read(i,j,2)
                        density_quat(i,j,3) = density_quat_read(i,j,3)
                        density_quat(i,j,4) = density_quat_read(i,j,4)
                     end do
                  end do

                  nvirt = NORBT
                  nocc  = N_CORR_OCC
                  call read_from_dirac(ao_basis,qorbital,'scf',ierr)
                  call convert_mo_to_quaternion(qorbital)
                  allocate(HF_AO(qorbital%nao,nvirt,4))
                  HF_AO = qorbital%coeff_q(:,nocc+1:nocc+nvirt,:)
                  call print_date('Read scf orbital coefficients')

                  call diagonal_matrix_quater(density_quat,N_OCCUPATION,NOs_MO)    
                  call print_date('Diagonalized density matrix (quaternion basis)')       

                  call basis_transform_vector(NOs_MO,HF_AO,NOs_AO)
                  call print_date('Get new orbital in AO basis')

                  trace = 0.0
                  compact = 0.0 
                  do i=1, size(N_OCCUPATION)
                        if (N_OCCUPATION(i) < 0.0) then 
                              N_OCCUPATION(i) = N_OCCUPATION(i) + 2.0
                              compact = compact + N_OCCUPATION(i)
                        end if 
                        trace = trace + N_OCCUPATION(i)
                  end do


                  sz_print = size(N_OCCUPATION)
                  do i = 1, sz_print
                     do j =1, sz_print-i
                        if (N_OCCUPATION(j) < N_OCCUPATION(j+1)) then
                              temp = N_OCCUPATION(j)
                              N_OCCUPATION(j) = N_OCCUPATION(j+1)
                              N_OCCUPATION(j+1) = temp
                        end if
                     end do 
                  end do


                  print *, " Natural orbital occupation numbers"

                  k = size(N_OCCUPATION)/5
                  q = mod(size(N_OCCUPATION),5)
                  do i=1, k, 1
                        j= (i-1)*5+1
                        write(*,"(F10.7,3X,F10.7,3X,F10.7,3X,F10.7,3X,F10.7)") N_OCCUPATION(j), N_OCCUPATION(j+1), &
                        N_OCCUPATION(j+2),N_OCCUPATION(j+3), N_OCCUPATION(j+4)
                  end do

                  j = k*5+1
                  select case (q)
                     case (1)
                        write(*,"(F10.7)") N_OCCUPATION(j)

                     case (2)
                        write(*,"(F10.7,3X,F10.7)") N_OCCUPATION(j), N_OCCUPATION(j+1)

                     case (3)
                        write(*,"(F10.7,3X,F10.7,3X,F10.7)") N_OCCUPATION(j), N_OCCUPATION(j+1), N_OCCUPATION(j+2)
                                                                                       
                     case (4)
                        write(*,"(F10.7,3X,F10.7,3X,F10.7,3X,F10.7)") N_OCCUPATION(j), N_OCCUPATION(j+1), N_OCCUPATION(j+2), &
                                                                               N_OCCUPATION(j+3)                                                                                          
                  end select                                                                                                                                                



                  print *, "Sum of occupation numbers (1) :",trace
                  print *, "Sum of occupation numbers of occupied orbitals (2)", compact
                  print *, "Note : (1) should equal the total number or correlated electrons"

                  deallocate(density_quat_read)
                        
      end subroutine get_transform_NO

end module mp2no_driver

