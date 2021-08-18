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
!       This module is used to construct the natural orbital based on quaternion density matrix.
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The natural orbitals and occupation numbers are written into 'DFNOSAO' by Module3.
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
            call write_no_to_newfile (exa_input,NOs_AO,N_OCCUPATION,'DFCOEF','DFNOSAO')
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

                  integer                            :: i,j,k,q,sz_print
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

                  flb = 0                                    
                  sz = NORBT
                  sz_occ = N_CORR_OCC
                  call quater_read_from_dirac_orbital(HF_AO,sz,sz_occ,flb)
                  call print_date('Read HF orbital coefficient from DFCOEF')

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

!-----------------------------------------------------------------------------------------------------------------------------
      subroutine write_no_to_newfile (exa_input,NOs_AO,N_OCCUPATION,file_old,file_new)
            !      write new orbital information into new file 'DFNOSAO', which is similar to DFCOEF 
            !      copy the basis information from DFCOEF and use new occupation number and orbital coefficient 
            !      to replace the old orbital. 
            !      This file only includes the orbital which is correlated in CC. 
            !      revise subroutine write_mos_to_dirac of exacorr_dirac.F
            
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
                   
                   
                   do i = 1, nmo
                      qorbital%energy(i)=N_OCCUPATION(i) !we use occupation number to replace the orbital energy
                   end do
 
                             
                   do  j = 1, nmo
                      do  i = 1, qorbital%nao
                            qorbital%coeff_q(i,j,1)=NOs_AO(i,j,1)
                            qorbital%coeff_q(i,j,2)=NOs_AO(i,j,2)
                            qorbital%coeff_q(i,j,3)=NOs_AO(i,j,3)
                            qorbital%coeff_q(i,j,4)=NOs_AO(i,j,4)
                      enddo
                   enddo
          
          
                   call quater_write_mos_to_dirac(ao_basis,qorbital,file_new)
                   call print_date('Write natural orbitals into DFNOSAO file like DFCOEF')
          
      end subroutine write_no_to_newfile
!----------------------------------------------------------------------------------------------------------

end module mp2no_driver

