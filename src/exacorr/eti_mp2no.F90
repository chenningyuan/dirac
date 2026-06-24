module eti_mp2no
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        !use exacorr_utils
        use eti_common_routines
        use exacorr_mo
        use exacorr_global
        use exacorr_printout
        use exacorr_scaling_parameters
        use mp2no_matrixtransforms
        use mp2no_fileinterfaces
        use eti_cis
        use eti_ccrsp_common
        use, intrinsic:: ISO_C_BINDING

        implicit none

        type(tensor_t) :: one_tensor

        private

        public eti_mp2no_driver, get_new_orbital, get_transform_Fock

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine eti_mp2no_driver (exa_input,n_cano)
!       This module is used to construct the MP2 density matrix
!       The complex density matrix are saved on the file 'DM_complex'.
!       The natural orbital can then be generated in Module2 by a quaternion diagomalization.
!       The full natural orbital and Fock matrix in natural orbital basis are saved on file 'NOs_MO' and 'FM_in_NatOrb', respectively.
!       The new reduced canonical orbitals are written into 'MP2NOs_AO' by Module3.
!       Written by Xiang Yuan, modified by Johann Pototschnig

        use eti_mp2lap

        type(exacc_input), intent(in) :: exa_input
        integer, intent(inout)        :: N_cano                  ! number of new canonical orbital

!       Point connected with tensors
        integer                :: ierr
        real(8),pointer        :: NOs_AO(:,:,:)           ! new quaternion orbital in AO basis
        real(kind=8), pointer  :: N_OCCUPATION(:)         ! occupation number from qdia of DM
        real(kind=8), pointer  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
        logical                :: alive
        real(kind=8)           :: threshold
        integer                :: i, nocc
        character(len=10)      :: fn

        call print_date("Entering MP2_NO Module")

        threshold = exa_input%t_mp2no
        !LV: the line below should be improved, will not work if a user does not give the HOMO as the last orbital
        nocc = exa_input%mokr_occ(size(exa_input%mokr_occ))

        inquire(file="NOs_MO", exist=alive)

        if (alive) then
            call print_date('File NOs_MO exists')
        else
            call print_date('File NOs_MO does not exist. Start constructing MP2 density matrix')
            ! Generate the file DM_complex that stores the density matrix
            ! The file name DM_complex is hard-coded -> may need to put it in a config file
            if (exa_input%mp2no<9) then
                  call eti_generate_MP2_VNO_MObasis(exa_input)
            else
                  call mp2lap_get_DM(exa_input)
            endif
            call get_transform_Fock(exa_input)
        end if

        call get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION, &
                              New_orbital_energy,'NOs_MO',N_cano)
        call write_nos(nocc,NOs_AO,New_orbital_energy,'scf','mbpt')

        call print_date('Leaving MP2_NO Module')
      end subroutine eti_mp2no_driver


!------------------------------------------------------------------------------------------------------------------------------------
      subroutine eti_generate_MP2_VNO_MObasis (exa_input)
!       This subroutine obtains the complex mp2 density matrix in MO basis
!       The non-diagonal complex DM is written into file DM_complex.
!       Written by Xiang Yuan

         use checkpoint
         use eti_integrals
         use eti_cc, only : eti_init_fock

         type(exacc_input), intent(in)          :: exa_input
         integer                                :: nocc, nvir              ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable                  :: mo_vir(:)               ! the list index of virtual spinors
         integer,  allocatable                  :: mo_occ(:)               ! the list index of occupied spinors
         real(8), allocatable                   :: eps_occ(:),eps_vir(:)   ! HF orbital energies

!        Intermediate tensors (same notation as in relccsd)
         type(tensor_t)          :: MP2_vv_density
         type(tensor_t)          :: CIS_vv_density
         type(tensor_t)          :: result_tensor
         integer(INTD), dimension(4) :: vvoo_dims, oovv_dims
         integer(INTD), dimension(2) :: vv_dims
         integer(C_INT)              :: one_dims(1)

!        Common tensors
         type(eti_comm_tens_t) :: comm_t

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_intg_tens_t)       :: int_t

!        property integrals for PNO
!         type(eti_tens_shape_t)            :: tens_shape

!        Point connected with tensors
         integer(C_INT)              :: ierr
         complex(8), pointer         :: result_tens(:)
         real(8)                     :: mp2_energy,scf_energy

         type(C_PTR):: body_p

!        Scalars (need to be defined as tensor types)
         real(kind=8)  :: energy_shift(2) = 0.0d0
         integer       :: i, j, k, a, c

         integer :: num_gpus = 0, igpu
         character(len=20)                  :: fn, istr

         logical :: fock_flag ! used for FNO++ orbitals

         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir


!        Make some noise so that we know we are here
         call print_date('Entered mp2_driver routine')

         call eti_executor_initialize(exa_input%talsh_buff,get_default_executor())
!        Initialize libraries

         call print_date('Initializing interest library')
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         vv_dims(1) = nvir
         vv_dims(2) = nvir

         ierr=construct_tensor_t(one_tensor,exa_input%eti_precision%storage_t,one_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(MP2_vv_density,exa_input%eti_precision%storage_t,vv_dims,init_val=ZERO)

!        Get HF orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("retrieved orbital energies")

!        Get oovv integrals
         fock_flag = exa_input%mp2no.eq.3
         if (fock_flag) then
            int_t%oovv_active = .true.
            int_t%oooo_active = .true.
            int_t%ooov_active = .false.
            int_t%vovo_active = .true.
            int_t%vovv_active = .false.
            int_t%vvvv_active = .false.
            call get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t,1,exa_input%print_level,exa_input%eti_precision)
            call eti_init_fock(1, mo_occ, mo_vir, eps_occ, eps_vir, int_t, scf_energy, exa_input%print_level)
         else
            int_t%oovv_active = .true.
            int_t%oooo_active = .false.
            int_t%ooov_active = .false.
            int_t%vovo_active = .false.
            int_t%vovv_active = .false.
            int_t%vvvv_active = .false.
            call get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t,1,exa_input%print_level,exa_input%eti_precision)
         end if

!        Make some noise so that we know we are proceeding
         call print_date('started forming VV block of MP2 density matrix')
!        Initialize scale with denominators of oovv integrals
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc

         ierr=construct_tensor_t(comm_t%t2,exa_input%eti_precision%storage_t,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(result_tensor,exa_input%eti_precision%storage_t,one_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,one_dims)
!        Initialize t2 amplitudes with the MP2 values
         ierr=initialize_tensor_t(comm_t%t2)
         ierr=t_einsum("T(a,b,i,j)+=V+(i,j,a,b)",comm_t%t2,int_t%oovv,one_tensor)
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
!        Calculate and print the MP2 energy
!        Note: here we don't consider the t1 amplitude that is different with excorr_mp2no case. 
         ierr=initialize_tensor_t(result_tensor)
         ierr=t_einsum("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv) ! o(2)v(2)
         mp2_energy=result_tens(1)*ONE_QUARTER_R
         write(*,*) ""
         write(*,*) "    MP2 energy = ", mp2_energy
         write(*,*) ""

         call checkpoint_write ('/result/wavefunctions/mbpt/method',sdata="MP2")
         call checkpoint_write ('/result/wavefunctions/mbpt/e_corr',rdata=mp2_energy)
         call checkpoint_read ('/result/wavefunctions/scf/energy',rdata=scf_energy)
         call checkpoint_write ('/result/wavefunctions/mbpt/energy',rdata=scf_energy+mp2_energy)
         call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=int_t%oovv)

!        Calculate and print the MP2 density matrix of virtual-virtual block Dvv 
         ierr=initialize_tensor_t(MP2_vv_density)
         if(exa_input%mp2no.eq.3) then
            call eti_generate_FNOpp_density_matrix(exa_input, MP2_vv_density, comm_t%t2, int_t)
            call print_date('FNO++: Finished MP2_DM_VV calculation and save it to file DM_complex!')
         else
            ierr=t_einsum("D(a,b)+=V+(i,j,c,a)*V(i,j,c,b)",MP2_vv_density,&
                                        int_t%oovv,int_t%oovv,scale=ONE_HALF) ! o(2)v(2)
         endif

!        Make some noise so that we know we are proceeding
         call print_date('finished MP2_DM_VV calculation and save it to file DM_complex')

         if (exa_input%do_neo) then
            call generate_effec_densive(exa_input,MP2_vv_density,CIS_vv_density)
         end if

         call write_2dtensor_as_matrix(nvir,MP2_vv_density,'DM_complex')

         call eti_executor_shutdown

        end subroutine eti_generate_MP2_VNO_MObasis

!------------------------------------------------------------------------------------------------------------------------------------
      subroutine get_transform_Fock (exa_input)
         ! This routine is used for doing transform of fock matrix from MO basis to NO basis.
         ! The new basis set comes from the digonalization of density matrix.
         ! The all natural orbitals and occupation numbers are saved in file NOs_MO
         ! The transformed Fock matrix are saved in file FM_in_NatOrb
         ! In principle, this routine is general. Any type complex matrix could be used for getting the new basis and transformed Fock matrix.
         ! This routine doesn't truncate the new basis set.
         ! Written by Xiang Yuan

                  type(exacc_input), intent(in)      :: exa_input

                  type(cmo)                          :: qorbital
                  type(basis_set_info_t)             :: ao_basis
                  real(kind=8), pointer              :: N_OCCUPATION(:)         ! occupation number from qdia of DM
                  real(kind=8), pointer              :: N_OCCUPATION_tmp(:)     ! occupation number from qdia of DM
                  complex(kind=8), pointer           :: MP2_density_comp(:,:)   ! MP2 density matrix from DM_complex
                  real(8), pointer                   :: MP2_density_quat(:,:,:) ! quaternion MP2 density matrix
                  real(kind=8), pointer              :: NOs_MO(:,:,:)           ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: NOs_MO_tmp(:,:,:)       ! natural orbital from qdia of DM
                  real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix
                  real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
                  !integer                            :: sz,flb,sz_occ           ! size of matrix and flag for block
                  integer                            :: sz1, sz2                ! size of matrix and flag for block
                  real(kind=8)                       :: Sum_N_OCCUPATION_total, Sum_N_OCCUPATION_total_tmp
                  integer                            :: i, ierr, nocc, nvirt
                  integer                            :: i_mp2fno_prop, counter
                  integer                            :: N_reduced, N_reduced_tmp
                  character(len=20)                  :: fn, istr

                  if((exa_input%mp2no==3).and..not.exa_input%mp2fno_weighted_sum) then
                     ! First loop to generate the FNO and save them in NOs_MO files and get the size of the full NOs_NO space
                     sz1 = 0
                     sz2 = 0
                     do i_mp2fno_prop = 1,exa_input%mp2fno_n_prop
                        call print_date('Diagonalizing perturbed density...')
                        call read_matrix_comp(MP2_density_comp,&
                            trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))
                        call print_date('read density matrix (DM) in MO basis (spinor) from DM_complex')

                        call complex_to_quaternion(MP2_density_comp,MP2_density_quat)
                        call print_date('changed MO DM representation (spinor->quaternion)')

                        call diagonal_matrix_quater(MP2_density_quat,N_OCCUPATION_tmp,NOs_MO_tmp)
                        call print_date('diagonalized MO DM (quaternion)')

                        ! sort NOs_MO so that orbitals with higher occupation are first
                        call sort_MOs_by_occupation_number_v1 (N_OCCUPATION_tmp,NOs_MO_tmp)

                        call quater_write(NOs_MO_tmp,N_OCCUPATION_tmp,&
                            trim(get_filename_with_index('NOs_MO',i_mp2fno_prop)))
                        call print_date('wrote all natural orbital and occupation number in file NOs_MO')

                        ! Count orbitals in the final file
                        sz1 = size(NOs_MO_tmp,1)
                        do i = 1,size(N_OCCUPATION_tmp)
                           if(N_OCCUPATION_tmp(i).gt.exa_input%t_mp2no) then
                              sz2 = sz2 + 1
                           endif
                        enddo

                        deallocate(MP2_density_quat) ! allocation done in complex_to_quaternion
                        deallocate(NOs_MO_tmp)       ! allocation done in diagonal_matrix_quater
                        deallocate(N_OCCUPATION_tmp) ! allocation done in diagonal_matrix_quater
                     enddo

                     if(sz2.gt.sz1) then
                        print*, "------------------"
                        print*, "- MP2FNO++ ERROR -"
                        print*, "------------------"
                        print*, "Too many natural orbitals. Please reduce the MP2 occupation number threshold."
                        print*, "Maximum allowed number of orbitals = ",sz1
                        print*, "Number of generated nataural orbitals = ",sz2
                        print*, "Occupation number threshold = ",exa_input%t_mp2no
                        call quit ('Error in eti_mp2no.F90')
                     endif

                     allocate(NOs_MO(sz1,sz2,4))
                     allocate(N_OCCUPATION(sz2))

                     Sum_N_OCCUPATION_total = 0.0
                     N_reduced = 0
                     counter = 0
                     ! Read the FNO from files and concatenate them in arrays 
                     do i_mp2fno_prop = 1,exa_input%mp2fno_n_prop
                        call quater_read_orbital(NOs_MO_tmp,N_OCCUPATION_tmp,Sum_N_OCCUPATION_total_tmp,exa_input%t_mp2no,&
                            N_reduced_tmp,trim(get_filename_with_index('NOs_MO',i_mp2fno_prop)))

                        forall (i=1:size(N_OCCUPATION_tmp)) N_OCCUPATION(i+counter) = N_OCCUPATION_tmp(i)
                        forall (i=1:size(N_OCCUPATION_tmp)) NOs_MO(:,i+counter,:) = NOs_MO_tmp(:,i,:)

                        Sum_N_OCCUPATION_total = Sum_N_OCCUPATION_total + Sum_N_OCCUPATION_total_tmp
                        N_reduced = N_reduced + N_reduced_tmp
                        counter = counter + size(N_OCCUPATION_tmp)

                        deallocate(NOs_MO_tmp)       ! allocation done in quater_read_orbital
                        deallocate(N_OCCUPATION_tmp) ! allocation done in quater_read_orbital
                     enddo

                     ! sort the entire set of natural orbitals generated with all perturbations
                     call sort_MOs_by_occupation_number_v2(N_OCCUPATION,NOs_MO,size(NOs_MO,1),size(NOs_MO,2),4)
                     call print_mp2fvno_occupation_number_details(exa_input%t_mp2no, N_reduced, N_OCCUPATION,& 
                        Sum_N_OCCUPATION_total)
                  else
                     call read_matrix_comp(MP2_density_comp,'DM_complex')
                     call print_date('read density matrix (DM) in MO basis (spinor) from DM_complex')

                     call complex_to_quaternion(MP2_density_comp,MP2_density_quat)
                     call print_date('changed MO DM representation (spinor->quaternion)')

                     call diagonal_matrix_quater(MP2_density_quat,N_OCCUPATION,NOs_MO)
                     call print_date('diagonalized MO DM (quaternion)')

                     call sort_MOs_by_occupation_number_v1 (N_OCCUPATION,NOs_MO)
                     call print_date('Sorted natural orbitals by occupation number')

                  endif

                  if(exa_input%print_occupation_number.eq.1) then
                      print*, "-------------------------------"
                      print*, "Printing all occupation numbers"
                      print*, "-------------------------------"
                      do i=1,size(N_OCCUPATION)
                         write(*,'(I6, A2, E15.4)') i,') ', N_OCCUPATION(i)
                      enddo
                  endif

                  call quater_write(NOs_MO,N_OCCUPATION,'NOs_MO')
                  call print_date('wrote all natural orbital and occupation number in file NOs_MO')

                  ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
                  call read_from_dirac(ao_basis,qorbital,'scf',ierr)
                  call convert_mo_to_quaternion(qorbital)
                  nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
                  nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
                  allocate(Fock_HF(nvirt,nvirt,4))
                  Fock_HF = 0.0
                  forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
                  call print_date('read Fock matrix from checkpoint')

                  call basis_transform_matrix(Fock_HF,NOs_MO,Fock_NO)
                  call print_date('got full VV Fock matrix in natural orbital basis F_VV_NO')

                  call quater_write(Fock_NO,'FM_inNatOrb')
                  call print_date('wrote F_VV_NO in file FM_in_NatOrb')

      end subroutine get_transform_Fock


!------------------------------------------------------------------------------------------------------------------------------
      subroutine get_new_orbital(exa_input,threshold,NOs_AO,N_OCCUPATION_reduced,New_orbital_energy,file_orb,N_reduced)
      ! This routine is used for getting new truncated canonical orbital in AO basis
      ! and the corresponding orbital energy.


            type(exacc_input), intent(in)      :: exa_input
            real(kind=8),intent(in)            :: threshold
            real(kind=8), pointer,intent(out)  :: NOs_AO(:,:,:)           ! new canonical orbital in AO basis
            real(kind=8), pointer,intent(out)  :: N_OCCUPATION_reduced(:) ! occupation number from qdia of DM
            real(kind=8), pointer,intent(out)  :: New_orbital_energy(:)   ! new orbital energy from qdia of Fock matrix in new orbital basis
            character(len=*),intent(in)        :: file_orb
            integer,intent(inout)              :: N_reduced               ! number of reduced natural orbital

            type(cmo)                          :: qorbital
            type(basis_set_info_t)             :: ao_basis
            real(kind=8)                       :: Sum_N_OCCUPATION_total
            real(kind=8)                       :: Sum_N_OCCUPATION_reduced
            real(kind=8)                       :: Perc_OCCUPATION_reduced
            real(kind=8), pointer              :: NOs_MO_reduced(:,:,:)   ! natural orbitals which have occupation number larger than threshold
            real(kind=8), pointer              :: New_cano_MO(:,:,:)      ! new canonical orbital in AO basis
            real(kind=8), pointer              :: Fock_NO(:,:,:)          ! v-v block of Fock matrix in natural orbital basis
            real(kind=8), pointer              :: Fock_HF(:,:,:)          ! v-v block of Fock matrix
            real(kind=8), pointer              :: HF_AO(:,:,:)            ! Hartree-Fock orbital
            integer                            :: nocc, nvirt             ! number of occupied and virtuals
            integer                            :: ierr,i,j,k

            call quater_read(NOs_MO_reduced,N_OCCUPATION_reduced,Sum_N_OCCUPATION_total,&
                              threshold,N_reduced,file_orb)
            call print_date('got the truncated NO and occupation number from NOs_MO')

            if((exa_input%mp2no.ne.3).or.exa_input%mp2fno_weighted_sum) then
               call print_mp2fvno_occupation_number_details(threshold, N_reduced, N_OCCUPATION_reduced, Sum_N_OCCUPATION_total)
            endif

            nvirt = exa_input%nkr_vir        ! the number of virtual orbitals
            nocc  = exa_input%mokr_occ(size(exa_input%mokr_occ))
            call read_from_dirac(ao_basis,qorbital,'scf',ierr)
            call convert_mo_to_quaternion(qorbital)
            allocate(HF_AO(qorbital%nao,nvirt,4))
            HF_AO = qorbital%coeff_q(:,nocc+1:nocc+nvirt,:)
            call print_date('read HF orbital coefficient from checkpoint')

            ! we initialize the (diagonal) vv block of the Fock matrix with the orbital energies
            allocate(Fock_HF(nvirt,nvirt,4))
            Fock_HF = 0.0
            forall (i=1:nvirt) Fock_HF(i,i,1) = qorbital%energy(nocc+i)
            call print_date('read Fock matrix from checkpoint')

            call basis_transform_matrix(Fock_HF,NOs_MO_reduced,Fock_NO)
            call print_date('transformed Fock matrix into reduced NO basis')

            call get_Cano_Orbital(NOs_MO_reduced,New_cano_MO,Fock_NO,New_orbital_energy)
            call print_date('got new Canonical orbital in MO basis')

            call basis_transform_vector(New_cano_MO,HF_AO,NOs_AO)
            call print_date('got new orbitals in AO basis')

            call print_rMP2FVNO(New_orbital_energy,N_OCCUPATION_reduced)
      end subroutine

!------------------------------------------------------------------------------------------------------------------------------
      subroutine generate_effec_densive(exa_input,MP2_vv_density,CIS_vv_effden)

            type(exacc_input), intent(in) :: exa_input
            type(tensor_t)          :: MP2_vv_density
            type(tensor_t)          :: CIS_vv_effden
            integer(INTD)               :: ierr

            call print_date('Start computing CIS effective density matrix')
            call eti_cis_driver(exa_input,CIS_vv_effden,initialize_tensor_library=.false.,initialize_interest=.false.)
            ierr=t_einsum('D(a,b)+=DC(a,b)',MP2_vv_density,CIS_vv_effden,scale=ONE)
            call print_date('Add CIS effective density matrix into MP2 desnity matrix')

      end subroutine

      subroutine get_prop_bar_ov(prop_bar_ov,prop_ov,t2)
         type(tensor_t), intent(inout)       :: prop_bar_ov, prop_ov, t2
         integer(INTD)                           :: ierr

         ierr=t_einsum("A(i,a)+=B(i,a)",prop_bar_ov,prop_ov)
         ierr=t_einsum("A(i,a)+=T(a,b,i,j)*B(j,b)",prop_bar_ov,t2,prop_ov)
      end subroutine

      subroutine get_prop_bar_oovv(prop_bar_oovv,orb_int_tensor1,t2)
         type(tensor_t), intent(inout)               :: prop_bar_oovv, t2
         type(eti_orb_prop_integ_tens_t), intent(inout):: orb_int_tensor1
         integer(INTD)                                   :: ierr

         ! First term with P(a,c)
         ierr=t_einsum("A(i,j,a,c)+=T(e,c,i,j)*B(a,e)",prop_bar_oovv,t2,orb_int_tensor1%prop_vv)
         ierr=t_einsum("A(i,j,a,c)+=T(e,a,i,j)*B(c,e)",prop_bar_oovv,t2,orb_int_tensor1%prop_vv,scale=MINUS_ONE)

         ! Second term P(i,j)
         ierr=t_einsum("A(i,j,a,c)+=T(a,c,m,j)*B(m,i)",prop_bar_oovv,t2,orb_int_tensor1%prop_oo,scale=MINUS_ONE)
         ierr=t_einsum("A(i,j,a,c)+=T(a,c,m,i)*B(m,j)",prop_bar_oovv,t2,orb_int_tensor1%prop_oo)
      end subroutine

      subroutine eti_generate_FNOpp_density_matrix(exa_input, MP2_vv_density, t2_MP2, int_t)
         type(exacc_input), intent(in) :: exa_input
         type(tensor_t), intent(inout)      :: MP2_vv_density
         type(eti_intg_tens_t), intent(inout) :: int_t
         type(tensor_t), intent(inout) :: t2_MP2

         type(eti_comm_tens_t) :: comm_t
         type(tensor_t) :: MP2_vv_density_tmp
         type(eti_orb_prop_integ_tens_t)   :: orb_int_tensor1
         type(eti_rsp_interm_tens_t)       :: interm_t_prop1

         complex(8), pointer   :: hbar_oo_ptr(:,:), hbar_vv_ptr(:,:)
         complex(8), pointer   :: t1_ptr(:,:), t2_ptr(:,:,:,:)
         complex(8), pointer   :: MP2_vv_density_array_ptr(:,:)
         complex(8), pointer   :: MP2_vv_density_array_tmp(:,:)
         type(C_PTR):: body_p
         integer(INTD), dimension(4) :: vvoo_dims, oovv_dims
         integer(INTD), dimension(2) :: vv_dims, ov_dims, oo_dims, vo_dims
         integer                                :: nocc, nvir              ! the size of the mo basis for occupied and virtual spinors
         integer :: a, c, i, i_mp2fno_prop, j
         integer(INTD)                                   :: ierr

         nocc = exa_input%nocc
         nvir = exa_input%nvir

         call print_date("Computing MP2FNO++ density matrix with config:")
         write(*,*) "FNO++ omega = ", exa_input%mp2fno_prop_omega
         ! May want to print the perturbation name instead of its PRPTRA index
         write(*,*) "#perturbations = ", exa_input%mp2fno_n_prop
         write(*,*) "perturbations = ", exa_input%mp2fno_prop
         write(*,*) "Weighted sum of the perturbations = ", exa_input%mp2fno_weighted_sum
         if(exa_input%mp2fno_weighted_sum) then
            write(*,*) "weights = ", exa_input%mp2fno_prop_weights
         endif

         ! FNO++
         ! Notations from:
         ! Reduced-scaling coupled cluster response theory: Challenges and opportunities (2018)

         ! dimensions
         vv_dims(:) = nvir

         ov_dims(1) = nocc
         ov_dims(2) = nvir

         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir

         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc

         vo_dims(1) = nvir
         vo_dims(2) = nocc

         oo_dims = nocc

         ! Hbar in comm_t%hoo and comm_t%hvv (not to confuse with actual hoo and hvv)
         ierr=construct_tensor_t(comm_t%hoo,exa_input%eti_precision%storage_t,oo_dims,init_val=ZERO)
         ierr=t_einsum("H(i,j)+=F(i,j)",comm_t%hoo,int_t%foo)
         ierr=t_einsum("H(i,j)+=T(e,f,i,n)*V(j,n,e,f)",comm_t%hoo,t2_MP2,int_t%oovv,scale=ONE_HALF)

         ierr=construct_tensor_t(comm_t%hvv,exa_input%eti_precision%storage_t,vv_dims,init_val=ZERO)
         ierr=t_einsum("H(a,b)+=F(a,b)",comm_t%hvv,int_t%fvv)
         ierr=t_einsum("H(a,b)+=T(f,a,m,n)*V(m,n,f,b)",comm_t%hvv,t2_MP2,int_t%oovv,scale=MINUS_ONE_HALF)

         ierr=construct_tensor_t(comm_t%t1,exa_input%eti_precision%storage_t,vo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%t2,exa_input%eti_precision%storage_t,vvoo_dims,init_val=ZERO)

         ! pointer access
         ierr=direct_data_pointer(comm_t%t1,body_p)
         call c_f_pointer(body_p,t1_ptr,vo_dims)
         ierr=direct_data_pointer(comm_t%hoo,body_p)
         call c_f_pointer(body_p,hbar_oo_ptr,oo_dims)
         ierr=direct_data_pointer(comm_t%hvv,body_p)
         call c_f_pointer(body_p,hbar_vv_ptr,vv_dims)
         ierr=direct_data_pointer(comm_t%t2,body_p)
         call c_f_pointer(body_p,t2_ptr,vvoo_dims)

         !---------------
         ! Generate density matrix for each perturbation
         !---------------
         do i_mp2fno_prop=1,exa_input%mp2fno_n_prop

            write(*,*) "Generating density matrix of perturbation ", exa_input%mp2fno_prop(i_mp2fno_prop)

            ! Abar ai
            call get_prop_integral(exa_input,orb_int_tensor1,exa_input%mp2fno_prop(i_mp2fno_prop))
            interm_t_prop1%property_label = orb_int_tensor1%property_label

            ierr=construct_tensor_t(interm_t_prop1%prop_bar_ov,exa_input%eti_precision%storage_t,ov_dims,init_val=ZERO)
            call get_prop_bar_ov(interm_t_prop1%prop_bar_ov,orb_int_tensor1%prop_ov,t2_MP2)

            ! Abar acij
            ierr=construct_tensor_t(interm_t_prop1%prop_bar_oovv,exa_input%eti_precision%storage_t,oovv_dims,init_val=ZERO)
            call get_prop_bar_oovv(interm_t_prop1%prop_bar_oovv,orb_int_tensor1,t2_MP2)

            ! t1
            ierr=initialize_tensor_t(comm_t%t1)
            ierr=t_einsum("T(a,i)+=A(i,a)",comm_t%t1,interm_t_prop1%prop_bar_ov)

            do i=1,nocc
               do a=1,nvir
                   t1_ptr(a,i) = t1_ptr(a,i)/(hbar_vv_ptr(a,a) - hbar_oo_ptr(i,i) + exa_input%mp2fno_prop_omega)
               enddo
            enddo


            ! t2
            ierr=initialize_tensor_t(comm_t%t2)
            ierr=t_einsum("T(a,c,i,j)+=A(i,j,a,c)",comm_t%t2,interm_t_prop1%prop_bar_oovv)

            do i=1,nocc
               do j=1,nocc
                  do a=1,nvir
                     do c=1,nvir
                        t2_ptr(a,c,i,j) = t2_ptr(a,c,i,j)/(hbar_vv_ptr(a,a) + hbar_vv_ptr(c,c)&
                           - hbar_oo_ptr(i,i) - hbar_oo_ptr(j,j) + exa_input%mp2fno_prop_omega)
                     enddo
                  enddo
               enddo
            enddo

            ! MP2 density matrix
            ierr=initialize_tensor_t(MP2_vv_density)
            ierr=t_einsum("D(a,b)+=T(a,c,i,j)*T(b,c,i,j)",MP2_vv_density,&
                                        comm_t%t2, comm_t%t2,scale=ONE_HALF)
            ierr=t_einsum("D(a,b)+=T(a,i)*T(b,i)",MP2_vv_density,&
                                        comm_t%t1, comm_t%t1)

            ! Generate one file per perturbation -> use of HDF5?
            write(*,*) "Saving density matrix in  ", trim(get_filename_with_index('DM_complex',i_mp2fno_prop))
            call write_2dtensor_as_matrix(nvir,MP2_vv_density,&
                trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))

            ! destroy tensors
            ierr=destroy_tensor_t(orb_int_tensor1%prop_oo)
            ierr=destroy_tensor_t(orb_int_tensor1%prop_vv)
            ierr=destroy_tensor_t(orb_int_tensor1%prop_vo)
            ierr=destroy_tensor_t(orb_int_tensor1%prop_ov)
            ierr=destroy_tensor_t(interm_t_prop1%prop_bar_ov)
            ierr=destroy_tensor_t(interm_t_prop1%prop_bar_oovv)

         enddo

         ! Generate weighted density matrix
         if(exa_input%mp2fno_weighted_sum) then
            ! initialize arrays
            ierr=initialize_tensor_t(MP2_vv_density)
            ierr=construct_tensor_t(MP2_vv_density_tmp,exa_input%eti_precision%storage_t,vv_dims,init_val=ZERO)

            do i_mp2fno_prop=1,exa_input%mp2fno_n_prop
               ! read the matrix
               call read_matrix_comp(MP2_vv_density_array_tmp,&
                   trim(get_filename_with_index('DM_complex',i_mp2fno_prop)))

               ! assign the values of an array to a tensor
               ierr=direct_data_pointer(MP2_vv_density_tmp,body_p)
               call c_f_pointer(body_p,MP2_vv_density_array_ptr,vv_dims)
               do a=1,nvir
                  do c=1,nvir
                     MP2_vv_density_array_ptr(a,c) = MP2_vv_density_array_tmp(a,c)
                  enddo
               enddo

               ! contract
               ierr=t_einsum("D(a,b)+=D_tmp(a,b)",MP2_vv_density,MP2_vv_density_tmp,&
                   scale=dcmplx(exa_input%mp2fno_prop_weights(i_mp2fno_prop),0.0))

               ! deallocate/destroy arrays
               deallocate(MP2_vv_density_array_tmp)
            enddo
            ierr=direct_data_pointer(MP2_vv_density,body_p)

            ! destroy arrays
            ierr=destroy_tensor_t(MP2_vv_density_tmp)
         endif

      endsubroutine eti_generate_FNOpp_density_matrix

end module eti_mp2no

