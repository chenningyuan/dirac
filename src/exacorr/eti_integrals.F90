module eti_integrals
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_cfg, only: exacc_input
        use exacorr_utils
        use eti_common_routines
        use exacorr_printout
        use exacorr_scaling_parameters
        use, intrinsic:: ISO_C_BINDING

        implicit none

        private

        public get_CC_integrals

       contains

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t,moint_scheme,print_level)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list.

!        Written by Lucas Visscher, summer 2017

         use exacorr_global
         use exacorr_ao
         use eti_ao_to_mo

         integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(eti_intg_tens_t),intent(inout)  :: int_t
         integer, intent(in)                  :: print_level
         integer, intent(in)                  :: moint_scheme
         integer(INTD) :: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4),ooov_dims(1:4),vovv_dims(1:4)

         integer(INTD)                      :: ierr
         complex(8), pointer, contiguous    :: ao_tens(:,:,:,:)
         type(C_PTR)                        :: body_p

         integer            :: nao, nshells  ! number of AOs and number of AO shells
         type(basis_func_info_t), allocatable :: gto(:)   ! array with basis function information
         integer :: first_ao(4), last_ao(4)      ! first and last ao for each index
         integer :: basis_angular, koff, loff, ksh, lsh, lk, ll, nk, nl
      
!        auxilliary integral tensors (only needed inside this routine)
         type(tensor_t) :: aoint_tensor, a2int_tensor, ovvo_tensor
         type(tensor_t) :: oooo_tensor, ooov_tensor, ovov_tensor, vvoo_tensor, vvov_tensor, vvvv_tensor
         integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4), a2int_dims(1:4)
         
!        arrays needed to communicate with the MO integral generator
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)
         type(ao2mo_info) :: occ_space, vir_space, aomo_spaces(4)

!        scalars (need to be defined as tensor types)
         type(tensor_t) :: one_tensor, minusone_tensor
         integer(INTD)      :: minusone_dims(1)

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Store information about the occupied and virtual spaces that we transform (better pass this type already)
         occ_space%nao = nao
         occ_space%nmo = nocc
         occ_space%first_ao = 1
         occ_space%last_ao = nao
         allocate(occ_space%mo_list(nocc))
         occ_space%mo_list = mo_occ

         vir_space%nao = nao
         vir_space%nmo = nvir
         vir_space%first_ao = 1
         vir_space%last_ao = nao
         allocate(vir_space%mo_list(nvir))
         vir_space%mo_list = mo_vir

         if (moint_scheme<3) then ! schemes 1 and 2 take all AO integrals in memory
            ! Get tensor with AO integrals
            aoint_dims      = nao
            ierr=construct_tensor_t(aoint_tensor,C8,aoint_dims,init_val=ZERO)
            if (ierr /= 0) then
               print*, "Could not construct ao integral tensor, likely not enough memory available"
               call quit ('error constructing aoint_tensor')
            end if
            call print_date('allocated AO integral tensor')
!           Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
            ierr=direct_data_pointer(aoint_tensor,body_p)
            call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
            call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
            call print_date('finished ao integral calculation')

!           Get antisymmetric oovv tensor
!           Note that the dimensions for the tensor are given in physicists notation <12||34>
            oovv_dims(1:2) = nocc
            oovv_dims(3:4) = nvir
            ierr=construct_tensor_t(int_t%oovv,C8,oovv_dims,init_val=ZERO)
!           whereas the spaces for the index transform are in mullikens notation (13||24)
            aomo_spaces(1) = occ_space
            aomo_spaces(2) = vir_space
            aomo_spaces(3) = occ_space
            aomo_spaces(4) = vir_space
            call get_integral_tensor (aoint_tensor,int_t%oovv,12,aomo_spaces)
            call print_date('finished oovv integrals')

            oooo_dims = nocc
            ierr=construct_tensor_t(int_t%oooo,C8,oooo_dims,init_val=ZERO)
            aomo_spaces(1) = occ_space
            aomo_spaces(2) = occ_space
            aomo_spaces(3) = occ_space
            aomo_spaces(4) = occ_space
            call get_integral_tensor (aoint_tensor,int_t%oooo,12,aomo_spaces)
            call print_date('finished oooo integrals')

!           Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!           Get vovo tensor
            vovo_dims(1) = nvir
            vovo_dims(2) = nocc
            vovo_dims(3) = nvir
            vovo_dims(4) = nocc
            ierr=construct_tensor_t(int_t%vovo,C8,vovo_dims,init_val=ZERO)
            aomo_spaces(1) = vir_space
            aomo_spaces(2) = vir_space
            aomo_spaces(3) = occ_space
            aomo_spaces(4) = occ_space
            call get_integral_tensor (aoint_tensor,int_t%vovo,0,aomo_spaces)
!           Get ovvo tensor
            ovvo_dims(1) = nocc
            ovvo_dims(2) = nvir
            ovvo_dims(3) = nvir
            ovvo_dims(4) = nocc
            ierr=construct_tensor_t(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
            aomo_spaces(1) = occ_space
            aomo_spaces(2) = vir_space
            aomo_spaces(3) = vir_space
            aomo_spaces(4) = occ_space
            call get_integral_tensor (aoint_tensor,ovvo_tensor,0,aomo_spaces)
!           Make anti-symmetrized vovo tensor
            minusone_dims(1) = 1
            ierr=construct_tensor_t(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
            ierr=t_einsum("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,minusone_tensor)
            ierr=destroy_tensor_t(ovvo_tensor)
            call print_date('finished vovo integrals')

!           Get vvvv tensor
            vvvv_dims = nvir
            ierr=construct_tensor_t(int_t%vvvv,C8,vvvv_dims,init_val=ZERO)
            aomo_spaces(1) = vir_space
            aomo_spaces(2) = vir_space
            aomo_spaces(3) = vir_space
            aomo_spaces(4) = vir_space
            call get_integral_tensor (aoint_tensor,int_t%vvvv,12,aomo_spaces)
            call print_date('finished vvvv integrals')

!           Get ooov tensor
            ooov_dims(1:3) = nocc
            ooov_dims(4)   = nvir
            ierr=construct_tensor_t(int_t%ooov,C8,ooov_dims,init_val=ZERO)
            aomo_spaces(1) = occ_space
            aomo_spaces(2) = occ_space
            aomo_spaces(3) = occ_space
            aomo_spaces(4) = vir_space
            call get_integral_tensor (aoint_tensor,int_t%ooov,12,aomo_spaces)
            call print_date('finished ooov integrals')

!           Get vovv tensor
            vovv_dims(1)   = nvir
            vovv_dims(2)   = nocc
            vovv_dims(3:4) = nvir
            ierr=construct_tensor_t(int_t%vovv,C8,vovv_dims,init_val=ZERO)
            aomo_spaces(1) = vir_space
            aomo_spaces(2) = vir_space
            aomo_spaces(3) = occ_space
            aomo_spaces(4) = vir_space
            call get_integral_tensor (aoint_tensor,int_t%vovv,34,aomo_spaces)
            call print_date('finished vovv integrals')
            ! Clean-up
            ierr=destroy_tensor_t(aoint_tensor)

         else if (moint_scheme==4) then

           ! allocate the seven auxilliary tensors (Mulliken order) into which contributions get accumulated
           ierr=construct_tensor_t(oooo_tensor,C8,[nocc,nocc,nocc,nocc],init_val=ZERO)
           ierr=construct_tensor_t(ooov_tensor,C8,[nocc,nocc,nocc,nvir],init_val=ZERO)
           ierr=construct_tensor_t(ovov_tensor,C8,[nocc,nvir,nocc,nvir],init_val=ZERO)
           ierr=construct_tensor_t(ovvo_tensor,C8,[nocc,nvir,nvir,nocc],init_val=ZERO)
           ierr=construct_tensor_t(vvoo_tensor,C8,[nvir,nvir,nocc,nocc],init_val=ZERO)
           ierr=construct_tensor_t(vvov_tensor,C8,[nvir,nvir,nocc,nvir],init_val=ZERO)
           ierr=construct_tensor_t(vvvv_tensor,C8,[nvir,nvir,nvir,nvir],init_val=ZERO)

           allocate (mo_list(2*max(nocc,nvir))) ! allocate the mo list with the maximum size

           call print_date('starting MO-integral generation')
           ! in this algorithm we only keep part of the AO and HT integrals in memory
           basis_angular = get_basis_angular()
           call get_gtos(1,nao,gto,nshells)
           first_ao(1:2) = 1
           last_ao(1:2) = nao
           aoint_dims(1:2) = nao
           loff = 1
           do lsh = 1, nshells
               ll   =  gto(lsh)%orb_momentum
               nl   =  nfunctions(ll,basis_angular)
               first_ao(4)   = loff
               last_ao(4)    = loff + nl - 1
               aoint_dims(4) = nl
               koff = 1
               do ksh = 1, nshells
                  lk   =  gto(ksh)%orb_momentum
                  nk   =  nfunctions(lk,basis_angular)
                  ! adjust ranges of AOs that can be transformed
                  first_ao(3)   = koff
                  last_ao(3)    = koff + nk - 1
                  aoint_dims(3) = nk
                  ierr=construct_tensor_t(aoint_tensor,C8,aoint_dims,init_val=ZERO)
                  if (ierr /= 0) then
                     print*, "Could not construct ao integral tensor"
                     call quit ('error constructing aoint_tensor')
                  end if
                  ierr=direct_data_pointer(aoint_tensor,body_p)
                  call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
     
                  call compute_ao_integrals (first_ao,aoint_dims,ao_tens_complex=ao_tens)

                  ! Setting dimensions for the half transformed integral tensor, the last two are the untransformed AO spaces
                  a2int_dims(3:4)   = aoint_dims(3:4)

                  ! get the integrals with vv as first transform
                  a2int_dims(1:2)   = nvir
                  ierr=construct_tensor_t(a2int_tensor,C8,a2int_dims,init_val=ZERO)
                  if(ierr /= 0) call quit('error constructing a2int_tensor for vv')
                  aomo_spaces(1) = vir_space
                  aomo_spaces(2) = vir_space
                  call do_ht1 (aoint_tensor,a2int_tensor,aomo_spaces(1:2))
                  
                  ! vvvv transformation
                  aomo_spaces(3) = vir_space
                  aomo_spaces(4) = vir_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,vvvv_tensor,aomo_spaces(3:4))

                  ! vvov transformation
                  aomo_spaces(3) = occ_space
                  aomo_spaces(4) = vir_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,vvov_tensor,aomo_spaces(3:4))

                  ! vvoo transformation
                  aomo_spaces(3) = occ_space
                  aomo_spaces(4) = occ_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,vvoo_tensor,aomo_spaces(3:4))
                  
                  ierr=destroy_tensor_t(a2int_tensor)

                  ! get the integrals with ov as first transform
                  a2int_dims(1)   = nocc
                  a2int_dims(2)   = nvir
                  ierr=construct_tensor_t(a2int_tensor,C8,a2int_dims,init_val=ZERO)
                  if(ierr /= 0) call quit('error constructing a2int_tensor for ov')
                  aomo_spaces(1) = occ_space
                  aomo_spaces(2) = vir_space
                  call do_ht1 (aoint_tensor,a2int_tensor,aomo_spaces(1:2))

                  ! ovov transformation
                  aomo_spaces(3) = occ_space
                  aomo_spaces(4) = vir_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,ovov_tensor,aomo_spaces(3:4))
                  
                  ! ovvo transformation
                  aomo_spaces(3) = vir_space
                  aomo_spaces(4) = occ_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,ovvo_tensor,aomo_spaces(3:4))                

                  ierr=destroy_tensor_t(a2int_tensor)

                  ! get the integrals with oo as first transform
                  a2int_dims(1:2) = nocc
                  ierr=construct_tensor_t(a2int_tensor,C8,a2int_dims,init_val=ZERO)
                  if(ierr /= 0) call quit('error constructing a2int_tensor for oo')
                  aomo_spaces(1) = occ_space
                  aomo_spaces(2) = occ_space              
                  call do_ht1 (aoint_tensor,a2int_tensor,aomo_spaces(1:2))

                  ! oooo transformation
                  aomo_spaces(3) = occ_space
                  aomo_spaces(4) = occ_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)             
                  call do_ht2 (a2int_tensor,oooo_tensor,aomo_spaces(3:4))

                  ! ooov transformation
                  aomo_spaces(3) = occ_space
                  aomo_spaces(4) = vir_space
                  aomo_spaces(3:4)%first_ao = first_ao(3:4)
                  aomo_spaces(3:4)%last_ao  = last_ao(3:4)
                  call do_ht2 (a2int_tensor,ooov_tensor,aomo_spaces(3:4))
                  ierr=destroy_tensor_t(a2int_tensor)
                  
                  ierr=destroy_tensor_t(aoint_tensor)
                  koff = koff + nk
                  end do ! ksh
               loff = loff + nl
            end do ! lsh

            ! Sort and antisymmetrize the completed integrals

            ! The vovo tensor has contributions from two different Mulliken tensors, make this one first
            ierr=construct_tensor_t(int_t%vovo,C8,[nvir,nocc,nvir,nocc],init_val=ZERO)
            call mulliken_to_dirac_sort (vvoo_tensor,int_t%vovo,0)
            ierr=destroy_tensor_t(vvoo_tensor)
            minusone_dims(1) = 1
            ierr=construct_tensor_t(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
            ierr=t_einsum("V(a,i,b,j)+=V(i,b,a,j)*C()",int_t%vovo,ovvo_tensor,minusone_tensor)
            ! Note that this einsum is different than the one for scheme 1,2 because the tensor 
            ! is in Mulliken order here, while it is in Dirac order for scheme 1,2.
            ierr=destroy_tensor_t(ovvo_tensor)
            ierr=destroy_tensor_t(minusone_tensor)

            ierr=construct_tensor_t(int_t%oovv,C8,[nocc,nocc,nvir,nvir],init_val=ZERO)
            call mulliken_to_dirac_sort (ovov_tensor,int_t%oovv,12)
            ierr=destroy_tensor_t(ovov_tensor)

            ierr=construct_tensor_t(int_t%oooo,C8,[nocc,nocc,nocc,nocc],init_val=ZERO)
            call mulliken_to_dirac_sort (oooo_tensor,int_t%oooo,12)
            ierr=destroy_tensor_t(oooo_tensor)

            ierr=construct_tensor_t(int_t%ooov,C8,[nocc,nocc,nocc,nvir],init_val=ZERO)
            call mulliken_to_dirac_sort (ooov_tensor,int_t%ooov,12)
            ierr=destroy_tensor_t(ooov_tensor)

            ierr=construct_tensor_t(int_t%vvvv,C8,[nvir,nvir,nvir,nvir],init_val=ZERO)
            call mulliken_to_dirac_sort (vvvv_tensor,int_t%vvvv,12)
            ierr=destroy_tensor_t(vvvv_tensor)

            ierr=construct_tensor_t(int_t%vovv,C8,[nvir,nocc,nvir,nvir],init_val=ZERO)
            call mulliken_to_dirac_sort (vvov_tensor,int_t%vovv,34)
            ierr=destroy_tensor_t(vvov_tensor)

            deallocate (mo_list)

            call print_date('end of MO-integral generation')

         else
            call print_date('unknown index transformation algorithm specified')
            stop 'error in get_cc_integrals'
         end if ! branching between different schemes

       end subroutine get_CC_integrals

end module 
