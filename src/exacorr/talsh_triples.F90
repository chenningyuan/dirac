module eti_triples
!This module contains the routines needed to calculate perturbative triples
!written by Johann V. Pototschnig
!standard routine by L. Visscher, A. Papadopolus

        use tensor_algebra
        use exacorr_tensor_interface 
        use exacorr_datatypes
        use exacorr_utils
        use eti_common_routines
        use, intrinsic:: ISO_C_BINDING

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),TWO=(2.D0,0.D0), THREE=(3.0,0.0),   &
                                 SIX=(6.0,0.0), MINUS_ONE=(-1.D0,0.D0), MINUS_TWO=(-2.D0,0.D0),       &
                                 ONE_HALF=(0.5D0,0.D0), ONE_QUARTER=(0.25D0,0.D0),                    &
                                 ONE_TWELFTH=(1.D0,0.D0)/12, MINUS_ONE_HALF=(-0.5D0,0.D0),            &
                                 MINUS_ONE_QUARTER=(-0.25D0,0.D0), MINUS_ONE_TWELFTH=(-1.D0,0.D0)/12
        real(8), parameter    :: ONE_QUARTER_R=0.25D0, ONE_HALF_R=0.5D0

        private

        public t3_energy
        public t3_energy_laplace

        interface get_subtensor
           module procedure get_subtensor
        end interface

       contains

       subroutine t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor, &
                             ooov_tensor,vovv_tensor,oovv_tensor,print_level)

!        Evaluate +T, -T and (T) energy terms.
!        Definitions of terms: L. Visscher, T.J. Lee, K.G. Dyall, J Chem Phys. 105 (1996) 8769.
         
         type(tensor_t), intent(inout)  :: t1_tensor, t2_tensor
         type(tensor_t), intent(inout)  :: fov_tensor, ooov_tensor, vovv_tensor, oovv_tensor
         real(8), intent(in)                :: eps_occ(:),eps_vir(:)
         real(8), intent(out)               :: t_energy(3)
         integer, intent(in)                :: print_level

         integer(INTD)     :: t2_dims(4), abc_dims(3), ab_dims(2), ia_dims(2), abi_dims(3)
         type(tensor_t):: w_tensor, y_tensor, we_tensor
         type(tensor_t):: vovv_i_tensor, vovv_j_tensor, vovv_k_tensor
         type(tensor_t):: ooov_ij_tensor, ooov_jk_tensor, ooov_ki_tensor
         type(tensor_t):: t2_ij_tensor, t2_jk_tensor, t2_ki_tensor
         type(tensor_t):: t2_i_tensor, t2_j_tensor, t2_k_tensor

!        scalars (need to be defined as tensor types)
         type(C_PTR)         :: body_p
         type(tensor_t)  :: result_tensor
         integer(INTD)       :: result_dims(1)
         complex(8), pointer :: result_tens(:)

         real(8)        :: eps_ijk
         integer(INTD)  :: ierr, tens_rank
         integer        :: nocc, nvir, i, j, k
         complex(8), allocatable :: y(:,:,:) ! used in form_y_intermediate, but allocated here to have fewer alloc/dealloc calls

         integer(INTD)      :: one_dims(1)
         type(tensor_t) :: one_tripl
         !construct one tensor
         ierr=construct_tensor_t(one_tripl,C8,one_dims(1:0),init_val=ONE)

!        Get nocc and nvir from tensor information
         ierr = get_nmodes(t2_tensor,tens_rank)
         ierr = get_extents_static(t2_tensor,t2_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'error: t2 tensor corrupted'

         if (t2_dims(1).eq.t2_dims(2)) then
            nvir = t2_dims(1)
         else
            stop 'error: asymmetric tensor in t3_energy'
         end if

         if (t2_dims(3).eq.t2_dims(4)) then
            nocc = t2_dims(3)
         else
            stop 'error: asymmetric tensor in t3_energy'
         end if

!        Initialize scalars that are to be used as tensors in contractions
         result_dims(1) = 1
         ierr=construct_tensor_t(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)

!        Allocate tensors of dimension nvir**3 to hold intermediates
         abc_dims = nvir
         ierr=construct_tensor_t(w_tensor,C8,abc_dims,init_val=ZERO)
         ierr=construct_tensor_t(y_tensor,C8,abc_dims,init_val=ZERO)
         ierr=construct_tensor_t(we_tensor,C8,abc_dims,init_val=ZERO)
         allocate(y(nvir,nvir,nvir))

!        The following is a work-around, we need to implement contraction with restricted range for "inactive" indices
         ab_dims  = nvir
         ia_dims(1) = nocc
         ia_dims(2) = nvir
         abi_dims(1:2) = nvir
         abi_dims(3)   = nocc
         ierr=construct_tensor_t(vovv_i_tensor,C8,abc_dims,init_val=ZERO)
         ierr=construct_tensor_t(vovv_j_tensor,C8,abc_dims,init_val=ZERO)
         ierr=construct_tensor_t(vovv_k_tensor,C8,abc_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_jk_tensor,C8,ab_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_ki_tensor,C8,ab_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_ij_tensor,C8,ab_dims,init_val=ZERO)
         ierr=construct_tensor_t(ooov_ij_tensor,C8,ia_dims,init_val=ZERO)
         ierr=construct_tensor_t(ooov_jk_tensor,C8,ia_dims,init_val=ZERO)
         ierr=construct_tensor_t(ooov_ki_tensor,C8,ia_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_i_tensor,C8,abi_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_j_tensor,C8,abi_dims,init_val=ZERO)
         ierr=construct_tensor_t(t2_k_tensor,C8,abi_dims,init_val=ZERO)

!        Get access to result
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)

         t_energy = ZERO
         do k = 3, nocc
            call get_subtensor (vovv_tensor,vovv_k_tensor,432,k)
            call get_subtensor (t2_tensor,t2_k_tensor,433,k)
            do j = 2, k-1
               call get_subtensor (vovv_tensor,vovv_j_tensor,432,j)
               call get_subtensor (t2_tensor,t2_j_tensor,433,j)
               call get_subtensor (t2_tensor,t2_jk_tensor,4334,j,k)
               call get_subtensor (ooov_tensor,ooov_jk_tensor,4312,j,k)
               do i = 1, j-1
                  call get_subtensor (vovv_tensor,vovv_i_tensor,432,i)
                  call get_subtensor (t2_tensor,t2_i_tensor,433,i)
                  call get_subtensor (t2_tensor,t2_ij_tensor,4334,i,j)
                  call get_subtensor (t2_tensor,t2_ki_tensor,4334,k,i)
                  call get_subtensor (ooov_tensor,ooov_ij_tensor,4312,i,j)
                  call get_subtensor (ooov_tensor,ooov_ki_tensor,4312,k,i)

                  ! Form w-intermediate, first without symmetrization and minus sign in the auxilliary tensor array we
                  ierr=initialize_tensor_t(we_tensor)
                  ! For this term we use the permutation symmetry <vv||ov> = - <vo||vv>*, with vvov(abe:i) = vovv*(eab:i)
                  ierr=t_einsum("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_i_tensor,t2_jk_tensor) ! o(3)v(4)
                  ierr=t_einsum("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_j_tensor,t2_ki_tensor) ! o(3)v(4)
                  ierr=t_einsum("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_k_tensor,t2_ij_tensor) ! o(3)v(4)
                  ! For this term we use the permutation symmetry <vo||oo> = - <oo||ov>*, with vooo(al:ij) = ooov*(la:ij)
                  ierr=t_einsum("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_ij_tensor,t2_k_tensor) ! o(4)v(3)
                  ierr=t_einsum("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_jk_tensor,t2_i_tensor) ! o(4)v(3)
                  ierr=t_einsum("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_ki_tensor,t2_j_tensor) ! o(4)v(3)
                  if(print_level > 1) write(*,*) "W[sym] =",tensor_norm2(we_tensor)

                  ! Copy into w array to symmetrize and add the minus sign omitted in the above contractions
                  ierr=initialize_tensor_t(w_tensor)
                  ! symmetrization via contraction as addition with permutation is not implemented
                  ierr=t_einsum("W(a,b,c)+=W(a,b,c)",w_tensor,we_tensor,one_tripl,scale=MINUS_ONE)
                  ierr=t_einsum("W(a,b,c)+=W(b,c,a)",w_tensor,we_tensor,one_tripl,scale=MINUS_ONE)
                  ierr=t_einsum("W(a,b,c)+=W(c,a,b)",w_tensor,we_tensor,one_tripl,scale=MINUS_ONE)
                  if(print_level > 1) write(*,*) "W[+T ] =",tensor_norm2(w_tensor)

                  ! Form Y-intermediate and calculate 5th order -T energy contribution
                  call form_y3_intermediate (i,j,k,t1_tensor,t2_tensor,y_tensor,y)
                  if(print_level > 1) write(*,*) "Y[-T ] =",tensor_norm2(y_tensor)
                  ierr=initialize_tensor_t(result_tensor)
                  ierr=t_einsum("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(3)=t_energy(3) + result_tens(1)

                  ! Form Y-intermediate and calculate 5th order (T) energy contribution
                  call form_y2_intermediate (i,j,k,t1_tensor,t2_tensor,fov_tensor,oovv_tensor,y_tensor,y)
                  if(print_level > 1) write(*,*) "Y[(T)] =",tensor_norm2(y_tensor)
                  eps_ijk = eps_occ(i) + eps_occ(j) + eps_occ(k)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  ierr=initialize_tensor_t(result_tensor)
                  ierr=t_einsum("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(2)=t_energy(2) + result_tens(1)

                  ! Copy w into y array and calculate 4th order +T energy contribution
                  ierr=initialize_tensor_t(y_tensor)
                  ierr=t_einsum("W(a,b,c)+=W(a,b,c)",y_tensor,w_tensor)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  if(print_level > 1) write(*,*) " W[/D ] =",tensor_norm2(w_tensor)
                  ierr=initialize_tensor_t(result_tensor)
                  ierr=t_einsum("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(1)=t_energy(1) + result_tens(1)
               end do
            end do
            if (print_level >1) write (*,*)  't3_energy:',t_energy
         end do

!        We looped over all (abc) combinations, scale with the factor of 1/6 needed to go to the unique ones.
         t_energy = t_energy / SIX

!        Destruct intermediate tensors
         deallocate(y)
         ierr=destroy_tensor_t(t2_i_tensor)
         ierr=destroy_tensor_t(t2_j_tensor)
         ierr=destroy_tensor_t(t2_k_tensor)
         ierr=destroy_tensor_t(t2_ij_tensor)
         ierr=destroy_tensor_t(t2_jk_tensor)
         ierr=destroy_tensor_t(t2_ki_tensor)
         ierr=destroy_tensor_t(vovv_i_tensor)
         ierr=destroy_tensor_t(vovv_j_tensor)
         ierr=destroy_tensor_t(vovv_k_tensor)
         ierr=destroy_tensor_t(ooov_ij_tensor)
         ierr=destroy_tensor_t(ooov_jk_tensor)
         ierr=destroy_tensor_t(ooov_ki_tensor)
         ierr=destroy_tensor_t(we_tensor)
         ierr=destroy_tensor_t(w_tensor)
         ierr=destroy_tensor_t(y_tensor)
         ierr=destroy_tensor_t(one_tripl)

       end subroutine t3_energy

       subroutine form_y2_intermediate (i,j,k,t1_tensor,t2_tensor,fov_tensor,oovv_tensor,y_tensor,y)

!        Get y intermediate for -T fifth order energy correction.

         type(tensor_t), intent(inout) :: t1_tensor, t2_tensor, y_tensor
         type(tensor_t), intent(inout) :: fov_tensor, oovv_tensor
         integer, intent(in)               :: i, j, k
         complex(8), intent(inout)         :: y(:,:,:)

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank
         complex(8), pointer  :: t1_tens(:,:), t2_tens(:,:,:,:), y_tens(:,:,:)
         complex(8), pointer  :: fov_tens(:,:),oovv_tens(:,:,:,:)
         integer(INTD)        :: t1_dims(2), t2_dims(4), y_dims(3), fov_dims(2), oovv_dims(4)
         integer              :: nvir, a, b, c

!        Determine source dimensions and get direct access to the tensor bodies

         ierr = get_nmodes(t1_tensor,tens_rank)
         if (tens_rank .ne. 2) stop 'Error: t1 tensor corrupted in form_y2_intermediate'
         ierr = get_extents_static(t1_tensor,t1_dims)
         ierr=direct_data_pointer(t1_tensor,body_p)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,t1_tens,t1_dims)

         ierr = get_nmodes(t2_tensor,tens_rank)
         if (tens_rank .ne. 4) stop 'Error: t2 tensor corrupted in form_y2_intermediate'
         ierr = get_extents_static(t2_tensor,t2_dims)
         ierr=direct_data_pointer(t2_tensor,body_p)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,t2_tens,t2_dims)

         ierr = get_nmodes(y_tensor,tens_rank)
         if (tens_rank .ne. 3) stop 'Error: y tensor corrupted in form_y2_intermediate'
         ierr = get_extents_static(y_tensor,y_dims)
         ierr=direct_data_pointer(y_tensor,body_p)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,y_tens,y_dims)

         ierr = get_nmodes(fov_tensor,tens_rank)
         if (tens_rank .ne. 2) stop 'Error: fo tensor corrupted in form_y2_intermediate'
         ierr = get_extents_static(fov_tensor,fov_dims)
         ierr=direct_data_pointer(fov_tensor,body_p)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,fov_tens,fov_dims)

         ierr = get_nmodes(oovv_tensor,tens_rank)
         if (tens_rank .ne. 4) stop 'Error: oovv tensor corrupted in form_y2_intermediate'
         ierr = get_extents_static(oovv_tensor,oovv_dims)
         ierr=direct_data_pointer(oovv_tensor,body_p)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,oovv_tens,oovv_dims)

!        Form the intermediate
         nvir = t1_dims(1)

!        Formula with only the P_ijk permutation operators
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y(a,b,c) =            dconjg(oovv_tens(i,j,a,b))*t1_tens(c,k) + t2_tens(a,b,i,j)*dconjg(fov_tens(k,c))
                  y(a,b,c) = y(a,b,c) + dconjg(oovv_tens(j,k,a,b))*t1_tens(c,i) + t2_tens(a,b,j,k)*dconjg(fov_tens(i,c))
                  y(a,b,c) = y(a,b,c) + dconjg(oovv_tens(k,i,a,b))*t1_tens(c,j) + t2_tens(a,b,k,i)*dconjg(fov_tens(j,c))
               end do
            end do
         end do

!        Do the P_abc permutation to fully symmetrize the expression
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y_tens(a,b,c)  = y(a,b,c) + y(c,a,b) + y(b,c,a)
               end do
            end do
         end do

       end subroutine form_y2_intermediate

       subroutine form_y3_intermediate (i,j,k,t1_tensor,t2_tensor,y_tensor,y)

!        Get y intermediate for -T fifth order energy correction.

         type(tensor_t), intent(inout) :: t1_tensor, t2_tensor, y_tensor
         integer, intent(in)               :: i, j, k
         complex(8), intent(inout)         :: y(:,:,:)

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank
         complex(8), pointer  :: t1_tens(:,:), t2_tens(:,:,:,:), y_tens(:,:,:)
         integer(INTD)        :: t1_dims(2), t2_dims(4), y_dims(3)
         integer              :: nvir, a, b, c

!        Determine source dimensions and get direct access to the tensor bodies

         ierr = get_nmodes(t1_tensor,tens_rank)
         if (tens_rank .ne. 2) stop 'Error: t1 tensor corrupted in form_y_intermediate'
         ierr = get_extents_static(t1_tensor,t1_dims)
         ierr=direct_data_pointer(t1_tensor,body_p)
         call c_f_pointer(body_p,t1_tens,t1_dims)

         ierr = get_nmodes(t2_tensor,tens_rank)
         if (tens_rank .ne. 4) stop 'Error: t2 tensor corrupted in form_y_intermediate'
         ierr = get_extents_static(t2_tensor,t2_dims)
         ierr=direct_data_pointer(t2_tensor,body_p)
         call c_f_pointer(body_p,t2_tens,t2_dims)

         ierr = get_nmodes(y_tensor,tens_rank)
         if (tens_rank .ne. 3) stop 'Error: y tensor corrupted in form_y_intermediate'
         ierr = get_extents_static(y_tensor,y_dims)
         ierr=direct_data_pointer(y_tensor,body_p)
         call c_f_pointer(body_p,y_tens,y_dims)

!        Form the intermediate
         nvir = t1_dims(1)

!        Formula with only the P_ijk and P_ij permutation operators
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y(a,b,c) =            t1_tens(a,i)*t1_tens(b,j)*t1_tens(c,k) - t1_tens(a,j)*t1_tens(b,i)*t1_tens(c,k)
                  y(a,b,c) = y(a,b,c) + t1_tens(a,j)*t1_tens(b,k)*t1_tens(c,i) - t1_tens(a,k)*t1_tens(b,j)*t1_tens(c,i)
                  y(a,b,c) = y(a,b,c) + t1_tens(a,k)*t1_tens(b,i)*t1_tens(c,j) - t1_tens(a,i)*t1_tens(b,k)*t1_tens(c,j)
                  y(a,b,c) = y(a,b,c) / THREE
                  y(a,b,c) = y(a,b,c) + t2_tens(a,b,i,j)*t1_tens(c,k)+t2_tens(a,b,j,k)*t1_tens(c,i)+t2_tens(a,b,k,i)*t1_tens(c,j)
               end do
            end do
         end do

!        Do the P_abc permutation to fully symmetrize the expression
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y_tens(a,b,c)  = y(a,b,c) + y(c,a,b) + y(b,c,a)
               end do
            end do
         end do

       end subroutine form_y3_intermediate

       subroutine get_subtensor (tensor,subtensor,mode,i,j)

!        Severely limited implementation of extraction of a subtensor, as needed in t3_energy routine.
!        No shuffling of indices possible, only 4->3 (mode=1) and 4->2 (mode=2) extraction.
!        Type of extraction defined by 3 or 4 digit integer mode:
!        1st integer: dimension of source tensor
!        2nd integer: dimension of target tensor
!        3rd integer: passive index #1
!        4th integer: passive index #2

         type(tensor_t), intent(inout) :: tensor, subtensor
         integer, intent(in)               :: mode,i
         integer, intent(in), optional     :: j

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank, subtens_rank
         complex(8), pointer  :: tens(:,:,:,:), subtens3(:,:,:), subtens2(:,:)
         integer(INTD)        :: tens_dims(4), subtens_dims(3)

!        Determine source dimensions and get direct access to the tensor body to enable copying
         ierr = get_nmodes(tensor,tens_rank)
         if (tens_rank .ne. 4) stop 'Error: Only extraction from 4-dimensional tensors possible'
         ierr = get_extents_static(tensor,tens_dims)
         ierr=direct_data_pointer(tensor,body_p)
         call c_f_pointer(body_p,tens,tens_dims)

!        Determine target  dimensions and get direct access to the tensor body to enable copying
         ierr = get_nmodes(subtensor,subtens_rank)
         ierr = get_extents_static(subtensor,subtens_dims)
         if (subtens_rank .eq. 3) then
            ierr=direct_data_pointer(subtensor,body_p)
            call c_f_pointer(body_p,subtens3,subtens_dims)
         elseif (subtens_rank .eq. 2) then
            ierr=direct_data_pointer(subtensor,body_p)
            call c_f_pointer(body_p,subtens2,subtens_dims(1:2))
         else
            stop 'Error: Only insertion into 3- or 2-dimensional tensors possible'
         end if

!        Copy the data into the subtensor
         if (mode .eq. 433) then
            subtens3(:,:,:) = tens(:,:,i,:)
         elseif (mode .eq. 432) then
            subtens3(:,:,:) = tens(:,i,:,:)
         elseif (mode .eq. 4312) then
            subtens2(:,:)   = tens(i,j,:,:)
         elseif (mode .eq. 4334) then
            subtens2(:,:)   = tens(:,:,i,j)
         else
            write (*,*) mode,' type of tensor extraction is not implemented'
            stop 'Error in get_subtensor'
         end if

       end subroutine get_subtensor

       subroutine t3_energy_laplace (t_energy,comm_t,int_t,eps_occ,eps_vir,exa_input)

!       Evaluate +T, -T and (T) energy terms.
!       Definitions of terms: L. Visscher, T.J. Lee, K.G. Dyall, J Chem Phys. 105 (1996) 8769.
!       Implementation by Johann Pototschnig, Summer 2021
        use type_laplace, only : laplace_info, deallocate_laplace
        use exacorr_cfg, only: exacc_input
        real(8), intent(out)                 :: t_energy(3)
        type(eti_comm_tens_t), intent(inout) :: comm_t
        type(eti_intg_tens_t), intent(inout) :: int_t
        real(8), intent(in)                  :: eps_occ(:),eps_vir(:)
        type(exacc_input), intent(in)        :: exa_input

        !  variables to get laplace points
        type(laplace_info)    :: laplace

        real(8), allocatable  :: zero_occ(:),zero_vir(:),eps_occ_neg(:)
        complex(8)            :: T(40),Res(40)
        integer(INTD)         :: ierr, tens_rank
        integer(INTD)         :: t2_dims(4)
        integer               :: nocc, nvir, i

!       Get nocc and nvir from tensor information
        ierr = get_nmodes(comm_t%t2,tens_rank)
        ierr = get_extents_static(comm_t%t2,t2_dims)
        if (ierr.ne.0 .or. tens_rank.ne.4) stop 'error: t2 tensor corrupted'

        if (t2_dims(1).eq.t2_dims(2)) then
          nvir = t2_dims(1)
        else
          stop 'error: asymmetric tensor in t3_energy'
        end if
        if (t2_dims(3).eq.t2_dims(4)) then
          nocc = t2_dims(3)
        else
          stop 'error: asymmetric tensor in t3_energy'
        end if

!       get a list of positive occupied orbitals
        allocate(eps_occ_neg(nocc))
        eps_occ_neg(1:nocc)=-eps_occ(1:nocc)
!       create list of empty orbital
        allocate(zero_occ(nocc))
        allocate(zero_vir(nvir))
        zero_occ=0.0
        zero_vir=0.0
!       determine Laplace quadrature parameters with the minimax algorithm
        laplace%bounds(1) = 3.0d0*(MINVAL(eps_vir)-MAXVAL(eps_occ))
        laplace%bounds(2) = 3.0d0*(MAXVAL(eps_vir)-MINVAL(eps_occ))

        if (exa_input%print_level>5) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- occ orb. :', eps_occ
          write(*,*) '--- vir orb. :', eps_vir
          write(*,*) '--- limits :',MINVAL(eps_vir),MAXVAL(eps_occ),MAXVAL(eps_vir),MINVAL(eps_occ)
        end if

        
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

!       start computation of triples
        t_energy=0

        Res=DCMPLX(0.0D0)
        do i=1,laplace%num
          T=DCMPLX(0.0D0)
          call laplace_get_T1(T,comm_t%t2,int_t%ooov,int_t%vovv,t2_dims,  &
               eps_occ_neg,eps_vir,zero_occ,zero_vir,laplace%val(i)/2,exa_input%print_level)
          t_energy(1)=t_energy(1)+laplace%wei(i)*real(sum(T))

          write(*,*) 'Terms*Lapl:'
          write(*,'(10e10.2)') laplace%wei(i)*T

          Res=Res+laplace%wei(i)*T

        end do
        
        write(*,*) 'Sums of Res 1: '
        write(*,'(6e16.8)') sum(Res(1:9)),sum(Res(10:19)),sum(Res(20:29)), sum(Res(30:39)), sum(Res)

        Res=DCMPLX(0.0D0)
        do i=1,laplace%num
          T=DCMPLX(0.0D0)
          call laplace_get_T2(T,comm_t%t2,comm_t%t1,comm_t%fov,int_t%ooov,int_t%vovv,int_t%oovv,t2_dims, &
               eps_occ_neg,eps_vir,zero_occ,zero_vir,laplace%val(i)/2,exa_input%print_level)
          t_energy(2)=t_energy(2)+laplace%wei(i)*real(sum(T))
          
          write(*,*) 'Terms*Lapl:'
          write(*,'(10e10.2)') laplace%wei(i)*T

          Res=Res+laplace%wei(i)*T

        end do 

        write(*,*) 'Sums of Res 2: '
        write(*,'(6e16.8)') sum(Res(1:9)),sum(Res(10:19)),sum(Res(20:29)), sum(Res(30:39)), sum(Res)

        T=DCMPLX(0.0D0)
        call laplace_get_T3(T,comm_t%t2,comm_t%t1,int_t%ooov,int_t%vovv,t2_dims,exa_input%print_level)
        t_energy(3)=t_energy(3)-real(sum(T))

        Res=T
        write(*,*) 'Sums of Res 3: '
        write(*,'(6e16.8)') sum(Res(1:8)), sum(Res(9:20)), sum(Res(21:24)), &
                              sum(Res(25:34)), sum(Res(35:38)), sum(Res)

        t_energy=-t_energy  
        write(*,*) 'Energy: ',t_energy

!       clean up
        call deallocate_laplace(laplace)
        deallocate(zero_occ)
        deallocate(zero_vir)
        deallocate(eps_occ_neg)
    end subroutine t3_energy_laplace

    subroutine laplace_get_T1(T,i_t2,i_ooov,i_vovv,t2_dims, &
               eps_occ,eps_vir,zero_occ,zero_vir,s,print_level)

      complex(8), intent(out)           :: T(:)
      type(tensor_t), intent(inout) :: i_t2,i_ooov,i_vovv
      real(8), intent(in)               :: eps_occ(:),eps_vir(:),zero_occ(:),zero_vir(:)
      real(8), intent(in)               :: s
      integer(INTD), intent(in)         :: t2_dims(4)
      integer, intent(in)               :: print_level

      type(tensor_t)  :: t2_o,t2_v,vooo_o,vvov_v, I1, I2
      integer(INTD)       :: ierr
      integer(INTD)       :: dims4(4), dims2(2)
      integer             :: nvir, nocc

      ! scalars (need to be defined as tensor types)
      type(C_PTR)         :: body_p
      type(tensor_t)  :: T_tensor
      integer(INTD)       :: T_dims(1)
      complex(8), pointer :: T_tens(:)

      integer(INTD)      :: one_dims(1)
      type(tensor_t) :: one_tripl
      !construct one tensor
      ierr=construct_tensor_t(one_tripl,C8,one_dims(1:0),init_val=ONE)

      if (print_level>1) call print_date('Entered laplace_get_T1 routine')
      if (print_level>1) write(*,*) '-------------'

      !Initialize scalars that are to be used as tensors in contractions
      T_dims(1) = 1
      ierr=construct_tensor_t(T_tensor,C8,T_dims(1:0),init_val=ZERO)
      ierr=direct_data_pointer(T_tensor,body_p)
      call c_f_pointer(body_p,T_tens,T_dims)
     
      nvir = t2_dims(1)
      nocc = t2_dims(3)
      T_tens(1)=ZERO

      ! ------------------- compute first transformed integrals and amplitudes
      !ierr=construct_tensor_t(t2_v,C8,t2_dims,init_val=ONE/SQRT(DCMPLX(nvir)))  !+++++ testing
      ierr=construct_tensor_t(t2_v,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(e,a,i,j)+=O(e,a,i,j)",t2_v,i_t2)
      call laplace_mult(t2_v,s,zero_vir,eps_vir,eps_occ,eps_occ)

      dims4=nvir
      dims4(3)=nocc
      !ierr=construct_tensor_t(vvov_v,C8,dims4,init_val=ONE/SQRT(DCMPLX(nvir))) !+++++ testing
      ierr=construct_tensor_t(vvov_v,C8,dims4,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,e)+=O+(e,i,b,a)",vvov_v,i_vovv,one_tripl)
      call laplace_mult(vvov_v,s,eps_vir,eps_vir,eps_occ,zero_vir)

      ! ------------------- compute first transformed integrals and amplitudes

      ! 0.25 T2^{ef}_{}G2^{ef}_{}
      dims2=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,f)+=Y+(e,c,j,k)*X(f,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T2 wrong')
      ierr=t_einsum("G(e,f)+=Y+(a,b,i,e)*X(a,b,i,f)",I2,vvov_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G2 wrong')
      ierr=t_einsum("R()+=T(e,f)*G(e,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 1 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(1)=T_tens(1)
      

      ! 0.25 T5^{ef}_{ki}G4^{}_{iekf}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,f,k,i)+=Y+(e,c,j,k)*X(f,c,i,j)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T5 wrong')
      ierr=t_einsum("G(e,f,k,i)+=Y+(a,b,i,e)*X(a,b,k,f)",I2,vvov_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G4 wrong')
      ierr=t_einsum("R()+=T(e,f,k,i)*G(e,f,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 2 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(2)=T_tens(1)

      ! 0.25 T6^{ef}_{ji}G4^{}_{}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,f,j,i)+=Y+(e,c,j,k)*X(f,c,k,i)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T6 wrong')
      ierr=t_einsum("G(e,f,j,i)+=Y+(a,b,i,e)*X(a,b,j,f)",I2,vvov_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G4 wrong')
      ierr=t_einsum("R()+=T(e,f,j,i)*G(e,f,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 3 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(3)=T_tens(1)

      ! 0.25 T4^{ecfb}_{}G7^{bc}_{ef}
      T_tens(1)=ZERO
      dims4=nvir
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(3:4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("G(e,c,f,b)+=Y+(a,b,i,e)*X(c,a,i,f)",I1,vvov_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G7 wrong')
      ierr=t_einsum("T(f,b,j,k)+=Y+(e,c,j,k)*X(e,c,f,b)",I2,t2_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T4 wrong')
      ierr=t_einsum("R()+=T(f,b,j,k)*G(f,b,j,k)",T_tensor,t2_v,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 4 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(4)=T_tens(1)

      ! 0.25 T4^{e,c,f,a}_{}G8^{ac}_{ef}
      T_tens(1)=ZERO
      dims4=nvir
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(3:4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("G(e,c,f,a)+=Y+(a,b,i,e)*X(b,c,i,f)",I1,vvov_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G8 wrong')
      ierr=t_einsum("T(f,a,j,k)+=Y+(e,c,j,k)*X(e,c,f,a)",I2,t2_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T4 wrong')
      ierr=t_einsum("R()+=T(f,a,j,k)*G(f,a,j,k)",T_tensor,t2_v,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 5 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(5)=T_tens(1)

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I9^{ea}_{kf}J13^{af}_{ek}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,a,k,f)+=X+(e,c,j,k)*Y(c,a,j,f)",I1,t2_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I9 wrong')
      ierr=t_einsum("N(f,b,k,i)+=X+(a,b,i,e)*Y(e,a,k,f)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J13 wrong')
      ierr=t_einsum("R()+=M(f,b,k,i)*N(f,b,k,i)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 6 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(6)=T_tens(1)

      ! 0.25 I13^{eb}_{jf}J9^{bf}_{ej}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,b,j,f)+=X+(e,c,j,k)*Y(b,c,k,f)",I1,t2_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I13 wrong')
      ierr=t_einsum("N(f,a,i,j)+=X+(a,b,i,e)*Y(e,b,j,f)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J9 wrong')
      ierr=t_einsum("R()+=M(f,a,i,j)*N(f,a,i,j)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 7 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(7)=T_tens(1)      

      ! 0.25 I14^{ea}_{jf}J6^{af}_{ej}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,a,j,f)+=X+(e,c,j,k)*Y(c,a,k,f)",I1,t2_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I14 wrong')
      ierr=t_einsum("N(f,b,i,j)+=X+(a,b,i,e)*Y(e,a,j,f)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J6 wrong')
      ierr=t_einsum("R()+=M(f,b,i,j)*N(f,b,i,j)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 8 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(8)=T_tens(1)

      ! 0.25 I6^{eb}_{kf}J14^{bf}_{ek}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,b,k,f)+=X+(e,c,j,k)*Y(b,c,j,f)",I1,t2_v,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I6 wrong')
      ierr=t_einsum("N(f,a,k,i)+=X+(a,b,i,e)*Y(e,b,k,f)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J14 wrong')
      ierr=t_einsum("R()+=M(f,a,k,i)*N(f,a,k,i)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 9 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(9)=T_tens(1)

      ! ------------------- compute next transformed integrals and amplitudes
      !ierr=construct_tensor_t(t2_o,C8,t2_dims,init_val=ONE/SQRT(DCMPLX(nocc))) !+++++ testing
      ierr=construct_tensor_t(t2_o,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,m)+=O(a,b,i,m)",t2_o,i_t2)
      call laplace_mult(t2_o,s,eps_vir,eps_vir,eps_occ,zero_occ)

      dims4=nocc
      dims4(1)=nvir
      !ierr=construct_tensor_t(ooov_o,C8,dims4,init_val=ONE/SQRT(DCMPLX(nocc)))    !+++++ testing
      ierr=construct_tensor_t(vooo_o,C8,dims4,init_val=ZERO)
      ierr=t_einsum("R(a,m,i,j)+=O+(j,i,m,a)",vooo_o,i_ooov,one_tripl)
      call laplace_mult(vooo_o,s,eps_vir,zero_occ,eps_occ,eps_occ)
      ! ------------------- compute next transformed integrals and amplitudes

      ! 0.25 T11^{cf}_{mi}G9^{mc}_{if}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(c,f,m,i)+=X+(b,c,k,m)*Y(f,b,k,i)",I1,t2_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T11 wrong')
      ierr=t_einsum("G(c,f,m,i)+=X+(a,m,i,j)*Y(c,a,j,f)",I2,vooo_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G9 wrong')
      ierr=t_einsum("R()+=T(c,f,m,i)*G(c,f,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 21 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(21)=T_tens(1)

      ! 0.25 T10^{cf}_{mj}G10^{mc}_{jf}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(c,f,m,j)+=X+(b,c,k,m)*Y(f,b,j,k)",I1,t2_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T10 wrong')
      ierr=t_einsum("G(c,f,m,j)+=X+(a,m,i,j)*Y(c,a,i,f)",I2,vooo_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G10 wrong')
      ierr=t_einsum("R()+=T(c,f,m,j)*G(c,f,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 22 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(22)=T_tens(1)

      ! 0.25 T12^{ea}_{kn}G11^{an}_{ek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,k,n)+=X+(e,c,j,k)*Y(c,a,j,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T12 wrong')
      ierr=t_einsum("G(e,a,k,n)+=X+(a,b,i,e)*Y(b,n,k,i)",I2,vvov_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G11 wrong')
      ierr=t_einsum("R()+=T(e,a,k,n)*G(e,a,k,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 23 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(23)=T_tens(1)

      ! 0.25 T9^{ea}_{jn}G12^{an}_{ej}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,j,n)+=X+(e,c,j,k)*Y(c,a,k,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T9 wrong')
      ierr=t_einsum("G(e,a,j,n)+=X+(a,b,i,e)*Y(b,n,i,j)",I2,vvov_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G12 wrong')
      ierr=t_einsum("R()+=T(e,a,j,n)*G(e,a,j,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 24 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(24)=T_tens(1)

      ! 0.25 T13^{bf}_{mj}G13^{nb}_{je}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,f,m,j)+=X+(b,c,k,m)*Y(f,c,j,k)",I1,t2_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T13 wrong')
      ierr=t_einsum("G(b,f,m,j)+=X+(a,m,i,j)*Y(a,b,i,f)",I2,vooo_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G13 wrong')
      ierr=t_einsum("R()+=T(b,f,m,j)*G(b,f,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 25 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(25)=T_tens(1)

      ! 0.25 T15^{bf}_{mi}G14^{mb}_{if}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,f,m,i)+=X+(b,c,k,m)*Y(f,c,k,i)",I1,t2_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T15 wrong')
      ierr=t_einsum("G(b,f,m,i)+=X+(a,m,i,j)*Y(a,b,j,f)",I2,vooo_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G14 wrong')
      ierr=t_einsum("R()+=T(b,f,m,i)*G(b,f,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 26 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(26)=T_tens(1)

      ! 0.25 T16^{eb}_{kn}G15^{bn}_{ek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,b,k,n)+=X+(e,c,j,k)*Y(b,c,j,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T16 wrong')
      ierr=t_einsum("G(e,b,k,n)+=X+(a,b,i,e)*Y(a,n,k,i)",I2,vvov_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G15 wrong')
      ierr=t_einsum("R()+=T(e,b,k,n)*G(e,b,k,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 27 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(27)=T_tens(1)

      ! 0.25 T14^{eb}_{jn}G16^{bn}_{ej}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,b,j,n)+=X+(e,c,j,k)*Y(b,c,k,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T14 wrong')
      ierr=t_einsum("G(e,b,j,n)+=X+(a,b,i,e)*Y(a,n,i,j)",I2,vvov_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G16 wrong')
      ierr=t_einsum("R()+=T(e,b,j,n)*G(e,b,j,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 28 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(28)=T_tens(1)

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I3^{}_{en}J3^{}_{en}
      dims2=nvir
      dims2(2)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("M(e,n)+=X+(e,c,j,k)*Y(c,n,j,k)",I1,t2_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I3 wrong')
      ierr=t_einsum("N(e,n)+=X+(a,b,i,e)*Y(a,b,i,n)",I2,vvov_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J3 wrong')
      ierr=t_einsum("R()+=M(e,n)*N(e,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 30 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(30)=T_tens(1)

      ! 0.25 I4^{}_{mf}I4^{}_{mf}
      dims2=nvir
      dims2(1)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("M(m,f)+=X+(b,c,k,m)*Y(b,c,k,f)",I1,t2_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I4 wrong')
      ierr=t_einsum("N(m,f)+=X+(a,m,i,j)*Y(f,a,i,j)",I2,vooo_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J4 wrong')
      ierr=t_einsum("R()+=M(m,f)*N(m,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 31 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(31)=T_tens(1)

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I1^{}_{kmif}J16^{mf}_{ik}
      dims4=nocc
      dims4(4)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(k,m,i,f)+=X+(b,c,k,m)*Y(b,c,i,f)",I1,t2_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=t_einsum("N(k,m,i,f)+=X+(a,m,i,j)*Y(f,a,j,k)",I2,vooo_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J16 wrong')
      ierr=t_einsum("F()+=M(k,m,i,f)*N(k,m,i,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 32 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(32)=T_tens(1)

      ! 0.25 I16^{en}_{ki}J1^{}_{iekn}
      dims4=nocc
      dims4(1)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,n,k,i)+=X+(e,c,j,k)*Y(c,n,i,j)",I1,t2_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I16 wrong')
      ierr=t_einsum("N(e,n,k,i)+=X+(a,b,i,e)*Y(a,b,k,n)",I2,vvov_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J1 wrong')
      ierr=t_einsum("F()+=M(e,n,k,i)*N(e,n,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 33 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(33)=T_tens(1)


      ! 0.25 I1^{}_{kmjf}J7^{mf}_{jk}
      dims4=nocc
      dims4(4)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(k,m,j,f)+=X+(b,c,k,m)*Y(b,c,j,f)",I1,t2_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=t_einsum("N(k,m,j,f)+=X+(a,m,i,j)*Y(f,a,k,i)",I2,vooo_o,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J7 wrong')
      ierr=t_einsum("F()+=M(k,m,j,f)*N(k,m,j,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 34 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(34)=T_tens(1)

      ! 0.25 I7^{en}_{ji}J1^{}_{iejn}
      dims4=nocc
      dims4(1)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,n,j,i)+=X+(e,c,j,k)*Y(c,n,k,i)",I1,t2_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=t_einsum("N(e,n,j,i)+=X+(a,b,i,e)*Y(a,b,j,n)",I2,vvov_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J7 wrong')
      ierr=t_einsum("F()+=M(e,n,j,i)*N(e,n,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 35 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(35)=T_tens(1)

      ! 0.25 I10^{ca}_{mf}J2^{ca}_{ne}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(c,a,m,f)+=X+(b,c,k,m)*Y(a,b,k,f)",I1,t2_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I10 wrong')
      ierr=t_einsum("N(f,c,i,j)+=X+(a,m,i,j)*Y(c,a,m,f)",I2,vooo_o,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J2 wrong')
      ierr=t_einsum("F()+=M(f,c,i,j)*N(f,c,i,j)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 36 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(36)=T_tens(1)

      ! 0.25 I2^{ecan}_{}J10^{ac}_{en}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(4)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(3)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,c,a,n)+=X+(e,c,j,k)*Y(a,n,j,k)",I1,t2_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I2 wrong')
      ierr=t_einsum("N(b,c,i,n)+=X+(a,b,i,e)*Y(e,c,a,n)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J10 wrong')
      ierr=t_einsum("F()+=M(b,c,i,n)*N(b,c,i,n)",T_tensor,I2,t2_o,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 37 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(37)=T_tens(1)

      ! 0.25 I11^{ba}_{mf}J2^{amfb}_{}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(4)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(b,a,m,f)+=X+(b,c,k,m)*Y(c,a,k,f)",I1,t2_o,vvov_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I11 wrong')
      ierr=t_einsum("N(f,b,i,j)+=X+(a,m,i,j)*Y(b,a,m,f)",I2,vooo_o,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J2 wrong')
      ierr=t_einsum("F()+=M(f,b,i,j)*N(f,b,i,j)",T_tensor,I2,t2_v,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 38 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(38)=T_tens(1)

      ! 0.25 I2^{ecbn}_{}J11^{bc}_{en}
      T_tens(1)=ZERO
      dims4=nvir
      dims4(4)=nocc
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      dims4(3)=nocc
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(e,c,b,n)+=X+(e,c,j,k)*Y(b,n,j,k)",I1,t2_v,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I11 wrong')
      ierr=t_einsum("N(c,a,i,n)+=X+(a,b,i,e)*Y(e,c,b,n)",I2,vvov_v,I1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J2 wrong')
      ierr=t_einsum("F()+=M(c,a,i,n)*N(c,a,i,n)",T_tensor,I2,t2_o,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 39 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(39)=T_tens(1)

      ! ------------------- delete first transformed integrals and amplitudes
      ierr=destroy_tensor_t(t2_v)
      ierr=destroy_tensor_t(vvov_v)
      ! ------------------- delete first transformed integrals and amplitudes

      ! 0.25 T1^{}_{mn}G1^{}_{nm}
      dims2=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,n)+=X+(b,c,k,m)*Y(b,c,k,n)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T1 wrong')
      ierr=t_einsum("G(m,n)+=X+(a,m,i,j)*Y(a,n,i,j)",I2,vooo_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G1 wrong')
      ierr=t_einsum("R()+=T(m,n)*G(m,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 11 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(11)=T_tens(1)


      ! 0.25 T3^{}_{kmjn}G5^{mn}_{jk}
      dims4=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(k,m,j,n)+=X+(b,c,k,m)*Y(b,c,j,n)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T3 wrong')
      ierr=t_einsum("G(k,m,j,n)+=X+(a,m,i,j)*Y(a,n,k,i)",I2,vooo_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G5 wrong')
      ierr=t_einsum("R()+=T(k,m,j,n)*G(k,m,j,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 12 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(12)=T_tens(1)

      ! 0.25 T3^{}_{kmin}G6^{mn}_{ik}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(k,m,i,n)+=X+(b,c,k,m)*Y(b,c,i,n)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T3 wrong')
      ierr=t_einsum("G(k,m,i,n)+=X+(a,m,i,j)*Y(a,n,j,k)",I2,vooo_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G6 wrong')
      ierr=t_einsum("R()+=T(k,m,i,n)*G(k,m,i,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 13 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(13)=T_tens(1)

      ! 0.25 T7^{ca}_{mn}G3^{amcn}_{}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(c,a,m,n)+=X+(b,c,k,m)*Y(a,b,k,n)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T7 wrong')
      ierr=t_einsum("G(c,a,m,n)+=X+(a,m,i,j)*Y(c,n,i,j)",I2,vooo_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G3 wrong')
      ierr=t_einsum("R()+=T(c,a,m,n)*G(c,a,m,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 14 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(14)=T_tens(1)

      ! 0.25 T8^{ba}_{mn}G3^{amcn}_{}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,a,m,n)+=X+(b,c,k,m)*Y(c,a,k,n)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T8 wrong')
      ierr=t_einsum("G(b,a,m,n)+=X+(a,m,i,j)*Y(b,n,i,j)",I2,vooo_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G3 wrong')
      ierr=t_einsum("R()+=T(b,a,m,n)*G(b,a,m,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 15 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(15)=T_tens(1)

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I5^{cn}_{mi}J15^{mc}_{in}
      dims4=nocc
      dims4(1)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(c,n,m,i)+=X+(b,c,k,m)*Y(b,n,k,i)",I1,t2_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I5 wrong')
      ierr=t_einsum("N(c,n,m,i)+=X+(a,m,i,j)*Y(c,a,j,n)",I2,vooo_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J15 wrong')
      ierr=t_einsum("R()+=M(c,n,m,i)*N(c,n,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 16 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(16)=T_tens(1)

      ! 0.25 I15^{bn}_{mj}J5^{mb}_{jn}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(b,n,m,j)+=X+(b,c,k,m)*Y(c,n,j,k)",I1,t2_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I15 wrong')
      ierr=t_einsum("N(b,n,m,j)+=X+(a,m,i,j)*Y(a,b,i,n)",I2,vooo_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J5 wrong')
      ierr=t_einsum("R()+=M(b,n,m,j)*N(b,n,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 17 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(17)=T_tens(1)

      ! 0.25 I8^{cn}_{mj}J12^{mc}_{jn}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(c,n,m,j)+=X+(b,c,k,m)*Y(b,n,j,k)",I1,t2_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I8 wrong')
      ierr=t_einsum("N(c,n,m,j)+=X+(a,m,i,j)*Y(c,a,i,n)",I2,vooo_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J12 wrong')
      ierr=t_einsum("R()+=M(c,n,m,j)*N(c,n,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 18 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(18)=T_tens(1)

      ! 0.25 I12^{bn}_{mi}J8^{mb}_{in}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("M(b,n,m,i)+=X+(b,c,k,m)*Y(c,n,k,i)",I1,t2_o,vooo_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I12 wrong')
      ierr=t_einsum("N(b,n,m,i)+=X+(a,m,i,j)*Y(a,b,j,n)",I2,vooo_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction J8 wrong')
      ierr=t_einsum("R()+=M(b,n,m,i)*N(b,n,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 19 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(19)=T_tens(1)

      !clean up
      ierr=destroy_tensor_t(T_tensor)
      ierr=destroy_tensor_t(t2_o)
      ierr=destroy_tensor_t(vooo_o)
      ierr=destroy_tensor_t(one_tripl)

    end subroutine laplace_get_T1


    subroutine laplace_get_T2(T,i_t2,i_t1,i_fov,i_ooov,i_vovv,i_oovv,t2_dims, &
               eps_occ,eps_vir,zero_occ,zero_vir,s,print_level)

      complex(8), intent(out)           :: T(:)
      type(tensor_t), intent(inout) :: i_t2,i_t1,i_fov
      type(tensor_t), intent(inout) :: i_ooov,i_vovv,i_oovv
      real(8), intent(in)               :: eps_occ(:),eps_vir(:),zero_occ(:),zero_vir(:)
      real(8), intent(in)               :: s
      integer(INTD), intent(in)         :: t2_dims(4)
      integer, intent(in)               :: print_level

      type(tensor_t)  :: t2_o, t2_v, t2, t1
      type(tensor_t)  :: vvov_v, vooo_o, vvoo, fvo
      type(tensor_t)  :: I1, I2
      integer(INTD)       :: ierr
      integer(INTD)       :: dims4(4), dims2(2)
      integer             :: nvir, nocc

      ! scalars (need to be defined as tensor types)
      type(C_PTR)         :: body_p
      type(tensor_t)  :: T_tensor
      integer(INTD)       :: T_dims(1)
      complex(8), pointer :: T_tens(:)

      integer(INTD)      :: one_dims(1)
      type(tensor_t) :: one_tripl
      !construct one tensor
      ierr=construct_tensor_t(one_tripl,C8,one_dims(1:0),init_val=ONE)

      if (print_level>1) call print_date('Entered laplace_get_T2 routine')
      if (print_level>1) write(*,*) '-------------'

      !Initialize scalars that are to be used as tensors in contractions
      T_dims(1) = 1
      ierr=construct_tensor_t(T_tensor,C8,T_dims(1:0),init_val=ZERO)
      ierr=direct_data_pointer(T_tensor,body_p)
      call c_f_pointer(body_p,T_tens,T_dims)
     
      nvir = t2_dims(1)
      nocc = t2_dims(3)
      T_tens(1)=ZERO

      ! ------------------- compute first transformed integrals and amplitudes


      ierr=construct_tensor_t(t2_v,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(e,a,i,j)+=O+(e,a,i,j)*I()",t2_v,i_t2,one_tripl)
      call laplace_mult(t2_v,s,zero_vir,eps_vir,eps_occ,eps_occ)

      ierr=construct_tensor_t(t2_o,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,m)+=O+(a,b,i,m)*I()",t2_o,i_t2,one_tripl)
      call laplace_mult(t2_o,s,eps_vir,eps_vir,eps_occ,zero_occ)

      dims4=nvir
      dims4(3)=nocc
      ierr=construct_tensor_t(vvov_v,C8,dims4,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,e)+=O(e,i,b,a)*I()",vvov_v,i_vovv,one_tripl)
      call laplace_mult(vvov_v,s,eps_vir,eps_vir,eps_occ,zero_vir)

      dims4=nocc
      dims4(1)=nvir
      ierr=construct_tensor_t(vooo_o,C8,dims4,init_val=ZERO)
      ierr=t_einsum("R(a,m,i,j)+=O(j,i,m,a)*I()",vooo_o,i_ooov,one_tripl)
      call laplace_mult(vooo_o,s,eps_vir,zero_occ,eps_occ,eps_occ)

      ! ----------------------------------------------------------------------------------------- compute terms with t1
      dims2=nvir
      dims2(2)=nocc
      ierr=construct_tensor_t(t1,C8,dims2,init_val=ZERO)
      ierr=t_einsum("R(a,i)+=O(a,i)",t1,i_t1)
      call laplace_mult(t1,s,eps_vir,eps_occ)

      dims4=nocc
      dims4(1:2)=nvir
      ierr=construct_tensor_t(vvoo,C8,dims4,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,j)+=O+(i,j,a,b)*I()",vvoo,i_oovv,one_tripl)
      call laplace_mult(vvoo,s,eps_vir,eps_vir,eps_occ,eps_occ)

      ! 0.25 H1^{mb}_{}U1^{b}_{m}
      dims2=nvir
      dims2(1)=nocc
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,b)+=X(a,m,i,j)*Y(a,b,i,j)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H1 wrong')
      ierr=t_einsum("G(m,b)+=X(b,c,k,m)*Y(c,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U1 wrong')
      ierr=t_einsum("R()+=T(m,b)*G(m,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 1 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(1)=T_tens(1)

      ! 0.25 H2^{mc}_{}U2^{c}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,c)+=X(a,m,i,j)*Y(c,a,i,j)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H2 wrong')
      ierr=t_einsum("G(m,c)+=X(b,c,k,m)*Y(b,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U2 wrong')
      ierr=t_einsum("R()+=T(m,c)*G(m,c)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 2 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(2)=T_tens(1)

      dims2=nocc
      dims2(1)=nvir
      ! 0.25 H3^{}_{ej}U3^{e}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,j)+=X(a,b,i,e)*Y(a,b,i,j)",I1,vvov_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H3 wrong')
      ierr=t_einsum("G(e,j)+=X(e,c,j,k)*Y(c,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U3 wrong')
      ierr=t_einsum("R()+=T(e,j)*G(e,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 3 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(3)=T_tens(1)  

      ! 0.25 H4^{}_{ek}U4^{e}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,k)+=X(a,b,i,e)*Y(a,b,k,i)",I1,vvov_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H4 wrong')
      ierr=t_einsum("G(e,k)+=X(e,c,j,k)*Y(c,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U4 wrong')
      ierr=t_einsum("R()+=T(e,k)*G(e,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 4 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(4)=T_tens(1)     

      ! 0.25 H5^{}_{iejk}U5^{e}_{jki}
      dims4=nocc
      dims4(2)=nvir
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(i,e,j,k)+=X(a,b,i,e)*Y(a,b,j,k)",I1,vvov_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H5 wrong')
      ierr=t_einsum("G(i,e,j,k)+=X(e,c,j,k)*Y(c,i)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U5 wrong')
      ierr=t_einsum("R()+=T(i,e,j,k)*G(i,e,j,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 5 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(5)=T_tens(1)

      ! 0.25 H6^{mb}_{ik}U6^{b}_{kmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(m,b,i,k)+=X(a,m,i,j)*Y(a,b,j,k)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H6 wrong')
      ierr=t_einsum("G(m,b,i,k)+=X(b,c,k,m)*Y(c,i)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U6 wrong')
      ierr=t_einsum("R()+=T(m,b,i,k)*G(m,b,i,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 6 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(6)=T_tens(1)

      ! 0.25 H7^{mc}_{ik}U7^{c}_{kmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(m,c,i,k)+=X(a,m,i,j)*Y(c,a,j,k)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H7 wrong')
      ierr=t_einsum("G(m,c,i,k)+=X(b,c,k,m)*Y(b,i)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U7 wrong')
      ierr=t_einsum("R()+=T(m,c,i,k)*G(m,c,i,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 7 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(7)=T_tens(1)

      ! 0.25 H8^{mc}_{jk}U8^{c}_{kmj}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(m,c,j,k)+=X(a,m,i,j)*Y(c,a,k,i)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H8 wrong')
      ierr=t_einsum("G(m,c,j,k)+=X(b,c,k,m)*Y(b,j)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U8 wrong')
      ierr=t_einsum("R()+=T(m,c,j,k)*G(m,c,j,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 8 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(8)=T_tens(1)

      ! 0.25 H9^{mb}_{jk}U9^{b}_{kmj}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(m,b,j,k)+=X(a,m,i,j)*Y(a,b,k,i)",I1,vooo_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H9 wrong')
      ierr=t_einsum("G(m,b,j,k)+=X(b,c,k,m)*Y(c,j)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U9 wrong')
      ierr=t_einsum("R()+=T(m,b,j,k)*G(m,b,j,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 9 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(9)=T_tens(1)

      dims4=nocc
      ! 0.25 H10^{ambc}_{}U10^{bca}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(k,m,i,j)+=X(b,c,k,m)*Y(b,c,i,j)",I1,t2_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H10 wrong')
      ierr=t_einsum("G(k,m,i,j)+=X(a,m,i,j)*Y(a,k)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U10 wrong')
      ierr=t_einsum("R()+=T(k,m,i,j)*G(k,m,i,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 10 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(10)=T_tens(1)

      dims4=nvir
      dims4(3:4)=nocc
      ! 0.25 H11^{acej}_{}U11^{ace}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,b,k,i)+=X(e,c,j,k)*Y(b,c,i,j)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H11 wrong')
      ierr=t_einsum("G(e,b,k,i)+=X(a,b,i,e)*Y(a,k)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U11 wrong')
      ierr=t_einsum("R()+=T(e,b,k,i)*G(e,b,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 11 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(11)=T_tens(1)

      ! 0.25 H12^{bcej}_{}U12^{ecb}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,k,i)+=X(e,c,j,k)*Y(c,a,i,j)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H12 wrong')
      ierr=t_einsum("G(e,a,k,i)+=X(a,b,i,e)*Y(b,k)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U12 wrong')
      ierr=t_einsum("R()+=T(e,a,k,i)*G(e,a,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 12 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(12)=T_tens(1)

      ! 0.25 H13^{acek}_{}U13^{eca}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,e,j,i)+=X(e,c,j,k)*Y(b,c,k,i)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H13 wrong')
      ierr=t_einsum("G(b,e,j,i)+=X(a,b,i,e)*Y(a,j)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U13 wrong')
      ierr=t_einsum("R()+=T(b,e,j,i)*G(b,e,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 13 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(13)=T_tens(1)

      ! 0.25 H14^{acek}_{}U14^{eca}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,j,i)+=X(e,c,j,k)*Y(c,a,k,i)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction H14 wrong')
      ierr=t_einsum("G(e,a,j,i)+=X(a,b,i,e)*Y(b,j)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction U14 wrong')
      ierr=t_einsum("R()+=T(e,a,j,i)*G(e,a,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 14 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(14)=T_tens(1)


      dims2=nvir
      ! 0.25 O1^{be}_{}V1^{b}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(b,e)+=X(e,c,j,k)*Y(b,c,j,k)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction O1 wrong')
      ierr=t_einsum("G(b,e)+=X(a,b,i,e)*Y(a,i)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction V1 wrong')
      ierr=t_einsum("R()+=T(b,e)*G(b,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 15 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(15)=T_tens(1) 

      ! 0.25 O2^{ea}_{}V2^{a}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(a,e)+=X(e,c,j,k)*Y(c,a,j,k)",I1,t2_v,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction O2 wrong')
      ierr=t_einsum("G(a,e)+=X(a,b,i,e)*Y(b,i)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction V2 wrong')
      ierr=t_einsum("R()+=T(a,e)*G(a,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 16 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(16)=T_tens(1) 

      dims2=nocc
      ! 0.25 O3^{}_{mj}V3^{m}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,j)+=X(b,c,k,m)*Y(b,c,j,k)",I1,t2_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction O3 wrong')
      ierr=t_einsum("G(m,j)+=X(a,m,i,j)*Y(a,i)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction V3 wrong')
      ierr=t_einsum("R()+=T(m,j)*G(m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 17 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(17)=T_tens(1) 

      ! 0.25 O4^{}_{mi}V4^{m}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,i)+=X(b,c,k,m)*Y(b,c,k,i)",I1,t2_o,vvoo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction O4 wrong')
      ierr=t_einsum("G(m,i)+=X(a,m,i,j)*Y(a,j)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction V4 wrong')
      ierr=t_einsum("R()+=T(m,i)*G(m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 18 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(18)=T_tens(1) 

      ierr=destroy_tensor_t(t1)
      ierr=destroy_tensor_t(vvoo)
      ! --------------------------------------------------------------------------------------- compute terms with fvo
      dims2=nvir
      dims2(2)=nocc
      ierr=construct_tensor_t(fvo,C8,dims2,init_val=ZERO)
      ierr=t_einsum("R(a,i)+=O+(i,a)",fvo,i_fov,one_tripl)
      call laplace_mult(fvo,s,eps_vir,eps_occ)

      ierr=construct_tensor_t(t2,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,j)+=O(a,b,i,j)*I()",t2,i_t2,one_tripl)
      call laplace_mult(t2,s,eps_vir,eps_vir,eps_occ,eps_occ)

      dims2=nvir
      ! 0.25 E1^{eb}_{}F1^{b}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,b)+=X(e,c,j,k)*Y(b,c,j,k)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction É1 wrong')
      ierr=t_einsum("G(e,b)+=X(a,b,i,e)*Y(a,i)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F1 wrong')
      ierr=t_einsum("R()+=T(e,b)*G(e,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 21 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(21)=T_tens(1) 

      ! 0.25 E2^{ea}_{}F2^{a}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,a)+=X(e,c,j,k)*Y(c,a,j,k)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E2 wrong')
      ierr=t_einsum("G(e,a)+=X(a,b,i,e)*Y(b,i)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F2 wrong')
      ierr=t_einsum("R()+=T(e,a)*G(e,a)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 22 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(22)=T_tens(1)

      dims2=nocc
      ! 0.25 E3^{}_{mj}F3^{m}_{j}  
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,j)+=X(b,c,k,m)*Y(b,c,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E3 wrong')
      ierr=t_einsum("G(m,j)+=X(a,m,i,j)*Y(a,i)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F3 wrong')
      ierr=t_einsum("R()+=T(m,j)*G(m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 23 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(23)=T_tens(1)

      ! 0.25 E4^{}_{mj}F4^{m}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,i)+=X(b,c,k,m)*Y(b,c,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E4 wrong')
      ierr=t_einsum("G(m,i)+=X(a,m,i,j)*Y(a,j)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F4 wrong')
      ierr=t_einsum("R()+=T(m,i)*G(m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 24 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(24)=T_tens(1)

      dims4=nocc
      ! 0.25 E5^{}_{kmij}F5^{m}_{ijk}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(k,m,i,j)+=X(b,c,k,m)*Y(b,c,i,j)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E5 wrong')
      ierr=t_einsum("G(k,m,i,j)+=X(a,m,i,j)*Y(a,k)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F5 wrong')
      ierr=t_einsum("R()+=T(k,m,i,j)*G(k,m,i,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 25 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(25)=T_tens(1)

      dims4(1)=nvir
      ! 0.25 E6^{ecab}_{}F6^{abc}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,i,j,k)+=X(a,b,i,e)*Y(a,b,j,k)",I1,vvov_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E6 wrong')
      ierr=t_einsum("G(e,i,j,k)+=X(e,c,j,k)*Y(c,i)",I2,t2_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F6 wrong')
      ierr=t_einsum("R()+=T(e,i,j,k)*G(e,i,j,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 26 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(26)=T_tens(1)

      dims4=nvir
      dims4(3:4)=nocc
      ! 0.25 E7^{eb}_{ki}F7^{b}_{iek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,b,k,i)+=X(e,c,j,k)*Y(b,c,i,j)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E7 wrong')
      ierr=t_einsum("G(e,b,k,i)+=X(a,b,i,e)*Y(a,k)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F7 wrong')
      ierr=t_einsum("R()+=T(e,b,k,i)*G(e,b,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 27 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(27)=T_tens(1)

      ! 0.25 E8^{eb}_{ki}F8^{b}_{iek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,b,j,i)+=X(e,c,j,k)*Y(b,c,k,i)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E8 wrong')
      ierr=t_einsum("G(e,b,j,i)+=X(a,b,i,e)*Y(a,j)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F8 wrong')
      ierr=t_einsum("R()+=T(e,b,j,i)*G(e,b,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 28 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(28)=T_tens(1)

      ! 0.25 E9^{ea}_{ki}F9^{a}_{iek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,k,i)+=X(e,c,j,k)*Y(c,a,i,j)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E9 wrong')
      ierr=t_einsum("G(e,a,k,i)+=X(a,b,i,e)*Y(b,k)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F9 wrong')
      ierr=t_einsum("R()+=T(e,a,k,i)*G(e,a,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 29 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(29)=T_tens(1)

      ! 0.25 E10^{ea}_{ji}F10^{a}_{iej}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(e,a,j,i)+=X(e,c,j,k)*Y(c,a,k,i)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E10 wrong')
      ierr=t_einsum("G(e,a,j,i)+=X(a,b,i,e)*Y(b,j)",I2,vvov_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F10 wrong')
      ierr=t_einsum("R()+=T(e,a,j,i)*G(e,a,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 30 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(30)=T_tens(1)

      ! 0.25 E11^{ba}_{mj}F11^{amb}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,a,m,j)+=X(b,c,k,m)*Y(c,a,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E11 wrong')
      ierr=t_einsum("G(b,a,m,j)+=X(a,m,i,j)*Y(b,i)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F11 wrong')
      ierr=t_einsum("R()+=T(b,a,m,j)*G(b,a,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 31 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(31)=T_tens(1)

      ! 0.25 E12^{ba}_{mi}F12^{amb}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(b,a,m,i)+=X(b,c,k,m)*Y(c,a,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E12 wrong')
      ierr=t_einsum("G(b,a,m,i)+=X(a,m,i,j)*Y(b,j)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F12 wrong')
      ierr=t_einsum("R()+=T(b,a,m,i)*G(b,a,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 32 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(32)=T_tens(1)

      ! 0.25 E13^{ca}_{mj}F13^{amc}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(c,a,m,j)+=X(b,c,k,m)*Y(a,b,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E13 wrong')
      ierr=t_einsum("G(c,a,m,j)+=X(a,m,i,j)*Y(c,i)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F13 wrong')
      ierr=t_einsum("R()+=T(c,a,m,j)*G(c,a,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 33 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(33)=T_tens(1)

      ! 0.25 E14^{ca}_{mi}F14^{amc}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims4,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims4,init_val=ZERO)
      ierr=t_einsum("T(c,a,m,i)+=X(b,c,k,m)*Y(a,b,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E14 wrong')
      ierr=t_einsum("G(c,a,m,i)+=X(a,m,i,j)*Y(c,j)",I2,vooo_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction F14 wrong')
      ierr=t_einsum("R()+=T(c,a,m,i)*G(c,a,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 34 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(34)=T_tens(1)

      dims2=nvir
      dims2(2)=nocc
      ! 0.25 A1^{}_{ej}B1^{e}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,j)+=X(a,b,i,e)*Y(a,b,i,j)",I1,vvov_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction A1 wrong')
      ierr=t_einsum("G(e,j)+=X(e,c,j,k)*Y(c,k)",I2,t2_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction B1 wrong')
      ierr=t_einsum("R()+=T(e,j)*G(e,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 35 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(35)=T_tens(1)

      ! 0.25 A2^{}_{ek}B2^{e}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(e,k)+=X(a,b,i,e)*Y(a,b,k,i)",I1,vvov_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction A2 wrong')
      ierr=t_einsum("G(e,k)+=X(e,c,j,k)*Y(c,j)",I2,t2_v,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction B2 wrong')
      ierr=t_einsum("R()+=T(e,k)*G(e,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 36 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(36)=T_tens(1)      

      dims2=nocc
      dims2(2)=nvir
      ! 0.25 A3^{mb}_{}B3^{b}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,b)+=X(a,m,i,j)*Y(a,b,i,j)",I1,vooo_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction A3 wrong')
      ierr=t_einsum("G(m,b)+=X(b,c,k,m)*Y(c,k)",I2,t2_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction B3 wrong')
      ierr=t_einsum("R()+=T(m,b)*G(m,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 37 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(37)=T_tens(1)

      ! 0.25 A4^{mc}_{}B4^{c}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims2,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims2,init_val=ZERO)
      ierr=t_einsum("T(m,c)+=X(a,m,i,j)*Y(c,a,i,j)",I1,vooo_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction A4 wrong')
      ierr=t_einsum("G(m,c)+=X(b,c,k,m)*Y(b,k)",I2,t2_o,fvo)
      if (ierr.ne.0) call quit('laplace_get_T: contraction B4 wrong')
      ierr=t_einsum("R()+=T(m,c)*G(m,c)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 38 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(38)=T_tens(1)


      ierr=destroy_tensor_t(t2)
      ierr=destroy_tensor_t(fvo)
      !------------------------------------------------------------------------------------------
      
      !clean up
      ierr=destroy_tensor_t(T_tensor)
      ierr=destroy_tensor_t(t2_o)
      ierr=destroy_tensor_t(t2_v)
      ierr=destroy_tensor_t(vvov_v)
      ierr=destroy_tensor_t(vooo_o)
      ierr=destroy_tensor_t(one_tripl)
      
    end subroutine laplace_get_T2

    subroutine laplace_get_T3(T,i_t2,i_t1,i_ooov,i_vovv,t2_dims, print_level)

      complex(8), intent(out)           :: T(:)
      type(tensor_t), intent(inout) :: i_t2,i_t1
      type(tensor_t), intent(inout) :: i_ooov,i_vovv
      integer(INTD), intent(in)         :: t2_dims(4)
      integer, intent(in)               :: print_level

      type(tensor_t)  :: t2_o, t2_v, t2, t1
      type(tensor_t)  :: vvov_v, vooo_o, vvoo, fvo
      type(tensor_t)  :: I1, I2, I3
      integer(INTD)       :: ierr
      integer(INTD)       :: dims_vvvv(4),dims_vvov(4),dims_vvoo(4),dims_vvvo(4)
      integer(INTD)       :: dims_oooo(4),dims_vooo(4),dims_vovo(4)
      integer(INTD)       :: dims_oo(2),dims_ov(2),dims_vo(2),dims_vv(2)
      integer             :: nvir, nocc

      ! scalars (need to be defined as tensor types)
      type(C_PTR)         :: body_p
      type(tensor_t)  :: T_tensor
      integer(INTD)       :: T_dims(1)
      complex(8), pointer :: T_tens(:)

      integer(INTD)      :: one_dims(1)
      type(tensor_t) :: one_tripl
      !construct one tensor
      ierr=construct_tensor_t(one_tripl,C8,one_dims(1:0),init_val=ONE)

      if (print_level>1) call print_date('Entered laplace_get_T3 routine')
      if (print_level>1) write(*,*) '-------------'

      !Initialize scalars that are to be used as tensors in contractions
      T_dims(1) = 1
      ierr=construct_tensor_t(T_tensor,C8,T_dims(1:0),init_val=ZERO)
      ierr=direct_data_pointer(T_tensor,body_p)
      call c_f_pointer(body_p,T_tens,T_dims)
     
      nvir = t2_dims(1)
      nocc = t2_dims(3)
      T_tens(1)=ZERO

      !-------------------- setting up indices
      dims_vv=nvir
      dims_vo=nvir
      dims_vo(2)=nocc

      dims_oo=nocc
      dims_ov=nocc
      dims_ov(2)=nvir
      
      dims_vvvv=nvir
      dims_vvov=nvir
      dims_vvov(3)=nocc
      dims_vvvo=nvir
      dims_vvvo(4)=nocc

      dims_oooo=nocc
      dims_vooo=nocc
      dims_vooo(1)=nvir
      

      dims_vvoo=nocc
      dims_vvoo(1:2)=nvir
      dims_vovo=nocc
      dims_vovo(1)=nvir
      dims_vovo(3)=nvir

      ! ------------------- compute first transformed integrals and amplitudes

      ierr=construct_tensor_t(t2_v,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(e,a,i,j)+=O+(e,a,i,j)*I()",t2_v,i_t2,one_tripl)

      ierr=construct_tensor_t(t2_o,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,m)+=O+(a,b,i,m)*I()",t2_o,i_t2,one_tripl)

      ierr=construct_tensor_t(vvov_v,C8,dims_vvov,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,e)+=O(e,i,b,a)*I()",vvov_v,i_vovv,one_tripl)

      ierr=construct_tensor_t(vooo_o,C8,dims_vooo,init_val=ZERO)
      ierr=t_einsum("R(a,m,i,j)+=O(j,i,m,a)*I()",vooo_o,i_ooov,one_tripl)

      ierr=construct_tensor_t(t1,C8,dims_vo,init_val=ZERO)
      ierr=t_einsum("R(a,i)+=O(a,i)",t1,i_t1)
      ! ----------------------------------------------------------------------------------------- compute terms with t1

      ! 0.25 M1^{b}_{e}N1^{e}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(b,e)+=X(a,b,i,e)*Y(a,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M1 wrong')
      ierr=t_einsum("G(e,j)+=X(e,c,j,k)*Y(c,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N1 wrong')
      ierr=t_einsum("K(b,e)+=G(e,j)*Y(b,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X1 wrong')
      ierr=t_einsum("R()+=T(b,e)*K(b,e)",T_tensor,I1,I3,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 1 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(1)=T_tens(1) 

      ! -0.25 M2^{a}_{e}N2^{e}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(a,e)+=X(a,b,i,e)*Y(b,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M2 wrong')
      ierr=t_einsum("G(e,j)+=X(e,c,j,k)*Y(c,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N2 wrong')
      ierr=t_einsum("K(a,e)+=G(e,j)*Y(a,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X2 wrong')
      ierr=t_einsum("R()+=T(a,e)*K(a,e)",T_tensor,I1,I3,scale=MINUS_ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 2 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(2)=T_tens(1)

      ! 0.25 M3^{a}_{e}N3^{e}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(a,e)+=X(a,b,i,e)*Y(b,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M3 wrong')
      ierr=t_einsum("G(e,k)+=X(e,c,j,k)*Y(c,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N3 wrong')
      ierr=t_einsum("K(a,e)+=G(e,k)*Y(a,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X3 wrong')
      ierr=t_einsum("R()+=T(a,e)*K(a,e)",T_tensor,I1,I3,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 3 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(3)=T_tens(1)

      ! -0.25 M4^{b}_{e}N4^{e}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(b,e)+=X(a,b,i,e)*Y(a,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M4 wrong')
      ierr=t_einsum("G(e,k)+=X(e,c,j,k)*Y(c,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N4 wrong')
      ierr=t_einsum("K(b,e)+=G(e,k)*Y(b,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X4 wrong')
      ierr=t_einsum("R()+=T(b,e)*K(b,e)",T_tensor,I1,I3,scale=MINUS_ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 4 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(4)=T_tens(1)

      ! 0.25 M4^{m}_{j}N4^{a}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,j)+=X(a,m,i,j)*Y(a,i)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M5 wrong')
      ierr=t_einsum("G(b,m)+=X(b,c,k,m)*Y(c,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N5 wrong')
      ierr=t_einsum("K(m,j)+=G(b,m)*Y(b,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X5 wrong')
      ierr=t_einsum("R()+=T(m,j)*K(m,j)",T_tensor,I1,I3,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 5 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(5)=T_tens(1)

      ! -0.25 M5^{m}_{i}N5^{b}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,i)+=X(a,m,i,j)*Y(a,j)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M5 wrong')
      ierr=t_einsum("G(b,m)+=X(b,c,k,m)*Y(c,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N5 wrong')
      ierr=t_einsum("K(m,i)+=G(b,m)*Y(b,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X5 wrong')
      ierr=t_einsum("R()+=T(m,i)*K(m,i)",T_tensor,I1,I3,scale=MINUS_ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 6 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(6)=T_tens(1)

      ! 0.25 M6^{m}_{i}N6^{c}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,i)+=X(a,m,i,j)*Y(a,j)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M6 wrong')
      ierr=t_einsum("G(c,m)+=X(b,c,k,m)*Y(b,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N6 wrong')
      ierr=t_einsum("K(m,i)+=G(c,m)*Y(c,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X6 wrong')
      ierr=t_einsum("R()+=T(m,i)*K(m,i)",T_tensor,I1,I3,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 7 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(7)=T_tens(1)

      ! 0.25 M7^{m}_{j}N7^{c}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,j)+=X(a,m,i,j)*Y(a,i)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M7 wrong')
      ierr=t_einsum("G(c,m)+=X(b,c,k,m)*Y(b,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N7 wrong')
      ierr=t_einsum("K(m,j)+=G(c,m)*Y(c,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X7 wrong')
      ierr=t_einsum("R()+=T(m,j)*K(m,j)",T_tensor,I1,I3,scale=MINUS_ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 8 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(8)=T_tens(1)

      ! 0.08 M8^{b}_{iej}N8^{ecb}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(b,i,e,j)+=X(a,b,i,e)*Y(a,j)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M8 wrong')
      ierr=t_einsum("G(e,c,b,j)+=X(e,c,j,k)*Y(b,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N8 wrong')
      ierr=t_einsum("K(b,i,e,j)+=G(e,c,b,j)*Y(c,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X8 wrong')
      ierr=t_einsum("R()+=T(b,i,e,j)*K(b,i,e,j)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 9 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(9)=T_tens(1)

      ! -0.08 M9^{b}_{iek}N9^{ecb}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(b,i,e,k)+=X(a,b,i,e)*Y(a,k)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M9 wrong')
      ierr=t_einsum("G(e,c,b,k)+=X(e,c,j,k)*Y(b,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N9 wrong')
      ierr=t_einsum("K(b,i,e,k)+=G(e,c,b,k)*Y(c,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X9 wrong')
      ierr=t_einsum("R()+=T(b,i,e,k)*K(b,i,e,k)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 10 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(10)=T_tens(1)

      ! 0.08 M10^{abc}_{e}N10^{eca}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvvv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vvvv,init_val=ZERO)
      ierr=t_einsum("T(a,b,c,e)+=X(a,b,i,e)*Y(c,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M10e wrong')
      ierr=t_einsum("G(e,c,a,k)+=X(e,c,j,k)*Y(a,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N10e wrong')
      ierr=t_einsum("K(a,b,c,e)+=G(e,c,a,k)*Y(b,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X10e wrong')
      ierr=t_einsum("R()+=T(a,b,c,e)*K(a,b,c,e)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 11 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(11)=T_tens(1)

      ! -0.08 M10^{amc}_{j}N10^{cba}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,m,c,j)+=X(a,m,i,j)*Y(c,i)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M10m wrong')
      ierr=t_einsum("G(b,c,a,m)+=X(b,c,k,m)*Y(a,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N10m wrong')
      ierr=t_einsum("K(a,m,c,j)+=G(b,c,a,m)*Y(b,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X10m wrong')
      ierr=t_einsum("R()+=T(a,m,c,j)*K(a,m,c,j)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 12 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(12)=T_tens(1)

      ! -0.08 M11^{abc}_{e}N11^{eca}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvvv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vvvv,init_val=ZERO)
      ierr=t_einsum("T(a,b,c,e)+=X(a,b,i,e)*Y(c,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M11 wrong')
      ierr=t_einsum("G(e,c,a,j)+=X(e,c,j,k)*Y(a,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N11 wrong')
      ierr=t_einsum("K(a,b,c,e)+=G(e,c,a,j)*Y(b,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X11 wrong')
      ierr=t_einsum("R()+=T(a,b,c,e)*K(a,b,c,e)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 13 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(13)=T_tens(1)

      ! -0.08 M12^{a}_{iej}N12^{e}_{jki}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,i,e,j)+=X(a,b,i,e)*Y(b,j)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M12 wrong')
      ierr=t_einsum("G(e,j,k,i)+=X(e,c,j,k)*Y(c,i)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N12 wrong')
      ierr=t_einsum("K(a,i,e,j)+=G(e,j,k,i)*Y(a,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X12 wrong')
      ierr=t_einsum("R()+=T(a,i,e,j)*K(a,i,e,j)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 14 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(14)=T_tens(1)

      ! 0.08 M12^{a}_{jei}N12^{e}_{kij}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,j,e,i)+=X(a,b,j,e)*Y(b,i)",I1,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M12a wrong')
      ierr=t_einsum("G(e,k,i,j)+=X(e,c,k,i)*Y(c,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N12a wrong')
      ierr=t_einsum("K(a,j,e,i)+=G(e,k,i,j)*Y(a,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X12a wrong')
      ierr=t_einsum("R()+=T(a,j,e,i)*K(a,j,e,i)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 15 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(15)=T_tens(1)

      ! -0.08 M13^{m}_{ijk}N13^{c}_{kmj}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oooo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oooo,init_val=ZERO)
      ierr=t_einsum("T(m,i,j,k)+=X(a,m,i,j)*Y(a,k)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M13a wrong')
      ierr=t_einsum("G(c,k,m,j)+=X(b,c,k,m)*Y(b,j)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N13a wrong')
      ierr=t_einsum("K(m,i,j,k)+=G(c,k,m,j)*Y(c,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X13a wrong')
      ierr=t_einsum("R()+=T(m,i,j,k)*K(m,i,j,k)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 16 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(16)=T_tens(1)

      ! 0.08 M13^{m}_{ijk}N13^{a}_{kmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oooo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_oooo,init_val=ZERO)
      ierr=t_einsum("T(m,i,j,k)+=X(a,m,i,j)*Y(a,k)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M13b wrong')
      ierr=t_einsum("G(c,k,m,i)+=X(b,c,k,m)*Y(b,i)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N13b wrong')
      ierr=t_einsum("K(m,i,j,k)+=G(c,k,m,i)*Y(c,j)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X13b wrong')
      ierr=t_einsum("R()+=T(m,i,j,k)*K(m,i,j,k)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 17 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(17)=T_tens(1)

      ! 0.08 M14^{amb}_{i}N14^{b}_{jmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,m,b,j)+=X(a,m,i,j)*Y(b,i)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M14 wrong')
      ierr=t_einsum("G(b,k,m,j)+=X(b,c,k,m)*Y(c,j)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N14 wrong')
      ierr=t_einsum("K(a,m,b,j)+=G(b,k,m,j)*Y(a,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X14 wrong')
      ierr=t_einsum("R()+=T(m,i,j,k)*K(m,i,j,k)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 18 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(18)=T_tens(1)

      ! -0.08 M15^{amb}_{i}N15^{b}_{kmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vooo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,m,b,i)+=X(a,m,i,j)*Y(b,j)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M15 wrong')
      ierr=t_einsum("G(b,k,m,i)+=X(b,c,k,m)*Y(c,i)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N15 wrong')
      ierr=t_einsum("K(a,m,b,i)+=G(b,k,m,i)*Y(a,k)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X15 wrong')
      ierr=t_einsum("R()+=T(a,m,b,i)*K(a,m,b,i)",T_tensor,I1,I3,scale=MINUS_ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 19 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(19)=T_tens(1)

      ! 0.08 M16^{amc}_{i}N16^{b}_{kmi}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vovo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvo,init_val=ZERO)
      ierr=construct_tensor_t(I3,C8,dims_vovo,init_val=ZERO)
      ierr=t_einsum("T(a,m,c,i)+=X(a,m,i,j)*Y(c,j)",I1,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction M16 wrong')
      ierr=t_einsum("G(b,c,a,m)+=X(b,c,k,m)*Y(a,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction N16 wrong')
      ierr=t_einsum("K(a,m,c,i)+=G(b,c,a,m)*Y(b,i)",I3,I2,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction X16 wrong')
      ierr=t_einsum("R()+=T(a,m,c,i)*K(a,m,c,i)",T_tensor,I1,I3,scale=ONE_TWELFTH)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 20 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      ierr=destroy_tensor_t(I3)
      T(20)=T_tens(1)


      ! ----------------------------------------------------------------------------------------- compute terms with t2
      ierr=construct_tensor_t(t2,C8,t2_dims,init_val=ZERO)
      ierr=t_einsum("R(a,b,i,j)+=O(a,b,i,j)",t2,i_t2)

      ! 0.25 K1^{eb}_{}L1^{b}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(e,b)+=X(e,c,j,k)*Y(b,c,j,k)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E1 wrong')
      ierr=t_einsum("G(e,b)+=X(a,b,i,e)*Y(a,i)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R1 wrong')
      ierr=t_einsum("R()+=T(e,b)*G(e,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 21 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(21)=T_tens(1)  

      ! 0.25 K2^{ea}_{}L2^{a}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vv,init_val=ZERO)
      ierr=t_einsum("T(e,a)+=X(e,c,j,k)*Y(c,a,j,k)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E2 wrong')
      ierr=t_einsum("G(e,a)+=X(a,b,i,e)*Y(b,i)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R2 wrong')
      ierr=t_einsum("R()+=T(e,a)*G(e,a)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 22 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(22)=T_tens(1)  

      ! 0.25 K3^{}_{mj}L3^{m}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,j)+=X(b,c,k,m)*Y(b,c,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E3 wrong')
      ierr=t_einsum("G(m,j)+=X(a,m,i,j)*Y(a,i)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R3 wrong')
      ierr=t_einsum("R()+=T(m,j)*G(m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 23 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(23)=T_tens(1)  

      ! 0.25 K4^{}_{mi}L4^{m}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_oo,init_val=ZERO)
      ierr=t_einsum("T(m,i)+=X(b,c,k,m)*Y(b,c,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E4 wrong')
      ierr=t_einsum("G(m,i)+=X(a,m,i,j)*Y(a,j)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R4 wrong')
      ierr=t_einsum("R()+=T(m,i)*G(m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 24 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(24)=T_tens(1)  

      ! 0.25 K5^{}_{kmij}L5^{m}_{ijk}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_oooo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_oooo,init_val=ZERO)
      ierr=t_einsum("T(k,m,i,j)+=X(b,c,k,m)*Y(b,c,i,j)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E5 wrong')
      ierr=t_einsum("G(k,m,i,j)+=X(a,m,i,j)*Y(a,k)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R5 wrong')
      ierr=t_einsum("R()+=T(k,m,i,j)*G(k,m,i,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 25 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(25)=T_tens(1)  

      ! 0.25 K6^{ecab}_{}L6^{abc}_{e}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvvv,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvvv,init_val=ZERO)
      ierr=t_einsum("T(e,c,a,b)+=X(e,c,j,k)*Y(a,b,j,k)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E6 wrong')
      ierr=t_einsum("G(e,c,a,b)+=X(a,b,i,e)*Y(c,i)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R6 wrong')
      ierr=t_einsum("R()+=T(e,c,a,b)*G(e,c,a,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 26 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(26)=T_tens(1)      

      ! 0.25 K7^{eb}_{ki}L7^{b}_{iek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(e,b,k,i)+=X(e,c,j,k)*Y(b,c,i,j)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E7 wrong')
      ierr=t_einsum("G(e,b,k,i)+=X(a,b,i,e)*Y(a,k)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R7 wrong')
      ierr=t_einsum("R()+=T(e,b,k,i)*G(e,b,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 27 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(27)=T_tens(1)

      ! 0.25 K8^{eb}_{ji}L8^{b}_{iej}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(e,b,j,i)+=X(e,c,j,k)*Y(b,c,k,i)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E8 wrong')
      ierr=t_einsum("G(e,b,j,i)+=X(a,b,i,e)*Y(a,j)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R8 wrong')
      ierr=t_einsum("R()+=T(e,b,j,i)*G(e,b,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 28 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(28)=T_tens(1)

      ! 0.25 K9^{ea}_{ki}L9^{a}_{iek}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(e,a,k,i)+=X(e,c,j,k)*Y(c,a,i,j)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E9 wrong')
      ierr=t_einsum("G(e,a,k,i)+=X(a,b,i,e)*Y(b,k)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R9 wrong')
      ierr=t_einsum("R()+=T(e,a,k,i)*G(e,a,k,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 29 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(29)=T_tens(1)

      ! 0.25 K10^{ea}_{ji}L10^{a}_{iej}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(e,a,j,i)+=X(e,c,j,k)*Y(c,a,k,i)",I1,t2_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E10 wrong')
      ierr=t_einsum("G(e,a,j,i)+=X(a,b,i,e)*Y(b,j)",I2,vvov_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R10 wrong')
      ierr=t_einsum("R()+=T(e,a,j,i)*G(e,a,j,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 30 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(30)=T_tens(1)

      ! 0.25 K11^{ba}_{mj}L11^{amb}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(b,a,m,j)+=X(b,c,k,m)*Y(c,a,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E11 wrong')
      ierr=t_einsum("G(b,a,m,j)+=X(a,m,i,j)*Y(b,i)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R11 wrong')
      ierr=t_einsum("R()+=T(b,a,m,j)*G(b,a,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 31 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(31)=T_tens(1)

      ! 0.25 K12^{ba}_{mi}L12^{amb}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(b,a,m,i)+=X(b,c,k,m)*Y(c,a,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E12 wrong')
      ierr=t_einsum("G(b,a,m,i)+=X(a,m,i,j)*Y(b,j)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R12 wrong')
      ierr=t_einsum("R()+=T(b,a,m,i)*G(b,a,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 32 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(32)=T_tens(1)

      ! 0.25 K13^{ca}_{mj}L13^{amc}_{j}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(c,a,m,j)+=X(b,c,k,m)*Y(a,b,j,k)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E13 wrong')
      ierr=t_einsum("G(c,a,m,j)+=X(a,m,i,j)*Y(c,i)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R13 wrong')
      ierr=t_einsum("R()+=T(c,a,m,j)*G(c,a,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 33 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(33)=T_tens(1)

      ! 0.25 E14^{ca}_{mi}R14^{amc}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vvoo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vvoo,init_val=ZERO)
      ierr=t_einsum("T(c,a,m,i)+=X(b,c,k,m)*Y(a,b,k,i)",I1,t2_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E14 wrong')
      ierr=t_einsum("G(c,a,m,i)+=X(a,m,i,j)*Y(c,j)",I2,vooo_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R14 wrong')
      ierr=t_einsum("R()+=T(c,a,m,i)*G(c,a,m,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 34 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(34)=T_tens(1)

      ! 0.25 P1^{}_{ej}Q1^{e}_{i}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=t_einsum("T(e,j)+=X(a,b,i,e)*Y(a,b,i,j)",I1,vvov_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E15 wrong')
      ierr=t_einsum("G(e,j)+=X(e,c,j,k)*Y(c,k)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R15 wrong')
      ierr=t_einsum("R()+=T(e,j)*G(e,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 35 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(35)=T_tens(1)

      ! 0.25 P2^{}_{ek}Q2^{e}_{k}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_vo,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_vo,init_val=ZERO)
      ierr=t_einsum("T(e,k)+=X(a,b,i,e)*Y(a,b,k,i)",I1,vvov_v,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E16 wrong')
      ierr=t_einsum("G(e,k)+=X(e,c,j,k)*Y(c,j)",I2,t2_v,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R16 wrong')
      ierr=t_einsum("R()+=T(e,k)*G(e,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 36 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(36)=T_tens(1) 

      ! 0.25 P3^{mb}_{}Q3^{b}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_ov,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_ov,init_val=ZERO)
      ierr=t_einsum("T(m,b)+=X(a,m,i,j)*Y(a,b,i,j)",I1,vooo_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E17 wrong')
      ierr=t_einsum("G(m,b)+=X(b,c,k,m)*Y(c,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R17 wrong')
      ierr=t_einsum("R()+=T(m,b)*G(m,b)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 37 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(37)=T_tens(1)

      ! 0.25 P4^{mc}_{}Q4^{c}_{m}
      T_tens(1)=ZERO
      ierr=construct_tensor_t(I1,C8,dims_ov,init_val=ZERO)
      ierr=construct_tensor_t(I2,C8,dims_ov,init_val=ZERO)
      ierr=t_einsum("T(m,c)+=X(a,m,i,j)*Y(c,a,i,j)",I1,vooo_o,t2)
      if (ierr.ne.0) call quit('laplace_get_T: contraction E18 wrong')
      ierr=t_einsum("G(m,c)+=X(b,c,k,m)*Y(b,k)",I2,t2_o,t1)
      if (ierr.ne.0) call quit('laplace_get_T: contraction R18 wrong')
      ierr=t_einsum("R()+=T(m,c)*G(m,c)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 38 wrong')
      ierr=destroy_tensor_t(I1)
      ierr=destroy_tensor_t(I2)
      T(38)=T_tens(1) 

      ierr=destroy_tensor_t(t2)
      !------------------------------------------------------------------------------------------
      
      !clean up
      ierr=destroy_tensor_t(t1)
      ierr=destroy_tensor_t(T_tensor)
      ierr=destroy_tensor_t(t2_o)
      ierr=destroy_tensor_t(t2_v)
      ierr=destroy_tensor_t(vvov_v)
      ierr=destroy_tensor_t(vooo_o)
      ierr=destroy_tensor_t(one_tripl)

    end subroutine laplace_get_T3

end module eti_triples
