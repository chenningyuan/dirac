module eti_diis
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tapp_top
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_cfg, only: exacc_input
        use exacorr_utils
        use intermediates_cc
        use intermediates_lambda
        use intermediates_eom_ccsd
        use eti_common_routines
        use mfsolver_dirac_input
        use exacorr_printout
        use exacorr_scaling_parameters
        use, intrinsic:: ISO_C_BINDING

        implicit none

        type(tensor_t) :: one_tensor

        private

        public diis

!        interface scale_with_denominators
!           module procedure scale_with_denominators
!        end interface
!
        interface diis
           module procedure diis
        end interface

       contains
!-------------------------------------------------------------------------------------------

       subroutine diis (t5_tensor,t4_tensor,t3_tensor,t2_tensor,t1_tensor,t_convergence,mode)
         logical,parameter  :: crop=.true. ! if switching off, also increase dimension of maxdim (used to be 5)
         integer, parameter :: maxdim = 3 ! maximum number of tensors that are stored
         integer, parameter :: mindim = 2 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tensor_t), intent(inout)  :: t2_tensor
         type(tensor_t), intent(inout), optional  :: t1_tensor
         type(tensor_t), intent(inout), optional  :: t3_tensor
         type(tensor_t), intent(inout), optional  :: t4_tensor
         type(tensor_t), intent(inout), optional  :: t5_tensor
         integer             :: i
         integer(INTD)       :: ierr, tens_rank
         integer(INTD), allocatable, save :: t1_dims(:)
         integer(INTD), allocatable, save :: t2_dims(:)
         integer(INTD), allocatable, save :: t3_dims(:)
         integer(INTD), allocatable, save :: t4_dims(:)
         integer(INTD), allocatable, save :: t5_dims(:)
         logical             :: init, clean_up, CCD
         logical             :: CCSDitT, CCSDitTQ, CCSDitTQP
         integer, save       :: itens, ntens
         complex(8), save    :: b_diis(maxdim,maxdim)
         complex(8)          :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         real(8)             :: b_max, b_scale
         type(tensor_t), save :: t1_inp_tensor
         type(tensor_t), save :: t2_inp_tensor
         type(tensor_t), save :: t3_inp_tensor
         type(tensor_t), save :: t4_inp_tensor
         type(tensor_t), save :: t5_inp_tensor
         type(tensor_t), save :: t1_res_tensor(maxdim), t1_err_tensor(maxdim)
         type(tensor_t), save :: t2_res_tensor(maxdim), t2_err_tensor(maxdim)
         type(tensor_t), save :: t3_res_tensor(maxdim), t3_err_tensor(maxdim)
         type(tensor_t), save :: t4_res_tensor(maxdim), t4_err_tensor(maxdim)
         type(tensor_t), save :: t5_res_tensor(maxdim), t5_err_tensor(maxdim)
         integer             :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)             :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tensor_t)  :: scalar_tensor
         integer(INTD)       :: scalar_dims(1)
         complex(8), pointer :: scalar_tens(:)
         type(C_PTR):: body_p

         complex(4), pointer :: scalar_tens_c8(:)
         real(8), pointer :: scalar_tens_r8(:)
         real(4), pointer :: scalar_tens_r4(:)

         integer(kind(TAPP_datatype)) :: eti_storage_datatype

         init       = .false.
         clean_up   = .false.
         CCD = .false.
         CCSDitT = .false.
         CCSDitTQ = .false.
         CCSDitTQP = .false.

         eti_storage_datatype =  t2_tensor%datatype

         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.
         if (present(t3_tensor)) CCSDitT = .true.
         if (present(t4_tensor)) CCSDitTQ = .true.
         if (present(t5_tensor)) CCSDitTQP = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                ierr = get_nmodes(t1_tensor,tens_rank)
                ierr = get_extents(t1_tensor,t1_dims)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'

                if (CCSDitT) then 
                   ierr = get_nmodes(t3_tensor,tens_rank)
                   ierr = get_extents(t3_tensor,t3_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.6) stop 'input t3-tensor corrupted in diis'
                end if

                if (CCSDitTQ) then 
                   ierr = get_nmodes(t4_tensor,tens_rank)
                   ierr = get_extents(t4_tensor,t4_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.8) stop 'input t4-tensor corrupted in diis'
                end if

                if (CCSDitTQP) then
                   ierr = get_nmodes(t5_tensor,tens_rank)
                   ierr = get_extents(t5_tensor,t5_dims) 
                   if (ierr.ne.0 .or. tens_rank.ne.10) stop 'input t5-tensor corrupted in diis'
                end if
            end if
            ierr = get_nmodes(t2_tensor,tens_rank)
            ierr = get_extents(t2_tensor,t2_dims) 
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis '
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) then
               ierr=construct_tensor_t(t1_inp_tensor,eti_storage_datatype,t1_dims,init_val=ZERO)
               if (CCSDitT) ierr=construct_tensor_t(t3_inp_tensor,eti_storage_datatype,t3_dims,init_val=ZERO)
               if (CCSDitTQ) ierr=construct_tensor_t(t4_inp_tensor,eti_storage_datatype,t4_dims,init_val=ZERO)
               if (CCSDitTQP) ierr=construct_tensor_t(t5_inp_tensor,eti_storage_datatype,t5_dims,init_val=ZERO)
            end if
            ierr=construct_tensor_t(t2_inp_tensor,eti_storage_datatype,t2_dims,init_val=ZERO)
            do i = 1, maxdim
               if (.not.CCD) then
                    ierr=construct_tensor_t(t1_res_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    ierr=construct_tensor_t(t1_err_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    if (CCSDitT) then
                       ierr=construct_tensor_t(t3_res_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t3_err_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQ) then
                       ierr=construct_tensor_t(t4_res_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t4_err_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQP) then
                       ierr=construct_tensor_t(t5_res_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t5_err_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                    end if
               end if
               ierr=construct_tensor_t(t2_res_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               ierr=construct_tensor_t(t2_err_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               if (ierr.ne.0) stop 'not enough space for diis'
            end do
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=initialize_tensor_t(t1_inp_tensor)
                ierr=t_einsum("T1O(a,i)+=T1N(a,i)",t1_inp_tensor,t1_tensor)
                if (CCSDitT) then
                   ierr=initialize_tensor_t(t3_inp_tensor)
                   ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3N(a,b,c,i,j,k)",t3_inp_tensor,t3_tensor)
                end if
                if (CCSDitTQ) then
                   ierr=initialize_tensor_t(t4_inp_tensor)
                   ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4N(a,b,c,d,i,j,k,l)",t4_inp_tensor,t4_tensor)
                end if
                if (CCSDitTQP) then
                   ierr=initialize_tensor_t(t5_inp_tensor)
                   ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5N(a,b,c,d,e,i,j,k,l,m)",t5_inp_tensor,t5_tensor)
                end if
            end if
            ierr=initialize_tensor_t(t2_inp_tensor)
            ierr=t_einsum("T2O(a,b,i,j)+=T2N(a,b,i,j)",t2_inp_tensor,t2_tensor)
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) then
               ierr=destroy_tensor_t(t1_inp_tensor)
               if (CCSDitT) ierr=destroy_tensor_t(t3_inp_tensor)
               if (CCSDitTQ) ierr=destroy_tensor_t(t4_inp_tensor)
               if (CCSDitTQP) ierr=destroy_tensor_t(t5_inp_tensor)
            end if
            ierr=destroy_tensor_t(t2_inp_tensor)
            do i = 1, maxdim
               if (.not.CCD) then
                  ierr=destroy_tensor_t(t1_res_tensor(i))
                  ierr=destroy_tensor_t(t1_err_tensor(i))
                  if (CCSDitT) then
                     ierr=destroy_tensor_t(t3_res_tensor(i))
                     ierr=destroy_tensor_t(t3_err_tensor(i))
                  end if
                  if (CCSDitTQ) then
                     ierr=destroy_tensor_t(t4_res_tensor(i))
                     ierr=destroy_tensor_t(t4_err_tensor(i))
                  end if
                  if (CCSDitTQP) then
                     ierr=destroy_tensor_t(t5_res_tensor(i))
                     ierr=destroy_tensor_t(t5_err_tensor(i))
                  end if
               end if
               ierr=destroy_tensor_t(t2_res_tensor(i))
               ierr=destroy_tensor_t(t2_err_tensor(i))
            end do
            return
         end if

!        define scalar tensor needed to extract the error
         scalar_dims(1) = 1
         ierr=construct_tensor_t(scalar_tensor,eti_storage_datatype,scalar_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(scalar_tensor,body_p)

         if (eti_storage_datatype.eq.TAPP_F32) then
            call c_f_pointer(body_p,scalar_tens_r4,scalar_dims)
         else if (eti_storage_datatype.eq.TAPP_F64) then
            call c_f_pointer(body_p,scalar_tens_r8,scalar_dims)
         else if (eti_storage_datatype.eq.TAPP_C32) then
            call c_f_pointer(body_p,scalar_tens_c8,scalar_dims)
         else if (eti_storage_datatype.eq.TAPP_C64) then
            call c_f_pointer(body_p,scalar_tens,scalar_dims)
         else
            call quit ('unsupported storage datatype in tensor_norm2')
         end if

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_res_tensor(itens))
            ierr=t_einsum("T1R(a,i)+=T1(a,i)",t1_res_tensor(itens),t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_res_tensor(itens))
               ierr=t_einsum("T3R(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",t3_res_tensor(itens),t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_res_tensor(itens))
               ierr=t_einsum("T4R(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",t4_res_tensor(itens),t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_res_tensor(itens))
               ierr=t_einsum("T5R(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",t5_res_tensor(itens),t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_res_tensor(itens))
         ierr=t_einsum("T2R(a,b,i,j)+=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not. CCD) then
            ierr=initialize_tensor_t(t1_err_tensor(itens))
            ierr=t_einsum("T1E(a,i)+=T1R(a,i)",t1_err_tensor(itens),t1_tensor)
            ierr=t_einsum("T1E(a,i)+=T1I(a,i)",t1_err_tensor(itens),t1_inp_tensor,scale=MINUS_ONE)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_err_tensor(itens))
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3R(a,b,c,i,j,k)",t3_err_tensor(itens),t3_tensor)
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3I(a,b,c,i,j,k)",t3_err_tensor(itens),t3_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_err_tensor(itens))
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4R(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_tensor)
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4I(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_err_tensor(itens))
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5R(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_tensor)
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5I(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_inp_tensor,scale=MINUS_ONE)
            end if
         end if
         ierr=initialize_tensor_t(t2_err_tensor(itens))
         ierr=t_einsum("T2E(a,b,i,j)+=T2R(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
         ierr=t_einsum("T2E(a,b,i,j)+=T2I(a,b,i,j)",t2_err_tensor(itens),t2_inp_tensor,scale=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         ierr=initialize_tensor_t(scalar_tensor)
         do i = 1, ntens
            ierr=initialize_tensor_t(scalar_tensor)
            if (.not.CCD) then
               ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
               if (CCSDitT) then
                  ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)", &
                &     scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
               end if
               if (CCSDitTQ) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)", &
                &     scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
               end if
               if (CCSDitTQP) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)", &
                &     scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
               end if
            end if
            ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
! aspg, this association has already been done above             
!           ierr=direct_data_pointer(scalar_tensor,body_p)
!           call c_f_pointer(body_p,scalar_tens,scalar_dims)
            if (eti_storage_datatype.eq.TAPP_F32) then
               b_diis(itens,i) = cmplx(scalar_tens_r4(1),0.0D0)
               b_diis(i,itens) = b_diis(itens,i)
            else if (eti_storage_datatype.eq.TAPP_F64) then
               b_diis(itens,i) = cmplx(scalar_tens_r8(1),0.0D0)
               b_diis(i,itens) = b_diis(itens,i)
            else if (eti_storage_datatype.eq.TAPP_C32) then
               b_diis(itens,i) = cmplx(scalar_tens_c8(1)%re, scalar_tens_c8(1)%im)
               b_diis(i,itens) = cmplx(scalar_tens_c8(1)%re,-scalar_tens_c8(1)%im)
            else if (eti_storage_datatype.eq.TAPP_C64) then
               b_diis(itens,i) = scalar_tens(1)
               b_diis(i,itens) = dconjg(scalar_tens(1))
            else
               call quit ('unsupported storage datatype in diis()')
            end if
!           b_diis(itens,i) = scalar_tens(1)
!           b_diis(i,itens) = dconjg(scalar_tens(1))
            ! print *,'b_diis(',i,',',itens,'):',b_diis(i,itens)
! aspg, generalizing

         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             write(*,*) "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
         !   print *,'coeff:',coeff
           if (crop) then
              ! replace latest error vector by the optimized error vector using t1_tensor and t2_tensor as scratch space
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if 
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_err_tensor(i),scale=coeff(i))
                     if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)",&
                  &         t3_tensor,t3_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                  &         t4_tensor,t4_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                  &         t5_tensor,t5_err_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_err_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=copy_tensor_t(t1_err_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=copy_tensor_t(t3_err_tensor(itens),t3_tensor) 
                 if (CCSDitTQ) ierr=copy_tensor_t(t4_err_tensor(itens),t4_tensor) 
                 if (CCSDitTQP) ierr=copy_tensor_t(t5_err_tensor(itens),t5_tensor) 
              end if
              ierr=copy_tensor_t(t2_err_tensor(itens),t2_tensor)
              ! correct the latest row of the b_diis array for the fact that we use an optimized error vector
              ierr=initialize_tensor_t(scalar_tensor)
              do i = 1, ntens
                 ierr=initialize_tensor_t(scalar_tensor)
                 if (.not.CCD) then
                    ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
                    if (CCSDitT) ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)",&
                   &                                        scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
                    if (CCSDitTQ) ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)",&
                   &                                        scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
                    if (CCSDitTQP) ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)",&
                   &                                        scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
                 end if
                 ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))

!                ierr=direct_data_pointer(scalar_tensor,body_p)
!                call c_f_pointer(body_p,scalar_tens,scalar_dims)
!                b_diis(itens,i) = scalar_tens(1)
!                b_diis(i,itens) = dconjg(scalar_tens(1))
            if (eti_storage_datatype.eq.TAPP_F32) then
               b_diis(itens,i) = cmplx(scalar_tens_r4(1),0.0D0)
               b_diis(i,itens) = b_diis(itens,i)
            else if (eti_storage_datatype.eq.TAPP_F64) then
               b_diis(itens,i) = cmplx(scalar_tens_r8(1),0.0D0)
               b_diis(i,itens) = b_diis(itens,i)
            else if (eti_storage_datatype.eq.TAPP_C32) then
               b_diis(itens,i) = cmplx(scalar_tens_c8(1)%re, scalar_tens_c8(1)%im)
               b_diis(i,itens) = cmplx(scalar_tens_c8(1)%re,-scalar_tens_c8(1)%im)
            else if (eti_storage_datatype.eq.TAPP_C64) then
               b_diis(itens,i) = scalar_tens(1)
               b_diis(i,itens) = dconjg(scalar_tens(1))
            else
               call quit ('unsupported storage datatype in diis()')
            end if

              end do
              ! replace latest result vector by the optimized result vector, this is also the tensor to be used in the next iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=copy_tensor_t(t1_res_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=copy_tensor_t(t3_res_tensor(itens),t3_tensor)
                 if (CCSDitTQ) ierr=copy_tensor_t(t4_res_tensor(itens),t4_tensor)
                 if (CCSDitTQP) ierr=copy_tensor_t(t5_res_tensor(itens),t5_tensor)
              end if
              ierr=copy_tensor_t(t2_res_tensor(itens),t2_tensor)
           else ! not using crop
              ! make the amplitudes for the new iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
           end if ! end of crop if 
         end if

!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_inp_tensor)
            ierr=t_einsum("T1O(a,i)+=T1(a,i)",t1_inp_tensor,t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_inp_tensor)
               ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",&
                                      t3_inp_tensor,t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_inp_tensor)
               ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",&
                                      t4_inp_tensor,t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_inp_tensor)
               ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",&
              &                       t5_inp_tensor,t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_inp_tensor)
         ierr=t_einsum("T2O(a,b,i,j)+=T2(a,b,i,j)",t2_inp_tensor,t2_tensor)

! ierr=deallocate_temporary_full_data_tensor_t(scalar_tensor) ! todo like in ctf branch
         ierr=destroy_tensor_t(scalar_tensor)
       end subroutine diis

#if 0
       subroutine diis (t5_tensor,t4_tensor,t3_tensor,t2_tensor,t1_tensor,t_convergence,mode)
         logical,parameter  :: crop=.true. ! if switching off, also increase dimension of maxdim (used to be 5)
         integer, parameter :: maxdim = 3 ! maximum number of tensors that are stored
         integer, parameter :: mindim = 2 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tensor_t), intent(inout)  :: t2_tensor
         type(tensor_t), intent(inout), optional  :: t1_tensor
         type(tensor_t), intent(inout), optional  :: t3_tensor
         type(tensor_t), intent(inout), optional  :: t4_tensor
         type(tensor_t), intent(inout), optional  :: t5_tensor
         integer             :: i
         integer(INTD)       :: ierr, tens_rank
         integer(INTD), save :: t1_dims(1:2)
         integer(INTD), save :: t2_dims(1:4)
         integer(INTD), save :: t3_dims(1:6)
         integer(INTD), save :: t4_dims(1:8)
         integer(INTD), save :: t5_dims(1:10)
         logical             :: init, clean_up, CCD
         logical             :: CCSDitT, CCSDitTQ, CCSDitTQP
         integer, save       :: itens, ntens
         complex(8), save    :: b_diis(maxdim,maxdim)
         complex(8)          :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         real(8)             :: b_max, b_scale
         type(tensor_t), save :: t1_inp_tensor
         type(tensor_t), save :: t2_inp_tensor
         type(tensor_t), save :: t3_inp_tensor
         type(tensor_t), save :: t4_inp_tensor
         type(tensor_t), save :: t5_inp_tensor
         type(tensor_t), save :: t1_res_tensor(maxdim), t1_err_tensor(maxdim)
         type(tensor_t), save :: t2_res_tensor(maxdim), t2_err_tensor(maxdim)
         type(tensor_t), save :: t3_res_tensor(maxdim), t3_err_tensor(maxdim)
         type(tensor_t), save :: t4_res_tensor(maxdim), t4_err_tensor(maxdim)
         type(tensor_t), save :: t5_res_tensor(maxdim), t5_err_tensor(maxdim)
         integer             :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)             :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tensor_t)  :: scalar_tensor
         integer(INTD)       :: scalar_dims(1)
         complex(8), pointer :: scalar_tens(:)
         type(C_PTR):: body_p

         integer(kind(TAPP_datatype)) :: eti_storage_datatype

         init       = .false.
         clean_up   = .false.
         CCD = .false.
         CCSDitT = .false.
         CCSDitTQ = .false.
         CCSDitTQP = .false.

         eti_storage_datatype =  t2_tensor%datatype


         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.
         if (present(t3_tensor)) CCSDitT = .true.
         if (present(t4_tensor)) CCSDitTQ = .true.
         if (present(t5_tensor)) CCSDitTQP = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'

                if (CCSDitT) then
                   ierr = talsh_tensor_dimensions(t3_tensor,tens_rank,t3_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.6) stop 'input t3-tensor corrupted in diis'
                end if

                if (CCSDitTQ) then
                   ierr = talsh_tensor_dimensions(t4_tensor,tens_rank,t4_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.8) stop 'input t4-tensor corrupted in diis'
                end if

                if (CCSDitTQP) then
                   ierr = talsh_tensor_dimensions(t5_tensor,tens_rank,t5_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.10) stop 'input t5-tensor corrupted in diis'
                end if
            end if
            ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis '
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) then
               ierr=construct_tensor_t(t1_inp_tensor,eti_storage_datatype,t1_dims,init_val=ZERO)
               if (CCSDitT) ierr=construct_tensor_t(t3_inp_tensor,eti_storage_datatype,t3_dims,init_val=ZERO)
               if (CCSDitTQ) ierr=construct_tensor_t(t4_inp_tensor,eti_storage_datatype,t4_dims,init_val=ZERO)
               if (CCSDitTQP) ierr=construct_tensor_t(t5_inp_tensor,eti_storage_datatype,t5_dims,init_val=ZERO)
            end if
            ierr=construct_tensor_t(t2_inp_tensor,eti_storage_datatype,t2_dims,init_val=ZERO)
            do i = 1, maxdim
               if (.not.CCD) then
                    ierr=construct_tensor_t(t1_res_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    ierr=construct_tensor_t(t1_err_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    if (CCSDitT) then
                       ierr=construct_tensor_t(t3_res_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t3_err_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQ) then
                       ierr=construct_tensor_t(t4_res_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t4_err_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQP) then
                       ierr=construct_tensor_t(t5_res_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t5_err_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                    end if
               end if
               ierr=construct_tensor_t(t2_res_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               ierr=construct_tensor_t(t2_err_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               if (ierr.ne.0) stop 'not enough space for diis'
            end do
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=initialize_tensor_t(t1_inp_tensor)
                ierr=t_einsum("T1O(a,i)+=T1N(a,i)",t1_inp_tensor,t1_tensor)
                if (CCSDitT) then
                   ierr=initialize_tensor_t(t3_inp_tensor)
                   ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3N(a,b,c,i,j,k)",t3_inp_tensor,t3_tensor)
                end if
                if (CCSDitTQ) then
                   ierr=initialize_tensor_t(t4_inp_tensor)
                   ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4N(a,b,c,d,i,j,k,l)",t4_inp_tensor,t4_tensor)
                end if
                if (CCSDitTQP) then
                   ierr=initialize_tensor_t(t5_inp_tensor)
                   ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5N(a,b,c,d,e,i,j,k,l,m)",t5_inp_tensor,t5_tensor)
                end if
            end if
            ierr=initialize_tensor_t(t2_inp_tensor)
            ierr=t_einsum("T2O(a,b,i,j)+=T2N(a,b,i,j)",t2_inp_tensor,t2_tensor)
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) then
               ierr=destroy_tensor_t(t1_inp_tensor)
               if (CCSDitT) ierr=destroy_tensor_t(t3_inp_tensor)
               if (CCSDitTQ) ierr=destroy_tensor_t(t4_inp_tensor)
               if (CCSDitTQP) ierr=destroy_tensor_t(t5_inp_tensor)
            end if
            ierr=destroy_tensor_t(t2_inp_tensor)
            do i = 1, maxdim
               if (.not.CCD) then
                  ierr=destroy_tensor_t(t1_res_tensor(i))
                  ierr=destroy_tensor_t(t1_err_tensor(i))
                  if (CCSDitT) then
                     ierr=destroy_tensor_t(t3_res_tensor(i))
                     ierr=destroy_tensor_t(t3_err_tensor(i))
                  end if
                  if (CCSDitTQ) then
                     ierr=destroy_tensor_t(t4_res_tensor(i))
                     ierr=destroy_tensor_t(t4_err_tensor(i))
                  end if
                  if (CCSDitTQP) then
                     ierr=destroy_tensor_t(t5_res_tensor(i))
                     ierr=destroy_tensor_t(t5_err_tensor(i))
                  end if
               end if
               ierr=destroy_tensor_t(t2_res_tensor(i))
               ierr=destroy_tensor_t(t2_err_tensor(i))
            end do
            return
         end if

!        define scalar tensor needed to extract the error
         scalar_dims(1) = 1
         ierr=construct_tensor_t(scalar_tensor,eti_storage_datatype,scalar_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(scalar_tensor,body_p)
         call c_f_pointer(body_p,scalar_tens,scalar_dims)

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_res_tensor(itens))
            ierr=t_einsum("T1R(a,i)+=T1(a,i)",t1_res_tensor(itens),t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_res_tensor(itens))
               ierr=t_einsum("T3R(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",t3_res_tensor(itens),t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_res_tensor(itens))
               ierr=t_einsum("T4R(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",t4_res_tensor(itens),t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_res_tensor(itens))
               ierr=t_einsum("T5R(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",t5_res_tensor(itens),t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_res_tensor(itens))
         ierr=t_einsum("T2R(a,b,i,j)+=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not. CCD) then
            ierr=initialize_tensor_t(t1_err_tensor(itens))
            ierr=t_einsum("T1E(a,i)+=T1R(a,i)",t1_err_tensor(itens),t1_tensor)
            ierr=t_einsum("T1E(a,i)+=T1I(a,i)",t1_err_tensor(itens),t1_inp_tensor,scale=MINUS_ONE)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_err_tensor(itens))
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3R(a,b,c,i,j,k)",t3_err_tensor(itens),t3_tensor)
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3I(a,b,c,i,j,k)",t3_err_tensor(itens),t3_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_err_tensor(itens))
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4R(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_tensor)
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4I(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_err_tensor(itens))
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5R(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_tensor)
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5I(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_inp_tensor,scale=MINUS_ONE)
            end if
         end if
         ierr=initialize_tensor_t(t2_err_tensor(itens))
         ierr=t_einsum("T2E(a,b,i,j)+=T2R(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
         ierr=t_einsum("T2E(a,b,i,j)+=T2I(a,b,i,j)",t2_err_tensor(itens),t2_inp_tensor,scale=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         ierr=initialize_tensor_t(scalar_tensor)
         do i = 1, ntens
            ierr=initialize_tensor_t(scalar_tensor)
            if (.not.CCD) then
               ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
               if (CCSDitT) then
                  ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)", &
                &     scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
               end if
               if (CCSDitTQ) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)", &
                &     scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
               end if
               if (CCSDitTQP) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)", &
                &     scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
               end if
            end if
            ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
            b_diis(itens,i) = scalar_tens(1)
            b_diis(i,itens) = dconjg(scalar_tens(1))
            ! print *,'b_diis(',i,',',itens,'):',b_diis(i,itens)
         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             write(*,*) "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
         !   print *,'coeff:',coeff
           if (crop) then
              ! replace latest error vector by the optimized error vector using t1_tensor and t2_tensor as scratch space
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if 
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_err_tensor(i),scale=coeff(i))
                     if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)",&
                  &         t3_tensor,t3_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                  &         t4_tensor,t4_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                  &         t5_tensor,t5_err_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_err_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=t_einsum("T1E(a,i)=T1(a,i)",t1_err_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=t_einsum("T3E(a,b,c,i,j,k)=T3(a,b,c,i,j,k)",&
                                                           t3_err_tensor(itens),t3_tensor) 
                 if (CCSDitTQ) ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)=T4(a,b,c,d,i,j,k,l)",&
                                                           t4_err_tensor(itens),t4_tensor) 
                 if (CCSDitTQP) ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)=T5(a,b,c,d,e,i,j,k,l,m)",&
                                                           t5_err_tensor(itens),t5_tensor) 
              end if
              ierr=t_einsum("T2E(a,b,i,j)=T2(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
              ! correct the latest row of the b_diis array for the fact that we use an optimized error vector
              ierr=initialize_tensor_t(scalar_tensor)
              do i = 1, ntens
                 ierr=initialize_tensor_t(scalar_tensor)
                 if (.not.CCD) then
                    ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
                    if (CCSDitT) ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)",&
                   &                                        scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
                    if (CCSDitTQ) ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)",&
                   &                                        scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
                    if (CCSDitTQP) ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)",&
                   &                                        scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
                 end if
                 ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
                 b_diis(itens,i) = scalar_tens(1)
                 b_diis(i,itens) = dconjg(scalar_tens(1))
              end do
              ! replace latest result vector by the optimized result vector, this is also the tensor to be used in the next iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=t_einsum("T1R(a,i)=T1(a,i)",t1_res_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=t_einsum("T3R(a,b,c,i,j,k)=T3(a,b,c,i,j,k)",&
                &                                     t3_res_tensor(itens),t3_tensor)
                 if (CCSDitTQ) ierr=t_einsum("T4R(a,b,c,d,i,j,k,l)=T4(a,b,c,d,i,j,k,l)",&
                &                                      t4_res_tensor(itens),t4_tensor)
                 if (CCSDitTQP) ierr=t_einsum("T5R(a,b,c,d,e,i,j,k,l,m)=T5(a,b,c,d,e,i,j,k,l,m)",&
                &                                       t5_res_tensor(itens),t5_tensor)
              end if
              ierr=t_einsum("T2R(a,b,i,j)=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)
           else ! not using crop
              ! make the amplitudes for the new iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
           end if ! end of crop if
         end if

!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_inp_tensor)
            ierr=t_einsum("T1O(a,i)+=T1(a,i)",t1_inp_tensor,t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_inp_tensor)
               ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",&
                                      t3_inp_tensor,t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_inp_tensor)
               ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",&
                                      t4_inp_tensor,t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_inp_tensor)
               ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",&
              &                       t5_inp_tensor,t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_inp_tensor)
         ierr=t_einsum("T2O(a,b,i,j)+=T2(a,b,i,j)",t2_inp_tensor,t2_tensor)

         ierr=destroy_tensor_t(scalar_tensor)

       end subroutine diis

       subroutine diis (t5_tensor,t4_tensor,t3_tensor,t2_tensor,t1_tensor,t_convergence,mode)
         logical,parameter  :: crop=.true. ! if switching off, also increase dimension of maxdim (used to be 5)
         integer, parameter :: maxdim = 3 ! maximum number of tensors that are stored
         integer, parameter :: mindim = 2 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(tensor_t), intent(inout)  :: t2_tensor
         type(tensor_t), intent(inout), optional  :: t1_tensor
         type(tensor_t), intent(inout), optional  :: t3_tensor
         type(tensor_t), intent(inout), optional  :: t4_tensor
         type(tensor_t), intent(inout), optional  :: t5_tensor
         integer             :: i
         integer(INTD)       :: ierr, tens_rank
         integer(INTD), allocatable, save :: t1_dims(:)
         integer(INTD), allocatable, save :: t2_dims(:)
         integer(INTD), allocatable, save :: t3_dims(:)
         integer(INTD), allocatable, save :: t4_dims(:)
         integer(INTD), allocatable, save :: t5_dims(:)
         logical             :: init, clean_up, CCD
         logical             :: CCSDitT, CCSDitTQ, CCSDitTQP
         integer, save       :: itens, ntens
         complex(8), save    :: b_diis(maxdim,maxdim)
         complex(8)          :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         real(8)             :: b_max, b_scale
         type(tensor_t), save :: t1_inp_tensor
         type(tensor_t), save :: t2_inp_tensor
         type(tensor_t), save :: t3_inp_tensor
         type(tensor_t), save :: t4_inp_tensor
         type(tensor_t), save :: t5_inp_tensor
         type(tensor_t), save :: t1_res_tensor(maxdim), t1_err_tensor(maxdim)
         type(tensor_t), save :: t2_res_tensor(maxdim), t2_err_tensor(maxdim)
         type(tensor_t), save :: t3_res_tensor(maxdim), t3_err_tensor(maxdim)
         type(tensor_t), save :: t4_res_tensor(maxdim), t4_err_tensor(maxdim)
         type(tensor_t), save :: t5_res_tensor(maxdim), t5_err_tensor(maxdim)
         integer             :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)             :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(tensor_t)  :: scalar_tensor
         integer(INTD)       :: scalar_dims(1)
         complex(8), pointer :: scalar_tens(:)
         type(C_PTR):: body_p

         integer(kind(TAPP_datatype)) :: eti_storage_datatype

         init       = .false.
         clean_up   = .false.
         CCD = .false.
         CCSDitT = .false.
         CCSDitTQ = .false.
         CCSDitTQP = .false.

         eti_storage_datatype =  t2_tensor%datatype

         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.
         if (present(t3_tensor)) CCSDitT = .true.
         if (present(t4_tensor)) CCSDitTQ = .true.
         if (present(t5_tensor)) CCSDitTQP = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                ierr = get_nmodes(t1_tensor,tens_rank)
                ierr = get_extents(t1_tensor,t1_dims)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'

                if (CCSDitT) then 
                   ierr = get_nmodes(t3_tensor,tens_rank)
                   ierr = get_extents(t3_tensor,t3_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.6) stop 'input t3-tensor corrupted in diis'
                end if

                if (CCSDitTQ) then 
                   ierr = get_nmodes(t4_tensor,tens_rank)
                   ierr = get_extents(t4_tensor,t4_dims)
                   if (ierr.ne.0 .or. tens_rank.ne.8) stop 'input t4-tensor corrupted in diis'
                end if

                if (CCSDitTQP) then
                   ierr = get_nmodes(t5_tensor,tens_rank)
                   ierr = get_extents(t5_tensor,t5_dims) 
                   if (ierr.ne.0 .or. tens_rank.ne.10) stop 'input t5-tensor corrupted in diis'
                end if
            end if
            ierr = get_nmodes(t2_tensor,tens_rank)
            ierr = get_extents(t2_tensor,t2_dims) 
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis '
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) then
               ierr=construct_tensor_t(t1_inp_tensor,eti_storage_datatype,t1_dims,init_val=ZERO)
               if (CCSDitT) ierr=construct_tensor_t(t3_inp_tensor,eti_storage_datatype,t3_dims,init_val=ZERO)
               if (CCSDitTQ) ierr=construct_tensor_t(t4_inp_tensor,eti_storage_datatype,t4_dims,init_val=ZERO)
               if (CCSDitTQP) ierr=construct_tensor_t(t5_inp_tensor,eti_storage_datatype,t5_dims,init_val=ZERO)
            end if
            ierr=construct_tensor_t(t2_inp_tensor,eti_storage_datatype,t2_dims,init_val=ZERO)
            do i = 1, maxdim
               if (.not.CCD) then
                    ierr=construct_tensor_t(t1_res_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    ierr=construct_tensor_t(t1_err_tensor(i),eti_storage_datatype,t1_dims,init_val=ZERO)
                    if (CCSDitT) then
                       ierr=construct_tensor_t(t3_res_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t3_err_tensor(i),eti_storage_datatype,t3_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQ) then
                       ierr=construct_tensor_t(t4_res_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t4_err_tensor(i),eti_storage_datatype,t4_dims,init_val=ZERO)
                    end if
                    if (CCSDitTQP) then
                       ierr=construct_tensor_t(t5_res_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                       ierr=construct_tensor_t(t5_err_tensor(i),eti_storage_datatype,t5_dims,init_val=ZERO)
                    end if
               end if
               ierr=construct_tensor_t(t2_res_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               ierr=construct_tensor_t(t2_err_tensor(i),eti_storage_datatype,t2_dims,init_val=ZERO)
               if (ierr.ne.0) stop 'not enough space for diis'
            end do
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=initialize_tensor_t(t1_inp_tensor)
                ierr=t_einsum("T1O(a,i)+=T1N(a,i)",t1_inp_tensor,t1_tensor)
                if (CCSDitT) then
                   ierr=initialize_tensor_t(t3_inp_tensor)
                   ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3N(a,b,c,i,j,k)",t3_inp_tensor,t3_tensor)
                end if
                if (CCSDitTQ) then
                   ierr=initialize_tensor_t(t4_inp_tensor)
                   ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4N(a,b,c,d,i,j,k,l)",t4_inp_tensor,t4_tensor)
                end if
                if (CCSDitTQP) then
                   ierr=initialize_tensor_t(t5_inp_tensor)
                   ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5N(a,b,c,d,e,i,j,k,l,m)",t5_inp_tensor,t5_tensor)
                end if
            end if
            ierr=initialize_tensor_t(t2_inp_tensor)
            ierr=t_einsum("T2O(a,b,i,j)+=T2N(a,b,i,j)",t2_inp_tensor,t2_tensor)
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) then
               ierr=destroy_tensor_t(t1_inp_tensor)
               if (CCSDitT) ierr=destroy_tensor_t(t3_inp_tensor)
               if (CCSDitTQ) ierr=destroy_tensor_t(t4_inp_tensor)
               if (CCSDitTQP) ierr=destroy_tensor_t(t5_inp_tensor)
            end if
            ierr=destroy_tensor_t(t2_inp_tensor)
            do i = 1, maxdim
               if (.not.CCD) then
                  ierr=destroy_tensor_t(t1_res_tensor(i))
                  ierr=destroy_tensor_t(t1_err_tensor(i))
                  if (CCSDitT) then
                     ierr=destroy_tensor_t(t3_res_tensor(i))
                     ierr=destroy_tensor_t(t3_err_tensor(i))
                  end if
                  if (CCSDitTQ) then
                     ierr=destroy_tensor_t(t4_res_tensor(i))
                     ierr=destroy_tensor_t(t4_err_tensor(i))
                  end if
                  if (CCSDitTQP) then
                     ierr=destroy_tensor_t(t5_res_tensor(i))
                     ierr=destroy_tensor_t(t5_err_tensor(i))
                  end if
               end if
               ierr=destroy_tensor_t(t2_res_tensor(i))
               ierr=destroy_tensor_t(t2_err_tensor(i))
            end do
            return
         end if

!        define scalar tensor needed to extract the error
         scalar_dims(1) = 1
         ierr=construct_tensor_t(scalar_tensor,eti_storage_datatype,scalar_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(scalar_tensor,body_p)
         call c_f_pointer(body_p,scalar_tens,scalar_dims)

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_res_tensor(itens))
            ierr=t_einsum("T1R(a,i)+=T1(a,i)",t1_res_tensor(itens),t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_res_tensor(itens))
               ierr=t_einsum("T3R(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",t3_res_tensor(itens),t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_res_tensor(itens))
               ierr=t_einsum("T4R(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",t4_res_tensor(itens),t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_res_tensor(itens))
               ierr=t_einsum("T5R(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",t5_res_tensor(itens),t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_res_tensor(itens))
         ierr=t_einsum("T2R(a,b,i,j)+=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not. CCD) then
            ierr=initialize_tensor_t(t1_err_tensor(itens))
            ierr=t_einsum("T1E(a,i)+=T1R(a,i)",t1_err_tensor(itens),t1_tensor)
            ierr=t_einsum("T1E(a,i)+=T1I(a,i)",t1_err_tensor(itens),t1_inp_tensor,scale=MINUS_ONE)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_err_tensor(itens))
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3R(a,b,c,i,j,k)",t3_err_tensor(itens),t3_tensor)
               ierr=t_einsum("T3E(a,b,c,i,j,k)+=T3I(a,b,c,i,j,k)",t3_err_tensor(itens),t3_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_err_tensor(itens))
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4R(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_tensor)
               ierr=t_einsum("T4E(a,b,c,d,i,j,k,l)+=T4I(a,b,c,d,i,j,k,l)",&
              &                      t4_err_tensor(itens),t4_inp_tensor,scale=MINUS_ONE)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_err_tensor(itens))
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5R(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_tensor)
               ierr=t_einsum("T5E(a,b,c,d,e,i,j,k,l,m)+=T5I(a,b,c,d,e,i,j,k,l,m)",&
              &                      t5_err_tensor(itens),t5_inp_tensor,scale=MINUS_ONE)
            end if
         end if
         ierr=initialize_tensor_t(t2_err_tensor(itens))
         ierr=t_einsum("T2E(a,b,i,j)+=T2R(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
         ierr=t_einsum("T2E(a,b,i,j)+=T2I(a,b,i,j)",t2_err_tensor(itens),t2_inp_tensor,scale=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         ierr=initialize_tensor_t(scalar_tensor)
         do i = 1, ntens
            ierr=initialize_tensor_t(scalar_tensor)
            if (.not.CCD) then
               ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
               if (CCSDitT) then
                  ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)", &
                &     scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
               end if
               if (CCSDitTQ) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)", &
                &     scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
               end if
               if (CCSDitTQP) then
                  ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)", &
                &     scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
               end if
            end if
            ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
            b_diis(itens,i) = scalar_tens(1)
            b_diis(i,itens) = dconjg(scalar_tens(1))
            ! print *,'b_diis(',i,',',itens,'):',b_diis(i,itens)
         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             write(*,*) "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
         !   print *,'coeff:',coeff
           if (crop) then
              ! replace latest error vector by the optimized error vector using t1_tensor and t2_tensor as scratch space
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if 
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_err_tensor(i),scale=coeff(i))
                     if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)",&
                  &         t3_tensor,t3_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                  &         t4_tensor,t4_err_tensor(i),scale=coeff(i))
                     if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                  &         t5_tensor,t5_err_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_err_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=copy_tensor_t(t1_err_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=copy_tensor_t(t3_err_tensor(itens),t3_tensor) 
                 if (CCSDitTQ) ierr=copy_tensor_t(t4_err_tensor(itens),t4_tensor) 
                 if (CCSDitTQP) ierr=copy_tensor_t(t5_err_tensor(itens),t5_tensor) 
              end if
              ierr=copy_tensor_t(t2_err_tensor(itens),t2_tensor)
              ! correct the latest row of the b_diis array for the fact that we use an optimized error vector
              ierr=initialize_tensor_t(scalar_tensor)
              do i = 1, ntens
                 ierr=initialize_tensor_t(scalar_tensor)
                 if (.not.CCD) then
                    ierr=t_einsum("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
                    if (CCSDitT) ierr=t_einsum("E()+=T+(a,b,c,i,j,k)*T(a,b,c,i,j,k)",&
                   &                                        scalar_tensor,t3_err_tensor(itens),t3_err_tensor(i))
                    if (CCSDitTQ) ierr=t_einsum("E()+=T+(a,b,c,d,i,j,k,l)*T(a,b,c,d,i,j,k,l)",&
                   &                                        scalar_tensor,t4_err_tensor(itens),t4_err_tensor(i))
                    if (CCSDitTQP) ierr=t_einsum("E()+=T+(a,b,c,d,e,i,j,k,l,m)*T(a,b,c,d,e,i,j,k,l,m)",&
                   &                                        scalar_tensor,t5_err_tensor(itens),t5_err_tensor(i))
                 end if
                 ierr=t_einsum("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
                 b_diis(itens,i) = scalar_tens(1)
                 b_diis(i,itens) = dconjg(scalar_tens(1))
              end do
              ! replace latest result vector by the optimized result vector, this is also the tensor to be used in the next iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
              if (.not.CCD) then
                 ierr=copy_tensor_t(t1_res_tensor(itens),t1_tensor)
                 if (CCSDitT) ierr=copy_tensor_t(t3_res_tensor(itens),t3_tensor)
                 if (CCSDitTQ) ierr=copy_tensor_t(t4_res_tensor(itens),t4_tensor)
                 if (CCSDitTQP) ierr=copy_tensor_t(t5_res_tensor(itens),t5_tensor)
              end if
              ierr=copy_tensor_t(t2_res_tensor(itens),t2_tensor)
           else ! not using crop
              ! make the amplitudes for the new iteration
              if (.not.CCD) then
                 ierr=initialize_tensor_t(t1_tensor)
                 if (CCSDitT) ierr=initialize_tensor_t(t3_tensor)
                 if (CCSDitTQ) ierr=initialize_tensor_t(t4_tensor)
                 if (CCSDitTQP) ierr=initialize_tensor_t(t5_tensor)
              end if
              ierr=initialize_tensor_t(t2_tensor)
              do i = 1, ntens
                 if (.not.CCD) then
                    ierr=t_einsum("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
                    if (CCSDitT) ierr=t_einsum("T3(a,b,c,i,j,k)+=T3O(a,b,c,i,j,k)", &
                   &                                    t3_tensor,t3_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQ) ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=T4O(a,b,c,d,i,j,k,l)",&
                   &                                     t4_tensor,t4_res_tensor(i),scale=coeff(i))
                    if (CCSDitTQP) ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=T5O(a,b,c,d,e,i,j,k,l,m)",&
                   &                                     t5_tensor,t5_res_tensor(i),scale=coeff(i))
                 end if
                 ierr=t_einsum("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
              end do
           end if ! end of crop if 
         end if

!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=initialize_tensor_t(t1_inp_tensor)
            ierr=t_einsum("T1O(a,i)+=T1(a,i)",t1_inp_tensor,t1_tensor)
            if (CCSDitT) then
               ierr=initialize_tensor_t(t3_inp_tensor)
               ierr=t_einsum("T3O(a,b,c,i,j,k)+=T3(a,b,c,i,j,k)",&
                                      t3_inp_tensor,t3_tensor)
            end if
            if (CCSDitTQ) then
               ierr=initialize_tensor_t(t4_inp_tensor)
               ierr=t_einsum("T4O(a,b,c,d,i,j,k,l)+=T4(a,b,c,d,i,j,k,l)",&
                                      t4_inp_tensor,t4_tensor)
            end if
            if (CCSDitTQP) then
               ierr=initialize_tensor_t(t5_inp_tensor)
               ierr=t_einsum("T5O(a,b,c,d,e,i,j,k,l,m)+=T5(a,b,c,d,e,i,j,k,l,m)",&
              &                       t5_inp_tensor,t5_tensor)
            end if
         end if
         ierr=initialize_tensor_t(t2_inp_tensor)
         ierr=t_einsum("T2O(a,b,i,j)+=T2(a,b,i,j)",t2_inp_tensor,t2_tensor)

         ierr=destroy_tensor_t(scalar_tensor)

       end subroutine diis
#endif

end module eti_diis
