module talsh_common_routines

!This module contains subroutines both needed in solving t equations as well as Lambda equations.
!get_orbital_energies,diis,scale_with_denominators

    use talsh
    use tensor_algebra

    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)
    private

    public get_orbital_energies
    public scale_with_denominators
    public diis

    contains

      subroutine diis (t1_tensor,t2_tensor,t_convergence,mode)

       integer, parameter :: maxdim = 5 ! maximum number of tensors that will be stored 
       integer, parameter :: mindim = 3 ! minimum number of tensors to allow extrapolation
       integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
       real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

       type(talsh_tens_t), intent(inout)  :: t1_tensor, t2_tensor
       integer            :: ierr, i, tens_rank
       integer, save      :: t1_dims(1:2), t2_dims(1:4)
       logical            :: init, clean_up
       integer, save      :: itens, jtens, ntens
       complex(8), save   :: b_diis(maxdim,maxdim)
       complex(8)         :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
       real(8)            :: b_max, b_scale
       type(talsh_tens_t), save :: t1_inp_tensor, t2_inp_tensor
       type(talsh_tens_t), save :: t1_res_tensor(maxdim), t2_res_tensor(maxdim)
       type(talsh_tens_t), save :: t1_err_tensor(maxdim), t2_err_tensor(maxdim)
       integer            :: ipiv(0:maxdim), info  ! needed for ZHESV solver
       real(8)            :: work(maxdim,maxdim)   ! needed for ZHESV solver
       type(talsh_tens_t), save  :: scalar_tensor
       integer:: scalar_dims(1)
       complex(8), pointer :: scalar_tens(:)
       type(C_PTR):: body_p

       init       = .false.
       clean_up   = .false.

       if (present(mode)) then
          if (mode.eq.1) init       = .true.
          if (mode.eq.2) clean_up   = .true.
       end if 

!      define scalar tensor needed to extract the error
       scalar_dims(1) = 1
       ierr=talsh_tensor_construct(scalar_tensor,C8,scalar_dims(1:0),init_val=ZERO)
       ierr=talsh_tensor_get_body_access(scalar_tensor,body_p,C8,0,DEV_HOST)
       call c_f_pointer(body_p,scalar_tens,scalar_dims)

       if (init) then
!         get dimensions of the tensors
          ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
          ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!         initialize DIIS matrix
          b_diis      = ZERO
!         create storage for input, result and error tensors (in memory for the time being)
          ierr=talsh_tensor_construct(t1_inp_tensor,C8,t1_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(t2_inp_tensor,C8,t2_dims,init_val=ZERO)
          do i = 1, maxdim
             ierr=talsh_tensor_construct(t1_res_tensor(i),C8,t1_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(t1_err_tensor(i),C8,t1_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(t2_res_tensor(i),C8,t2_dims,init_val=ZERO)
             ierr=talsh_tensor_construct(t2_err_tensor(i),C8,t2_dims,init_val=ZERO)
          end do
          itens = 0
          ntens = 0
!         store the tensors in memory as current input tensors, needed to calculate the error vector later on
          ierr=talsh_tensor_init(t1_inp_tensor)
          ierr=talsh_tensor_add("T1O(a,i)+=T1N(a,i)",t1_inp_tensor,t1_tensor)
          ierr=talsh_tensor_init(t2_inp_tensor)
          ierr=talsh_tensor_add("T2O(a,b,i,j)+=T2N(a,b,i,j)",t2_inp_tensor,t2_tensor)
          return ! no extrapolation possible
       end if

       if (clean_up) then
!         de-allocate tensors
          ierr=talsh_tensor_destruct(t1_inp_tensor)
          ierr=talsh_tensor_destruct(t2_inp_tensor)
          do i = 1, maxdim
             ierr=talsh_tensor_destruct(t1_res_tensor(i))
             ierr=talsh_tensor_destruct(t1_err_tensor(i))
             ierr=talsh_tensor_destruct(t2_res_tensor(i))
             ierr=talsh_tensor_destruct(t2_err_tensor(i))
          end do
          return
       end if

!      update counters
       itens = itens + 1                                                                                                                    
       if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
       ntens = ntens + 1                    ! number of tensors stored after returning
       if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!      store the result tensors
       ierr=talsh_tensor_init(t1_res_tensor(itens))
       ierr=talsh_tensor_add("T1R(a,i)+=T1(a,i)",t1_res_tensor(itens),t1_tensor)
       ierr=talsh_tensor_init(t2_res_tensor(itens))
       ierr=talsh_tensor_add("T2R(a,b,i,j)+=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)

!      calculate the error tensors by comparing the result of the t-equations with the tensors used as input
       ierr=talsh_tensor_init(t1_err_tensor(itens))
       ierr=talsh_tensor_add("T1E(a,i)+=T1R(a,i)",t1_err_tensor(itens),t1_tensor)
       ierr=talsh_tensor_add("T1E(a,i)+=T1I(a,i)",t1_err_tensor(itens),t1_inp_tensor,scale=MINUS_ONE)
       ierr=talsh_tensor_init(t2_err_tensor(itens))
       ierr=talsh_tensor_add("T2E(a,b,i,j)+=T2R(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
       ierr=talsh_tensor_add("T2E(a,b,i,j)+=T2I(a,b,i,j)",t2_err_tensor(itens),t2_inp_tensor,scale=MINUS_ONE)

!      update the b_diis array that holds the scalar products of the error tensors
       ierr=talsh_tensor_init(scalar_tensor)
       do i = 1, ntens
          ierr=talsh_tensor_init(scalar_tensor)
          ierr=talsh_tensor_contract("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
          ierr=talsh_tensor_contract("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
          b_diis(itens,i) = scalar_tens(1)
          b_diis(i,itens) = dconjg(scalar_tens(1))
       end do

       if (present(t_convergence)) then
          t_convergence = dsqrt(dble(b_diis(itens,itens)))
       end if 

!      transfer to a work array in the way requested by the solver and normalize the matrix
       b_max = maxval (abs(b_diis))
       if (b_max .gt. 1.D-30) then
           b_scale = 1.D0 / b_max
       else
           print*, "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
           return
       end if
       bb = ZERO
       bb(1:ntens,0) = MINUS_ONE
       bb(0,1:ntens) = MINUS_ONE
       bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!      solve the linear system of equations using the LAPACK solver ZHESV
       coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
       coeff(1:ntens) = ZERO 
       call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!      make the extrapolated amplitudes if the space is large enough
       if (ntens.ge.mindim) then
         ierr=talsh_tensor_init(t1_tensor)
         ierr=talsh_tensor_init(t2_tensor)
          do i = 1, ntens
             ierr=talsh_tensor_add("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
             ierr=talsh_tensor_add("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
          end do
       end if

!      the extrapolated amplitudes are the new input tensors
       ierr=talsh_tensor_init(t1_inp_tensor)
       ierr=talsh_tensor_add("T1O(a,i)+=T1(a,i)",t1_inp_tensor,t1_tensor)
       ierr=talsh_tensor_init(t2_inp_tensor)
       ierr=talsh_tensor_add("T2O(a,b,i,j)+=T2(a,b,i,j)",t2_inp_tensor,t2_tensor)

       ierr=talsh_tensor_destruct(scalar_tensor)

      end subroutine diis

      subroutine scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor,t2_tensor,t3_tensor,eps_ijk)

       type(talsh_tens_t), optional, intent(inout)  :: t1_tensor, t2_tensor, t3_tensor
       real(8), intent(in)                          :: eps_occ(:),eps_vir(:)
       real(8), intent(in), optional                :: eps_ijk

       real(8)  :: denominator
       complex(8), pointer, contiguous:: t1_tens(:,:), t2_tens(:,:,:,:), t3_tens(:,:,:)
       type(C_PTR):: body_p

       integer:: t1_dims(1:2), t2_dims(1:4), t3_dims(1:3)

       integer   :: ierr, tens_rank, a, b, c, i, j

       integer, optional, intent(in) :: nocc

       if (present(t1_tensor)) then

          ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
          if (ierr.ne.0 .or. tens_rank.ne.2) stop 'program error: t1 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t1_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t1_tens,t1_dims(1:2)) ! to use <t1_tens> as a regular Fortran 2d array

          if (t1_dims(2) .eq. nocc*2) then
           !t amplitudes
           do i = 1, t1_dims(2)
             do a = 1, t1_dims(1)
                denominator = eps_occ(i) - eps_vir(a)
                t1_tens(a,i) = t1_tens(a,i) / denominator
             end do
           end do
          else
           !Lambda
           do i = 1, t1_dims(1)
             do a = 1, t1_dims(2)
                denominator = eps_occ(i) - eps_vir(a)
                t1_tens(i,a) = t1_tens(i,a) / denominator
             end do
           end do
          end if

       end if

       if (present(t2_tensor)) then

          ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
          if (ierr.ne.0 .or. tens_rank.ne.4) stop 'program error: t2 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t2_tens,t2_dims(1:4)) ! to use <t2_tens> as a regular Fortran 4d array

          if (t2_dims(3) .eq. nocc*2) then
           !t amplitudes
           do j = 1, t2_dims(4)
             do i = 1, t2_dims(3)
               do b = 1, t2_dims(2)
                 do a = 1, t2_dims(1)
                    denominator = eps_occ(i) + eps_occ(j) - eps_vir(a) - eps_vir(b)
                    t2_tens(a,b,i,j) = t2_tens(a,b,i,j) / denominator
                 end do
               end do
             end do
           end do
          else
           !Lambda
           do i = 1, t2_dims(1)
             do j = 1, t2_dims(2)
               do a = 1, t2_dims(3)
                 do b = 1, t2_dims(4)
                    denominator = eps_occ(i) + eps_occ(j) - eps_vir(a) - eps_vir(b)
                    t2_tens(i,j,a,b) = t2_tens(i,j,a,b) / denominator
                 end do
               end do
             end do
           end do
          end if

       end if

       if (present(t3_tensor)) then

          ierr = talsh_tensor_dimensions(t3_tensor,tens_rank,t3_dims)
          if (ierr.ne.0 .or. tens_rank.ne.3) stop 'program error: t3 tensor corrupted'
          ierr=talsh_tensor_get_body_access(t3_tensor,body_p,C8,0,DEV_HOST)
          call c_f_pointer(body_p,t3_tens,t3_dims(1:3)) ! to use <t3_tens> as a regular Fortran 3d array

          do c = 1, t3_dims(3)
            do b = 1, t3_dims(2)
              do a = 1, t3_dims(1)
                 denominator = eps_ijk - eps_vir(a) - eps_vir(b) - eps_vir(c)
                 t3_tens(a,b,c) = t3_tens(a,b,c) / denominator
              end do
            end do
          end do

       end if

      end subroutine scale_with_denominators

      subroutine get_orbital_energies (nmo, mo_list, eps)

!This subroutine reads the desired subset of mo energies from DIRAC

!      Written by Lucas Visscher, winter 2016/2017 (but in Goa, India, temperature about 25 C)

       use dirac_MO
       use exacorr_global

       integer, intent(in) :: nmo          ! the length of the mo basis
       integer, intent(in) :: mo_list(:)   ! and their indices
       real(8) :: eps(:)

       integer             :: k
       type(MO)            :: qmo

       call get_mo_coefficients (qmo,mo_list,nmo)

       do k = 1, nmo
          eps(k)     = qmo%energy(k)
          eps(k+nmo) = qmo%energy(k)
       end do
       call delete_qmo (qmo)

      end subroutine get_orbital_energies

end module talsh_common_routines
