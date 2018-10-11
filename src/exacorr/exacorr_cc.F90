module exacorr_cc
!This module contains the routines needed to calculate ae coupled cluster wave function

!       use tensor_algebra
        use exatensor

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0),ONE_QUARTER=(0.25D0,0.D0)!, &
!                                THREE=(3.0,0.0), SIX=(6.0,0.0)
        private

        public cc_driver

       contains
!-------------------------------------------------------------------------------------------
        subroutine cc_driver (nocc,nvir,mo_occ,mo_vir)

!This subroutine drives the ccsd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code
 
!        Written by Lucas Visscher, winter 2016/2017
!        ExaTensor conversion by Anastasios Papadopoulos, October 2018
 
         use exacorr_global

         integer, intent(in) :: nocc,nvir           ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in) :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer                        :: ierr, i
         complex(8), pointer            :: energy(:),result_tens(:)
         integer(C_SIZE_T), parameter   :: BUF_SIZE=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_SIZE_T)              :: host_buf_size
         integer(C_INT)                 :: host_arg_max

!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_tens_t) :: foo_tensor, fvv_tensor, fov_tensor
         type(talsh_tens_t) :: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor, ooov_tensor, vovv_tensor
!        intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t) :: goint_tensor, gvint_tensor, hoint_tensor, hvint_tensor, hovint_tensor
         type(talsh_tens_t) :: aint_tensor, bint_tensor, cint_tensor, hint_tensor
         type(talsh_tens_t) :: tau_tensor, tau2_tensor
         integer            :: goint_dims(1:2),gvint_dims(1:2),hoint_dims(1:2),hvint_dims(1:2),hovint_dims(1:2)
         integer            :: aint_dims(1:4),bint_dims(1:4),cint_dims(1:4),hint_dims(1:4)
!        auxilliary tensors for t1 equation
         type(talsh_tens_t) :: hovaux_tensor, hoaux_tensor
!        solution tensors
         type(talsh_tens_t) :: t1_tensor, t2_tensor, s1_tensor, s2_tensor
         integer            :: t1_dims(1:2), t2_dims(1:4), s1_dims(1:2), s2_dims(1:4)
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor, one_tensor
         integer            :: result_dims(1),one_dims(1)
!        CCSD control and result variables
         integer            :: iccsd, nccsd
         real(8)            :: mp2_energy, ccsd_energy!, t_energy(3), t_convergence  ! triples not implemented
!        real(8), parameter :: t_target_precision = 1.D-13

!        Make some noise so that we know we are here
         print*, " Arrived in cc_driver routine"

!        Initialize libraries
         host_buf_size=BUF_SIZE
         ierr=talsh_init(host_buf_size,host_arg_max)
         print*, " Initialized talsh library"
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
         call interest_initialize(.true.)
         print*, " Initialized interest library"
         call initialize_global_data
         print*, " Initialized global data"

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,0,DEV_HOST)

!        Get fixed one-body tensors ! set to 0 now at the beginning
!                                 fo_tensor,fov_tensor,fv_tensor)
         foo_tensor = 0
         fov_tensor = 0
         fvv_tensor = 0   

!        Get fixed two-body tensors
         call get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
                                oooo_tensor,ooov_tensor,oovv_tensor, &
                                vovo_tensor,vovv_tensor,vvvv_tensor)

!        Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
         t1_dims(1) = nvir * 2
         t1_dims(2) = nocc * 2
         ierr=talsh_tensor_construct(t1_tensor,C8,t1_dims,init_val=ZERO)

!        Initialize t2 amplitudes with the MP2 values
         t2_dims(1:2) = nvir * 2
         t2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(t2_tensor,C8,t2_dims,init_val=ZERO)
!        ierr=talsh_tensor_add("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor,one_tensor) ! o(2)v(2)
!        end of work-around
!        call scale_with_denominators (eps_occ,eps_vir,t1_tensor,t2_tensor)

!        Calculate and print the MP2 energy
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,t2_tensor,oovv_tensor) ! o(2)v(2)
         mp2_energy=result_tens(1) * ONE_QUARTER
         print*," MP2 energy = ", mp2_energy

!        Start CCSD iterations

         nccsd = 30 ! hard wired for the moment

         goint_dims = nocc  * 2
         ierr=talsh_tensor_construct(goint_tensor,C8,goint_dims,init_val=ZERO)
         gvint_dims = nvir * 2
         ierr=talsh_tensor_construct(gvint_tensor,C8,gvint_dims,init_val=ZERO)
         hoint_dims = nocc  * 2
         ierr=talsh_tensor_construct(hoint_tensor,C8,hoint_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(hoaux_tensor,C8,hoint_dims,init_val=ZERO)
         hvint_dims = nvir * 2
         ierr=talsh_tensor_construct(hvint_tensor,C8,hvint_dims,init_val=ZERO)
         hovint_dims(1) = nocc * 2
         hovint_dims(2) = nvir * 2
         ierr=talsh_tensor_construct(hovint_tensor,C8,hovint_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(hovaux_tensor,C8,hovint_dims,init_val=ZERO)
         aint_dims  = nocc  * 2
         ierr=talsh_tensor_construct(aint_tensor,C8,aint_dims,init_val=ZERO)
         cint_dims(1) = nocc * 2
         cint_dims(2) = nvir * 2
         cint_dims(3) = nocc * 2
         cint_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(cint_tensor,C8,cint_dims,init_val=ZERO)
         hint_dims(1) = nvir * 2
         hint_dims(2) = nocc * 2
         hint_dims(3) = nvir * 2
         hint_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(hint_tensor,C8,hint_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(tau_tensor, C8,t2_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(tau2_tensor,C8,t2_dims,init_val=ZERO)
         s1_dims(1) = nvir * 2
         s1_dims(2) = nocc * 2
         ierr=talsh_tensor_construct(s1_tensor,C8,s1_dims,init_val=ZERO)
         s2_dims(1:2) = nvir * 2
         s2_dims(3:4) = nocc * 2
         ierr=talsh_tensor_construct(s2_tensor,C8,s2_dims,init_val=ZERO)
!        call diis (t1_tensor,t2_tensor,mode=1)
         if (ierr.eq.TALSH_SUCCESS) then
            print*," ccdriver: allocated all tensors"
         else
            print*," ccdriver: error in tensor allocation, code=",ierr
            return
         end if

         do iccsd = 1, nccsd

!        make tau ! validated
         ierr=talsh_tensor_init(tau_tensor)
         ierr=talsh_tensor_add("T(a,b,i,j)+=R(a,b,i,j)",tau_tensor,t2_tensor)
         ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",tau_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
         ierr=talsh_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",tau_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE) ! o(2)v(2)

!NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!        calculate intermediaries

         ierr=talsh_tensor_init(s2_tensor)

!        The vvvv contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,vvvv_tensor,tau_tensor,scale=ONE_HALF) !validated o(2)v(4)

!        vv ! validated
         ierr=talsh_tensor_init(hvint_tensor)
         ierr=talsh_tensor_add("H(a,c)+=F(a,c)",hvint_tensor,fv_tensor)
         ierr=talsh_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",hvint_tensor,oovv_tensor,tau_tensor, &
                                     scale=ONE_HALF) ! o(2)v(3)
         ierr=talsh_tensor_init(gvint_tensor)
         ierr=talsh_tensor_add("G(a,c)+=H(a,c)",gvint_tensor,hvint_tensor)
         ierr=talsh_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",gvint_tensor,vovv_tensor,t1_tensor) ! o(1)v(3)
         ierr=talsh_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",gvint_tensor,fov_tensor,t1_tensor,scale=MINUS_ONE) ! o(1)v(2)

!        vovo ! moretests
         ierr=talsh_tensor_init(tau2_tensor)
         ierr=talsh_tensor_add("T(a,d,l,i)+=R(a,d,l,i)",tau2_tensor,t2_tensor,scale=ONE_HALF)
         ierr=talsh_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,t1_tensor,t1_tensor) ! o(2)v(2)
         ierr=talsh_tensor_init(hint_tensor)
         ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",hint_tensor,oovv_tensor,tau2_tensor)! o(3)v(3)
         ierr=talsh_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",hint_tensor,vovv_tensor,t1_tensor) ! o(2)v(3)
         ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",hint_tensor,ooov_tensor,t1_tensor,scale=MINUS_ONE) ! o(3)v(2)
         ierr=talsh_tensor_add("H(a,k,c,i)+=V(a,k,c,i)",hint_tensor,vovo_tensor)

!        This term appears four times as we need to permute both i and j and a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,hint_tensor,t2_tensor)  ! validated o(3)v(3)

!        ovoo ! validated
         ierr=talsh_tensor_init(cint_tensor)
         ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,vovv_tensor,tau_tensor,scale=MINUS_ONE_HALF) ! o(3)v(3)
         ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,vovo_tensor,t1_tensor) ! o(3)v(2)
         ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,vovo_tensor,t1_tensor,scale=MINUS_ONE) ! o(3)v(2)
         ierr=talsh_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,ooov_tensor,one_tensor) ! o(3)v(1)

!        oo ! validated
         ierr=talsh_tensor_init(hoint_tensor)
         ierr=talsh_tensor_add("H(k,j)+=F(k,j)",hoint_tensor,fo_tensor)
         ierr=talsh_tensor_contract("H(k,j)+=V(k,l,c,d)*T(c,d,l,j)",hoint_tensor,oovv_tensor,tau_tensor, &
                                     scale=ONE_HALF) ! o(3)v(2)
         ierr=talsh_tensor_init(goint_tensor)
         ierr=talsh_tensor_add("G(k,j)+=H(k,j)",goint_tensor,hoint_tensor)
         ierr=talsh_tensor_contract("G(k,j)+=V(k,l,j,c)*T(c,l)",goint_tensor,ooov_tensor,t1_tensor,scale=MINUS_ONE) ! o(3)v(1)
         ierr=talsh_tensor_contract("G(k,j)+=F(k,c)*T(c,j)",goint_tensor,fov_tensor,t1_tensor) ! o(2)v(1)

!        ov ! validated
         ierr=talsh_tensor_init(hovint_tensor)
         ierr=talsh_tensor_init(hovaux_tensor)
         ierr=talsh_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",hovint_tensor,oovv_tensor,t1_tensor) ! o(2)v(2)
         ierr=talsh_tensor_add("L(k,c)+=H(k,c)",hovaux_tensor,hovint_tensor)                                  ! copy common term
         ierr=talsh_tensor_add("H(k,c)+=F(k,c)",hovint_tensor,fov_tensor)                                     ! add fock matrix
         ierr=talsh_tensor_contract("H(k,c)+=F(k,c)*K()",hovaux_tensor,fov_tensor,one_tensor,scale=MINUS_ONE) ! subtract fock matrix

!        oooo ! validated
         ierr=talsh_tensor_init(aint_tensor)
         ierr=talsh_tensor_add("A(k,l,i,j)+=V(k,l,i,j)",aint_tensor,oooo_tensor)
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",aint_tensor,ooov_tensor,t1_tensor) ! o(4)v(1)
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",aint_tensor,ooov_tensor,t1_tensor, &
                                     scale=MINUS_ONE) ! o(4)v(1)
         ierr=talsh_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",aint_tensor,oovv_tensor,tau_tensor, &
                                     scale=ONE_HALF) ! o(4)v(2)

         if (ierr.ne.TALSH_SUCCESS) then
            print*," ccdriver: error in calculation of intermediates, code=",ierr
            return
         end if

         if (debug) then
!        For debugging: print quasi-norm (NOT using complex conjugation so that we can also see whether tensors are imaginary)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=T(a,i)*T(a,i)",result_tensor,t1_tensor,t1_tensor)
            print*, "T1amp = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,t2_tensor,t2_tensor,scale=ONE_QUARTER)
            print*, "T2amp = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=T(a,b,i,j)*T(a,b,i,j)",result_tensor,tau_tensor,tau_tensor,scale=ONE_QUARTER)
            print*, "Tau   = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=H(k,j)*H(k,j)",result_tensor,hoint_tensor,hoint_tensor)
            print*, "hoint = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=G(k,j)*G(k,j)",result_tensor,goint_tensor,goint_tensor)
            print*, "goint = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=H(a,c)*H(a,c)",result_tensor,hvint_tensor,hvint_tensor)
            print*, "hvint = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=G(a,c)*G(a,c)",result_tensor,gvint_tensor,gvint_tensor)
            print*, "gvint = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=H(k,c)*H(k,c)",result_tensor,hovint_tensor,hovint_tensor)
            print*, "hovint = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=A(k,l,i,j)*A(k,l,i,j)",result_tensor,aint_tensor,aint_tensor)
            print*, "a int = ",result_tens(1)/4.D0
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=H(a,i,c,k)*H(a,i,c,k)",result_tensor,hint_tensor,hint_tensor)
            print*, "h int = ",result_tens(1)
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("S()+=C(k,b,i,j)*C(k,b,i,j)",result_tensor,cint_tensor,cint_tensor)
            print*, "c int = ",result_tens(1)
         end if

!        This term appears twice as we need to permute i and j
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,goint_tensor,t2_tensor)                   ! validated o(3)v(2)
!        Again two times, also note the complex conjugation as we use vovv instead of vvvo
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,vovv_tensor,t1_tensor)                   ! validated o(2)v(3)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,vovv_tensor,t1_tensor,scale=MINUS_ONE)   ! validated o(2)v(3)

!        This term appears twice as we need to permute a and b
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,gvint_tensor,t2_tensor)  ! validated o(2)v(3)
!        Again two times, with the C intermediate (note: in RELCCSD this is done differently)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,t1_tensor,scale=MINUS_ONE)    ! o(3)v(2)
         ierr=talsh_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,t1_tensor)                    ! o(3)v(2)

!        The oooo contribution
         ierr=talsh_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,aint_tensor,tau_tensor,scale=ONE_HALF) !validated o(4)v(2)

!        ierr=talsh_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,oovv_tensor,one_tensor)    ! validated o(2)v(2)
!        end of work-around

!        calculate s1 tensor
         ierr=talsh_tensor_init(s1_tensor)
!        ierr=talsh_tensor_add("S(a,i)+=F+(i,a)",s1_tensor,fov_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,fov_tensor,one_tensor)
!        end of work-around
         ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,vovv_tensor,tau_tensor,scale=ONE_HALF) !validated o(2)v(3)
         ierr=talsh_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,hvint_tensor,t1_tensor)                        !validated o(1)v(2)
         ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoint_tensor,t1_tensor)                        !validated o(2)v(1)
         ierr=talsh_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,hovint_tensor,t2_tensor)                   !validated o(2)v(2)
         ierr=talsh_tensor_init(hoaux_tensor)                                                                        !moretests
         ierr=talsh_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoaux_tensor,hovaux_tensor,t1_tensor)                    !moretests o(2)v(1)
         ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoaux_tensor,t1_tensor)                        !moretests o(2)v(1)
         ierr=talsh_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,ooov_tensor,tau_tensor,scale=ONE_HALF) !validated o(3)v(2)
         ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,vovo_tensor,t1_tensor,scale=MINUS_ONE)     !validated o(2)v(2)

         if (ierr.ne.TALSH_SUCCESS) then
            print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
            return
         end if

!        prepare for next iteration
         ierr=talsh_tensor_init(t1_tensor)
         ierr=talsh_tensor_add("T1(a,i)+=S1(a,i)",t1_tensor,s1_tensor)
         ierr=talsh_tensor_init(t2_tensor)
         ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",t2_tensor,s2_tensor)

!        call diis (t1_tensor,t2_tensor,t_convergence=t_convergence)

         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,i)*F(i,a)",result_tensor,t1_tensor,fov_tensor) ! o(1)v(1)
         ccsd_energy=result_tens(1)
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,tau_tensor,oovv_tensor) ! o(2)v(2)
         ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER
         write (*,'(A,I4,A,F20.15,A,E10.1)') "  Iteration ",iccsd," CCSD energy = ",ccsd_energy, " convergence =",t_convergence
         if (t_convergence .lt. t_target_precision) exit

         end do

!        Calculate (T) energy correction, not implemented yet
!        call t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor,ooov_tensor,vovv_tensor,oovv_tensor)
!        write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
!        write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
!        write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
!        write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",ccsd_energy+t_energy(1)
!        write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(2)
!        write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(3)

!        Destruct tensors and shut down library
!        call diis (t1_tensor,t2_tensor,mode=2)
         ierr=talsh_tensor_destruct(oooo_tensor)
         ierr=talsh_tensor_destruct(oovv_tensor)
         ierr=talsh_tensor_destruct(vovo_tensor)
         ierr=talsh_tensor_destruct(vvvv_tensor)
         ierr=talsh_tensor_destruct(ooov_tensor)
         ierr=talsh_tensor_destruct(vovv_tensor)
         ierr=talsh_tensor_destruct(s1_tensor)
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(t1_tensor)
         ierr=talsh_tensor_destruct(t2_tensor)
         ierr=talsh_tensor_destruct(hoint_tensor)
         ierr=talsh_tensor_destruct(hoaux_tensor)
         ierr=talsh_tensor_destruct(hvint_tensor)
         ierr=talsh_tensor_destruct(hovint_tensor)
         ierr=talsh_tensor_destruct(hovaux_tensor)
         ierr=talsh_tensor_destruct(goint_tensor)
         ierr=talsh_tensor_destruct(gvint_tensor)
         ierr=talsh_tensor_destruct(aint_tensor)
         ierr=talsh_tensor_destruct(cint_tensor)
         ierr=talsh_tensor_destruct(hint_tensor)
!        Get statistics 
         ierr=talsh_stats()
         ierr=talsh_shutdown()

!        Free global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         print*," Final CCSD energy    = ",ccsd_energy
!        print*," Final CCSD(T) energy = ",ccsd_energy+t_energy(1)+t_energy(2)  ! no (T) implementation yet
         print*, " Leaving cc_driver routine"
         return
        end subroutine cc_driver

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir, &
                                     oooo_tensor,ooov_tensor,oovv_tensor, &
                                     vovo_tensor,vovv_tensor,vvvv_tensor)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Lucas Visscher, summer 2017
!        ExaTensor conversion by Anastasios Papadopoulos, October 2018   
 
         use exacorr_datatypes
         use exacorr_ao_to_mo
         use exacorr_global

         integer, intent(in)               :: nocc,nvir             ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)               :: mo_occ(:),mo_vir(:)   ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout) :: oovv_tensor, vvvv_tensor, oooo_tensor, vovo_tensor, ooov_tensor, vovv_tensor
         integer                           :: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4),ooov_dims(1:4),vovv_dims(1:4)

         integer                            :: ierr, i
         complex(8), pointer, contiguous    :: ao_tens(:,:,:,:)

         integer    :: nao

!        auxilliary integral tensors (only needed inside this routine)
         type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
         integer            :: aoint_dims(1:4), ovvo_dims(1:4)

!        arrays needed to communicate with the MO integral generator
         integer                :: nmo(4)
         integer, allocatable   :: mo_list(:)

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: minusone_tensor
         integer            :: minusone_dims(1)

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,0,DEV_HOST)
         call compute_ao_integrals (nao,ao_tens)

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oovv_dims(1:2) = nocc * 2
         oovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(oovv_tensor,C8,oovv_dims,init_val=ZERO)
!        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
         allocate (mo_list(2*nocc+2*nvir))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nocc
         nmo(4) = nvir
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nocc+nvir+1:2*nocc+nvir)     = mo_occ(1:nocc)
         mo_list(2*nocc+nvir+1:2*nocc+2*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,oovv_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get antisymmetric oooo tensor
         oooo_dims = nocc * 2
         ierr=talsh_tensor_construct(oooo_tensor,C8,oooo_dims,init_val=ZERO)
         allocate (mo_list(4*nocc))
         nmo = nocc 
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,oooo_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         vovo_dims(1) = nvir * 2
         vovo_dims(2) = nocc * 2
         vovo_dims(3) = nvir * 2
         vovo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(vovo_tensor,C8,vovo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,vovo_tensor,0,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
         ovvo_dims(1) = nocc * 2
         ovvo_dims(2) = nvir * 2
         ovvo_dims(3) = nvir * 2
         ovvo_dims(4) = nocc * 2
         ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1) = nocc
         nmo(2) = nvir
         nmo(3) = nvir
         nmo(4) = nocc
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         mo_list(nvir+nocc+1:2*nvir+nocc)     = mo_vir(1:nvir)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,ovvo_tensor,0,nmo,mo_list)
         deallocate (mo_list)
!        Make anti-symmetrized vovo tensor
         minusone_dims(1) = 1
         ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",vovo_tensor,ovvo_tensor,minusone_tensor)
         ierr=talsh_tensor_destruct(ovvo_tensor)

!        Get vvvv tensor
         vvvv_dims = nvir * 2
         ierr=talsh_tensor_construct(vvvv_tensor,C8,vvvv_dims,init_val=ZERO)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vvvv_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get ooov tensor
         ooov_dims(1:3) = nocc * 2
         ooov_dims(4)   = nvir * 2
         ierr=talsh_tensor_construct(ooov_tensor,C8,ooov_dims,init_val=ZERO)
         allocate (mo_list(3*nocc+nvir))
         nmo(1:3) = nocc
         nmo(4)   = nvir
         mo_list(1:nocc)               = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)        = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc)      = mo_occ(1:nocc)
         mo_list(3*nocc+1:3*nocc+nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,ooov_tensor,12,nmo,mo_list)
         deallocate (mo_list)

!        Get vovv tensor
         vovv_dims(1)   = nvir * 2
         vovv_dims(2)   = nocc * 2
         vovv_dims(3:4) = nvir * 2
         ierr=talsh_tensor_construct(vovv_tensor,C8,vovv_dims,init_val=ZERO)
         allocate (mo_list(3*nvir+nocc))
         nmo(1)   = nvir
         nmo(2)   = nvir
         nmo(3)   = nocc
         nmo(4)   = nvir
         mo_list(1:nvir)                    = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)             = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)      = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:3*nvir+nocc) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,vovv_tensor,34,nmo,mo_list)
         deallocate (mo_list)

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)

       end subroutine get_CC_integrals

end module exacorr_cc
