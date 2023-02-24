module talsh_1rdm_2rdm_ccsd_kernels

! This module contains routines to calculate density matrices
! See Shee, A.; Visscher, L.; Saue, T. J. Chem. Phys. 2016, 145, 184107.
! and Halbert, L and Gomes, ASP, in preparation (ground-state 2RDM, to come EOMCC 1RDMs) 
! Written by Stan Papadopoulos, February 2019
! Refactored by Andre Gomes and Loic Halbert summer 2022

!  use talsh
!  use tensor_algebra
!  use talsh_common_routines, only : tensor_norm2
!  use exacorr_datatypes, only : talsh_dm_tens, exacc_input
!  use exacorr_utils, only     : print_date
  use, intrinsic:: ISO_C_BINDING

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_utils
  use intermediates
  use talsh_common_routines

  use exacorr_global

  implicit none

  complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0), &
                             MINUS_ONE_QUARTER=(-0.25D0,0.D0), ONE_QUARTER=(0.25D0,0.D0), &
                             ONE_EIGHT=(0.125D0,0.D0), MINUS_ONE_EIGHT=(-0.125D0,0.D0)

  private
   
    public talsh_calc_dm_1rdm_ccsd_GS_mo
    public talsh_calc_dm_1rdm_ccsd_EE_mo
    public talsh_calc_dm_1rdm_ccsd_IP_mo
    public talsh_calc_dm_1rdm_ccsd_EA_mo

    public calc_dmSandD_CCD_mo_bis

    public talsh_calc_dm_2rdm_ccsd_GS_mo
    public talsh_calc_dm_2rdm_ccsd_EE_mo
    public talsh_calc_dm_2rdm_ccsd_IP_mo
    public talsh_calc_dm_2rdm_ccsd_EA_mo


  contains

    subroutine talsh_calc_dm_1rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,print_level,t1_tensor,l1_tensor)

!     input variables
      type(talsh_dm_tens), intent(inout) :: gamma
      type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
      type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor
      integer,             intent(in   ) :: nocc, nvir, print_level
!     auxiliary tensors
      type(talsh_tens_t) :: one_tensor
      type(talsh_tens_t) :: tau2_tensor
      type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor
!     tensor dimensions
      integer(C_INT)              :: one_dims(1)=1
!     CCD switch
      logical :: CCD = .false.
!     error code
      integer :: ierr

      if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

!*********************************
!     construct auxiliary tensor *
!*********************************
 
      ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
      if (.not.CCD) then
        ierr=talsh_tensor_construct(tau2_tensor,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
        ierr=talsh_tensor_construct(vv_aux_tensor,C8,(/nvir,nvir/),init_val=ZERO)
        ierr=talsh_tensor_construct(oo_aux_tensor,C8,(/nocc,nocc/),init_val=ZERO)
      end if

      if (print_level.gt.8) then
        print *,'T1',tensor_norm2(t1_tensor)
        print *,'T2',tensor_norm2(t2_tensor)*0.25d0
        print *,'L1',tensor_norm2(l1_tensor)
        print *,'L2',tensor_norm2(l2_tensor)*0.25d0
      end if

!********************
!     calculate DMs *
!********************

!----------------------------------------------------------------
!     oo: g'(j,i) = - 1/2 t(ef,im) * l(jm,ef) - t(e,i) * l(j,e) |
!----------------------------------------------------------------
      ierr=talsh_tensor_init(gamma%prime_oo,ZERO)
      ierr=talsh_tensor_contract("G(j,i)+=T(e,f,i,m)*L(j,m,e,f)",gamma%prime_oo,t2_tensor,l2_tensor,scale=MINUS_ONE_HALF)
      if (print_level.gt.8) print *,'G(j,i) a',tensor_norm2(gamma%prime_oo)
      if (.not.CCD) ierr=talsh_tensor_contract("G(j,i)+=T(e,i)*L(j,e)",gamma%prime_oo,t1_tensor,l1_tensor,scale=MINUS_ONE)
      if (print_level.gt.8) print *,'G(j,i) b',tensor_norm2(gamma%prime_oo)

!--------------------------------------------------------------
!     vv: g'(a,b) = 1/2 t(ae,mn) * l(mn,be) + t(a,m) * l(m,b) |
!--------------------------------------------------------------
      ierr=talsh_tensor_init(gamma%prime_vv,ZERO)
      ierr=talsh_tensor_contract("G(a,b)+=T(a,e,m,n)*L(m,n,b,e)",gamma%prime_vv,t2_tensor,l2_tensor,scale=ONE_HALF)
      if (print_level.gt.8) print *,'G(a,b) a',tensor_norm2(gamma%prime_vv)
      if (.not.CCD) ierr=talsh_tensor_contract("G(a,b)+=T(a,m)*L(m,b)",gamma%prime_vv,t1_tensor,l1_tensor)
      if (print_level.gt.8) print *,'G(a,b) b',tensor_norm2(gamma%prime_vv)

!     no off-diagonal terms in CCD
      if (.not.CCD) then
!----------------------------
!     vo: g'(a,i) += t(a,i) |
!----------------------------
        ierr=talsh_tensor_init(gamma%prime_vo,ZERO)
        ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
        if (print_level.gt.8) print *,'G(a,i) a',tensor_norm2(gamma%prime_vo)

!-----------------------------------------------------------
!     vo: g'(a,i) += l(m,e) * (t(ae,im) - t(e,i) * t(a,m)) |
!-----------------------------------------------------------
        ierr=talsh_tensor_init(tau2_tensor,ZERO)
        ierr=talsh_tensor_add("T(a,e,i,m)+=R(a,e,i,m)",tau2_tensor,t2_tensor)
        if (print_level.gt.8) print *,'G(a,i) tau2 a',tensor_norm2(tau2_tensor)
        ierr=talsh_tensor_contract("T(a,e,i,m)+=T(e,i)*T(a,m)",tau2_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)
        if (print_level.gt.8) print *,'G(a,i) tau2 b',tensor_norm2(tau2_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,tau2_tensor)
        if (print_level.gt.8) print *,'G(a,i) b',tensor_norm2(gamma%prime_vo)

!----------------------------------------------------------------------------
!     vo: g'(a,i) -= 1/2 l(mn,ef) * (t(ef,in) * t(a,m) + t(e,i) * t(af,mn)) |
!----------------------------------------------------------------------------
        ierr=talsh_tensor_init(oo_aux_tensor,ZERO)
        ierr=talsh_tensor_init(vv_aux_tensor,ZERO)
        ierr=talsh_tensor_contract("X(m,i)+=L(m,n,e,f)*T(e,f,i,n)",oo_aux_tensor,l2_tensor,t2_tensor)
        if (print_level.gt.8) print *,'G(a,i) aux oo',tensor_norm2(oo_aux_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF) 
        if (print_level.gt.8) print *,'G(a,i) term c',tensor_norm2(gamma%prime_vo)
        ierr=talsh_tensor_contract("X(a,e)+=L(m,n,e,f)*T(a,f,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
        if (print_level.gt.8) print *,'G(a,i) aux vv',tensor_norm2(vv_aux_tensor)
        ierr=talsh_tensor_contract("G(a,i)+=X(a,e)*T(e,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
        if (print_level.gt.8) print *,'G(a,i) term d',tensor_norm2(gamma%prime_vo)

!---------------------------
!     ov: g'(i,a) = l(i,a) |
!---------------------------
        ierr=talsh_tensor_init(gamma%prime_ov,ZERO)
        ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gamma%prime_ov,l1_tensor)
        if (print_level.gt.8) print *,'G(i,a) term a',tensor_norm2(gamma%prime_ov)

!--------------------------------------------------
!       ov symm: g(i,a) = 1/2 (g'(i,a) + g'(a,i)) |
!--------------------------------------------------
        ierr=talsh_tensor_init(gamma%symm_ov,ZERO)
        ierr=talsh_tensor_contract("G(i,a)+=X(i,a)*K()",gamma%symm_ov,gamma%prime_ov,one_tensor,scale=ONE_HALF)
        ierr=talsh_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
        if (print_level.gt.8) print *,'G(i,a) symmetrized',tensor_norm2(gamma%symm_ov)
        

!--------------------------------------------------
!       vo symm: g(a,i) = symmetrized c.c. g(i,a)
!--------------------------------------------------
         ierr=talsh_tensor_init(gamma%symm_vo,ZERO)
         ierr=talsh_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%symm_ov,one_tensor)
         if (print_level.gt.8) print *,'G(a,i) symmetrized',tensor_norm2(gamma%symm_vo)

      end if

!------------------------------------------------
!     vv symm: g(b,a) = 1/2 (g'(b,a) + g'(a,b)) |
!------------------------------------------------
      ierr=talsh_tensor_init(gamma%symm_vv,ZERO)
      ierr=talsh_tensor_contract("G(b,a)+=X(b,a)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(b,a)+=X+(a,b)*K()",gamma%symm_vv,gamma%prime_vv,one_tensor,scale=ONE_HALF)
      if (print_level.gt.8) print *,'G(b,a) symmetrized',tensor_norm2(gamma%symm_vv)

!------------------------------------------------
!     oo symm: g(j,i) = 1/2 (g'(j,i) + g'(i,j)) |
!------------------------------------------------
      ierr=talsh_tensor_init(gamma%symm_oo,ZERO)
      ierr=talsh_tensor_contract("G(j,i)+=X(j,i)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)
      ierr=talsh_tensor_contract("G(j,i)+=X+(i,j)*K()",gamma%symm_oo,gamma%prime_oo,one_tensor,scale=ONE_HALF)
      if (print_level.gt.8) print *,'G(j,i) symmetrized',tensor_norm2(gamma%symm_oo)


!     debug
      if (print_level.gt.8) then
        if (.not.CCD) then
          print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
          print*, "gsymm_ov   = ", tensor_norm2(gamma%symm_vo)
          print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
          print*, "gsymm_vo   = ", tensor_norm2(gamma%symm_vo)
        end if
        print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)
        print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
        print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
        print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)

        if (print_level.gt.9) then
          if (.not.CCD) then
            print *,'>>>>> Loic_ In my branch no  relccsd order Allowed'
            print *,'>>>>> elements, DOV not symmetrized, relccsd order'
!            call print_tensor_elements(gamma%prime_ov)!,relccsd_order=.true.)
            print *,'<<<<< elements, DOV not symmetrized, relccsd order'
            print *,'>>>>> elements, DVO not symmetrized, relccsd order'
!            call print_tensor_elements(gamma%prime_vo)!,relccsd_order=.true.)
            print *,'<<<<< elements, DVO not symmetrized, relccsd_order'
          end if
          print *,'>>>>> elements, DVV not symmetrized, relccsd order'
!          call print_tensor_elements(gamma%prime_vv)!,relccsd_order=.true.)
          print *,'<<<<< elements, DVV not symmetrized, relccsd order'
          print *,'>>>>> elements, DOO not symmetrized, relccsd order'
!          call print_tensor_elements(gamma%prime_oo)!,relccsd_order=.true.)
          print *,'<<<<< elements, DOO not symmetrized, relccsd order'

          if (.not.CCD) then
            print *,'>>>>> elements, DOV symmetrized'
!            call print_tensor_elements(gamma%symm_ov)!,relccsd_order=.true.)
            print *,'<<<<< elements, DOV symmetrized'
          end if
          print *,'>>>>> elements, DOO symmetrized'
!          call print_tensor_elements(gamma%symm_oo)!,relccsd_order=.true.)
          print *,'<<<<< elements, DOO symmetrized'
          print *,'>>>>> elements, DVV symmetrized'
!          call print_tensor_elements(gamma%symm_vv)!,relccsd_order=.true.)
          print *,'<<<<< elements, DVV symmetrized'
          if (.not.CCD) then
            print *,'>>>>> elements, DVO symmetrized'
!            call print_tensor_elements(gamma%symm_vo)!,relccsd_order=.true.)
            print *,'<<<<< elements, DVO symmetrized'
          end if
        end if ! extended print level, 9
      end if ! extended print level,8


!**************
!     cleanup *
!**************

      ierr=talsh_tensor_destruct(one_tensor)
      if (.not.CCD) then
        ierr=talsh_tensor_destruct(tau2_tensor)
        ierr=talsh_tensor_destruct(vv_aux_tensor)
        ierr=talsh_tensor_destruct(oo_aux_tensor)
      end if

    end subroutine talsh_calc_dm_1rdm_ccsd_GS_mo


    subroutine talsh_calc_dm_2rdm_ccsd_GS_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,print_level,t1_tensor,l1_tensor)

      type(talsh_dm_tens), intent(inout) :: gamma

      type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
      type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor
      integer,             intent(in   ) :: nocc, nvir, print_level

      type(talsh_tens_t) :: one_tensor
!     tensor dimensions
      integer(C_INT)              :: one_dims(1)=1

      type(talsh_tens_t) :: Hinter_vovv
      type(talsh_tens_t) :: Hinter_ooov, Hinter_ovvo
      type(talsh_tens_t) :: Hinter_oooo, Hinter_vvoo


!     CCD switch
      logical :: CCD = .false.
!     error code
      integer :: ierr

      logical :: debug=.true.
      ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)


      ierr=talsh_tensor_construct(Hinter_vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_ovvo,C8,(/nocc,nvir,nvir,nocc/),init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_vvoo,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
      ! gamma oovv
      ! j1
      ierr=talsh_tensor_add("G(i,j,a,b)+=L(i,j,a,b)",gamma%prime_oovv,l2_tensor)
      if (ierr.ne.0)  print*, " j1, ierr  ", ierr

      ! gamma vovv
      ! h2
      ierr=talsh_tensor_contract("G(a,i,b,c)+=L(m,i,b,c)*T(a,m)",gamma%prime_vovv,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " h2, ierr  ", ierr

      !gamma ooov
      ! e2
      ierr=talsh_tensor_contract("G(i,j,k,a)+=L(i,j,e,a)*T(e,k)",gamma%prime_ooov,l2_tensor,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0) print*, " e2, ierr  ", ierr

      ! gamma vvvv
      ! i3
      ierr=talsh_tensor_contract("G(a,b,c,d)+=L(m,n,c,d)*T(a,b,m,n)",gamma%prime_vvvv,l2_tensor,t2_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " i3, ierr  ", ierr
      ! i4
      ierr=talsh_tensor_init(Hinter_vovv)
      ierr=talsh_tensor_contract("H(a,n,c,d)+=L(m,n,c,d)*T(a,m)",Hinter_vovv,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " i4-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,d)+=H(a,n,c,d)*T(b,n)",gamma%prime_vvvv,Hinter_vovv,t1_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " i4-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,d)+=H(b,n,c,d)*T(a,n)",gamma%prime_vvvv,Hinter_vovv,t1_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " i4-c, ierr  ", ierr

      ! gamma oooo
      ! n3
      ierr=talsh_tensor_contract("G(i,j,k,l)+=L(i,j,e,f)*T(e,f,k,l)",gamma%prime_oooo,l2_tensor,t2_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " n3, ierr  ", ierr
      ! n4
      ierr=talsh_tensor_init(Hinter_ooov)
      ierr=talsh_tensor_contract("H(i,j,k,f)+=L(i,j,e,f)*T(e,k)",Hinter_ooov,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " n4-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(i,j,k,l)+=H(i,j,k,f)*T(f,l)",gamma%prime_oooo,Hinter_ooov,t1_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " n4-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(i,j,k,l)+=H(i,j,l,f)*T(f,k)",gamma%prime_oooo,Hinter_ooov,t1_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " n4-c, ierr  ", ierr

      ! Start ovov f7 intermediate for g16 and k16
      ! f7
      ierr=talsh_tensor_contract("G(i,a,j,b)+=L(i,m,e,b)*T(e,a,j,m)",gamma%prime_ovov,l2_tensor,t2_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " f7, ierr  ", ierr
      ! gamma vvvo
      ! g16 sign of f7 and inversion
      ierr=talsh_tensor_contract("G(a,b,c,i)+=H(j,b,i,c)*T(a,j)",gamma%prime_vvvo,gamma%prime_ovov,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " g16-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,i)+=H(j,a,i,c)*T(b,j)",gamma%prime_vvvo,gamma%prime_ovov,t1_tensor)
      if (ierr.ne.0)  print*, " g16-b, ierr  ", ierr

      ! g15 + g18
      ierr=talsh_tensor_contract("G(a,b,c,i)+=H(a,c)*T(b,i)",gamma%prime_vvvo,gamma%prime_vv,t1_tensor)
      if (ierr.ne.0)  print*, "g15+g18-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,i)+=H(b,c)*T(a,i)",gamma%prime_vvvo,gamma%prime_vv,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, "g15+g18-b, ierr  ", ierr

      ! g17 + g19
      ierr=talsh_tensor_contract("G(a,b,c,i)+=H(a,b,c,d)*T(d,i)",gamma%prime_vvvo,gamma%prime_vvvv,t1_tensor,scale=MINUS_ONE)

      ! g14
      ierr=talsh_tensor_contract("G(a,b,c,i)+=L(m,c)*T(a,b,m,i)",gamma%prime_vvvo,l1_tensor,t2_tensor)
      if (ierr.ne.0)  print*, " g14, ierr  ", ierr

      ! gamma ovoo
      ! k16  sign of f7
      ierr=talsh_tensor_contract("G(i,a,j,k)+=H(i,a,j,e)*T(e,k)",gamma%prime_ovoo,gamma%prime_ovov,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " k16-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(i,a,j,k)+=H(i,a,k,e)*T(e,j)",gamma%prime_ovoo,gamma%prime_ovov,t1_tensor)
      if (ierr.ne.0)  print*, " k16-b, ierr  ", ierr

      ! k15+k18
      ierr=talsh_tensor_contract("G(i,a,j,k)+=H(i,j)*T(a,k)",gamma%prime_ovoo,gamma%prime_oo,t1_tensor)
      if (ierr.ne.0)  print*, " k15+k18-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(i,a,j,k)+=H(i,k)*T(a,j)",gamma%prime_ovoo,gamma%prime_oo,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " k15+k18-b, ierr  ", ierr

      ! k17+k19
      ierr=talsh_tensor_contract("G(i,a,j,k)+=H(i,l,j,k)*T(a,l)",gamma%prime_ovoo,gamma%prime_oooo,t1_tensor)
      if (ierr.ne.0)  print*, " k19-d, ierr  ", ierr

      ! k14
      ierr=talsh_tensor_contract("G(i,a,j,k)+=L(i,e)*T(e,a,j,k)",gamma%prime_ovoo,l1_tensor,t2_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " k14, ierr  ", ierr

      ! End of gamma ovov as f7 was earlier construct
      ! f6
      ierr=talsh_tensor_contract("G(i,a,j,b)+=L(i,b)*T(a,j)",gamma%prime_ovov,l1_tensor,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " f6, ierr  ", ierr
      ! f8
      ierr=talsh_tensor_init(Hinter_ooov)
      ierr=talsh_tensor_contract("H(i,m,j,b)+=L(i,m,e,b)*T(e,j)",Hinter_ooov,l2_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " f8-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(i,a,j,b)+=H(i,m,j,b)*T(a,m)",gamma%prime_ovov,Hinter_ooov,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " f8-b, ierr  ", ierr

      ! gamma vvoo
      ! l1
      ierr=talsh_tensor_add("G(a,b,i,j)+=T(a,b,i,j)",gamma%prime_vvoo,t2_tensor)
      if (ierr.ne.0)  print*, " l1-a, ierr  ", ierr
 
      !m5+m9
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,i)*T(a,b,n,j)",gamma%prime_vvoo,gamma%prime_oo,t2_tensor)
      if (ierr.ne.0)  print*, " m9-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,j)*T(a,b,n,i)",gamma%prime_vvoo,gamma%prime_oo,t2_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m9-c, ierr  ", ierr

      !m7+m11
      if (ierr.ne.0)  print*, " m11-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(a,e)*T(e,b,i,j)",gamma%prime_vvoo,gamma%prime_vv,t2_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m11-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(b,e)*T(e,a,i,j)",gamma%prime_vvoo,gamma%prime_vv,t2_tensor)
      if (ierr.ne.0)  print*, " m11-c, ierr  ", ierr

      ! m20 
      ierr=talsh_tensor_init(Hinter_ovvo)
      ierr=talsh_tensor_contract("H(n,a,f,i)+=L(m,n,e,f)*T(a,e,i,m)",Hinter_ovvo,l2_tensor,t2_tensor)
      if (ierr.ne.0)  print*, " m20-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,i)*T(b,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " m20-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,j)*T(b,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " m20-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,i)*T(a,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " m20-d, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,j)*T(a,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " m20-e, ierr  ", ierr
      ! m21
      ierr=talsh_tensor_contract("A(f,b,j,n)+=T(f,j)*T(b,n)",Hinter_vvoo,t1_tensor,t1_tensor)
      if (ierr.ne.0)  print*, " m21-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,i)*A(f,b,j,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m21-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,j)*A(f,b,i,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo)
      if (ierr.ne.0)  print*, " m21-d, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,i)*A(f,a,j,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo)
      if (ierr.ne.0)  print*, " m21-e, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,j)*A(f,a,i,n)",gamma%prime_vvoo,Hinter_ovvo,Hinter_vvoo,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m21-f, ierr  ", ierr

      ! m3+m4+m6+m8 = Pij Pab \ts bj \gammas vo -l2
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(a,i)*T(b,j)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor)
      if (ierr.ne.0)  print*, " m8-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(a,j)*T(b,i)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m8-d, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(b,i)*T(a,j)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, " m8-e, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(b,j)*T(a,i)",gamma%prime_vvoo,gamma%prime_vo,t1_tensor)
      if (ierr.ne.0)  print*, " m8-f, ierr  ", ierr

      ! MINUS l2
      ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " l2-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " l2-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " l2-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,j)*T(a,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
      if (ierr.ne.0)  print*, " l2-d, ierr  ", ierr
      !m12+m14+m15+m24
      ierr=talsh_Tensor_init(Hinter_vvoo)
      ierr=talsh_tensor_contract("H(a,b,i,j)+=T(a,b,i,j)*K()",Hinter_vvoo,t2_tensor,one_tensor)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -a, ierr  ", ierr
      ierr=talsh_tensor_contract("H(a,b,i,j)+=T(a,i)*T(b,j)",Hinter_vvoo,t1_tensor,t1_tensor)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -b, ierr  ", ierr
      ierr=talsh_tensor_contract("H(a,b,i,j)+=T(a,j)*T(b,i)",Hinter_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -c, ierr  ", ierr

      ierr=talsh_Tensor_init(Hinter_oooo)
      ierr=talsh_tensor_contract("H(m,n,i,j)+=L(m,n,e,f)*T(e,f,i,j)",Hinter_oooo,l2_tensor,Hinter_vvoo)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -d, ierr  ", ierr

      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(m,n,i,j)*A(a,b,m,n)",gamma%prime_vvoo,Hinter_oooo,Hinter_vvoo,scale=ONE_QUARTER)
      if (ierr.ne.0)  print*, "m12+m14+m15+m24 -e, ierr  ", ierr

      !! gamma%prime construction


      !! Start Symetrisation
      ! gamma%symm_ooov 
      ierr=talsh_tensor_contract("G(i,j,k,a)+=G(i,j,k,a)*K()",gamma%symm_ooov,gamma%prime_ooov,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ooov 1-a  ", ierr
      ierr=talsh_tensor_contract("G(i,j,k,a)+=G+(k,a,i,j)*K()",gamma%symm_ooov,gamma%prime_ovoo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ooov 1-b  ", ierr
      ! gamma%symm_oovo = -gamma%symm_ooov
      ierr=talsh_tensor_contract("G(i,j,a,k)+=G(i,j,k,a)*K()",gamma%symm_oovo,gamma%symm_ooov,one_tensor,scale=MINUS_ONE)

      ! gamma%symm_ovoo 
      ierr=talsh_tensor_contract("G(i,a,k,l)+=G(i,a,k,l)*K()",gamma%symm_ovoo,gamma%prime_ovoo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovoo 2-a  ", ierr
      ierr=talsh_tensor_contract("G(i,a,k,l)+=G+(k,l,i,a)*K()",gamma%symm_ovoo,gamma%prime_ooov,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovoo 2-b  ", ierr
      ! gamma%symm_vooo = -gamma%symm_ovoo
      ierr=talsh_tensor_contract("G(a,i,k,l)+=G(i,a,k,l)*K()",gamma%symm_vooo,gamma%symm_ovoo,one_tensor,scale=MINUS_ONE)

      ! gamma%symm_ovov 
      ierr=talsh_tensor_contract("G(i,a,j,b)+=G(i,a,j,b)*K()",gamma%symm_ovov,gamma%prime_ovov,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovov 3-a  ", ierr
      ierr=talsh_tensor_contract("G(i,a,j,b)+=G+(j,b,i,a)*K()",gamma%symm_ovov,gamma%prime_ovov,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm ovov 3-b  ", ierr
      ! gamma%symm_vovo = gamma%symm_ovov  
      ierr=talsh_tensor_contract("G(a,i,b,j)+=G(i,a,j,b)*K()",gamma%symm_vovo,gamma%symm_ovov,one_tensor)
      ! gamma%symm_voov = -gamma%symm_ovov  
      ierr=talsh_tensor_contract("G(a,i,j,b)+=G(i,a,j,b)*K()",gamma%symm_voov,gamma%symm_ovov,one_tensor,scale=MINUS_ONE)
      ! gamma%symm_ovvo = -gamma%symm_ovov  
      ierr=talsh_tensor_contract("G(i,a,b,j)+=G(i,a,j,b)*K()",gamma%symm_ovvo,gamma%symm_ovov,one_tensor,scale=MINUS_ONE)

      ! gamma%symm_vvvo 
      ierr=talsh_tensor_contract("G(a,b,c,i)+=G(a,b,c,i)*K()",gamma%symm_vvvo,gamma%prime_vvvo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvo 4-a  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,i)+=G+(c,i,a,b)*K()",gamma%symm_vvvo,gamma%prime_vovv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvo 4-b  ", ierr
      ! gamma%symm_vvov = -gamma%symm_vvvo
      ierr=talsh_tensor_contract("G(a,b,i,c)+=G(a,b,c,i)*K()",gamma%symm_vvov,gamma%symm_vvvo,one_tensor,scale=MINUS_ONE)

      ! gamma%symm_vovv
      ierr=talsh_tensor_contract("G(a,i,b,c)+=G(a,i,b,c)*K()",gamma%symm_vovv,gamma%prime_vovv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vovv 5-a  ", ierr
      ierr=talsh_tensor_contract("G(a,i,b,c)+=G+(b,c,a,i)*K()",gamma%symm_vovv,gamma%prime_vvvo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vovv 5-b  ", ierr
      ! gamma%symm_ovvv = -gamma%symm_vovv
      ierr=talsh_tensor_contract("G(i,a,b,c)+=G(a,i,b,c)*K()",gamma%symm_ovvv,gamma%symm_vovv,one_tensor,scale=MINUS_ONE)

      ! gamma%symm_vvvv
      ierr=talsh_tensor_contract("G(a,b,c,d)+=G(a,b,c,d)*K()",gamma%symm_vvvv,gamma%prime_vvvv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvv 6-a  ", ierr
      ierr=talsh_tensor_contract("G(a,b,c,d)+=G+(c,d,a,b)*K()",gamma%symm_vvvv,gamma%prime_vvvv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvvv 6-b  ", ierr

      ! gamma%symm_oooo
      ierr=talsh_tensor_contract("G(i,j,k,l)+=G(i,j,k,l)*K()",gamma%symm_oooo,gamma%prime_oooo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm oooo 7-a  ", ierr
      ierr=talsh_tensor_contract("G(i,j,k,l)+=G+(k,l,i,j)*K()",gamma%symm_oooo,gamma%prime_oooo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm oooo 7-b  ", ierr


      ! gamma%symm_oovv 
      ierr=talsh_tensor_contract("G(i,j,a,b)+=G(i,j,a,b)*K()",gamma%symm_oovv,gamma%prime_oovv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm oovv 8-a  ", ierr
      ierr=talsh_tensor_contract("G(i,j,a,b)+=G+(a,b,i,j)*K()",gamma%symm_oovv,gamma%prime_vvoo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm oovv 8-b  ", ierr

      ! gamma%symm_vvoo 
      ierr=talsh_tensor_contract("G(a,b,i,j)+=G(a,b,i,j)*K()",gamma%symm_vvoo,gamma%prime_vvoo,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvoo 9-a  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=G+(i,j,a,b)*K()",gamma%symm_vvoo,gamma%prime_oovv,one_tensor,scale=ONE_HALF)
      if (ierr.ne.0)  print*, " symm vvoo 9-b  ", ierr



      ierr=talsh_tensor_destruct(Hinter_vovv)
      ierr=talsh_tensor_destruct(Hinter_ooov)
      ierr=talsh_tensor_destruct(Hinter_ovvo)
      ierr=talsh_tensor_destruct(Hinter_oooo)
      ierr=talsh_tensor_destruct(Hinter_vvoo)





    end subroutine talsh_calc_dm_2rdm_ccsd_GS_mo  


    subroutine calc_dmSandD_CCD_mo_bis(gamma,t1_tensor,t2_tensor,l1_tensor,l2_tensor,nocc,nvir,print_level,int_t)

        use talsh
        use exacorr_datatypes !, only :talsh_intg_tens,talsh_comm_tens,talsh_lambda_tens
 !       use exacorr_global

!        Common tensors
!      type(talsh_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron
!        integrals)
      type(talsh_intg_tens_t),intent(inout) :: int_t


      type(talsh_dm_tens), intent(inout) :: gamma

      type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
      type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor
      integer,             intent(in   ) :: nocc, nvir, print_level

      logical  :: debug = .true.

!     tensor dimensions
      integer(C_INT)              :: one_dims(1)=1

      integer                     :: i,j  

      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

!      integer(INTD), dimension(4) :: vvoo_dims
      type(talsh_tens_t) :: aux_tensor_vvoo

      type(talsh_tens_t) :: Hinter_oooo, Hinter_ovvo
      type(talsh_tens_t) :: Hinter_vv, Hinter_oo

      integer :: ierr

      type(C_PTR):: body_p
      complex(8), pointer :: result_tens(:)  
      type(talsh_tens_t) :: result_tensor,one_tensor
      integer(INTD)      :: result_dims(1)

!        Initialize scalars that are to be used as tensors in contractions
      ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)


!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)

      ierr=talsh_tensor_construct(Hinter_oo,C8,oo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_vv,C8,vv_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_oooo,C8,oooo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(Hinter_ovvo,C8,ovvo_dims,init_val=ZERO)
      ierr=talsh_tensor_construct(aux_tensor_vvoo,C8,vvoo_dims,init_val=ZERO)

!      Eq1-1 Gij=-\undemi \ld jmef \td efim 
      ierr=talsh_tensor_contract("G(i,j)+=L(j,m,e,f)*T(e,f,i,m)",gamma%prime_oo,l2_tensor,t2_tensor,scale=MINUS_ONE_HALF)
      if (debug)  print*, " Eq1-1, ierr  ", ierr

!      Eq1-2 Gab=+\undemi \ld mnae \td bemn
      ierr=talsh_tensor_contract("G(a,b)+=L(m,n,a,e)*T(b,e,m,n)",gamma%prime_vv,l2_tensor,t2_tensor,scale=ONE_HALF)
      if (debug)  print*, " Eq1-2, ierr  ", ierr

!      Eq2-1 Gijkl=+ ?\unhuit? \ld klef \td efij
     ! gamma oooo
      ! n3
      ierr=talsh_tensor_contract("G(i,j,k,l)+=L(i,j,e,f)*T(e,f,k,l)",gamma%prime_oooo,l2_tensor,t2_tensor,scale=ONE_HALF)
      if (debug)  print*, " n3, ierr  ", ierr

!      Eq2-2 Gabcd = + ?\unhuit? \ld mnab \td cdmn
      ! gamma vvvv
      ! i3
      ierr=talsh_tensor_contract("G(a,b,c,d)+=L(m,n,c,d)*T(a,b,m,n)",gamma%prime_vvvv,l2_tensor,t2_tensor,scale=ONE_HALF)
      if (debug)  print*, " i3, ierr  ", ierr

!      Eq2-3 Gajib = + ?\unquart? \ld imae \td bejm
      ! f7
      ierr=talsh_tensor_contract("G(i,a,j,b)+=L(i,m,e,b)*T(e,a,j,m)",gamma%prime_ovov,l2_tensor,t2_tensor,scale=MINUS_ONE)
      if (debug)  print*, " f7, ierr  ", ierr
  
!      Eq2-4 Gabij = + ?\unquart? \ld ijab
      ! j1
      ierr=talsh_tensor_add("G(i,j,a,b)+=L(i,j,a,b)",gamma%prime_oovv,l2_tensor)
      if (debug)  print*, " j1, ierr  ", ierr

!      Eq2-5 Goovv
!      Eq2-5-a : Gijab = + ?\unquart? \td abij
      ! l1
      ierr=talsh_tensor_add("G(a,b,i,j)+=T(a,b,i,j)",gamma%prime_vvoo,t2_tensor)
      if (debug)  print*, " l1-a, ierr  ", ierr
!      Eq2-5-b : Gijab = + ?\unseize? \ld mnef \td efij \td abmn         
      ! m12
      ierr=talsh_tensor_init(Hinter_oooo)
      ierr=talsh_tensor_contract("H(m,n,i,j)+=L(m,n,e,f)*T(e,f,i,j)",Hinter_oooo,l2_tensor,t2_tensor)
      if (debug)  print*, " m12-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(m,n,i,j)*T(a,b,m,n)",gamma%prime_vvoo,Hinter_oooo,t2_tensor,scale=ONE_QUARTER)
      if (debug)  print*, " m12-b, ierr  ", ierr
!      Eq2-5-c : Gijab = ?\undemi? \ld mnef \td aeim \td bfjm
      ! m20
      ierr=talsh_tensor_init(Hinter_ovvo)
      ierr=talsh_tensor_contract("H(n,a,f,i)+=L(m,n,e,f)*T(a,e,i,m)",Hinter_ovvo,l2_tensor,t2_tensor)
      if (debug)  print*, " m20-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,i)*T(b,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor)
      if (debug)  print*, " m20-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,a,f,j)*T(b,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=MINUS_ONE)
      if (debug)  print*, " m20-c, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,i)*T(a,f,j,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor,scale=MINUS_ONE)
      if (debug)  print*, " m20-d, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,b,f,j)*T(a,f,i,n)",gamma%prime_vvoo,Hinter_ovvo,t2_tensor)
      if (debug)  print*, " m20-e, ierr  ", ierr
!      Eq2-5-d : Gjiab = ?-\unquart? \ld menf \td ebmn \td afij       
      ! m7
      ierr=talsh_tensor_init(Hinter_vv)
      ierr=talsh_tensor_contract("H(a,f)+=L(m,n,e,f)*T(e,a,m,n)",Hinter_vv,l2_tensor,t2_tensor)
      if (debug)  print*, " m7-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(a,f)*T(f,b,i,j)",gamma%prime_vvoo,Hinter_vv,t2_tensor,scale=MINUS_ONE_HALF)
      if (debug)  print*, " m7-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(b,f)*T(f,a,i,j)",gamma%prime_vvoo,Hinter_vv,t2_tensor,scale=ONE_HALF)
      if (debug)  print*, " m7-c, ierr  ", ierr
!      Eq2-5-e : Gijab = ?-\unquart? \ld mnef \td efmj \td abin
      ! m9
      ierr=talsh_tensor_init(Hinter_oo)  !Loic_ Take care 
      ierr=talsh_tensor_contract("H(n,i)+=L(m,n,e,f)*T(e,f,m,i)",Hinter_oo,l2_tensor,t2_tensor)
      if (debug)  print*, " m9-a, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,i)*T(a,b,n,j)",gamma%prime_vvoo,Hinter_oo,t2_tensor,scale=MINUS_ONE_HALF)
      if (debug)  print*, " m9-b, ierr  ", ierr
      ierr=talsh_tensor_contract("G(a,b,i,j)+=H(n,j)*T(a,b,n,i)",gamma%prime_vvoo,Hinter_oo,t2_tensor,scale=ONE_HALF)
      if (debug)  print*, " m9-c, ierr  ", ierr


      print *,"CCD Results"   
!      call print_tensor_elements(gamma%prime_oo)  
!      call print_tensor_elements(int_t%foo)

      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_contract("E()+=G(i,j)*F(j,i)",result_tensor,gamma%prime_oo,int_t%foo)
      if (debug) print*, "result_tensor_oo ,ierr ", ierr
      if (debug) print*, "result_tensor_oo ,result_tensor_oo ", result_tens(1)
      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_contract("E()+=G(a,b)*F(b,a)",result_tensor,gamma%prime_vv,int_t%fvv)
      if (debug) print*, "result_tensor_aa ,ierr ", ierr
      if (debug) print*, "result_tensor_aa ,result_tensor_aa ", result_tens(1)
      ierr=talsh_tensor_init(result_tensor)

      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_contract("E()+=G(a,b,c,d)*V(c,d,a,b)",result_tensor,gamma%prime_vvvv,int_t%vvvv)
      if (debug) print*, "result_tensor_vvvv ,ierr ", ierr
      if (debug) print*, "result_tensor_vvvv ,result_tensor_vvvv ", result_tens(1)

      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_contract("E()+=G(i,j,k,l)*V(k,l,i,j)",result_tensor,gamma%prime_oooo,int_t%oooo)
      if (debug) print*, "result_tensor_oooo ,ierr ", ierr
      if (debug) print*, "result_tensor_oooo ,result_tensor_oooo ", result_tens(1)


      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_init(aux_tensor_vvoo)
      ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",aux_tensor_vvoo,int_t%oovv,one_tensor) 
      ierr=talsh_tensor_contract("E()+=G(i,j,a,b)*V(a,b,i,j)",result_tensor,gamma%prime_oovv,aux_tensor_vvoo)
      if (debug) print*, "result_tensor_oovv ,ierr ", ierr
      if (debug) print*, "result_tensor_oovv ,result_tensor_oovv ", result_tens(1)

      ierr=talsh_tensor_init(result_tensor)
      ierr=talsh_tensor_contract("E()+=G(a,b,i,j)*V(i,j,a,b)",result_tensor,gamma%prime_vvoo,int_t%oovv)
      if (debug) print*, "result_tensor_vvoo ,ierr ", ierr
      if (debug) print*, "result_tensor_vvoo ,result_tensor_vvoo ", result_tens(1)

      ierr=talsh_tensor_init(result_tensor)
!      ierr=talsh_tensor_init(aux_tensor_vvoo)
!      ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",aux_tensor_vvoo,int_t%oovv,one_tensor) 
      ierr=talsh_tensor_contract("E()+=G(i,a,j,b)*V(b,j,a,i)",result_tensor,gamma%prime_ovov,int_t%vovo)
      if (debug) print*, "result_tensor_ovov ,ierr ", ierr
      if (debug) print*, "result_tensor_ovov ,result_tensor_ovov ",result_tens(1)


      ierr=talsh_tensor_destruct(Hinter_oooo)
      ierr=talsh_tensor_destruct(Hinter_ovvo)
      ierr=talsh_tensor_destruct(Hinter_vv)
      ierr=talsh_tensor_destruct(Hinter_oo)
      ierr=talsh_tensor_destruct(aux_tensor_vvoo)

    end subroutine calc_dmSandD_CCD_mo_bis

        subroutine talsh_calc_dm_1rdm_ccsd_EE_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)


          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor

    !     tensor dimensions
          integer :: one_dims(1)=1
          integer, dimension(2) :: oo_dims, vv_dims
          integer, dimension(2) :: ov_dims,vo_dims
          integer, dimension(4) :: vvoo_dims, ovvo_dims, ovvv_dims, oovo_dims
          integer, dimension(4) :: oooo_dims,vvvv_dims,oovv_dims !Loic_
          integer, dimension(4) :: vvvo_dims, ooov_dims, vooo_dims, vovv_dims
          integer, dimension(4) :: vovo_dims, ovoo_dims, ovov_dims, voov_dims
          integer, dimension(4) :: vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)

    !*********************************
    !     construct auxiliary tensor *
    !*********************************
          ierr=talsh_tensor_construct(gamma%prime_oo,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_ov,C8,ov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vo,C8,vo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vv,C8,vv_dims,init_val=ZERO)

          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
          ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)

    !********************
    ! Intermediates     *
    !********************

            ierr=talsh_tensor_construct(MATone_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauUP_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauDown_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauNo_tensor,C8,vvoo_dims,init_val=ZERO)

            !  Loic_Tenseurs_Intermediaires


           !  MATone_oo
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,j)*L(i,e)",MATone_oo_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,f,j,m)*L(i,m,e,f)",MATone_oo_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,j)*L(i,e)",MATtwo_oo_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,f,j,m)*L(i,m,e,f)",MATtwo_oo_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(a,m)*L(m,b)",MATone_vv_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(e,a,m,n)*L(m,n,e,b)",MATone_vv_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(a,m)*L(m,b)",MATtwo_vv_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(e,a,m,n)*L(m,n,e,b)",MATtwo_vv_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           !  Tau
           
           ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauUp(a,b,i,j)+=T(a,b,i,j)",TauUp_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(a,i)*T(b,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(b,i)*T(a,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

    !********************
    !     calculate DMs *
    !********************


    ! rho ij
           ! eq a1 + a3 :: - \MATone ij   
           ierr=talsh_tensor_add("G(i,j)+=X(i,j)",gamma%prime_oo,MATone_oo_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr1 =",ierr
           ! eq a5      :: - \ld imef \rs fm \ts ej
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(i,e)+=L(i,m,e,f)*R(f,m)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr2 =",ierr
           ierr=talsh_tensor_contract("G(i,j)+=X(i,e)*T(e,j)",gamma%prime_oo,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr3 =",ierr
           ! eq a2 + a4 :: - \ro \MATtwo ij
           ierr=talsh_tensor_add("G(i,j)+=X(i,j)",gamma%prime_oo,MATtwo_oo_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr4 =",ierr

    !
    !   TAKE CARE no \kronecker_{ij}
    !


    ! rho ab
           ! eq b3 + b4 :: +MATone ab
           ierr=talsh_tensor_add("G(a,b)+=X(a,b)",gamma%prime_vv,MATone_vv_tensor)
           if(ierr.ne.0) write(*,*) "ierr5 =",ierr
           ! eq b5      :: +\ld mneb \rs em \ts an
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(n,b)+=L(m,n,e,b)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr6 =",ierr
           ierr=talsh_tensor_contract("G(a,b)+=X(n,b)*T(a,n)",gamma%prime_vv,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr7 =",ierr
           ! eq b1 + b2 :: + \ro \MATtwo ab
           ierr=talsh_tensor_add("G(a,b)+=X(a,b)",gamma%prime_vv,MATtwo_vv_tensor)
           if(ierr.ne.0) write(*,*) "ierr8 =",ierr


    ! rho ai
           ! eq c1     :: +\ro \ls me \td aeim
           ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*T(a,e,i,m)",gamma%prime_vo,l1_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr9 =",ierr
           ! eq c2 + c3   :: -\ro \ts am \MATtwo mi
           ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,MATtwo_oo_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr10 =",ierr
           ! eq c4        :: -\ro 1/2 \ld mnef \td eamn \ts fi
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr11 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*T(f,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr12 =",ierr
           ! eq c7 + c11  :: -\MATone ni \ts an 
           ierr=talsh_tensor_contract("G(a,i)+=X(n,i)*T(a,n)",gamma%prime_vo,MATone_oo_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr13 =",ierr
           ! eq c8 + c12  :: -\MATone af \ts fi
           ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*T(f,i)",gamma%prime_vo,MATone_vv_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr14 =",ierr
           ! eq c13 + c14 :: - \ld mnef \rs em \tauNo fain
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(n,f)+=L(m,n,e,f)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr15 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(n,f)*TauNo(f,a,i,n)",gamma%prime_vo,ov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr16 =",ierr
           ! eq c5        :: +\ts ai
           ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr17 =",ierr
           ! eq c6        :: +\ls me \rd aeim
           ierr=talsh_tensor_contract("G(a,i)+=L(m,e)*R(a,e,i,m)",gamma%prime_vo,l1_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr18 =",ierr
           ! eq c9        :: -1/2 \ld mnef \rs fi \td eamn
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr19 =",ierr
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*R(f,i)",gamma%prime_vo,vv_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr20 =",ierr
           ! eq c10       :: - 1/2 \ld mnef \rs n \td efmi
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("X(n,i)+=L(m,n,e,f)*T(e,f,m,i)",oo_aux_tensor,l2_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr21 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(n,i)*R(a,n)",gamma%prime_vo,oo_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr22 =",ierr


    ! rho_ia
           ! eq d1        :: +\ro \ls ia
           ierr=talsh_tensor_add("G(i,a)+=L(i,a)",gamma%prime_ov,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr23 =",ierr
           ! eq d2        :: +\ld imae \rs em
           ierr=talsh_tensor_contract("G(i,a)+=L(i,m,a,e)*R(e,m)",gamma%prime_ov,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr24 =",ierr


    !  Symmetrization


           ierr=talsh_tensor_add("G(i,a)+=X(i,a)",gamma%symm_ov,gamma%prime_ov,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr25 =",ierr
           ierr=talsh_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr26 =",ierr
           ierr=talsh_tensor_add("G(a,i)+=X(a,i)",gamma%symm_vo,gamma%prime_vo,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr27 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%prime_ov,one_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr28 =",ierr


    !     debug
          if (.true.) then
            print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)
            if (.not.CCD) then
              print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
              print*, "gsymm_ov   = ", tensor_norm2(gamma%symm_ov)
              print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
              print*, "gsymm_vo   = ", tensor_norm2(gamma%symm_vo)
            end if
            print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
            print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
            print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)
          end if

    !**************
    !     cleanup *
    !**************

          ierr=talsh_tensor_destruct(one_tensor)
          if (.not.CCD) then
            ierr=talsh_tensor_destruct(tau2_tensor)
            ierr=talsh_tensor_destruct(vv_aux_tensor)
            ierr=talsh_tensor_destruct(oo_aux_tensor)
          end if

    end subroutine talsh_calc_dm_1rdm_ccsd_EE_mo



        subroutine talsh_calc_dm_2rdm_ccsd_EE_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)


          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor
          type(talsh_tens_t) :: ovvo_aux_tensor, ovov_aux_tensor,voov_aux_tensor
          type(talsh_tens_t) :: oovo_aux_tensor

    !     tensor dimensions
          integer :: one_dims(1)=1

      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

    !*********************************
    !     construct auxiliary tensor *
    !*********************************

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)

          ierr=talsh_tensor_construct(gamma%prime_oooo,C8,oooo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_ovoo,C8,ovoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_voov,C8,voov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvvo,C8,vvvo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvvv,C8,vvvv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_oovo,C8,oovo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vovv,C8,vovv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_oovv,C8,oovv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvoo,C8,vvoo_dims,init_val=ZERO)

          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
          ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vo_aux_tensor,C8,vo_dims,init_val=ZERO)

    !********************
    ! Intermediates     *
    !********************
            ierr=talsh_tensor_construct(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(vvoo_aux_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)

            ! Loic_Taille

            ierr=talsh_tensor_construct(MATone_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ovvv_aux_tensor,C8,ovvv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ovov_aux_tensor,C8,ovov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(voov_aux_tensor,C8,voov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauUP_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauDown_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauNo_tensor,C8,vvoo_dims,init_val=ZERO)


            !  Loic_Tenseurs_Intermediaires


           !  MATone_oo
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,j)*L(i,e)",MATone_oo_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,f,j,m)*L(i,m,e,f)",MATone_oo_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,j)*L(i,e)",MATtwo_oo_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,f,j,m)*L(i,m,e,f)",MATtwo_oo_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(a,m)*L(m,b)",MATone_vv_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(e,a,m,n)*L(m,n,e,b)",MATone_vv_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(a,m)*L(m,b)",MATtwo_vv_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(e,a,m,n)*L(m,n,e,b)",MATtwo_vv_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           !  Tau
           
           ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauUp(a,b,i,j)+=T(a,b,i,j)",TauUp_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(a,i)*T(b,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(b,i)*T(a,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

    !********************
    !     calculate DMs *
    !********************


    ! rho ijkl
           ! eq 1              :: +1/2 \ld ijef \rd efkl 
    ierr=talsh_tensor_contract("G(i,j,k,l)+=L(i,j,e,f)*R(e,f,k,l)",gamma%prime_oooo,l2_tensor,r2_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr1 =",ierr
    ! eq 2              :: +\Pm kl \ld ijef \rs ek \ts fl
    ierr=talsh_tensor_contract("A(e,f,k,l)+=R(e,k)*R(f,l)",vvoo_aux_tensor,r1_tensor,r1_tensor)
    ierr=talsh_tensor_contract("A(e,f,k,l)+=R(e,l)*R(f,k)",vvoo_aux_tensor,r1_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr2 =",ierr
    ierr=talsh_tensor_contract("G(i,j,k,l)+=A(e,f,k,l)*L(i,j,e,f)",gamma%prime_oooo,vvoo_aux_tensor,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr3 =",ierr
    ! eq 3              :: +\ro 1/2 \ld ijef \TauDown efkl 
    ierr=talsh_tensor_contract("G(i,j,k,l)+=L(i,j,e,f)*TauDown(e,f,k,l)",gamma%prime_oooo,l2_tensor,TauDown_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr4 =",ierr
      
        !!! Pas de Kronecker

    ! rho iajk
    ! eq j2 + j10        :: - (1/2) \Pm jk \MATone ij \ts ak
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,j)*T(a,k)",gamma%prime_ovoo,MATone_oo_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,k)*T(a,j)",gamma%prime_ovoo,MATone_oo_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr5 =",ierr

    ! eq j3 + j9         :: - (1/2) \Pm jk \rs ak \MATtwo ij
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,j)*R(a,k)",gamma%prime_ovoo,MATtwo_oo_tensor,r1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,k)*R(a,j)",gamma%prime_ovoo,MATtwo_oo_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr6 =",ierr

    ! eq j6 + j12        :: + (1/2) \Pm jk \ld imfe \rs ek \TauNo fajm
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(i,m,k,f)+=L(i,m,e,f)*R(e,k)",ooov_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr7 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,m,k,f)*TauNo(f,a,j,m)",gamma%prime_ovoo,&
         ooov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,m,j,f)*TauNo(f,a,k,m)",gamma%prime_ovoo,ooov_aux_tensor,TauNo_tensor)
    if(ierr.ne.0) write(*,*) "ierr8 =",ierr

    ! eq j7 + j11        :: + (1/4) \ld imfe \rs am \TauDown efjk
    ierr=talsh_tensor_init(ovvv_aux_tensor)
    ierr=talsh_tensor_contract("A(i,a,f,e)+=L(i,m,f,e)*R(a,m)",ovvv_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr9 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,a,f,e)*TauDown(e,f,j,k)",gamma%prime_ovoo,&
         ovvv_aux_tensor,TauDown_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr10 =",ierr

    ! eq j8 + j13        :: - (1/2) \ld imfe \rs em \TauDown fajk
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(i,f)+=L(i,m,f,e)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr11 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,f)*TauDown(f,a,j,k)",gamma%prime_ovoo,ov_aux_tensor,TauDown_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr12 =",ierr

    ! eq j1              :: - (1/2) \ls ie \rd eajk
    ierr=talsh_tensor_contract("G(i,a,j,k)+=L(i,e)*R(e,a,j,k)",gamma%prime_ovoo,l1_tensor,r2_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr13 =",ierr

    ! eq j4              :: + (1/2) \Pm jk \ld imfe \rd fajm \ts ek
    ierr=talsh_tensor_init(ovvo_aux_tensor)
    ierr=talsh_tensor_contract("A(i,a,e,j)+=L(i,m,e,f)*R(f,a,m,j)",ovvo_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr14 =",ierr
    ierr= talsh_tensor_contract("G(i,a,j,k)+=A(i,a,e,j)*T(e,k)",gamma%prime_ovoo,ovvo_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,a,e,k)*T(e,j)",gamma%prime_ovoo,ovvo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr15 =",ierr

    ! eq j5              :: + (1/4) \ld imfe \rd fejk \ts am
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_contract("A(i,m,j,k)+=L(i,m,f,e)*R(f,e,j,k)",oooo_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr16 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,m,j,k)*T(a,m)",gamma%prime_ovoo,oooo_aux_tensor,t1_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr17 =",ierr

    ! eq j15 + j18       :: - \ro (1/2) \Pm jk \MATtwo ij \ts ak
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,j)*T(a,k)",gamma%prime_ovoo,MATtwo_oo_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr18 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=M(i,k)*T(a,j)",gamma%prime_ovoo,MATtwo_oo_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr19 =",ierr

    ! eq j17+j19         :: + \ro (1/4) \ld imfe \ts am \TauUp fejk
    ierr=talsh_tensor_init(ovvv_aux_tensor)
    ierr=talsh_tensor_contract("A(i,a,f,e)+=L(i,m,f,e)*T(a,m)",ovvv_aux_tensor,l2_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr20 =",ierr
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,a,f,e)*TauDown(e,f,j,k)",gamma%prime_ovoo,&
         ovvv_aux_tensor,TauDown_tensor,scale=ONE_HALF) 
    if(ierr.ne.0) write(*,*) "ierr21 =",ierr

    ! eq j14             :: - \ro (1/2) \ls ie \td eajk
    ierr=talsh_tensor_contract("G(i,a,j,k)+=L(i,e)*T(e,a,j,k)",gamma%prime_ovoo,l1_tensor,t2_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr22 =",ierr

    ! eq j16             :: + (1/2) \Pm jk \ld imfe \td fajm \ts ek
    ierr=talsh_tensor_init(ovvo_aux_tensor)
    ierr=talsh_tensor_contract("A(i,a,e,j)+=L(i,m,e,f)*T(f,a,m,j)",ovvo_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr23 =",ierr
    ierr= talsh_tensor_contract("G(i,a,j,k)+=A(i,a,e,j)*T(e,k)",gamma%prime_ovoo,ovvo_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,a,e,k)*T(e,j)",gamma%prime_ovoo,ovvo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr24 =",ierr

  ! rho abci

    ! eq g2+g10          :: + (1/2) \Pm ab \MATone ac \ts bi
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(a,c)*T(b,i)",gamma%prime_vvvo,MATone_vv_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(b,c)*T(a,i)",gamma%prime_vvvo,MATone_vv_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr26 =",ierr

    ! eq g3+g9           :: + (1/2) \Pm ab \MATtwo ac \rs bi
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(a,c)*R(b,i)",gamma%prime_vvvo,MATtwo_vv_tensor,r1_tensor)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(b,c)*R(a,i)",gamma%prime_vvvo,MATtwo_vv_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr27 =",ierr

    ! eq g7+g11          :: - (1/4) \ld mnec \rs ei \TauUp abmn        
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,c,i)+=L(m,n,c,e)*R(e,i)",oovo_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr28 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,n,c,i)*TauUp(a,b,m,n)",gamma%prime_vvvo,oovo_aux_tensor,TauUp_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr29 =",ierr

    ! eq g8+g13          :: + (1/2) \ld mnec \rs em \TauUp abni  
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(n,c)+=L(m,n,e,c)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr30 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(n,c)*TauUp(a,b,n,i)",gamma%prime_vvvo,ov_aux_tensor,TauUp_tensor)
    if(ierr.ne.0) write(*,*) "ierr31 =",ierr

    ! eq g6+g12          :: -(1/2) \Pm ab \ld mnec \rs bm \TauNo eain
    ierr=talsh_tensor_init(ovvv_aux_tensor)
    ierr=talsh_tensor_contract("A(n,b,e,c)+=L(n,m,e,c)*R(b,m)",ovvv_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr32 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(n,b,e,c)*TauNo(e,a,i,n)",gamma%prime_vvvo,ovvv_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(n,a,e,c)*TauNo(e,b,i,n)",gamma%prime_vvvo,&
         ovvv_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr33 =",ierr

    ! eq g1              :: +(1/2) \ls mc \rd ba im
    ierr=talsh_tensor_contract("G(a,b,c,i)+=L(m,c)*R(b,a,i,m)",gamma%prime_vvvo,l1_tensor,r2_tensor) 
    if(ierr.ne.0) write(*,*) "ierr34 =",ierr

    ! eq g4              :: -(1/2) \Pm ab \ld mnec \rd eain \ts bm
    ierr=talsh_tensor_init(ovvo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,c,a,i)+=L(m,n,c,e)*R(e,a,n,i)",ovvo_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr35 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,a,c,i)*T(b,m)",gamma%prime_vvvo,&
        ovvo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,b,c,i)*T(a,m)",gamma%prime_vvvo,ovvo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr36 =",ierr

    ! eq g5              :: -(1/4) \ld mnec \rd bamn \ts ei
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,c)+=L(m,n,e,c)*T(e,i)",ooov_aux_tensor,l2_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr37 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,n,i,c)*R(b,a,m,n)",gamma%prime_vvvo,&
        ooov_aux_tensor,r2_tensor,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr38 =",ierr

    ! eq g15+g18         :: + \ro (1/2) \Pm ab \ts bi \MATtwo ac
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(a,c)*T(b,i)",gamma%prime_vvvo,MATtwo_vv_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(b,c)*T(a,i)",gamma%prime_vvvo,&
         MATtwo_vv_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr39 =",ierr

    ! eq g17+g19         :: - \ro (1/4) \ld mnec \ts ei \TauUp bamn                      XX
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,c)+=L(m,n,e,c)*T(e,i)",ooov_aux_tensor,l2_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr40 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,n,i,c)*TauUp(b,a,m,n)",gamma%prime_vvvo,&
         ooov_aux_tensor,TauUp_tensor,scale=MINUS_ONE_HALF) 
    if(ierr.ne.0) write(*,*) "ierr41 =",ierr

    ! eq g14             :: + \ro (1/2) \ls mc \td baim    XX
    ierr=talsh_tensor_contract("G(a,b,c,i)+=L(m,c)*T(a,b,m,i)",gamma%prime_vvvo,l1_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr42 =",ierr

    ! eq g16             :: - \ro \Pm ab \ld mnec \ts bm \td eain
    ierr=talsh_tensor_init(ovov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,a,i,c)+=L(m,n,e,c)*T(e,a,i,n)",ovov_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr43 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,a,i,c)*T(b,m)",gamma%prime_vvvo,ovov_aux_tensor,t1_tensor,scale=MINUS_ONE) 
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,b,i,c)*T(a,m)",gamma%prime_vvvo,ovov_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr44 =",ierr

    ! rho ijak
    ! eq                 :: -(1/2) \ld ijea \rs ek
    ierr=talsh_tensor_contract("G(i,j,a,k)+=L(i,j,e,a)*R(e,k)",gamma%prime_oovo,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr45 =",ierr
    ! eq                 :: - \ro (1/2) \ld ijea \ts ek
    ierr=talsh_tensor_contract("G(i,j,a,k)+=L(i,j,e,a)*T(e,k)",gamma%prime_oovo,l2_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr46 =",ierr


    ! rho aibj XX
    ! eq  f1             :: - 1/4 \Pm ij \Pm ab \ls ib \rs aj
    ierr=talsh_tensor_contract("G(a,i,j,b)+=R(a,j)*L(i,b)",gamma%prime_voov,r1_tensor,l1_tensor,scale=MINUS_ONE_QUARTER)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=R(b,j)*L(i,a)",gamma%prime_voov,r1_tensor,l1_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=R(a,i)*L(j,b)",gamma%prime_voov,r1_tensor,l1_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=R(b,i)*L(j,a)",gamma%prime_voov,r1_tensor,l1_tensor,scale=MINUS_ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr47 =",ierr
    ! eq f2              :: - 1/4 \Pm ij \pm ab \ld imeb \rd eajm
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,b)*R(e,a,j,m)",gamma%prime_voov,l2_tensor,r2_tensor, scale= MINUS_ONE_QUARTER)    
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,b)*R(e,a,i,m)",gamma%prime_voov,l2_tensor,r2_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,a)*R(e,b,j,m)",gamma%prime_voov,l2_tensor,r2_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,a)*R(e,b,i,m)",gamma%prime_voov,l2_tensor,r2_tensor, scale= MINUS_ONE_QUARTER)    
    if(ierr.ne.0) write(*,*) "ierr48 =",ierr
    ! eq f3              :: -1/4 \Pm ij \Pm ab \ld iimeb \rs am \ts ej
    ierr=talsh_tensor_init(ovvv_aux_tensor) 
    ierr=talsh_tensor_contract("A(i,a,e,b)+=L(i,m,e,b)*R(a,m)",ovvv_aux_tensor,l2_tensor,r1_tensor)
    ierr=talsh_tensor_contract("A(i,a,e,b)+=L(i,m,e,a)*R(b,m)",ovvv_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr49 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,a,e,b)*T(e,j)",gamma%prime_voov,ovvv_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,a,e,b)*T(e,i)",gamma%prime_voov,ovvv_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr50 =",ierr
    ! eq f4              :: -1/4 \Pm ij \Pm ab \ld imeb \rs ej \ts am
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(i,m,j,b)+=L(i,m,e,b)*R(e,j)",ooov_aux_tensor,l2_tensor,r1_tensor)
    ierr=talsh_tensor_contract("A(i,m,j,b)+=L(j,m,e,b)*R(e,i)",ooov_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr51 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,m,j,b)*T(a,m)",gamma%prime_voov,ooov_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,m,j,a)*T(b,m)",gamma%prime_voov,ooov_aux_tensor,t1_tensor,scale=ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr52 =",ierr
    ! eq f5              :: - \Pm ij \Pm ab \ld imbe \rs em \ts aj
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(i,b)+=L(i,m,e,b)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr53 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,b)*,T(a,j)",gamma%prime_voov,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,a)*,T(b,j)",gamma%prime_voov,ov_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,b)*,T(a,i)",gamma%prime_voov,ov_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,a)*,T(b,i)",gamma%prime_voov,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr54 =",ierr


    ! eq  f6             :: \ro \Pm ij \Pm ab \ls ib \ts aj
    ierr=talsh_tensor_contract("G(a,i,j,b)+=T(a,j)*L(i,b)",gamma%prime_voov,t1_tensor,l1_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=T(b,j)*L(i,a)",gamma%prime_voov,t1_tensor,l1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=T(a,i)*L(j,b)",gamma%prime_voov,t1_tensor,l1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=T(b,i)*L(j,a)",gamma%prime_voov,t1_tensor,l1_tensor)
    if(ierr.ne.0) write(*,*) "ierr55 =",ierr
    ! eq f7              :: - \ro 1/4 \Pm ij \pm ab \ld imeb \td eajm
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,b)*T(e,a,j,m)",gamma%prime_voov,l2_tensor,t2_tensor,scale= MINUS_ONE_QUARTER)    
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,b)*T(e,a,i,m)",gamma%prime_voov,l2_tensor,t2_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,a)*T(e,b,j,m)",gamma%prime_voov,l2_tensor,t2_tensor)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,a)*T(e,b,i,m)",gamma%prime_voov,l2_tensor,t2_tensor,scale= MINUS_ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr56 =",ierr
    ! eq f8              :: - \ro 1/4 \Pm ij \Pm ab \ld imeb \ts am \ts ej
    ierr=talsh_tensor_init(ovvv_aux_tensor)
    ierr=talsh_tensor_contract("A(i,a,e,b)+=L(i,m,e,b)*T(a,m)",ovvv_aux_tensor,l2_tensor,t1_tensor)
    ierr=talsh_tensor_contract("A(i,a,e,b)+=L(i,m,e,a)*T(b,m)",ovvv_aux_tensor,l2_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr57 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,a,e,b)*T(e,j)",gamma%prime_voov,ovvv_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,a,e,b)*T(e,i)",gamma%prime_voov,ovvv_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr58 =",ierr

    ! rho aibc
    ! eq                 :: + (1/2) \ld mibc \rs am
    ierr=talsh_tensor_contract("G(a,i,b,c)+=R(a,m)*L(m,i,b,c)",gamma%prime_vovv,r1_tensor,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr59 =",ierr
    ! eq                 :: + \ro (1/2) \ld mibc \ts am
    ierr=talsh_tensor_contract("G(a,i,b,c)+=T(a,m)*L(m,i,b,c)",gamma%prime_vovv,t1_tensor,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr60 =",ierr

    ! rho abcd     XX
    ! eq i1              :: +(1/4) 1/2 \ld mncd \rd abmn 
    ierr=talsh_tensor_contract("G(a,b,c,d)+=L(m,n,c,d)*R(a,b,m,n)",gamma%prime_vvvv,l2_tensor,r2_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr61 =",ierr
    ! eq i2              :: +(1/4) \Pm ab \ld mncd \rs am \ts bn
    ierr=talsh_tensor_init(ovvv_aux_tensor)
    ierr=talsh_tensor_contract("A(n,a,c,d)+=R(a,m)*L(m,n,c,d)",ovvv_aux_tensor,r1_tensor,l2_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr62 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,d)+=A(n,a,c,d)*T(b,n)",gamma%prime_vvvv,ovvv_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,c,d)+=A(n,b,c,d)*T(a,n)",gamma%prime_vvvv,ovvv_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr63 =",ierr
    ! eq i3+i4           :: \ro (1/4) 1/2 \ld mncd \TauUp abmn
    ierr=talsh_tensor_contract("G(a,b,c,d)+=L(m,n,c,d)*TauUp(a,b,m,n)",gamma%prime_vvvv,l2_tensor,TauUp_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr64 =",ierr


    ! rho ijab    XX
    ! eq 1               :: +\ro \ld ijab 
    ierr=talsh_tensor_add("G(i,j,a,b)+=L(i,j,a,b)",gamma%prime_oovv,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr65 =",ierr

    ! rho abij
    ! eq h4+h5+h6+h7     :: - (1/4) \Pm ij  \MATone ni \TauUp abnj
    ierr=talsh_tensor_contract("G(a,b,i,j)+=M(n,i)*TauUp(a,b,n,j)",gamma%prime_vvoo,MATone_oo_tensor,TauUp_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=M(n,j)*TauUp(a,b,n,i)",gamma%prime_vvoo,MATone_oo_tensor,TauUp_tensor)
    if(ierr.ne.0) write(*,*) "ierr66 =",ierr
    ! eq h6+h7+h10+h11   :: - (1/4) \Pm ab \MATone ae \TauDown ebij
    ierr=talsh_tensor_contract("G(a,b,i,j)+=M(a,e)*TauDown(e,b,i,j)",gamma%prime_vvoo,&
         MATone_vv_tensor,TauDown_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=M(b,e)*TauDown(e,a,i,j)",gamma%prime_vvoo,MATone_vv_tensor,TauDown_tensor)
    if(ierr.ne.0) write(*,*) "ierr67 =",ierr
    ! eq h12+h14         :: + (1/4) 1/4 \ld mnef \rd ef ij \TauUp abmn
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,j)+=L(m,n,e,f)*R(e,f,i,j)",oooo_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr68 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,n,i,j)*TauUp(a,b,m,n)",gamma%prime_vvoo,&
        oooo_aux_tensor,TauUp_tensor,scale=ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr69 =",ierr
    ! eq j13+j15         :: + (1/4) 1/4 \ld mnef \rd abmn \TauDown efij
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,j)+=L(m,n,e,f)*TauDown(e,f,i,j)",oooo_aux_tensor,l2_tensor,TauDown_tensor)
    if(ierr.ne.0) write(*,*) "ierr70 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,n,i,j)*R(a,b,m,n)",gamma%prime_vvoo,&
        oooo_aux_tensor,r2_tensor,scale=ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr71 =",ierr
    ! eq j24+j26         :: + (1/4) 1/2 \Pm ij \ld mnef \rs ei \ts fj \TauUp abmn
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_init(ooov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,f)+=L(m,n,e,f)*R(e,i)",ooov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr72 =",ierr
    ierr=talsh_tensor_contract("B(m,n,i,j)+=A(m,n,i,f)*T(f,j)",oooo_aux_tensor,ooov_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("B(m,n,i,j)+=A(m,n,j,f)*T(f,i)",oooo_aux_tensor,ooov_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr73 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(m,n,i,j)*TauUp(a,b,m,n)",gamma%prime_vvoo,&
        oooo_aux_tensor,TauUp_tensor,scale=ONE_HALF) 
    if(ierr.ne.0) write(*,*) "ierr74 =",ierr
    ! eq j25+j27         :: + (1/4) 1/2 \Pm ab \ld mnef \rs am \ts bn \TauDown efij
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_init(vvoo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,j)+=L(m,n,e,f)*TauDown(e,f,i,j)",oooo_aux_tensor,l2_tensor,TauDown_tensor)
    if(ierr.ne.0) write(*,*) "ierr75 =",ierr
    ierr=talsh_tensor_contract("B(a,b,m,n)+=R(a,m)*T(b,n)",vvoo_aux_tensor,r1_tensor,t1_tensor)
    ierr=talsh_tensor_contract("B(a,b,m,n)+=R(b,m)*T(a,n)",vvoo_aux_tensor,r1_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr76 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,n,i,j)*B(a,b,m,n)",gamma%prime_vvoo,&
         oooo_aux_tensor,vvoo_aux_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr77 =",ierr
    ! eq j35+j36         :: - (1/4) \Pm ij \ld mnef \rs em \ts fi \TauUp abnj
    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(n,f)+=L(m,n,e,f)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr78 =",ierr
    ierr=talsh_tensor_contract("B(n,i)+=A(n,f)*T(f,i)",oo_aux_tensor,ov_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr79 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,i)*TauUp(a,b,n,j)",gamma%prime_vvoo,&
        oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,j)*TauUp(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor)
    if(ierr.ne.0) write(*,*) "ierr80 =",ierr

    ! eq j18+19          :: - (1/4) \Pm ij \Pm ab \ls me \rs bj \TauNo aemi
    ierr=talsh_tensor_init(ovvo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,b,e,j)+=L(m,e)*R(b,j)",ovvo_aux_tensor,l1_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr81 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,b,e,j)*TauNo(a,e,m,i)",gamma%prime_vvoo,&
         ovvo_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,a,e,j)*TauNo(b,e,m,i)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,b,e,i)*TauNo(a,e,m,j)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,a,e,i)*TauNo(b,e,m,j)",gamma%prime_vvoo,&
        ovvo_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr82 =",ierr
    ! eq j20+j21         :: - (1/4) \Pm ij \Pm ab \ld mnef \rd aeim \TauNo bfnj
    ierr=talsh_tensor_init(voov_aux_tensor)
    ierr=talsh_tensor_contract("A(a,n,i,f)+=R(a,e,i,m)*L(m,n,e,f)",voov_aux_tensor,r2_tensor,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr83 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,i,f)*TauNo(b,f,n,j)",gamma%prime_vvoo,&
         voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,i,f)*TauNo(a,f,n,j)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,j,f)*TauNo(b,f,n,i)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,j,f)*TauNo(a,f,n,i)",gamma%prime_vvoo,&
         voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr84 =",ierr
    ! eq j16+j30         :: - (1/4) \Pm ij \rd abin \MATtwo nj
    ierr=talsh_tensor_contract("G(a,b,i,j)+=R(a,b,i,n)*M(n,j)",gamma%prime_vvoo,r2_tensor,MATtwo_oo_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=R(a,b,j,n)*M(n,i)",gamma%prime_vvoo,r2_tensor,MATtwo_oo_tensor)
    if(ierr.ne.0) write(*,*) "ierr85 =",ierr
    ! eq j17+j31         :: - (1/4) \Pm ab \rd afij \MATtwo bn
    ierr=talsh_tensor_contract("G(a,b,i,j)+=R(a,f,i,j)*M(b,f)",gamma%prime_vvoo,r2_tensor,MATtwo_vv_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=R(b,f,i,j)*M(a,f)",gamma%prime_vvoo,r2_tensor,MATtwo_vv_tensor)
    if(ierr.ne.0) write(*,*) "ierr86 =",ierr

    ! eq j1              :: + (1/4) \td abij
    ierr=talsh_tensor_add("G(a,b,i,j)+=T(a,b,i,j)",gamma%prime_vvoo,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr87 =",ierr
    ! eq j2              :: + (1/4) 1/2 \Pm ab \Pm ij \ts ai \ts bj
    ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,j)*T(a,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr88 =",ierr
    ! eq j3              :: + (1/4) \ls me \rd aeim \ts bj
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(a,i)+=R(a,e,i,m)*L(m,e)",vo_aux_tensor,r2_tensor,l1_tensor)
    if(ierr.ne.0) write(*,*) "ierr89 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr90 =",ierr
    ! eq j22             :: - (1/4) 1/2 \Pm ij \Pm ab \ld mnef \rs fi \td eamn \ts bj
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr91 =",ierr
    ierr=talsh_tensor_contract("B(a,i)+=A(a,f)*R(f,i)",vo_aux_tensor,vv_aux_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr90 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr92 =",ierr
! Signe XX ! eq j23             :: - (1/4) 1/2 \Pm ij \Pm ab \ld mnef \rs an \td efmi \ts bj
    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(n,i)+=L(m,n,e,f)*T(e,f,m,i)",oo_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr93 =",ierr
    ierr=talsh_tensor_contract("B(a,i)+=X(n,i)*R(a,n)",vo_aux_tensor,oo_aux_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr94 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr95 =",ierr

    ! eq j28             :: - (1/4)  \Pm ij \Pm ab \ld mnef \rs an \td bejm \ts fi
    ierr=talsh_tensor_init(voov_aux_tensor)
    ierr=talsh_tensor_init(vvoo_aux_tensor)
    ierr=talsh_tensor_contract("A(b,n,j,f)+=L(m,n,e,f)*T(b,e,j,m)",voov_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr96 =",ierr
    ierr=talsh_tensor_contract("B(a,f,n,i)+=R(a,n)*T(f,i)",vvoo_aux_tensor,r1_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr97 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,j,f)*B(a,f,n,i)",gamma%prime_vvoo,&
         voov_aux_tensor,vvoo_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,j,f)*B(b,f,n,i)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,i,f)*B(a,f,n,j)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,i,f)*B(b,f,n,j)",gamma%prime_vvoo,&
         voov_aux_tensor,vvoo_aux_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr98 =",ierr
!! I keep  voov_aux - Not usefull :  ierr=talsh_tensor_contract("X(n,b,f,j)+=L(m,n,e,f)*T(b,e,j,m)",ovvo_aux_tensor,l2_tensor,t2_tensor)
    ! eq j29              :: - (1/4) \Pm ij \Pm ab \ld mnef \rs fi \tdbejm \ts an
    ierr=talsh_tensor_init(vvoo_aux_tensor) 
    ierr=talsh_tensor_contract("C(a,f,n,i)+=T(a,n)*R(f,i)",vvoo_aux_tensor,t1_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr99 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=C(b,n,j,f)*B(a,f,n,i)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=C(a,n,j,f)*B(b,f,n,i)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=C(b,n,i,f)*B(a,f,n,j)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=C(a,n,i,f)*B(b,f,n,j)",gamma%prime_vvoo,voov_aux_tensor,vvoo_aux_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr100 =",ierr

 ! XX          ! eq j32             :: - (1/4) 1/2 \Pm ij \Pm ab \ld mnef \rs bj \td eamn \ts fi
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr101 =",ierr
    ierr=talsh_tensor_contract("B(a,i)+=A(a,f)*T(f,i)",vo_aux_tensor,vv_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr102 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*R(b,j)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*R(b,i)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*R(a,j)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*R(a,i)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr103 =",ierr

    ! eq j33             :: + (1/4) \Pm ij \Pm ab \ld mnef \rs em \td afin \ts bj 
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(n,f)+=L(m,n,e,f)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr104 =",ierr
    ierr=talsh_tensor_contract("B(a,i)+=A(n,f)*T(a,f,i,n)",vo_aux_tensor,ov_aux_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr105 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr106 =",ierr

    ! eq 34              :: - (1/4) 1/2 \Pm ij \Pm ab \ld mnef \rs bj \td efim \ts an
    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(n,i)+=L(m,n,e,f)*T(e,f,i,m)",oo_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr107 =",ierr
    ierr=talsh_tensor_contract("B(a,i)+=A(n,i)*T(a,n)",vo_aux_tensor,oo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr108 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*R(b,j)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*R(b,i)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*R(a,j)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*R(a,i)",gamma%prime_vvoo,vo_aux_tensor,r1_tensor,scale=MINUS_ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr109 =",ierr

    ! eq 37              :: - (1/4) \Pm ab \ld mnef \rs em \td fbij \ts an
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(n,f)+=L(m,n,e,f)*R(e,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr110 =",ierr
    ierr=talsh_tensor_contract("B(a,f)+=A(n,f)*T(a,n)",vv_aux_tensor,ov_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr111 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr112 =",ierr

    !!!! rho abij la partie \ro
    ! eq j3              :: + (1/4) \ro  \Pm ij \Pm ab \ls me \td aeim \ts bj
    ierr=talsh_tensor_init(vo_aux_tensor)
    ierr=talsh_tensor_contract("A(a,i)+=L(m,e)*T(a,e,i,m)",vo_aux_tensor,l1_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr113 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr114 =",ierr

    ! eq h4+h5           :: - (1/4) \ro \Pm ij \ls me \ts ei \TauUp abmj
    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,i)+=L(m,e)*T(e,i)",oo_aux_tensor,l1_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr115 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,i)*TauUp(a,b,m,j)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,j)*TauUp(a,b,m,i)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr116 =",ierr

!XX Erreur sur \Pm ij           ! eq h6+h7           :: - (1/4) \ro 1/2 \Pm ij \Pm ab \ld mnef \td ea mn \TauNo fbij
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("A(a,f)+=L(m,n,e,f)*T(e,a,m,n)",vv_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr117 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,f)*X(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,TauNo_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,f)*X(f,a,i,j)",gamma%prime_vvoo,vv_aux_tensor,TauNo_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr118 =",ierr

!XX Erreur sur \Pm ij           ! eq h8+h9           :: - (1/4) \ro 1/2 \Pm ij \Pm ab \ld mnef  \td feim \TauNo abnj
    ierr=talsh_tensor_init(oo_aux_tensor)
    ierr=talsh_tensor_contract("A(n,i)+=L(m,n,e,f)*T(f,e,i,m)",oo_aux_tensor,l2_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr119 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,i)*X(a,b,n,j)",gamma%prime_vvoo,oo_aux_tensor,TauNo_tensor,scale=MINUS_ONE_HALF)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,j)*X(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,TauNo_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr120 =",ierr

    ! eq j11             :: - (1/4) \ro \Pm ab \s me \ts am \td ebij
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("A(a,e)+=L(m,e)*T(a,m)",vv_aux_tensor,l1_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr121 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,e)*T(e,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,e)*T(e,a,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
    if(ierr.ne.0) write(*,*) "ierr122 =",ierr

    ! eq h20+h21         :: - (1/4) \ro \Pm ij \Pm ab \ld mnef \td aeim \TauNo fbjn
    ierr=talsh_tensor_init(voov_aux_tensor)
    ierr=talsh_tensor_contract("A(a,n,i,f)+=T(a,e,i,m)*L(m,n,e,f)",voov_aux_tensor,t2_tensor,l2_tensor)
    if(ierr.ne.0) write(*,*) "ierr123 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,i,f)*X(f,b,j,n)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,i,f)*X(f,a,j,n)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,j,f)*X(f,b,i,n)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,j,f)*X(f,a,i,n)",gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr124 =",ierr

    ! eq h12+h14+h15+h24 :: + (1/4) \ro 1/4 \ld mnef \TauUp ab mn \TauUp efij
    ierr=talsh_tensor_init(oooo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,i,j)+=L(m,n,e,f)*TauUp(e,f,i,j)",oooo_aux_tensor,l2_tensor,TauUp_tensor)
    if(ierr.ne.0) write(*,*) "ierr125 =",ierr
    ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,n,i,j)*TauUp(a,b,m,n)",gamma%prime_vvoo,&
         oooo_aux_tensor,TauUp_tensor,scale=ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr126 =",ierr


    !     debug
    print *, "gamma%prime_ovoo  = ", tensor_norm2(gamma%prime_ovoo)
    print *, "gamma%prime_oooo  = ", tensor_norm2(gamma%prime_oooo) 
    print *, "gamma%prime_vvvo  = ", tensor_norm2(gamma%prime_vvvo) 
    print *, "gamma%prime_oovo  = ", tensor_norm2(gamma%prime_oovo)  
    print *, "gamma%prime_voov  = ", tensor_norm2(gamma%prime_voov)  
    print *, "gamma%prime_vovv  = ", tensor_norm2(gamma%prime_vovv)  
    print *, "gamma%prime_vvvv  = ", tensor_norm2(gamma%prime_vvvv)  
    print *, "gamma%prime_oovv  = ", tensor_norm2(gamma%prime_oovv)  
    print *, "gamma%prime_vvoo  = ", tensor_norm2(gamma%prime_vvoo)  


    !**************
    !     cleanup *
    !**************

          ierr=talsh_tensor_destruct(one_tensor)
          if (.not.CCD) then
            ierr=talsh_tensor_destruct(tau2_tensor)
            ierr=talsh_tensor_destruct(vv_aux_tensor)
            ierr=talsh_tensor_destruct(oo_aux_tensor)
          end if

    end subroutine talsh_calc_dm_2rdm_ccsd_EE_mo



        subroutine talsh_calc_dm_1rdm_ccsd_IP_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)


          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor
          type(talsh_tens_t) :: oovo_aux_tensor

    !     tensor dimensions
          integer :: one_dims(1)=1
      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

    !*********************************
    !     construct auxiliary tensor *
    !*********************************

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)


          ierr=talsh_tensor_construct(gamma%prime_oo,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_ov,C8,ov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vo,C8,vo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vv,C8,vv_dims,init_val=ZERO)


          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

          ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vvoo_aux_tensor,C8,vvoo_dims,init_val=ZERO)

          ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
 
    ! rho ij
           ! eq a3              :: - \ls iX \rs Xj   
           ierr=talsh_tensor_contract("G(i,j)+=L(i,X)*R(X,j)",gamma%prime_oo,l1_tensor,r1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr1 =",ierr
           ! eq a4              :: - \ld mieX \rd eXmj   
           ierr=talsh_tensor_contract("G(i,j)+=L(m,i,e,X)*R(e,X,m,j)",gamma%prime_oo,l2_tensor,r2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr2 =",ierr
           ! eq a5      :: - \ld imeX \rs Xm \ts ej
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("A(i,e)+=L(i,m,e,X)*R(X,m)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr3 =",ierr
           ierr=talsh_tensor_contract("G(i,j)+=A(i,e)*T(e,j)",gamma%prime_oo,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr4 =",ierr

    ! rho ab
           ! eq         :: + 1/2 \ld mnbX \rd aXmn 
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("G(a,b)+=L(m,n,b,X)*R(a,X,m,n)",gamma%prime_vv,l2_tensor,r2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr5 =",ierr
           ! eq         :: + \ld mnbX \rs Xn \ts am
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("A(m,b)+=L(m,n,b,X)*R(X,n)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr6 =",ierr
           ierr=talsh_tensor_contract("G(a,b)+=A(m,b)*T(a,m)",gamma%prime_vv,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr7 =",ierr


    ! rho ai
           ! eq c6        :: +\ls me \rd aeim
           ierr=talsh_tensor_contract("G(a,i)+=L(m,X)*R(a,X,i,m)",gamma%prime_vo,l1_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr8 =",ierr
           ! eq c7        :: - \ls nx \rs xi \ts am
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("A(n,i)+=L(n,X)*R(X,i)",oo_aux_tensor,l1_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr9 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=A(n,i)*T(a,n)",gamma%prime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr10 =",ierr
           ! eq c13 + c14 :: + \ld mneX \rs Xn ( \td aeim - \ts ei \ts am)
           ierr=talsh_tensor_init(vvoo_aux_tensor)
           ierr=talsh_tensor_contract("Tau(a,e,i,m)+=T(a,e,i,m)*K()",vvoo_aux_tensor,t2_tensor,one_tensor)
           if(ierr.ne.0) write(*,*) "ierr11-a =",ierr
           ierr=talsh_tensor_contract("Tau(e,a,i,m)+=T(e,i)*T(a,m)",vvoo_aux_tensor,t1_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr11-b =",ierr
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(m,e)+=L(m,n,e,X)*R(X,n)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr12 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(m,e)*Tau(e,a,i,m)",gamma%prime_vo,ov_aux_tensor,vvoo_aux_tensor)
           if(ierr.ne.0) write(*,*) "ierr13 =",ierr

           ! eq c10       :: - 1/2 \ld mneX \rs Xi \td eamn
           ierr=talsh_tensor_init(oovo_aux_tensor)
           ierr=talsh_tensor_contract("A(m,n,e,i)+=L(m,n,e,X)*R(X,i)",oovo_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr14 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=A(m,n,e,i)*T(e,a,m,n)",gamma%prime_vo,oovo_aux_tensor,t2_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr15 =",ierr

           ! eq c5        :: +\ts ai
           ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr16 =",ierr

           ! eq c12       :: - 1/2 \ld mneX \rd aXmn \ts ei
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,e)+=L(m,n,e,X)*R(a,X,m,n)",vv_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr17 =",ierr
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("G(a,i)+=X(a,e)*T(e,i)",gamma%prime_vo,vv_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr18 =",ierr

           ! eq c11       :: -  \ld mneX \rd eXin \ts am
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("X(m,i)+=L(m,n,e,X)*R(e,X,i,n)",oo_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr19 =",ierr
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_contract("G(a,i)+=X(m,i)*T(a,m)",gamma%prime_vo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr20 =",ierr


    ! rho_ia
           ! eq d2        :: +\ld imaX \rs Xm
           ierr=talsh_tensor_contract("G(i,a)+=L(i,m,a,X)*R(X,m)",gamma%prime_ov,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr21 =",ierr

        end subroutine talsh_calc_dm_1rdm_ccsd_IP_mo



        subroutine talsh_calc_dm_2rdm_ccsd_IP_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)


          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor
          type(talsh_tens_t) :: ovvo_aux_tensor, ovov_aux_tensor,voov_aux_tensor
          type(talsh_tens_t) :: oovo_aux_tensor, vvov_aux_tensor

          type(talsh_tens_t) :: vX_aux_tensor


    !     tensor dimensions
          integer :: one_dims(1)=1
          integer, dimension(2) :: vX_dims
      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD =.true.

    !*********************************
    !     construct auxiliary tensor *
    !*********************************

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)


          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
          ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vo_aux_tensor,C8,vo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vovo_aux_tensor,C8,vovo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vooo_aux_tensor,C8,vooo_dims,init_val=ZERO)

          vX_dims = nvir
          vX_dims(2) = nvir
          ierr=talsh_tensor_construct(vX_aux_tensor,C8,vX_dims,init_val=ZERO)

    !********************
    ! Intermediates     *
    !********************
            ierr=talsh_tensor_construct(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)

            ! Loic_Taille

            ierr=talsh_tensor_construct(MATone_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ovvv_aux_tensor,C8,ovvv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ovov_aux_tensor,C8,ovov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(voov_aux_tensor,C8,voov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(vvov_aux_tensor,C8,vvov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ooov_aux_tensor,C8,ooov_dims,init_val=ZERO)

            ierr=talsh_tensor_construct(TauUP_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauDown_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauNo_tensor,C8,vvoo_dims,init_val=ZERO)

           ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr


    !!!!!!!!!!!!
    ! rho ijak !
    !!!!!!!!!!!!
    ! eq                 :: -(1/2) \ld ijaX \rs Xk
    ierr=talsh_tensor_contract("G(i,j,a,k)+=L(i,j,a,X)*R(X,k)",gamma%prime_oovo,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr1 =",ierr

    !!!!!!!!!!!!
    ! rho aibj !
    !!!!!!!!!!!!

!Loic_ Why no permutation \Pm ij \Pm ab ?           ! eq f2              :: - 1/4  \ld mibX \rd aXmj
    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,b)*R(e,a,j,m)",gamma%prime_voov,l2_tensor,r2_tensor,scale= MINUS_ONE_QUARTER)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,b)*R(e,a,i,m)",gamma%prime_voov,l2_tensor,r2_tensor)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(i,m,e,a)*R(e,b,j,m)",gamma%prime_voov,l2_tensor,r2_tensor)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=L(j,m,e,a)*R(e,b,i,m)",gamma%prime_voov,l2_tensor,r2_tensor,scale= MINUS_ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr2 =",ierr

!Loic_ Why no permutation \Pm ab ?          ! eq f4              :: -1/4 \ld mibX \rs Xj \ts am
    ierr=talsh_tensor_init(oovo_aux_tensor) !Loic_ 090622 Sign Problem
    ierr=talsh_tensor_contract("A(m,i,b,j)+=L(m,i,b,X)*R(X,j)",oovo_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr3 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(m,i,b,j)*T(a,m)",gamma%prime_voov,oovo_aux_tensor,t1_tensor,scale=MINUS_ONE_QUARTER)
!   ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,m,a,j)*T(b,m)",gamma%prime_voov,oovo_aux_tensor,t1_tensor,scale=ONE_QUARTER)
    if(ierr.ne.0) write(*,*) "ierr4 =",ierr

!Loic_ Why no permutation \Pm ij \Pm ab ?          ! eq f5              :: - 1/4  \ld imbX \rs Xm \ts aj
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(i,b)+=L(i,m,b,X)*R(X,m)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr5 =",ierr
    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,b)*,T(a,j)",gamma%prime_voov,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(i,a)*,T(b,j)",gamma%prime_voov,ov_aux_tensor,t1_tensor)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,b)*,T(a,i)",gamma%prime_voov,ov_aux_tensor,t1_tensor)
!    ierr=talsh_tensor_contract("G(a,i,j,b)+=A(j,a)*,T(b,i)",gamma%prime_voov,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr6 =",ierr

  !!!!!!!!!!!!
  ! rho abci !
  !!!!!!!!!!!!

! \Pm ab ? ! eq g4              ::  -  (\Pm ab) \ld mneXc \rd Xain \ts bm
    ierr=talsh_tensor_init(ovvo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,c,a,i)+=L(m,n,X,c)*R(X,a,n,i)",ovvo_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr7 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,a,c,i)*T(b,m)",gamma%prime_vvvo,ovvo_aux_tensor,t1_tensor,scale=MINUS_ONE)
!    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,b,c,i)*T(a,m)",gamma%prime_vvvo,ovvo_aux_tensor,t1_tensor)
    if(ierr.ne.0) write(*,*) "ierr8 =",ierr
    ! eq g8+g13          :: + \ld mncX \rs Xn \TauUp abmi  
    ierr=talsh_tensor_init(ov_aux_tensor)
    ierr=talsh_tensor_contract("A(m,c)+=L(m,n,c,X)*R(X,n)",ov_aux_tensor,l2_tensor,r1_tensor)
    if(ierr.ne.0) write(*,*) "ierr9 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,c)*TauUp(a,b,m,i)",gamma%prime_vvvo,ov_aux_tensor,TauUp_tensor)
    if(ierr.ne.0) write(*,*) "ierr10 =",ierr
    ! eq g7+g11          :: - \ld mncX \rs Xi \TauUp abmn        
    ierr=talsh_tensor_init(oovo_aux_tensor)
    ierr=talsh_tensor_contract("A(m,n,c,i)+=L(m,n,c,X)*R(X,i)",oovo_aux_tensor,l2_tensor,r1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr11 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(m,n,c,i)*TauUp(a,b,m,n)",gamma%prime_vvvo,oovo_aux_tensor,TauUp_tensor,scale=ONE_HALF)
    if(ierr.ne.0) write(*,*) "ierr12 =",ierr


! \Pm ab ? ! eq g10          :: + (\Pm ab) \ld mncX \rd aXmn \ts bi
    ierr=talsh_tensor_init(vv_aux_tensor)
    ierr=talsh_tensor_contract("A(a,c)+=L(m,n,c,X)*R(a,X,m,n)",vv_aux_tensor,l2_tensor,r2_tensor)
    if(ierr.ne.0) write(*,*) "ierr13 =",ierr
    ierr=talsh_tensor_contract("G(a,b,c,i)+=A(a,c)*T(b,i)",gamma%prime_vvvo,vv_aux_tensor,t1_tensor)
    ierr=talsh_tensor_contract("G(a,b,c,i)+=M(b,c)*T(a,i)",gamma%prime_vvvo,vv_aux_tensor,t1_tensor,scale=MINUS_ONE)
    if(ierr.ne.0) write(*,*) "ierr14 =",ierr

 !!!!!!!!!!!!
 ! rho iajk !
 !!!!!!!!!!!!
     ! eq k1              :: - (1/4) \ls iX \rd aXkj
     ierr=talsh_tensor_contract("G(i,a,j,k)+=L(i,X)*R(a,X,k,j)",gamma%prime_ovoo,l1_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr15 =",ierr
     ! eq k2              :: - (1/4) \Pm jk \ls iX \rs Xj \ts ak
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("A(i,j)+=L(i,X)*R(X,j)",oo_aux_tensor,l1_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr16 =",ierr
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,j)*T(a,k)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr17 =",ierr
     ! eq k4              :: + (1/4) \Pm jk \ld mieX \rd aXmj \ts ek
     ierr=talsh_tensor_init(vovo_aux_tensor)
     ierr=talsh_tensor_contract("A(a,i,e,j)+=L(m,i,e,X)*R(a,X,m,j)",vovo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr18 =",ierr
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(a,i,e,j)*T(e,k)",gamma%prime_ovoo,vovo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(a,i,e,k)*T(e,j)",gamma%prime_ovoo,vovo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr19 =",ierr
     ! eq k10             :: - (1/4) \Pm \ld mieX \rd eXmj \ts ak
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("A(i,j)+=L(m,i,e,X)*R(e,X,m,j)",oo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr20 =",ierr
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,j)*T(a,k)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr21 =",ierr
     
     ! eq k13 + k8             :: - (1/4) \ld imeX \rs Xm \TauDown eajk
     ierr=talsh_tensor_init(ov_aux_tensor)
     ierr=talsh_tensor_contract("A(i,e)+=L(i,m,e,X)*R(X,m)",ov_aux_tensor,l2_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr22 =",ierr
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,e)*T(e,a,j,k)",gamma%prime_ovoo,ov_aux_tensor,TauDown_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr23 =",ierr

     ! eq k5             :: + (1/4) \ld imeX \rd eXkj \ts am
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("A(m,i,k,j)+=L(m,i,e,X)*R(e,X,k,j)",oooo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr24 =",ierr
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(m,i,k,j)*T(a,m)",gamma%prime_ovoo,oooo_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr25 =",ierr

     ! eq k16+k12         :: + (1/4) \Pm kj \ld imeX \rs Xk \TauNo aemj
     ierr=talsh_tensor_init(oovo_aux_tensor)
     ierr=talsh_tensor_contract("A(i,m,e,k)+=L(m,i,e,X)*R(X,k)",oovo_aux_tensor,l2_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr26 =",ierr  
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,m,e,k)*T(a,e,m,j)",gamma%prime_ovoo,oovo_aux_tensor,TauNo_tensor)
     ierr=talsh_tensor_contract("G(i,a,j,k)+=A(i,m,e,j)*T(a,e,m,k)",gamma%prime_ovoo,oovo_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr27 =",ierr         

        !!!!!!!!!!!!
        ! rho ijkl !
        !!!!!!!!!!!!

     ! eq n1              :: + (1/4) \ld ijeX \rd eXkl
     ierr=talsh_tensor_contract("G(i,j,k,l)+=L(i,j,e,X)*R(e,X,k,l)",gamma%prime_oooo,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr28 =",ierr
     ! eq n2              :: + (1/4) \Pm kl \ld ijeX \rs Xl \ts ek
     ierr=talsh_tensor_init(oovo_aux_tensor)
     ierr=talsh_tensor_contract("A(i,j,e,l)+=L(i,j,e,X)*R(X,l)",oovo_aux_tensor,l2_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr29 =",ierr
     ierr=talsh_tensor_contract("G(i,j,k,l)+=A(i,j,e,l)*T(e,k)",gamma%prime_oooo,oovo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(i,j,k,l)+=A(i,j,e,k)*T(e,l)",gamma%prime_oooo,oovo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr30 =",ierr


        !!!!!!!!!!!!
        ! rho abij !
        !!!!!!!!!!!!

     ! eq l1+l2           :: + (1/4) \ls mX \rs Xm (\td abij + 1/2 \Pm ab \Pm ij \ts ai \ts bj)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,b,i,j)*K()",gamma%prime_vvoo,t2_tensor,one_tensor)
     if(ierr.ne.0) write(*,*) "ierr000 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr000 =",ierr

     ! eq l0+l00          :: + (1/4) 1/2 \ld nmeX \rd eXnm (\td abij + 1/2 \Pm  ab \Pm ij \ts ai \ts bj)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,b,i,j)*K()",gamma%prime_vvoo,t2_tensor,one_tensor)
     if(ierr.ne.0) write(*,*) "ierr0000 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr0000 =",ierr

!!!!!!!!!
!!!!!!!!!!!!!                            Take Care with TauUP
!!!!!!!!!
     ! eq l4+l5           :: - (1/4) \Pm ij \ls mX \rs Xj (\ts abim + \Pm ab  \ts ai \ts bm)
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("A(m,j)+=L(m,X)*T(X,j)",oo_aux_tensor,l1_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr31 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,j)*T(a,b,i,m)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(m,i)*T(a,b,j,m)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor)
     if(ierr.ne.0) write(*,*) "ierr32 =",ierr

     ! eq l3              :: + (1/4) \Pm ij \Pm ab \ls mX \rd aXim \ts bj
     ierr=talsh_tensor_init(vo_aux_tensor)
     ierr=talsh_tensor_contract("A(a,i)+=L(m,X)*R(a,X,i,m)",vo_aux_tensor,l1_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr33 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr34 =",ierr

     ! eq l17             :: - (1/4) \Pm ab \ls mX \rd bXji \ts am 
     ierr=talsh_tensor_init(vooo_aux_tensor)
     ierr=talsh_tensor_contract("A(b,m,j,i)+=L(m,X)*R(b,X,j,i)",vooo_aux_tensor,l1_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr35 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,m,j,i)*T(a,m)",gamma%prime_vvoo,vooo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,m,j,i)*T(b,m)",gamma%prime_vvoo,vooo_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr36 =",ierr

     ! eq l33+l35         :: - (1/4) \Pm ij \Pm ab \ld mneX \rs Xm (\td aeni+ \ts an \ts ei ) \ts bj
     ierr=talsh_tensor_init(ov_aux_tensor)
     ierr=talsh_tensor_contract("A(n,e)+=L(m,n,e,X)*R(X,m)",ov_aux_tensor,l2_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr37 =",ierr
     ierr=talsh_tensor_init(vo_aux_tensor)
     ierr=talsh_tensor_contract("B(a,i)+=A(n,e)*TauNo(a,e,n,i)",vo_aux_tensor,ov_aux_tensor,TauNo_tensor)
     if(ierr.ne.0) write(*,*) "ierr38 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr39 =",ierr

     ! eq l36+l37         :: - (1/4) \ld mneX \rs Xm (\P ab \ts an \td ebij + \Pmij \ts ei \td abnj)
!! DO NOT CHANGE contraction l2_tensor*r1_tensor
!           ierr=talsh_tensor_init(ov_aux_tensor)
!           ierr=talsh_tensor_contract("A(n,e)+=L(m,n,e,X)*R(X,m)",ov_aux_tensor,l2_tensor,r1_tensor)
!           if(ierr.ne.0) write(*,*) "ierr37 =",ierr
     ierr=talsh_tensor_init(vv_aux_tensor)          
     ierr=talsh_tensor_contract("B(a,e)+=A(n,e)*T(a,n)",vv_aux_tensor,ov_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr40 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,e)*T(e,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,e)*T(e,a,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
     if(ierr.ne.0) write(*,*) "ierr41 =",ierr
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("B(n,i)+=A(n,e)*T(e,i)",oo_aux_tensor,ov_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr42 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,i)*T(a,b,n,j)",gamma%prime_vvoo,oo_aux_tensor,t2_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,j)*T(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,t2_tensor)
     if(ierr.ne.0) write(*,*) "ierr43 =",ierr


     ! eq l24+l26         :: + (1/4) 1/2 \Pm ij \ld nmeX \rs Xi (\td abmn + \Pm ab \ts am \ts bn ) \ts ej
     ierr=talsh_tensor_init(oovo_aux_tensor)
     ierr=talsh_tensor_contract("A(n,m,e,i)+=L(m,n,e,X)*R(X,i)",oovo_aux_tensor,l2_tensor,r1_tensor)
     if(ierr.ne.0) write(*,*) "ierr44 =",ierr
     ierr=talsh_tensor_init(vvov_aux_tensor)
     ierr=talsh_tensor_contract("B(a,b,i,e)+=A(n,m,e,i)*TauUp(a,b,m,n)",vvov_aux_tensor,oovo_aux_tensor,TauUp_tensor)
     if(ierr.ne.0) write(*,*) "ierr45 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,b,i,e)*T(e,j)",gamma%prime_vvoo,vvov_aux_tensor,t1_tensor,scale=ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,b,j,e)*T(e,i)",gamma%prime_vvoo,vvov_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr46 =",ierr

     ! eq l22+l29         :: - (1/4) \Pm ab \Pm ij  \ld mneX \rs Xi ( \ts am \td bejn + 1/2 \ts bj \td eanm ) 
!!! DO NOT CHANGE contraction l2_tensor*r1_tensor
!!!           ierr=talsh_tensor_init(oovo_aux_tensor)
!!!           ierr=talsh_tensor_contract("A(n,m,e,i)+=L(m,n,e,X)*R(X,i)",oovo_aux_tensor,l2_tensor,r1_tensor)
     ierr=talsh_tensor_contract("B(n,a,e,i)+=A(n,m,e,i)*T(a,m)",ovvo_aux_tensor,oovo_aux_tensor,t1_tensor)
     if(ierr.ne.0) write(*,*) "ierr47 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,a,e,i)*T(b,e,j,n)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,b,e,i)*T(a,e,j,n)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,a,e,j)*T(b,e,i,n)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,b,e,j)*T(a,e,i,n)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr48 =",ierr
     ierr=talsh_tensor_init(vo_aux_tensor)
     ierr=talsh_tensor_contract("B(a,i)+=A(n,m,e,i)*T(e,a,n,m)",vo_aux_tensor,oovo_aux_tensor,t2_tensor)
     if(ierr.ne.0) write(*,*) "ierr49 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr50 =",ierr

     ! eq l6+l7           :: - (1/4) 1/2 \Pm ab \ld nmeX \rd aXnm (\td  ebij + \Pm ij \ts ei \ts bj )
     ierr=talsh_tensor_init(vv_aux_tensor)
     ierr=talsh_tensor_contract("A(a,e)+=L(n,m,e,X)*R(a,X,n,m)",vv_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr51 =",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,e)*T(e,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,TauDown_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,e)*T(e,a,i,j)",gamma%prime_vvoo,vv_aux_tensor,TauDown_tensor,scale=ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr52 =",ierr

     ! eq l9+l8           :: - (1/4) \Pm ij \ld nmeX \rd eXim (\td abnj + \Pm ab \ts ai \ts bj)
     ierr=talsh_tensor_init(oo_aux_tensor)
     ierr=talsh_tensor_contract("A(n,i)+=L(n,m,e,X)*R(e,X,i,m)",oo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr53 = ",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,i)*T(a,b,n,j)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,j)*T(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor)
     if(ierr.ne.0) write(*,*) "ierr54 = ",ierr

     ! eq l20, l21        :: - (1/4) \Pm ij \Pm ab \ld nmeX \rd aXim (\td bejn + \ts bj \ts en)
     ierr=talsh_tensor_init(ovvo_aux_tensor)
     ierr=talsh_tensor_contract("A(n,a,e,i)+=L(m,n,e,X)*R(a,X,i,m)",ovvo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr55 = ",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,a,e,i)*T(b,e,j,n)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,b,e,i)*T(a,e,j,n)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,a,e,j)*T(b,e,i,n)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,b,e,j)*T(a,e,i,n)",gamma%prime_vvoo,ovvo_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
     if(ierr.ne.0) write(*,*) "ierr56 = ",ierr

     ! eq l31             :: - (1/4) 1/2 \Pm ab \ld nmeX \rd bXji \td aemn
     ierr=talsh_tensor_init(vX_aux_tensor)
     ierr=talsh_tensor_contract("A(a,X)+=L(n,m,e,X)*T(a,e,m,n)",vX_aux_tensor,l2_tensor,t2_tensor)
     if(ierr.ne.0) write(*,*) "ierr57 = ",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,X)*R(b,X,j,i)",gamma%prime_vvoo,vX_aux_tensor,r2_tensor,scale=MINUS_ONE_HALF)
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,X)*R(a,X,j,i)",gamma%prime_vvoo,vX_aux_tensor,r2_tensor,scale=ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr58 = ",ierr

     ! eq l14+l12         :: + (1/4) 1/2 \ld nmeX \rd eXij (\td abmn +\Pm  ab \ts an \ts bm) 
     ierr=talsh_tensor_init(oooo_aux_tensor)
     ierr=talsh_tensor_contract("A(n,m,i,j)+=L(n,m,e,X)*R(e,X,i,j)",oooo_aux_tensor,l2_tensor,r2_tensor)
     if(ierr.ne.0) write(*,*) "ierr59 = ",ierr
     ierr=talsh_tensor_contract("G(a,b,i,j)+=A(n,m,i,j)*T(a,b,m,n)",gamma%prime_vvoo,oooo_aux_tensor,TauUp_tensor,scale=ONE_HALF)
     if(ierr.ne.0) write(*,*) "ierr60 = ",ierr


         end subroutine talsh_calc_dm_2rdm_ccsd_IP_mo

        subroutine talsh_calc_dm_1rdm_ccsd_EA_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)

          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_EA_vv_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor
          type(talsh_tens_t) :: vovv_aux_tensor

    !     tensor dimensions
          integer :: one_dims(1)=1
      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

    !*********************************
    !     construct auxiliary tensor *
    !*********************************


          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
            vvoo_dims(1:2) = nvir
            vvoo_dims(3:4) = nocc
            ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
            vv_dims = nvir
            ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
            oo_dims = nocc
            ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
           
            ov_dims    = nvir
            ov_dims(1) = nocc
            ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)

    !********************
    ! Intermediates     *
    !********************

            ! Loic_Taille

            vv_dims = nvir
            ierr=talsh_tensor_construct(MATone_EA_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
            oo_dims = nocc
            ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)

            vovv_dims = nvir
            vovv_dims(2) = nocc
            ierr=talsh_tensor_construct(vovv_aux_tensor,C8,vovv_dims,init_val=ZERO)

            vvoo_dims(1:2) = nvir
            vvoo_dims(3:4) = nocc
            ierr=talsh_tensor_construct(TauUP_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauDown_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauNo_tensor,C8,vvoo_dims,init_val=ZERO)

            !  Loic_Tenseurs_Intermediaires


           !  MATone_oo
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,j)*L(i,e)",MATone_oo_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,f,j,m)*L(i,m,e,f)",MATone_oo_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,j)*L(i,e)",MATtwo_oo_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(i,j)+=T(e,f,j,m)*L(i,m,e,f)",MATtwo_oo_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ! MATone_EA_vv
           ierr=talsh_tensor_contract("MATone(a,b)+=R(a,X)*L(X,b)",MATone_EA_vv_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(e,a,X,n)*L(X,n,e,b)",MATone_EA_vv_tensor,r2_tensor,l2_tensor)!,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(a,m)*L(m,b)",MATtwo_vv_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(e,a,m,n)*L(m,n,e,b)",MATtwo_vv_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           !  Tau
           
           ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauUp(a,b,i,j)+=T(a,b,i,j)",TauUp_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(a,i)*T(b,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(b,i)*T(a,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

    !********************
    !     calculate DMs *
    !********************


    ! rho ij
           ! eq a3 :: - \ld iXef \rd efjX   
           ierr=talsh_tensor_contract("G(i,j)+=L(i,X,e,f)*R(e,f,j,X)",gamma%prime_oo,l2_tensor,r2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr1 =",ierr
           ! eq a5      :: - \ld iXef \rs fX \ts ej
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(i,e)+=L(i,X,e,f)*R(f,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr2 =",ierr
           ierr=talsh_tensor_contract("G(i,j)+=X(i,e)*T(e,j)",gamma%prime_oo,ov_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr3 =",ierr

    ! rho ab
           ! eq b3 + b4 :: +MATone ab
           ierr=talsh_tensor_add("G(a,b)+=X(a,b)",gamma%prime_vv,MATone_EA_vv_tensor)
           if(ierr.ne.0) write(*,*) "ierr4 =",ierr
           ! eq b5      :: +\ld mneb \rs em \ts an
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(n,b)+=L(X,n,e,b)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr5 =",ierr
           ierr=talsh_tensor_contract("G(a,b)+=X(n,b)*T(a,n)",gamma%prime_vv,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr6 =",ierr


    ! rho ai
           ! eq c11  :: -\ld Xnef \rd efXi \ts an 
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("X(n,i)+=L(X,n,e,f)*R(e,f,X,i)",oo_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr7 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(n,i)*T(a,n)",gamma%prime_vo,MATone_oo_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr8 =",ierr



           ! eq c8 + c12  :: -\MATone af \ts fi
           ierr=talsh_tensor_contract("G(a,i)+=X(a,f)*T(f,i)",gamma%prime_vo,MATone_EA_vv_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr9 =",ierr

           ! eq c13 + c14 :: - \ld Xnef \rs eX \tauNo fain
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(n,f)+=L(X,n,e,f)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr10 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(n,f)*TauNo(f,a,i,n)",gamma%prime_vo,ov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr11 =",ierr
           ! eq c5        :: +\ts ai
           ierr=talsh_tensor_add("G(a,i)+=T(a,i)",gamma%prime_vo,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr12 =",ierr
           ! eq c6        :: +\ls Xe \rd aeiX
           ierr=talsh_tensor_contract("G(a,i)+=L(X,e)*R(a,e,i,X)",gamma%prime_vo,l1_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr13 =",ierr
           ! eq c10       :: - 1/2 \ld Xnef \rs aX \td efin
           ierr=talsh_tensor_init(vovv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,n,e,f)+=L(X,n,e,f)*R(a,X)",vovv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr14 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X(a,n,e,f)*T(e,f,i,n)",gamma%prime_vo,vovv_aux_tensor,t2_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr15 =",ierr


    ! rho_ia
           ! eq d2        :: +\ld iXae \rs eX
           ierr=talsh_tensor_contract("G(i,a)+=L(i,X,a,e)*R(e,X)",gamma%prime_ov,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr16 =",ierr


    !  Symmetrization


           ierr=talsh_tensor_add("G(i,a)+=X(i,a)",gamma%symm_ov,gamma%prime_ov,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr17 =",ierr
           ierr=talsh_tensor_contract("G(i,a)+=X+(a,i)*K()",gamma%symm_ov,gamma%prime_vo,one_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr18 =",ierr
           ierr=talsh_tensor_add("G(a,i)+=X(a,i)",gamma%symm_vo,gamma%prime_vo,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr19 =",ierr
           ierr=talsh_tensor_contract("G(a,i)+=X+(i,a)*K()",gamma%symm_vo,gamma%prime_ov,one_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr20 =",ierr


    !     debug
          if (.true.) then
            print*, "gprime_oo  = ", tensor_norm2(gamma%prime_oo)
            if (.not.CCD) then
              print*, "gprime_ov  = ", tensor_norm2(gamma%prime_ov)
              print*, "gsymm_ov   = ", tensor_norm2(gamma%symm_ov)
              print*, "gprime_vo  = ", tensor_norm2(gamma%prime_vo)
              print*, "gsymm_vo   = ", tensor_norm2(gamma%symm_vo)
            end if
            print*, "gprime_vv  = ", tensor_norm2(gamma%prime_vv)
            print*, "gsymm_oo   = ", tensor_norm2(gamma%symm_oo)
            print*, "gsymm_vv   = ", tensor_norm2(gamma%symm_vv)
          end if

    !**************
    !     cleanup *
    !**************

          ierr=talsh_tensor_destruct(one_tensor)
          if (.not.CCD) then
            ierr=talsh_tensor_destruct(tau2_tensor)
            ierr=talsh_tensor_destruct(vv_aux_tensor)
            ierr=talsh_tensor_destruct(oo_aux_tensor)
          end if


        end subroutine talsh_calc_dm_1rdm_ccsd_EA_mo


        subroutine talsh_calc_dm_2rdm_ccsd_EA_mo(gamma,t2_tensor,l2_tensor,nocc,nvir,t1_tensor,l1_tensor,r1_tensor,r2_tensor)


          real(8)                        :: print_level=2

    !     input variables
          type(talsh_dm_tens), intent(inout) :: gamma
          type(talsh_tens_t),  intent(inout) :: t2_tensor, l2_tensor
          type(talsh_tens_t),  intent(inout), optional :: t1_tensor, l1_tensor

          type(talsh_tens_t),  intent(inout), optional :: r1_tensor, r2_tensor
    !Loic_

          integer,             intent(in   ) :: nocc, nvir
    !     auxiliary tensors
          type(talsh_tens_t) :: one_tensor
          type(talsh_tens_t) :: tau2_tensor
          type(talsh_tens_t) :: oo_aux_tensor, vv_aux_tensor

          ! Loic_tensor 
          type(talsh_tens_t) :: TauUp_tensor,TauDown_tensor, TauNo_tensor
          type(talsh_tens_t) :: MATone_oo_tensor, MATone_vv_tensor
          type(talsh_tens_t) :: MATone_EA_oo_tensor, MATone_EA_vv_tensor
          type(talsh_tens_t) :: MATtwo_EA_oo_tensor
          type(talsh_tens_t) :: MATtwo_oo_tensor, MATtwo_vv_tensor

          type(talsh_tens_t) :: ov_aux_tensor, vo_aux_tensor
          type(talsh_tens_t) :: vvvv_aux_tensor, oooo_aux_tensor,oovv_aux_tensor
          type(talsh_tens_t) :: vvoo_aux_tensor, vvvo_aux_tensor,ooov_aux_tensor
          type(talsh_tens_t) :: vooo_aux_tensor, ovvv_aux_tensor,vovo_aux_tensor
          type(talsh_tens_t) :: ovvo_aux_tensor, ovov_aux_tensor,voov_aux_tensor
          type(talsh_tens_t) :: vovv_aux_tensor

    !     tensor dimensions
          integer :: one_dims(1)=1
          integer, dimension(2) :: Xo_dims, Xv_dims
      integer(INTD), dimension(2) :: oo_dims, vv_dims, vo_dims, ov_dims
      integer(INTD), dimension(4) :: vvoo_dims

      integer(INTD), dimension(4) :: ooov_dims, ovov_dims, vvvo_dims
      integer(INTD), dimension(4) :: vovv_dims, vvvv_dims
      integer(INTD), dimension(4) :: voov_dims, oovv_dims
      integer(INTD), dimension(4) :: ovoo_dims, ovvo_dims, vooo_dims
      integer(INTD), dimension(4) :: oooo_dims, vovo_dims, ovvv_dims
      integer(INTD), dimension(4) :: oovo_dims, vvov_dims

    !     CCD switch
          logical :: CCD = .false.
    !     error code
          integer :: ierr

          if (.not.present(t1_tensor) .and. .not.present(l1_tensor)) CCD = .true.

    !*********************************
    !     construct auxiliary tensor *
    !*********************************

      oo_dims        = (/nocc,nocc/)
      vv_dims        = (/nvir,nvir/)
      vo_dims        = (/nvir,nocc/)
      ov_dims        = (/nocc,nvir/)
      ooov_dims      = (/nocc,nocc,nocc,nvir/)
      ovov_dims      = (/nocc,nvir,nocc,nvir/)
      vvvv_dims      = (/nvir,nvir,nvir,nvir/)
      vvvo_dims      = (/nvir,nvir,nvir,nocc/)
      vovv_dims      = (/nvir,nocc,nvir,nvir/)
      voov_dims      = (/nvir,nocc,nocc,nvir/)
      oovv_dims      = (/nocc,nocc,nvir,nvir/)
      ovoo_dims      = (/nocc,nvir,nocc,nocc/)
      ovvo_dims      = (/nocc,nvir,nvir,nocc/)
      oooo_dims      = (/nocc,nocc,nocc,nocc/)
      vvoo_dims      = (/nvir,nvir,nocc,nocc/)
      ovvv_dims      = (/nocc,nvir,nvir,nvir/)
      vooo_dims      = (/nvir,nocc,nocc,nocc/)
      oovo_dims      = (/nocc,nocc,nvir,nocc/)
      vovo_dims      = (/nvir,nocc,nvir,nocc/)
      vvov_dims      = (/nvir,nvir,nocc,nvir/)


          ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
          ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vv_aux_tensor,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(oo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(ov_aux_tensor,C8,ov_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(vo_aux_tensor,C8,vo_dims,init_val=ZERO)

    !********************
    ! Intermediates     *
    !******************** 
            ierr=talsh_tensor_construct(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(oooo_aux_tensor,C8,oooo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(vovv_aux_tensor,C8,vovv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(vooo_aux_tensor,C8,vooo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(vvvv_aux_tensor,C8,vvvv_dims,init_val=ZERO)

!          ierr=talsh_tensor_construct(gamma%prime_oo,C8,oo_dims,init_val=ZERO)
!          ierr=talsh_tensor_construct(gamma%prime_ov,C8,ov_dims,init_val=ZERO)
!          ierr=talsh_tensor_construct(gamma%prime_vo,C8,vo_dims,init_val=ZERO)
!          ierr=talsh_tensor_construct(gamma%prime_vv,C8,vv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vovo,C8,vovo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvvo,C8,vvvo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vovv,C8,vovv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvvv,C8,vvvv_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_ovoo,C8,ovoo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(gamma%prime_vvoo,C8,vvoo_dims,init_val=ZERO)


            ! Loic_Taille

            ierr=talsh_tensor_construct(MATone_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_vv_tensor,C8,vv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATone_oo_tensor,C8,oo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(MATtwo_oo_tensor,C8,oo_dims,init_val=ZERO)

            Xo_dims = nocc
            ierr=talsh_tensor_construct(MATtwo_EA_oo_tensor,C8,Xo_dims,init_val=ZERO)
            Xv_dims = nvir
            ierr=talsh_tensor_construct(MATone_EA_vv_tensor,C8,Xv_dims,init_val=ZERO)

            ierr=talsh_tensor_construct(ovvv_aux_tensor,C8,ovvv_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(ovov_aux_tensor,C8,ovov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(voov_aux_tensor,C8,voov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauUP_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauDown_tensor,C8,vvoo_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(TauNo_tensor,C8,vvoo_dims,init_val=ZERO)

            !  Loic_Tenseurs_Intermediaires


           !  MATone_oo
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,j)*L(i,e)",MATone_oo_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(i,j)+=R(e,f,j,m)*L(i,m,e,f)",MATone_oo_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           ierr=talsh_tensor_contract("MATtwo(X,j)+=T(e,j)*L(X,e)",MATtwo_EA_oo_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(X,j)+=T(e,f,j,m)*L(X,m,e,f)",MATtwo_EA_oo_tensor,t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           ierr=talsh_tensor_contract("MATone(a,b)+=R(a,m)*L(m,b)",MATone_vv_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATone(a,b)+=R(e,a,m,n)*L(m,n,e,b)",MATone_vv_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(a,m)*L(m,b)",MATtwo_vv_tensor,t1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("MATtwo(a,b)+=T(e,a,m,n)*L(m,n,e,b)",MATtwo_vv_tensor,&
              t2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           ierr=talsh_tensor_contract("M(a,b)+=R(a,X)*L(X,b)",MATone_EA_vv_tensor,r1_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("M(a,b)+=R(e,a,X,n)*L(X,n,e,b)",MATone_EA_vv_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr

           !  Tau
           
           ierr=talsh_tensor_add("TauNo(a,b,i,j)+=T(a,b,i,j)",TauNo_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauNo(a,b,i,j)+=T(a,i)*T(b,j)",TauNo_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauUp(a,b,i,j)+=T(a,b,i,j)",TauUp_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(a,i)*T(b,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauUp(a,b,i,j)+=T(b,i)*T(a,j)",TauUp_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

           ierr=talsh_tensor_add("TauDown(a,b,i,j)+=T(a,b,i,j)",TauDown_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,i)*T(b,j)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr
           ierr=talsh_tensor_contract("TauDown(a,b,i,j)+=T(a,j)*T(b,i)",TauDown_tensor,t1_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr00 =",ierr

    !********************
    !     calculate DMs *
    !********************

        !!!!!!!!!!!!
        ! rho ijak      : No ijak for EA !
        !!!!!!!!!!!!

        !!!!!!!!!!!!
        ! rho aibj !
        !!!!!!!!!!!!

! \Pm ab \Pm ij           !  eq f2             :: - (1/4) \ld Xibe \rd eaXj
           ierr=talsh_tensor_contract("G(a,i,b,j)+=L(X,i,b,e)*R(e,a,X,j)",gamma%prime_vovo,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr1 =",ierr           

! \Pm ab \Pm ij           ! eq  f3             :: - (1/4) \ld imeb \rs am \ts ej
           ierr=talsh_tensor_init(ovvv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,i,b,e)+=L(i,X,e,b)*R(a,X)",vovv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr2 =",ierr
           ierr=talsh_tensor_contract("G(a,i,b,j)+=X(a,i,b,e)*T(e,j)",gamma%prime_vovo,vovv_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr3 =",ierr

! \Pm ab \Pm ij           ! eq f5              :: +/- (1/4) \ld imbe \rs em \ts aj
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(i,b)+=L(i,X,e,b)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr4 =",ierr
           ierr=talsh_tensor_contract("G(a,i,b,j)+=X(i,b)*T(a,j)",gamma%prime_vovo,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr5 =",ierr
           
        !!!!!!!!!!!!
        ! rho abci !
        !!!!!!!!!!!!

           ! eq g1              ::  + (1/2) \s mc \rd baim
           ierr=talsh_tensor_contract("G(a,b,c,i)+=L(X,c)*R(b,a,i,X)",gamma%prime_vvvo,l1_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr6 =",ierr

           ! eq g2 + g10        :: + (1/2) \Pm ab \MATone_EA_vv \ts bi
           ierr=talsh_tensor_contract("G(a,b,c,i)+=M(a,c)*T(b,i)",gamma%prime_vvvo,MATone_EA_vv_tensor,t1_tensor)
           ierr=talsh_tensor_contract("G(a,b,c,i)+=M(b,c)*T(a,i)",gamma%prime_vvvo,MATone_EA_vv_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr7 =",ierr

           ! eq g4              :: - (1/2) \Pm ab \ld mXec \rd eaiX \ts bm
           ierr=talsh_tensor_init(voov_aux_tensor)
           ierr=talsh_tensor_contract("X(m,a,i,c)+=L(m,X,e,c)*R(e,a,i,X)",ovov_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr8 =",ierr
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(m,a,i,c)*T(b,m)",gamma%prime_vvvo,ovov_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(m,b,i,c)*T(a,m)",gamma%prime_vvvo,ovov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr9 =",ierr

           ! eq g5              :: - (1/2) 1/2 \ld Xnec \rd baXn \ts ei
           ierr=talsh_tensor_init(vvvv_aux_tensor)
           ierr=talsh_tensor_contract("X(b,a,e,c)+=L(X,n,e,c)*R(b,a,X,n)",vvvv_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr10 =",ierr 
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(b,a,e,c)*T(e,i)",gamma%prime_vvvo,&
               vvvv_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(b,a,e,c)*T(e,i)",gamma%prime_vvvo,vvvv_aux_tensor,t1_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr11 =",ierr 

           ! eq g6 + g12        :: - (1/2) \Pm ab \ld Xnec \rs bX \TauNo eain
           ierr=talsh_tensor_init(vovv_aux_tensor)
           ierr=talsh_tensor_contract("X(b,n,e,c)+=L(X,n,e,c)*R(b,X)",vovv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr12 =",ierr  
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(b,n,e,c)*T(e,a,i,n)",gamma%prime_vvvo,&
               vovv_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(a,n,e,c)*T(e,b,i,n)",gamma%prime_vvvo,vovv_aux_tensor,TauNo_tensor)
           if(ierr.ne.0) write(*,*) "ierr13 =",ierr  

           ! eq g8 + g13        :: + (1/2) \ld Xnec \rs eX \TauUp abni
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(n,c)+=L(X,n,e,c)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr14 =",ierr
           ierr=talsh_tensor_contract("G(a,b,c,i)+=X(n,c)*T(a,b,n,i)",gamma%prime_vvvo,ov_aux_tensor,TauUp_tensor)
           if(ierr.ne.0) write(*,*) "ierr14-b =",ierr


        !!!!!!!!!!!!
        ! rho aibc !
        !!!!!!!!!!!!

           ! eq h1         :: + 1/2 \ld Xibc \rs aX
           ierr=talsh_tensor_contract("G(a,i,b,c)+=L(X,i,b,c)*R(a,X)",gamma%prime_vovv,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr15 =",ierr

        !!!!!!!!!!!!
        ! rho abcd !
        !!!!!!!!!!!!

           ! eq i1         :: + (1/4) \ld Xnec \rd abXn  i!!!!!! PROBLEM
           ierr=talsh_tensor_contract("G(a,b,c,d)+=L(X,n,c,d)*R(a,b,X,n)",gamma%prime_vvvv,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr16 =",ierr
           ! eq i2         :: + (1/4) \Pm ab \ld mncd \rs am \ts bn
           ierr=talsh_tensor_init(vovv_aux_tensor)
           ierr=talsh_tensor_contract("X(a,n,c,d)+=L(X,n,c,d)*R(a,X)",vovv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr17 =",ierr   
           ierr=talsh_tensor_contract("G(a,b,c,d)+=X(a,n,c,d)*T(b,n)",gamma%prime_vvvv,vovv_aux_tensor,t1_tensor)
           ierr=talsh_tensor_contract("G(a,b,c,d)+=X(b,n,c,d)*T(a,n)",gamma%prime_vvvv,vovv_aux_tensor,t1_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr18 =",ierr 

        !!!!!!!!!!!!
        ! rho ijab              : No rho !
        !!!!!!!!!!!!

        !!!!!!!!!!!!
        ! rho iajk !
        !!!!!!!!!!!!

           ! eq k4              :: +\Pm ij  \ld iXfe \rd fajX (i,a,j,e) \ts ek

           ierr=talsh_tensor_init(ovov_aux_tensor)
           ierr=talsh_tensor_contract("X(i,a,j,e)+=L(i,X,f,e)*R(f,a,j,X)",ovov_aux_tensor,l2_tensor,r2_tensor)
           ierr=talsh_tensor_contract("X(i,a,j,e)+=L(j,X,f,e)*R(f,a,i,X)",ovov_aux_tensor,l2_tensor,r2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr19 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=X(i,a,j,e)*T(e,k)",gamma%prime_ovoo,ovov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr20 =",ierr

           ! eq k7              :: + 1/2 \ld iXfe \rs aX (iafe) \td fejk
           ierr=talsh_tensor_init(ovvv_aux_tensor)
           ierr=talsh_tensor_contract("X(i,a,f,e)+=L(i,X,f,e)*R(a,X)",ovvv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr21 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=X(i,a,f,e)*T(f,e,j,k)",gamma%prime_ovoo,ovvv_aux_tensor,t2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr22 =",ierr

           ! eq k8              :: - \Pm jk \ld iXfe \rs eX \ts fj \ts ak
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("X(i,f)+=L(i,X,f,e)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr23 =",ierr
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("Y(i,j)+=X(i,f)*R(f,j)",oo_aux_tensor,ov_aux_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr24 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=Y(i,j)*T(a,k)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(i,a,j,k)+=Y(i,k)*T(a,j)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr25 =",ierr

           ! eq k10             :: - \Pm jk \ld iXfe \ts ak \rs fejX 
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("X(i,j)+=L(i,m,f,e)*R(f,e,j,m)",oo_aux_tensor,l2_Tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr26 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=X(i,j)*T(a,k)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(i,a,j,k)+=X(i,k)*T(a,j)",gamma%prime_ovoo,oo_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr27 =",ierr

           ! eq k11             :: + 1/2 \Pm jk \ld iXfe \ts fj \ts ek \rs aX
           ierr=talsh_tensor_init(ovvv_aux_tensor)
           ierr=talsh_tensor_contract("X(i,a,f,e)+=L(i,X,f,e)*R(a,X)",ovvv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr28 =",ierr
           ierr=talsh_tensor_init(ovov_aux_tensor)                                   
           ierr=talsh_tensor_contract("Y(i,a,j,e)+=X(i,a,f,e)*T(f,j)",ovov_aux_tensor,ovvv_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr29 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=Y(i,a,j,e)*T(e,k)",gamma%prime_ovoo,ovov_aux_tensor,t1_tensor,scale=ONE_HALF)
           ierr=talsh_tensor_contract("G(i,a,j,k)+=Y(i,a,k,e)*T(e,j)",gamma%prime_ovoo,&
               ovov_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr30 =",ierr   
   
           ! eq k13             :: - \ld iXfe \rs eX \td fajk
           ierr=talsh_tensor_init(ov_aux_tensor)                                   
           ierr=talsh_tensor_contract("X(i,f)+=L(i,X,f,e)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr31 =",ierr
           ierr=talsh_tensor_contract("G(i,a,j,k)+=Y(i,f)*T(f,a,j,k)",gamma%prime_ovoo,&
               ov_aux_tensor,t2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr32 =",ierr


        !!!!!!!!!!!!
        ! rho abij !
        !!!!!!!!!!!!
!XXX idem           ! eq h4+h5+h8+h9     :: - (1/4) \Pm ij  \MATone ni \TauUp abnj
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_contract("X(i,j)+=R(e,f,j,m)*L(i,m,e,f)",oo_aux_tensor,r2_tensor,l2_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr0 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=X(n,i)*TauUp(a,b,n,j)",gamma%prime_vvoo,&
               oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=X(n,j)*TauUp(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor)
           if(ierr.ne.0) write(*,*) "ierr66 =",ierr
!XXX idem          ! eq h6+h7+h10+h11   :: - (1/4) \Pm ab \MATone ae \TauDown ebij
           ierr=talsh_tensor_contract("G(a,b,i,j)+=M(a,e)*TauDown(e,b,i,j)",gamma%prime_vvoo,&
               MATone_EA_vv_tensor,TauDown_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=M(b,e)*TauDown(e,a,i,j)",gamma%prime_vvoo,MATone_EA_vv_tensor,TauDown_tensor)
           if(ierr.ne.0) write(*,*) "ierr67 =",ierr
!ok           ! eq j13+j15         :: + (1/4) 1/4 \ld Xnef \rd abXn \TauDown efij
           ierr=talsh_tensor_init(vvvv_aux_tensor)
           ierr=talsh_tensor_contract("A(a,b,e,f)+=L(X,n,e,f)*R(a,b,X,n)",vvvv_aux_tensor,l2_tensor,r2_tensor)
           if(ierr.ne.0) write(*,*) "ierr70 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,b,e,f)*TauDawn(e,f,i,j)",gamma%prime_vvoo,&
               vvvv_aux_tensor,TauDown_tensor,scale=ONE_QUARTER)
           if(ierr.ne.0) write(*,*) "ierr71 =",ierr
!ok           ! eq j25+j27         :: + (1/4) 1/2 \Pm ab \ld Xnef \rs aX \ts bn \TauDown efij
           ierr=talsh_tensor_init(vovv_aux_tensor)
           ierr=talsh_tensor_contract("A(a,n,e,f)+=L(X,n,e,f)*R(a,X)",vovv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr75 =",ierr
           ierr=talsh_tensor_init(vooo_aux_tensor)
           ierr=talsh_tensor_contract("B(a,n,i,j)+=A(a,n,e,f)*TauDown(e,f,i,j)",vooo_aux_tensor,vovv_aux_tensor,TauDown_tensor)
           if(ierr.ne.0) write(*,*) "ierr76 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,n,i,j)*T(b,n)",gamma%prime_vvoo,vooo_aux_tensor,t1_tensor,scale=ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,n,i,j)*T(a,n)",gamma%prime_vvoo,&
               vooo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr77 =",ierr
!ok           ! eq j35+j36         :: - (1/4) \Pm ij \ld Xnef \rs eX \ts fi \TauUp abnj
           ierr=talsh_tensor_init(oo_aux_tensor)
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("A(n,f)+=L(X,n,e,f)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr78 =",ierr
           ierr=talsh_tensor_contract("B(n,i)+=A(n,f)*T(f,i)",oo_aux_tensor,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr79 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,i)*TauUp(a,b,n,j)",gamma%prime_vvoo,&
               oo_aux_tensor,TauUp_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(n,j)*TauUp(a,b,n,i)",gamma%prime_vvoo,oo_aux_tensor,TauUp_tensor)
           if(ierr.ne.0) write(*,*) "ierr80 =",ierr
!ok           ! eq j20+j21         :: - (1/4) \Pm ij \Pm ab \ld Xnef \rd aeiX \TauNo bfnj
           ierr=talsh_tensor_init(voov_aux_tensor)
           ierr=talsh_tensor_contract("A(a,n,i,f)+=R(a,e,i,X)*L(X,n,e,f)",voov_aux_tensor,r2_tensor,l2_tensor)
           if(ierr.ne.0) write(*,*) "ierr83 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,i,f)*TauNo(b,f,n,j)",gamma%prime_vvoo,&
                voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,i,f)*TauNo(a,f,n,j)",&
                &gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,n,j,f)*TauNo(b,f,n,i)",&
                &gamma%prime_vvoo,voov_aux_tensor,TauNo_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,n,j,f)*TauNo(a,f,n,i)",gamma%prime_vvoo,&
               voov_aux_tensor,TauNo_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr84 =",ierr

!ok           ! eq j16+j30         :: - (1/4) \Pm ij \rd abin \MATtwo nj
           ierr=talsh_tensor_contract("G(a,b,i,j)+=R(a,b,i,X)*M(X,j)",gamma%prime_vvoo,&
               r2_tensor,MATtwo_EA_oo_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=R(a,b,j,X)*M(X,i)",&
               &gamma%prime_vvoo,r2_tensor,MATtwo_EA_oo_tensor)
           if(ierr.ne.0) write(*,*) "ierr85 =",ierr

!ok           ! eq j1              :: + (1/4) \td abij
           ierr=talsh_tensor_add("G(a,b,i,j)+=T(a,b,i,j)",gamma%prime_vvoo,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr87 =",ierr
!ok           ! eq j2              :: + (1/4) 1/2 \Pm ab \Pm ij \ts ai \ts bj
           ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,i)*T(b,j)",&
               &gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,i)*T(a,j)",&
               &gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=T(a,j)*T(b,i)",&
               &gamma%prime_vvoo,t1_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=T(b,j)*T(a,i)",&
               &gamma%prime_vvoo,t1_tensor,t1_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr88 =",ierr

!ok           ! eq j3              :: + (1/4) \ls Xe \rd aeiX \ts bj
           ierr=talsh_tensor_init(vo_aux_tensor)
           ierr=talsh_tensor_contract("A(a,i)+=R(a,e,i,X)*L(X,e)",vo_aux_tensor,r2_tensor,l1_tensor)
           if(ierr.ne.0) write(*,*) "ierr89 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=A(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr90 =",ierr

!ok ! Signe XX ! eq j23             :: - (1/4) 1/2 \Pm ij \Pm ab \ld mXef \rs aX \td efmi \ts bj
           ierr=talsh_tensor_init(ovvv_aux_tensor)
           ierr=talsh_tensor_contract("A(m,a,e,f)+=L(m,X,e,f)*R(a,X)",ovvv_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr93 =",ierr

           ierr=talsh_tensor_init(vo_aux_tensor)
           ierr=talsh_tensor_contract("B(a,i)+=A(m,a,e,f)*T(e,f,m,i)",vo_aux_tensor,ovvv_aux_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr93 =",ierr


           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",&
                & gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",&
                &gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",&
                &gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE_HALF)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",&
                &gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=ONE_HALF)
           if(ierr.ne.0) write(*,*) "ierr95 =",ierr

!ok           ! eq j28             :: - (1/4)  \Pm ij \Pm ab \ld mXef \rs aX \td bejm \ts fi
!!!!!!!! Don't init ovvv_aux_Tensor, keep the same value as before
!           ierr=talsh_tensor_init(ovvv_aux_tensor)
!           ierr=talsh_tensor_contract("A(m,a,e,f)+=L(m,X,e,f)*R(a,X)",ovvv_aux_tensor,l2_tensor,r1_tensor)

           ierr=talsh_tensor_init(ovvo_aux_tensor)
           ierr=talsh_tensor_contract("B(m,a,e,i)+=A(m,a,e,f)*T(f,i)",ovvo_aux_tensor,ovvv_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr97 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(m,a,e,i)*T(b,e,j,m)",gamma%prime_vvoo,&
               ovvo_aux_tensor,t2_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(m,b,e,i)*T(a,e,j,m)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(m,a,e,j)*T(b,e,i,m)",gamma%prime_vvoo,ovvo_aux_tensor,t2_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(m,b,e,j)*T(a,e,i,m)",gamma%prime_vvoo,&
               ovvo_aux_tensor,t2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr98 =",ierr


!ok           ! eq j33             :: + (1/4) \Pm ij \Pm ab \ld Xnef \rs eX \td afin \ts bj 
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_init(vo_aux_tensor)
           ierr=talsh_tensor_contract("A(n,f)+=L(X,n,e,f)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr104 =",ierr
           ierr=talsh_tensor_contract("B(a,i)+=A(n,f)*T(a,f,i,n)",vo_aux_tensor,ov_aux_tensor,t2_tensor)
           if(ierr.ne.0) write(*,*) "ierr105 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,i)*T(b,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,j)*T(b,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,i)*T(a,j)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(b,j)*T(a,i)",gamma%prime_vvoo,vo_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr106 =",ierr


!ok           ! eq 37              :: - (1/4) \Pm ab \ld Xnef \rs eX \td fbij \ts an
           ierr=talsh_tensor_init(vv_aux_tensor)
           ierr=talsh_tensor_init(ov_aux_tensor)
           ierr=talsh_tensor_contract("A(n,f)+=L(X,n,e,f)*R(e,X)",ov_aux_tensor,l2_tensor,r1_tensor)
           if(ierr.ne.0) write(*,*) "ierr110 =",ierr
           ierr=talsh_tensor_contract("B(a,f)+=A(n,f)*T(a,n)",vv_aux_tensor,ov_aux_tensor,t1_tensor)
           if(ierr.ne.0) write(*,*) "ierr111 =",ierr
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor)
           ierr=talsh_tensor_contract("G(a,b,i,j)+=B(a,f)*T(f,b,i,j)",gamma%prime_vvoo,vv_aux_tensor,t2_tensor,scale=MINUS_ONE)
           if(ierr.ne.0) write(*,*) "ierr112 =",ierr



        !!!!!!!!!!!!
        ! rho ijkl              : No rho ijkl for EA  !
        !!!!!!!!!!!!
                

            print *, "gamma%prime_vovo  = ", tensor_norm2(gamma%prime_vovo)
            print *, "gamma%prime_vvvo  = ", tensor_norm2(gamma%prime_vvvo)
            print *, "gamma%prime_vovv  = ", tensor_norm2(gamma%prime_vovv)
            print *, "gamma%prime_vvvv  = ", tensor_norm2(gamma%prime_vvvv)

            print *, "gamma%prime_ovoo  = ", tensor_norm2(gamma%prime_ovoo)
            print *, "gamma%prime_vvoo  = ", tensor_norm2(gamma%prime_vvoo)

        end subroutine talsh_calc_dm_2rdm_ccsd_EA_mo


end module 


