!   
!    Outputting TAL-SH Fortran code generated with :
!      Codegen, a toolset to process tensor contraction DSLs (SIAL etc)
!      Andre Gomes (CNRS UMR8523, Lille) and Dmitry Lyakh (OCLF, Oak Ridge)
!   
subroutine eti_CCD_codegen_call(nocc,nvir,S23,H13,H3,H2,H14,H1,H9,H10,Z26,Z27)

    use exacorr_tensor_interface
    use tensor_algebra
    use, intrinsic:: ISO_C_BINDING

    implicit none

    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE_HALF=(0.5D0,0.D0), &
                             MINUS_ONE=(-1.D0,0.D0),ONE=(1.0D0,0.D0), MINUS_ONE_HALF=(-0.5D0,0.D0), &
                             MINUS_ONE_QUARTER=(-0.25D0,0.D0), ONE_QUARTER=(0.25D0,0.D0), &
                             MINUS_ONE_EIGHT=(-0.125D0,0.D0)


    integer, intent(in) :: nocc
    integer, intent(in) :: nvir
    type(tensor_t), intent(inout) :: S23    !       VVOO
    type(tensor_t), intent(inout) :: H13    !       OOVV
    type(tensor_t), intent(inout) :: H3    !       VVOO
    type(tensor_t), intent(inout) :: H2    !       VVVV
    type(tensor_t), intent(inout) :: H14    !       OOOO
    type(tensor_t), intent(inout) :: H1    !       VV
    type(tensor_t), intent(inout) :: H9    !       OO
    type(tensor_t), intent(inout) :: H10    !       VOOV
    type(tensor_t), intent(inout) :: Z26    !       UnitTensor
    type(tensor_t), intent(inout) :: Z27    !       VVOO
    type(tensor_t) :: Z33_1    !       OOOO
    type(tensor_t) :: Z30_1    !       VV
    type(tensor_t) :: Z32_1    !       OO
    type(tensor_t) :: Z37_1    !       VOOV

    integer :: ierr

!   original expression: Z26(|)+=S23(d1bd2b|l1bl2b)*H13(l1bl2b|d1bd2b)
    ierr=t_einsum("Z26()+=S23(d1,d2,l1,l2)*H13(l1,l2,d1,d2)",Z26,S23,H13,scale=(0.25d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=H3(e1be2b|m1bm2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=H3(e1,e2,m1,m2)",Z27,H3,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=H2(e1be2b|d1bd2b)*S23(d1bd2b|m1bm2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=H2(e1,e2,d1,d2)*S23(d1,d2,m1,m2)",Z27,H2,S23,scale=(0.5d0,0.0d0))
    ierr=construct_tensor_t(Z33_1, H1%datatype, (/nocc,nocc,nocc,nocc/), init_val=ZERO)
!   original expression: Z33_1(l1bl2b|m1bm2b)+=S23(d1bd2b|m1bm2b)*H13(l1bl2b|d1bd2b)
    ierr=t_einsum("Z33_1(l1,l2,m1,m2)+=S23(d1,d2,m1,m2)*H13(l1,l2,d1,d2)",Z33_1,S23,H13,scale=(0.5d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|l1bl2b)*Z33_1(l1bl2b|m1bm2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,l1,l2)*Z33_1(l1,l2,m1,m2)",Z27,S23,Z33_1,scale=(0.5d0,0.0d0))
    ierr=destroy_tensor_t(Z33_1)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|l1bl2b)*H14(l1bl2b|m1bm2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,l1,l2)*H14(l1,l2,m1,m2)",Z27,S23,H14,scale=(0.5d0,0.0d0))
    ierr=construct_tensor_t(Z30_1, H1%datatype, (/nvir,nvir/), init_val=ZERO)
!   original expression: Z30_1(e1b|d1b)+=S23(e1b,d2b|l1bl2b)*H13(l1bl2b|d1b,d2b)*-1.
    ierr=t_einsum("Z30_1(e1,d1)+=S23(e1,d2,l1,l2)*H13(l1,l2,d1,d2)",Z30_1,S23,H13,scale=(-0.5d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m1bm2b)*Z30_1(e2b|d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m1,m2)*Z30_1(e2,d1)",Z27,S23,Z30_1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m1bm2b)*Z30_1(e1b|d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m1,m2)*Z30_1(e1,d1)",Z27,S23,Z30_1,scale=(-1.0d0,0.0d0))
    ierr=destroy_tensor_t(Z30_1)
    ierr=construct_tensor_t(Z30_1, H1%datatype, (/nvir,nvir/), init_val=ZERO)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m1bm2b)*Z30_1(e2b|d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m1,m2)*Z30_1(e2,d1)",Z27,S23,Z30_1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m1bm2b)*Z30_1(e1b|d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m1,m2)*Z30_1(e1,d1)",Z27,S23,Z30_1,scale=(-1.0d0,0.0d0))
    ierr=destroy_tensor_t(Z30_1)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m1bm2b)*H1(e2b|d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m1,m2)*H1(e2,d1)",Z27,S23,H1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m1bm2b)*H1(e1b|d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m1,m2)*H1(e1,d1)",Z27,S23,H1,scale=(-1.0d0,0.0d0))
    ierr=construct_tensor_t(Z32_1, H1%datatype, (/nocc,nocc/), init_val=ZERO)
!   original expression: Z32_1(l1b|m1b)+=S23(d1bd2b|m1b,l2b)*H13(l1b,l2b|d1bd2b)*-1.
    ierr=t_einsum("Z32_1(l1,m1)+=S23(d1,d2,m1,l2)*H13(l1,l2,d1,d2)",Z32_1,S23,H13,scale=(-0.5d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m1b,l1b)*Z32_1(l1b|m2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m1,l1)*Z32_1(l1,m2)",Z27,S23,Z32_1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m2b,l1b)*Z32_1(l1b|m1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m2,l1)*Z32_1(l1,m1)",Z27,S23,Z32_1,scale=(-1.0d0,0.0d0))
    ierr=destroy_tensor_t(Z32_1)
    ierr=construct_tensor_t(Z32_1, H1%datatype, (/nocc,nocc/), init_val=ZERO)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m1b,l1b)*Z32_1(l1b|m2b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m1,l1)*Z32_1(l1,m2)",Z27,S23,Z32_1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m2b,l1b)*Z32_1(l1b|m1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m2,l1)*Z32_1(l1,m1)",Z27,S23,Z32_1,scale=(-1.0d0,0.0d0))
    ierr=destroy_tensor_t(Z32_1)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m1b,l1b)*H9(l1b|m2b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m1,l1)*H9(l1,m2)",Z27,S23,H9,scale=(-1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1be2b|m2b,l1b)*H9(l1b|m1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,e2,m2,l1)*H9(l1,m1)",Z27,S23,H9,scale=(1.0d0,0.0d0))
    ierr=construct_tensor_t(Z37_1, H1%datatype, (/nvir,nocc,nocc,nvir/), init_val=ZERO)
!   original expression: Z37_1(e1b,l1b|m1b,d1b)+=S23(e1b,d2b|m1b,l2b)*H13(l1b,l2b|d1b,d2b)*0.5
    ierr=t_einsum("Z37_1(e1,l1,m1,d1)+=S23(e1,d2,m1,l2)*H13(l1,l2,d1,d2)",Z37_1,S23,H13,scale=(0.5d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m1b,l1b)*Z37_1(e2b,l1b|m2b,d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m1,l1)*Z37_1(e2,l1,m2,d1)",Z27,S23,Z37_1,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m1b,l1b)*Z37_1(e1b,l1b|m2b,d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m1,l1)*Z37_1(e1,l1,m2,d1)",Z27,S23,Z37_1,scale=(-1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m2b,l1b)*Z37_1(e2b,l1b|m1b,d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m2,l1)*Z37_1(e2,l1,m1,d1)",Z27,S23,Z37_1,scale=(-1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m2b,l1b)*Z37_1(e1b,l1b|m1b,d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m2,l1)*Z37_1(e1,l1,m1,d1)",Z27,S23,Z37_1,scale=(1.0d0,0.0d0))
    ierr=destroy_tensor_t(Z37_1)
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m1b,l1b)*H10(e2b,l1b|m2b,d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m1,l1)*H10(e2,l1,m2,d1)",Z27,S23,H10,scale=(1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m1b,l1b)*H10(e1b,l1b|m2b,d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m1,l1)*H10(e1,l1,m2,d1)",Z27,S23,H10,scale=(-1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e1b,d1b|m2b,l1b)*H10(e2b,l1b|m1b,d1b)*-1.
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e1,d1,m2,l1)*H10(e2,l1,m1,d1)",Z27,S23,H10,scale=(-1.0d0,0.0d0))
!   original expression: Z27(e1be2b|m1bm2b)+=S23(e2b,d1b|m2b,l1b)*H10(e1b,l1b|m1b,d1b)
    ierr=t_einsum("Z27(e1,e2,m1,m2)+=S23(e2,d1,m2,l1)*H10(e1,l1,m1,d1)",Z27,S23,H10,scale=(1.0d0,0.0d0))


end subroutine eti_CCD_codegen_call
