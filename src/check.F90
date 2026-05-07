module check

contains

 subroutine check_normalization(M)
       implicit none
       real(8), intent(in)     :: M(:,:,:)
       integer                 :: ierr,i,j,k,rcq
       real(8)                 :: sum_column, absMji

       rcq = size(M,3)
       ierr=0
       do i = 1, size(M,2)
         sum_column = 0.D0
         do j = 1, size(M,1)
           if (rcq .eq. 1) then
              absMji = M(j,i,1)**2
           else if (rcq .eq. 2) then
              absMji = M(j,i,1)**2 + M(j,i,2)**2
           else if (rcq .eq. 4) then
              absMji = M(j,i,1)**2 + M(j,i,2)**2 + M(j,i,3)**2 + M(j,i,4)**2
           end if
           sum_column = sum_column + absMji
         end do
         if (sum_column .lt. 0.999D0 .or. sum_column .gt. 1.001) then
            ierr=1
         endif
       end do

       if (ierr.eq.1) then
        write(*,*)"**WARNING** The matrix is not normalized."
       else if (ierr.eq.0) then
        write(*,*)"The matrix is normalized."
       endif
 end subroutine check_normalization

 subroutine check_hermiticity(M)

       use linear_algebra

       implicit none
       real(8), intent(in)     :: M(:,:,:)
       real(8), allocatable    :: transconjgM(:,:,:)
       integer                 :: i,j,rcq
       real(8)                 :: norm

       rcq = size(M,3)
       if (size(M,1) .ne. size(M,2)) stop 'Hermiticity check: matrix is not square'

       allocate(transconjgM(size(M,2),size(M,1),size(M,3)))
       call transpose_conjg(M,transconjgM)
       call diff_norm2_rcq(M,transconjgM,norm)
       deallocate(transconjgM)
       print*,"Deviation from hermicity: ",norm

 end subroutine check_hermiticity


 subroutine check_SymOrth(C,S)

       use linear_algebra

       implicit none
       real(8), intent(in)  :: C(:,:,:), S(:,:,:)
       real(8), allocatable :: check(:,:,:),Id(:,:,:), transconjgC(:,:,:)
       integer              :: i,j
       real(8)              :: norm

       if (size(C,3) .ne. size(S,3)) then
           write(*,*)"SymOrth: two different matrix types (real, complex or quaternion)."
           stop
       endif

       allocate(Id(size(C,2),size(C,2),size(S,3)))
       allocate(check(size(C,2),size(C,2),size(S,3)))
       allocate(transconjgC(size(C,2),size(C,1),size(S,3)))

       Id = 0.D0
       forall (i=1:size(C,2)) Id(i,i,1)=1.D0
       
       call transpose_conjg(C,transconjgC)
       call matrix_mul(transconjgC,S,C,check)
       call diff_norm2_rcq(check,Id,norm)

       write(*,*)"Deviation from orthogonality      :", norm
       deallocate(Id,check,transconjgC)

 end subroutine check_SymOrth

 subroutine deviation_nelec(rdm,nelec)
       implicit none
       integer              :: rcq, i, k
       integer, intent(in)  :: nelec
       real(8), intent(in)  :: rdm(:,:,:)
       real(8)              :: sum_diag
     
       rcq = size(rdm,3)
       if (size(rdm,1) .ne. size(rdm,2)) then
          print*,"deviation_nelec : rdm is not square."
          stop
       end if
       sum_diag = 0.D0
       do i=1,size(rdm,1)
          sum_diag = sum_diag + rdm(i,i,1)
       end do

       write(*,*)"Deviation in number of electrons  :",sum_diag-nelec
          
 end subroutine deviation_nelec

 subroutine deviation_idempotency(rdm)

       use linear_algebra

       integer              :: rcq
       real(8), intent(in)  :: rdm(:,:,:)
       real(8), allocatable :: rdm_squared(:,:,:)
       real(8)              :: norm

       rcq = size(rdm,3)
       if (size(rdm,1) .ne. size(rdm,2)) then
          print*,"deviation_idempotency : rdm is not square."
          stop
       end if

       allocate(rdm_squared(size(rdm,1),size(rdm,2),size(rdm,3)))
       call matrix_mul(rdm,rdm,rdm_squared)
       call diff_norm2_rcq(rdm,rdm_squared,norm)
       write(*,*)"Deviation in idempotency  :", norm
       deallocate(rdm_squared)

 end subroutine deviation_idempotency

 subroutine purity_indicator(rdm,nelec)

       use linear_algebra

       integer              :: rcq,i,nelec
       real(8), intent(in)  :: rdm(:,:,:)
       real(8), allocatable :: rdm_squared(:,:,:)
       real(8)              :: purity_indic,norm

       rcq = size(rdm,3)
       if (size(rdm,1) .ne. size(rdm,2)) then
          print*,"deviation_idempotence : rdm is not square."
          stop
       end if

       allocate(rdm_squared(size(rdm,1),size(rdm,2),size(rdm,3)))
       call matrix_mul(rdm,rdm,rdm_squared)
       call diff_norm2_rcq(rdm,rdm_squared,norm)
       write(*,*)"Deviation in idempotency  :", norm
       purity_indic = norm/dble(nelec)
       write(*,*)"Purity indicator : ", purity_indic
       deallocate(rdm_squared)
 
 end subroutine purity_indicator

 subroutine test_IAO_and_IBO(IAO,Cocc_IAO,Cocc_IBO_moB1,spinor_original,ncore,restricted)
        use rose_ao
        use rose_mo
        use rose_utils
        use linear_algebra

        implicit none
        real(8), intent(in)        :: IAO(:,:,:),Cocc_IAO(:,:,:), Cocc_IBO_moB1(:,:,:)
        real(8), allocatable       :: S11(:,:,:)
        logical, intent(in)        :: restricted
        type(mo_basis), intent(in) :: spinor_original
        integer, intent(in)        :: ncore
        logical                    :: tobe
        integer                    :: rcq, nmo1, nmo2, nocc, nelec, i, j, ierr, luxyz, nelec_val
        real(8), allocatable       :: HF_Fock(:,:,:), HF_RDM(:,:,:), IAO_Hcore(:,:,:), IAO_Fock(:,:,:), IAO_RDM(:,:,:)
        real(8), allocatable       :: IBO_RDM(:,:,:), core_RDM(:,:,:),Hcore(:,:,:), transconjgIAO(:,:,:), transconjgIBO(:,:,:)
        real(8)                    :: E_core, E_HF, E_HF_IAO, E_HF_IBO, norm

        rcq  = size(IAO,3)
        nmo1 = size(IAO,1) + ncore ! The IAO and NBO unitaries do not contain the core
        nmo2 = size(IAO,2)
        nocc = size(Cocc_IAO,2)

        write(*,*)
        write(*,*) "####################################"
        write(*,*) "#         Starting TESTS           #"
        write(*,*) "####################################"
        write(*,*)

        ! S11 (is unity as should be tested long before this stage):
        allocate(S11(size(IAO,1),size(IAO,1),rcq))
        S11 = 0.D0
        forall (i=1:size(IAO,1)) S11(i,i,1) = 1.D0

        ! Full Fock matrix:
        allocate(HF_Fock(nmo1,nmo1,rcq))
        HF_Fock = 0.D0
        forall (i=1:nmo1) HF_Fock(i,i,1) = spinor_original%energy(i)

        ! Full (core+valence) reference 1RDM:
        allocate(HF_RDM(nmo1,nmo1,rcq)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal (only for nocc)
        HF_RDM = 0.D0
        forall (i=1:ncore+nocc) HF_RDM(i,i,1) = 1.D0
        if (restricted) HF_RDM = 2*HF_RDM

        ! Core-only reference 1RDM:
        allocate(core_RDM(nmo1,nmo1,rcq)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal (only for core)
        core_RDM = 0.D0
        forall (i=1:ncore) core_RDM(i,i,1) = 1.D0
        if (restricted) core_RDM = 2*core_RDM

        ! Valence 1RDM in the IAO basis
        allocate(transconjgIAO(nmo2,nmo1-ncore,rcq))
        allocate(IAO_RDM(nmo2,nmo2,rcq))
        call transpose_conjg(IAO,transconjgIAO)
        call matrix_mul(transconjgIAO,HF_RDM(ncore+1:nmo1,ncore+1:nmo1,:),IAO,IAO_RDM)
        deallocate(transconjgIAO)

        ! Valence 1RDM in the IBO basis
        allocate(IBO_RDM(nmo1,nmo1,rcq))
        IBO_RDM = 0.D0
        forall (i=1:ncore) IBO_RDM(i,i,1) = 1.D0
        allocate(transconjgIBO(nocc,nmo1-ncore,rcq))
        call transpose_conjg(Cocc_IBO_moB1,transconjgIBO)
        call matrix_mul(Cocc_IBO_moB1,transconjgIBO,IBO_RDM(ncore+1:nmo1,ncore+1:nmo1,:))
        deallocate(transconjgIBO)

        nelec = nocc + ncore
        nelec_val = nocc
        if (restricted) then
           IBO_RDM(:,:,1) = 2*IBO_RDM(:,:,1)
           nelec   = 2*(nocc+ncore)
           nelec_val = 2*nocc
        end if

        call check_normalization(IAO)
        call check_normalization(Cocc_IAO)
        call check_normalization(Cocc_IBO_moB1)
        call diff_norm2_rcq(IBO_RDM,HF_RDM,norm)
        write(*,*) "Deviation of IBO_RDM from HF_RDM  :", norm
        call check_SymOrth(IAO,S11)
        call check_SymOrth(Cocc_IBO_moB1,S11)
        call deviation_nelec(IAO_RDM,nelec_val)
        call deviation_nelec(IBO_RDM,nelec)

     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.
 end subroutine test_IAO_and_IBO

end module
