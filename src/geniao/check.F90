module check

contains

 subroutine check_normalization(matrix,row,col)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer, intent(in)     :: row,col
       complex(8), intent(in)  :: matrix(row,col)
       integer                 :: ierr,i,j
       complex(8)              :: sum_column

       ierr=0
       do i=1,col
         sum_column = ZERO
         do j=1,row
           sum_column = sum_column + conjg(matrix(j,i))*matrix(j,i)
         end do
         if (REAL(sum_column) .lt. 0.999D0 .or. AIMAG(sum_column) .gt. 1.0D-3) then
            ierr=1
         endif
       end do

       if (ierr.eq.1) then
        write(*,*)"**Warning** The matrix is not normalized."
       else if (ierr.eq.0) then
        write(*,*)"The matrix is normalized."
       endif
 end subroutine check_normalization

 subroutine check_SymOrth(C,S,row,col)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer, intent(in)     :: row,col
       complex(8), intent(in)  :: C(row,col), S(row,row)
       complex(8)              :: check(col,col),Id(col,col)
       integer                 :: i,j

       Id = ZERO
       forall (i=1:col) Id(i,i)=ONE
       check = matmul(conjg(transpose(C)),matmul(S,C))
       write(*,*)"Deviation from orthogonality      :", &
                 norm2([norm2(real(check-Id)),norm2(aimag(check-Id))])
 end subroutine check_SymOrth

 subroutine deviation_nelec(rdm,n,nelec)
       complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
       integer,    intent(in)  :: nelec,n
       complex(8), intent(in)  :: rdm(n,n)
       complex(8)              :: sum_diag
     
       sum_diag = ZERO
       do i=1,n
         sum_diag = sum_diag + rdm(i,i)
       enddo
       write(*,*)"Deviation in number of electrons  :", real(sum_diag - nelec)
 end subroutine deviation_nelec

 subroutine test_IAO_and_IBO(IAO,Cocc_IAO,Cocc_IBO_moB1,S11,cspinor_original,nmo1,nmo2,nocc)
        use exacorr_datatypes
        use exacorr_mo
        use energy

        complex(8), parameter   :: ZERO=(0.D0,0.D0), ONE=(1.D0,0.D0)
        integer, intent(in)     :: nmo1,nmo2,nocc
        complex(8), intent(in)  :: IAO(nmo1,nmo2),Cocc_IAO(nocc,nmo2), Cocc_IBO_moB1(nmo1,nocc), S11(nmo1,nmo1)
        type(cmo), intent(in)   :: cspinor_original
        complex(8), allocatable :: HF_Fock(:,:), HF_RDM(:,:), IAO_CoreH(:,:), IAO_Fock(:,:), IAO_RDM(:,:)
        complex(8), allocatable :: IBO_RDM(:,:)
        type(one_el_t)          :: int_1el
        real(8)                 :: E_HF, E_HF_IAO, E_HF_IBO

        ! HF energy = 0.5* Tr[(h + F)*gamma] + nuclear_repulsion
        ! (If it was only the Fock matrix, it would be the some of orbital energies, which is obviously only correct for a non-interacting case.)

        ! Therefore, we need to read the one-electron integrals, in the MRCONEE file.
        ! Read the integrals to get the H core matrix (int_1el%h_core)
        call read_from_mrconee (int_1el, ierr)
        ! Fock matrix:
        allocate(HF_Fock(nmo1,nmo1))
        HF_Fock = ZERO
        forall (i=1:nmo1) HF_Fock(i,i) = cspinor_original%energy(i)

        ! 1RDM:
        allocate(HF_RDM(nmo1,nmo1)) ! In the MO basis, should be 0 everywhere except 1 on part of the diagonal (only for nocc)
        HF_RDM = ZERO
        forall (i=1:nocc) HF_RDM(i,i) = ONE

        allocate(IAO_CoreH(nmo2,nmo2))
        IAO_CoreH = matmul(conjg(transpose(IAO)),matmul(int_1el%h_core,IAO))

        ! Fock matrix in the IAO basis
        allocate(IAO_Fock(nmo2,nmo2))
        IAO_Fock = matmul(conjg(transpose(IAO)),matmul(HF_Fock,IAO))

        ! 1RDM in the IAO basis
        allocate(IAO_RDM(nmo2,nmo2))
        IAO_RDM = matmul(conjg(transpose(IAO)),matmul(S11,matmul(HF_RDM,matmul(S11,IAO))))

        ! 1RDM in the IBO basis
        allocate(IBO_RDM(nmo1,nmo1))
        IBO_RDM = matmul(Cocc_IBO_moB1,conjg(transpose(Cocc_IBO_moB1)))

        call check_normalization(IAO,nmo1,nmo2)
        call check_normalization(Cocc_IAO,nmo2,nocc)
        call check_normalization(Cocc_IBO_moB1,nmo1,nocc)
        write(*,*) "Deviation of IBO_RDM from HF_RDM  :", norm2([norm2(real(IBO_RDM-HF_RDM)),norm2(aimag(IBO_RDM-HF_RDM))])
        call check_SymOrth(IAO,S11,nmo1,nmo2)
        call check_SymOrth(Cocc_IBO_moB1,S11,nmo1,nocc)
        call deviation_nelec(IAO_RDM,nmo2,nocc)
        call deviation_nelec(IBO_RDM,nmo1,nocc)
        ! Energy:
        if (int_1el%n_spinor .ne. nmo1) then
          print*,
          print*,'**WARNING** number of spinor when reading MRCONEE is not the same as nmo1 of MO basis B1.'
          print*,"            Cannot compute the HF energy."
        else
          call HF_energy(E_HF,int_1el%h_core,HF_Fock,HF_RDM,nmo1,int_1el%e_core)
          call HF_energy(E_HF_IAO,IAO_CoreH,IAO_Fock,IAO_RDM,nmo2,int_1el%e_core)
          call HF_energy(E_HF_IBO,int_1el%h_core,HF_Fock,IBO_RDM,nmo1,int_1el%e_core)
          print*,"HF energy (moB1)   : ",E_HF
          print*,"HF energy (IAO)    : ",E_HF_IAO
          print*,"HF energy (IBO/MO) : ",E_HF_IBO
          print*,"Deviation from the HF energy in the original basis and in the IAO basis:",REAL(E_HF_IAO - E_HF)
          print*,"Deviation from the HF energy in the original basis and in the IBO basis:",REAL(E_HF_IBO - E_HF)
        endif

        deallocate(IAO_Fock,IAO_RDM,HF_RDM,HF_Fock,IBO_RDM)
     55 FORMAT( 100(:,1X,'(',G10.2,',',G10.2,')') ) ! Useful format to write complex matrices.
 end subroutine test_IAO_and_IBO
end module
