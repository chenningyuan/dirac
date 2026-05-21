module eti_drivers_gs_automatic
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_cfg, only: exacc_input
        use exacorr_utils
        use intermediates_cc
        use intermediates_lambda
        use eti_common_routines
        use eti_diis
        use mfsolver_dirac_input
        use exacorr_tensor_interface
        use exacorr_printout
        use exacorr_scaling_parameters
        use, intrinsic:: ISO_C_BINDING

        implicit none

        type(tensor_t) :: one_tensor

        private

        public solve_t_equations_automatic
        public solve_lambda_equations_automatic

       contains
!-------------------------------------------------------------------------------------------

! aspg
        subroutine solve_t_equations_automatic(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision,&
                                       ccsd_energy,CCD,iterative_level,iterative_scheme,auto_scheme,&
                                       ncycles,print_level,benchmarking_skip_diis)

!        Routine for solving the t amplitude equations

         logical, intent(in), optional :: benchmarking_skip_diis
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_intg_tens_t), intent(inout) :: int_t
!        intermediate tensors (same notation as in relccsd)
!        temporarily using tensors with permutations required by the automatically generated code
!        the next update to the code generator is expected to handle index permutations/conjugation 
         type(tensor_t) :: vvoo_aux_tensor
         type(tensor_t) :: voov_aux_tensor 
         type(tensor_t) :: vooo_aux_tensor 
         type(tensor_t) :: vvov_aux_tensor 
         type(tensor_t) :: hvo_aux_tensor
!        solution tensors
         type(tensor_t) :: s1_tensor
         type(tensor_t) :: s2_tensor
         type(tensor_t) :: s3_tensor
         type(tensor_t) :: s4_tensor
         type(tensor_t) :: s5_tensor
!        scalars (need to be defined as tensor types)
         type(tensor_t) :: result_tensor
         integer(INTD)      :: result_dims(1)
         !Tensor dimensions
         integer(C_INT)     :: one_dims(1)
         !Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims
         integer(INTD), dimension(2) :: ov_dims
         integer(INTD), dimension(2) :: vo_dims
         integer(INTD), dimension(4) :: vvoo_dims
         integer(INTD), dimension(4) :: ovoo_dims
         integer(INTD), dimension(4) :: oovv_dims
         integer(INTD), dimension(4) :: voov_dims
         integer(INTD), dimension(4) :: vvov_dims
         integer(INTD), dimension(4) :: vooo_dims
         integer(INTD), dimension(6) :: vvvooo_dims
         integer(INTD), dimension(8) :: vvvvoooo_dims
         integer(INTD), dimension(10):: vvvvvooooo_dims
         !Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: iccsd
         integer, intent(in) :: ncycles
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         real(8), intent(out) :: ccsd_energy
         logical, intent(in)  :: CCD
         integer, intent(in)  :: iterative_level
         integer, intent(in)  :: iterative_scheme
         integer, intent(in)  :: auto_scheme
         integer, intent(in)  :: print_level
         integer (INTD)       :: ierr
         logical :: benchmarking_skip_diis_local

         if (present(benchmarking_skip_diis)) then
            benchmarking_skip_diis_local = benchmarking_skip_diis
         else
            benchmarking_skip_diis_local = .false.
         end if
         print *,''
         print *,' implementation with automatically generated residual equations'
         print *,''

         result_dims = 1
         ierr=construct_tensor_t(one_tensor,C8,result_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)

         if (.not.CCD) then
            oo_dims = (/nocc, nocc/)
            vo_dims = (/nvir, nocc/)
            ov_dims = (/nocc, nvir/)
            vo_dims = (/nvir, nocc/)

            ierr=construct_tensor_t(hvo_aux_tensor,C8,vo_dims,init_val=ZERO)
            ierr=t_einsum("h(a,i)+=h+(i,a)",hvo_aux_tensor,int_t%fov)
            ierr=construct_tensor_t(s1_tensor,C8,vo_dims,init_val=ZERO)

            vvov_dims = (/nvir, nvir, nocc, nvir/)
            vooo_dims = (/nvir, nocc, nocc, nocc/)

            ierr=construct_tensor_t(vooo_aux_tensor,C8,vooo_dims,init_val=ZERO)
            ierr=construct_tensor_t(vvov_aux_tensor,C8,vvov_dims,init_val=ZERO)
            ierr=initialize_tensor_t(vvov_aux_tensor)
            ierr=initialize_tensor_t(vooo_aux_tensor)
            ierr=t_einsum("T(a,k,i,j)+=V+(i,j,k,a)*K()", vooo_aux_tensor, int_t%ooov,one_tensor,scale=MINUS_ONE)
            ierr=t_einsum("T(a,b,i,c)+=V+(c,i,a,b)*K()", vvov_aux_tensor, int_t%vovv,one_tensor,scale=MINUS_ONE)
         end if
         vvoo_dims = (/nvir, nvir, nocc, nocc/)
         oovv_dims = (/nocc, nocc, nvir, nvir/)
         voov_dims = (/nvir, nocc, nocc, nvir/)

         ierr=construct_tensor_t(s2_tensor,C8,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(vvoo_aux_tensor,C8,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(voov_aux_tensor,C8,voov_dims,init_val=ZERO)

         ierr=initialize_tensor_t(vvoo_aux_tensor)
         ierr=initialize_tensor_t(voov_aux_tensor)
         ierr=t_einsum("T(a,b,i,j)+=V+(i,j,a,b)*K()",vvoo_aux_tensor, int_t%oovv,one_tensor)
         ierr=t_einsum("T(a,i,j,b)+=V(a,i,b,j)*K()", voov_aux_tensor, int_t%vovo,one_tensor,scale=MINUS_ONE)

         if (iterative_level.ge.3) then
            vvvooo_dims = (/nvir, nvir, nvir, nocc, nocc, nocc/)
            ierr=construct_tensor_t(s3_tensor,C8,vvvooo_dims,init_val=ZERO)
         end if

         if (iterative_level.ge.4) then
            vvvvoooo_dims = (/nvir, nvir, nvir, nvir, nocc, nocc, nocc, nocc/)
            ierr=construct_tensor_t(s4_tensor,C8,vvvvoooo_dims,init_val=ZERO)
         end if

         if (iterative_level.ge.5) then
            vvvvvooooo_dims = (/nvir,nvir, nvir, nvir, nvir, nocc,nocc, nocc, nocc, nocc/)
            ierr=construct_tensor_t(s5_tensor,C8,vvvvvooooo_dims,init_val=ZERO)
         end if

         do iccsd = 1, ncycles

          if (print_level.gt.2) then
            if (.not.CCD) then 
               write(*,*) "t1    = ",tensor_norm2(comm_t%t1)
               if (iterative_level.ge.3) write(*,*) "t3    = ",tensor_norm2(comm_t%t3)
               if (iterative_level.ge.4) write(*,*) "t4    = ",tensor_norm2(comm_t%t4)
               if (iterative_level.ge.5) write(*,*) "t5    = ",tensor_norm2(comm_t%t5)
            end if
            write(*,*) "t2    = ",tensor_norm2(comm_t%t2)
          end if

          if (print_level.gt.2) then
            if (.not.CCD) then
               write(*,*) "s1    = ",tensor_norm2(s1_tensor)
               if (iterative_level.ge.3) write(*,*) "s3    = ",tensor_norm2(s3_tensor)
               if (iterative_level.ge.4) write(*,*) "s4    = ",tensor_norm2(s4_tensor)
               if (iterative_level.ge.5) write(*,*) "s5    = ",tensor_norm2(s5_tensor)
            end if
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)
          end if

          if (ierr.ne.0) then
             write(*,*) " ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

          ierr=initialize_tensor_t(result_tensor)
          if (.not.CCD) then
             ierr=initialize_tensor_t(s1_tensor)
             if (iterative_level.ge.3) ierr=initialize_tensor_t(s3_tensor)
             if (iterative_level.ge.4) ierr=initialize_tensor_t(s4_tensor)
             if (iterative_level.ge.5) ierr=initialize_tensor_t(s5_tensor)
          end if
          ierr=initialize_tensor_t(s2_tensor)

          write(*,'(A,I5)') " auto_scheme :      ", auto_scheme
          if ((iterative_level.eq.2).and.(iterative_scheme.eq.21)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
              call  eti_CCSD_codegen(nocc,nvir,          &
                                     int_t%oovv,         &   ! H24 OOVV 
                                     comm_t%t1,          &   ! S43 VO
                                     comm_t%fov,         &   ! H15 OV
                                     comm_t%t2,          &   ! S44 VVOO
                                     hvo_aux_tensor,     &   ! H2  VO
                                     comm_t%fvv,         &   ! H1  VV
                                     int_t%vovv,         &   ! H17 VOVV
                                     voov_aux_tensor,    &   ! H18 VOOV
                                     int_t%ooov,         &   ! H25 OOOV
                                     comm_t%foo,         &   ! H16 OO
                                     vvoo_aux_tensor,    &   ! H5  VVOO
                                     int_t%vvvv,         &   ! H3  VVVV
                                     vooo_aux_tensor,    &   ! H19 VOOO
                                     vvov_aux_tensor,    &   ! H4  VVOV
                                     int_t%oooo,         &   ! H26 OOOO
                                     result_tensor,      &   ! Z48 UnitTensor
                                     s1_tensor,          &   ! Z49 VO
                                     s2_tensor)              ! Z50 VVOO
            else if (auto_scheme.eq.2) then  
            ! tenpi
                 call   eti_tenpi_ccsd(nocc,nvir,   &
                                      comm_t%foo,         &    !       OO
                                      comm_t%fov,         &    !       OV
                                      hvo_aux_tensor,     &    !       f_VO
                                      comm_t%fvv,         &    !       VV
                                      comm_t%t1,          &    !       t1_VO
                                      comm_t%t2,          &    !       t2_VVOO
                                      int_t%oooo,         &    !       OOOO
                                      int_t%ooov,         &    !       OOOV
                                      int_t%oovv,         &    !       OOVV
                                      vooo_aux_tensor,    &    !       VOOO
                                      voov_aux_tensor,    &    !       VOOV
                                      int_t%vovv,         &    !       VOVV
                                      vvoo_aux_tensor,    &    !       v_VVOO
                                      vvov_aux_tensor,    &    !       VVOV
                                      int_t%vvvv,         &    !       VVVV
                                      result_tensor,      &    !       UnitTensor
                                      s1_tensor,          &    !       r_VO
                                      s2_tensor)               !       r_VVOO
            else if (auto_scheme.eq.3) then  
            ! tenpi
               call tapp_tenpi_ccsd(nocc,nvir,   &
                                    comm_t%foo,         &    !       OO
                                    comm_t%fov,         &    !       OV
                                    hvo_aux_tensor,     &    !       f_VO
                                    comm_t%fvv,         &    !       VV
                                    comm_t%t1,          &    !       t1_VO
                                    comm_t%t2,          &    !       t2_VVOO
                                    int_t%oooo,         &    !       OOOO
                                    int_t%ooov,         &    !       OOOV
                                    int_t%oovv,         &    !       OOVV
                                    vooo_aux_tensor,    &    !       VOOO
                                    voov_aux_tensor,    &    !       VOOV
                                    int_t%vovv,         &    !       VOVV
                                    vvoo_aux_tensor,    &    !       v_VVOO
                                    vvov_aux_tensor,    &    !       VVOV
                                    int_t%vvvv,         &    !       VVVV
                                    result_tensor,      &    !       UnitTensor
                                    s1_tensor,          &    !       r_VO
                                    s2_tensor)               !       r_VVOO
            end if
          else if ((iterative_level.eq.2).and.(iterative_scheme.eq.20)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
              call   eti_CCD_codegen_call(nocc,nvir,      &
                                          comm_t%t2,      &   ! S23 
                                          int_t%oovv,     &   ! H13   OOVV
                                          vvoo_aux_tensor,&   ! H3    VVOO
                                          int_t%vvvv,     &   ! H2    VVVV
                                          int_t%oooo,     &   ! H14   OOOO 
                                          comm_t%fvv,     &   ! H1    VV
                                          comm_t%foo,     &   ! H9    OO
                                          voov_aux_tensor,&   ! H10   VOOV
                                          result_tensor,  &   ! Z26   unit tensor
                                          s2_tensor)          ! Z27   VVOO
            else if (auto_scheme.eq.2) then  
            ! tenpi
                 call   eti_tenpi_ccd(nocc,nvir,            &
                                         comm_t%foo,         &    !       OO
                                         comm_t%fvv,         &    !       VV
                                         comm_t%t2,          &    !       VVOO
                                         int_t%oooo,         &    !       OOOO
                                         int_t%oovv,         &    !       OOVV
                                         voov_aux_tensor,    &    !       VOOV
                                         vvoo_aux_tensor,    &    !       v_VVOO
                                         int_t%vvvv,         &    !       VVVV
                                         result_tensor,      &    !       UnitTensor
                                         s2_tensor)               !       r_VVOO
            else if (auto_scheme.eq.3) then  
            ! tenpi
               call tapp_tenpi_ccd(nocc,nvir,            &
                                       comm_t%foo,         &    !       OO
                                       comm_t%fvv,         &    !       VV
                                       comm_t%t2,          &    !       VVOO
                                       int_t%oooo,         &    !       OOOO
                                       int_t%oovv,         &    !       OOVV
                                       voov_aux_tensor,    &    !       VOOV
                                       vvoo_aux_tensor,    &    !       v_VVOO
                                       int_t%vvvv,         &    !       VVVV
                                       result_tensor,      &    !       UnitTensor
                                       s2_tensor)               !       r_VVOO
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.30)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDT_codegen_call(nocc,nvir,          &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor,  &  !  VO
                                             comm_t%fvv,      &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor, &  !  VVOO
                                             int_t%vvvv,      &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor, &  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor)         !  VVVOOO
            else if (auto_scheme.eq.2) then 
            ! tenpi
                 call   eti_tenpi_ccsdt(nocc,nvir,   &
                                            comm_t%foo,         &    ! F2     OO
                                            comm_t%fov,         &    ! F1     OV
                                            hvo_aux_tensor,     &    ! F4      VO
                                            comm_t%fvv,         &    ! F3      VV
                                            comm_t%t1,          &    ! T1      VO
                                            comm_t%t2,          &    ! T2      VVOO
                                            comm_t%t3,          &    ! T3      VVVOOO
                                            int_t%oooo,         &    ! V7      OOOO
                                            int_t%ooov,         &    ! V3      OOOV
                                            int_t%oovv,         &    ! V1      OOVV
                                            vooo_aux_tensor,    &    ! V8      VOOO
                                            voov_aux_tensor,    &    ! V4      VOOV
                                            int_t%vovv,         &    ! V2      VOVV
                                            vvoo_aux_tensor,    &    ! V9      v_VVOO
                                            vvov_aux_tensor,    &    ! V6      VVOV
                                            int_t%vvvv,         &    ! V5      VVVV
                                            result_tensor,      &    !       UnitTensor
                                            s1_tensor,          &    !       r_VO
                                            s2_tensor,          &    !       r_VVOO
                                            s3_tensor)               !       r_VVVOOO
            else if (auto_scheme.eq.3) then 
            ! tenpi
               call tapp_tenpi_ccsdt(nocc,nvir,   &
                                          comm_t%foo,         &    ! F2     OO
                                          comm_t%fov,         &    ! F1     OV
                                          hvo_aux_tensor,     &    ! F4      VO
                                          comm_t%fvv,         &    ! F3      VV
                                          comm_t%t1,          &    ! T1      VO
                                          comm_t%t2,          &    ! T2      VVOO
                                          comm_t%t3,          &    ! T3      VVVOOO
                                          int_t%oooo,         &    ! V7      OOOO
                                          int_t%ooov,         &    ! V3      OOOV
                                          int_t%oovv,         &    ! V1      OOVV
                                          vooo_aux_tensor,    &    ! V8      VOOO
                                          voov_aux_tensor,    &    ! V4      VOOV
                                          int_t%vovv,         &    ! V2      VOVV
                                          vvoo_aux_tensor,    &    ! V9      v_VVOO
                                          vvov_aux_tensor,    &    ! V6      VVOV
                                          int_t%vvvv,         &    ! V5      VVVV
                                          result_tensor,      &    !       UnitTensor
                                          s1_tensor,          &    !       r_VO
                                          s2_tensor,          &    !       r_VVOO
                                          s3_tensor)               !       r_VVVOOO
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.31)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDT_1a_codegen_call(nocc,nvir,          &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor,  &  !  VO
                                             comm_t%fvv,      &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor, &  !  VVOO
                                             int_t%vvvv,      &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor, &  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor)         !  VVVOOO
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.32)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDT_1b_codegen_call(nocc,nvir,          &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor,  &  !  VO
                                             comm_t%fvv,      &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor, &  !  VVOO
                                             int_t%vvvv,      &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor, &  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor)         !  VVVOOO
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.33)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDT_2_codegen_call(nocc,nvir,          &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor,  &  !  VO
                                             comm_t%fvv,      &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor, &  !  VVOO
                                             int_t%vvvv,      &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor, &  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor)         !  VVVOOO
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.34)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDT_3_codegen_call(nocc,nvir,          &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor,  &  !  VO
                                             comm_t%fvv,      &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor, &  !  VVOO
                                             int_t%vvvv,      &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor, &  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor)         !  VVVOOO
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.4).and.(iterative_scheme.eq.40)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDTQ_codegen_call(nocc,nvir,     &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor, &  !  VO
                                             comm_t%fvv,     &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor,&  !  VVOO
                                             int_t%vvvv,     &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor,&  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             comm_t%t4,      &  !  VVVVOOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor,      &  !  VVVOOO
                                             s4_tensor)         !  VVVVOOOO
            else if (auto_scheme.eq.2) then
            ! tenpi
                 call   eti_tenpi_ccsdtq(nocc,nvir,   &
                                            comm_t%foo,         &    ! F2     OO
                                            comm_t%fov,         &    ! F1     OV
                                            hvo_aux_tensor,     &    ! F4      VO
                                            comm_t%fvv,         &    ! F3      VV
                                            comm_t%t1,          &    ! T1      VO
                                            comm_t%t2,          &    ! T2      VVOO
                                            comm_t%t3,          &    ! T3      VVVOOO
                                            comm_t%t4,          &    ! T4      VVVOOO
                                            int_t%oooo,         &    ! V7      OOOO
                                            int_t%ooov,         &    ! V3      OOOV
                                            int_t%oovv,         &    ! V1      OOVV
                                            vooo_aux_tensor,    &    ! V5      VOOO
                                            voov_aux_tensor,    &    ! V2      VOOV
                                            int_t%vovv,         &    ! V4      VOVV
                                            vvoo_aux_tensor,    &    ! V9      v_VVOO
                                            vvov_aux_tensor,    &    ! V8      VVOV
                                            int_t%vvvv,         &    ! V6      VVVV
                                            result_tensor,      &    !       UnitTensor
                                            s1_tensor,          &    !       r_VO
                                            s2_tensor,          &    !       r_VVOO
                                            s3_tensor,          &    !       r_VVVOOO
                                            s4_tensor)               !       r_VVVVOOOO
            else if (auto_scheme.eq.3) then
            ! tenpi
               call tapp_tenpi_ccsdtq(nocc,nvir,   &
                                          comm_t%foo,         &    ! F2     OO
                                          comm_t%fov,         &    ! F1     OV
                                          hvo_aux_tensor,     &    ! F4      VO
                                          comm_t%fvv,         &    ! F3      VV
                                          comm_t%t1,          &    ! T1      VO
                                          comm_t%t2,          &    ! T2      VVOO
                                          comm_t%t3,          &    ! T3      VVVOOO
                                          comm_t%t4,          &    ! T4      VVVOOO
                                          int_t%oooo,         &    ! V7      OOOO
                                          int_t%ooov,         &    ! V3      OOOV
                                          int_t%oovv,         &    ! V1      OOVV
                                          vooo_aux_tensor,    &    ! V5      VOOO
                                          voov_aux_tensor,    &    ! V2      VOOV
                                          int_t%vovv,         &    ! V4      VOVV
                                          vvoo_aux_tensor,    &    ! V9      v_VVOO
                                          vvov_aux_tensor,    &    ! V8      VVOV
                                          int_t%vvvv,         &    ! V6      VVVV
                                          result_tensor,      &    !       UnitTensor
                                          s1_tensor,          &    !       r_VO
                                          s2_tensor,          &    !       r_VVOO
                                          s3_tensor,          &    !       r_VVVOOO
                                          s4_tensor)               !       r_VVVVOOOO
            end if
          else if ((iterative_level.eq.5).and.(iterative_scheme.eq.50)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
               call   eti_CCSDTQP_codegen_call(nocc,nvir,    &
                                             int_t%oovv,     &  !  OOVV
                                             comm_t%t1,      &  !  VO
                                             comm_t%fov,     &  !  OV
                                             comm_t%t2,      &  !  VVOO
                                             hvo_aux_tensor, &  !  VO
                                             comm_t%fvv,     &  !  VV
                                             int_t%vovv,     &  !  VOVV
                                             voov_aux_tensor,&  !  VOOV
                                             int_t%ooov,     &  !  OOOV
                                             comm_t%foo,     &  !  OO
                                             comm_t%t3,      &  !  VVVOOO
                                             vvoo_aux_tensor,&  !  VVOO
                                             int_t%vvvv,     &  !  VVVV
                                             vooo_aux_tensor,&  !  VOOO
                                             vvov_aux_tensor,&  !  VVOV
                                             int_t%oooo,     &  !  OOOO
                                             comm_t%t4,      &  !  VVVVOOOO
                                             comm_t%t5,      &  !  VVVVVOOOOO
                                             result_tensor,  &  !  UnitTensor
                                             s1_tensor,      &  !  VO
                                             s2_tensor,      &  !  VVOO
                                             s3_tensor,      &  !  VVVOOO
                                             s4_tensor,      &  !  VVVVOOOO
                                             s5_tensor)         !  VVVVVOOOOO 
            else if (auto_scheme.eq.2) then
            ! tenpi
               call quit ('Error: CCSDTQP not yet implemented in tenpi.')
            end if
          else
             write(*,*) "Implementation not available for this method"
          end if

          if (print_level.gt.2) then
            if (.not.CCD) then
               write(*,*) "s1    = ",tensor_norm2(s1_tensor)
               if (iterative_level.ge.3)  write(*,*) "s3    = ",tensor_norm2(s3_tensor)
               if (iterative_level.ge.4)  write(*,*) "s4    = ",tensor_norm2(s4_tensor)
               if (iterative_level.ge.5)  write(*,*) "s5    = ",tensor_norm2(s5_tensor)
            end if
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)
          end if

!         prepare for next iteration
          if (.not.CCD) then
             ierr=initialize_tensor_t(comm_t%t1)
             ierr=t_einsum("T1(a,i)+=S1(a,i)",comm_t%t1,s1_tensor)
             if (iterative_level.ge.3) then
                ierr=initialize_tensor_t(comm_t%t3)
                ierr=t_einsum("T3(a,b,c,i,j,k)+=S3(a,b,c,i,j,k)",comm_t%t3,s3_tensor)
             end if
             if (iterative_level.ge.4) then
                ierr=initialize_tensor_t(comm_t%t4)
                ierr=t_einsum("T4(a,b,c,d,i,j,k,l)+=S4(a,b,c,d,i,j,k,l)",comm_t%t4,s4_tensor)
             end if
             if (iterative_level.ge.5) then
                ierr=initialize_tensor_t(comm_t%t5)
                ierr=t_einsum("T5(a,b,c,d,e,i,j,k,l,m)+=S5(a,b,c,d,e,i,j,k,l,m)",comm_t%t5,s5_tensor)
             end if
          end if
          ierr=initialize_tensor_t(comm_t%t2)
          ierr=t_einsum("T2(a,b,i,j)+=S2(a,b,i,j)",comm_t%t2,s2_tensor)

          if (.not.CCD) then
            if (iterative_level.ge.5) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%t1,&
                              t2_tensor=comm_t%t2,full_t3_tensor=comm_t%t3,           &
                              full_t4_tensor=comm_t%t4,full_t5_tensor=comm_t%t5)
            else if (iterative_level.ge.4) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%t1,&
                              t2_tensor=comm_t%t2,full_t3_tensor=comm_t%t3,full_t4_tensor=comm_t%t4)
            else if (iterative_level.ge.3) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%t1,&
                              t2_tensor=comm_t%t2,full_t3_tensor=comm_t%t3)
            else
               call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)
            end if
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
          end if

! aspg: for testing/debugging purposes we can choose to skip diis (optional argument) 
      if (benchmarking_skip_diis_local) then
          if (.not.CCD) then
            if (iterative_level.ge.5) then
               call diis (t5_tensor=comm_t%t5,t4_tensor=comm_t%t4,t3_tensor=comm_t%t3,&
                         t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,t_convergence=t_convergence)
            else if (iterative_level.ge.4) then
               call diis (t4_tensor=comm_t%t4,t3_tensor=comm_t%t3,t2_tensor=comm_t%t2,&
                       & t1_tensor=comm_t%t1,t_convergence=t_convergence)
            else if (iterative_level.ge.3) then
               call diis (t3_tensor=comm_t%t3,t2_tensor=comm_t%t2,&
                       & t1_tensor=comm_t%t1,t_convergence=t_convergence)
            else 
               call diis (t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,t_convergence=t_convergence)
            end if
          else
            call diis (t2_tensor=comm_t%t2,t_convergence=t_convergence)
          end if
      end if

!!aspg, checking energy calculation/printout
          if (.true.) then
          if (.not.CCD) then
            ierr=initialize_tensor_t(result_tensor)
            ierr=t_einsum("E()+=T(a,i)*F(i,a)",result_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
            ccsd_energy=result_tens(1)
          else
            ccsd_energy = 0.D0
          end if

!         update tau
          call get_tau(comm_t,CCD)

          ierr=initialize_tensor_t(result_tensor)
          ierr=t_einsum("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER_R
          endif
!!aspg, checking energy calculation/printout

          call print_iteration(iccsd,t_convergence,ccsd_energy,print_level)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         write(*,*) ""


         if (.not.CCD) then
            ierr=destroy_tensor_t(s1_tensor)
            ierr=destroy_tensor_t(hvo_aux_tensor)
            ierr=destroy_tensor_t(vooo_aux_tensor)
            ierr=destroy_tensor_t(vvov_aux_tensor)
            if (iterative_level.ge.3) ierr=destroy_tensor_t(s3_tensor)
            if (iterative_level.ge.4) ierr=destroy_tensor_t(s4_tensor)
            if (iterative_level.ge.5) ierr=destroy_tensor_t(s5_tensor)
         end if
         ierr=destroy_tensor_t(vvoo_aux_tensor)
         ierr=destroy_tensor_t(voov_aux_tensor)
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(result_tensor)
         ierr=destroy_tensor_t(one_tensor)

        end subroutine
! aspg
        subroutine solve_lambda_equations_automatic(comm_t, int_t, nocc, nvir,eps_occ, eps_vir, t_target_precision, &
                                           CCD, iterative_level,iterative_scheme, auto_scheme, &
                                           ncycles, print_level)
        !Routine for solving the lambda equations

         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron
!        integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_intg_tens_t), intent(inout) :: int_t
         !type(talsh_lambda_tens_t),intent(inout) :: lambda_t
!        intermediate tensors (same notation as in relccsd)
!        temporarily using tensors with permutations required by the
!        automatically generated code
!        the next update to the code generator is expected to handle index
!        permutations/conjugation 
         type(tensor_t) :: vvoo_aux_tensor
         type(tensor_t) :: voov_aux_tensor 
         type(tensor_t) :: vooo_aux_tensor 
         type(tensor_t) :: vvov_aux_tensor 
         type(tensor_t) :: ovvo_aux_tensor
         type(tensor_t) :: vvvo_aux_tensor
         type(tensor_t) :: oovo_aux_tensor
         type(tensor_t) :: ovoo_aux_tensor
         type(tensor_t) :: ovov_aux_tensor
         type(tensor_t) :: ovvv_aux_tensor

         type(tensor_t) :: hvo_aux_tensor
!        solution tensors
         type(tensor_t) :: s1_tensor
         type(tensor_t) :: s2_tensor
         type(tensor_t) :: s3_tensor
         type(tensor_t) :: s4_tensor
         type(tensor_t) :: s5_tensor
         type(tensor_t) :: result_tensor
         integer(INTD)      :: result_dims(1)
         !Tensor dimensions
         integer(C_INT)     :: one_dims(1)
         !Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims
         integer(INTD), dimension(2) :: ov_dims
         integer(INTD), dimension(2) :: vo_dims
         integer(INTD), dimension(2) :: vv_dims
         integer(INTD), dimension(4) :: oovo_dims
         integer(INTD), dimension(4) :: vvoo_dims
         integer(INTD), dimension(4) :: ovoo_dims
         integer(INTD), dimension(4) :: oovv_dims
         integer(INTD), dimension(4) :: ovvo_dims
         integer(INTD), dimension(4) :: voov_dims
         integer(INTD), dimension(4) :: ovov_dims
         integer(INTD), dimension(4) :: vvov_dims
         integer(INTD), dimension(4) :: ovvv_dims
         integer(INTD), dimension(4) :: vooo_dims
         integer(INTD), dimension(4) :: ooov_dims
         integer(INTD), dimension(4) :: vovo_dims
         integer(INTD), dimension(4) :: vovv_dims
         integer(INTD), dimension(4) :: vvvv_dims
         integer(INTD), dimension(4) :: vvvo_dims
         integer(INTD), dimension(4) :: oooo_dims
         integer(INTD), dimension(6) :: ooovvv_dims
         integer(INTD), dimension(8) :: oooovvvv_dims
         integer(INTD), dimension(10):: ooooovvvvv_dims
         !Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: ilambda
         integer, intent(in) :: ncycles
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         !real(8), intent(out) :: ccsd_energy
         logical, intent(in)  :: CCD
         integer, intent(in)  :: iterative_level
         integer, intent(in)  :: iterative_scheme
         integer, intent(in)  :: auto_scheme
         integer, intent(in)  :: print_level
         integer (INTD)       :: ierr
          

         print*,'**************************************************'
         print*,' Automatically generated lambda residual equations'
         print*,'**************************************************'
        
         one_dims=1
         ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
 

         if (.not.CCD) then
            oo_dims = (/nocc, nocc/)
            vo_dims = (/nvir, nocc/)
            ov_dims = (/nocc, nvir/)
            vv_dims = (/nvir, nvir/) 

            ierr=construct_tensor_t(s1_tensor,C8,ov_dims,init_val=ZERO)
         
            oovo_dims = (/nocc, nocc, nvir, nocc/)
            ierr=construct_tensor_t(oovo_aux_tensor,C8,oovo_dims,init_val=ZERO)
            ierr=initialize_tensor_t(oovo_aux_tensor)
            ierr=t_einsum("T(i,j,a,k)+=V(i,j,k,a)",oovo_aux_tensor,int_t%ooov,one_tensor,scale=MINUS_ONE)

           ! write(*,*) 'oovo_aux_tensor=', tensor_norm2(oovo_aux_tensor)
            ovoo_dims = (/nocc, nvir, nocc, nocc/)
            ierr=construct_tensor_t(ovoo_aux_tensor,C8,ovoo_dims,init_val=ZERO)
            ierr=initialize_tensor_t(ovoo_aux_tensor)
            ierr=t_einsum("T(k,a,i,j)+=V+(i,j,k,a)",ovoo_aux_tensor,int_t%ooov,one_tensor)
           ! write(*,*) 'ovoo_aux_tensor=', tensor_norm2(ovoo_aux_tensor)

            ovov_dims = (/nocc, nvir, nocc, nvir/)
            ierr=construct_tensor_t(ovov_aux_tensor,C8,ovov_dims,init_val=ZERO)
            ierr=initialize_tensor_t(ovov_aux_tensor)
            ierr=t_einsum("T(i,a,j,b)+=V(a,i,b,j)",ovov_aux_tensor,int_t%vovo,one_tensor)
           ! write(*,*) 'ovov_aux_tensor=', tensor_norm2(ovov_aux_tensor)

            ovvv_dims = (/nocc, nvir, nvir, nvir/)
            ierr=construct_tensor_t(ovvv_aux_tensor,C8,ovvv_dims,init_val=ZERO)
            ierr=initialize_tensor_t(ovvv_aux_tensor)
            ierr=t_einsum("T(i,a,b,c)+=V(a,i,b,c)",ovvv_aux_tensor,int_t%vovv,one_tensor,scale=MINUS_ONE)
           ! write(*,*) 'ovvv_aux_tensor=', tensor_norm2(ovvv_aux_tensor)
            vvvo_dims = (/nvir, nvir, nvir, nocc/)
            ierr=construct_tensor_t(vvvo_aux_tensor,C8,vvvo_dims,init_val=ZERO)
            ierr=initialize_tensor_t(vvvo_aux_tensor)
            ierr=t_einsum("T(a,b,c,i)+=V+(c,i,a,b)",vvvo_aux_tensor,int_t%vovv,one_tensor)
           ! write(*,*) 'vvvo_aux_tensor=', tensor_norm2(vvvo_aux_tensor)
         end if


         oovv_dims = (/nocc, nocc, nvir, nvir/)
         ierr=construct_tensor_t(s2_tensor,C8,oovv_dims,init_val=ZERO)

         ovvo_dims = (/nocc, nvir, nvir, nocc/)
         ierr=construct_tensor_t(ovvo_aux_tensor,C8,ovvo_dims,init_val=ZERO)
         ierr=initialize_tensor_t(ovvo_aux_tensor)
         ierr=t_einsum("T(i,a,b,j)+=V(a,i,b,j)", ovvo_aux_tensor,int_t%vovo,one_tensor,scale=MINUS_ONE)
         write(*,*) 'ovvo_aux_tensor=', tensor_norm2(ovvo_aux_tensor)

         if (iterative_level.ge.3) then
            ooovvv_dims = (/nocc, nocc, nocc, nvir, nvir, nvir/)
            ierr=construct_tensor_t(s3_tensor,C8,ooovvv_dims,init_val=ZERO)
         end if
         if (iterative_level.ge.4) then 
            oooovvvv_dims = (/nocc, nocc, nocc, nocc, nvir, nvir, nvir, nvir/)
            ierr=construct_tensor_t(s4_tensor,C8,oooovvvv_dims,init_val=ZERO)
         end if
         if (iterative_level.ge.5) then 
            ooooovvvvv_dims = (/nocc,nocc, nocc, nocc, nocc, nvir,nvir, nvir, nvir, nvir/)
            ierr=construct_tensor_t(s5_tensor,C8,ooooovvvvv_dims,init_val=ZERO)
         end if
       
           if (.not.CCD) then
            if (iterative_level.eq.5) then
               call diis (t5_tensor=comm_t%l5,t4_tensor=comm_t%l4,t3_tensor=comm_t%l3, &
                          t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=1)
            else if (iterative_level.eq.4) then
               call diis (t4_tensor=comm_t%l4,t3_tensor=comm_t%l3,t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,&
                           mode=1)
            else if (iterative_level.eq.3) then
               call diis (t3_tensor=comm_t%l3,t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=1)
            else
               call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=1)
            end if
         else
            call diis (t2_tensor=comm_t%l2,mode=1)
         end if

         do ilambda = 1, ncycles

          if (print_level.gt.2) then
!         if (.true.) then
            if (.not.CCD) then 
               write(*,*) "l1    = ",tensor_norm2(comm_t%l1)
               if (iterative_level.ge.3) write(*,*) "l3    = ",tensor_norm2(comm_t%l3)
               if (iterative_level.ge.4) write(*,*) "l4    = ",tensor_norm2(comm_t%l4)
               if (iterative_level.ge.5) write(*,*) "l5    = ",tensor_norm2(comm_t%l5)
            end if
            write(*,*) "l2    = ",tensor_norm2(comm_t%l2)
          end if
          if (print_level.gt.2) then
            if (.not.CCD) then
               write(*,*) "s1    = ",tensor_norm2(s1_tensor)
               if (iterative_level.ge.3) write(*,*) "s3    = ",tensor_norm2(s3_tensor)
               if (iterative_level.ge.4) write(*,*) "s4    = ",tensor_norm2(s4_tensor)
               if (iterative_level.ge.5) write(*,*) "s5    = ",tensor_norm2(s5_tensor)
            end if
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)
          end if
          if (ierr.ne.0) then
             write(*,*) " ccdriver: error in calculation of common blocks,code=",ierr
             return
          end if

       !   ierr=initialize_tensor_t(result_tensor)
          if (.not.CCD) then
             ierr=initialize_tensor_t(s1_tensor)
             if (iterative_level.ge.3) ierr=initialize_tensor_t(s3_tensor)
             if (iterative_level.ge.4) ierr=initialize_tensor_t(s4_tensor)
             if (iterative_level.ge.5) ierr=initialize_tensor_t(s5_tensor)
          end if
          ierr=initialize_tensor_t(s2_tensor)
          if ((iterative_level.eq.2).and.(iterative_scheme.eq.21)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            !call talsh_CCSD_codegen()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then  
            ! tenpi
               call eti_tenpi_lambda_ccsd(nocc,nvir,&
                                            comm_t%foo, &       !OO
                                            comm_t%fov, &       !OV
                                            comm_t%fvv, &       !VV
                                            comm_t%l2, &        !l2_OOVV
                                            comm_t%l1, &        !l1_OV
                                            comm_t%t1, &        !t1_VO
                                            comm_t%t2, &        !t2_VVOO
                                            int_t%oooo,&        !OOOO
                                            int_t%ooov,&        !OOOV
                                            oovo_aux_tensor,&   !OOVO
                                            int_t%oovv,&        !OOVV
                                            ovoo_aux_tensor,&   !OVOO 
                                            ovov_aux_tensor,&   !OVOV
                                            ovvo_aux_tensor, &  !OVVO
                                            ovvv_aux_tensor, &  !OVVV
                                            int_t%vovo, &       !VOVO
                                            int_t%vovv, &       !VOVV
                                            vvvo_aux_tensor, &  !VVVO
                                            int_t%vvvv, &       !VVVV
                                            s2_tensor, &        !r_OOVV
                                            s1_tensor &         !r_OV    
                                                                   )
           else if (auto_scheme.eq.3) then 
             call tapp_tenpi_lambda_ccsd(nocc,nvir,&
                                            comm_t%foo, &       !OO
                                            comm_t%fov, &       !OV
                                            comm_t%fvv, &       !VV
                                            comm_t%l2, &        !l2_OOVV
                                            comm_t%l1, &        !l1_OV
                                            comm_t%t1, &        !t1_VO
                                            comm_t%t2, &        !t2_VVOO
                                            int_t%oooo,&        !OOOO
                                            int_t%ooov,&        !OOOV
                                            oovo_aux_tensor,&   !OOVO
                                            int_t%oovv,&        !OOVV
                                            ovoo_aux_tensor,&   !OVOO 
                                            ovov_aux_tensor,&   !OVOV
                                            ovvo_aux_tensor, &  !OVVO
                                            ovvv_aux_tensor, &  !OVVV
                                            int_t%vovo, &       !VOVO
                                            int_t%vovv, &       !VOVV
                                            vvvo_aux_tensor, &  !VVVO
                                            int_t%vvvv, &       !VVVV
                                            s2_tensor, &        !r_OOVV
                                            s1_tensor &         !r_OV    
                                                                   )
            end if
          else if ((iterative_level.eq.2).and.(iterative_scheme.eq.20)) then
            if (auto_scheme.eq.1) then       
            ! mb-codegen
            !call talsh_CCD_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then 
            ! tenpi
               call eti_tenpi_lambda_ccd(nocc,nvir,&
                                                     comm_t%foo, & !OO
                                                     comm_t%fvv, & !VV
                                                     comm_t%l2, & !l2_OOVV
                                                     comm_t%t2, & ! t2_VVOO
                                                     int_t%oooo,& ! OOOO
                                                     int_t%oovv,& ! OOVV
                                                     ovvo_aux_tensor, & !OVVO
                                                     int_t%vvvv, &      !VVVV
                                                     s2_tensor &    !r_OOVV
                                                        )
            else if (auto_scheme.eq.3) then 
             call tapp_tenpi_lambda_ccd(nocc,nvir,&
                                                     comm_t%foo, & !OO
                                                     comm_t%fvv, & !VV
                                                     comm_t%l2, & !l2_OOVV
                                                     comm_t%t2, & ! t2_VVOO
                                                     int_t%oooo,& ! OOOO
                                                     int_t%oovv,& ! OOVV
                                                     ovvo_aux_tensor, & !OVVO
                                                     int_t%vvvv, &      !VVVV
                                                     s2_tensor &    !r_OOVV
                                                        )

            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.30)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            ! call talsh_CCSDT_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then 
            ! tenpi
               call eti_tenpi_lambda_ccsdt(nocc,nvir,&
                                                     comm_t%foo, &       !OO
                                                     comm_t%fov, &       !OV
                                                     comm_t%fvv, &       !VV
                                                     comm_t%l3, &        !l3_OOOVVV
                                                     comm_t%l2, &        !l2_OOVV
                                                     comm_t%l1, &        !l1_OV
                                                     comm_t%t1, &        !t1_VO
                                                     comm_t%t2, &        !t2_VVOO
                                                     comm_t%t3, &        !t3_VVVOOO
                                                     int_t%oooo,&        !OOOO
                                                     int_t%ooov,&        !OOOV
                                                     oovo_aux_tensor,&   !OOVO
                                                     int_t%oovv,&        !OOVV
                                                     ovoo_aux_tensor,&   !OVOO 
                                                     ovov_aux_tensor,&   !OVOV                   
                                                     ovvo_aux_tensor, &  !OVVO
                                                     ovvv_aux_tensor, &  !OVVV
                                                     int_t%vovo, &       !VOVO
                                                     int_t%vovv, &       !VOVV
                                                     vvvo_aux_tensor, &  !VVVO
                                                     int_t%vvvv, &       !VVVV
                                                     s3_tensor, &        !r_OOOVVV
                                                     s2_tensor, &        !r_OOVV  
                                                     s1_tensor  )        !r_OV
                                                
                                                 
            else if (auto_scheme.eq.3) then 
             call tapp_tenpi_lambda_ccsdt(nocc,nvir,&
                                                     comm_t%foo, &       !OO
                                                     comm_t%fov, &       !OV
                                                     comm_t%fvv, &       !VV
                                                     comm_t%l3, &        !l3_OOOVVV
                                                     comm_t%l2, &        !l2_OOVV
                                                     comm_t%l1, &        !l1_OV
                                                     comm_t%t1, &        !t1_VO
                                                     comm_t%t2, &        !t2_VVOO
                                                     comm_t%t3, &        !t3_VVVOOO
                                                     int_t%oooo,&        !OOOO
                                                     int_t%ooov,&        !OOOV
                                                     oovo_aux_tensor,&   !OOVO
                                                     int_t%oovv,&        !OOVV
                                                     ovoo_aux_tensor,&   !OVOO 
                                                     ovov_aux_tensor,&   !OVOV                   
                                                     ovvo_aux_tensor, &  !OVVO
                                                     ovvv_aux_tensor, &  !OVVV
                                                     int_t%vovo, &       !VOVO
                                                     int_t%vovv, &       !VOVV
                                                     vvvo_aux_tensor, &  !VVVO
                                                     int_t%vvvv, &       !VVVV
                                                     s3_tensor, &        !r_OOOVVV
                                                     s2_tensor, &        !r_OOVV  
                                                     s1_tensor  )        !r_OV
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.31)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            ! call talsh_CCSDT_1a_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.32)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
            ! call talsh_CCSDT_1b_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.33)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
           !  call talsh_CCSDT_2_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.3).and.(iterative_scheme.eq.34)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
           !  call talsh_CCSDT_3_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
            end if
          else if ((iterative_level.eq.4).and.(iterative_scheme.eq.40)) then 
            if (auto_scheme.eq.1) then
            ! mb-codegen
          !   call talsh_CCSDTQ_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
               call eti_tenpi_lambda_ccsdtq(nocc,nvir,&
                                                     comm_t%foo, &       !OO
                                                     comm_t%fov, &       !OV
                                                     comm_t%fvv, &       !VV
                                                     comm_t%l4, & !l4_OOOOVVVV
                                                     comm_t%l3, & !l3_OOOVVV
                                                     comm_t%l2, & !l2_OOVV
                                                     comm_t%l1, &        !l1_OV
                                                     comm_t%t1, &        !t1_VO
                                                     comm_t%t2, & !t2_VVOO
                                                     comm_t%t3, & !t3_VVVOOO
                                                     comm_t%t4, & !t4_VVVVOOOO
                                                     int_t%oooo,&        !OOOO
                                                     int_t%ooov,&        !OOOV
                                                     oovo_aux_tensor,&   !OOVO
                                                     int_t%oovv,&        !OOVV
                                                     ovoo_aux_tensor,&   !OVOO 
                                                     ovov_aux_tensor,&   !OVOV                   
                                                     ovvo_aux_tensor, &  !OVVO
                                                     ovvv_aux_tensor, &  !OVVV
                                                     int_t%vovo, &       !VOVO
                                                     int_t%vovv, &       !VOVV
                                                     vvvo_aux_tensor, &  !VVVO
                                                     int_t%vvvv, &       !VVVV
                                                     s4_tensor, & !r_OOOOVVVV
                                                     s3_tensor, & !r_OOOVVV
                                                     s2_tensor, &        !r_OOVV  
                                                     s1_tensor  )        !r_OV
            else if (auto_scheme.eq.3) then 
              call tapp_tenpi_lambda_ccsdtq(nocc,nvir,&
                                                     comm_t%foo, &       !OO
                                                     comm_t%fov, &       !OV
                                                     comm_t%fvv, &       !VV
                                                     comm_t%l4, & !l4_OOOOVVVV
                                                     comm_t%l3, & !l3_OOOVVV
                                                     comm_t%l2, & !l2_OOVV
                                                     comm_t%l1, &        !l1_OV
                                                     comm_t%t1, &        !t1_VO
                                                     comm_t%t2, & !t2_VVOO
                                                     comm_t%t3, & !t3_VVVOOO
                                                     comm_t%t4, & !t4_VVVVOOOO
                                                     int_t%oooo,&        !OOOO
                                                     int_t%ooov,&        !OOOV
                                                     oovo_aux_tensor,&   !OOVO
                                                     int_t%oovv,&        !OOVV
                                                     ovoo_aux_tensor,&   !OVOO 
                                                     ovov_aux_tensor,&   !OVOV                   
                                                     ovvo_aux_tensor, &  !OVVO
                                                     ovvv_aux_tensor, &  !OVVV
                                                     int_t%vovo, &       !VOVO
                                                     int_t%vovv, &       !VOVV
                                                     vvvo_aux_tensor, &  !VVVO
                                                     int_t%vvvv, &       !VVVV
                                                     s4_tensor, & !r_OOOOVVVV
                                                     s3_tensor, & !r_OOOVVV
                                                     s2_tensor, &        !r_OOVV  
                                                     s1_tensor  )        !r_OV

            end if
          else if ((iterative_level.eq.5).and.(iterative_scheme.eq.50)) then
            if (auto_scheme.eq.1) then
            ! mb-codegen
          !   call talsh_CCSDTQP_codegen_call()
               call quit('Lambda equations not implemented with mbauto')
            else if (auto_scheme.eq.2) then
            ! tenpi
               call quit ('Error: CCSDTQP not yet implemented in tenpi.')
            end if
          else 
             write(*,*) "Implementation not available for this method"
          end if

          if (print_level.gt.2) then
!         if (.true.) then
            if (.not.CCD) then
               write(*,*) "s1    = ",tensor_norm2(s1_tensor)
               if (iterative_level.ge.3)  write(*,*) "s3    = ",tensor_norm2(s3_tensor)
               if (iterative_level.ge.4)  write(*,*) "s4    = ",tensor_norm2(s4_tensor)
               if (iterative_level.ge.5)  write(*,*) "s5    = ",tensor_norm2(s5_tensor)
            end if
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)
          end if

!         prepare for next iteration
                

          if (.not.CCD) then
             ierr=initialize_tensor_t(comm_t%l1)
             ierr=t_einsum("L1(i,a)+=S1(i,a)",comm_t%l1,s1_tensor)
             if (iterative_level.ge.3) then
                ierr=initialize_tensor_t(comm_t%l3)
                ierr=t_einsum("L3(i,j,k,a,b,c)+=S3(i,j,k,a,b,c)",comm_t%l3,s3_tensor)
             end if
             if (iterative_level.ge.4) then
                ierr=initialize_tensor_t(comm_t%l4)
                ierr=t_einsum("L4(i,j,k,l,a,b,c,d)+=S4(i,j,k,l,a,b,c,d)",comm_t%l4,s4_tensor)
             end if
             if (iterative_level.ge.5) then
                ierr=initialize_tensor_t(comm_t%l5)
                ierr=t_einsum("L5(i,j,k,l,m,a,b,c,d,e)+=S5(i,j,k,l,m,a,b,c,d,e)",comm_t%l5,s5_tensor)
             end if
          end if
          ierr=initialize_tensor_t(comm_t%l2)
          ierr=t_einsum("L2(i,j,a,b)+=S2(i,j,a,b)",comm_t%l2,s2_tensor)
       
          if (.not.CCD) then
            if (iterative_level.ge.5) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%l1,&
                              t2_tensor=comm_t%l2,full_t3_tensor=comm_t%l3,&
                              full_t4_tensor=comm_t%l4,full_t5_tensor=comm_t%l5)
            else if (iterative_level.ge.4) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%l1,&
                              t2_tensor=comm_t%l2,full_t3_tensor=comm_t%l3,full_t4_tensor=comm_t%l4)
            else if (iterative_level.ge.3) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%l1,&
                              t2_tensor=comm_t%l2,full_t3_tensor=comm_t%l3)
            else
               call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%l1,&
                              t2_tensor=comm_t%l2)
            end if
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%l2)
          end if


        if (.not.CCD) then
            if (iterative_level.ge.5) then
               call diis (t5_tensor=comm_t%l5,t4_tensor=comm_t%l4,t3_tensor=comm_t%l3,&
                          t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,t_convergence=t_convergence)
            else if (iterative_level.ge.4) then
               call diis (t4_tensor=comm_t%l4,t3_tensor=comm_t%l3,t2_tensor=comm_t%l2,&
                        & t1_tensor=comm_t%l1,t_convergence=t_convergence)
            else if (iterative_level.ge.3) then
               call diis (t3_tensor=comm_t%l3,t2_tensor=comm_t%l2,&
                        t1_tensor=comm_t%l1,t_convergence=t_convergence)
            else
               call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,t_convergence=t_convergence)
            end if
          else
            call diis (t2_tensor=comm_t%l2,t_convergence=t_convergence)
          end if


        call print_iteration(ilambda,t_convergence)

          if (t_convergence .lt. t_target_precision) then 
                write(*,*) 'CONVERGENCE' 

                exit
          end if
                
         end do

        write(*,*) "------------------------"
         write(*,*) ""



         
        
         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(comm_t%l1)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(comm_t%l2)*ONE_QUARTER_R

         write(*,*) ""

         if (.not.CCD) then
            ierr=destroy_tensor_t(s1_tensor)
            ierr=destroy_tensor_t(hvo_aux_tensor)
            ierr=destroy_tensor_t(vooo_aux_tensor)
            ierr=destroy_tensor_t(vvov_aux_tensor)
            if (iterative_level.ge.3) ierr=destroy_tensor_t(s3_tensor)
            if (iterative_level.ge.4) ierr=destroy_tensor_t(s4_tensor)
            if (iterative_level.ge.5) ierr=destroy_tensor_t(s5_tensor)
         end if
         ierr=destroy_tensor_t(vvoo_aux_tensor)
         ierr=destroy_tensor_t(voov_aux_tensor)
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(result_tensor)
         ierr=destroy_tensor_t(one_tensor)

        end subroutine 
!
end module 
