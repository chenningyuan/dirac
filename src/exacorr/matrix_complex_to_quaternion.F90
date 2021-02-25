module matrix_complex_to_quaternion
!This module contains the routines needed to calculate the mp2 natural orbitals 

        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use exacorr_mo
        use exacorr_global,only: nfunctions
        use, intrinsic:: ISO_C_BINDING

        implicit none
        private

        interface basis_transform
         module procedure basis_transform_vector
         module procedure basis_transform_matrix
        end interface

        public read_matrix_comp,sort_values_vectors,basis_transform,diagoanl_matrix_quater, &
               complex_to_quaternion,get_Cano_Orbital

       contains
!-------------------------------------------------------------------------------------------------------------------------------------
 
       subroutine read_matrix_comp(Matrix, Matrix_name)
!       Written by Xiang Yuan, Dec 2020.

        complex(kind=8), pointer, intent(inout)  :: Matrix(:,:)   ! complex matrix from DM_complex
        character(len=*), intent(in)             :: Matrix_name
        logical :: debug=.true.

        integer :: i, j, k, ierr
        integer, parameter :: cm = 10
        real(kind=8)       :: a,b

        ! first round for reading matrix size
        open (cm,file=Matrix_name,status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
        if (debug) then
          print *,'<<<<< Complex Matrix Existed '
          ierr = 0
          do while(ierr==0)
                read(cm,'(2x,2F28.16,2x,2I4)',iostat=ierr) a,b,i,j
          end do
           print *,'>>>>> Read Matrix Index'
        end if
        close(cm,status='KEEP')

        allocate(Matrix(i,j))

        open (cm,file=Matrix_name,status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
          ierr = 0
          do while(ierr==0)
                read(cm,'(2x,2F28.16,2x,2I4)',iostat=ierr) a,b,i,j
                Matrix(i,j) = cmplx(a,b)
          end do
           print *,'>>>>> Read Matrix'
        close(cm,status='KEEP')



      end subroutine read_matrix_comp


 !-----------------------------------------------------------------------------------------------------------------------------------------------

      subroutine sort_values_vectors (NOs_OCC_quat,NOs_MO_quat)

!       This subroutine sorts the natural orbital based on occupation number 
!       Written by Xiang Yuan, July 2020

         !input variables
         real(kind=8), pointer, intent(inout)       :: NOs_OCC_quat(:)          ! occupation number from qdia of DM
         real(kind=8), pointer, intent(inout)       :: NOs_MO_quat(:,:,:)      ! natural orbital from qdia of DM


         real(kind=8)      :: temp
         real, allocatable :: temp_arr(:,:,:)
         integer           :: i,j,sz


         sz = size(NOs_OCC_quat,1)
         allocate(temp_arr(sz,1,4))

         ! sorting occupation number in descending order 
         do i = 1, sz/2
            temp = NOs_OCC_quat(i)
            NOs_OCC_quat(i) = NOs_OCC_quat(sz+1-i)
            NOs_OCC_quat(sz+1-i) = temp
         end do

         ! sorting natural orbital in descending order 
         do i = 1, sz/2
            temp_arr(:,1,:) = NOs_MO_quat(:,i,:)
            NOs_MO_quat(:,i,:) = NOs_MO_quat(:,sz+1-i,:)
            NOs_MO_quat(:,sz+1-i,:) = temp_arr(:,1,:)
         end do
         
        end subroutine sort_values_vectors

!-----------------------------------------------------------------------------------------------------------------------------------------

        subroutine basis_transform_vector (mat_basis1_in,basis2_in,mat_basis2_out,fl)
!        Routine to do basis transform by quaternion matrix multiplication(vectors).
!        We follows V' = U * V.
!        V': vectors in new basis;  V: vectors in old basis; U: transform matrix;        
                
           !input variables
           real(kind=8), pointer, intent(inout)   :: mat_basis1_in(:,:,:)    ! vectors in old basis
           real(kind=8), pointer, intent(inout)   :: basis2_in(:,:,:)        ! transform matrix
           real(kind=8), pointer, intent(inout)   :: mat_basis2_out(:,:,:)   ! tvectors in new basis
           integer, intent(inout)                 :: fl



           integer                   :: LRA, LCA, NZA
           integer                   :: LRB, LCB, NZB
           integer                   :: LRC, LCC, NZC

           integer                   :: re



           LRA = size(basis2_in,1)
           LCA = size(basis2_in,2)
           LRB = size(mat_basis1_in,1)
           LCB = size(mat_basis1_in,2)
           LRC = size(basis2_in,1)
           LCC = size(mat_basis1_in,2)


           NZA = 4
           NZB = 4
           NZC = 4

           allocate(mat_basis2_out(LRC,LCC,4))

           call qua_mat_mulp(basis2_in,mat_basis1_in,mat_basis2_out,fl,LRA,LCA,NZA,LRB,LCB,NZB,LRC,LCC,NZC)
                      
        end subroutine basis_transform_vector


!-----------------------------------------------------------------------------------------------------------------------------------

        subroutine basis_transform_matrix (mat_basis1_in,basis2_in,mat_basis2_out,mat)
!        Routine to do basis transform by quaternion matrix multiplication(matrix).
!        We follows C = A^{-1}BA
!        C: matrix in new basis;  B: matrix in old basis; A: change of basis matrix;       
                         
         !input variables
         real(kind=8), pointer, intent(inout)   :: mat_basis1_in(:,:,:)    ! matrix in old basis(matrix_B)
         real(kind=8), pointer, intent(inout)   :: basis2_in(:,:,:)        ! change of basis matrix(matrix_A)
         real(kind=8), pointer, intent(inout)   :: mat_basis2_out(:,:,:)   ! matrix in new basis(matrix_C)
         character                              :: mat                     ! flag used for matrix basis transform

         integer                                :: re,fl
         real(kind=8), pointer                  :: aux(:,:,:)              ! matrix used in qgemm
         integer                                :: LRA, LCA, NZA
         integer                                :: LRB, LCB, NZB
         integer                                :: LR_au, LC_au, NZ_au
         integer                                :: LRC, LCC, NZC


         LRA = size(basis2_in,1)
         LCA = size(basis2_in,2)
         LRB = size(mat_basis1_in,1)
         LCB = size(mat_basis1_in,2)
         LR_au = size(basis2_in,2)
         LC_au = size(mat_basis1_in,2)
         LRC = size(basis2_in,2)
         LCC = size(basis2_in,2)

         NZA =   4
         NZB =   4
         NZ_au = 4 
         NZC =   4

         allocate(aux(LR_au,LC_au,NZ_au))
         allocate(mat_basis2_out(LRC,LCC,4))
         
         fl = 1
         call qua_mat_mulp(basis2_in,mat_basis1_in,aux,fl,LRA,LCA,NZA,LRB,LCB,NZB,LR_au,LC_au,NZ_au)
         fl = 0
         call qua_mat_mulp(aux,basis2_in,mat_basis2_out,fl,LR_au,LC_au,NZ_au,LRA,LCA,NZA,LRC,LCC,NZC) 


         deallocate(aux)

        end subroutine basis_transform_matrix


!-----------------------------------------------------------------------------------------------------------------------------------
        
      subroutine qua_mat_mulp (A,B,C,flag,LRA,LCA,NZA,LRB,LCB,NZB,LRC,LCC,NZC)
!        Routine to do quaternion matrix multiplication. We use qgemm routine from qalebra.F.
                         
         !input variables         
         real(kind=8), pointer, intent(inout)  :: A(:,:,:), B(:,:,:), C(:,:,:)               ! matrix used in qgemm
         integer, intent(in)                   :: flag
         integer, intent(in)                   :: LRA, LCA, NZA
         integer, intent(in)                   :: LRB, LCB, NZB
         integer, intent(in)                   :: LRC, LCC, NZC 


         integer                   :: M, N, K                              !paramter for blas routine dgemm

         integer, allocatable      :: IQPA(:), IQPB(:), IQPC(:)        
         real(kind=8)              :: ALPHA, BETA
         character                 :: FA, FB ,TA ,TB
         integer                   :: i, tem
            
         M = LRA
         N = LCB
         K = LCA


         FA = 'N'
         FB = 'N'
         TA = 'N'
         TB = 'N'
         
         ALPHA = 1.0
         BETA = 0.0
         
         
         allocate(IQPA(NZA))
         allocate(IQPB(NZB))
         allocate(IQPC(NZC)) 

         do i = 1, NZA
            IQPA(i) = i
            IQPB(i) = i
            IQPC(i) = i
         end do            

!       There three options: case1 does C=A^{dagger}*B; case2 does C=A*B
         if (flag == 1) then

            FA = 'H'
            tem = M
            M = K
            K = tem

            call qgemm(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA, & 
                     FB,TB,IQPB,B,LRB,LCB,NZB,BETA, &
                     IQPC,C,LRC,LCC,NZC)      

         else if (flag == -1) then            
            FB = 'H'
            tem = N
            N = K
            K = tem
         
            call qgemm(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA, & 
                        FB,TB,IQPB,B,LRB,LCB,NZB,BETA, &
                        IQPC,C,LRC,LCC,NZC)     

         else 
            call qgemm(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA, & 
                       FB,TB,IQPB,B,LRB,LCB,NZB,BETA, &
                       IQPC,C,LRC,LCC,NZC)    
         end if


         deallocate(IQPA)
         deallocate(IQPB)
         deallocate(IQPC)
                               

      end subroutine qua_mat_mulp
         
!-----------------------------------------------------------------------------------------------------------------------------------                                                   
        subroutine diagoanl_matrix_quater(A,EIG,VEC)

         ! Set parameter needed in qdiag
                 
                 real(kind=8), intent(in), pointer :: A(:,:,:)
                 real(kind=8), pointer :: EIG(:)
                 real(kind=8), pointer :: VEC(:,:,:)                 
                 integer :: LRA, LCA, LRV, LCV
                 integer :: NZ, N, MATZ, LWORK, IERR
                 integer :: i, j
                 real(kind=8), allocatable :: WORK(:)
                 
               
                 LRA = size(A,1)
                 LCA = size(A,2)
                 LRV = LRA
                 LCV = LCA
                 !lwork = 1000000*LRA
                 lwork = 10*LRA
                 NZ = 4
                 N = LRA
                 MATZ = 1
                 allocate(EIG(N))
                 allocate(VEC(LRV,LCV,NZ))
                 allocate(work(MAX(1,LWORK)))


                 call qdiag (NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,WORK,LWORK,IERR)


               deallocate(work)


      end subroutine diagoanl_matrix_quater
                
               
!---------------------------------------------------------------------------------------------------------------------------------------------------
            
      subroutine complex_to_quaternion(exa_input,dm_mo_complex,dm_mo_quater)
!               use exacorr_global, only : get_nsolutions

               !input variables
               type(exacc_input),   intent(in)     :: exa_input
               complex(8), pointer, intent(inout)  :: dm_mo_complex(:,:)
               real(8),pointer,     intent(inout)  :: dm_mo_quater(:,:,:)

               !matrix dimensions
               integer :: nocc, nvir, nkr_occ, nkr_vir, nesh
               !complete mo list
               integer, allocatable :: mo_list(:)
               !error code
               integer :: ierr
               !loop variables
               integer :: i, j, norbt
               integer :: i_unbar, j_unbar, j_bar
               !debug
               logical :: debug=.true.
               integer :: k


!***************
!     set dims *
!***************
 
!      nocc    = exa_input%nocc
!      nvir    = exa_input%nvir
!      nkr_occ = exa_input%nkr_occ
!      nkr_vir = exa_input%nkr_vir
!      norbt = nkr_occ + nkr_vir

      norbt = size(dm_mo_complex,1)/2
      allocate(mo_list(norbt))
      do i = 1, norbt
            mo_list(i) = i
      end do
      
!      nesh = get_nsolutions() / 2
      nesh = norbt
      allocate (dm_mo_quater(nesh,nesh,4))
      dm_mo_quater = 0
      do j = 1, norbt
        j_unbar = 2*j-1
        j_bar   = 2*j
        do i = 1, norbt
          ! We copy the unbar-unbar part of the density matrix, as well as the unbar-bar part
          ! This follows from eq. 27 of Shee, Visscher, Saue JCP 145 184107, 2016
          ! it would be good to test whether the Kramers-restricted assumption holds
          ! Or even to average over the two density matrices and enforce a Kramers-restricted density
          i_unbar = 2*i-1
          dm_mo_quater(mo_list(i),mo_list(j),1) = real(dm_mo_complex(i_unbar,j_unbar))
          dm_mo_quater(mo_list(i),mo_list(j),2) = imag(dm_mo_complex(i_unbar,j_unbar))
          dm_mo_quater(mo_list(i),mo_list(j),3) = real(dm_mo_complex(i_unbar,j_bar))
          dm_mo_quater(mo_list(i),mo_list(j),4) = imag(dm_mo_complex(i_unbar,j_bar))
        end do
      end do

!      deallocate (dm_mo_quater)
      deallocate (mo_list)
!      deallocate (dm_mo_comple)
      end subroutine complex_to_quaternion

!------------------------------------------------------------------------------------------------------------------------------------------------------------------      
      subroutine get_Cano_Orbital(orbital_noncano,orbital_cano,Fock_NO,new_energy)
         ! This subroutine is used for constructing the canonical orbital

         !input variables
         real(kind=8), pointer,intent(inout)   :: orbital_noncano(:,:,:)             ! non-canonical orbital in HF mo basis
         real(kind=8), pointer,intent(out)     :: orbital_cano(:,:,:)                ! canonical orbital in HF mo basis
         real(kind=8), pointer,intent(out)     :: New_energy(:)                      ! new orbital energy from qdia of Fock matrix in non-canonical orbital basis
         real(kind=8), pointer,intent(inout)   :: Fock_NO(:,:,:)                     ! v-v block of Fock matrix in non-canonical orbital basis

         real(kind=8), pointer              :: orbital_cano_basis_noncano(:,:,:)  ! canonical orbital in non-canonical orbital basis
         integer                            :: fl, i, j

         ! Diagonalize new FM to obtain the new canonical orbital 
         call diagoanl_matrix_quater(Fock_NO, new_energy, orbital_cano_basis_noncano)
         ! Transform the new canonical orbital from non-canonical orbital basis into HF mo basis
         fl = 0
         call basis_transform(orbital_cano_basis_noncano, orbital_noncano, orbital_cano, fl)         

      end subroutine get_Cano_Orbital
 
!-----------------------------------------------------------------------------------------------------------------------------------   

end module matrix_complex_to_quaternion

