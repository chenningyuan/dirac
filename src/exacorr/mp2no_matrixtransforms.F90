module mp2no_matrixtransforms
! This module contains the routines needed to calculate the mp2 natural orbitals 
! Including quaternion matrix digonalization and multiplication
! Written by Xiang Yuan.

        implicit none

        private

        public basis_transform_vector, basis_transform_matrix
        public complex_to_quaternion, diagonal_matrix_quater, get_cano_orbital, read_matrix_comp!, sort_values_vectors
        public sort_MOs_by_occupation_number_v2, sort_MOs_by_occupation_number_v1

       contains

!------------------------------------------------------------------------------------------------------------------------------------- 
       subroutine read_matrix_comp(Matrix, Matrix_file)
!       read complex density matrix from 'DM_complex' 
        use exacorr_utils

        complex(kind=8), pointer, intent(inout)  :: Matrix(:,:)   ! complex matrix from DM_complex
        character(len=*), intent(in)             :: Matrix_file
        logical :: debug=.false.

        integer :: sz_all
        integer :: cm 

        call get_free_fileunit(cm)
        ! first round for reading matrix size
        open (cm,file=Matrix_file,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
        ! read dimension of density matrix 
        read(cm) sz_all
        allocate(Matrix(sz_all,sz_all))
        ! read density matrix element 
        read(cm) Matrix

        close(cm,status='KEEP')

      end subroutine read_matrix_comp

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine sort_MOs_by_occupation_number_v1(N_OCCUPATION,NOs_MO)
!     This subroutine sorts the natural orbital based on occupation number in descending order.  

         !input variables
         real(kind=8), pointer, intent(inout)       :: N_OCCUPATION(:)    ! occupation number from qdia of DM
         real(kind=8), pointer, intent(inout)       :: NOs_MO(:,:,:)      ! natural orbital from qdia of DM


         real(kind=8)      :: temp
         real, allocatable :: temp_arr(:,:,:)
         integer           :: i,j,sz1, sz2


         sz1 = size(NOs_MO,1)
         sz2 = size(NOs_MO,2)
         allocate(temp_arr(sz1,1,4))

         ! sorting occupation number in descending order 
         do i = 1, sz2/2
            temp = N_OCCUPATION(i)
            N_OCCUPATION(i) = N_OCCUPATION(sz2+1-i)
            N_OCCUPATION(sz2+1-i) = temp
         end do

         ! sorting natural orbital in descending order 
         do i = 1, sz2/2
            temp_arr(:,1,:) = NOs_MO(:,i,:)
            NOs_MO(:,i,:) = NOs_MO(:,sz2+1-i,:)
            NOs_MO(:,sz2+1-i,:) = temp_arr(:,1,:)
         end do
         
        end subroutine sort_MOs_by_occupation_number_v1

        !subroutine sort_MOs_by_occupation_number_v2(N_OCCUPATION, NOs_MO, nao, nto, nz)
        subroutine sort_MOs_by_occupation_number_v2(N_OCCUPATION, NOs_MO)
           use exacorr_utils
           ! Using merge sort algorithm for sorting FNO++
           !integer, intent(in) :: nao, nto, nz
           !real(kind=8), dimension(nao,nto,nz), intent(inout)  :: NOs_MO
           !real(kind=8), dimension(nto),        intent(inout)  :: N_OCCUPATION

           !real(kind=8), dimension(nao,nto,nz)      :: NOs_MO_tmp
           !real(kind=8), dimension(nto)             :: N_OCCUPATION_tmp !, N_OCCUPATION_original(:)
           !integer, dimension(nto)                  :: indices, indices_tmp

           !integer :: i, j, width, n, iLeft, iRight, iEnd
           !!!!!!!!!!!!!
           real(kind=8), pointer, intent(inout)  :: NOs_MO(:,:,:)
           real(kind=8), pointer,        intent(inout)  :: N_OCCUPATION(:)

           real(kind=8), dimension(:,:,:),allocatable      :: NOs_MO_tmp
           real(kind=8), dimension(:),allocatable             :: N_OCCUPATION_tmp !, N_OCCUPATION_original(:)
           integer, dimension(:), allocatable                  :: indices, indices_tmp

           integer :: i, j, width, n, iLeft, iRight, iEnd
           integer :: nao, nto, nz


           nao = size(NOs_MO,1)
           nto = size(NOs_MO,2)
           n = nto
           allocate(NOs_MO_tmp(nao,nto,4))
           allocate(indices(n))
           allocate(N_OCCUPATION_tmp(n))
           allocate(indices_tmp(n))

           n = nto

           N_OCCUPATION_tmp = N_OCCUPATION
           NOs_MO_tmp = NOs_MO

           indices = (/ (i, i=1,n) /)
           indices_tmp = (/ (i, i=1,n) /)

           ! sort N_OCCUPATION
           width = 1
           do while(width<n)

              do i = 1, n, 2*width
                 call bottomup_merge(N_OCCUPATION, i, min(i+width, n+1), min(i+width*2,n+1),N_OCCUPATION_tmp, indices, indices_tmp)
              enddo

              N_OCCUPATION = N_OCCUPATION_tmp
              indices = indices_tmp
              width = width * 2
           enddo

           ! check for debug
           !do i=1,n-1
           !   write(*,*) i, N_OCCUPATION(i),N_OCCUPATION(i+1)
           !   if(N_OCCUPATION(i)<N_OCCUPATION(i+1)) then
           !      call print_date("Not sorted!")
           !      STOP
           !   endif
           !enddo

           ! sort NOs_MO with the indices
           do i=1,n
              NOs_MO(:,i,:) = NOs_MO_tmp(:,indices(i),:)
           enddo
           
           deallocate(indices)
           deallocate(indices_tmp)
           deallocate(NOs_MO_tmp)
           deallocate(N_OCCUPATION_tmp)

           call print_date("Finished sorting NOs!")
           return
        end subroutine sort_MOs_by_occupation_number_v2

        subroutine bottomup_merge(arr, iLeft, iRight, iEnd, arr_tmp, indices, indices_tmp)
           real(kind=8), dimension(:), intent(in) :: arr
           integer,      dimension(:), intent(in) :: indices
           real(kind=8), dimension(:), intent(inout) :: arr_tmp
           integer,      dimension(:), intent(inout) :: indices_tmp
           integer, intent(in) :: iLeft, iRight, iEnd
           integer :: k, i, j

           i = iLeft
           j = iRight
           do k = iLeft, iEnd-1
               if(i<iRight.and.(j>=iEnd.or.arr(i)>=arr(j))) then
                   arr_tmp(k) = arr(i)
                   indices_tmp(k) = indices(i)
                   i = i + 1
               else
                   arr_tmp(k) = arr(j)
                   indices_tmp(k) = indices(j)
                   j = j + 1
               endif
           enddo
        endsubroutine



!-------------------------------------------------------------------------------------------------------------------------------------
        subroutine basis_transform_vector (mat_basis1_in,basis2_in,mat_basis2_out)
!        Routine to do basis transform by quaternion matrix multiplication(vectors).
!        We follow V' = U * V.
!        V': vectors in new basis;  V: vectors in old basis; U: transformation matrix;        

           !input variables
           real(kind=8), intent(in)               :: mat_basis1_in(:,:,:)    ! vectors in old basis
           real(kind=8), intent(in)               :: basis2_in(:,:,:)        ! transformation matrix
           real(kind=8), intent(inout), pointer   :: mat_basis2_out(:,:,:)   ! vectors in new basis

           integer                   :: LRA, LCA
           integer                   :: LRB, LCB
           integer                   :: LRC, LCC

           LRA = size(basis2_in,1)
           LCA = size(basis2_in,2)
           LRB = size(mat_basis1_in,1)
           LCB = size(mat_basis1_in,2)
           LRC = size(basis2_in,1)
           LCC = size(mat_basis1_in,2)
           allocate(mat_basis2_out(LRC,LCC,4))

           ! We do not use symmetry so NZ is always 4 and the active quaternion units are 1 to 4.
           call qgemm(LRA,LCB,LCA,1.D0,'N','N',(/1,2,3,4/),basis2_in,LRA,LCA,4, & 
                      'N','N',(/1,2,3,4/),mat_basis1_in,LRB,LCB,4,0.D0, &
                      (/1,2,3,4/),mat_basis2_out,LRC,LCC,4)    
                      
        end subroutine basis_transform_vector


!-----------------------------------------------------------------------------------------------------------------------------------
        subroutine basis_transform_matrix (mat_basis1_in,basis2_in,mat_basis2_out)
!        Routine to do basis transform by quaternion matrix multiplication(matrix). 
!        We do similar Transform follows C = A^{-1}BA
!        C: matrix in new basis;  B: matrix in old basis; A: change of basis matrix, A should be a unitary matrix;       

         !input variables
         real(kind=8), pointer, intent(inout)   :: basis2_in(:,:,:)        ! transformation matrix (matrix_A)
         real(kind=8), pointer, intent(in)      :: mat_basis1_in(:,:,:)    ! matrix in old basis(matrix_B)
         real(kind=8), pointer, intent(inout)   :: mat_basis2_out(:,:,:)   ! matrix in new basis (matrix_C)

         integer                                :: LRA, LCA
         integer                                :: LRB, LCB
         integer                                :: LRC, LCC
         integer                                :: iprint

         LRA = size(basis2_in,1)
         LCA = size(basis2_in,2)
         LRB = size(mat_basis1_in,1)
         LCB = size(mat_basis1_in,2)
         LRC = size(basis2_in,2)
         LCC = size(basis2_in,2)
         allocate(mat_basis2_out(LRC,LCC,4))

         ! We do not use symmetry so NZ is always 4 and the active quaternion units are 1 to 4.
         iprint = 0
         call qtrans90('AOMO','S',0.D0,LRB,LCB,LRC,LCC,&
                        mat_basis1_in,LRB,LCB,4,(/1,2,3,4/), &    
                        mat_basis2_out,LRC,LCC,4,(/1,2,3,4/), &
                        basis2_in,LRA,LCA,4,(/1,2,3,4/), & 
                        basis2_in,LRA,LCA,4,(/1,2,3,4/), & 
                        iprint)
         
        end subroutine basis_transform_matrix


!-------------------------------------------------------------------------------------------------------------------------------------                                                   
        subroutine diagonal_matrix_quater(A,EIG,VEC)
! Set parameter needed in qdiag

            real(kind=8), intent(in) :: A(:,:,:)
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
            lwork = 10*LRA
            NZ = 4
            N = LRA
            MATZ = 1
            allocate(EIG(N))
            allocate(VEC(LRV,LCV,NZ))
            allocate(work(MAX(1,LWORK)))

            call qdiag (NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,WORK,LWORK,IERR)

            deallocate(work)


      end subroutine diagonal_matrix_quater
                
               
!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine complex_to_quaternion(dm_mo_complex,dm_mo_quater)
!    refer to routine 'talsh_dm_form_ao_diracformat_from_dm_mo' of talsh_gradient.F90 

      !input variables
      complex(8), pointer, intent(inout)  :: dm_mo_complex(:,:)
      real(8),pointer,     intent(inout)  :: dm_mo_quater(:,:,:)

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
 
      norbt = size(dm_mo_complex,1)/2
      allocate(mo_list(norbt))
      do i = 1, norbt
            mo_list(i) = i
      end do
      
      allocate (dm_mo_quater(norbt,norbt,4))
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
      deallocate (dm_mo_complex)
      end subroutine complex_to_quaternion

!-------------------------------------------------------------------------------------------------------------------------------------     
      subroutine get_Cano_Orbital(orbital_noncano,orbital_cano,Fock_NO,New_orbital_energy)
! This subroutine is used for constructing the canonical orbitals

         !input variables
         real(kind=8), pointer,intent(inout)   :: orbital_noncano(:,:,:)          ! non-canonical orbital in HF mo basis
         real(kind=8), pointer,intent(out)     :: orbital_cano(:,:,:)             ! canonical orbital in HF mo basis
         real(kind=8), pointer,intent(out)     :: New_orbital_energy(:)           ! new orbital energy from qdia of Fock matrix in non-canonical orbital basis
         real(kind=8), pointer,intent(inout)   :: Fock_NO(:,:,:)                  ! v-v block of Fock matrix in non-canonical orbital basis

         real(kind=8),allocatable     :: orbital_cano_basis_noncano(:,:,:)  ! canonical orbital in non-canonical orbital basis
         integer                      :: n, ierr

         ! Diagonalize new FM to obtain the new canonical orbital 
         N = size(Fock_NO,1)
         if (size(Fock_NO,2) /= N) stop " non-square matrix in MP2NO "
         allocate (New_orbital_energy(N))
         allocate(orbital_cano_basis_noncano(N,N,4))
         call qdiag90(4,N,Fock_NO,N,N,New_orbital_energy,1,orbital_cano_basis_noncano,N,N,ierr)

         ! Transform the new canonical orbital from non-canonical orbital basis into HF mo basis
         call basis_transform_vector(orbital_cano_basis_noncano, orbital_noncano, orbital_cano)         
         deallocate (orbital_cano_basis_noncano)

      end subroutine get_Cano_Orbital
 
end module mp2no_matrixtransforms
