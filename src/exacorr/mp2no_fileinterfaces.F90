module mp2no_fileinterfaces
!This module contains the routines needed to do read and write.
!In principle, one could add more routines to read or write the information that they want
!Written by Xiang Yuan.

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_global,only: nfunctions
        use exacorr_utils
        use mp2no_matrixtransforms
        use, intrinsic:: ISO_C_BINDING

        implicit none

        interface quater_write
            module procedure quater_write_orbital
            module procedure quater_write_fock
        end interface

        interface quater_read
            module procedure quater_read_orbital
            module procedure quater_read_matrix
        end interface

        public quater_write, quater_read
        public generate_unified_NOs_MO_file
        public print_rMP2FVNO, write_nos

       contains

!-------------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_write_orbital(orbital_coeff,orbital_eigenval,file_name)
!     routine for writing new orbital coefficient and eigenvalue into the file.

!     input variables
         real(kind=8), pointer              :: orbital_coeff(:,:,:)      ! new orbital coefficient
         real(kind=8), pointer,intent(in)   :: orbital_eigenval(:)       ! new eigenvalue, it could be occupation number or orbital energy
         character(len=*),intent(in)        :: file_name
 
         integer :: sz1, sz2

         integer :: write_orb
         logical :: debug=.true.

         call get_free_fileunit(write_orb)
         sz1 = size(orbital_coeff,1)
         sz2 = size(orbital_coeff,2)
         open (write_orb,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
!           print *,'<<<<< Quaternion orbital, Writing Begin.'
           ! write dimension of the matrix 
           write(write_orb) sz1, sz2
           ! write eigenvalues 
           write(write_orb) orbital_eigenval
           ! write orbitals coeffcient 
           write(write_orb) orbital_coeff
!            print *,'>>>>>  Quaternion orbital, Writing Completed.'
         end if
        close(write_orb,status='KEEP')

      end subroutine quater_write_orbital

      subroutine  quater_write_orbital2(orbital_coeff,orbital_eigenval,file_name,nao,nto,nz)
!     routine for writing new orbital coefficient and eigenvalue into the file.

!     input variables
         integer, intent(in) :: nao,nto,nz
         real(kind=8), dimension(nao,nto,nz), intent(in) :: orbital_coeff      ! new orbital coefficient
         real(kind=8), dimension(nto),intent(in)      :: orbital_eigenval       ! new eigenvalue, it could be occupation number or orbital energy
         character(len=*),intent(in)        :: file_name
 
         integer :: sz1, sz2

         integer :: write_orb
         logical :: debug=.true.

         call get_free_fileunit(write_orb)
         sz1 = size(orbital_coeff,1)
         sz2 = size(orbital_coeff,2)
         open (write_orb,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
!           print *,'<<<<< Quaternion orbital, Writing Begin.'
           ! write dimension of the matrix 
           write(write_orb) sz1, sz2
           print *,"Quater write size = ",sz1,sz2
           ! write eigenvalues 
           write(write_orb) orbital_eigenval
           ! write orbitals coeffcient 
           write(write_orb) orbital_coeff
!            print *,'>>>>>  Quaternion orbital, Writing Completed.'
         end if
        close(write_orb,status='KEEP')

      end subroutine quater_write_orbital2
!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_write_fock(matrix,file_name)
   !     routine for writing fock matrix into the file.

   !     input variables
         real(kind=8), pointer              :: matrix(:,:,:)      ! in principle any quaternion matrix is possible. 
         character(len=*),intent(in)        :: file_name

   !     loop variables
         integer :: sz

         integer :: write_mat
         logical :: debug=.true.

         call get_free_fileunit(write_mat)
         sz = size(matrix,1)
         open (write_mat,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')

         if (debug) then
!           print *,'<<<<< Quaternion Fock Matrix, Writing Begin.'
           ! write dimension of the matrix 
           write(write_mat) sz
           ! write eigenvalues 
           write(write_mat) matrix
!            print *,'>>>>>  Quaternion Fock Matrix, Writing Completed'
         end if
        close(write_mat,status='KEEP')      

      end subroutine quater_write_fock

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_read_orbital(orbital_coeff,orbital_eigenval,Sum_orbital_eigenval_total,threshold,n_orb_trunca,file_name)
!     routine for reading the orbital coefficient and eigenvalue from the file according to the threshold.

!     input variables
         real(kind=8), pointer              :: orbital_coeff(:,:,:)      
         real(kind=8), pointer,intent(out)  :: orbital_eigenval(:)   
         real(kind=8),intent(inout)         :: Sum_orbital_eigenval_total   
         integer,intent(inout)              :: n_orb_trunca  ! number of truncated orbital 
         real(kind=8),intent(in)            :: threshold 
         character(len=*),intent(in)        :: file_name


         real(kind=8), pointer              :: orbital_coeff_all(:,:,:)      
         real(kind=8), pointer              :: orbital_eigenval_all(:)         
   !     loop variables
         integer :: i, j

         !integer :: sz_all ! sz_all -> sz1
         integer :: sz1, sz2, sz
         integer :: read_orb, ierr=0

         logical :: debug=.true.

         call get_free_fileunit(read_orb)
         open (read_orb,file=file_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')

         if (debug) then
!            print *,'<<<<< Start reading all Natural Orbital and Occupation number from the NOs_MO.'
           ! read orbital matrix size
           read(read_orb) sz1, sz2
           print *,"Quater read size = ",sz1,sz2

           allocate(orbital_eigenval_all(sz2))
           allocate(orbital_coeff_all(sz1,sz2,4))
           
           ! read eigenvalue
           read(read_orb) orbital_eigenval_all
           ! read orbitals 
           read(read_orb) orbital_coeff_all

!           print *,'>>>>> All Orbitals, Read completed'
         end if
        close(read_orb,status='KEEP')

        n_orb_trunca = 0
        do i=1,sz2
         if (orbital_eigenval_all(i) >= threshold) then
            n_orb_trunca = n_orb_trunca + 1
         end if
        end do


        sz = n_orb_trunca
        allocate(orbital_eigenval(sz))
        allocate(orbital_coeff(sz1,sz,4))


        do i=1, sz
         orbital_eigenval(i) = orbital_eigenval_all(i)
        end do


        do j=1, sz
            do i=1, sz1
               orbital_coeff(i,j,:) = orbital_coeff_all(i,j,:)
            end do
        end do

        Sum_orbital_eigenval_total = 0.0
        do i=1, sz1
           Sum_orbital_eigenval_total = Sum_orbital_eigenval_total + orbital_eigenval_all(i)
        end do

        deallocate(orbital_coeff_all)
        deallocate(orbital_eigenval_all)
      end subroutine quater_read_orbital

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_read_matrix(NORBT,matrix)
         !     routine for reading the fock matrix from the file.
      
         !     input variables
               real(kind=8), pointer              :: matrix(:,:,:)      ! natural orbital from qdia of DM
               integer,intent(in)                 :: NORBT
               character(8)           :: label(4)
      
               
         !     loop variables
               integer :: sz
      
               integer :: read_mat 
               logical :: debug=.true.
      
               call get_free_fileunit(read_mat)
               open (read_mat,file='CCDENS',status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
               read (read_mat) label
               read (read_mat) matrix
               close(read_mat,status='KEEP')      
      
            end subroutine quater_read_matrix

!--------------------------------------------------------------------------------------------------------------      
      subroutine print_rMP2FVNO(number_energy,number_occupation)
            real(kind=8), pointer,intent(inout)  :: number_occupation(:) ! occupation number from qdia of DM
            real(kind=8), pointer,intent(inout)  :: number_energy(:)
            integer                            :: i, j, k, q

            print *, "  "
            print *,"**********************************************************"
            print *,"***** Recanonized MP2 frozen virtual natural spinors *****"
            print *,"**********************************************************"
            print *, "  "
            print *, "                  ------------------"
            print *, "                  Occupation numbers"
            print *, "                  ------------------"
            print *, "  "
            k = size(number_occupation)/5
            q = mod(size(number_occupation),5)
            do i=1, k, 1
                  j= (i-1)*5+1
                  write(*,"(F20.10,3X,F20.10,3X,F20.10,3X,F20.10,3X,F20.10)") number_occupation(j), number_occupation(j+1), &
                  number_occupation(j+2),number_occupation(j+3), number_occupation(j+4)
!                  print *, ' Orbital',i,'  :Occupation number: ', N_OCCUPATION(i)
            end do

            j = k*5+1
            select case (q)
               case (1)
                  write(*,"(F20.10)") number_occupation(j)

               case (2)
                  write(*,"(F20.10,3X,F20.10)") number_occupation(j), number_occupation(j+1)

               case (3)
                  write(*,"(F20.10,3X,F20.10,3X,F20.10)") number_occupation(j), number_occupation(j+1), number_occupation(j+2)

               case (4)
                  write(*,"(F20.10,3X,F20.10,3X,F20.10,3X,F20.10)") number_occupation(j), number_occupation(j+1), &
                  number_occupation(j+2), number_occupation(j+3)
            end select



            print *, "  "
            print *, "                  ----------------"
            print *, "                  Orbital energies"
            print *, "                  ----------------"
            print *, "  "
            k = size(number_energy)/5
            q = mod(size(number_energy),5)
            do i=1, k, 1
                  j= (i-1)*5+1
                  write(*,"(F20.10,3X,F20.10,3X,F20.10,3X,F20.10,3X,F20.10)") number_energy(j), number_energy(j+1), &
                  number_energy(j+2),number_energy(j+3), number_energy(j+4)
!                  print *, ' Orbital',i,'  :Occupation number: ', N_OCCUPATION(i)
            end do

            j = k*5+1
            select case (q)
               case (1)
                  write(*,"(F20.10)") number_energy(j)

               case (2)
                  write(*,"(F20.10,3X,F20.10)") number_energy(j), number_energy(j+1)

               case (3)
                  write(*,"(F20.10,3X,F20.10,3X,F20.10)") number_energy(j), number_energy(j+1), number_energy(j+2)

               case (4)
                  write(*,"(F20.10,3X,F20.10,3X,F20.10,3X,F20.10)") number_energy(j), number_energy(j+1), number_energy(j+2), &
                                                            number_energy(j+3)
            end select

            print *,"**********************************************************"
            print *, " "


      end subroutine 

      subroutine write_nos (nocc,NOs_AO,New_orbital_energy,reference_orbitals,new_orbitals)
     !      write new orbital information into the appropriate section of checkpoint 

            implicit none

            character(len=*),intent(in)           :: reference_orbitals,new_orbitals
            real(8),pointer,intent(in)            :: NOs_AO(:,:,:)          ! new canonical orbital in AO basis
            real(kind=8),pointer,intent(in)       :: New_orbital_energy(:)  ! new orbital energy 

            type(basis_set_info_t)   :: ao_basis
            type(cmo)                :: qorbital        ! the new canonical orbital 
            type(cmo)                :: qorbital_hf     ! the original HF canonical orbital 

            integer                  :: nocc            ! occupied orbital 
            integer                  :: nao, nmo, nno   ! number of atomic orbital; molecular orbital; natural orbital
            integer                  :: i,j,iz
            integer                  :: ierr
            integer                  :: nmo_hf


   !        read the reference orbitals
            call read_from_dirac(ao_basis,qorbital_hf,reference_orbitals,ierr)
            call convert_mo_to_quaternion(qorbital_hf)

            nao = qorbital_hf%nao
            nno = size(NOs_AO,2)
            nmo = nocc+nno
            nmo_hf = qorbital_hf%nmo

            ! make a copy of the reference (Hartree-Fock) orbitals
            call copy_mo(qorbital_hf,qorbital)

            ! replace the orbital energies by the new (recanonized) ones
            do i = nocc+1, nmo
               qorbital%energy(i)=New_orbital_energy(i-nocc)
            end do

            ! we did not recanonize the deleted virtuals, give them a sufficiently high energy to not interfere
            if (nmo < nmo_hf) then
                  do i=nmo+1, nmo_hf
                     qorbital%energy(i)=qorbital_hf%energy(nmo_hf)*100
                  end do
            end if

            ! replace the coefficients by the recanonized natural orbital selection
            do iz = 1, 4
               do  j = nocc+1, nmo
                  do  i = 1, qorbital%nao
                     qorbital%coeff_q(i,j,iz)=NOs_AO(i,j-nocc,iz)
                  enddo
               enddo
            enddo

            call write_mos_to_checkpoint(2,qorbital,new_orbitals)
            call print_date('Written natural orbitals to '//trim(new_orbitals)//' section of checkpoint')

      end subroutine write_nos

      subroutine generate_unified_NOs_MO_file(exa_input)
         ! Routine specific to MP2FNO++ to regroup into the same NOs_MO file
         ! the NOs_MO generated for different perturbation.
         ! Also print some information about the FVO size reduction.
         implicit none

         type(exacc_input), intent(in) :: exa_input

         real(kind=8), allocatable  :: NOs_MO(:,:,:), NOs_MO_tmp(:,:,:)
         real(kind=8), allocatable  :: N_OCCUPATION(:), N_OCCUPATION_tmp(:)
         real(kind=8)               :: N_OCCUPATION_tmp_j
         real(kind=8)               :: t_mp2no, Sum_N_OCCUPATION_total
         character(len=10)          :: fn
         integer :: unit
         integer :: nao, nto, nz
         integer :: i, j, iMO
         

         t_mp2no = exa_input%t_mp2no
         Sum_N_OCCUPATION_total = 0.0

         ! read nao from checkpoint file
         nz = 4 
         nto = 0

         ! first read all files to get nto and nao
         call print_date('First reading every NOs_MO file...')
         do i=1,exa_input%n_mp2fno_prop
            write(fn, '(A,I0)') 'NOs_MO', i

            call get_free_fileunit(unit)
            open (unit,file=fn,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')

            read(unit) nao, j ! j is a dummy variable

            allocate(N_OCCUPATION_tmp(nao))
            read(unit) N_OCCUPATION_tmp

            close(unit,status='KEEP')

            ! number of orbitals whose occupation number is high enough
            do j =1,nao
               N_OCCUPATION_tmp_j = N_OCCUPATION_tmp(j)
               Sum_N_OCCUPATION_total = Sum_N_OCCUPATION_total + N_OCCUPATION_tmp_j
               if(N_OCCUPATION_tmp_j.ge.t_mp2no) then
                  nto = nto + 1
               endif
            enddo

            deallocate(N_OCCUPATION_tmp)
         enddo

         write(*,*) "nto = ",nto

         ! allocate arrays
         allocate(N_OCCUPATION(nto))
         allocate(NOs_MO(nao,nto,nz))
         allocate(N_OCCUPATION_tmp(nao))
         allocate(NOs_MO_tmp(nao,nao,nz))

         ! read all coefficients + occupation numbers
         call print_date('Second read of NOs_MO files...')
         iMO = 1
         do i=1,exa_input%n_mp2fno_prop
            write(fn, '(A,I0)') 'NOs_MO', i

            call get_free_fileunit(unit)
            open (unit,file=fn,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')

            read(unit) nao, j ! j is a dummy variable
            read(unit) N_OCCUPATION_tmp
            read(unit) NOs_MO_tmp

            close(unit,status='KEEP')

            do j = 1,nao
               if(N_OCCUPATION_tmp(j).ge.t_mp2no) then
                  NOs_MO(:,iMO,:) = NOs_MO_tmp(:,j,:)
                  N_OCCUPATION(iMO) = N_OCCUPATION_TMP(j)
                  iMO = iMO + 1
               endif
            enddo

         enddo
         write(*,*) "iMO = ",iMO
         call print_date('Finished!')
         
         !deallocate(N_OCCUPATION_tmp)
         !deallocate(NOs_MO_tmp)

         ! sort by occupation number
         call sort_MOs_by_occupation_number(N_OCCUPATION,NOs_MO,nao,nto,nz)
         !do i=1,size(N_OCCUPATION)
         !    write(*,*) "DEBUG Nocc", i,"  = ",N_OCCUPATION(i) 
         !enddo

         ! save into the final NOs_MO file
         call quater_write_orbital2(NOs_MO,N_OCCUPATION,'NOs_MO',nao,nto,nz)
         call print_date('MP2FNO++: Wrote all natural orbital and occupation number in file NOs_MO.')

         ! print info
         call print_mp2fvno_occupation_number_details(t_mp2no, nto, N_OCCUPATION, Sum_N_OCCUPATION_total)

         !deallocate(N_OCCUPATION)
         !deallocate(NOs_MO)
      end subroutine generate_unified_NOs_MO_file

end module mp2no_fileinterfaces

