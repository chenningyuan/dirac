module mp2no_fileinterfaces
!This module contains the routines needed to do read and write.
!In principle, one could add more routines to read or write the information that they want
!Written by Xiang Yuan.

        use exacorr_datatypes
        use exacorr_mo
        use exacorr_global,only: nfunctions
        use exacorr_utils
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

        public quater_write_mos_to_dirac,quater_read_from_dirac_fock, &
               quater_read_from_dirac_orbital, quater_write, quater_read

        public print_rMP2FVNO, write_to_newfile

       contains

!-------------------------------------------------------------------------------------------------------------------------------------------
         subroutine quater_write_mos_to_dirac (ao_basis,qorbital,file_out_name)
     !      write mo coefficients in dirac format
     !      refer to subroutine write_mos_to_dirac of exacorr_dirac.F

     
            implicit none
     
            character(len=*),intent(in)           :: file_out_name
            type(basis_set_info_t), intent(inout) :: ao_basis
            type(cmo), intent(inout)              :: qorbital

            integer                :: ierr
            character              :: text*74      ! holds text read from DFCOEF, only used for printing this
            integer                :: nsym
            real(8)                :: total_energy
            logical                :: tobe
            integer                :: lucoef 
            integer                :: ishell
            integer                :: nshells, iao, ibas
            integer                :: npo, nto, npriexp, ncnt, nzbuf
            integer                :: nao, nmo, nno ! number of atomic orbital; molecular orbital; natural orbital
            integer, allocatable   :: numcf(:),jstrt(:)
            integer,allocatable    :: index(:)
            real(8), allocatable   :: priexp(:),priccf(:)
            integer                :: i,j
            integer                :: mp2dm
            character(8)           :: label(4)
            character(len=10) :: datex
            character(len=8)  :: timex
     
   
     !      switch to Kramers-restricted format
            text = "MO coefficients produced by geniao "
            npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
            nsym = 1 ! we do not use symmetry, so nsym is always 1
            nzbuf = 4
   

   !        print warning in case we overwrite a file and open the file
            call get_free_fileunit(lucoef)
            inquire (file=file_out_name,exist=tobe)
            if (tobe) then
               print *, "warning: overwriting file "//file_out_name
               open (lucoef,file=file_out_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
            else
               open (lucoef,file=file_out_name,status='NEW',FORM='UNFORMATTED',access='SEQUENTIAL')
            end if
     
            ! convert basis set info to 1-d lists used by DIRAC (assuming uncontracted basis sets)
            ! TODO: change storage on DIRAC side to use the basis set type
   
            ! we always have uncontracted basis sets, so can hardwire the next variable and equate the number of exponents to the number of shells
            nshells = ao_basis%nshells
            ncnt = 1
            npriexp =  nshells
            allocate(jstrt(nshells))
            allocate(numcf(nshells))
            allocate(priexp(nshells))
            allocate(priccf(nshells))
            do ishell = 1, nshells
               jstrt(ishell) = ishell-1
               numcf(ishell) = 1
               priexp(ishell) = ao_basis%gtos(ishell)%exponent(1)
               priccf(ishell) = ao_basis%gtos(ishell)%coefficient(1)
            end do
      
   
            ! write the header
            label(1) = '********'
            call daytime(datex,timex)
            label(2) = datex
            label(3) = timex
            label(4) = 'INFO    '
            write (lucoef) label
            write (lucoef) text,nsym,nzbuf,npo,qorbital%nmo,qorbital%nao,qorbital%total_energy

   
            ! write the mo-coefficients
            call daytime(datex,timex)
            label(2) = datex
            label(3) = timex
            label(4) = 'COEFS   '
            write (lucoef) label
            write (lucoef) qorbital%coeff_q


            ! write the orbital energies
            call daytime(datex,timex)
            label(2) = datex
            label(3) = timex
            label(4) = 'EVALS   '           
            write (lucoef) label
            write (lucoef) qorbital%energy


            ! write the symnmetry information
            call daytime(datex,timex)
            label(2) = datex
            label(3) = timex
            label(4) = 'SUPERSYM' 
            write (lucoef) label
            write (lucoef) qorbital%symmetry


            ! write the ao basis information header
            call daytime(datex,timex)
            label(2) = datex
            label(3) = timex
            label(4) = 'BASINFO ' 
            write (lucoef) label
            write (lucoef) ao_basis%nao,ao_basis%nshells,npriexp,ncnt 



     !      orb_momentum, atom_number, n_primives for all shells
            write (lucoef) ao_basis%gtos(:)%orb_momentum,ao_basis%gtos(:)%atom_number,ao_basis%gtos(:)%n_primitives
     !      start of exponent values, number of contracted for all shells
            write (lucoef) jstrt,numcf
     !       write (lucoef) (ao_basis%gtos(ishell)%coord(:),ishell=1,nshells)
            write (lucoef) (ao_basis%gtos(:)%coord(i),i=1,3)
     !      values of the exponents for all shells (packed into a 1-d array)
            write (lucoef) priexp
     !      values of the contraction coefficients for all shells (packed into a 1-d array)
            write (lucoef) priccf
   
            deallocate(jstrt)
            deallocate(numcf)
            deallocate(priexp)
            deallocate(priccf)
           
            close (lucoef,status='KEEP')   
            return
   
           end subroutine quater_write_mos_to_dirac
   

!----------------------------------------------------------------------------------------------------------------------
      subroutine quater_read_from_dirac_fock(Fock_HF,sz,k,flb)
!     routine for only reading orbital energy from DFCOEF and construting the diagonal Fock matrix 
!     flb =1 we construt the oo block otherwise we build vv block. 

            real(kind=8), pointer, intent(inout)    :: Fock_HF(:,:,:)          ! Fock matrix 
            integer, intent(inout)                  :: sz,flb,k                ! size of matrix and flag for block
            type(cmo)                               :: qorbital
            type(basis_set_info_t)                  :: ao_basis
            integer                                 :: ierr
            integer                                 :: i,sz_total_energy


!     read information from DFCOEF
            call read_from_dirac(ao_basis,qorbital,'DFCOEF',ierr)
            call convert_mo_to_quaternion(qorbital)
            sz_total_energy = size(qorbital%energy,1)

            allocate(Fock_HF(sz,sz,4))

!     construct the FM 
            Fock_HF = 0.0
            if (flb == 1) then
               do i = 1, sz
                  Fock_HF(i,i,1) = qorbital%energy(i)
               end do
            else
   !            k = sz_total_energy - sz
               do i = 1, sz
                  Fock_HF(i,i,1) = qorbital%energy(k+i)
               end do
            end if 
         
      end subroutine quater_read_from_dirac_fock


!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine quater_read_from_dirac_orbital(HF_AO,sz,k,flb)
!     routine for only reading HF orbital coefficient from DFCOEF.
!     flb =1 we read occupied orbital otherwise we read virtual orbital 
 

         real(kind=8), pointer, intent(inout)    :: HF_AO(:,:,:)            ! new canonical orbital 
         integer, intent(inout)                  :: sz,flb,k                ! size of matrix and flag for block
         type(cmo)                               :: qorbital
         type(basis_set_info_t)                  :: ao_basis
         integer                                 :: ierr
         integer                                 :: i, j, sz_total_mo, sz_total_ao


!     read information of from DFCOEF
         call read_from_dirac(ao_basis,qorbital,'DFCOEF',ierr)
         call convert_mo_to_quaternion(qorbital)
         sz_total_ao = size(qorbital%coeff_q,1)
         sz_total_mo = size(qorbital%coeff_q,2)

         allocate(HF_AO(sz_total_ao,sz,4))
      
!     constrcut HF orbital matrix
         if (flb == 1) then
            do j = 1, sz
               do  i = 1, sz_total_ao
                  HF_AO(i,j,1) = qorbital%coeff_q(i,j,1)
                  HF_AO(i,j,2) = qorbital%coeff_q(i,j,2)
                  HF_AO(i,j,3) = qorbital%coeff_q(i,j,3)
                  HF_AO(i,j,4) = qorbital%coeff_q(i,j,4)
               end do
            end do
         else
!            k = sz_total_mo - sz
            do j = 1, sz
               do i = 1, sz_total_ao
                  HF_AO(i,j,1) = qorbital%coeff_q(i,k+j,1)
                  HF_AO(i,j,2) = qorbital%coeff_q(i,k+j,2)
                  HF_AO(i,j,3) = qorbital%coeff_q(i,k+j,3)
                  HF_AO(i,j,4) = qorbital%coeff_q(i,k+j,4)
               end do
            end do
         end if 

      end subroutine quater_read_from_dirac_orbital

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_write_orbital(orbital_coeff,orbital_eigenval,file_name)
!     routine for writing new orbital coefficient and eigenvalue into the file.

!     input variables
         real(kind=8), pointer              :: orbital_coeff(:,:,:)      ! new orbital coefficient
         real(kind=8), pointer,intent(in)   :: orbital_eigenval(:)       ! new eigenvalue, it could be occupation number or orbital energy
         character(len=*),intent(in)        :: file_name
 
         integer :: sz

         integer :: write_orb
         logical :: debug=.true.

         call get_free_fileunit(write_orb)
         sz = size(orbital_coeff,2)
         open (write_orb,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
!           print *,'<<<<< Quaternion orbital, Writing Begin.'
           ! write dimension of the matrix 
           write(write_orb) sz
           ! write eigenvalues 
           write(write_orb) orbital_eigenval
           ! write orbitals coeffcient 
           write(write_orb) orbital_coeff
!            print *,'>>>>>  Quaternion orbital, Writing Completed.'
         end if
        close(write_orb,status='KEEP')

      end subroutine quater_write_orbital
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

         integer :: sz_all, sz
         integer :: read_orb, ierr=0

         logical :: debug=.true.

         call get_free_fileunit(read_orb)
         open (read_orb,file=file_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')

         if (debug) then
!            print *,'<<<<< Start reading all Natural Orbital and Occupation number from the NOs_MO.'
           ! read orbital matrix size
           read(read_orb) sz_all

           allocate(orbital_eigenval_all(sz_all))
           allocate(orbital_coeff_all(sz_all,sz_all,4))
           
           ! read eigenvalue
           read(read_orb) orbital_eigenval_all
           ! read orbitals 
           read(read_orb) orbital_coeff_all

!           print *,'>>>>> All Orbitals, Read completed'
         end if
        close(read_orb,status='KEEP')

        n_orb_trunca = 0
        do i=1,sz_all
         if (orbital_eigenval_all(i) >= threshold) then
            n_orb_trunca = n_orb_trunca + 1
         end if
        end do


        sz = n_orb_trunca
        allocate(orbital_eigenval(sz))
        allocate(orbital_coeff(sz_all,sz,4))


        do i=1, sz
         orbital_eigenval(i) = orbital_eigenval_all(i)
        end do


        do j=1, sz
            do i=1, sz_all
               orbital_coeff(i,j,:) = orbital_coeff_all(i,j,:)
            end do
        end do

        Sum_orbital_eigenval_total = 0.0
        do i=1, sz_all
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

      subroutine write_to_newfile (nocc,NOs_AO,N_OCCUPATION,New_orbital_energy,file_old,file_new)
     !      write new orbital information into new file 'MP2NOs_AO', which is similar to DFCOEF 
     !      copy the basis information from DFCOEF and use new orbital energy and coefficient to replace the old orbital.   
     !      revise subroutine write_mos_to_dirac of exacorr_dirac.F

            implicit none

            character(len=*),intent(in)           :: file_old,file_new
            real(8),pointer,intent(in)            :: NOs_AO(:,:,:)          ! new canonical orbital in AO basis
            real(kind=8),pointer,intent(in)       :: N_OCCUPATION(:)        ! occupation number
            real(kind=8),pointer,intent(in)       :: New_orbital_energy(:)  ! new orbital energy 


            type(basis_set_info_t)   :: ao_basis
            type(cmo)                :: qorbital        ! save the new canonical orbital 
            type(cmo)                :: qorbital_hf     ! save the original HF canonical orbital 

            integer                  :: nocc            ! occupied orbital 
            integer                  :: nao, nmo, nno   ! number of atomic orbital; molecular orbital; natural orbital
            integer                  :: i,j,iz
            integer                  :: ierr
            integer                  :: nmo_hf


   !        read information of from DFCOEF
            call read_from_dirac(ao_basis,qorbital_hf,file_old,ierr)
            call convert_mo_to_quaternion(qorbital_hf)

            nao = qorbital_hf%nao
            nno = size(NOs_AO,2)
            nmo = nocc+nno
            nmo_hf = qorbital_hf%nmo

            ! make a copy of the Hartree-Fock orbitals
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

            call quater_write_mos_to_dirac(ao_basis,qorbital,file_new)
            call print_date('Write new orbital into MP2NOs_AO file like DFCOEF')

      end subroutine write_to_newfile

end module mp2no_fileinterfaces

