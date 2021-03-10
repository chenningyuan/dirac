module write_orbital_to_file
!This module contains the routines needed to read and write.
!In principle, people could add more routines to read or write the information that they need
!At the moment, there are three routines:
! 1) read all information from DFCOEF
! 2) write new energy and orbital into the new file like DFCOEF.
! 3) only read HF energy and orbital coefficient from DFCOEF

        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use exacorr_mo
        use exacorr_global,only: nfunctions
        use, intrinsic:: ISO_C_BINDING

        implicit none

        interface quater_write
            module procedure quater_write_orbital
            module procedure quater_write_fock
        end interface

        interface quater_read
            module procedure quater_read_orbital
            module procedure quater_read_fock
        end interface

        public quater_read_from_dirac_all,quater_write_mos_to_dirac,quater_read_from_dirac_fock, &
               quater_read_from_dirac_orbital, quater_write, quater_read

       contains

!-------------------------------------------------------------------------------------------------------------------------------------------

        subroutine quater_read_from_dirac_all (ao_basis,qorbital,file_name,ierr)

         use exacorr_mo
         use exacorr_datatypes
         use exacorr_global,only: nfunctions
         use exacorr_utils
  
         implicit none

  !      we need information of ao_basis
  !      read all positive energy mo coefficients present on DFCOEF file
  !      check for imcompatibiklities (use of symmetry) as well and continue
  !      by reading also the basis set information into gtos
  !      this subroutine similar to routine read_from_dirac of exacorr_dirac.F 
  
         character(len=*),intent(in) :: file_name
         type(basis_set_info_t), intent(inout) :: ao_basis
         type(qumo), intent(inout):: qorbital
         integer, intent(out)   :: ierr

         integer, parameter     :: NOERROR  = 0 ! no error
         integer, parameter     :: NOTFOUND = 1 ! file not found
         integer, parameter     :: CORRUPT  = 2 ! file unreadable
         integer, parameter     :: SYMMETRY = 3 ! calculation used symmetry
         integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old DFCOEF?)
         integer, parameter     :: CONTRAC  = 5 ! basis set is contracted (Interest can not handle this)
         integer, parameter     :: INCONSIS = 6 ! DFCOEF does not contain consistent information wrt number of basis functions
         character              :: text*74      ! holds text read from DFCOEF, only used for printing this
         integer                :: nsym
         integer                :: idim(3,2)
         real(8)                :: total_energy
         logical                :: tobe
         integer                :: lucoef
         integer                :: i, j, k, ishell, iprimstart, iprimend
         integer                :: nao, nao_bas, nshells, iao, ibas
         integer                :: nmo, npo, nto, npriexp, ncnt
         integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
         integer, allocatable   :: jstrt(:), numcf(:)
         real(8), allocatable   :: coeff(:,:,:),eig(:)
         integer, allocatable   :: ibeig(:),index(:)
         integer, allocatable   :: ncent(:),nhkt(:),nuco(:)
         real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)
  
         ierr = NOERROR
  
         ! open file and read dimension info
         inquire (file=file_name,exist=tobe)
         if (.not.tobe) then
            print*, "error: MO coefficient file"//file_name//" not found"
            ierr = NOTFOUND
            return
         else
            call get_free_fileunit(lucoef)
            open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
            read (lucoef,end=10,err=10) text,nsym,((idim(i,j),i=1,3),j=1,nsym),total_energy
         end if
  
         ! check symmetry and get dimensions of the coefficient array
         if (nsym > 1) then
            ierr = SYMMETRY
            close (lucoef,status='KEEP')
            return
         else
            nao = idim(3,1) ! number of basis functions
            npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
            nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
            nto = npo+nmo   ! total number of solutions on this file
            if (npo > 0) print*, "warning: reading from 4C DFCOEF file" ! this warning will go once we have the full index transform implemented
         end if
  
!         print *,'print nto', nto
         allocate(coeff(nao,nto,nz))
         allocate(eig(nto))
         allocate(ibeig(nto))
         allocate(index(nmo))
               
         ! read the mo coefficients and orbital eigenvalues
         read (lucoef,end=10,err=10) coeff
         read (lucoef,end=10,err=10) eig
         read (lucoef,end=10,err=10) ibeig
         ! fill in index values
         do i = 1, nmo
            index(i) = i ! we set the index compatible with the fact that we store all mo's
         end do

  
         ! put in qumo
         call alloc_mo(qorbital,nao,nmo,4)
         qorbital%total_energy = total_energy
         qorbital%index        = index
         qorbital%coeff        = coeff
         qorbital%energy       = eig
         qorbital%boson_irrep  = 0 ! because nosym. if there is, include ibeig in converter
  
         deallocate(coeff)
         deallocate(eig)
         deallocate(ibeig)
         deallocate(index)


         ! read the basis set information and check for consistency
         read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
         if (nao < nao_bas) then
            ierr = INCONSIS
            close (lucoef,status='KEEP')
            print *,"ierr", ierr     
            return
         elseif (nao > nao_bas) then
            print*,"WARNING: assuming Levy-Leblond Hamiltonian"
            print*,"Reading fewer basis functions than MO-coefficients"
         end if

         
   !     allocate the arrays to read the basis set information in DIRAC format
         allocate(ncent(nshells))
         allocate(nhkt(nshells))
         allocate(nuco(nshells))
         allocate(numcf(nshells))
         allocate(jstrt(nshells))
         allocate(cent(nshells,3))
         allocate(priexp(npriexp))
         allocate(priccf(npriexp,ncnt))
  
  !      orb_momentum, atom_number, n_primives for all shells
         read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
  !      start of exponent values, number of contracted for all shells
         read (LUCOEF,err=20,end=20) jstrt,numcf
  !      coordinates (x, y, z) of center for all shells
         read (LUCOEF,err=20,end=20) cent
  !      values of the exponents for all shells (packed into a 1-d array)
         read (LUCOEF,err=20,end=20) priexp
  !      values of the contraction coefficients for all shells (packed into a 1-d array)
         read (LUCOEF,err=20,end=20) priccf
  
         ! put in gtos
         ao_basis%nshells       = nshells
         ao_basis%nao           = nao_bas
         ao_basis%basis_angular = 1
         nullify(ao_basis%gtos)
         allocate(ao_basis%gtos(nshells))
         do ishell = 1, nshells
            ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
            ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
            ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
            nullify  (ao_basis%gtos(ishell)%exponent) 
            allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
            nullify  (ao_basis%gtos(ishell)%coefficient)
            allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
            iPrimStart = jstrt(ishell) + 1
            iPrimEnd   = jstrt(ishell) + nuco(ishell)
            k          = numcf(ishell)
            ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
            ao_basis%gtos(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
            ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
  !         In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
           call compress_unc_shell (ao_basis%gtos(ishell),ierr)
            if (ierr == 1) ierr = CONTRAC            
         end do

         deallocate(ncent)
         deallocate(nhkt)
         deallocate(nuco)
         deallocate(numcf)
         deallocate(jstrt)
         deallocate(cent)
         deallocate(priexp)
         deallocate(priccf)
  
   !      Make array with pointers from a basis function index to its shell
         nullify(ao_basis%shell_indices)
         allocate(ao_basis%shell_indices(nao))
         iao = 0
         do ishell = 1, nshells
            do ibas = 1, &
                   nfunctions(ao_basis%gtos(ishell)%orb_momentum,ao_basis%basis_angular)
               iao = iao + 1
               if (iao.gt.nao) then
                  ierr = INCONSIS
                  call quit ('Error in read_from_dirac')
               endif
               ao_basis%shell_indices(iao) = ishell
            end do
         end do
         if (iao.gt.nao) then
            ierr = INCONSIS
            close (lucoef,status='KEEP')
            call quit ('Error in read_from_dirac')
         endif
  
         close (lucoef,status='KEEP')
  
         call print_date("read from DFCOEF:"//text)
         return
  
    10   ierr = CORRUPT
         close (lucoef,status='KEEP')
         return
  
    20   ierr = NOBASIS
         close (lucoef,status='KEEP')
         return
  
         end subroutine quater_read_from_dirac_all

                
!-----------------------------------------------------------------------------------------------------------------------------------------------

         subroutine quater_write_mos_to_dirac (ao_basis,qorbital,file_out_name)

            use exacorr_mo
            use exacorr_datatypes
            use exacorr_utils
     
            implicit none
     
     !      write mo coefficients in dirac format
     !      revise subroutine write_mos_to_dirac from exacorr_dirac.F
     
            character(len=*),intent(in)           :: file_out_name
            type(basis_set_info_t), intent(inout) :: ao_basis
            type(qumo), intent(inout)             :: qorbital

            integer                :: ierr
            character              :: text*74      ! holds text read from DFCOEF, only used for printing this
            integer                :: nsym
            real(8)                :: total_energy
            logical                :: tobe
            integer                :: lucoef 
            integer                :: ishell
            integer                :: nshells, iao, ibas
            integer                :: npo, nto, npriexp, ncnt
            integer                :: nao, nmo, nno ! number of atomic orbital; molecular orbital; natural orbital
            integer, allocatable   :: numcf(:),jstrt(:)
            integer,allocatable    :: index(:)
            real(8), allocatable   :: priexp(:),priccf(:)
            integer                :: i,j
            integer                :: mp2dm
     
   
     !      switch to Kramers-restricted format
            text = "MO coefficients produced by geniao "
            npo = 0  ! number of "positronic orbitals is hardwired to 0 for now
            nsym = 1 ! we do not use symmetry, so nsym is always 1
   

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
            write (lucoef) text,nsym,npo,qorbital%nmo,qorbital%nao,qorbital%total_energy
            ! write the mo-coefficients
            write (lucoef) qorbital%coeff
            ! write the orbital energies
            write (lucoef) qorbital%energy
            ! write the boson symnmetry information
            write (lucoef) qorbital%boson_irrep
            ! write the ao basis information header
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
   
   
!            open (mp2dm,file='Orbital_check',status='REPLACE',FORM='FORMATTED',access='SEQUENTIAL')
!            Print *, "write HF orbital cofficient into file Orbital_check."
!              do j=1,size(qorbital%coeff,2)
!                 do i=1,size(qorbital%coeff,1)
!                    write(mp2dm,'(2x,4F28.10,2x,4I4)') qorbital%coeff(i,j,:),i,j
!                 end do
!               end do
!           close(mp2dm,status='KEEP')


            deallocate(jstrt)
            deallocate(numcf)
            deallocate(priexp)
            deallocate(priccf)

            
            close (lucoef,status='KEEP')
   
            return
   
           end subroutine quater_write_mos_to_dirac
   
   
   
   !---------------------------------------------------------------------------------------------------------------------------------------

        subroutine quater_read_from_dirac_fock(Fock_HF,sz,flb)
         use exacorr_mo
         use exacorr_datatypes
         use exacorr_global,only: nfunctions
         use exacorr_utils
      !     routine for reading orbital energy and construting the diagonal Fock matrix 

            real(kind=8), pointer, intent(inout)    :: Fock_HF(:,:,:)          ! Fock matrix 
            integer, intent(inout)                  :: sz,flb                  ! size of matrix and flag for block
            type(qumo)                              :: qorbital
            type(basis_set_info_t)                  :: ao_basis
            integer                                 :: ierr
            integer                                 :: i, k, sz_total_energy


!           read information from DFCOEF
            call quater_read_from_dirac_all(ao_basis,qorbital,'DFCOEF',ierr)
            sz_total_energy = size(qorbital%energy,1)

            allocate(Fock_HF(sz,sz,4))

!           construct the FM 
            Fock_HF = 0.0

            if (flb == 1) then
               do i = 1, sz
                  Fock_HF(i,i,1) = qorbital%energy(i)
               end do
            else
               k = sz_total_energy - sz
               do i = 1, sz
                  Fock_HF(i,i,1) = qorbital%energy(k+i)
               end do
            end if 
         
        end subroutine quater_read_from_dirac_fock


!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine quater_read_from_dirac_orbital(HF_AO,sz,flb)
         use exacorr_mo
         use exacorr_datatypes
         use exacorr_global,only: nfunctions
         use exacorr_utils

         real(kind=8), pointer, intent(inout)    :: HF_AO(:,:,:)            ! new canonical orbital 
         integer, intent(inout)                  :: sz,flb                  ! size of matrix and flag for block
         type(qumo)                              :: qorbital
         type(basis_set_info_t)                  :: ao_basis
         integer                                 :: ierr
         integer                                 :: i, j, k, sz_total_mo, sz_total_ao


!           read information of from DFCOEF
         call quater_read_from_dirac_all(ao_basis,qorbital,'DFCOEF',ierr)
         sz_total_ao = size(qorbital%coeff,1)
         sz_total_mo = size(qorbital%coeff,2)

         allocate(HF_AO(sz_total_ao,sz,4))
      

!          constrcut HF orbital matrix
         if (flb == 1) then
            do j = 1, sz
               do  i = 1, sz_total_ao
                  HF_AO(i,j,1) = qorbital%coeff(i,j,1)
                  HF_AO(i,j,2) = qorbital%coeff(i,j,2)
                  HF_AO(i,j,3) = qorbital%coeff(i,j,3)
                  HF_AO(i,j,4) = qorbital%coeff(i,j,4)
               end do
            end do
         else
            k = sz_total_mo - sz
            do j = 1, sz
               do i = 1, sz_total_ao
                  HF_AO(i,j,1) = qorbital%coeff(i,k+j,1)
                  HF_AO(i,j,2) = qorbital%coeff(i,k+j,2)
                  HF_AO(i,j,3) = qorbital%coeff(i,k+j,3)
                  HF_AO(i,j,4) = qorbital%coeff(i,k+j,4)
               end do
            end do
         end if 


      end subroutine quater_read_from_dirac_orbital

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_write_orbital(orbital,eigenval,file_name)
   !     input variables
         real(kind=8), pointer              :: orbital(:,:,:)      ! natural orbital from qdia of DM
         real(kind=8), pointer,intent(out)  :: eigenval(:)         ! occupation number from qdia of DM
         character(len=*),intent(in)        :: file_name
 

   !     loop variables
         integer :: i, j, k, sz

         integer, parameter :: write_orb = 10
         logical :: debug=.true.

         sz = size(orbital,1)
         open (write_orb,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
           print *,'<<<<< Quaternion orbital, Writing Begin.'
           ! write dimension of the matrix 
           write(write_orb) sz
           ! write eigenvalues 
           write(write_orb) eigenval
           ! write orbitals 
           write(write_orb) orbital

!           k = 1
!           do j=1,sz
!              do i=1,sz
!               if(j == i) then
!                  write(write_orb,'(1X,I8,2x,4F28.16,2x,2I4,2x,F28.16)') k, orbital(i,j,:),i,j,eigenval(i)
!                  k = k + 1
!               else 
!                  write(write_orb,'(1X,I8,2x,4F28.16,2x,2I4)') k, orbital(i,j,:),i,j
!                  k = k + 1
!               end if 
!              end do
!            end do

            print *,'>>>>>  Quaternion orbital, Writing Completed.'
         end if
        close(write_orb,status='KEEP')

      end subroutine quater_write_orbital
!-------------------------------------------------------------------------------------------------------------------------------------

      subroutine  quater_write_fock(matrix,file_name)
   !     input variables
         real(kind=8), pointer              :: matrix(:,:,:)      ! natural orbital from qdia of DM
         character(len=*),intent(in)        :: file_name

   !     loop variables
         integer :: i, j, k, sz

         integer, parameter :: write_mat = 10
         logical :: debug=.true.

         sz = size(matrix,1)
         open (write_mat,file=file_name,status='REPLACE',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
           print *,'<<<<< Quaternion Fock Matrix, Writing Begin.'
           ! write dimension of the matrix 
           write(write_mat) sz
           ! write eigenvalues 
           write(write_mat) matrix

!           k = 1
!           do j=1,sz
!              do i=1,sz
!                  write(write_mat,'(1X,I8,2x,4F28.16,2x,2I4)') k,matrix(i,j,:),i,j
!                  k = k + 1
!              end do
!            end do

            print *,'>>>>>  Quaternion Fock Matrix, Writing Completed'
         end if
        close(write_mat,status='KEEP')      

      end subroutine quater_write_fock

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_read_orbital(orbital,eigenval,threshold,n_orb,file_name)

   !     input variables
         real(kind=8), pointer              :: orbital(:,:,:)      ! natural orbital from qdia of DM
         real(kind=8), pointer,intent(out)  :: eigenval(:)         ! occupation number from qdia of DM
         integer,intent(inout)              :: n_orb
         real(kind=8),intent(in)            :: threshold 
         character(len=*),intent(in)        :: file_name


         real(kind=8), pointer              :: orbital_all(:,:,:)      
         real(kind=8), pointer              :: eigenval_all(:)         
   !     loop variables
         integer :: i, j, k, sz_all, sz, p, q, ierr=0
         real(kind=8) :: k_re

         integer, parameter :: read_orb = 10
         logical :: debug=.true.

         open (read_orb,file=file_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')

         if (debug) then
            print *,'<<<<< Start reading all Natural Orbital and Occupation number from the NOs_MO.'
           ! read orbital matrix size
           read(read_orb) sz_all

           allocate(orbital_all(sz_all,sz_all,4))
           allocate(eigenval_all(sz_all))

           ! read eigenvalue
           read(read_orb) eigenval_all
           ! read orbitals 
           read(read_orb) orbital_all

!           do while(ierr==0)
!            read(read_orb,'(3X,I8)',iostat=ierr) k
!           end do

           print *,'>>>>> All Orbitals, Read completed'
         end if
        close(read_orb,status='KEEP')

!         k_re = real(k)
!         sz_all = int(sqrt(k_re))

!         allocate(orbital_all(sz_all,sz_all,4))
!         allocate(eigenval_all(sz_all))
!         open (read_orb,file=file_name,status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
!         if (debug) then
!           print *,'<<<<< Start reading all Natural Orbital and Occupation number from the NOs_MO.'
!           do j=1,sz_all
!              do i=1,sz_all
!               if(j == i) then
!                  read(read_orb,'(1X,I8,2x,4F28.16,2x,2I4,2x,F28.16)') k,orbital_all(i,j,:),p,q,eigenval_all(i)
!               else 
!                  read(read_orb,'(1X,I8,2x,4F28.16,2x,2I4)') k,orbital_all(i,j,:),p,q
!               end if 
!              end do
!            end do
!            print *,'>>>>> All Orbitals, Read completed'
!         end if
!        close(read_orb,status='KEEP')


        n_orb = 0
        do i=1,sz_all
         if (eigenval_all(i) >= threshold) then
            n_orb = n_orb + 1
         end if
        end do


        sz = n_orb
        allocate(eigenval(sz))
        allocate(orbital(sz_all,sz,4))



        do i=1, sz
         eigenval(i) = eigenval_all(i)
        end do


        do j=1, sz
            do i=1, sz_all
               orbital(i,j,:) = orbital_all(i,j,:)
            end do
        end do

        print *,' Reduce Orbital space.'
        print *,' Now', sz, 'Orbitals'

        deallocate(orbital_all)
        deallocate(eigenval_all)
      end subroutine quater_read_orbital

!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine  quater_read_fock(matrix,n_orb,file_name)

   !     input variables
         real(kind=8), pointer              :: matrix(:,:,:)      ! natural orbital from qdia of DM
         integer,intent(inout)              :: n_orb
         character(len=*),intent(in)        :: file_name

         
   !     loop variables
         integer :: i, j, k, sz, p, q

         integer, parameter :: read_mat = 10
         logical :: debug=.true.

!         sz = n_orb
         allocate(matrix(sz,sz,4))


         open (read_mat,file=file_name,status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
         if (debug) then
           print *,'<<<<< Start reading fock matrix into the file.'
           ! read dimension of the matrix 
           read(read_mat) sz
           allocate(matrix(sz,sz,4))
           ! read matrix
           read(read_mat) matrix

!          do j=1,sz
!              do i=1,sz
!                  read(read_mat,'(2x,4F28.16,2x,2I4)') matrix(i,j,:),p,q
!              end do
!            end do

            print *,'>>>>>  quater fock matrix, reading completed'
         end if
        close(read_mat,status='KEEP')      

      end subroutine quater_read_fock

!-------------------------------------------------------------------------------------------------------------------------------------        

end module write_orbital_to_file

