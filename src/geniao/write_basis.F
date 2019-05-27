!
!   Subroutines (modified from exacorr_dirac.F) needed to 
!   extract basis set info
!   end write it in the end of the .dfcoef file.
!
!
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

       subroutine write_basis_dfcoef
     &  (ao_basis,file_name,iunit)
!      iunit should be the integer corresponding to the .dfcoef file,
!      still opened, in which we have already written the info about
!      CMO coefficients. (CMO, EIG, IBEIG)

       use exacorr_mo
       use exacorr_datatypes
       use exacorr_global
       use exacorr_utils

       implicit none

!      read all positive energy mo coefficients present on DFCOEF file
!      check for imcompatibiklities (use of symmetry) as well and continue
!      by reading also the basis set information into gtos

       character(len=*),intent(in) :: file_name
       type(basis_set_info_t), intent(inout) :: ao_basis
       integer, intent(in)    :: iunit
       integer                :: ierr 
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
       real(8), allocatable   :: coeff(:,:,:),eig(:)
       integer, allocatable   :: ibeig(:),index(:)
       integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
       integer, allocatable   :: jstrt(:)
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
          open (lucoef,file=file_name,status='OLD',FORM='UNFORMATTED',
     &     access='SEQUENTIAL')
          read (lucoef,end=10,err=10) text,nsym,
     &         ((idim(i,j),i=1,3),j=1,nsym),total_energy
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

       allocate(coeff(nao,nto,nz))
       allocate(eig(nto))
       allocate(ibeig(nto))
       allocate(index(2*nmo))
        
       ! read the mo coefficients and orbital eigenvalues
       read (lucoef,end=10,err=10) coeff
       read (lucoef,end=10,err=10) eig
       read (lucoef,end=10,err=10) ibeig

       ! read the basis set information and check for consistency
       read (LUCOEF,err=20,end=20) nao_bas,nshells,npriexp,ncnt 
       if (nao /= nao_bas) then
          ierr = INCONSIS
          print*,"Inconsistency: use NONREL, X2C, or X2CMMF Hamiltonian"
          close (lucoef,status='KEEP')
          return
       end if

       ! allocate the arrays to read the basis set information in DIRAC format
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

       close (lucoef,status='KEEP')

       write (iunit) nao_bas,nshells,npriexp,ncnt 
       write (iunit) nhkt,ncent,nuco
       write (iunit) jstrt,numcf
       write (iunit) cent
       write (iunit) priexp
       write (iunit) priccf

       deallocate(ncent)
       deallocate(nhkt)
       deallocate(nuco)
       deallocate(numcf)
       deallocate(jstrt)
       deallocate(cent)
       deallocate(priexp)
       deallocate(priccf)

       return

  10   ierr = CORRUPT
       close (lucoef,status='KEEP')
       return

  20   ierr = NOBASIS
       close (lucoef,status='KEEP')
       return

       end subroutine write_basis_dfcoef
