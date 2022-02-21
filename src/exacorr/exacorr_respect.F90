module exacorr_respect

! Module for interfacing the ReSpect SCF with exacorr
! Written by Stan Papadopoulos, January/February 2019

  use exacorr_utils

  implicit none

  private

  public read_from_respect

  contains

    subroutine read_from_respect (ao_basis,cspinor,ierr)
      use labeled_storage
!     Routine to read coefficients and ao_basis data from ReSpect

      use exacorr_mo
      use exacorr_datatypes, only : basis_set_info_t

!     input variables
      type(basis_set_info_t), intent(inout) :: ao_basis
      type(cmo),intent(inout) :: cspinor
      integer,  intent(  out) :: ierr
!     error messages
      integer, parameter      :: NOERROR  = 0 ! no error
      integer, parameter      :: NOTFOUND = 1 ! file not found
      integer, parameter      :: CORRUPT  = 2 ! file unreadable
      integer, parameter      :: NOBASIS  = 3 ! no basis set information could be read (old RSD_MOS?)
      integer, parameter      :: INCONSIS = 4 ! RSD_MOS does not contain consistent information wrt number of basis functions
!     file reading
      logical                 :: tobe
      integer                 :: lucoef
!     data to be read/stored
      complex(8), allocatable :: coeff(:,:)
      complex(8), allocatable :: sphcoeff(:), overlap(:), &
                                 diplenx(:), dipleny(:), diplenz(:)
      integer                 :: nao, nelms, nshells, iao, ibas
      integer                 :: npriexp, ncnt
      integer                 :: ncentr, ncoeff,nvec
      integer, allocatable    :: index(:)
      integer, allocatable    :: ncent(:), nhkt(:)
      integer, allocatable    :: natomc(:)
      real(8)                 :: total_energy
      real(8), allocatable    :: cent(:,:), cent_rsd(:,:)
      real(8), allocatable    :: priexp(:), priccf(:,:)
      real(8), allocatable    :: eig(:)
!     summation/count variables
      integer                 :: i, j, ishell, nfunctions
!     ReSpect HDF5 data variables
      logical                 :: doreex=.true.
      character(6)            :: reex(6)
      type(file_info_t), save :: reex_file
  
      ierr = NOERROR

!****************
!     open file *
!****************

      inquire (file='RSD_MOS',exist=tobe)
      if (.not.tobe) then
        print*, "error: MO coefficient file RSD_MOS not found"
        ierr = NOTFOUND
        return
      else
        call get_free_fileunit(lucoef)
        open (lucoef,file='RSD_MOS',status='unknown',form='formatted',position='rewind',access='sequential')
      end if

!***********************
!     read atomic data *
!***********************

!     number of atoms
      read (lucoef,*,end=10,err=10) ncentr
      allocate(natomc(ncentr))
      allocate(cent_rsd(3,ncentr))

!     number of protons per atom
      read (lucoef,*,end=10,err=10) natomc
!     coordinates of atoms
      read (lucoef,*,end=10,err=10) cent_rsd

      deallocate(natomc) ! not needed

!**************************
!     read basis set data *
!**************************

      ncnt = 1 ! read_from_respect only supported for uncontracted basis sets

!     number of shells
      read (lucoef,*,end=20,err=20) nshells

      npriexp = nshells ! uncontracted basis set

      allocate(nhkt(       nshells))
      allocate(ncent(      nshells))
      allocate(priexp(     npriexp))
      allocate(priccf(ncnt,npriexp))

!     orbital momentum
      read (lucoef,*,end=20,err=20) nhkt
!     atom number
      read (lucoef,*,end=20,err=20) ncent
!     primitive exponents
      read (lucoef,*,end=20,err=20) priexp
!     contraction coefficients
      read (lucoef,*,end=20,err=20) priccf

!     calculate number of functions
      nao = 0
      do i = 1, nshells
!       number of Cartesian GTOs / "electronic" (=positive energy) solutions"
        nao = nao + nhkt(i) * ( nhkt(i) + 1 ) / 2
      end do

!*************************************************************
!     read total SCF energy, eigenvalues and MO coefficients *
!*************************************************************

!     SCF energy
      read(lucoef,*,err=10,end=10) total_energy

!     read eigenvalues
      read (lucoef,*,end=10,err=10) nelms
      allocate(eig(nelms))
      read (lucoef,*,end=10,err=10) eig

!     we set the index compatible with the fact that we store all mo's
      allocate(index(nelms))
      do i = 1, nelms
        index(i) = i
      end do

!     read MO coefficients
      read (lucoef,*,end=10,err=10) reex(1), ncoeff ! ncoeff = 2*nao * nelms = 2*nao * 2*nmo

      if (2*nao*nelms /= ncoeff) goto 10 ! sanity check

      allocate(coeff(2*nao,nelms))
      read (lucoef,*,end=10,err=10) coeff

!************************************************************************
!     read MO coefficients and property matrices in the spherical basis *
!************************************************************************
      if (doreex) then

!     read MO coefficients in the spherical basis
      read (lucoef,*,end=10,err=10) reex(2), nvec !nvec = no. array elements
      allocate(sphcoeff(nvec))
      read (lucoef,*,end=10,err=10) sphcoeff

!     read overlap matrix in the spherical basis
      read (lucoef,*,end=10,err=10) reex(3), nvec !nvec = no. array elements
      allocate(overlap(nvec))
      read (lucoef,*,end=10,err=10) overlap

!     read dipole length x matrix in the spherical basis
      read (lucoef,*,end=10,err=10) reex(4), nvec !nvec = no. array elements
      allocate(diplenx(nvec))
      read (lucoef,*,end=10,err=10) diplenx

!     read dipole length y matrix in the spherical basis
      read (lucoef,*,end=10,err=10) reex(5), nvec !nvec = no. array elements
      allocate(dipleny(nvec))
      read (lucoef,*,end=10,err=10) dipleny

!     read dipole length z matrix in the spherical basis
      read (lucoef,*,end=10,err=10) reex(6), nvec !nvec = no. array elements
      allocate(diplenz(nvec))
      read (lucoef,*,end=10,err=10) diplenz

      endif
  
!*********************************************
!       sort and create all necessary arrays *
!*********************************************

!     sort ReSpect order into exacorr order
      if (ncentr.gt.1) call sort_data(ncent,nhkt,priexp,priccf,coeff,nshells)

!     Create array containing centers per shell, instead of per atom
      allocate(cent(nshells,3))
      j = 1
      do i = 1, nshells-1
        cent(i,1:3) = cent_rsd(1:3,j)
        if (ncent(i) < ncent(i+1)) then
          j = j + 1
        else if (i == nshells-1) then
          cent(nshells,1:3) = cent_rsd(1:3,j)
        end if
      end do
      deallocate(cent_rsd)

!***********************************************
!     store all data in derived type arguments *
!***********************************************
!     CHIMA STUFF
!     Write ExaCorr sorted (w.r.t. AO basis) coefficients to (HDF5) file!
      write(*,*) "Writing ReSpect Coefficients to file!"
      !write(*,*) "nao:", 2*nao
      !write(*,*) "nmo:", nelms

      if (reex_file%status == -1) then
         write(*,*) " Opening HDF5 REEX file"
         reex_file%type = 2
         reex_file%name = "REEX.h5"
      end if

      write(*,*) " Writing REEX data to REEX file"
      call lab_write (reex_file,'nmo',idata=nelms)
      call lab_write (reex_file,'nsph',idata=int(sqrt(real(nvec))))
      call lab_write (reex_file,reex(1)//'_RE',rdata=reshape(transpose(real(coeff,8)), shape=(/2*nao*nelms/) ) )
      call lab_write (reex_file,reex(1)//'_IM',rdata=reshape(transpose(dimag(coeff )), shape=(/2*nao*nelms/) ) )

      if (doreex) then
!*****************************************************************************
!     store MO coefficients and property matrices in spherical basis on file *
!*****************************************************************************

!     coefficients
      call lab_write (reex_file,reex(2)//'_RE',rdata=reshape(real( sphcoeff,8), shape=(/nvec/) ) )
      call lab_write (reex_file,reex(2)//'_IM',rdata=reshape(dimag(sphcoeff  ), shape=(/nvec/) ) )
!     overlap
      call lab_write (reex_file,reex(3)//'_RE',rdata=reshape(real( overlap, 8), shape=(/nvec/) ) )
      call lab_write (reex_file,reex(3)//'_IM',rdata=reshape(dimag(overlap   ), shape=(/nvec/) ) )
!     dipole x
      call lab_write (reex_file,reex(4)//'_RE',rdata=reshape(real( diplenx, 8), shape=(/nvec/) ) )
      call lab_write (reex_file,reex(4)//'_IM',rdata=reshape(dimag(diplenx   ), shape=(/nvec/) ) )
!     dipole y
      call lab_write (reex_file,reex(5)//'_RE',rdata=reshape(real( dipleny, 8), shape=(/nvec/) ) )
      call lab_write (reex_file,reex(5)//'_IM',rdata=reshape(dimag(dipleny   ), shape=(/nvec/) ) )
!     dipole z      
      call lab_write (reex_file,reex(6)//'_RE',rdata=reshape(real( diplenz, 8), shape=(/nvec/) ) )
      call lab_write (reex_file,reex(6)//'_IM',rdata=reshape(dimag(diplenz   ), shape=(/nvec/) ) )
      
      endif
  
!     put in cmo
      call alloc_mo(cspinor,nao,nelms)
      cspinor%total_energy = total_energy
      cspinor%index        = index
      cspinor%coeff(:,:,1) = coeff(1:nao,:)
      cspinor%coeff(:,:,2) = coeff(1+nao:2*nao,:)
      cspinor%energy       = eig(:)
      cspinor%symmetry     = 0 

      deallocate(coeff)
      deallocate(eig  )
      deallocate(index)

!     put in gtos
      ao_basis%nshells       = nshells
      ao_basis%nao           = nao
      ao_basis%basis_angular = 1
      nullify(ao_basis%gtos)
      allocate(ao_basis%gtos(nshells))
      do ishell = 1, nshells
        ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
        ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
        ao_basis%gtos(ishell)%n_primitives = 1
        nullify  (ao_basis%gtos(ishell)%exponent)
        allocate (ao_basis%gtos(ishell)%exponent(1))
        nullify  (ao_basis%gtos(ishell)%coefficient)
        allocate (ao_basis%gtos(ishell)%coefficient(1))
        ao_basis%gtos(ishell)%exponent = priexp(ishell)
        ao_basis%gtos(ishell)%coefficient = priccf(1,ishell)
        ao_basis%gtos(ishell)%coord = cent(ishell,1:3)         
      end do

      deallocate(nhkt  )
      deallocate(ncent )
      deallocate(priexp)
      deallocate(priccf)
      deallocate(cent  )

!     make array with pointers from a basis function index to its shell
      nullify(ao_basis%shell_indices)
      allocate(ao_basis%shell_indices(nao))
      iao = 0
      do ishell = 1, nshells
        nfunctions = ao_basis%gtos(ishell)%orb_momentum * ( ao_basis%gtos(ishell)%orb_momentum + 1 ) / 2
        do ibas = 1, nfunctions
          iao = iao + 1
          if (iao.gt.nao) then
            ierr = INCONSIS
            call quit ('Error in read_from_respect')
          endif
          ao_basis%shell_indices(iao) = ishell
        end do
      end do

      if (iao.gt.nao) then
        ierr = INCONSIS
        call quit ('Error in read_from_respect')
      endif

      close (lucoef,status='KEEP')

      call print_date("Finished reading data from RSD_MOS")
      return

10    ierr = CORRUPT
      return

20    ierr = NOBASIS
      return

      end subroutine read_from_respect

      subroutine sort_data(ncent,nhkt,priexp,priccf,coeff,nshells)

!       Sort of MO coefficients and basis set data

!       input variables
        integer,    intent(inout) :: ncent(:), nhkt(:)
        integer,    intent(in   ) :: nshells
        real(8),    intent(inout) :: priexp(:), priccf(:,:)
        complex(8), intent(inout) :: coeff(:,:)
!       local variables
        integer :: i, j, index, index2
        integer :: tmp_ncent, tmp_nhkt
        integer :: dims(2)
        integer :: nao, nfunctions, ncentr
        real(8) :: tmp_priexp, tmp_priccf
        complex(8), allocatable :: tmp_coeff(:,:)

        !  extract variables from arguments
        dims   = shape(coeff)
        nao    = dims(1) / 2
        ncentr = maxval(ncent)

        ! make tmp copy of coefficients
        allocate(tmp_coeff(dims(1),dims(2)))
        tmp_coeff = coeff

        ! MO coefficient matrix ReSpect
        !       2*nmo
        !     ---------
        ! nao | alpha |
        !     ---------
        ! nao | beta  |
        !     ---------

        ! sort MO coefficients
        i = 1
        index = 1
        do while (i.le.ncentr)
          index2 = 1
          do j = 1, nshells
            nfunctions = nhkt(j) * ( nhkt(j) + 1 ) / 2
            if (i == ncent(j)) then
              ! alpha
              coeff(index:index+nfunctions-1,:) = tmp_coeff(index2:index2+nfunctions-1,:)
              ! beta
              coeff(index+nao:index+nao+nfunctions-1,:) = tmp_coeff(index2+nao:index2+nao+nfunctions-1,:)
              index = index + nfunctions
            end if
            index2 = index2 + nfunctions
          end do
          i = i + 1
        end do

        deallocate(tmp_coeff)

        do i = 2, nshells
          tmp_ncent  = ncent (  i)
          tmp_nhkt   = nhkt  (  i)
          tmp_priexp = priexp(  i)
          tmp_priccf = priccf(1,i)
          j = i - 1
          do while (j >= 1)
            if (ncent(j) <= tmp_ncent) exit
            ncent (  j+1) = ncent (  j)
            nhkt  (  j+1) = nhkt  (  j)
            priexp(  j+1) = priexp(  j)
            priccf(1,j+1) = priccf(1,j)
            j = j - 1
          end do
          ncent (  j+1) = tmp_ncent
          nhkt  (  j+1) = tmp_nhkt
          priexp(  j+1) = tmp_priexp
          priccf(1,j+1) = tmp_priccf
        end do

      end subroutine sort_data

end module exacorr_respect
