module exacorr_respect

    implicit none

    private

    public read_from_respect

    contains

        subroutine read_from_respect (ao_basis,qmo,ierr)

        use dirac_mo
        use exacorr_datatypes, only : basis_set_info_t
        use exacorr_global, only    : ncartesian

!       read all positive energy mo coefficients present on DFCOEF file
!       check for imcompatibiklities (use of symmetry) as well and continue
!       by reading also the basis set information into gtos
!       ReSpect does not use symmetry

        type(basis_set_info_t), intent(inout) :: ao_basis
        type(MO), intent(inout):: qmo
        integer, intent(out)   :: ierr 
        integer, parameter     :: NOERROR  = 0 ! no error
        integer, parameter     :: NOTFOUND = 1 ! file not found
        integer, parameter     :: CORRUPT  = 2 ! file unreadable
        integer, parameter     :: NOBASIS  = 4 ! no basis set information could be read (old DFCOEF?)
        integer, parameter     :: INCONSIS = 6 ! DFCOEF does not contain consistent information wrt number of basis functions
        !character              :: text*74      ! holds text read from DFCOEF, only used for printing this
        integer                :: idim(3,2)
        real(8)                :: total_energy
        logical                :: tobe
        integer                :: lucoef
        integer                :: i, j, k, ishell, iprimstart, iprimend
        integer                :: nao, nshells, iao, ibas
        integer                :: nmo, npriexp, ncnt
        integer, parameter     :: nz = 4 ! we use no symmetry, so always have this
        real(8), allocatable   :: coeff(:,:,:),eig(:)
        integer, allocatable   :: index(:)
        integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
        integer, allocatable   :: jstrt(:)
        real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)
!--------------------------------------------------------------------------------------
        integer                 :: ncentr
        integer                 :: nelms
        integer                 :: ncoeff
        integer,    allocatable :: natomc(:)
        real(8),    allocatable :: cent_rsd(:,:)
        real(8),    allocatable :: eig_rsd(:)
        complex(8), allocatable :: coeff_rsd(:)

        ierr = NOERROR

!       open file and read dimension info
        inquire (file='RSD_MOS',exist=tobe)
        if (.not.tobe) then
           print*, "error: MO coefficient file RSD_MOS not found"
           ierr = NOTFOUND
           return
        else
           call get_free_fileunit(lucoef)
           open (lucoef,file='RSD_MOS',status='unknown',form='formatted',position='rewind',access='sequential')
        end if

!-------Read atomic data

!       number of atoms
        read (lucoef,*,end=10,err=10) ncentr

        allocate(natomc(ncentr))
        allocate(cent_rsd(3,ncentr))

!       number of protons per atom
        read (lucoef,*,end=10,err=10) natomc
!       coordinates of atoms
        read (lucoef,*,end=10,err=10) cent_rsd

        deallocate(natomc)

        ncnt = 1 ! read_from_respect only supported for uncontracted basis sets

!-------Read basis set data

!       number of shells
        read (lucoef,*,end=20,err=20) nshells

        allocate(nhkt(nshells))
        allocate(ncent(nshells))
        allocate(priexp(nshells))
        allocate(priccf(ncnt,nshells))

!       orbital momentum
        read (lucoef,*,end=20,err=20) nhkt
!       atom number
        read (lucoef,*,end=20,err=20) ncent
!       primitive exponents
        read (lucoef,*,end=20,err=20) priexp
!       contraction coefficients
        read (lucoef,*,end=20,err=20) priccf

!       sort ReSpect order into exacorr order
        call sort_basis_data(ncent,nhkt,priexp,priccf,nshells)

!       Create array containing centers per shell, instead of per atom
        allocate(cent(3,nshells))
        j = 1
        do i = 1, nshells - 1
            cent(1:3,i) = cent_rsd(1:3,j)
            if (ncent(i) < ncent(i+1)) then
                j = j + 1
            else if (i == nshells-1) then
                cent(1:3,nshells) = cent_rsd(1:3,j)
            end if
        end do
        deallocate(cent_rsd)

!       calculate number of functions
        nao = 0
        nmo = 0
        do i = 1, nshells
!           number of Cartesian GTOs / "electronic" (=positive energy) solutions"
            nao = nao + nhkt(i) * (nhkt(i) + 1) / 2
!           number of sperical GTOs
            nmo = nmo + 2 * (nhkt(i) - 1) + 1
        end do

!       we set the index compatible with the fact that we store all mo's
        allocate(index(nmo))
        do i = 1, nmo
           index(i) = i
        end do

        print*, "nao: ", nao, "nmo: ", nmo
        print*, "          nhkt      ncent   priexp              priccf"
        do i = 1, nshells
            print*, nhkt(i), ncent(i), priexp(i), priccf(1,i)
        end do

        print*, "cent"
        do i = 1, nshells
            print*, cent(1,i), cent(2,i), cent(3,i)
        end do

        return

!-------read Eigenvalues

        read (lucoef,*,end=10,err=10) nelms ! nelms = 2 * nmo
        allocate(eig_rsd(nelms))
        read (lucoef,*,end=10,err=10) eig_rsd

        return

!-------read MO coefficients

        read (lucoef,*,end=10,err=10) ncoeff ! ncoeff = 2*nao * 2*nmo
        allocate(coeff_rsd(nao*2,nelms))
        read (lucoef,*,end=10,err=10) coeff_rsd

!-------DONT FORGET TO DEALLOC EVERYTHING----------------------

!------------------------------------------------------------------------

!        !allocate(coeff(nao,nto,nz))
!        allocate(coeff_a(nao,nto))
!        allocate(coeff_b(nao,nto))
!        allocate(eig(nto))
!         
!!       read the mo coefficients and orbital eigenvalues
!        read (lucoef,end=10,err=10) coeff_a
!        read (lucoef,end=10,err=10) coeff_b
!        read (lucoef,end=10,err=10) eig
!!       fill in index values
!
!        allocate(coeff(nao,nto,nz))
!
!!       Transformation of complex coefficients to quaternion units
!        call coeff_comp_to_quaternion(coeff,coeff_a,nao,nmo,1)
!        call coeff_comp_to_quaternion(coeff,coeff_b,nao,nmo,2)
!
!        deallocate(coeff_a)
!        deallocate(coeff_b)
!
!!       put in qmo
!        call create_qmo (qmo,nao,nmo,nz)
!        qmo%total_energy = total_energy
!        qmo%index        = index
!        qmo%coeff        = coeff(:,npo+1:nto,1:nz)
!        qmo%energy       = eig(npo+1:nto)
!        qmo%boson_irrep  = 0 
!
!        deallocate(coeff)
!        deallocate(eig)
!        deallocate(ibeig)
!        deallocate(index)
!
!!       read the basis set information and check for consistency
!        read (LUCOEF,err=20,end=20) nshells,npriexp,ncnt 
!
!!       allocate the arrays to read the basis set information in DIRAC format
!        allocate(ncent(nshells))
!        allocate(nhkt(nshells))
!        allocate(nuco(nshells))
!        allocate(numcf(nshells))
!        allocate(jstrt(nshells))
!        allocate(cent(nshells,3))
!        allocate(priexp(npriexp))
!        allocate(priccf(npriexp,ncnt))
!
!!       orb_momentum, atom_number, n_primives for all shells
!        read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
!!       start of exponent values, number of contracted for all shells
!        read (LUCOEF,err=20,end=20) jstrt,numcf
!!       coordinates (x, y, z) of center for all shells
!        read (LUCOEF,err=20,end=20) cent
!!       values of the exponents for all shells (packed into a 1-d array)
!        read (LUCOEF,err=20,end=20) priexp
!!       values of the contraction coefficients for all shells (packed into a 1-d array)
!        read (LUCOEF,err=20,end=20) priccf
!
!!       put in gtos
!        ao_basis%nshells = nshells
!        ao_basis%nao     = nao
!        nullify(ao_basis%gtos)
!        allocate(ao_basis%gtos(nshells))
!        do ishell = 1, nshells
!           ao_basis%gtos(ishell)%orb_momentum = nhkt(ishell)
!           ao_basis%gtos(ishell)%atom_number  = ncent(ishell)
!           ao_basis%gtos(ishell)%n_primitives = nuco(ishell)
!           nullify  (ao_basis%gtos(ishell)%exponent) 
!           allocate (ao_basis%gtos(ishell)%exponent(nuco(ishell)))
!           nullify  (ao_basis%gtos(ishell)%coefficient)
!           allocate (ao_basis%gtos(ishell)%coefficient(nuco(ishell)))
!           iPrimStart = jstrt(ishell) + 1
!           iPrimEnd   = jstrt(ishell) + nuco(ishell)
!           k          = numcf(ishell)
!           ao_basis%gtos(ishell)%exponent = priexp(iPrimStart:iPrimEnd)
!           ao_basis%gtos(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
!           ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
!!          In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
!           call compress_unc_shell (ao_basis%gtos(ishell),ierr)
!           if (ierr == 1) ierr = CONTRAC
!        end do
!
!        deallocate(ncent)
!        deallocate(nhkt)
!        deallocate(nuco)
!        deallocate(numcf)
!        deallocate(jstrt)
!        deallocate(cent)
!        deallocate(priexp)
!        deallocate(priccf)
!
!!       make array with pointers from a basis function index to its shell
!        nullify(ao_basis%shell_indices)
!        allocate(ao_basis%shell_indices(nao))
!        iao = 0
!        do ishell = 1, nshells
!           do ibas = 1, ncartesian(ao_basis%gtos(ishell)%orb_momentum)
!              iao = iao + 1
!              if (iao.gt.nao) then
!                 ierr = INCONSIS
!                 call quit ('Error in read_from_respect')
!              endif
!              ao_basis%shell_indices(iao) = ishell
!           end do
!        end do
!        if (iao.gt.nao) then
!           ierr = INCONSIS
!           call quit ('Error in read_from_respect')
!        endif

        close (lucoef,status='KEEP')

        call print_date("Finished reading data from RSD_MOS")
        return

  10    ierr = CORRUPT
        return

  20    ierr = NOBASIS
        return

        end subroutine read_from_respect

        subroutine sort_basis_data(ncent,nhkt,priexp,priccf,nshells)

!        Insertion sort (so far the only sort which gives the desired ordening)

         integer, intent(inout) :: ncent(:), nhkt(:)
         real(8), intent(inout) :: priexp(:), priccf(:,:)
         integer, intent(in)    :: nshells
!        local variables
         integer :: tmp_ncent, tmp_nhkt
         real(8) :: tmp_priexp, tmp_priccf
         integer :: i, j

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

        end subroutine sort_basis_data

end module exacorr_respect
