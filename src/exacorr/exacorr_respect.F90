module exacorr_respect

!   Module for interfacing the ReSpect SCF with exacorr

    implicit none

    private

    public read_from_respect

    contains

        subroutine read_from_respect (ao_basis,como,ierr)

        use exacorr_mo
        use exacorr_datatypes, only : basis_set_info_t
        use exacorr_global, only    : ncartesian

!       read all positive energy mo coefficients present on DFCOEF file
!       check for imcompatibiklities (use of symmetry) as well and continue
!       by reading also the basis set information into gtos
!       ReSpect does not use symmetry

        type(basis_set_info_t), intent(inout) :: ao_basis
        type(cmo),intent(inout):: como
        integer, intent(out)   :: ierr 
        integer, parameter     :: NOERROR  = 0 ! no error
        integer, parameter     :: NOTFOUND = 1 ! file not found
        integer, parameter     :: CORRUPT  = 2 ! file unreadable
        integer, parameter     :: NOBASIS  = 3 ! no basis set information could be read (old DFCOEF?)
        integer, parameter     :: INCONSIS = 4 ! RSD_MOS does not contain consistent information wrt number of basis functions
        real(8)                :: total_energy
        logical                :: tobe
        integer                :: lucoef
        integer                :: i, j, ishell
        integer                :: nao, nmo, nshells, iao, ibas
        integer                :: npriexp, ncnt
        integer                :: ncentr, ncoeff
        complex(8), allocatable:: coeff(:,:)
        integer, allocatable   :: index(:)
        integer, allocatable   :: ncent(:),nhkt(:)
        integer, allocatable   :: natomc(:)
        real(8), allocatable   :: cent(:,:),cent_rsd(:,:)
        real(8), allocatable   :: priexp(:),priccf(:,:)
        real(8), allocatable   :: eig(:)

        ierr = NOERROR

!-------open file

        inquire (file='RSD_MOS',exist=tobe)
        if (.not.tobe) then
           print*, "error: MO coefficient file RSD_MOS not found"
           ierr = NOTFOUND
           return
        else
           call get_free_fileunit(lucoef)
           open (lucoef,file='RSD_MOS',status='unknown',form='formatted',position='rewind',access='sequential')
        end if

!-------read atomic data

!       number of atoms
        read (lucoef,*,end=10,err=10) ncentr
        allocate(natomc(ncentr))
        allocate(cent_rsd(3,ncentr))

!       number of protons per atom
        read (lucoef,*,end=10,err=10) natomc
!       coordinates of atoms
        read (lucoef,*,end=10,err=10) cent_rsd

        deallocate(natomc) ! not needed

!-------read basis set data

        ncnt = 1 ! read_from_respect only supported for uncontracted basis sets

!       number of shells
        read (lucoef,*,end=20,err=20) nshells

        npriexp = nshells ! uncontracted basis set

        allocate(nhkt(       nshells))
        allocate(ncent(      nshells))
        allocate(priexp(     npriexp))
        allocate(priccf(ncnt,npriexp))

!       orbital momentum
        read (lucoef,*,end=20,err=20) nhkt
!       atom number
        read (lucoef,*,end=20,err=20) ncent
!       primitive exponents
        read (lucoef,*,end=20,err=20) priexp
!       contraction coefficients
        read (lucoef,*,end=20,err=20) priccf

!       calculate number of functions
        nao = 0
        do i = 1, nshells
!           number of sperical GTOs / "electronic" (=positive energy) solutions"
            nao = nao + 2 * (nhkt(i) - 1) + 1
        end do
        nao = nao * 2

!-------sort and create all necessary arrays from available data

!       sort ReSpect order into exacorr order
        call sort_basis_data(ncent,nhkt,priexp,priccf,nshells)

!       Create array containing centers per shell, instead of per atom
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

!-------read total SCF energy

        read(lucoef,*,err=10,end=10) total_energy

!-------read eigenvalues

        read (lucoef,*,end=10,err=10) nmo
        allocate(eig(nmo))
        read (lucoef,*,end=10,err=10) eig

!       we set the index compatible with the fact that we store all mo's
        allocate(index(nmo))
        do i = 1, nmo
           index(i) = i
        end do

!-------read MO coefficients

        read (lucoef,*,end=10,err=10) ncoeff ! ncoeff = nao * nmo

        if (nao*nmo /= ncoeff) goto 10 ! sanity check

        allocate(coeff(nao,nmo))
        read (lucoef,*,end=10,err=10) coeff

!-------store all data in derived type arguments

!       put in cmo
        call alloc_mo(como,nao,nmo)
        como%total_energy = total_energy
        como%index        = index
        como%coeff        = coeff(:,:)
        como%energy       = eig(:)
        como%boson_irrep  = 0 

        deallocate(index)

!       put in gtos
        ao_basis%nshells = nshells
        ao_basis%nao     = nao
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
        deallocate(eig   )
        deallocate(coeff )

!       make array with pointers from a basis function index to its shell
        nullify(ao_basis%shell_indices)
        allocate(ao_basis%shell_indices(nao))
        iao = 0
        do ishell = 1, nshells
           do ibas = 1, ncartesian(ao_basis%gtos(ishell)%orb_momentum)
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

        subroutine sort_basis_data(ncent,nhkt,priexp,priccf,nshells)

!        Insertion sort (so far the only sort which gives the desired ordening)

!        input variables
         integer, intent(inout) :: ncent(:), nhkt(:)
         integer, intent(in   ) :: nshells
         real(8), intent(inout) :: priexp(:), priccf(:,:)
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
