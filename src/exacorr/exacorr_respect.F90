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
        integer, parameter     :: INCONSIS = 5 ! RSD_MOS does not contain consistent information wrt number of basis functions
        real(8)                :: total_energy
        logical                :: tobe
        integer                :: lucoef
        integer                :: i, j, k, ishell, iprimstart, iprimend
        integer                :: counter, first
        integer                :: nao, nshells, iao, ibas
        integer                :: nmo, npriexp, ncnt
        integer                :: ncentr, nelms, ncoeff
        complex(8), allocatable:: coeff(:,:)
        integer, allocatable   :: index(:)
        integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
        integer, allocatable   :: jstrt(:)
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

!-------sort and create all necessary arrays from available data

!       sort ReSpect order into exacorr order
        call sort_basis_data(ncent,nhkt,priexp,priccf,nshells)

!       Create array containing centers per shell, instead of per atom
        allocate(cent(3,nshells))
        j = 1
        do i = 1, nshells-1
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

!!       we set the index compatible with the fact that we store all mo's
!        allocate(index(nmo))
!        do i = 1, nmo
!           index(i) = i
!        end do

!       create nuco array -> number of uncontracted functions
        allocate(nuco(nshells))

        counter = 0
        first   = 1
        do i = 1, nshells-1
            if (nhkt(i) == nhkt(i+1)) then
                counter = counter + 1
                if (i == nshells-1) nuco(nshells) = counter
            else
                nuco(first:i) = counter + 1
                first         = i + 1
                counter       = 0
                if (i == nshells-1) nuco(nshells) = 1 
            end if
        end do

!       start of exponent values
        allocate(jstrt(nshells))
        i = 1
        do while (i <= nshells)
            if (i == 1) then
                counter =  nuco(i)
                jstrt(i:counter) = 0
            else if (i == nshells .and. nuco(i) == 1) then
                jstrt(nshells) = counter
            else
                jstrt(i:i+nuco(i)-1) = counter
                counter = counter + nuco(i)
            end if
            i = i + nuco(i)
        end do

!       number of contracted for all shells
        allocate(numcf(nshells))
        counter = 1
        do i = 1, nshells-1
            if (nhkt(i) == nhkt(i+1)) then
                numcf(i) = counter
                counter  = counter + 1
                if (i == nshells-1) numcf(nshells) = counter
            else
                numcf(i) = counter
                counter  = 1
                if (i == nshells-1) numcf(nshells) = counter
            end if
        end do

!-------read total_energy

!        read(lucoef,*,err=10,end=10) total_energy
!        print*, "total energy: ", total_energy

!        return

!-------read Eigenvalues

        read (lucoef,*,end=10,err=10) nelms ! nelms = 2 * nmo
        allocate(eig(nelms))
        read (lucoef,*,end=10,err=10) eig

!       we set the index compatible with the fact that we store all mo's
        allocate(index(nelms))
        do i = 1, nelms
           index(i) = i
        end do

!-------read MO coefficients

        read (lucoef,*,end=10,err=10) ncoeff ! ncoeff = 2*nao * 2*nmo
        allocate(coeff(nao*2,nelms))
        read (lucoef,*,end=10,err=10) coeff

        return

!-------store all data in derived type arguments

!REPLACE THIS WITH CMO
!!       put in qmo
!        call create_qmo (qmo,nao,nmo,nz)
!        qmo%total_energy = total_energy
!        qmo%index        = index
!        qmo%coeff        = coeff(:,npo+1:nto,1:nz)
!        qmo%energy       = eig(1:nmo)
!        qmo%energy       = eig(1:nelms)
!        qmo%boson_irrep  = 0 

        deallocate(index)

!        allocate(priccf(nshells,ncnt))
!
!       put in gtos
        ao_basis%nshells = nshells
        ao_basis%nao     = nao
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
!!!!!!!!!!           ao_basis%gtos(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
           ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
        end do

!-------DONT FORGET TO DEALLOC EVERYTHING----------------------
        deallocate(nhkt  )
        deallocate(ncent )
        deallocate(priexp)
        deallocate(priccf)
        deallocate(cent  )
        deallocate(nuco  )
        deallocate(jstrt )
        deallocate(numcf )
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
