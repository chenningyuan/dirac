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

        type(basis_set_info_t), intent(inout) :: ao_basis
        type(MO), intent(inout):: qmo
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
        complex(8), allocatable:: coeff_a(:,:),coeff_b(:,:)
        real(8), allocatable   :: coeff(:,:,:),eig(:)
        integer, allocatable   :: ibeig(:),index(:)
        integer, allocatable   :: ncent(:),nhkt(:),nuco(:),numcf(:)
        integer, allocatable   :: jstrt(:)
        real(8), allocatable   :: cent(:,:),priexp(:),priccf(:,:)

        ierr = NOERROR

!       open file and read dimension info
        inquire (file='DFCOEF',exist=tobe)
        if (.not.tobe) then
           print*, "error: MO coefficient file DFCOEF not found"
           ierr = NOTFOUND
           return
        else
           call get_free_fileunit(lucoef)
           open (lucoef,file='DFCOEF',status='OLD',FORM='UNFORMATTED',access='SEQUENTIAL')
           read (lucoef,end=10,err=10) text,nsym,((idim(i,j),i=1,3),j=1,nsym),total_energy
        end if

!       check symmetry and get dimensions of the coefficient array
        if (nsym > 1) then
           ierr = SYMMETRY
           return
        else
           npo = idim(1,1) ! "positronic" (= negative energy) solutions (can be thrown away for our purposes)
           nmo = idim(2,1) ! "electronic" (= positive energy) solutions (we only need these)
           nao = idim(3,1) ! number of basis functions
           nto = npo+nmo   ! total number of solutions on this file
           if (npo > 0) print*, "warning: reading from 4C DFCOEF file" ! this warning will go once we have the full index transform implementedd
        end if

        !allocate(coeff(nao,nto,nz))
        allocate(coeff_a(nao,nto))
        allocate(coeff_b(nao,nto))
        allocate(eig(nto))
        allocate(ibeig(nto))
        allocate(index(nto))
         
!       read the mo coefficients and orbital eigenvalues
        read (lucoef,end=10,err=10) coeff_a
        read (lucoef,end=10,err=10) coeff_b
        read (lucoef,end=10,err=10) eig
        read (lucoef,end=10,err=10) ibeig
!       fill in index values
        do i = 1, nmo
           index(i) = i ! we set the index compatible with the fact that we store all mo's
        end do

        allocate(coeff(nao,nto,nz))

!       Transformation of complex coefficients to quaternion units
        call coeff_comp_to_quaternion(coeff,coeff_a,nao,nmo,1)
        call coeff_comp_to_quaternion(coeff,coeff_b,nao,nmo,2)

        deallocate(coeff_a)
        deallocate(coeff_b)

!       put in qmo
        call create_qmo (qmo,nao,nmo,nz)
        qmo%total_energy = total_energy
        qmo%index        = index
        qmo%coeff        = coeff(:,npo+1:nto,1:nz)
        qmo%energy       = eig(npo+1:nto)
        qmo%boson_irrep  = ibeig(npo+1:nto)

        deallocate(coeff)
        deallocate(eig)
        deallocate(ibeig)
        deallocate(index)

!       read the basis set information and check for consistency
        read (LUCOEF,err=20,end=20) nshells,npriexp,ncnt 

!       allocate the arrays to read the basis set information in DIRAC format
        allocate(ncent(nshells))
        allocate(nhkt(nshells))
        allocate(nuco(nshells))
        allocate(numcf(nshells))
        allocate(jstrt(nshells))
        allocate(cent(nshells,3))
        allocate(priexp(npriexp))
        allocate(priccf(npriexp,ncnt))

!       orb_momentum, atom_number, n_primives for all shells
        read (LUCOEF,err=20,end=20) nhkt,ncent,nuco
!       start of exponent values, number of contracted for all shells
        read (LUCOEF,err=20,end=20) jstrt,numcf
!       coordinates (x, y, z) of center for all shells
        read (LUCOEF,err=20,end=20) cent
!       values of the exponents for all shells (packed into a 1-d array)
        read (LUCOEF,err=20,end=20) priexp
!       values of the contraction coefficients for all shells (packed into a 1-d array)
        read (LUCOEF,err=20,end=20) priccf

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
           ao_basis%gtos(ishell)%coefficient = priccf(iPrimStart:iPrimEnd,k)
           ao_basis%gtos(ishell)%coord = cent(ishell,1:3)
!          In DIRAC uncontracted shells are stored as quasi-contracted, fix this.
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

        call print_date("read from ReSpect:"//text)
        return

  10    ierr = CORRUPT
        return

  20    ierr = NOBASIS
        return

        end subroutine read_from_respect

        subroutine coeff_comp_to_quaternion(coeff,coeff_comp,nao,nmo,spin)

         real(8), intent(out)   :: coeff(:,:,:)
         complex(8), intent(in) :: coeff_comp(:,:) !???
         integer, intent(in)    :: nao, nmo, spin
         integer                :: j, k

!        Copy the coefficients into the quaternion array, using the quaternion relations as defined in
!        T. Saue, H. Jensen, J. Chem. Phys. 111 (1999) 6211–6222, Equation 16.

         if (spin == 1) then
            do k = 1, nmo
               do j = 1, nao
                !quaternion to complex
                !this%mocoef_all(j,k)     = dcmplx(qmo%coeff(j,k,1),qmo%coeff(j,k,2))

                !complex to quaternion
                coeff(j,k,1) = dreal(coeff_comp(j,k))
                coeff(j,k,2) = dimag(coeff_comp(j,k))
               end do
            end do
         elseif (spin == 2) then
            do k = 1, nmo
               do j = 1, nao
                !quaternion to complex
                !this%mocoef_all(j,k)     = dcmplx(-qmo%coeff(j,k,3),qmo%coeff(j,k,4))

                !complex to quaternion
                coeff(j,k,3) = -dreal(coeff_comp(j,k))
                coeff(j,k,4) =  dimag(coeff_comp(j,k))
               end do
            end do
         else
            call quit ('coeff_comp_to_quaternion called with unknown spin value')
         end if
 
        end subroutine coeff_comp_to_quaternion

end module exacorr_respect
