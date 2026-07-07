module rose_import_fragments

implicit none
private

public rose_define_fragments, rose_read_fragments, rose_read_mos, init_fragment

contains

   subroutine rose_define_fragments (molecule,fragments)

        use rose_input_data
        use read_xyzfile
        use read_xyzfile, only: atom

        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragment
        type(atom), pointer                :: current_atom ! linked list to store coordinate information
        integer                            :: i

        ! n_fragments == 0 signals that we use atomic fragments, otherwise molecular fragments indicated as frag1, frag2, etc.
        if (n_fragments == 0) then
            if (.not. file_for_MOs .eq. "adf" .and. .not. file_for_MOs .eq. "ams") atom_fragments = .true.
            n_fragments = molecule%natoms
            allocate(fragments(n_fragments))
            allocate(frag_ncor(n_fragments))
            allocate(frag_nval(n_fragments))
            allocate(frag_bias(n_fragments))
            frag_ncor = -1 ! signal that this is not set by the user
            frag_nval = -1 ! signal that this is not set by the user
            frag_bias = 0.0 ! by default we do not bias (except the first fragment, but this is done later)
        else
            allocate(fragments(n_fragments))
        end if

        if (atom_fragments) then
            ! We make the fragments using information from the supermolecule
           current_atom => molecule%firstatom
           do i = 1, n_fragments
              fragments(i)%natoms  = 1
              fragments(i)%name    = current_atom%name
              fragments(i)%Znumber = current_atom%bascharge
              allocate(fragments(i)%firstatom)
              fragments(i)%firstatom = current_atom
              nullify(fragments(i)%firstatom%next)
              current_atom => current_atom%next
           end do
        else
           do i = 1, n_fragments
              write(fragments(i)%name,'(a4,i0)') 'frag',i ! Used to find the fragment checkpoint file
              call rose_read_structure (fragments(i))
              call rose_print_structure (fragments(i))
           enddo
        end if

        ! Set the predefined information about core and valence
        call set_coreval_info (molecule,include_core)
        do i = 1, n_fragments
           call set_coreval_info (fragments(i),include_core)
        end do

        ! Adjust if necessary by filling in user-supplied information
        do i = 1, n_fragments
           ! Overwrite with user-defined information (if not set, this has a value of -1)
           ! We have minimal sanity checks, user should make sure ncor and nval make sense for these fragments
           if (allocated(frag_nval)) then
              if (frag_nval(i) > -1) fragments(i)%nvale = frag_nval(i)
           end if
           if (allocated(frag_ncor)) then
              if (frag_ncor(i) > -1) fragments(i)%ncore = frag_ncor(i)
           end if
        end do

        ! Check whether the user defined a favorite fragment and set this to fragment 1 in case not
        if (maxval(frag_bias) < 0.001 ) then
           ! bias fragment 1 with a default value of 0.1 to have a tie-breaker for non-polar bonds
           frag_bias(1) = 0.1
        end if

   end subroutine rose_define_fragments

   subroutine rose_read_fragments (mode,fragments,ao_basis_2,spinor_fragment)
        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyzfile
        use hdf5_interface
        use fchk_interface
        use ADF_interface
        use rose_input_data

        implicit none

        integer, intent(in)                :: mode ! 0: all; 1: valence; 2: hard virtuals
        type(molecule_info),intent(in)     :: fragments(:) ! molecule information for the fragments
        type(basis_set_info_t)             :: ao_basis_2(:) ! the ao basis sets of the fragments
        type(mo_basis)                     :: spinor_fragment(:) ! MO-coefficients of the fragments
        integer                            :: ierr,i, j, nval, ncor, nhard
        character(len=256)                 :: fragment_file
 
        ! Read the basis sets and relevant spinors of all fragments
        do i = 1, n_fragments

           ! Read the basis set and the full set of spinors for each fragment
           ! TODO: make it possible to mix atomic and molecular fragments
           fragment_file = fragments(i)%name
           if (atom_fragments) write(fragment_file,'(i3.3)') nint(fragments(i)%Znumber)
           call rose_read_mos (trim(fragment_file),file_for_molfrag,ao_basis_2(i),spinor_fragment(i))

           ! All atomic coefficients are generated with the atom in the origin, translate the basis to the actual position
           if (atom_fragments) then
              do j = 1, ao_basis_2(i)%nshells
                 ao_basis_2(i)%gtos(j)%coord = fragments(i)%firstatom%x
              end do
              if (ao_basis_2(i)%has_small_component) then
                 do j = 1, ao_basis_2(i)%nshells_small
                    ao_basis_2(i)%gtos_small(j)%coord = fragments(i)%firstatom%x
                 end do
              end if
           end if

           ! Keep only the spinors that we need
           if (restricted) then
              nval = fragments(i)%nvale/2
              ncor = fragments(i)%ncore/2
           else
              nval = fragments(i)%nvale
              ncor = fragments(i)%ncore
           endif
           
           select case (mode)
           case (1)
               call reduce_mo_basis (spinor_fragment(i),nval,ncor+1)
           case (2)
               nhard = 0
               do j = nval+ncor+1, spinor_fragment(i)%nmo
                  if (spinor_fragment(i)%energy(j) > fragment_cutoff) exit
                  nhard = j - nval - ncor
               end do
               call reduce_mo_basis (spinor_fragment(i),nhard,ncor+nval+1)
           case default
           end select
        end do
     
   end subroutine rose_read_fragments

   subroutine rose_read_mos (file_name,file_type,ao_basis,spinor)

! Read MOs from MO coefficient file

        use rose_ao
        use rose_mo
        use hdf5_interface
        use fchk_interface
        use adf_interface
        use rose_input_data
        use rose_x2c, only: rose_fw_transform

        implicit none

        character(len=*),intent(in)        :: file_name,file_type ! name and type of file from which to read
        type(basis_set_info_t),intent(out) :: ao_basis      ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(out)         :: spinor        ! all spinors as read from the checkpoint file
        integer                            :: ierr,spinor_type

        spinor_type = mo_type(restricted,spatial_orbitals)

        if (trim(file_type) .eq. "dirac") then
           call read_from_hdf5 (file_name,ierr,ao_basis,spinor)
           if (ierr .ne. 0) write(*,*)"**WARNING** Reading the DIRAC file has led to an error : ierr = ",ierr

        elseif (trim(file_type) .eq. "h5") then
           call read_from_hdf5 (file_name,ierr,ao_basis,spinor)

        else if (trim(file_type) .eq. "fchk") then
           spinor%type = spinor_type
           call read_from_fchk (file_name,ao_basis,spinor)

        else if (trim(file_type) .eq. "ams") then
           call read_in_adf (ao_basis,spinor)

        else if (trim(file_type) .eq. "adf") then
           call read_from_adf (file_name,ao_basis,spinor)

        else
           write(*,*)"**WARNING** wrong file extension: ",trim(file_type)
        end if

        if (spinor%has_parent .and. spinor%parent_active) &
           call rose_fw_transform(spinor,ao_basis)

        call convert_mo_basis(spinor,spinor_type)
 
   end subroutine rose_read_mos

   subroutine init_fragment (ao_basis,spinor)

! Initialize fragment with unit vectors, but allowing for projection from a Cartesian basis to a reduced basis consisting of solid harmonic (real) Gaussians

        use rose_ao
        use rose_mo
        use rose_utils

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis       ! the ao basis
        type(mo_basis),intent(out) :: spinor               ! set of spinors corresponding to complete filling of this basis
        real(8), allocatable :: cartesian_to_spherical(:,:)
        integer ishell, iao, imo, nao, nmo, orb_momentum, nao_l, nmo_l

        ! Initialize the spinors as Kramers restricted spinors
        nao = ao_basis%nao
        nmo = 0
        do ishell = 1, ao_basis%nshells
           nmo = nmo + nfunctions(ao_basis%gtos(ishell)%orb_momentum,2)
        end do
        call alloc_mo_basis (spinor,nao,nmo,4)

        spinor%energy   = 0.D0
        spinor%coeff_q  = 0.D0
        nao = 0 ! reset to use as an offset
        nmo = 0 ! reset to use as an offset
        do ishell = 1, ao_basis%nshells
           orb_momentum =  ao_basis%gtos(ishell)%orb_momentum
           nao_l = nfunctions(orb_momentum,ao_basis%basis_angular) ! number of basis functions for this l-value
           nmo_l = nfunctions(orb_momentum,2)                      ! number of spherical functions for this l-value
           allocate (cartesian_to_spherical(nao_l,nmo_l))
           call get_spherical_projection(orb_momentum,cartesian_to_spherical)   ! get projection matrix from cartesian to sphericals
           do imo = 1, nmo_l
              do iao = 1, nao_l
                 spinor%coeff_q(nao+iao,nmo+imo,1) = cartesian_to_spherical(iao,imo)
              end do
           end do
           nao = nao + nao_l
           nmo = nmo + nmo_l
           deallocate (cartesian_to_spherical)
        end do

   end subroutine init_fragment

   subroutine get_spherical_projection (orb_momentum,cartesian_to_spherical)

      integer, intent(in) :: orb_momentum
      real(8), intent(out) :: cartesian_to_spherical(:,:)

      integer :: l, m, lm, lx, ly, lz, lxyz

      l = orb_momentum - 1
      lm = 0
      do m = -l, l
         lm = lm+1
         lxyz = 0
         do lx = l, 0, -1
         do ly = l-lx, 0, -1
            lz = l-lx-ly
            lxyz = lxyz + 1
            cartesian_to_spherical(lxyz,lm) = rcsp(l,m,lx,ly,lz)
         end do
         end do
      end do

   end subroutine get_spherical_projection

   function rcsp (l,m,lx,ly,lz) result (coeff)

      ! returns coefficients for the transformation from Cartesians to the real solid harmonics

      integer, intent(in) :: l,m,lx,ly,lz
      integer    :: ma
      complex(8) :: fp, fm
      real(8)    :: coeff

      if (lx+ly+lz /= l .or. l<m .or. l<-m) then
         coeff = 0.D0
         return
      end if

      ma = abs(m)

      call csp (l,ma,lx,ly,lz,fp,fm)

      ! Use the real combinations
      if (m == 0) then
         coeff = real(fp,8)
      else if (m > 0) then
         coeff = dsqrt(0.5D0)*real(fp+fm,8)
      else
         coeff = dsqrt(0.5D0)*dimag(fp-fm)
      end if

      ! Dirac and Dalton work with unnormalized cartesian Gaussians
      if (lx>1) coeff = coeff / dsqrt(semifactorial(2*lx-1))
      if (ly>1) coeff = coeff / dsqrt(semifactorial(2*ly-1))
      if (lz>1) coeff = coeff / dsqrt(semifactorial(2*lz-1))

   end function rcsp

   subroutine csp (l,ma,lx,ly,lz,fp,fm)

      ! Equation (15) of H.B. Schlegel, M.J. Frisch, Int J Quantum Chem. 54 (1995) 83–87. doi:10.1002/qua.560540202.
      ! returns coefficients for the transformation from Cartesians to complex spherical harmonics
      ! NB: ma is the absolute value of m

      integer, intent(in) :: l,ma,lx,ly,lz
      complex(8), intent(out) :: fp,fm

      integer    :: i, j, k
      real(8)    :: f1, f2, f3, f4
      complex(8) :: f5

      fp  = dcmplx(0.D0,0.D0)
      fm  = dcmplx(0.D0,0.D0)
      f1 = factorial(l-ma)*factorial(l)/(factorial(l+ma)*factorial(2*l))
      f1 = f1 * factorial(2*lx)*factorial(2*ly)*factorial(2*lz)/(factorial(lx)*factorial(ly)*factorial(lz))
      f1 = dsqrt(f1)/(2**l*factorial(l))
      do i = 0, (l-ma)/2
         if (mod(lx+ly-ma,2) /= 0) cycle ! coefficients are zero for half-integer j
         j = (lx + ly - ma)/2
         f2 = binomial(l,i) * binomial(i,j)
         f3 = (-1)**i*factorial(2*l-2*i)/factorial(l-ma-2*i)
         do k = 0, j
            f4 = binomial(j,k) * binomial(ma,lx-2*k)
            f5 = dcmplx(-1.D0,0.D0)**((ma-lx+2*k)/2.D0)
            fp = fp + f2*f3*f4*f5
            f5 = dcmplx(-1.D0,0.D0)**(-(ma-lx+2*k)/2.D0)
            fm = fm + f2*f3*f4*f5
         end do
      end do
      fp = fp*f1
      fm = fm*f1

   end subroutine csp

   function factorial (n) result (fac)

    integer, intent (in) :: n
    real(8) :: fac
    integer :: i

    if (n<0) then
       fac = 0.D0
       return
    else
       fac = 1.D0
    end if
    do i = 2, n
       fac = fac * real(i,8)
    end do

   end function factorial

   function semifactorial (n) result (semifac)

    integer, intent (in) :: n
    real(8) :: semifac
    integer :: k

    if (n<0) then
       semifac = 0.D0
       return
    end if

    if (mod(n,2) == 0) then
       !even argument
       k = n/2
       semifac = 2.D0**k * factorial(k)
    else
       !odd argument
       k = (n-1)/2
       semifac = factorial(n) / (factorial(k) * 2.0**k)
    end if

   end function semifactorial

   function binomial (n, k) result (bin)

    integer, intent (in) :: n
    integer, intent (in) :: k
    real(8) :: bin
    integer :: i

    if (min(k,n-k) < 0) then
       bin = 0.D0
    else
       bin = 1.D0
       do i = 1, min(k,n-k)
          bin = bin * real(n-i+1,8) / real(i,8)
       end do
    end if

   end function binomial

   subroutine atom_minbas_size(atom_number,ncore,nvale)

!        returns size of minimal basis (NB: counting spinors, so there is a factor of 2 wrt orbitals)

         implicit none

         integer,intent(in)  ::  atom_number  ! the number that is used to search the basis set and ao-coefficient libs
         integer,intent(out) ::  ncore, nvale ! the number of core and valence spinors

!        Easiest is to simply hardwire this

         select case (atom_number)
            case (:0)        ! 0 is to be used for point charges that do not carry a basis set
                 ncore = 0
                 nvale = 0
            case (1:2)       ! H, He
                 ncore = 0
                 nvale = 2
            case (3:4)       ! Li,Be
                 ncore = 2
                 nvale = 2   ! better changed into 8, but keeping for the moment to let tests pass
            case (5:10)      ! B..Ne
                 ncore = 2
                 nvale = 8
            case (11:12)     ! Na,Mg
                 ncore = 10
                 nvale = 2
            case (13:18)     ! Al..Ar
                 ncore = 10
                 nvale = 8
            case (19:20)     ! K,Ca
                 ncore = 18
                 nvale = 2
            case (21:30)     ! Sc..Zn
                 ncore = 18
                 nvale = 12
            case (31:36)     ! Ga..Kr
                 ncore = 18
                 nvale = 18
            case (37:38)     ! Rb,Sr
                 ncore = 36
                 nvale = 2
            case (39:48)     ! Y..Cd
                 ncore = 36
                 nvale = 12
            case (49:54)     ! In..Xe
                 ncore = 36
                 nvale = 18
            case (55:56)     ! Cs,Ba
                 ncore = 54
                 nvale = 2
            case (57:70)     ! La..Yb
                 ncore = 54
                 nvale = 16  ! assuming that we will enforce 4f rather than 5d occupation in the atomic scf
            case (71:80)     ! Lu..Hg
                 ncore = 54
                 nvale = 26
            case (81:86)     ! Tl..Rn
                 ncore = 68  ! place the 4f in the core
                 nvale = 18
            case (87:88)     ! Fr,Ra
                 ncore = 86 
                 nvale = 2
            case (89:102)    ! Ac,RNo
                 ncore = 86 
                 nvale = 16  ! assuming that we will enforce 5f rather than 6d occupation in the atomic scf
            case (103:112)   ! Lr..Cn
                 ncore = 86
                 nvale = 26
            case (113:118)   ! Nh..Og
                 ncore = 100 ! place the 5f in the core
                 nvale = 18
            case (119:)      ! Venturing into the unknown superheavies, impossible to decide the valence
                 ncore = 118
                 nvale = atom_number - ncore
         end select
     
   end subroutine atom_minbas_size

   subroutine set_coreval_info (molecule,include_core)

! Loop over atoms and define what is core and valence. Works also for a single atom fragment.
! Always counting # of electrons, so needs to be divided by two to get restricted orbitals

        use read_xyzfile, only : molecule_info
        use read_xyzfile, only: atom

        implicit none

        type(molecule_info),intent(inout) :: molecule      ! the molecule
        logical, intent(in)               :: include_core  ! this implies no core, so all electrons are made valence

        type(atom),pointer :: current_atom
        integer            :: ncore, nvale

        molecule%ncore = 0
        molecule%nvale = 0
        current_atom => molecule%firstatom
        do while(associated(current_atom))
           call atom_minbas_size (nint(current_atom%bascharge),ncore,nvale)
           if (include_core) then
              nvale = ncore + nvale
              ncore = 0
           end if
           molecule%ncore = molecule%ncore + ncore
           molecule%nvale = molecule%nvale + nvale
           current_atom => current_atom%next
        enddo

    end subroutine set_coreval_info

end module
