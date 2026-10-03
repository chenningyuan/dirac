module rose_import_fragments

implicit none
private

public rose_define_fragments, rose_read_fragments, rose_read_mos, init_fragment

contains

   subroutine rose_define_fragments (molecule,fragments,ao_basis)

        use rose_input_data
        use read_xyzfile
        use read_xyzfile, only: atom
        use rose_ao, only: basis_set_info_t, collapse_ao_centers, match_center, &
                           shell_coord, shell_orb_momentum, set_shell_coord
        use rose_global, only: luout

        implicit none

        type(molecule_info)                :: molecule ! molecule information : name, charge, nelec, multiplicity, coordinates, ..
        type(molecule_info), allocatable   :: fragments(:) ! same information for the fragment
        type(basis_set_info_t), intent(in) :: ao_basis ! supermolecular AO basis; its shell centers define the atomic fragments
        type(atom), pointer                :: current_atom ! linked list to store coordinate information
        integer                            :: i, icen
        integer                            :: n_centers    ! number of distinct AO centers (centers carrying a basis)
        integer, allocatable               :: center_owner(:)            ! center -> fragment that claimed it (0 = unclaimed)
        integer                            :: n_lab, ilab, ig, jg, ifrag  ! labelled-xyz bookkeeping
        real(8), allocatable               :: label_coord(:,:)           ! join-key coordinates from the labelled xyz
        character(len=label_length), allocatable :: labels(:)            ! one group label per line of that file
        character(len=label_length), allocatable :: group_name(:)        ! distinct labels, in order of first appearance
        real(8), allocatable               :: center_coord(:,:)          ! coordinates of those centers
        integer, allocatable               :: shell_to_center(:)         ! large-component shell -> center
        integer, allocatable               :: shell_to_center_small(:)   ! small-component shell -> center
        integer, external                  :: adf_n_atomic_centers       ! AMS (Slater) center count

        ! n_fragments == 0 signals that we use atomic fragments, otherwise molecular fragments indicated as frag1, frag2, etc.
        if (n_fragments == 0) then
            if (.not. file_for_MOs .eq. "ams") atom_fragments = .true.
            if (atom_fragments .and. ao_basis%basis_type == 2 .and. &
                .not. associated(ao_basis%stos)) then
               ! Slater basis inside AMS, where the basis is an ADF object and
               ! stos(:) is null, so there are no shell coordinates to collapse and
               ! ADF is asked for the centers instead. A Slater basis read from file
               ! carries stos(:) and goes through collapse_ao_centers like any other.
               ! The size is asked for first and the array allocated here, because
               ! adf_atomic_centers is an external procedure without an explicit
               ! interface and so cannot take an allocatable dummy argument.
               n_centers = adf_n_atomic_centers(ao_basis)
               allocate(center_coord(3,n_centers))
               call adf_atomic_centers (ao_basis,n_centers,center_coord)
               allocate(shell_to_center(0),shell_to_center_small(0))
               n_fragments = n_centers
            else if (atom_fragments) then
               ! The number of fragments is the number of centers that actually carry a
               ! basis, taken from the AO basis itself rather than from the nucleus count.
               ! For a conventional molecule the two are equal; they differ exactly where
               ! ghosts or point nuclei are present. The ordering is fixed further below.
               call collapse_ao_centers (ao_basis,center_match_tol,n_centers,center_coord, &
                                         shell_to_center,shell_to_center_small)
               n_fragments = n_centers
            else
               ! The old ADF pathway supplies its own per-atom fragment files.
               n_fragments = molecule%natoms
            end if
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
            ! The nucleus list drives the ordering and the AO basis drives membership.
            !
            ! Ordering has to come from the nucleus list because add_atom groups that
            ! list by element (a new atom is inserted after the last one of the same
            ! name), and the AO shell ordering does not reproduce that grouping. Walking
            ! the nucleus list therefore leaves the fragment numbering exactly as it was.
            !
            ! Membership comes from the AO basis: an entry carrying no basis (a point
            ! nucleus) has no center and is skipped, and every AO center must be claimed
            ! by exactly one entry.
           allocate(center_owner(n_centers))
           center_owner = 0
           i = 0
           current_atom => molecule%firstatom
           do while (associated(current_atom))
              icen = match_center(current_atom%x,center_coord,n_centers,center_match_tol)
              if (icen > 0) then
                 if (center_owner(icen) /= 0) then
                    write(*,'(A,3F14.8)') ' Two entries of the molecule list share the AO center at ', &
                       center_coord(1:3,icen)
                    stop 'rose_define_fragments: center_match_tol is too loose for this geometry'
                 end if
                 if (current_atom%bascharge <= charge_match_tol) then
                    write(*,'(A,3F14.8)') ' The molecule list declares no basis at ',current_atom%x
                    write(*,'(A)') ' but the AO basis has shells on that center. The geometry and the basis'
                    write(*,'(A)') ' in the checkpoint file are inconsistent.'
                    stop 'rose_define_fragments: molecule list and AO basis disagree about this center'
                 end if
                 i = i + 1
                 center_owner(icen) = i
                 fragments(i)%natoms  = 1
                 if (abs(current_atom%charge) < ghost_threshold .and. current_atom%bascharge > ghost_threshold) then
                    fragments(i)%name = 'Gh(' // trim(current_atom%name) // ')'
                 else
                    fragments(i)%name = current_atom%name
                 end if
                 fragments(i)%Znumber = current_atom%charge
                 allocate(fragments(i)%firstatom)
                 fragments(i)%firstatom = current_atom
                 ! Place the fragment on its AO center rather than on the nucleus, so that
                 ! the translation of the atomic library basis in rose_read_fragments uses
                 ! the coordinate the supermolecular orbitals were actually built on.
                 fragments(i)%firstatom%x = center_coord(1:3,icen)
                 nullify(fragments(i)%firstatom%next)
              end if
              current_atom => current_atom%next
           end do

           ! Every center carrying a basis must have been claimed. One that was not has no
           ! entry in the molecule list, so its element -- and hence its reference atomic
           ! basis -- cannot be determined, while its functions still span part of the AO
           ! space. Continuing would silently leave them in no fragment at all.
           do icen = 1, n_centers
              if (center_owner(icen) == 0) then
                 write(*,'(A,3F14.8)') ' The AO basis has shells on a center at ',center_coord(1:3,icen)
                 write(*,'(A)') ' that has no entry in the molecule list, so its element is unknown and no'
                 write(*,'(A)') ' reference atomic basis can be found for it. Add it to the geometry as a'
                 write(*,'(A)') ' ghost center (zero nuclear charge, nonzero basis charge), or use molecular fragments.'
                 stop 'rose_define_fragments: AO center without a corresponding entry in the molecule list'
              end if
           end do

           ! ---------------------------------------------------------------------
           ! Optional: read the recanonization groups from a labelled xyz file.
           !
           ! The coordinates in that file are a join key only. Each line is matched
           ! against the AO centers found above and then discarded, so the file can
           ! never disagree with the orbitals about where a center is; it can only
           ! fail to match, which is an error. The label decides which group the
           ! fragment on that center belongs to, numbered by first appearance.
           ! ---------------------------------------------------------------------
           if (len_trim(fragment_xyz_file) > 0) then

              call read_xyz_labels (trim(fragment_xyz_file),n_lab,label_coord,labels)

              allocate(group_name(n_lab))
              allocate(frag_to_group(n_fragments))
              frag_to_group = 0
              n_groups      = 0

              do ilab = 1, n_lab
                 icen = match_center(label_coord(1:3,ilab),center_coord,n_centers,center_match_tol)
                 if (icen == 0) then
                    write(*,'(A,I0,A)') ' Line ',ilab+2,' of '//trim(fragment_xyz_file)//' is at'
                    write(*,'(A,3F14.8)') '   ',label_coord(1:3,ilab)
                    write(*,'(A)') ' which matches no center of the AO basis. The labelled file and the'
                    write(*,'(A)') ' orbitals describe different geometries, or center_match_tol is too tight.'
                    stop 'rose_define_fragments: labelled xyz line matches no AO center'
                 end if

                 ifrag = center_owner(icen)
                 if (frag_to_group(ifrag) /= 0) then
                    write(*,'(A,3F14.8)') ' Two lines of '//trim(fragment_xyz_file)//' claim the center at ', &
                       center_coord(1:3,icen)
                    stop 'rose_define_fragments: duplicate center in the labelled xyz file'
                 end if

                 ! Find this label among those seen so far, or open a new group for it.
                 ig = 0
                 do jg = 1, n_groups
                    if (trim(labels(ilab)) == trim(group_name(jg))) then
                       ig = jg
                       exit
                    end if
                 end do
                 if (ig == 0) then
                    n_groups = n_groups + 1
                    group_name(n_groups) = labels(ilab)
                    ig = n_groups
                 end if

                 frag_to_group(ifrag) = ig
              end do

              ! Every fragment must have been given a group.
              do i = 1, n_fragments
                 if (frag_to_group(i) == 0) then
                    write(*,'(A,I0,A)') ' Fragment ',i,' ('//trim(fragments(i)%name)//') at'
                    write(*,'(A,3F14.8)') '   ',fragments(i)%firstatom%x
                    write(*,'(A)') ' was not named in '//trim(fragment_xyz_file)//'.'
                    write(*,'(A)') ' Every center carrying a basis needs a line with a group label.'
                    stop 'rose_define_fragments: fragment left without a group'
                 end if
              end do

              ! Report the grouping, so it is visible in the output.
              write(luout,*)
              write(luout,*) "####################################"
              write(luout,*) "#     Recanonization  Groups       #"
              write(luout,*) "####################################"
              write(luout,*)
              write(luout,'(A,I0,A,I0,A)') ' ',n_fragments,' IAO fragments collected into ',n_groups,' groups'
              do ig = 1, n_groups
                 write(luout,'(A,I3,A,A)') ' group ',ig,' : ',trim(group_name(ig))
                 do i = 1, n_fragments
                    if (frag_to_group(i) == ig) write(luout,'(A,I4,A,A)') '      fragment ',i,'  ',trim(fragments(i)%name)
                 end do
              end do
              write(luout,*)

              deallocate(label_coord,labels,group_name)
           end if

           deallocate(center_owner,center_coord,shell_to_center,shell_to_center_small)
        else
           if (len_trim(fragment_xyz_file) > 0) then
              write(*,'(A)') ' .FRAGMENT_XYZ is currently supported only with atomic fragments,'
              write(*,'(A)') ' where the groups are built from the centers of the AO basis.'
              write(*,'(A)') ' Grouping molecular fragments is not implemented; it would be done by'
              write(*,'(A)') ' extending this keyword to that case rather than by a separate one.'
              stop 'rose_define_fragments: .FRAGMENT_XYZ with molecular fragments is not supported'
           end if

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

        ! Initialise recanonisation groups now that n_fragments is final.
        ! Installs the identity mapping when no grouping was requested.
        call rose_init_groups

   end subroutine rose_define_fragments

   subroutine rose_read_fragments (mode,fragments,ao_basis_2,spinor_fragment)
        use rose_ao
        use rose_mo
        use rose_utils
        use read_xyzfile
        use hdf5_interface
        use fchk_interface
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
           if (atom_fragments) write(fragment_file,'(i3.3)') nint(fragments(i)%firstatom%bascharge)
           call rose_read_mos (trim(fragment_file),file_for_molfrag,ao_basis_2(i),spinor_fragment(i))

           ! All atomic coefficients are generated with the atom in the origin, translate the basis to the actual position
           if (atom_fragments .and. ao_basis_2(i)%basis_type == 2 .and. &
               .not. associated(ao_basis_2(i)%stos)) then
              ! Slater basis inside AMS: the basis is an ADF object and
              ! SlaterOverlapMatrix takes the geometry alongside it, so the fragment
              ! is moved by setting the coordinates on that object. A Slater basis
              ! read from file carries stos(:) and is translated shell by shell
              ! below, exactly as a Gaussian one is.
              call adf_place_fragment (ao_basis_2(i),fragments(i)%firstatom%x)
           else if (atom_fragments) then
              do j = 1, ao_basis_2(i)%nshells
                 call set_shell_coord (ao_basis_2(i),j,fragments(i)%firstatom%x)
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
           call read_basis_in_adf (file_name,ao_basis,spinor)

        else
           write(*,*)"**WARNING** wrong file extension: ",trim(file_type)
        end if

        if (spinor%has_parent .and. spinor%parent_active) &
           call rose_fw_transform(spinor,ao_basis)

        call convert_mo_basis(spinor,spinor_type)
 
   end subroutine rose_read_mos

   subroutine init_fragment (ao_basis,spinor)

! Initialize fragment with unit vectors, but allowing for projection from a Cartesian
! basis to a reduced basis of solid harmonics. The reduction is the same one ADF and
! DIRAC both perform before diagonalising, so it applies to Slater functions as well
! as to Gaussians and the shell data is reached through shell_orb_momentum.

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
           nmo = nmo + nfunctions(shell_orb_momentum(ao_basis,ishell),2)
        end do
        call alloc_mo_basis (spinor,nao,nmo,4)

        spinor%energy   = 0.D0
        spinor%coeff_q  = 0.D0
        nao = 0 ! reset to use as an offset
        nmo = 0 ! reset to use as an offset
        do ishell = 1, ao_basis%nshells
           orb_momentum =  shell_orb_momentum(ao_basis,ishell)
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

        use read_xyzfile, only : molecule_info, ghost_threshold
        use read_xyzfile, only: atom
        use rose_input_data, only : file_for_MOs

        implicit none

        type(molecule_info),intent(inout) :: molecule      ! the molecule
        logical, intent(in)               :: include_core  ! this implies no core, so all electrons are made valence

        type(atom),pointer :: current_atom
        integer            :: ncore, nvale
        integer, external  :: adf_frozen_core_capacity

        molecule%ncore = 0
        molecule%nvale = 0
        current_atom => molecule%firstatom
        do while(associated(current_atom))
           if (abs(current_atom%charge) < ghost_threshold .and. current_atom%bascharge > ghost_threshold) then
              ncore = 0
              nvale = 0
           else
              call atom_minbas_size (nint(current_atom%bascharge),ncore,nvale)
           end if
           if (trim(file_for_MOs) .eq. "ams") then
              ! Follow the frozen core that AMS has chosen rather than imposing our own.
              ! ADF leaves the frozen orbitals out of the fragment MO list altogether, so
              ! there is nothing for ROSE to skip (ncore = 0) and the valence reference is
              ! what remains of the minimal set once the frozen shells are taken off it.
              ! With Core None the capacity is zero and this reduces exactly to the
              ! include_core branch below, which is what the ADF pathway did until now.
              nvale = ncore + nvale - adf_frozen_core_capacity(nint(current_atom%bascharge))
              ncore = 0
              if (nvale < 0) then
                 write(*,'(A,I0)') ' The frozen core of AMS leaves no valence shells for element ', &
                    nint(current_atom%bascharge)
                 stop 'set_coreval_info: frozen core larger than the minimal basis'
              end if
           else if (include_core) then
              nvale = ncore + nvale
              ncore = 0
           end if
           molecule%ncore = molecule%ncore + ncore
           molecule%nvale = molecule%nvale + nvale
           current_atom => current_atom%next
        enddo

    end subroutine set_coreval_info

end module
