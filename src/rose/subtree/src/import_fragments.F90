module import_fragments

implicit none
private

public get_atomic_fragment,init_fragment,get_molecular_fragment,atom_minbas_size, set_coreval_info

contains


   subroutine get_atomic_fragment (atom_number,atom_coord,ao_basis,spinor,spinor_type,file_for_MOs)

! Read atomic fragment from MO coefficient file

        use read_xyz, only : molecule_info
        use rose_ao
        use rose_mo
        use hdf5_interface
        use fchk_interface
        use adf_interface

        implicit none

        integer, intent(in)                :: atom_number   ! needed to find the file that is to be read
        real(8), intent(in)                :: atom_coord(3) ! the basis function center will be moved to this position
        integer, intent(in)                :: spinor_type   ! desired spinor type (may be different from the standard on the checkpoint file)
        type(basis_set_info_t),intent(out) :: ao_basis      ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(out)         :: spinor        ! all spinors as read from the checkpoint file
        character(len=10)                  :: fragfile
        integer                            :: ierr, ishell
        character(len=*),intent(in)        :: file_for_MOs
        logical                            :: restricted
        type(molecule_info)                :: frag          ! geometric information that is read (not used for these atomic fragments)

        if (spinor_type == 1 .or. spinor_type == 4) then
           restricted = .true.
        else
           restricted = .false.
        end if

        ! Should include path information, for now assume this is copied into the working dir
        write (frag%name,'(i3.3)') atom_number
        if (trim(file_for_MOs) .eq. "dirac") then
           call read_from_hdf5 (trim(frag%name),ierr,ao_basis,spinor,frag)
           if (ierr .ne. 0) write(*,*)"**WARNING** Reading the DIRAC file has led to an error : ierr = ",ierr
           call convert_mo_basis(spinor,spinor_type)

        elseif (trim(file_for_MOs) .eq. "h5") then
           call read_from_hdf5 (trim(frag%name),ierr,ao_basis,spinor,frag)

        else if (trim(file_for_MOs) .eq. "fchk") then
           call read_from_fchk (trim(frag%name),ao_basis,spinor,restricted,frag)

        else if (trim(file_for_MOs) .eq. "ams") then
           call read_in_adf (ao_basis,spinor,frag)

        else if (trim(file_for_MOs) .eq. "adf") then
           call read_from_adf (trim(frag%name),ao_basis,spinor,frag)

        else
           write(*,*)"**WARNING** wrong file extension: ",trim(file_for_MOs)
        end if

        ! All atomic coefficients are generated with the atom in the origin, translate the basis to the actual position
        do ishell = 1, ao_basis%nshells ! somehow the compiler wants an explicit loop for this
            ao_basis%gtos(ishell)%coord = atom_coord
        end do

   end subroutine get_atomic_fragment

   subroutine get_molecular_fragment (frag,ao_basis,spinor,spinor_type,file_for_MOs)

        use read_xyz, only : molecule_info
        use rose_ao
        use rose_mo
        use hdf5_interface
        use fchk_interface
        use adf_interface

        implicit none

        type(molecule_info),intent(inout)  :: frag          ! molecule information read from the checkpoint file
        type(basis_set_info_t),intent(inout) :: ao_basis      ! the ao basis used in the generation of the MO coefficients
        type(mo_basis),intent(inout)       :: spinor        ! all spinors as read from the checkpoint file
        integer, intent(in)                :: spinor_type   ! desired spinor type (may be different from the standard on the checkpoint file)
        integer                            :: ierr
        character(len=*),intent(in)        :: file_for_MOs
        logical                            :: restricted

        ! Should include path information, for now assume this is copied into the working dir

        if (spinor_type == 1 .or. spinor_type == 4) then
           restricted = .true.
        else
           restricted = .false.
        end if
        
        if (trim(file_for_MOs) .eq. "dirac") then
           call read_from_hdf5 (frag%name,ierr,ao_basis,spinor,frag)
           if (ierr .ne. 0) write(*,*)"**WARNING** Reading the DIRAC file has led to an error : ierr = ",ierr
           call convert_mo_basis(spinor,spinor_type)

        elseif (trim(file_for_MOs) .eq. "h5") then
           call read_from_hdf5 (frag%name,ierr,ao_basis,spinor,frag)

        else if (trim(file_for_MOs) .eq. "fchk") then
           call read_from_fchk (frag%name,ao_basis,spinor,restricted,frag)

        else if (trim(file_for_MOs) .eq. "ams") then
           call read_in_adf (ao_basis,spinor,frag)

        else if (trim(file_for_MOs) .eq. "adf") then
           call read_from_adf (frag%name,ao_basis,spinor,frag)

        else
           write(*,*)"**WARNING** wrong file extension: ",trim(file_for_MOs)

        end if

   end subroutine get_molecular_fragment

   subroutine init_fragment (ao_basis,spinor)

! Initialize fragment with unit vectors, but allowing for projection from a Cartesian basis to a reduced basis consisting of solid harmonic (real) Gaussians

        use rose_ao
        use rose_mo
        use rose_utils

        implicit none

        type(basis_set_info_t),intent(in) :: ao_basis       ! the ao basis
        type(mo_basis),intent(out) :: spinor               ! set of spinors corresponding to complete filling of this basis
        real(8), allocatable :: cartesian_to_spherical(:,:)
        integer i, ishell, iao, imo, nao, nmo, orb_momentum, nao_l, nmo_l

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
      real(8) :: value

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
                 nvale = 2
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

        use read_xyz, only : atom,molecule_info

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
