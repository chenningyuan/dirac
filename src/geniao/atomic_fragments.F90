module atomic_fragments

contains


   subroutine get_atomic_fragment (atom_number,atom_coord,ao_basis,cspinor)

        use exacorr_mo
        use exacorr_datatypes

        implicit none

        integer, intent(in) :: atom_number
        real(8), intent(in) :: atom_coord(3)
        type(basis_set_info_t),intent(out) :: ao_basis ! the ao basis used in the generation of the MO coefficients
        type(cmo),intent(out) :: cspinor               ! occupied spinors (minimal basis set)

        type(cmo):: cspinor_full  ! all spinors as they are read from file
        character(len=10) fragfile
        integer ierr, ishell, nao, nmo, ncore, nvale

        ! Should include path information, for now assume this is copied into the working dir
        write (fragfile,'(i3.3,a7)') atom_number,'.dfcoef'
        call read_from_dirac (ao_basis,cspinor_full,fragfile,ierr)

        if (ierr .ne. 0) then
          write(*,*)"**WARNING** Reading the DFCOEF file has led to an error : ierr = ",ierr
        endif
        ! All atomic coefficients are generated with the atom in the origin, translate the basis to the actual position
        do ishell = 1, ao_basis%nshells ! somehow the compiler wants an explicit loop for this
            ao_basis%gtos(ishell)%coord = atom_coord
        end do
        ! Query the size of the (mo) basis that should be made
        call atom_minbas_size (atom_number,ncore,nvale)
        nmo = ncore+nvale 
        nao = cspinor_full%nao
        ! Copy only the relevant MOs to the output array
        call alloc_mo (cspinor,nao,nmo)
        cspinor%energy = cspinor_full%energy(1:nmo)
        cspinor%coeff  = cspinor_full%coeff(:,1:nmo,:)
        cspinor%index  = cspinor_full%index(1:nmo)
        ! we don't need this coefficient set anymore, free up the memory
        call dealloc_mo(cspinor_full)

   end subroutine get_atomic_fragment

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

end module
