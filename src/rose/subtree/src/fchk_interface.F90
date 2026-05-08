module fchk_interface

!
!    Interface with formatted checkpoint files (fchk)
!
      
contains

      subroutine reorder_spinor_fchk(spinors,ao_basis,back)

       ! As shown in /psi4/src/psi4/libmints/writer.cc, the FCHK file has a different ordering than PSI4/DIRAC/PYSCF (which have the same between themselves).
       ! One has to sort the functions to work with INTEREST Library, which also has the same ordering as PSI4/DIRAC/PYSCF.

       use rose_mo
       use rose_ao
       use rose_utils
       use read_xyz
       implicit none
       type(mo_basis), intent(inout)         :: spinors
       type(basis_set_info_t), intent(inout) :: ao_basis
       logical, intent(in)                   :: back ! if true, it reorder back to normal, i.e. is will take the transpose and inverse of normalization factors.
       integer                               :: offset, nmo, k
       real(8)                               :: cartD(6,6), cartF(10,10), cartG(15,15)
       real(8)                               :: pf1,pf2,pf3,pf4
       real(8), allocatable                  :: coeff_ordered(:,:), coeff_to_order(:,:)
       
       ! For d functions:
       pf1 = 1.D0             ! aa
       pf2 = dsqrt(1.D0/3.D0) ! ab
       if (back) pf2 = 1.D0/pf2
       ! Psi4/DIRAC:  XX   XY   XZ   YY   YZ   ZZ
       !       FCHK:  XX   YY   ZZ   XY   XZ   YZ
                    ! pf1  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  pf1  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  pf1
                    ! 0.0  pf2  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  pf2  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  pf2  0.0
       cartD = 0.D0
       cartD(1,1) = pf1
       cartD(2,4) = pf1
       cartD(3,6) = pf1
       cartD(4,2) = pf2
       cartD(5,3) = pf2
       cartD(6,5) = pf2
       if (back) cartD = transpose(cartD)
       
       ! For f functions:
       pf1 = 1.D0              ! aaa
       pf2 = dsqrt(1.D0/5.D0)  ! aab
       pf3 = dsqrt(1.D0/15.D0) ! abc
       if (back) then
          pf2 = 1.D0/pf2
          pf3 = 1.D0/pf3
       endif
       ! Psi4/DIRAC:  XXX  XXY  XXZ  XYY  XYZ  XZZ  YYY  YYZ  YZZ  ZZZ
       !       FCHK:  XXX  YYY  ZZZ  XYY  XXY  XXZ  XZZ  YZZ  YYZ  XYZ
                    ! pf1  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  pf1  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf1
                    ! 0.0  0.0  0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  pf2  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  pf3  0.0  0.0  0.0  0.0  0.0
       cartF = 0.D0
       cartF(1,1) = pf1
       cartF(2,7) = pf1
       cartF(3,10)= pf1
       cartF(4,4) = pf2
       cartF(5,2) = pf2
       cartF(6,3) = pf2
       cartF(7,6) = pf2
       cartF(8,9) = pf2
       cartF(9,8) = pf2
       cartF(10,5)= pf3
       if (back) cartF = transpose(cartF)
       
       ! For g functions:
       pf1 = 1.D0              ! aaaa
       pf2 = dsqrt(1.D0/7.D0)  ! aaab
       pf3 = dsqrt(3.D0/35.D0) ! aabb
       pf4 = dsqrt(1.D0/35.D0) ! abcc
       if (back) then
          pf2 = 1.D0/pf2
          pf3 = 1.D0/pf3
          pf4 = 1.D0/pf4
       endif
       ! Psi4/DIRAC: XXXX XXXY XXXZ XXYY XXYZ XXZZ XYYY XYYZ XYZZ XZZZ YYYY YYYZ YYZZ YZZZ ZZZZ
       !       FCHK: XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ
                    ! pf1  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf1  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf1
                    ! 0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf2  0.0
                    ! 0.0  0.0  0.0  pf3  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  pf3  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf3  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  pf4  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf4  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                    ! 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  pf4  0.0  0.0  0.0  0.0  0.0  0.0
       cartG = 0.D0
       cartG(1,1)  = pf1
       cartG(2,11) = pf1
       cartG(3,15) = pf1
       cartG(4,2)  = pf2
       cartG(5,3)  = pf2
       cartG(6,7)  = pf2
       cartG(7,12) = pf2
       cartG(8,10) = pf2
       cartG(9,14) = pf2
       cartG(10,4) = pf3
       cartG(11,6) = pf3
       cartG(12,13)= pf3
       cartG(13,5) = pf4
       cartG(14,8) = pf4
       cartG(15,9) = pf4
       if (back) cartG = transpose(cartG)

       offset = 0
       nmo = spinors%nmo
       do k=1,ao_basis%nshells
         if (ao_basis%gtos(k)%orb_momentum .eq. 1) then
            offset = offset + 1
         else if (ao_basis%gtos(k)%orb_momentum .eq. 2) then
            offset = offset + 3
         else if (ao_basis%gtos(k)%orb_momentum .eq. 3) then
            allocate(coeff_ordered(6,nmo))
            allocate(coeff_to_order(6,nmo))
            coeff_to_order = spinors%coeff_r(offset+1:offset+6,:,1)
            coeff_ordered = matmul(cartD,coeff_to_order)
            spinors%coeff_r(offset+1:offset+6,:,1) = coeff_ordered
            deallocate(coeff_ordered,coeff_to_order)
            offset = offset + 6
         else if (ao_basis%gtos(k)%orb_momentum .eq. 4) then
            allocate(coeff_ordered(10,nmo))
            allocate(coeff_to_order(10,nmo))
            coeff_to_order = spinors%coeff_r(offset+1:offset+10,:,1)
            coeff_ordered = matmul(cartF,coeff_to_order)
            spinors%coeff_r(offset+1:offset+10,:,1) = coeff_ordered
            offset = offset + 10
            deallocate(coeff_ordered,coeff_to_order)
         else if (ao_basis%gtos(k)%orb_momentum .eq. 5) then
            allocate(coeff_ordered(15,nmo))
            allocate(coeff_to_order(15,nmo))
            coeff_to_order = spinors%coeff_r(offset+1:offset+15,:,1)
            coeff_ordered = matmul(cartG,coeff_to_order)
            spinors%coeff_r(offset+1:offset+15,:,1) = coeff_ordered
            offset = offset + 15
            deallocate(coeff_ordered,coeff_to_order)
         else
            stop "Cannot handle more than G functions to write FCHK file."
         end if
       enddo

      end subroutine reorder_spinor_fchk

      subroutine read_from_fchk (fname,ao_basis,spinors,molecule)

       use rose_mo
       use rose_ao
       use rose_utils
       use read_xyz

       implicit none

       character(len=*),intent(in)                     :: fname
       type(mo_basis), intent(inout), optional         :: spinors
       type(basis_set_info_t), intent(inout), optional :: ao_basis
       type(molecule_info),    intent(inout), optional :: molecule
       logical                :: restricted
       logical                :: tobe
       character(len=100)     :: ACHAR
       real(8)                :: total_energy,core_energy
       real(8)                :: coeff_norm
       integer                :: lucoef
       integer                :: i, j, k, ij
       integer                :: START, FIN, BUFF, DIFF
       integer                :: tmp, lineskip
       integer                :: nao, nmo
       integer                :: natoms, ncoord, ncoord_shells
       integer                :: nshells, nshells_cont
       integer                :: ncoeff, pureamd, pureamf
       integer, allocatable   :: atomic_numbers(:)
       real(8), allocatable   :: coord(:,:)
       integer, allocatable   :: orbital_momentum(:)
       integer, allocatable   :: shell_to_atom_map(:)
       real(8), allocatable   :: MOcoeff_alpha(:),MOcoeff_beta(:),one_elec_int(:)
       real(8), allocatable   :: MOenergy_alpha(:),MOenergy_beta(:)
       real(8), allocatable   :: all_exponents(:)
       real(8), allocatable   :: coord_shells(:)
       character(len=12)      :: fmt_tmp
       character(len=10)      :: format_size_single
       character(len=8)       :: single_precision
       character(len=20)      :: file_name

       format_size_single = '(A1,I1,A7)'
       single_precision   = 'ES16.8)'

       file_name = trim(fname)//".fchk"
       inquire (file=file_name,exist=tobe)

       if (.not.tobe) then
          print*,"error: MO coefficient file"//file_name//" not found"
          stop
       end if

       call get_free_fileunit(lucoef)
       open (lucoef,file=file_name,status='OLD',FORM='FORMATTED',access='SEQUENTIAL')

       ! read first comment line.
       read (lucoef,'(A100)') ACHAR
       ! read number of atoms
       read (lucoef,'(A49,I12)') ACHAR,natoms
       molecule%natoms = natoms
       read (lucoef,'(A49,I12)') ACHAR,molecule%charge
       read (lucoef,'(A49,I12)') ACHAR,molecule%multiplicity
       read (lucoef,'(A49,I12)') ACHAR,molecule%nelec
       read (lucoef,'(A49,I12)') ACHAR,molecule%nalpha
       read (lucoef,'(A49,I12)') ACHAR,molecule%nbeta
       ! get dimensions of the number of basis functions
       read(lucoef,'(A49,I12)') ACHAR,nao ! number of basis functions
       ! read atomic numbers
       read(lucoef,*) ! same as number of atoms
       START = 1
       FIN   = natoms
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(atomic_numbers(natoms))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          read(lucoef,trim(fmt_tmp)) atomic_numbers(START:START+DIFF-1)
          START = START + 6
       end do
       molecule%znumber = sum(atomic_numbers)

       ! skip nuclear charges (same as atomic number but real)
       read(lucoef,'(A49,I12)') ACHAR, tmp
       lineskip = ceiling(real(tmp)/5.D0)
       do i=1,lineskip
          read(lucoef,*)
       end do

       ! read the coordinates
       allocate(coord(3,natoms))
       read(lucoef,'(A49,I12)') ACHAR, tmp
       read(lucoef,*) coord

       ! process coordinates and atomic_numbers
       call get_atoms (natoms,coord,molecule%firstatom,atomic_numbers)

       if (.not.present(ao_basis).and..not.present(spinors)) return

       ! read the number of primitive shells
       read(lucoef,'(A49,I12)') ACHAR,nshells

       ! read if Cartesian functions are used
       read(lucoef,'(A49,I12)') ACHAR,pureamd
       read(lucoef,'(A49,I12)') ACHAR,pureamf
       if (pureamd .eq. 0 .or. pureamf .eq. 0) print*,"**Warning** cartesian functions should be used, not spherical ones."

       ! get the shell types (orbital momentum)
       read(lucoef,*) ! skip line that gives the number of shells.
       START = 1
       FIN   = nshells
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(orbital_momentum(nshells))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          read(lucoef,trim(fmt_tmp)) orbital_momentum(START:START+DIFF-1)
          START = START + 6
       end do

       ! get the number of primitives per shell (to skip)
       read(lucoef,'(A49,I12)') ACHAR, tmp
       lineskip = ceiling(real(tmp)/6.D0)
       do i=1,lineskip
          read(lucoef,*)
       end do

       ! get the shell to atom map
       read(lucoef,*) ! this line gives the number of shells which I already have
       START = 1
       FIN   = nshells
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(shell_to_atom_map(nshells))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          read(lucoef,trim(fmt_tmp)) shell_to_atom_map(START:START+DIFF-1)
          START = START + 6
       end do

       ! read primitive exponents
       read(lucoef,*) ! this line gives the number of shells which I already have
       START = 1
       FIN   = nshells
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(all_exponents(nshells))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          read(lucoef,trim(fmt_tmp)) all_exponents(START:START+DIFF-1)
          START = START + 5
       end do

       ! get the number of contraction coefficients (always equal to 1 for decontracted basis set)
       read(lucoef,'(A49,I12)') ACHAR, tmp
       lineskip = ceiling(real(tmp)/5.D0)
       do i=1,lineskip
          read(lucoef,*)
       end do

       ! get the coordinates of each shell
       read(lucoef,'(A49,I12)') ACHAR,ncoord_shells
       START = 1
       FIN   = ncoord_shells
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(coord_shells(ncoord_shells))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          read(lucoef,trim(fmt_tmp)) coord_shells(START:START+DIFF-1)
          START = START + 5
       end do

       ! read total energy
       read(lucoef,'(A49,ES22.15)') ACHAR,total_energy
       ! read alpha orbital energies
       read(lucoef,'(A49,I12)') ACHAR,nmo ! number of MOs
       START = 1
       FIN   = nmo
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(MOenergy_alpha(nmo))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          read(lucoef,trim(fmt_tmp)) MOenergy_alpha(START:START+DIFF-1)
          START = START + 5
       end do

       ! read alpha orbital coefficients
       read(lucoef,'(A49,I12)') ACHAR, ncoeff
       START = 1
       FIN   = ncoeff
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       allocate(MOcoeff_alpha(ncoeff))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          read(lucoef,trim(fmt_tmp)) MOcoeff_alpha(START:START+DIFF-1)
          START = START + 5
       end do

       if (.not. restricted) then
          ! read beta orbital energies
          read(lucoef,*) ! skip line which gives the number of nmo
          START = 1
          FIN   = nmo
          DIFF  = 5
          BUFF  = ceiling(real(FIN)/real(DIFF))
          allocate(MOenergy_beta(nmo))
          do i=1,BUFF
             if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
             write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
             read(lucoef,trim(fmt_tmp)) MOenergy_beta(START:START+DIFF-1)
             START = START + 5
          end do
          ! read beta orbital coefficients
          read(lucoef,'(A49,I12)') ACHAR, ncoeff
          START = 1
          FIN   = ncoeff
          DIFF  = 5
          BUFF  = ceiling(real(FIN)/real(DIFF))
          allocate(MOcoeff_beta(ncoeff))
          do i=1,BUFF
             if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
             write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
             read(lucoef,trim(fmt_tmp)) MOcoeff_beta(START:START+DIFF-1)
             START = START + 5
          end do
       end if

       ! everything important has been read, close the file.
       close(lucoef)

       if (present(ao_basis)) then
         ! put the result in the ao_basis type format
         ao_basis%nshells       = nshells
         ao_basis%nao           = nao
         ao_basis%basis_angular = 1
         nullify(ao_basis%gtos)
         allocate(ao_basis%gtos(nshells))
         do k=1,nshells
            nullify  (ao_basis%gtos(k)%exponent)
            allocate (ao_basis%gtos(k)%exponent(1))
            nullify  (ao_basis%gtos(k)%coefficient)
            allocate (ao_basis%gtos(k)%coefficient(1))
            coeff_norm = normalization_constant_psi4(orbital_momentum(k),all_exponents(k))
            ao_basis%gtos(k)%exponent     = all_exponents(k)
            ao_basis%gtos(k)%coefficient  = coeff_norm
            ao_basis%gtos(k)%coord        = coord_shells(3*(k-1)+1:3*(k-1)+3)
            ao_basis%gtos(k)%orb_momentum = orbital_momentum(k)+1 ! Following the DIRAC/DALTON convention of numbering momentum from 1 to l+1
            ao_basis%gtos(k)%atom_number  = shell_to_atom_map(k)
            ao_basis%gtos(k)%atom_element = atomic_numbers(ao_basis%gtos(k)%atom_number)
          enddo
       end if

       if (.not.present(spinors)) return

       if (spinors%type == 1 .or. spinors%type == 4) then
           restricted = .true.
       else
           restricted = .false.
       end if

       if (restricted) then
          call alloc_mo_basis (spinors,nao,nmo,1)
          spinors%total_energy = total_energy
          spinors%irrep        = 0 ! because no symmetry
          spinors%coeff_r      = 0.D0
          ij = 0
          do i=1,nmo
           spinors%energy(i) = MOenergy_alpha(i)
           do j=1,nao
             ij = ij+1
             spinors%coeff_r(j,i,1) = MOcoeff_alpha(ij)
           end do
          end do
       else
          call alloc_mo_basis (spinors,nao,2*nmo,2)
          spinors%total_energy = total_energy
          spinors%irrep        = 0 ! because no symmetry
          spinors%coeff_r      = 0.D0
          ij = 0
          do i=1,nmo
           spinors%energy(2*i-1) = MOenergy_alpha(i)
           spinors%energy(2*i) = MOenergy_beta(i)
           do j=1,nao
             ij = ij+1
             spinors%coeff_r(j,2*i-1,1) = MOcoeff_alpha(ij)
             spinors%coeff_r(j,2*i,2)   = MOcoeff_beta(ij)
           end do  
          end do
          deallocate(MOcoeff_beta)
          deallocate(MOenergy_beta)
       end if

       ! Reorder the coefficients back to normal (i.e. back to Psi4/Dirac/Interest standard.)
       call reorder_spinor_fchk(spinors,ao_basis,.true.)

       ! Deallocation
       deallocate(MOcoeff_alpha)
       deallocate(MOenergy_alpha)
       deallocate(all_exponents)
       deallocate(orbital_momentum)
       deallocate(coord_shells)
       deallocate(shell_to_atom_map)

      end subroutine read_from_fchk

      subroutine write_mos_to_fchk(spinors,ao_basis,molecule,fname)

       use rose_mo
       use rose_ao
       use rose_utils
       use read_xyz

       implicit none

       character(len=*),intent(in)           :: fname
       type(mo_basis), intent(inout)         :: spinors
       type(molecule_info), intent(in)       :: molecule
       type(basis_set_info_t), intent(inout) :: ao_basis
       character(len=80)      :: ACHAR
       logical                :: tobe
       integer                :: oldcoef,newcoef
       integer                :: i, j, k, ij
       integer                :: START, FIN, BUFF, DIFF
       integer                :: nao, nmo, natoms, nshells, pureamd, pureamf
       real(8), allocatable   :: array_of_one(:), atomic_number(:)
       real(8), allocatable   :: MOcoeff_alpha(:),MOcoeff_beta(:)
       real(8), allocatable   :: MOenergy_alpha(:),MOenergy_beta(:)
       real(8), allocatable   :: atom_coord(:), allcoord(:), allexponent(:)
       real(8)                :: pf1,pf2,pf3,pf4
       character(len=12)      :: fmt_tmp
       character(len=10)      :: format_size_single
       character(len=8)       :: single_precision
       character(len=20)      :: file_name

       if (spinors%type .gt. 2) then
          stop 'FCHK extension is only for real algebra'
       end if

       format_size_single = '(A1,I1,A7)'
       single_precision   = 'ES16.8)'

       file_name = trim(fname)//".fchk"
       inquire (file=file_name,exist=tobe)

       ! Reorder the coefficients to FCHK order
       call reorder_spinor_fchk(spinors,ao_basis,.false.)

       nao = spinors%nao
       nmo = spinors%nmo
       natoms = molecule%natoms
       nshells = ao_basis%nshells

       ! print warning in case we overwrite a file and open the file
       call get_free_fileunit(newcoef)
       inquire (file=file_name,exist=tobe)
       if (tobe) then
          print*, "warning: overwriting file "//file_name
          open (newcoef,file=file_name,status='OLD',FORM='FORMATTED',access='SEQUENTIAL')
       else
          open (newcoef,file=file_name,status='NEW',FORM='FORMATTED',access='SEQUENTIAL')
       end if

       write(newcoef,'(A27)')"Generated by GENIBO program"
       ! skip writing the second line that gives the method and the basis set name

       ! write Number of atoms
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of atoms                            ","I","  ",natoms
       ! write molecular informations:
       write(newcoef,'(A43,A1,A5,I12)') &
        "Charge                                     ","I","  ",molecule%charge
       write(newcoef,'(A43,A1,A5,I12)') &
        "Multiplicity                               ","I","  ",molecule%multiplicity
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of electrons                        ","I","  ",molecule%nelec
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of alpha electrons                  ","I","  ",molecule%nalpha
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of beta electrons                   ","I","  ",molecule%nbeta

       ! write number of basis functions
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of basis functions                  ","I","  ",nao

       ! write atomic number
       allocate(atomic_number(natoms))
       do i=1,natoms
        do k=1,nshells
          if (ao_basis%gtos(k)%atom_number .eq. i) then
             atomic_number(i) = ao_basis%gtos(k)%atom_element
          end if
        enddo
       enddo
       write(newcoef,'(A43,A1,A5,I12)') &
        "Atomic numbers                             ","I","N=",natoms
       START = 1
       FIN   = natoms
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          write(newcoef,fmt_tmp) nint(atomic_number(START:START+DIFF-1))
          START = START + 6
       end do

       ! same for nuclear charges (real of the atomic number...)
       write(newcoef,'(A43,A1,A5,I12)') &
        "Nuclear charges                            ","R","N=",natoms
       START = 1
       FIN   = natoms
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) atomic_number(START:START+DIFF-1)
          START = START + 5
       end do
       deallocate(atomic_number)

       ! LV: to be simplified (coordinates are now available in molecule variable)
       ! write current cartesian coordinates:
       allocate(atom_coord(3*natoms))
       do i=1,natoms
        do k=1,nshells
          if (ao_basis%gtos(k)%atom_number .eq. i) then
             atom_coord(3*(i-1)+1) = ao_basis%gtos(k)%coord(1)
             atom_coord(3*(i-1)+2) = ao_basis%gtos(k)%coord(2)
             atom_coord(3*(i-1)+3) = ao_basis%gtos(k)%coord(3)
          end if
        enddo
       enddo
       write(newcoef,'(A43,A1,A5,I12)') &
        "Current cartesian coordinates              ","R","N=",3*natoms
       START = 1
       FIN   = 3*natoms
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) atom_coord(START:START+DIFF-1)
          START = START + 5
       end do

       ! skip writing Integer atomic weights and Real atomic weights
       ! write number of shells:
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of primitive shells                 ","I","  ",nshells

       pureamd = 1
       pureamf = 1
       ! write that cartesian functions are used
       write(newcoef,'(A43,A1,A5,I12)') &
        "Pure/Cartesian d shells                    ","I","  ",pureamd
       write(newcoef,'(A43,A1,A5,I12)') &
        "Pure/Cartesian f shells                    ","I","  ",pureamf

       ! write shell types:
       write(newcoef,'(A43,A1,A5,I12)') &
        "Shell types                                ","I","N=",nshells
       START = 1
       FIN   = nshells
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          write(newcoef,fmt_tmp) (ao_basis%gtos(START:START+DIFF-1)%orb_momentum - 1)
          START = START + 6
       end do

       ! write number of primitives per shell:
       write(newcoef,'(A43,A1,A5,I12)') &
        "Number of primitives per shell             ","I","N=",nshells
       START = 1
       FIN   = nshells
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          allocate(array_of_one(DIFF))
          array_of_one = 1.D0
          write(newcoef,fmt_tmp) nint(array_of_one)
          deallocate(array_of_one)
          START = START + 6
       end do

       ! write shell to atom map
       write(newcoef,'(A43,A1,A5,I12)') &
        "Shell to atom map                          ","I","N=",nshells
       START = 1
       FIN   = nshells
       DIFF  = 6
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 6) DIFF = FIN - START + 1
          write(fmt_tmp,'(A1,I1,A4)') '(',DIFF,'I12)'
          write(newcoef,fmt_tmp) ao_basis%gtos(START:START+DIFF-1)%atom_number
          START = START + 6
       end do

       ! write primitive exponents
       allocate(allexponent(nshells))
       do i=1,nshells
          allexponent(i) = ao_basis%gtos(i)%exponent(1)
       enddo
       write(newcoef,'(A43,A1,A5,I12)') &
        "Primitive exponents                        ","R","N=",nshells
       START = 1
       FIN   = nshells
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) allexponent(START:START+DIFF-1)
          START = START + 5
       end do

       ! write contration coefficients
       write(newcoef,'(A43,A1,A5,I12)') &
        "Contraction coefficients                   ","R","N=",nshells
       START = 1
       FIN   = nshells
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          allocate(array_of_one(DIFF))
          array_of_one = 1.D0
          write(newcoef,fmt_tmp) array_of_one
          deallocate(array_of_one)
          START = START + 5
       end do

       allocate(allcoord(3*nshells))
       do i=1,nshells
         allcoord(3*(i-1)+1) = ao_basis%gtos(i)%coord(1)
         allcoord(3*(i-1)+2) = ao_basis%gtos(i)%coord(2)
         allcoord(3*(i-1)+3) = ao_basis%gtos(i)%coord(3)
       enddo
       ! write coordinates of each shell:
       write(newcoef,'(A43,A1,A5,I12)') &
        "Coordinates of each shell                  ","R","N=",3*nshells
       START = 1
       FIN   = 3*nshells
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) allcoord(START:START+DIFF-1)
          START = START + 5
       end do

       ! write the total energy:
       write(newcoef,'(A43,A1,A5,es22.15)') &
        "Total Energy                               ","R","",spinors%total_energy
       
       ! write MO energies and coefficients:
       ij = 0
       if (spinors%type .eq. 1) then
         allocate(MOenergy_alpha(nmo))
         allocate(MOcoeff_alpha(nmo*nao))
         do i=1,nmo
           MOenergy_alpha(i) = spinors%energy(i)
           do j=1,nao
             ij = ij+1
             MOcoeff_alpha(ij) = spinors%coeff_r(j,i,1)
           end do
         end do
       else if (spinors%type .eq. 2) then
         allocate(MOenergy_alpha(nmo/2))
         allocate(MOenergy_beta(nmo/2))
         allocate(MOcoeff_alpha(nmo*nao/2))
         allocate(MOcoeff_beta(nmo*nao/2))
         do i=1,nmo/2
           MOenergy_alpha(i) = spinors%energy(2*i-1)
           MOenergy_beta(i)  = spinors%energy(2*i)
           do j=1,nao
             ij = ij+1
             MOcoeff_alpha(ij) = spinors%coeff_r(j,2*i-1,1)
             MOcoeff_beta(ij)  = spinors%coeff_r(j,2*i,2)
           end do 
         end do
       end if


       ! write alpha MO energies:
       write(newcoef,'(A43,A1,A5,I12)') & 
        "Alpha Orbital Energies                     ","R","N=",size(MOenergy_alpha)
       START = 1
       FIN   = size(MOenergy_alpha)
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) MOenergy_alpha(START:START+DIFF-1)
          START = START + 5
       end do
       deallocate(MOenergy_alpha)      

       ! write alpha MO coefficients:
       do i=1,size(MOcoeff_alpha)
        if (abs(MOcoeff_alpha(i)) .lt. 1D-12) MOcoeff_alpha(i) = 0.D0
       enddo
       write(newcoef,'(A43,A1,A5,I12)') &
        "Alpha MO coefficients                      ","R","N=",size(MOcoeff_alpha)
       START = 1
       FIN   = size(MOcoeff_alpha)
       DIFF  = 5
       BUFF  = ceiling(real(FIN)/real(DIFF))
       do i=1,BUFF
          if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
          write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
          write(newcoef,fmt_tmp) MOcoeff_alpha(START:START+DIFF-1)
          START = START + 5
       end do
       deallocate(MOcoeff_alpha)

       if (spinors%type .eq. 2) then
          ! write beta MO energies:
          do i=1,size(MOcoeff_beta)
           if (abs(MOcoeff_beta(i)) .lt. 1D-12) MOcoeff_beta(i) = 0.D0
          enddo
          write(newcoef,'(A43,A1,A5,I12)') &
        "Beta Orbital Energies                      ","R","N=",size(MOenergy_beta)
          START = 1
          FIN   = size(MOenergy_beta)
          DIFF  = 5
          BUFF  = ceiling(real(FIN)/real(DIFF))
          do i=1,BUFF
             if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
             write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
             write(newcoef,fmt_tmp) MOenergy_beta(START:START+DIFF-1)
             START = START + 5
          end do
          deallocate(MOenergy_beta)
          ! write beta MO coefficients:
          write(newcoef,'(A43,A1,A5,I12)') &
        "Beta MO coefficients                       ","R","N=",size(MOcoeff_beta)
          START = 1
          FIN   = size(MOcoeff_beta)
          DIFF  = 5
          BUFF  = ceiling(real(FIN)/real(DIFF))
          do i=1,BUFF
             if ((FIN - START) .lt. 5) DIFF = FIN - START + 1
             write(fmt_tmp,format_size_single) '(',DIFF,trim(single_precision)
             write(newcoef,fmt_tmp) MOcoeff_beta(START:START+DIFF-1)
             START = START + 5
          end do
          deallocate(MOcoeff_beta)
       end if

       deallocate(allcoord,atom_coord,allexponent)

       close(newcoef)

       ! Reorder the coefficients back to normal (i.e. back to Psi4/Dirac/Interest standard.)
       call reorder_spinor_fchk(spinors,ao_basis,.true.)

      end subroutine write_mos_to_fchk

end module fchk_interface
