      MODULE READ_XYZFILE

!     Module to read xyz files with the format
!       Line 1: Number of atoms
!       Line 2: Comment
!       Other lines: atom name, coordinates 
        IMPLICIT NONE

        PRIVATE

        PUBLIC atom, firstatom
        PUBLIC get_atoms, add_atom, find_atom, count_atoms, find_charge
        PUBLIC print_atoms, del_atoms, read_xyz_atoms, read_atoms
        PUBLIC molecule_info, unit_conversion, unit_conversion_set, xtang
        PUBLIC open_file, filter_comments

        ! We define some global data
        INTEGER, PARAMETER :: namelen=12, linelen=192
        REAL(KIND=8), PARAMETER :: D0 = 0.0D0
        real(8), parameter :: xtang = 0.5291772083d0
        ! The unit conversion global variable is used when reading/writing xyz files
        REAL(KIND=8) :: UNIT_CONVERSION = 1.D0 / XTANG
        LOGICAL      :: UNIT_CONVERSION_SET=.FALSE.
        INTEGER            :: i,ios
        CHARACTER(LEN=4), PARAMETER :: NUCLABEL(1:118)=(/               &
     &  'H  ', 'He ', 'Li ', 'Be ', 'B  ',                              &
     &  'C  ', 'N  ', 'O  ', 'F  ', 'Ne ',                              &
     &  'Na ', 'Mg ', 'Al ', 'Si ', 'P  ',                              &
     &  'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ',                              &
     &  'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ',                              &
     &  'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ',                              &
     &  'Ga ', 'Ge ', 'As ', 'Se ', 'Br ',                              &
     &  'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ',                              &
     &  'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ',                              &
     &  'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ',                              &
     &  'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ',                              &
     &  'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ',                              &
     &  'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ',                              &
     &  'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ',                              &
     &  'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ',                              &
     &  'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ',                              &
     &  'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ',                              &
     &  'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ',                              &
     &  'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ',                              &
     &  'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ',                              &
     &  'Md ', 'No ', 'Lr ', 'Rf ', 'Db ',                              &
     &  'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ',                              &
     &  'Rg ', 'Cn ', 'Nh ', 'Fl ', 'Mc ',                              &
     &  'Lv ', 'Ts ', 'Og '/)
        CHARACTER(LEN=8), DIMENSION(2), PARAMETER :: filetypes=(/'MOL', &
     &                                                           'XYZ'/)

        INTERFACE open_file
          MODULE PROCEDURE open_file_name, open_file_num, open_scratch
        END INTERFACE

        TYPE :: atom
!     Defines a list for storing atomic data, like coordinates, charges and nuclear exponents.
          character(len=namelen)  :: name=' '
          real(KIND=8)            :: x(3)=(/(0.0,i=1, 3)/)
          real(KIND=8)            :: charge=0.0D0    ! this is the actual charge that may be taken zero or fractional
          real(KIND=8)            :: bascharge=0.0D0 ! this is the charge for the specified element that is used when searching the basis set library
          real(KIND=8)            :: nucexp=0.0D0
          TYPE(atom), POINTER     :: next
        END TYPE atom

        TYPE(atom), POINTER :: firstatom

    type :: molecule_info ! This type is used for fchk/psi4/pyscf interfaces
        character(len=120) :: name
        integer :: natoms
        integer :: nelec, nalpha, nbeta, ncore, nvale
        integer :: charge=0
        integer :: multiplicity=1 ! default is singlet.
        real(kind=8) :: znumber=0.0
        type(atom),pointer :: firstatom
    end type molecule_info

       INTERFACE get_atoms
           module procedure read_xyz_atoms
           module procedure get_atoms
       END INTERFACE

      CONTAINS

        INTEGER FUNCTION open_file_name(filename)
!     Opens a file and returns its identifying number (>0) if the operation was successful.
!     If the opening was unsuccessful a negative or 0 is returned.
!      0  -   File does not exist.
!     -1  -   Too many files opened.
!     -99 -   Error while executing the open statement.
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER, PARAMETER :: maxfile=99
          LOGICAL :: ex
          INTEGER :: n

          INQUIRE(file=filename, exist=ex)
          IF (ex) THEN
            DO n=1,maxfile
              INQUIRE(n, opened=ex)
              IF (.NOT.(ex)) THEN
                OPEN(n, file=filename, err=563)
                EXIT
              END IF
            END DO
 563        open_file_name=-99
            IF (n>maxfile) THEN
              open_file_name=-1
            ELSE
              open_file_name=n
            END IF
          ELSE
            open_file_name=0
          END IF

        END FUNCTION open_file_name

        INTEGER FUNCTION open_scratch()
!     Opens a scratch file and returns its identifying number (>0).
          INTEGER, PARAMETER :: maxfile=99
          LOGICAL :: ex
          INTEGER :: n

          DO n=1,maxfile
            INQUIRE(n, opened=ex)
            IF (.NOT.(ex)) THEN
              OPEN(n, status='SCRATCH')
              EXIT
            END IF
          END DO
 563      open_scratch=-99
          IF (n>maxfile) THEN
            open_scratch=-1
          ELSE
            open_scratch=n
          END IF

        END FUNCTION open_scratch

        INTEGER FUNCTION open_file_num(filename, filenum, formatted,    &
     &                                 lupri)
!     Opens a file on a requested unit number.
          CHARACTER(LEN=*), INTENT(IN)  :: filename
          INTEGER, INTENT(IN)           :: filenum
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: formatted
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          LOGICAL                       :: ex

          INQUIRE(file=filename, exist=ex)
          IF (.NOT. ex) THEN
            open_file_num=0
            RETURN
          END IF
          INQUIRE(filenum, exist=ex)
          IF (.NOT. ex) THEN
            open_file_num=-2
            RETURN
          END IF
          OPEN(filenum, file=filename, FORM=formatted, err=463)
          open_file_num=filenum

          RETURN
 463      open_file_num=-99

        END FUNCTION open_file_num

        FUNCTION get_extension(filename)
          CHARACTER(LEN=*), INTENT(IN) :: filename
          CHARACTER(LEN=LEN(filetypes(1))) :: get_extension
          INTEGER :: i, size_ft
          LOGICAL :: ex

          size_ft = SIZE(filetypes)
          DO i=1, size_ft
            INQUIRE(file=filename//'.'//TRIM(filetypes(i)), exist=ex)
            IF (ex) THEN
              get_extension=TRIM(filetypes(i))
              EXIT
            END IF
          END DO
          IF (i>size_ft) get_extension=' '
        END FUNCTION get_extension

        SUBROUTINE filter_comments(filein,fileout)
          INTEGER             :: filein,fileout
          CHARACTER(LEN=192)  :: str
          INTEGER             :: i,ios

          REWIND(filein)
          fileout=open_scratch()

          DO
            READ(filein,'(192a)',iostat=ios) str
            IF (ios.ne.0) EXIT
            DO i=1,LEN(str)
              SELECT CASE (str(i:i))
              CASE ('#','!')
                EXIT
              END SELECT
            END DO
            WRITE(fileout,*) str(1:i-1)
          END DO
          REWIND(filein)
          REWIND(fileout)

        END SUBROUTINE filter_comments

        SUBROUTINE add_atom(latom, fatom, name, x)
!
!       fatom: pointer to the first atom in the list
!       latom: pointer to the last atom in the list
!

          TYPE(atom), POINTER :: latom,fatom,newatom,sameatom
          CHARACTER(LEN=*), INTENT(IN)     :: name
          REAL(KIND=8),INTENT(IN),OPTIONAL :: x(1:3)

!         Make a new atom
          ALLOCATE(newatom)
          newatom%name=name
          IF (PRESENT(x)) newatom%x=x
          newatom%charge=find_charge(name)
          newatom%bascharge=find_charge(name)
          newatom%next=>NULL()
!
!         This might be the first call, if so: set fatom to this atom 

          IF (.not. associated(latom)) THEN
             fatom=>newatom
             latom=>newatom
             RETURN
          END IF

!         Find out whether we already have atoms of this type 
          sameatom=>find_atom(fatom, name, .true.)
          IF (associated(sameatom)) THEN
!            insert after the others of this type
             newatom%next=>sameatom%next
             sameatom%next=>newatom
!            important: update the latom pointer as well
             if (latom%name .eq. sameatom%name) latom=>newatom
          ELSE
!            we did not see this type before, put at the end of the list
             latom%next=>newatom
             latom=>newatom
          END IF


        END SUBROUTINE add_atom

        FUNCTION find_atom(fatom, name, last)
!     JHS050707 Returns a pointer to the first (or last)  atom in the list with this name 
!               or NULL() if is no such atom. 
           TYPE(atom), POINTER :: fatom 
           TYPE(atom), POINTER :: find_atom
           TYPE(atom), POINTER :: mol
           LOGICAL, OPTIONAL   :: last
           CHARACTER(LEN=*), INTENT(IN)    :: name

           find_atom=>NULL()

           mol=>fatom
           DO WHILE(ASSOCIATED(mol))
              IF (mol%name==name) THEN
                 find_atom=>mol
                 if (.not. present(last)) then
                    return
                 else
                    if (.not. last) return
                 end if
              END IF
              mol=>mol%next
           END DO
        END FUNCTION find_atom

        INTEGER FUNCTION count_atoms(fatom, name)
!     JHS050707 Counts the length of the atomlist.
         ! TYPE(atom), POINTER, INTENT(IN) :: fatom
          TYPE(atom), POINTER :: fatom
          TYPE(atom), POINTER             :: mol
          CHARACTER(LEN=*), OPTIONAL      :: name
          count_atoms=0
          mol=>fatom
          IF (PRESENT(name)) THEN
            DO WHILE (associated(mol))
              IF (name==mol%name) count_atoms=count_atoms+1
              mol=>mol%next
            END DO
          ELSE
            DO WHILE (associated(mol))
              count_atoms=count_atoms+1
              mol=>mol%next
            END DO
          END IF
        END FUNCTION count_atoms

        SUBROUTINE PRINT_ATOMS(firstatom, lupri)
          TYPE(atom), POINTER :: firstatom
          INTEGER, INTENT(IN) :: lupri
          TYPE(atom), POINTER :: mol

          WRITE(lupri,*)' XYZ coordinates read'
          WRITE(lupri,'(a,t13,a,t25,a,t37,a,t49,a,t61,a)') 'Label',     &
     &         'x','y','z','charge','nuc exp'
          mol=>firstatom
          DO
            IF (associated(mol)) THEN
              WRITE(lupri,'(a,t9, 5g13.6)')                             & 
     &          TRIM(mol%name), mol%x, mol%charge, mol%nucexp
              mol=>mol%next
            ELSE
              EXIT
            ENDIF
          END DO
        WRITE(lupri,'(a)') ' '

        END SUBROUTINE PRINT_ATOMS

        SUBROUTINE DEL_ATOMS(firstatom)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          TYPE(atom),POINTER :: firstatom
          TYPE(atom), POINTER :: mol

          IF (.NOT.associated(firstatom)) RETURN

          mol=>firstatom
          DO
            IF (associated(mol%next)) THEN
              firstatom=>mol%next              
              DEALLOCATE(mol)
              mol=>firstatom
            ELSE
              NULLIFY(firstatom)
              EXIT
            ENDIF
          END DO

        END SUBROUTINE DEL_ATOMS

        SUBROUTINE READ_ATOMS(filenum, linelen, fatom, lupri)
!       Reads atom name and coordinates from the .XYZ file
!       and writes it to the atoms list.
!       TODO: better checking for errors in the xyz file

          TYPE(atom), POINTER :: fatom
          INTEGER, INTENT(IN)                :: filenum, linelen, lupri
          TYPE(atom), POINTER                :: latom
          CHARACTER(LEN=namelen)             :: name
          REAL(KIND=8)                       :: x(3)
          INTEGER                            :: ios

          latom=>null()

          ios=0
          do while (ios.eq.0)

            name='init'
            x(1:3)=0.0

            READ(filenum, *,iostat=ios) name, x(1), x(2), x(3)

            if (ios.eq.0.AND.name.ne.'init') then
              CALL ADD_ATOM(latom,fatom,name,x)
            end if
          end do

        END SUBROUTINE READ_ATOMS

        subroutine get_atoms (natoms,coordinates,fatom,atomic_numbers,symbols)

!       Process list of coordinates and symbols and returns them as a linked list

        integer, intent(in)                 :: natoms
        real(8), intent(in)                 :: coordinates(3,*)
        type(atom), pointer                 :: fatom
        integer, intent(in) , optional      :: atomic_numbers(*)
        character*(*), intent(in), optional :: symbols(*)
        character(len=:), allocatable       :: symbol

        type(atom), pointer       :: latom
        integer                   :: iatom

        allocate(fatom)
        latom=>null()
        do iatom =1, natoms
           if (present(symbols)) then
              symbol = symbols(iatom)
           elseif (present(atomic_numbers)) then
              symbol = NUCLABEL(atomic_numbers(iatom))
           else
              print*,"symbols or atomic_numbers should be specified in get_atoms"
              stop "programming error"
           end if
           call add_atom(latom,fatom,trim(symbol),coordinates(:,iatom))
        end do

        end subroutine get_atoms

        SUBROUTINE READ_XYZ_ATOMS (molfile,fatom, natoms, znumber)

!       Reads the coordinates from an XYZ-FILE and returns them as a linked list

        INTEGER, INTENT(IN) :: molfile
        type(atom), pointer, intent(out) :: fatom
        integer, intent(out), optional :: natoms
        real(8), intent(out), optional :: znumber

        real(8)             :: xyz(3)
        type(atom), pointer :: latom
        INTEGER             :: lupri=6
        INTEGER             :: ios, iatom, natoms_read
        character(len=20)   :: symbol

!       In the XYZ format coordinates are supposed to start at the third line
        REWIND (molfile)
!       Check whether the first line has a number that can be used to check the number of atoms
        natoms_read = 0
        READ (molfile,*,iostat=ios) natoms_read
!       Skip the comment line
        read (molfile,*)
        allocate(fatom)
        latom=>null()
        if (present(znumber)) znumber = 0.0
        do iatom =1, natoms_read
          read(molfile, *,iostat=ios) symbol, xyz
          xyz = xyz/xtang  ! xyz files have coordinates in Angstroms, convert to au
          call add_atom(latom,fatom,trim(symbol),xyz)
          if (present(znumber)) znumber = znumber + find_charge(trim(symbol),6)
        end do
        if (present(natoms)) natoms = natoms_read

       END SUBROUTINE READ_XYZ_ATOMS

        INTEGER FUNCTION FIND_CHARGE(atomname,LUPRI)
          CHARACTER(LEN=*),INTENT(IN) :: atomname
          INTEGER,OPTIONAL,INTENT(IN) :: LUPRI
          INTEGER :: i


          find_charge=0
          DO i=1,SIZE(NUCLABEL)
            IF (atomname==NUCLABEL(i)) THEN
              find_charge=i
              RETURN
            END IF
          END DO
          IF (PRESENT(LUPRI)) WRITE(LUPRI,*)                            &
     &         'WARNING: No charge found for element: ', atomname

        END FUNCTION

      END MODULE READ_XYZFILE
