      MODULE READ_XYZ
        IMPLICIT NONE

        PRIVATE

        PUBLIC get_atoms, del_atoms
!     Module to read xyz files with the format
!       Line 1: Number of atoms
!       Line 2: Comment
!       Other lines: atom name, coordinates 

        INTEGER, PARAMETER :: namelen=8, linelen=192
        REAL(KIND=8), PARAMETER :: D0 = 0.0D0
        real(8), parameter, public :: xtang = 0.5291772083d0
        INTEGER            :: i,ios

        TYPE,public :: atom
!     Defines a list for storing atomic data, like coordinates, charges and nuclear exponents.
          CHARACTER(LEN=namelen)  :: name=' '
          REAL(KIND=8)            :: x(3)=(/(0.0,i=1, 3)/)
          REAL(KIND=8)            :: charge=0.0 ! this is the actual charge that may be taken zero or fractional
          REAL(KIND=8)            :: bascharge=0.0 ! this is the charge for the specified element that is used when searching the basis set library
          REAL(KIND=8)            :: nucexp=0.0
          TYPE(atom), POINTER     :: next
        END TYPE atom

!     Global variables
!     Start of linked lists of atoms, basistypes and atomtypes
        TYPE(atom), POINTER :: firstatom

!     Symmetry of the molecule
        LOGICAL :: ADDSYM=.TRUE.
        INTEGER :: NSYMOP=0
        INTEGER :: natoms_given

!     The unit conversion global variable is used when reading/writing xyz files
       REAL(KIND=8),PUBLIC  :: UNIT_CONVERSION = 1.D0 / XTANG

       INTERFACE GET_ATOMS
           module procedure READ_XYZ_ATOMS
       END INTERFACE

      CONTAINS

        SUBROUTINE ADD_ATOM(latom, fatom, name,x)
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
!         and allocate also the last atom pointer

          IF (.not. associated(latom)) THEN
             fatom=>newatom
             ALLOCATE(latom)
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


        END SUBROUTINE ADD_ATOM

        SUBROUTINE READ_ATOMS(filenum, linelen,fatom,lupri)
!       Reads atom name and coordinates from the .XYZ file 
!       and writes it to the atoms list.
!       TODO: better checking for errors in the xyz file
          TYPE(atom), POINTER :: fatom
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
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

        SUBROUTINE POS_FILE(infile,lupri,linelen,section,subsection,name)
!         Position a file to a specific keyword.
        
          INTEGER, INTENT(IN)   :: infile,lupri,linelen
          CHARACTER(LEN=*), INTENT(IN) ::section,subsection,name
          CHARACTER(LEN=linelen) :: str,pressect,pressub
          character(len=linelen) :: tmp_name
          INTEGER                :: ios

          REWIND(infile,iostat=ios)
          DO            
            read(infile,*,iostat=ios) str
            IF (ios.ne.0) THEN
             WRITE(lupri,*) '*** ERROR KEY NOT FOUND. ***'
             write(lupri,*) 'POS_FILE: error in read command ! ios=',ios
             write(lupri,*) 'infile(in)=',infile
             WRITE(lupri,*) "section=",section
             WRITE(lupri,*) "subsection=",subsection
             WRITE(lupri,*) "name=",name
             WRITE(lupri,*) "linelen=",linelen
             STOP '*** ERROR KEY NOT FOUND. ***'
            ELSE IF (str(1:2)=='**') THEN
              pressect=str              
            ELSE IF (str(1:1)=='*') THEN
              pressub=str
            ELSE IF (pressect==section.AND.pressub==subsection) THEN
              IF (str==name) EXIT
            END IF
          END DO
        END SUBROUTINE POS_FILE

        SUBROUTINE PRINT_ATOMS(firstatom, lupri)
          !TYPE(atom), POINTER, INTENT(IN) :: firstatom
          TYPE(atom), POINTER :: firstatom
          INTEGER, INTENT(IN) :: lupri
          TYPE(atom), POINTER :: mol

          WRITE(lupri,*)' XYZ coordinates read'
          WRITE(lupri,'(a,t13,a,t25,a,t37,a,t49,a,t61,a)') 'Label',     &
     &         'x','y','z','charge','nuc exp'
          mol=>firstatom
          DO
            IF (associated(mol)) THEN
              WRITE(lupri,'(a,t9, 5g12.6)')                             & 
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

        SUBROUTINE READ_XYZ_ATOMS (molfile,atoms,luprint)

!       Reads the coordinates from an XYZ-FILE and returns them as a linked list

        INTEGER, INTENT(IN) :: molfile
        INTEGER, INTENT(IN), optional  :: luprint
        type(atom), pointer, intent(out) :: atoms

        INTEGER             :: lupri
        INTEGER             :: ios, i

        if (present(luprint)) then
           lupri = luprint
        else
           lupri = 6
        end if 

!       In the XYZ format coordinates are supposed to start at the third line
        REWIND (molfile)
!       Check whether the first line has a number that can be used to check the number of atoms
        natoms_given = 0
        READ (molfile,*,iostat=ios) natoms_given
!       Skip the comment line
        read (molfile,*)
        CALL READ_ATOMS(molfile,linelen,firstatom,lupri)
!       Check number of atoms
        IF (natoms_given .NE. count_atoms(firstatom)) THEN
          WRITE(LUPRI,*) " Read ",count_atoms(firstatom)," atoms", &
                         " but expected ",natoms_given," atoms"
          STOP ' *** ERROR (READ_XYZ_COORDINATES) wrong number of atoms'
        END IF

        allocate(atoms)
        atoms=>firstatom

       END SUBROUTINE READ_XYZ_ATOMS

        INTEGER FUNCTION FIND_CHARGE(atomname,LUPRI)
          CHARACTER(LEN=*),INTENT(IN) :: atomname
          INTEGER,OPTIONAL,INTENT(IN) :: LUPRI
          INTEGER :: i

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

      END MODULE READ_XYZ
