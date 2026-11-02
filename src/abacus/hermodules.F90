!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

!     !!! MODULE FOR DIRAC XYZ/BASIS SETUP !!!
!     Module for reading DIRAC basis set information connected to xyz files.
!
!     In the long run the readingsection for the .inp file move to dirac/dirrdn.F
!     JHS 05-07-07

      MODULE HERMODULES
        USE READ_XYZFILE
        IMPLICIT NONE

        PRIVATE

        PUBLIC READ_XYZ
        PUBLIC basistype, atomtype, firstatomtype
        PUBLIC count_types
        PUBLIC UNIT_CONVERSION, UNIT_CONVERSION_SET, XTANG
#ifdef HAS_PELIB
        PUBLIC pelib_sym_xyz
#endif

        INTEGER, PARAMETER :: namelen=8, basislen=32,linelen=192
        REAL(KIND=8), PARAMETER :: D0 = 0.0D0
        INTEGER            :: i,ios

        TYPE :: basistype
!     Defines a list for storing the information about the basis sets.
          CHARACTER(LEN=namelen)     :: name=' '
          CHARACTER(LEN=basislen)    :: type=' '
          CHARACTER(LEN=basislen)    :: comp=' '
          CHARACTER(LEN=basislen)    :: file=' '
!         The following can be assigned for explicitly given basissets or Family type.
          INTEGER, POINTER           :: JCO(:) => NULL()
          REAL, POINTER              :: FAM(:) => NULL()
          TYPE(basistype), POINTER   :: next
        END TYPE basistype

        TYPE :: atomtype
             TYPE(atom), POINTER      :: first
             TYPE(basistype), POINTER :: basis
             INTEGER                  :: natom
             TYPE(atomtype), POINTER  :: next
        END TYPE atomtype

!     Global variables
!     Start of linked lists of basistypes and atomtypes
        TYPE(basistype), POINTER :: firstbas
        TYPE(atomtype), POINTER  :: firstatomtype

!     Symmetry of the molecule
        LOGICAL :: ADDSYM=.TRUE.
        INTEGER :: NSYMOP=0
        CHARACTER(LEN=1) :: KASYM(3,3)=' '
        INTEGER :: natoms_given
#ifdef HAS_PELIB
!     Variable for PElib (to ensure symmetry detection is disabled for xyz)
        LOGICAL :: pelib_sym_xyz = .false. ! only xyz input (READ_BASIS) sets it
#endif

      CONTAINS

        SUBROUTINE ADD_BASIS(prevbas, fbas,basdata,basnext)
!     JHS050707 Adds a new node to the basislist
          TYPE(basistype),POINTER :: prevbas,fbas
          TYPE(basistype),TARGET  :: basdata
          TYPE(basistype),POINTER :: bas
          TYPE(basistype),OPTIONAL,POINTER :: basnext

          ALLOCATE(bas)
          bas=basdata
          IF (PRESENT(basnext)) THEN
            bas%next=>basnext
          ELSE
            bas%next=>NULL()
          END IF
          IF (associated(prevbas)) THEN
            prevbas%next=>bas
          ELSE
            fbas=>bas
          END IF
          prevbas=>bas

        END SUBROUTINE ADD_BASIS

        SUBROUTINE ADD_JCO(bas,JCO)
!     JHS050707 Assigns an array with JCO components.
          TYPE(basistype)    :: bas
          INTEGER            :: JCO(:)
          INTEGER,TARGET,ALLOCATABLE,SAVE :: JCO2(:)

          ALLOCATE(bas%jco(SIZE(JCO)))
          bas%jco=JCO
        END SUBROUTINE

        SUBROUTINE ADD_FAM(filenum,linelen,bas)
!     JHS050707 Assigns for different types of family basissets an array FAM.
!               And reads the needed coefficients from file.
      use priunit_mod, only: LUPRI
          INTEGER,INTENT(IN) :: filenum,linelen
          TYPE(basistype)    :: bas
          REAL,TARGET,ALLOCATABLE,SAVE :: FAM(:)
          INTEGER               :: i,ios

! JHS test, iqm needs to be included as well!
          SELECT CASE(bas%type)
          CASE('WELLTEMP')
            ALLOCATE(FAM(5))
            read(filenum,*,iostat=ios) FAM(1:5)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ",                       &
     &          "Specify 5 coefficients for WELLTEMP basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('EVENTEMP')
            ALLOCATE(FAM(3))
            read(filenum,*,iostat=ios) FAM(1:3)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ",                       &
     &          "Specify 3 coefficients for EVENTEMP basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('GEOM')
            ALLOCATE(FAM(5))
            read(filenum,*,iostat=ios) FAM(1:2), FAM(5)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ",                       &
     &          "Specify 3 coefficients for GEOM basis."
              call exit(1)
            end if
            FAM(3)=D0
            FAM(4)=D0
            bas%fam=FAM
          CASE('FAMILY')
            ALLOCATE(FAM(1))
            read(filenum,*,iostat=ios) FAM(1)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ",                       &
     &          "Specify 1 coefficients for FAMILY basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('DUALFAMILY')
            ALLOCATE(FAM(2))
            read(filenum,*,iostat=ios) FAM(1:2)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ",                       &
     &          "Specify 2 coefficients for DUAL basis."
              call exit(1)
            end if
            bas%fam=FAM
          END SELECT
        END SUBROUTINE

        FUNCTION find_basis(fbas, name)
!     JHS050707 Returns a pointer to the next basisset with name in the list with basissets
!               or NULL() if there is no such basis.
           TYPE(basistype), POINTER :: fbas
           TYPE(basistype), POINTER :: find_basis
           TYPE(basistype), POINTER :: bas
           CHARACTER(LEN=*), INTENT(IN)       :: name

           bas=>fbas
           DO WHILE(ASSOCIATED(bas))
             IF (bas%name==name) THEN
               find_basis=>bas
               RETURN
             END IF
             bas=>bas%next
           END DO

           find_basis=>NULL()

        END FUNCTION find_basis

        INTEGER FUNCTION count_types(ftype)
          TYPE(atomtype), POINTER :: ftype
          TYPE(atomtype), POINTER :: typepointer

          count_types=0
          typepointer=>ftype
          DO WHILE (associated(typepointer))
             count_types=count_types+1
             typepointer=>typepointer%next
          END DO
        END FUNCTION count_types


        SUBROUTINE MAKE_ATOMTYPES(lupri)
!         Matches basis sets and atoms: an atomtype
!         consists of a set of identical nuclei (same charge, same basis)

          INTEGER, INTENT(IN),OPTIONAL       :: lupri
          TYPE(atom), POINTER                :: fatom=>NULL()
          TYPE(basistype), POINTER           :: basis=>NULL()
          TYPE(basistype), POINTER           :: defaultbasis=>NULL()
          TYPE(atomtype), POINTER            :: ltype=>NULL(),newtype
          CHARACTER(LEN=namelen)             :: name
          INTEGER                            :: i

!       Find default basis
        defaultbasis=>find_basis(firstbas,'DEFAULT')
        if (.not.(associated(defaultbasis)))                            &
     &     CALL QUIT('*** ERROR: no default basis found')
!       Set atom pointer to the first atom
        fatom=>firstatom
!       Loop over atoms
        DO WHILE (associated(fatom))
!         Make a new type
          ALLOCATE(newtype)
          name=fatom%name
          basis=>find_basis(firstbas,name)
          IF (.not. associated(basis)) basis=>defaultbasis
          newtype%first=>fatom
          newtype%basis=>basis
          newtype%natom=count_atoms(fatom,name)
          newtype%next=>NULL()
!         Check whether this is the first type
          IF (.not. associated(ltype)) THEN
             firstatomtype=>newtype
          ELSE
             ltype%next=>newtype
          END IF
          ltype=>newtype
!         Skip to the end of the list of atoms of this type
          DO  i = 1, ltype%natom
              fatom=>fatom%next
          END DO
        END DO
        NULLIFY(ltype)

        END SUBROUTINE MAKE_ATOMTYPES

        SUBROUTINE READ_BASIS(filenum,linelen,lupri)
        implicit none
#include "codata.h"
!       use codata_mod, only: XTANG
!      use dcbgen_mod, only: KCHARG
!      use dgroup_mod, only: NONATM
#include "dgroup.h" 

! using OPTIMI
#include "dcbgen.h"

!     JHS050707 Reads basis names and possibly the coefficients typed in the DIRAC.INP file.
!               Reads also Charges, nuclear exponents and the units of the XYZ file.
!     Stores it in the the basis list.
!     !!! TMP: TESTED FOR:
!                         - 'LARGE BASIS'
!     (Implemented also for the rest)
!         Needed to read the basis set information
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
          TYPE(basistype),POINTER            :: bas
          TYPE(basistype)                    :: tmpbas
          integer                            :: tmp_file,ios

!     FOR explicitly typed basissets.
          INTEGER,TARGET,ALLOCATABLE         :: JCO(:)

!     FOR UNIT CONVERSION and explicitly defined charges and/or nuclear exponents
          REAL(KIND=8)      :: ONE=1.0d0
          CHARACTER(LEN=32) :: UNITNAME='Angstrom'
          TYPE(atom), POINTER                :: mol

          CHARACTER(LEN=linelen)             :: str,section(3)
          INTEGER                            :: i,j,index,iqm
          CHARACTER(LEN=namelen)             :: name, prevname
          REAL(KIND=8)                       :: charge, bascharge, nucexp
          LOGICAL                            :: CONT=.FALSE.
          LOGICAL                            :: first_call=.TRUE.
          LOGICAL                            :: inMOLECULE
          CHARACTER(LEN=9) :: symtxt
          SAVE first_call

          UNIT_CONVERSION = 1.D0 / XTANG
          UNIT_CONVERSION_SET = .TRUE.

!     take out comments from the input file and initialize
         call filter_comments(filenum,tmp_file)

          mol=>null()
          bas=>null()
          index=1
          section=''

          NSYMOP=0
          ADDSYM=.TRUE.
          NONATM=.FALSE. ! Forcing linear symmetry not using atomic supersymmetry
          KASYM(:,:)=' '
          KCHARG=0
#ifdef HAS_PELIB 
          pelib_sym_xyz = .true.
#endif 
          inMOLECULE = .false.
          DO
            read(tmp_file,*,iostat=ios) str
            !write(6,*) ' --> next line is ',trim(str)

            IF (ios.ne.0) THEN
              EXIT
            ELSE IF (str(1:2)=='**') THEN
              section(1)=str
              !> make sure we process only sections related to XYZ
              if(inMOLECULE)then
                if(str(1:6) /= '**MOLE') exit
              end if

            ELSE IF (str(1:1)=='*') THEN
              section(2)=str
!             Reinitialize section after we read an *END
              IF (str=='*END')then
                section=''
                inMOLECULE = .false.
              END IF
              !> make sure we  process only section related to XYZ
              if(inMOLECULE)then
                select case(str(1:4))
                    case ('*CHA','*SYM','*BAS','*CEN','*COO')
                    case default
                        exit
                end select
              endif
            ELSE IF (str(1:1)=='.') THEN
               section(3)=str
            ELSE
              if(inMOLECULE)then
                write(lupri,*) "ERROR in READ_BASIS: unknown keyword ", &
     &          trim(str)
                call quit ("Error reading XYZ input in READ_BASIS")
              end if
            ENDIF

!
!           We only need to do something if we are reading a **MOLECULE block.
!           The options and suboptions are flattened out in one if then else if construction.
            IF(section(1)=='**MOLECULE') THEN

              inMOLECULE = .true.

!             add molecular charge. it is perhaps redundant to have *CHARGE and .CHARGE but i'm not
!             sure we will not want to add other keywords in the menu later on - plus this keeps in
!             line with the other fields that always place keywords (=section(3)) under single-star
!             menus
              IF (section(2)=='*CHARGE'.AND.                            &
     &            section(3)=='.CHARGE') THEN
                      read(tmp_file,*,iostat=ios) KCHARG

!            Symmetry specification (Only NOSYM supported, otherwise automatic recognition)
              ELSE IF (section(2)=='*SYMMETRY'.AND.                     &
     &                 section(3)=='.NOSYM') THEN
                       ADDSYM=.FALSE.
#ifdef HAS_PELIB
                       pelib_sym_xyz = .false.
#endif

!            Use linear symmetry even when single atom
              ELSE IF (section(2)=='*SYMMETRY'.AND.                     &
     &                 section(3)=='.LINEAR') THEN
                       NONATM=.TRUE.
!             Basis specification
              ELSE IF (section(2)=='*BASIS'.AND.                        &
     &                 section(3)=='.DEFAULT') THEN
!               In the subsection .DEFAULT we give only give the name of a 'LARGE BASIS' type
!               like cc-pVDZ etc. This will be the basisset for all atoms that are not
!               specified in the .SPECIAL section.
                read(tmp_file,*,iostat=ios) tmpbas%file
                if (ios.eq.0 .and. tmpbas%file>' ') then
                  tmpbas%name='DEFAULT'
                  tmpbas%comp='LARGE'
                  tmpbas%type='BASIS'
                  tmpbas%next=>NULL()
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  if( tmpbas%file(1:3) == 'ecp') read(tmp_file,*,iostat=ios) str
                  CYCLE
                else if (ios.ne.0) then
                  write(lupri,*) "ERROR(readbasis): ",                  &
     &                           "Expecting default basisset name."

                  exit
                end if
              ELSE IF(section(2)=='*BASIS'.AND.                         &
     &                section(3)=='.SPECIAL') THEN
!     JHS050707 For some or all atoms we can give an other basisset than the .DEFAULT
                tmpbas%name=' '
                tmpbas%comp='LARGE'
                tmpbas%type=' '
                tmpbas%file=' '
                tmpbas%next=>NULL()

                read(tmp_file,*,iostat=ios) tmpbas%name,tmpbas%type

!               Check whether the user has specified less than two keywords
                if (ios.ne.0) then
                     write(lupri,*) "ERROR(readbasis) :",            &
     &                 "Specify the atom name plus an option forr ",    &
     &                 ".SPECIAL basis."
                     CALL QUIT ('Input Error')
                end if
                if(tmpbas%name>' ' .and. tmpbas%type>' ' ) then
                  SELECT CASE(tmpbas%type)
                  CASE ('NOBASIS','POINTCHARGE')
                     CALL ADD_BASIS(bas,firstbas,tmpbas)
                  CASE ('BASIS','MOLFBAS')
                    backspace(tmp_file)
                    read(tmp_file,*,iostat=ios) tmpbas%name,&
     &                                      tmpbas%type,tmpbas%file
                    if (.not.(ios.ne.0) .and. tmpbas%file>' ') then
                      CALL ADD_BASIS(bas,firstbas,tmpbas)
                    else if (ios.ne.0) then
                     write(tmp_file,*) "ERROR(readbasis) :",            &
     &                 "Specify a filename for ",tmpbas%type
                    end if
                    if( tmpbas%file(1:3) == 'ecp') read(tmp_file,*,iostat=ios) str
                    cycle
                  CASE('INTGRL','EXPLICIT')
!     FIRST read rest of the line and then read the explicit coefficients.
!     TMP We read the coefficients in in the old way in write common block
                    read(tmp_file,*,iostat=ios) iqm
                    if ( ios.ne.0 .or. iqm.eq.0) then
                      write(lupri,*) "ERROR(readbasis) :",           &
     &                 "No iqm for ",tmpbas%type," specified."
                      cycle
                    end if
                    ALLOCATE(JCO(iqm))
                    backspace(tmp_file)
                    read(tmp_file,*,iostat=ios) iqm,JCO(1:iqm)
                    if (ios.ne.0) then
                      write(lupri,*) "ERROR(readbasis) :",           &
     &                 "Specify at least ",iqm," numbers"
                      cycle
                    end if
                    CALL ADD_JCO(tmpbas,JCO)
                    CALL ADD_BASIS(bas,firstbas,tmpbas)
                    DEALLOCATE(JCO)
!                   If we set the section back to empty, the read will skip over the
!                   explicitly typed basis set (the exponents are read in a different routine)
                    section(3)=' '
                    !> skip until we reach a new keyword/section
                    do
                        read(tmp_file,*,iostat=ios) str
                        if(ios .ne. 0) exit
                        SELECT CASE(trim(str(1:1)))
                          CASE ('*','.')
                            exit
                          CASE DEFAULT
                            CYCLE
                        END SELECT
                    end do
                    backspace(tmp_file)

                CASE('WELLTEMP','EVENTEMP','GEOM','FAMILY','DUALFAMILY')
                    CALL ADD_FAM(tmp_file,linelen,tmpbas)

                    read(tmp_file,*,iostat=ios) iqm
                    if (ios.ne.0 .or. iqm==0) then
                      write(lupri,*) "ERROR(readbasis) :",           &
     &                 "No iqm for ",tmpbas%type," specified."
                      cycle
                    end if
                    ALLOCATE(JCO(iqm))
                    read(tmp_file,*,iostat=ios) iqm,JCO(1:iqm)
                    if (ios.ne.0) then
                      write(lupri,*) "ERROR(readbasis) :",           &
     &                 "Specify at least ",iqm," numbers"
                      cycle
                    end if
                    CALL ADD_JCO(tmpbas,JCO)
                    CALL ADD_BASIS(bas,firstbas,tmpbas)
                    DEALLOCATE(JCO)

                  CASE DEFAULT
                    READ(tmp_file,*)
                    CYCLE
                  END SELECT
                end if

              ELSE IF (section(2)=='*CENTERS'.AND.                      &
     &                 section(3)=='.NUCLEUS') THEN
!               Explicit specification of charges for non-standard nuclei, defaults are taken from nucdata.h
                name=' '
                charge=0.0
                bascharge=0.0
                nucexp=0.0

                read(tmp_file,*,iostat=ios) name, charge
                if (ios.ne.0) then
                    write(lupri,*) "ERROR(readbasis) :",             &
     &          "Invalid name or charge specified in *CENTERS"
                    CALL QUIT ('ERROR(readbasis): .NUCLEUS')
                end if
!               Try to read also element number so that we can use the basis set library (e.g. in
!               counterpoise calculations)
                backspace(tmp_file)
                read(tmp_file,*,iostat=ios) name, charge, bascharge
                if (ios.ne.0) then
!                  User did not give two values, take the element number
!                  equal to the charge
                   bascharge=charge
                   backspace(tmp_file)
                end if

!               Set charge information for all atoms of this type
                mol=>firstatom
                DO
                  mol=>find_atom(mol,name)
                  IF (associated(mol)) THEN
                    mol%charge=charge
                    mol%bascharge=bascharge
                    mol%nucexp=nucexp
                    mol=>mol%next
                  ELSE
                    EXIT
                  END IF
                END DO

              ELSE IF (section(2)=='*COORDINATES'.AND.                  &
     &                 section(3)=='.UNITS') THEN
!               If this block is specified, the coordinates in the .XYZ file are given in these
!               units and we have to transform them to a.u.
!               The conversion factor may be given in Angstrom typed as 'ANGSTROM' or 'A'
!               Note: Angstroms are known and can be specified by 'A' or 'ANGSTROM'
                read(tmp_file,*,iostat=ios) str
                SELECT CASE(str)
                CASE('ANGSTROM', 'A')
                  UNIT_CONVERSION = ONE / XTANG
                  unitname='angstrom'
                CASE('AU','a.u.','bohr')
                  UNIT_CONVERSION = ONE
                  unitname='bohr'
                CASE('PM', 'pm', 'picometer')
                  UNIT_CONVERSION = ONE / (100.D0 * XTANG)
                  unitname='picometer'
                CASE DEFAULT
                  backspace (tmp_file)
                  read(tmp_file,*,iostat=ios) unitname,unit_conversion
                  if (ios.ne.0) then
                    write(lupri,*) "ERROR(readbasis) :",                &
     &    "Specify a unitname and the conversion factor to atomic units"
                    call Quit ("error reading units in .UNITS")
                  end if
                END SELECT

!             End of processing of the **MOLECULE BLOCK
              END IF
            END IF
          END DO

!         We want users to specify the default basis, quit if they don't
          IF(.NOT.(ASSOCIATED(firstbas))) THEN
            WRITE(LUPRI,*) 'No default basisset given in input file.    &
     &            Specify a basisset in  **MOLECULE *BASIS .DEFAULT'
            CALL QUIT('*** ERROR READ_BASIS, no defaultbasis given in   &
     &           input file.')
          END IF

!     We now have add all data that was explicitly added in the .INP file.

!     Convert units
          IF(.NOT.(unitname=='bohr')) THEN
!            WRITE(LUPRI,*) 'Coordinates were entered in ', unitname
!            WRITE(LUPRI,*) 'Converting to atomic units with',           &
!     &                     ' conversion factor=', UNIT_CONVERSION
            mol=>firstatom
            DO WHILE (ASSOCIATED(mol))
              mol%x=mol%x*UNIT_CONVERSION
              mol=>mol%next
            END DO
          END IF

        first_call = .FALSE.
        inMOLECULE = .false.
        close (tmp_file)

        END SUBROUTINE READ_BASIS


        SUBROUTINE WRITE_ATOMS(lupri)
!     JHS050707 Writes the data contained in the list of atoms to the common blocks.

!      use nuclei_mod, only: CHARGE, CORD, GAUNUC, GNUEXP, IATOMTYP, ISOTOP, IZATOM, NAMEX, NAMN, NUCIND
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      IMPLICIT NONE
!     CONSTANTS

!     COMMON BLOCKS

      INTEGER, INTENT(IN)             :: lupri

      INTEGER :: ncent,natom,natom_type,nq
      TYPE(atom), POINTER             :: mol
#include "nuclei.h"
!
      REAL(KIND=8), PARAMETER :: THRMIN = 1.D-15,DSM=1.0D-30
      REAL(KIND=8) :: XMASS

      CALL QENTER('WRITE_ATOMS')

!     Check that we do not exceed the maximum number of centers in the common block
      ncent=count_atoms(firstatom)
      IF (NCENT>MXCENT) THEN
         WRITE(lupri,*) 'Too many centers for this version of DIRAC'
         WRITE(lupri,*) 'Adjust maxorb.h'
         CALL QUIT('*** ERROR (READIN) ncent>mxcent')
      END IF

      NUCIND=0
      mol=>firstatom
      natom=0
      natom_type=1

      DO

        IF (.NOT.(associated(mol))) THEN
          EXIT
        ELSE
          natom=natom+1
!     write the atomic data.
!     If we didn't read the nuclear model somewhere read the nuclear radius.

!    WRITE NUCLEAR SHAPE. TMP To be changed, now always D0 is taken
          nq = NINT(mol%charge)
          IF (GAUNUC.AND.nq.GT.0) THEN

            IF (mol%nucexp.EQ.D0) THEN
              CALL NUCSIZ(nq,mol%nucexp,xmass)
            END IF
            IF (mol%nucexp .LE. 1.0D8) WRITE(lupri,*)                   &
     &         'WARNING!!! Nuclear gaussian exponent small: ',mol%nucexp
          ELSE
            mol%nucexp=D0
          END IF

!     WRITE COMMON BLOCKS
          NAMN(natom)=mol%name
          CORD(:,natom)=mol%x
          NAMEX(3*natom)=NAMN(natom)//' z'
          NAMEX(3*natom-1)=NAMN(natom)//' y'
          NAMEX(3*natom-2)=NAMN(natom)//' x'
          CHARGE(natom)=mol%charge
          GNUEXP(natom)=mol%nucexp
          if (natom .gt. 1) then
!            in .xyz input different basis sets are
!            specified with .SPECIAL names,
!            therefore we test on change in name rather
!            than change in charge.
             if ( NAMN(natom) .ne. NAMN(natom-1) ) then
                natom_type = natom_type + 1
             end if
          end if
          IATOMTYP(natom)=natom_type
          IZATOM(natom)  = nq ! number of protons in nucleus
          ISOTOP(natom)  = 1  ! assume most abundant isotope

          mol=>mol%next
        END IF
      END DO

      NUCIND=natom

      CALL QEXIT('WRITE_ATOMS')

      END SUBROUTINE WRITE_ATOMS

        SUBROUTINE POS_FILE(infile,lupri,linelen,section,subsection,name)
!         Position a file to a specific keyword.

          INTEGER, INTENT(IN)   :: infile,lupri,linelen
          CHARACTER(LEN=*), INTENT(IN) ::section,subsection,name
          CHARACTER(LEN=linelen) :: str,pressect,pressub
          character(len=linelen) :: tmp_name
          INTEGER                :: ios

          REWIND(infile,iostat=ios)
! miro: insert this control printout
          if (ios.ne.0) then
       write(lupri,*) 'POS_FILE: error in rewind command ! ios=',ios
       write(lupri,*) 'POS_FILE: infile(in)=',infile
       write(lupri,*) 'POS_FILE: linelen(in)=',linelen
       write(lupri,*) 'POS_FILE: section(in)=',section
       write(lupri,*) 'POS_FILE: subsection(in)=',subsection
       write(lupri,*) 'POS_FILE: name(in)=',name
           call quit('POS_FILE: ios<>0, error in rewind command !')
          endif
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
             CALL QUIT('*** ERROR KEY NOT FOUND. ***')
            ELSE IF (str(1:2)=='**') THEN
              pressect=str
            ELSE IF (str(1:1)=='*') THEN
              pressub=str
            ELSE IF (pressect==section.AND.pressub==subsection) THEN
              IF (str==name) EXIT
            END IF
          END DO
        END SUBROUTINE POS_FILE

        SUBROUTINE WRITE_ATOMTYPES(infile,lupri,WORK,LWORK,KANG,KSETS,KBLOCK,KPRIM,HERMIT)
!     JHS050707 Writes the contained in the basisset list to the common blocks.
!               Note the routine is higly similar to the subroutine READI1 for
!               .MOL files, the main difference is the split between reading the files
!               and writing them to the common blocks

!      use dcbgen_mod, only: TWOCOMP
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use nuclei_mod, only: NAMN, NONT, NONTYP, NUCDEP
!      use symmet_mod, only: MAXREP
      use maxorb_mod, only: MXCORB
!      use cbirea_mod, only: FAMPAR, INPTST, IPREAD, NFAMEXP, TOL_SYMADD, UNCONT
!      use ccom_mod, only: DOCART, DOSPHE, NHTYP, THRS
!      use dgroup_mod, only: ATOMIC, LINEAR, NONATM
      use maxaqn_mod, only: MXAQN, MXQN
!      use infpar_mod, only: MYTID
      use RECP_NTR
      IMPLICIT NONE
!     CONSTANTS
#include "symmet.h"

!     COMMON BLOCKS
! cbirea.h : UNCONT, IPREAD, TOL_SYMADD, ...
! infpar.h : mytid
#include "cbirea.h"
#include "infpar.h"
      REAL(KIND=8), PARAMETER :: THRMIN = 1.D-15,DSM=1.0D-30

      INTEGER,INTENT(IN) :: lupri,infile,KANG,KSETS,KBLOCK,KPRIM
      INTEGER,INTENT(IN) :: LWORK

      INTEGER      :: KLAST,i,j,IBSFLAG, KAOVEC, NBLOCK,                &
     &                katom,basindex,index, iqm_check

      REAL(KIND=8) :: QEXP=D0
!Miro: MXSPD is integer, not real, and that was causing problems later with parallel-GNU-compilers !!!
      INTEGER      :: MXSPD
      LOGICAL      :: HERMIT,DOOWN,ANG,NOORBTS,NOFITFS,                 &
     &                CNTBAS=.FALSE.,DIRCON(2)=.FALSE.,IFXYZ(3)
      CHARACTER(LEN=72) :: TTITLE(2)=' '
      CHARACTER(LEN=11) :: SYMTXT
      CHARACTER(LEN=15) :: CLASS
      CHARACTER(LEN=namelen) :: prevname

      TYPE(basistype),POINTER       :: baspointer
      TYPE(atom),POINTER            :: mol
      TYPE(atomtype),POINTER        :: typepointer

      REAL(KIND=8):: WORK(LWORK)
      REAL(KIND=8),ALLOCATABLE::                                        &
     &               ALPHA(:,:,:),                                      &
     &               CPRIM(:,:,:,:),                                    &
     &               CPRIMU(:,:,:,:)
      INTEGER,ALLOCATABLE :: JCO(:,:,:),JCO2(:,:),IBLOCK(:),            &
     &                       NUC(:,:),IQM(:,:),NRC(:,:),NBLCK(:,:),     &
     &                       ISGEN(:)
      LOGICAL,ALLOCATABLE :: SEG(:,:)
      CHARACTER(LEN=80),ALLOCATABLE :: BASREF(:,:,:)
      CHARACTER(LEN=80)             :: BASNAM

!     JHS050707 For type conversion purposes.
      INTEGER :: int8
      integer blk
!
#include "ccom.h"
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dgroup.h"
#else
!#include "gnrinf.h"
      LOGICAL LINEAR, ATOMIC
#endif

      katom=count_types(firstatomtype)

      ALLOCATE(JCO(KANG,KATOM,KSETS))
      ALLOCATE(JCO2(KANG,KATOM))
      ALLOCATE(IBLOCK(MXBSETS))
      ALLOCATE(NUC(KBLOCK,KSETS))
      ALLOCATE(IQM(KATOM,KSETS))
      ALLOCATE(NRC(KBLOCK,KSETS))
      ALLOCATE(NBLCK(KATOM,KSETS))
      ALLOCATE(ISGEN(KBLOCK))
      ALLOCATE(SEG(KBLOCK,KSETS))
      ALLOCATE(BASREF(10,KATOM,KSETS))
      ALLOCATE(ALPHA(KPRIM,KBLOCK,KSETS))
      ALLOCATE(CPRIM(KPRIM,KPRIM,KBLOCK,KSETS))
      ALLOCATE(CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS))

!     JHS050707 TMP Initialization (NEEDED?)

      NONT=0
      JCO=0
      JCO2=0
      IBLOCK=0
      NUC=0
      IQM=0
      NRC=0
      NBLCK=0
      ISGEN=0
      SEG=.FALSE.
      BASREF=' '
      IBLOCK(1:4)=1

!     ccom.h
      THRS=THRMIN
!     JHS050707 TMP change this for other coordinte systems
      DOCART=.TRUE.
      DOSPHE=.FALSE.
      DOOWN=.FALSE.
      NHTYP=0

!     nuclei.h
      CALL SPHINP(infile, WORK, LWORK, DOOWN, MXSPD)

      typepointer=>firstatomtype
      basindex=0
      index=0
      NONTYP=0

      DO
!     Loop over all atomtypes
        NOORBTS=.TRUE.
        NOFITFS=.TRUE.

        IF (.NOT.(associated(typepointer))) EXIT
        basindex=basindex+1
        NONTYP=NONTYP+1
        NONT(basindex)=typepointer%natom

!       easier to use the original pointers now
        baspointer=>typepointer%basis
        mol=>typepointer%first

        BASREF(1,basindex,1) = ' Reference not found in input'
        BASREF(2:10,basindex,1:3) = ' Not initialized'

!       Clear the input buffer (if any) of the RDLINE functionality
!       This makes sure a new line is read when RDLINE is called next
        CALL RDLINE(infile)
!       POSITION the input file.
        IF (baspointer%name == 'DEFAULT') THEN
           CALL POS_FILE(infile, lupri, linelen,'**MOLECULE', '*BASIS', '.DEFAULT')
           CALL RDLINE(infile)
        ELSE
           CALL POS_FILE(infile, lupri, linelen,'**MOLECULE', '*BASIS', baspointer%name)
        ENDIF

        SELECT CASE (baspointer%type)
        CASE ('NOBASIS','POINTCHARGE')
          IQM(basindex,1)=0
          JCO(1,basindex,1)=0
          BASNAM='Explicit'
          BASREF(1,basindex,1) =                                        &
     &          "No basis functions (point charge)"
        CASE ('INTGRL','EXPLICIT')
          IBSFLAG=1
          IQM(basindex,1)=SIZE(baspointer%jco)
          DO i=1,IQM(basindex,1)
            JCO(i,basindex,1)=baspointer%jco(i)
          END DO
          BASNAM='Explicit'
          KAOVEC = KBLOCK + 1 - IBLOCK(1)
!     Useful check : read iqm once more
            READ (infile,*) IQM_CHECK
            IF (IQM_CHECK .NE. IQM(basindex,1))                         &
     &         CALL QUIT ('Error in write_atomtypes')
          CALL GTOINP(infile,IQM(basindex,1),JCO(1,basindex,1),         &
     &         NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),                       &
     &         SEG(IBLOCK(1),1),                                        &
     &         ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),             &
     &         CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),                &
     &         NBLOCK,KAOVEC,KPRIM)
          NBLCK(basindex,1) = NBLOCK
          BASREF(1,basindex,1) =                                        &
     &          "Basis set typed explicitly in input file "
        CASE ('MOLFBAS')
          IBSFLAG=2
          BASNAM = baspointer%file
          KLAST = 1 + (1+KPRIM)*KPRIM
          IF (KLAST.GT.LWORK)                                           &
     &      CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
          DIRCON(1) = .TRUE.
          CALL BFGINP(basindex,IBLOCK(1),1,WORK(1),WORK(1+KPRIM),       &
     &         IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,                  &
     &         NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,                     &
     &         baspointer%file,BASREF(1,1,1))
        CASE ('BASIS')
          IBSFLAG=3
          BASNAM = baspointer%file
          KAOVEC = KBLOCK + 1 - IBLOCK(1)

          NBLOCK=basindex+1
          CALL BASLIB(IQM(basindex,1),JCO(1,basindex,1),                &
     &             NUC(IBLOCK(1),1),                                    &
     &             NRC(IBLOCK(1),1),SEG(IBLOCK(1),1),                   &
     &             ALPHA(1,IBLOCK(1),1),                                &
     &             CPRIM(1,1,IBLOCK(1),1),CPRIMU(1,1,IBLOCK(1),1),      &
     &             NBLOCK,                                              &
     &             KAOVEC,KPRIM,mol%bascharge,mol%bascharge,            &
     &             DSM,UNCONT,BASNAM,BASREF(1,basindex,1),IPREAD)
          NBLCK(basindex,1) = NBLOCK
        CASE ('WELLTEMP','EVENTEMP','GEOM')
          IBSFLAG=4
          CALL GENFAMEXP(infile,0)
          i=SIZE(baspointer%fam)-1
          FAMPAR(1:i)=baspointer%fam(1:i)
          NFAMEXP(1)=TRANSFER(baspointer%fam(i+1),int8)
          KAOVEC = KBLOCK + 1 - IBLOCK(1)
          CALL FAMBAS(infile,IQM(basindex,1),JCO(1,basindex,1),         &
     &         NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),                       &
     &         SEG(IBLOCK(1),1),                                        &
     &         ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),             &
     &         CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),                &
     &         NBLOCK,KAOVEC,KPRIM)
          NBLCK(basindex,1) = NBLOCK
          BASREF(1,basindex,1) =                                        &
     & "Well-tempered basis set typed explicitly in input file "
        CASE ('FAMILY','DUALFAMILY')
          IBSFLAG=4
          IF (baspointer%type=='FAMILY') i=1
          IF (baspointer%type=='DUALFAMILY') i=2
          CALL GENFAMEXP(infile,i)
          NFAMEXP(1:i)=TRANSFER(baspointer%fam(1:i),int8)
          KAOVEC = KBLOCK + 1 - IBLOCK(1)
          CALL FAMBAS(infile,IQM(basindex,1),JCO(1,basindex,1),         &
     &         NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),                       &
     &         SEG(IBLOCK(1),1),                                        &
     &         ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),             &
     &         CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),                &
     &         NBLOCK,KAOVEC,KPRIM)
          NBLCK(basindex,1) = NBLOCK
          BASREF(1,basindex,1) =                                        &
     & "Well-tempered basis set typed explicitly in input file "
        CASE DEFAULT
          WRITE(lupri,*) 'BASISSET specification not implemented'
          CALL QUIT('ERROR: READING BASISSET')
        END SELECT

        NOORBTS=NOORBTS.AND.IQM(basindex,1).EQ.0
        NHTYP  =   MAX(NHTYP,IQM(basindex,1))

        DO i = 1, NBLCK(basindex,1)
           IF ((NUC(i,1).NE.NRC(i,1)).AND.(NRC(i,1).GT.0))              &
     &          CNTBAS = .TRUE.
        END DO
        IF ( CNTBAS .AND. (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
           DIRCON(1) = .TRUE.
        END IF

!       Explicit SC will probably not work, worry about that later
        IF(.NOT.TWOCOMP) THEN
          IF (ASSOCIATED(baspointer%next)) THEN
            IF (baspointer%next%name==baspointer%name)                    &
     &            baspointer=>baspointer%next
          END IF

          IF (.NOT.(baspointer%comp(1:5)=='SMALL')) THEN
!       We should typically go via this routine
            IBSFLAG=4
!
!           When uncontracted we generate the small component basis
!           using the kinetic balance relation
!           ...or...
!           Make uncontracted small component basis
!           set from kinetic balance.
!           =======================================================
!
            CALL KINBAL(basindex,IBLOCK(1),IBLOCK(2),                     &
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,                  &
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,                     &
     &           INPTST,CNTBAS)

            BASREF(1,basindex,2) = ' Derived from large component'

          ELSE IF (baspointer%comp=='INTGRL'.OR.                          &
     &           baspointer%comp=='EXPLICIT') THEN
            IBSFLAG=1
            IQM(basindex,2)=SIZE(baspointer%jco)
            DO i=1,IQM(basindex,2)
              JCO(i,basindex,2)=baspointer%jco(i)
            END DO
            BASNAM=baspointer%file
            KAOVEC = KBLOCK + 1 - IBLOCK(2)
!     Additional check : read iqm once more
            READ (infile,*) IQM_CHECK
            IF (IQM_CHECK .NE. IQM(basindex,2))                           &
     &         CALL QUIT ('Error in write_atomtypes')
            CALL GTOINP(infile,IQM(basindex,2),JCO(1,basindex,2),         &
     &           NUC(IBLOCK(2),2),NRC(IBLOCK(2),2),                       &
     &           SEG(IBLOCK(2),2),                                        &
     &           ALPHA(1,IBLOCK(2),2),CPRIM(1,1,IBLOCK(2),2),             &
     &           CPRIMU(1,1,IBLOCK(2),2),ISGEN(IBLOCK(2)),                &
     &           NBLOCK,KAOVEC,KPRIM)
            NBLCK(basindex,2) = NBLOCK
            BASREF(1,basindex,2) =                                        &
     &            "Basis set typed explicitly in input file "
          ELSE IF (baspointer%type(1:7)=='MOLFBAS') THEN
            IBSFLAG=2
            BASNAM=baspointer%file
            KLAST = 1 + (1+KPRIM)*KPRIM
            IF (KLAST.GT.LWORK)                                           &
     &        CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
            DIRCON(1) = .TRUE.
            CALL BFGINP(basindex,IBLOCK(1),1,WORK(1),WORK(1+KPRIM),       &
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,                  &
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,                     &
     &           baspointer%file,BASREF(1,1,1))

          ELSE IF (.NOT.(baspointer%type(1:6)=='KINBAL')) THEN
            IBSFLAG=4
            CALL KINBAL(basindex,                                         &
     &           IBLOCK(1),IBLOCK(2),                                     &
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,                  &
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,                     &
     &           INPTST,CNTBAS)
            BASREF(1,basindex,2) = ' Derived from large component'
          END IF

          NOORBTS = NOORBTS.AND.IQM(basindex,2).EQ.0

          NHTYP = MAX(NHTYP,IQM(basindex,2))

          IF (IQM(basindex,2).GT.0) THEN
            NHTYP     = MAX(NHTYP,IQM(basindex,2))
          ENDIF

          DO i = 1, NBLCK(basindex,2)
            IF ((NUC(i,2).NE.NRC(i,2)).AND.(NRC(i,2).GT.0))               &
     &             CNTBAS = .TRUE.
          END DO
          IF ( CNTBAS .AND. IBSFLAG .NE. 4 ) DIRCON(2) = .TRUE.

        END IF

        IBLOCK(1:2)=IBLOCK(1:2)+NBLCK(basindex,1:2)
        NOORBT(index+1:NONT(basindex)+index) = NOORBTS
        DO i=3,4
!     JHS060707 TMP Include also the reading of fitsets.
          IBSFLAG=2
          IBLOCK(i)=IBLOCK(i)+NBLCK(basindex,i)
        END DO

!       This routine reads the next line to check whether an ECP is specified for this type
        CALL RECP_LNK_READCP(infile,basindex,NONT(basindex),IPREAD)

        typepointer=>typepointer%next
        index=index+NONT(basindex)

      END DO



      LINEAR = .FALSE.
      ATOMIC = .FALSE.
      CLASS = 'N/A'
      IF (ADDSYM) THEN
         CALL SYMADD(WORK,LWORK,NSYMOP,KATOM,KASYM,                     &
     &               CLASS,TOL_SYMADD,IPREAD,.false.)
         IF (CLASS(3:4).EQ.'oo') LINEAR = .TRUE.
      END IF
!
!     ****************************************************
!     ***** Adds an STO-3G basis if Extended Hueckel *****
!     *****        wanted as starting orbitals       *****
!     ****************************************************
!
!      IF (ADDSTO) CALL STOADD()
!
!     ***************************************
!     ***** Process symmetry input data *****
!     ***************************************
!
!     Fix a problem with ADDSYM for some D(2) cases
      CALL FIX_DUPLICATE_GENERATORS (NSYMOP,KASYM)
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)

!
!     ***************************************************
!     ***** Process orbital and geometry input data *****
!     ****************************************************
!
      CALL BASPRO(WORK,LWORK,NSYMOP,IQM,NBLCK,JCO,NUC,NRC,              &
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,KBLOCK,KPRIM, &
     &            DOOWN, blk)
      IF (LINEAR.and. (NUCDEP .Eq. 1) .and. (.NOT. NONATM) ) ATOMIC = .TRUE.
!
!     *************************************************
!     ***** Print orbital and geometry input data *****
!     *************************************************
!
#ifdef BUILD_GEN1INT
!     initialize basis sets (large and small components) used in Gen1Int interface
!     NB! workers obtain gen1int info from gen1int's own mpi calls
      if (mytid .eq. 0)                                                 &
     &call gen1int_host_init(2, NONTYP, KATOM, NONT, IQM,               &
     &                       NBLCK, KANG, JCO, KBLOCK,                  &
     &                       NUC, NRC, KPRIM, ALPHA, CPRIMU)
#endif

      CALL BASOUT(WORK,LWORK,IQM,NBLCK,JCO,BASREF,                      &
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,      &
     &            KBLOCK,KPRIM,HERMIT, blk)
!
!     TMP copy input to .out file.

!
!     ****************************************
!     ***** Output on LUONEL *****************
!     ****************************************
!

      IF (HERMIT) CALL WRONEL(TTITLE,IQM,IFXYZ,KATOM,JCO2,KANG)

          CALL RECP_LNK_CHECKMOL(DOCART)
          CALL RECP_LNK_IRREP(MAXREP)
          CALL RECP_LNK_RDGEO(SYMTXT,NSYMOP,KASYM,blk,                    &
     &            IQM,JCO,KATOM,KANG)


          DEALLOCATE(JCO)
          DEALLOCATE(JCO2)
          DEALLOCATE(IBLOCK)
          DEALLOCATE(NUC)
          DEALLOCATE(IQM)
          DEALLOCATE(NRC)
          DEALLOCATE(NBLCK)
          DEALLOCATE(ISGEN)
          DEALLOCATE(BASREF)
          DEALLOCATE(SEG)
          DEALLOCATE(ALPHA)
          DEALLOCATE(CPRIM)
          DEALLOCATE(CPRIMU)

        END SUBROUTINE WRITE_ATOMTYPES

        SUBROUTINE PRINT_BASIS(firstbas, lupri)
          !TYPE(basistype), POINTER, INTENT(IN) :: firstbas
          TYPE(basistype), POINTER :: firstbas
          INTEGER, INTENT(IN) :: lupri
          TYPE(basistype), POINTER :: bas

          CALL TITLER('Basisset read','*',118)
          WRITE(lupri,'(a,t9,a,t41,a,t73,a)')                           &
     &         'Label','component', 'input type','file'
          bas=>firstbas
          DO
            IF (associated(bas)) THEN
              WRITE(lupri,'(a,t9,a,a,a)')                               &
     &          TRIM(bas%name),bas%comp,bas%type,bas%file
              IF(associated(bas%jco).AND.associated(bas%fam))            &
     &           WRITE(lupri,*) 'Exponents:', bas%jco, bas%fam
              bas=>bas%next
            ELSE
              EXIT
            ENDIF
          END DO
        WRITE(lupri,'(a)') ' '

        END SUBROUTINE PRINT_BASIS

        SUBROUTINE DEL_ATOMTYPES(ftype)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.

          TYPE(atomtype),POINTER :: ftype, typepointer

          IF (.NOT.associated(ftype)) RETURN

          typepointer=>ftype
          DO
            IF (associated(typepointer%next)) THEN
              ftype=>typepointer%next
              DEALLOCATE(typepointer)
              typepointer=>ftype
            ELSE
              NULLIFY(ftype)
              EXIT
            ENDIF
          END DO

        END SUBROUTINE DEL_ATOMTYPES

        SUBROUTINE DEL_BASIS(firstbas)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.

          TYPE(basistype),POINTER :: firstbas
          TYPE(basistype), POINTER :: bas

          IF (.NOT.associated(firstbas)) RETURN

          bas=>firstbas
          DO
            IF (associated(bas%next)) THEN
              firstbas=>bas%next
              IF (ASSOCIATED(bas%jco)) DEALLOCATE(bas%jco)
              IF (ASSOCIATED(bas%fam)) DEALLOCATE(bas%fam)
              DEALLOCATE(bas)
              bas=>firstbas
            ELSE
              NULLIFY(firstbas)
              exit
            ENDIF
          END DO

        END SUBROUTINE DEL_BASIS

        SUBROUTINE CLEAR_XYZ
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.

          CALL DEL_ATOMS(firstatom)
          CALL DEL_BASIS(firstbas)
          CALL DEL_ATOMTYPES(firstatomtype)

        END SUBROUTINE CLEAR_XYZ


        SUBROUTINE READ_XYZ(molfile,cmdfile, lupri, &
     &                      work,lwork,kang,ksets,kblock,kprim,hermit)

        INTEGER, INTENT(IN) :: molfile,cmdfile,lupri
        INTEGER, INTENT(IN) :: lwork
        REAL(KIND=8)        :: WORK(LWORK)
        INTEGER, INTENT(IN) :: KANG,KSETS,KBLOCK,KPRIM
        LOGICAL, INTENT(IN) :: HERMIT

        INTEGER             :: inpfile
        INTEGER             :: ios

!       Driver for reading coordinates and basis set information from xyz file and DIRAC.INP

!       In the XYZ format coordinates are supposed to start at the third line
        REWIND (molfile)
!       Check whether the first line has a number that can be used to check the number of atoms
        natoms_given = 0
        READ (molfile,*,iostat=ios) natoms_given
!       Skip the comment line
        CALL RDLINE(molfile)
        CALL READ_ATOMS(molfile,linelen,firstatom,lupri)
!       Check number of atoms
        IF (natoms_given .NE. count_atoms(firstatom)) THEN
          WRITE(LUPRI,*) " Read ",count_atoms(firstatom)," atoms",      &
     &                   " but expected ",natoms_given," atoms"
          CALL QUIT(' *** ERROR (READ_XYZ) wrong number of atoms')
        END IF

!       Read the basis set information from DIRAC.INP
        inpfile=OPEN_FILE('DIRAC.INP',cmdfile,'FORMATTED',LUPRI)
        IF (.NOT.(inpfile>0)) THEN
          WRITE(LUPRI, '(//a)') 'ERROR, could not open DIRAC.INP'
          CALL QUIT(' *** ERROR (READ_XYZ) unable to open DIRAC.INP')
        END IF
        CALL READ_BASIS(inpfile, linelen,lupri)

!       Define the list of atomtypes
        CALL MAKE_ATOMTYPES(lupri)

!       Write the coordinates to the common blocks
        CALL WRITE_ATOMS(LUPRI)
!       Write the basisset data to the common blocks
        CALL WRITE_ATOMTYPES(inpfile, LUPRI,WORK,LWORK,                 &
     &                      KANG,KSETS,KBLOCK,KPRIM,HERMIT)

!       Free the reserved memory.
        CALL CLEAR_XYZ()
        CLOSE(inpfile)

        END SUBROUTINE READ_XYZ

      END MODULE HERMODULES
!  -- End of abacus/hermodules.F90 --
