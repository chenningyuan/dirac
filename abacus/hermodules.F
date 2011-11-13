!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

!     !!! MODULE FOR HANDLING FILES !!!
!     Contains
!       - OPEN_FILE(filename) or OPEN_FILE(filename, filenum,formatted,output)
!         Opens a file and returns the id. Which is a new number if only the
!         name is given, otherwise it is equal to filenum.
!         A negative or zero indicates an error.
!       - GET_EXTENSION(filename)
!         Returns the first extension in the array filetypes for which
!         filename.extension is an existing filename.
!    
!     JHS 05-07-07

      MODULE FILE
      IMPLICIT NONE

      PRIVATE
      PUBLIC OPEN_FILE, GET_EXTENSION, filter_comments

        CHARACTER(LEN=8), DIMENSION(2), PARAMETER :: filetypes=(/'MOL', &
     &                                                           'XYZ'/)
        INTEGER :: linewidth=0

        INTERFACE OPEN_FILE
          MODULE PROCEDURE OPEN_FILE_NAME, OPEN_FILE_NUM, OPEN_SCRATCH
        END INTERFACE
        
        
      CONTAINS

        INTEGER FUNCTION OPEN_FILE_NAME(filename)
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0  -   File does not exists.
!     -1  -   To many files opened.
!     -99  -   Error while executing the open statement.
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER, parameter :: maxfile=99
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
              open_file_name=-1 ! To many files opened
            ELSE
               open_file_name=n ! Open file
            END IF                  
          ELSE
            open_file_name=0 ! Filename does not exist
          END IF

        END FUNCTION OPEN_FILE_NAME

        INTEGER FUNCTION OPEN_SCRATCH()
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0  -   File does not exists.
!     -1  -   To many files opened.
!     -99  -   Error while executing the open statement.
          INTEGER, parameter :: maxfile=99
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
            open_scratch=-1 ! To many files opened
          ELSE
             open_scratch=n ! Open file
          END IF                  

        END FUNCTION OPEN_SCRATCH

        INTEGER FUNCTION OPEN_FILE_NUM(filename, filenum, formatted,    &
     &                                 lupri)
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0   -   File does not exists.
!     -1   -   To many files opened.
!     -2   -   Number already assigned
!     -3   -   File already opened with a different number.
!     -99  -   Error while executing the open statement.
          CHARACTER(LEN=*), INTENT(IN)  :: filename
          INTEGER, INTENT(IN)           :: filenum
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: formatted
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          INTEGER, PARAMETER            :: maxfile=99
          LOGICAL                       :: ex
          INTEGER                       :: n

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

        END FUNCTION OPEN_FILE_NUM

        FUNCTION GET_EXTENSION(filename)
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
          IF (i>size_ft)  get_extension=' '
        END FUNCTION
               
        subroutine filter_comments(filein,fileout)
          integer             :: filein,fileout
          character(len=192)  :: str
          character(len=1)    :: chr
          integer             :: i,ios

          rewind(filein)
          fileout=OPEN_SCRATCH()

          do 
            read(filein,'(192a)',iostat=ios) str
            if (ios.ne.0) exit
            do i=1,len(str)
              select case (str(i:i))
              case ('#','!')               
                exit
              end select
            end do      
            write(fileout,*) str(1:i-1)
          end do
          rewind(filein)
          rewind(fileout)
          

        end subroutine

      END MODULE FILE


!     !!! MODULE FOR READING XYZ FILES !!!
!     Module for reading xyz files. 
!     The coordinates are given in the .xyz file.
!     The basisset and possible charges and nuclear exponents are 
!     specified in the .inp file
!
!     In the long run the readingsection for the .inp file move to dirac/dirrdn.F
!     JHS 05-07-07

      MODULE READFILE
        USE FILE
        IMPLICIT NONE

        PRIVATE
        PUBLIC READ_XYZ,WRITE_ATOMS,WRITE_BASIS,CLEAR_XYZ
!     DETERMINES HOW DIFFERENT MOLECULE FILES SHOULD BE READ.
!     !!! PAM CAPITALIZES THE EXTENSIONS !!!

        INTEGER, PARAMETER :: namelen=8, basislen=32,linelen=192
        REAL(KIND=8), PARAMETER :: D0 = 0.0D0
        INTEGER            :: i,ios

        TYPE atom
!     Defines a list for storing atomic data, like coordinates, charges and nuclear exponents.
          CHARACTER(LEN=namelen)  :: name=' '
          REAL(KIND=8)            :: x(3)=(/(0.0,i=1, 3)/)
          REAL(KIND=8)            :: charge=0.0
          REAL(KIND=8)            :: nucexp=0.0
          TYPE(atom), POINTER     :: next
        END TYPE atom
        TYPE(atom), POINTER :: firstatom

        TYPE basistype
!     Defines a list for storing the information about the basis sets.
          CHARACTER(LEN=namelen)     :: name=' '
          CHARACTER(LEN=basislen)    :: type=' '
          CHARACTER(LEN=basislen)    :: comp=' '
          CHARACTER(LEN=basislen)    :: file=' '
!         The following can be assigned for explicitely given basissets or Family type.
          INTEGER, POINTER           :: JCO(:)
          REAL, POINTER              :: FAM(:)
          TYPE(basistype), POINTER   :: next
        END TYPE basistype
        TYPE(basistype), POINTER :: firstbas

!     Global variables
!     Symmetry of the molecule
        LOGICAL :: ADDSYM=.TRUE.
        INTEGER :: NSYMOP=0
        CHARACTER(LEN=1) :: KASYM(3,3)=' '



      CONTAINS


        SUBROUTINE ADD_ATOM(prevmol, fatom, name,x)
!     JHS050707 Adds a new node to the atom list.
          TYPE(atom), POINTER :: prevmol,fatom,mol
          CHARACTER(LEN=*), INTENT(IN)     :: name
          REAL(KIND=8),INTENT(IN),OPTIONAL :: x(1:3)

          ALLOCATE(mol)
          mol%name=name
          IF (PRESENT(x)) mol%x=x

          mol%next=>NULL()
          IF (associated(prevmol)) THEN
            prevmol%next=>mol
          ELSE
            fatom=>mol
          END IF
          prevmol=>mol
        END SUBROUTINE ADD_ATOM

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
          bas%jco=>NULL()
          bas%fam=>NULL()
        END SUBROUTINE ADD_BASIS

        SUBROUTINE ADD_JCO(bas,JCO)
!     JHS050707 Assigns an array with JCO components.
          TYPE(basistype)    :: bas
          INTEGER            :: JCO(:)
          INTEGER,TARGET,ALLOCATABLE,SAVE :: JCO2(:)
          
          ALLOCATE(JCO2(SIZE(JCO)))
          JCO2=JCO
          bas%jco=>JCO2
        END SUBROUTINE

        SUBROUTINE ADD_FAM(filenum,linelen,bas)
!     JHS050707 Assigns for different types of family basissets an array FAM. 
!               And reads the needed coefficients from file.
          INTEGER,INTENT(IN) :: filenum,linelen
          TYPE(basistype)    :: bas
          REAL,TARGET,ALLOCATABLE,SAVE :: FAM(:)
          INTEGER               :: i,ios
#include "priunit.h"

! JHS test, iqm needs to be included as well!
          SELECT CASE(bas%type)
          CASE('WELLTEMP')
            ALLOCATE(FAM(5))
            read(filenum,*,iostat=ios) FAM(1:5)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ", 
     &          "Specify 5 coefficients for WELLTEMP basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('EVENTEMP')
            ALLOCATE(FAM(3))
            read(filenum,*,iostat=ios) FAM(1:3)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ", 
     &          "Specify 3 coefficients for EVENTEMP basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('GEOM')
            ALLOCATE(FAM(5))
            read(filenum,*,iostat=ios) FAM(1:2), FAM(5)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ", 
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
              write(lupri,*) "ERROR(add_fam) : ", 
     &          "Specify 1 coefficients for FAMILY basis."
              call exit(1)
            end if
            bas%fam=FAM
          CASE('DUALFAMILY')
            ALLOCATE(FAM(2))
            read(filenum,*,iostat=ios) FAM(1:2)
            if (ios.ne.0) then
              write(lupri,*) "ERROR(add_fam) : ", 
     &          "Specify 2 coefficients for DUAL basis."
              call exit(1)
            end if
            bas%fam=FAM
          END SELECT
        END SUBROUTINE

        SUBROUTINE READ_ATOMS(filenum, linelen,fatom,lupri)
!     JHS050707 Reads atom name and coordinates from the .XYZ file 
!               and writes it to the atoms list.
#include "infpar.h"
#include "cbirea.h"
          TYPE(atom), POINTER :: fatom
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
          TYPE(atom), POINTER                :: mol
          CHARACTER(len=linelen)             :: str
          CHARACTER(LEN=namelen)             :: name
          REAL(KIND=8)                       :: x(3)
          INTEGER                            :: ios
 
          CALL QENTER('READ_ATOMS')

          mol=>null()
          IF (.NOT.(PARHER).AND.IPREAD>1)                               
     &         WRITE(LUPRI,*) 'Reading coordinates from xyz file...'

          ios=0
          do while(.not.(ios.ne.0))

            name='init'
            x(1:3)=0.0

            READ(filenum, *,iostat=ios) name, x(1), x(2), x(3)
            IF (IPREAD >= 4) THEN
              write(lupri,*) 'READ_ATOMS: read coordinates:',           
     &        x(1),x(2),x(3)," error code = ",ios
            ENDIF

c           todo : the current reading relies on ios giving an error for no-atom lines, improve upon this
            if (ios.eq.0.AND.name.ne.'init') then
              CALL ADD_ATOM(mol,fatom,name,x)  
            end if
          end do


          IF(.NOT.(PARHER).AND. IPREAD>1) WRITE(LUPRI,*) 'Done'

          CALL QEXIT('READ_ATOMS')

        END SUBROUTINE READ_ATOMS

      
        FUNCTION find_atom(fatom, name)
!     JHS050707 Returns a pointer to the next atom in the list with name 
!               or NULL() if is no such atom. 
          ! TYPE(atom), POINTER, INTENT(IN) :: fatom 
           TYPE(atom), POINTER :: fatom 
           TYPE(atom), POINTER             :: find_atom
          CHARACTER(LEN=*), INTENT(IN)     :: name
           TYPE(atom), POINTER             :: mol
           mol=>fatom
           DO 
             IF (associated(mol)) THEN
               IF (mol%name==name) THEN
                 find_atom=>mol
                 EXIT
               END IF
               mol=>mol%next
             ELSE
               find_atom=>NULL()
               EXIT
             ENDIF 
           END DO

        END FUNCTION find_atom


        FUNCTION find_basis(fbas, name)
!     JHS050707 Returns a pointer to the next basisset with name in the list with basissets
!               or NULL() if there is no such basis.
          ! TYPE(basistype), POINTER, INTENT(IN) :: fbas
           TYPE(basistype), POINTER :: fbas
           TYPE(basistype), POINTER                :: find_basis
           TYPE(basistype), POINTER                :: bas
           CHARACTER(LEN=*), INTENT(IN)        :: name

           bas=>fbas
           find_basis=>NULL()
           DO WHILE(ASSOCIATED(bas))
             IF (bas%name==name) THEN
               find_basis=>bas
               EXIT
             END IF
             bas=>bas%next
           END DO

        END FUNCTION find_basis

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

        SUBROUTINE SORT_ATOMS(fatom,fbas,lupri)
!     JHS050707 Sorts the atomlist in the same order as the basissets 
!               delets atoms that don't have a associated basisset
         ! TYPE(atom),     POINTER,INTENT(INOUT) :: fatom
          TYPE(atom),     POINTER :: fatom
         ! TYPE(basistype),POINTER,INTENT(INOUT) :: fbas
          TYPE(basistype),POINTER :: fbas
          INTEGER, INTENT(IN),OPTIONAL       :: lupri
          TYPE(atom), POINTER                :: mol=>NULL()
          TYPE(atom), POINTER                :: newmol=>NULL()
          TYPE(atom), POINTER                :: newfirstatom=>NULL()
          TYPE(basistype), POINTER           :: bas=>NULL(),tmpbas
          CHARACTER(LEN=namelen)             :: name,oldname
          
          mol=>fatom
          bas=>fbas
          newfirstatom=>NULL()
          newmol=>newfirstatom

          DO WHILE(ASSOCIATED(bas))
!     Look wether we have already ordered this species of atoms
            mol=>find_atom(newfirstatom,bas%name)
            IF (ASSOCIATED(mol)) THEN
              bas=>bas%next
              CYCLE
            END IF

            mol=>fatom               
              mol=>find_atom(mol, bas%name)
            DO WHILE (ASSOCIATED(mol))
              CALL ADD_ATOM(newmol,newfirstatom,mol%name,mol%x)
              newmol%charge=mol%charge
              newmol%nucexp=mol%nucexp
              mol=>mol%next
              mol=>find_atom(mol, bas%name)
            END DO
            bas=>bas%next
          END DO
          mol=>fatom
          CALL DEL_ATOMS(mol)
          fatom=>newfirstatom       

        END SUBROUTINE SORT_ATOMS


        SUBROUTINE SORT_BASIS(fatom,fbas,usecomp,lupri)
!     JHS050707 Sorts the basisinput in the same order as the atoms 
!               and issues warnings if no basis can be found
         ! TYPE(atom),     POINTER, INTENT(INOUT) :: fatom
          TYPE(atom),     POINTER :: fatom
         ! TYPE(basistype),POINTER, INTENT(INOUT) :: fbas
          TYPE(basistype),POINTER :: fbas
          INTEGER, INTENT(IN),OPTIONAL       :: lupri
          LOGICAL,INTENT(IN),OPTIONAL        :: usecomp
          LOGICAL                            :: first
          TYPE(atom), POINTER                :: mol=>NULL()
          TYPE(basistype), POINTER           :: bas=>NULL()
          TYPE(basistype), POINTER           :: tmpbas=>NULL()
          TYPE(basistype), POINTER           :: tmpbas2=>NULL()
          TYPE(basistype), POINTER           :: default=>NULL()
          TYPE(basistype), POINTER           :: prevbas=>NULL()
          TYPE(basistype), POINTER           :: newbas=>NULL()
          TYPE(basistype), POINTER           :: newfirstbas=>NULL()
          CHARACTER(LEN=namelen)             :: name,oldname
          CHARACTER(LEN=basislen)            :: comp
          
          mol=>fatom
          bas=>fbas
          newfirstbas=>NULL()
          newbas=>newfirstbas
          prevbas=>newfirstbas
          
!     FIND DEFAULT BASIS.
          default=>find_basis(fbas,'DEFAULT')
          IF (PRESENT(usecomp).AND.usecomp) THEN
            DO WHILE(ASSOCIATED(mol))
              first=.TRUE.
              bas=>fbas
              DO WHILE(ASSOCIATED(bas))
                bas=>find_basis(bas, mol%name)
                IF (first.AND..NOT.(ASSOCIATED(bas))) THEN
                  IF (.NOT.(ASSOCIATED(default)))                       &
     &              CALL QUIT('*** ERROR SORT_BASIS, no defaultbasis    &
     &                    found for atom')

                  CALL ADD_BASIS(newbas, newfirstbas,default)           
                  newbas%name=mol%name
                END IF
                IF (.NOT.(ASSOCIATED(bas))) EXIT
                first=.FALSE.
                tmpbas=>find_basis(newfirstbas,mol%name)
                tmpbas2=>prevbas
                DO
                  IF (.NOT.(ASSOCIATED(tmpbas))) THEN 
                    CALL ADD_BASIS(newbas, newfirstbas,bas)
                    EXIT
                  END IF
                  IF (bas%comp<tmpbas%comp) THEN
                    CALL ADD_BASIS(tmpbas2,newfirstbas,bas,tmpbas)
                    EXIT
                  ELSE IF(.NOT.(tmpbas%name==mol%name)) THEN
                    CALL ADD_BASIS(tmpbas2,newfirstbas,bas,tmpbas)
                    EXIT                          
                  END IF
                  tmpbas2=>tmpbas
                  tmpbas=>tmpbas%next
                END DO                    
                bas=>bas%next
              END DO
              prevbas=>newbas
              oldname=mol%name
              DO WHILE(ASSOCIATED(mol))
                IF (.NOT.(mol%name==oldname)) EXIT
                mol=>mol%next
              END DO
            END DO
          ELSE
            DO WHILE(ASSOCIATED(mol))
              first=.TRUE.
              DO 
                bas=>find_basis(fbas, mol%name)                  
                IF (first.AND..NOT.(ASSOCIATED(bas))) THEN
                  WRITE(lupri,*) 'Atom not found the given basissets.'
                  CALL QUIT('*** ERROR SORT_BASIS, basis                & 
     &                   not found for atom')                  
                END IF
                IF (.NOT.(ASSOCIATED(bas))) EXIT
                CALL ADD_BASIS(newbas,newfirstbas,bas)
                first=.FALSE.
                bas=>bas%next
              END DO                
              mol=>mol%next
            END DO
          END IF
          bas=>fbas
          CALL DEL_BASIS(bas)
          fbas=>newfirstbas

        END SUBROUTINE SORT_BASIS

        REAL FUNCTION FIND_CHARGE(atomname,LUPRI)
!     JHS050707 For a given atomname the corresponding charge is searched for in nucdata.h
          CHARACTER(LEN=*),INTENT(IN) :: atomname
          INTEGER,OPTIONAL,INTENT(IN) :: LUPRI
          INTEGER :: i
#include "nucdata.h"
          
          find_charge=0.0
          DO i=1,SIZE(NUCLABEL)
            IF (atomname==NUCLABEL(i)) THEN
              find_charge=REAL(i)
              RETURN
            END IF
          END DO
          IF (PRESENT(LUPRI)) WRITE(LUPRI,*) 'WARING:                   & 
     &         No charge found for element: ', atomname

        END FUNCTION


        SUBROUTINE READ_BASIS(filenum,linelen,lupri)
        use codata
#include "infpar.h"
#include "cbirea.h"
!     JHS050707 Reads basis names and possibly the coefficients typed in the DIRAC.INP file.
!               Reads also Charges, nuclear exponents and the units of the XYZ file.
!     Stores it in the the basis list.
!     !!! TMP: TESTED FOR: 
!                         - 'LARGE BASIS'
!     (Implemented also for the rest)
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
          TYPE(atom), POINTER                :: mol
          TYPE(atom)                         :: tmpmol
          TYPE(basistype),POINTER            :: bas
          TYPE(basistype)                    :: tmpbas
          integer                            :: tmp_file,ios

!     FOR UNIT CONVERSION
          REAL(KIND=8)      :: UNIT_CONVERSION,ONE=1.0
          CHARACTER(LEN=32) :: UNITNAME='Angstrom'
        
          CHARACTER(LEN=linelen)             :: str,section(3)
          INTEGER                            :: i,j,index,iqm
          CHARACTER(LEN=namelen)             :: name, prevname
          REAL(KIND=8)                       :: charge, nucexp
          LOGICAL                            :: CONT=.FALSE.
          CHARACTER(LEN=9) :: symtxt

!     FOR explicitely typed basissets.
          INTEGER,TARGET,ALLOCATABLE         :: JCO(:)

!     take out comments from the input file and initialize
         call filter_comments(filenum,tmp_file)
          
          mol=>null()
          bas=>null()
          index=1
          section=''

          NSYMOP=0
          ADDSYM=.TRUE.
          KASYM(:,:)=' '

!     Standard for xyz is angstrom, initialize with the AU/ANG conversion factor
          UNIT_CONVERSION = ONE / XTANG


          DO
            read(tmp_file,*,iostat=ios) str

            IF (ios.ne.0) THEN
              EXIT
            ELSE IF (str(1:2)=='**') THEN
              section(1)=str              
            ELSE IF (str(1:1)=='*') THEN
              section(2)=str
!             Reinitialize section after we read an *END
              IF (str=='*END') section=''
            ELSE IF (str(1:1)=='.') THEN
               section(3)=str
            ENDIF

!          
!           We only need to do something if we are reading a **MOLECULE block.
!           The options and suboptions are flattened out in one if then else if construction.
            IF(section(1)=='**MOLECULE') THEN

!            Symmetry specification
              IF (section(2)=='*SYMMETRY'.AND.                         
     &            section(3)=='.NOSYM') THEN
                      ADDSYM=.FALSE.
!             ELSE IF(section(2)=='*SYMMETRY') THEN
!                BACKSPACE(tmp_file)
!                CALL READ_VAR(tmp_file,linelen,NSYMOP,index,skip=.FALSE.)                
!                do i=1,3
!                  do j=1,3
!                    if (index<=0) then
!                      write(LUPRI,*)                                    & 
!     &                   'ERROR:hermodules.F reading symmetry.'
!            CALL QUIT('*** ERROR READ_BASIS, symmetry not               &
!     &                     correctly specified')
!                    CALL READ_VAR(tmp_file,linelen,kasym(i,j),           &
!     &           index,1,.FALSE.)
!                    end if
!                  end do
!                end do

!             Basis specification
              ELSE IF (section(2)=='*BASIS'.AND.
     &                 section(3)=='.DEFAULT') THEN
!               In the subsection .DEFAULT we give only give the name of a 'LARGE BASIS' type
!               like cc-pVDZ etc. This will be the basisset for all atoms that are not
!               specified in the .SPECIAL section.
                read(tmp_file,*,iostat=ios) tmpbas%file
                if (.not.(ios.ne.0) .and. tmpbas%file>' ') then         
                  tmpbas%name='DEFAULT'
                  tmpbas%comp='LARGE'
                  tmpbas%type='BASIS'
                  tmpbas%next=>NULL()
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  CYCLE
                else if (ios.ne.0) then
                  write(lupri,*) "ERROR(readbasis): ",                  
     &                           "Expecting default basisset name."

                  exit
                end if
              ELSE IF(section(2)=='*BASIS'.AND.                         
     &                section(3)=='.SPECIAL') THEN
!     JHS050707 For some or all atoms we can give an other basisset than the .DEFAULT
                tmpbas%name=' '
                tmpbas%comp=' '
                tmpbas%type=' '
                tmpbas%file=' '
                tmpbas%next=>NULL()
                
                read(tmp_file,*,iostat=ios) tmpbas%name,tmpbas%comp,
     &                                      tmpbas%type
                if (ios.ne.0) then
                  backspace(tmp_file)
                  read(tmp_file,*,iostat=ios) tmpbas%name
                  if(.not.(ios.ne.0) .and. tmpbas%name>' ') then
                    tmpbas%file=tmpbas%name
                    tmpbas%comp='LARGE'
                    tmpbas%type='BASIS'
                    CALL ADD_BASIS(bas,firstbas,tmpbas)
                    CYCLE
                  else if (ios.ne.0) then
                     write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "Specify a filename or full basisname for ",
     &                 ".SPECIAL basis."
                     exit
                  end if
                else if(tmpbas%name>' ' .and. tmpbas%type>' ' 
     &            .and. tmpbas%comp>' ') then
                  SELECT CASE(tmpbas%type)
                  CASE ('BASIS','MOLFBAS')
                    backspace(tmp_file)               
                    read(tmp_file,*,iostat=ios) tmpbas%name,tmpbas%comp,
     &                                      tmpbas%type,tmpbas%file
                    if (.not.(ios.ne.0) .and. tmpbas%file>' ') then
                      CALL ADD_BASIS(bas,firstbas,tmpbas)
                    else if (ios.ne.0) then
                     write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "Specify a filename for ",tmpbas%type
                    end if
                  CASE('INTGRL','EXPLICIT')
!     FIRST read rest of the line and then read the explicit coefficients.
!     TMP We read the coefficients in in the old way in write commonblock
                    read(tmp_file,*,iostat=ios) iqm
                    if ( ios.ne.0 .or. iqm.eq.0) then
                      write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "No iqm for ",tmpbas%type," specified."
                      cycle
                    end if
                    ALLOCATE(JCO(iqm))
                    read(tmp_file,*,iostat=ios) iqm,JCO(1:iqm)
                    if (ios.ne.0) then
                      write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "Specify at least ",iqm," coefficients"
                      cycle
                    end if
                    CALL ADD_JCO(tmpbas,JCO)
                    CALL ADD_BASIS(bas,firstbas,tmpbas)
                    DEALLOCATE(JCO)
  
                CASE('WELLTEMP','EVENTEMP','GEOM','FAMILY','DUALFAMILY')
                    CALL ADD_FAM(tmp_file,linelen,tmpbas)
                    
                    read(tmp_file,*,iostat=ios) iqm
                    if (ios.ne.0 .or. iqm==0) then
                      write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "No iqm for ",tmpbas%type," specified."
                      cycle
                    end if
                    ALLOCATE(JCO(iqm))
                    read(tmp_file,*,iostat=ios) iqm,JCO(1:iqm)
                    if (ios.ne.0) then
                      write(tmp_file,*) "ERROR(readbasis) :", 
     &                 "Specify at least ",iqm," coefficients"
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
                  
              ELSE IF (section(2)==''.AND.                         
     &                 section(3)=='.CENTERS') THEN
!               Explicit specification of charges for non-standard nuclei, defaults are taken from nucdata.h
                BACKSPACE(tmp_file)
                name=' '
                charge=0.0
                nucexp=0.0
  
                read(tmp_file,*,iostat=ios) name, charge, nucexp
                if (ios.ne.0) then
                  read(tmp_file,*,iostat=ios) name, charge
                  if (ios.ne.0) then
                    write(tmp_file,*) "ERROR(readbasis) :", 
     &          "Not a valid name and charge specified for *CENTERS"
                    cycle
                  end if
                end if

                mol=>firstatom
                DO 
                  mol=>FIND_ATOM(mol,name)
                  IF (associated(mol)) THEN 
                    mol%charge=charge
                    mol%nucexp=nucexp
                    mol=>mol%next
                  ELSE
                    EXIT
                  END IF
                END DO

              ELSE IF (section(2)==''.AND.                         
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
                CASE DEFAULT
                  backspace (tmp_file)
                  read(tmp_file,*,iostat=ios) unitname,unit_conversion
                  if (ios.ne.0) then
                    write(lupri,*) "ERROR(readbasis) :", 
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
            WRITE(LUPRI,*) 'No basisset given in input file.            & 
     &            Specify a basisset in  **MOLECULE *BASIS .DEFAULT'
            CALL QUIT('*** ERROR SORT_BASIS, no defaultbasis given in   & 
     &           input file.')
          END IF
          
!     We now have add all data that was explicitly added in the .INP file.
!     Add default basis.
          IF(.NOT.(PARHER).AND.IPREAD>1)                                &
     &         WRITE(LUPRI,*) 'Processing basissets...'
          CALL SORT_BASIS(firstatom,firstbas,.TRUE.,lupri)
!     Sort the atoms according to the bases.
          CALL SORT_ATOMS(firstatom,firstbas,lupri)
!     Add charges if they are missing
          mol=>firstatom
          prevname=' '
          DO WHILE(ASSOCIATED(mol))
            IF (mol%charge==0) THEN
              IF(.NOT.(mol%name)==prevname)                             & 
     &              charge=FIND_CHARGE(mol%name,lupri)                
              mol%charge=charge
            END IF
  
            prevname=mol%name
            mol=>mol%next
          END DO
          
!     Convert units
          IF(.NOT.(unitname=='bohr')) THEN
            WRITE(LUPRI,*) 'Coordinates were entered in ', unitname
            WRITE(LUPRI,*) 'Converting to atomic units with',
     &                     ' conversion factor=', UNIT_CONVERSION 
            mol=>firstatom
            DO WHILE (ASSOCIATED(mol))
              mol%x=mol%x*UNIT_CONVERSION
              mol=>mol%next
            END DO
          END IF

        END SUBROUTINE READ_BASIS


        SUBROUTINE WRITE_ATOMS(lupri)
!     JHS050707 Writes the data contained in the list of atoms to the common blocks.


      use codata
      IMPLICIT NONE
!     CONSTANTS
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"

!     COMMONBLOCKS
#include "clsfmm.h"
#include "cbirea.h"

      INTEGER, INTENT(IN)             :: lupri

      INTEGER :: KATOM,natom,count
      TYPE(atom), POINTER             :: mol

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "nucdata.h"

      REAL(KIND=8), PARAMETER :: THRMIN = 1.D-15,DSM=1.0D-30

      CALL QENTER('WRITE_ATOMS')

      KATOM=count_atoms(firstatom)
      IF (KATOM>MXCENT) THEN
         WRITE(lupri,*) 'Number of atoms specified greater max number   & 
     &        of atoms'
         WRITE(lupri,*) 'Adjust maxorb.h'
         CALL QUIT('*** ERROR (READIN) natom>mxcent')
      END IF

!     nuclei.h
      NUCIND=0
      mol=>firstatom
      natom=0

      DO

        IF (.NOT.(associated(mol))) THEN
          EXIT
        ELSE
          natom=natom+1
!     write the atomic data.
!     If we didn't read the nuclear model somewhere read the nuclear radius.

!    WRITE NUCLEAR SHAPE. TMP To be changed, now always D0 is taken
          IF (GAUNUC.AND.NINT(mol%charge).NE.0) THEN
               
            IF (mol%nucexp.EQ.D0) THEN
              CALL NUCSIZ(mol%charge,mol%nucexp)
            END IF
            IF (mol%nucexp .LE. 1.0D8) THEN
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'Nuclear gaussian exponent small: ',       &
     &                         mol%nucexp
            END IF
            WRITE(lupri,'(/A,F8.3,A,1P,D20.10)')                        &
     &       '  Nuclear Gaussian exponent for atom of charge',          &
     &       mol%charge,' :',mol%nucexp
          ELSE
            mol%nucexp=D0
          END IF

!     REALLY WRITE COMMONBLOCKS
          NAMN(natom)=mol%name
          CORD(:,natom)=mol%x
          NAMEX(3*natom)=TRIM(mol%name)//' z'
          NAMEX(3*natom-1)=TRIM(mol%name)//' y'
          NAMEX(3*natom-2)=TRIM(mol%name)//' x'
          CHARGE(natom)=mol%charge
          GNUEXP(natom)=mol%nucexp

!MI     ... control printout
          IF (IPREAD.GE.4) THEN
            write(LUPRI,*) 'NAMN(natom)=',NAMN(natom)
            write(LUPRI,*) 'mol%x:',mol%x(1),mol%x(2),mol%x(3)
            write(LUPRI,*) 'CORD(:,natom)',                             &
     &      CORD(1,natom),CORD(2,natom),CORD(3,natom)
          ENDIF
          
          mol=>mol%next
        END IF
      END DO

!      IF (ADDSTO) CALL STOADD()
!
!     ***************************************
!     ***** Process symmetry input data *****
!     ***************************************
!     TMP DETECT SYM?

!      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)

!

      NUCIND=natom

      CALL QEXIT('WRITE_ATOMS')

      END SUBROUTINE WRITE_ATOMS

        SUBROUTINE POS_FILE(infile,lupri,linelen,                       &
     &     section,subsection,name,comp)
!     JHS050707 A function for later use, positions a file to a specific keyword.
        
          INTEGER, INTENT(IN)   :: infile,lupri,linelen
          CHARACTER(LEN=*), INTENT(IN) ::section,subsection,name,comp
          CHARACTER(LEN=linelen) :: str,pressect,pressub
          character(len=linelen) :: tmp_name,tmp_comp
          INTEGER                :: index

          REWIND(infile)
          DO            
            read(infile,*,iostat=ios) str
            IF (ios.ne.0) THEN
              WRITE(lupri,*) '*** ERROR KEY NOT FOUND. ***'
              CALL QUIT('*** ERROR KEY NOT FOUND. ***')
            ELSE IF (str(1:2)=='**') THEN
              pressect=str              
            ELSE IF (str(1:1)=='*') THEN
              pressub=str
            ELSE IF (pressect==section.AND.pressub==subsection) THEN
              IF (str==name) THEN
                BACKSPACE(infile)
                tmp_name=' '
                tmp_comp=' '
                read(infile,*,iostat=ios) tmp_name,tmp_comp
                IF (.not.(ios.ne.0) .and. str==comp) EXIT
              END IF
            END IF
          END DO
        END SUBROUTINE POS_FILE

        SUBROUTINE WRITE_BASIS(infile,lupri,WORK,LWORK,                 &
     &                  KANG,KSETS,KBLOCK,KPRIM,HERMIT)
!     JHS050707 Writes the contained in the basisset list to the commonblocks.
!               Note the routine is higly similar tot the subroutine READI1 for
!               .MOL files, the main difference is the split between reading the files
!               and writing them to the commonblocks

      use codata
      IMPLICIT NONE
!     CONSTANTS
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"

!     COMMONBLOCKS
#include "clsfmm.h"
#include "cbirea.h"
      REAL(KIND=8), PARAMETER :: THRMIN = 1.D-15,DSM=1.0D-30

      INTEGER,INTENT(IN)   :: lupri,infile,KANG,KSETS,KBLOCK,KPRIM
      INTEGER,INTENT(INOUT)::LWORK

      INTEGER      :: KLAST,i,j,IBSFLAG, KAOVEC, NBLOCK,                &
     &                katom,basindex,index

      REAL(KIND=8) :: QEXP=D0
      REAL         :: MXSPD
      LOGICAL      :: HERMIT,DOOWN,ANG,NOORBTS,NOFITFS,                 & 
     &                CNTBAS=.FALSE.,DIRCON(2)=.FALSE.,IFXYZ(3)
      CHARACTER(LEN=72) :: TTITLE(2)=' '
      CHARACTER(LEN=11) :: SYMTXT
      CHARACTER(LEN=15) :: CLASS
      CHARACTER(LEN=namelen) :: prevname

      TYPE(basistype),POINTER       :: baspointer
      TYPE(atom),POINTER            :: mol

      REAL(KIND=8):: WORK(LWORK) 
      REAL(KIND=8),ALLOCATABLE::                                        & 
     &               ALPHA(:,:,:),                                      &
     &               CPRIM(:,:,:,:),                                    &
     &               CPRIMU(:,:,:,:)
      INTEGER,ALLOCATABLE :: NONT(:),JCO(:,:,:),JCO2(:,:),IBLOCK(:),    &
     &                       NUC(:,:),IQM(:,:),NRC(:,:),NBLCK(:,:),     &
     &                       ISGEN(:)
      LOGICAL,ALLOCATABLE :: SEG(:,:)
      CHARACTER(LEN=80),ALLOCATABLE :: BASREF(:,:,:)

!     JHS050707 For type conversion purposes.
      INTEGER :: int8
      integer blk

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
!#include "primit.h"
!#include "shells.h"
!#include "symmet.h"
!#include "aosotr.h"
!#include "infpar.h"
!#include "dcbrkb.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
      LOGICAL LINEAR
#endif
!#include "ibtfun.h"
      
      
      katom=count_atoms(firstatom)

      ALLOCATE(NONT(KATOM))
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

      baspointer=>firstbas
      mol=>firstatom
      basindex=0

      DO
!     Loop over all basissets
        NOORBTS=.TRUE.
        NOFITFS=.TRUE.

        IF (.NOT.(associated(baspointer))) EXIT
        basindex=basindex+1
        NONT(basindex)=count_atoms(firstatom,baspointer%name)
        BASREF(1,basindex,1) = ' Reference not found in input'
        BASREF(2:10,basindex,1:3) = ' Not initialized'

!     JHS060707 TMP We could introduce some additional error checks?
        SELECT CASE (baspointer%type)
        CASE ('INTGRL','EXPLICIT')
          IBSFLAG=1
          IQM(basindex,1)=SIZE(baspointer%jco)
          DO i=1,IQM(basindex,1)
            JCO(i,basindex,1)=baspointer%jco(i)
!MI         write(lupri,*) 'JCO(i,basindex,1)=',JCO(i,basindex,1)
!MI         call flshfo(lupri)
          END DO           
          BASNAM=baspointer%file
          KAOVEC = KBLOCK + 1 - IBLOCK(1)

!     POSITION the input file.
          CALL POS_FILE(infile, lupri, linelen,                         &
     &         '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
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
          mol=>find_atom(mol,baspointer%name)
          CALL BASLIB(IQM(basindex,1),JCO(1,basindex,1),                &
     &             NUC(IBLOCK(1),1),                                    &
     &             NRC(IBLOCK(1),1),SEG(IBLOCK(1),1),                   &
     &             ALPHA(1,IBLOCK(1),1),                                &
     &             CPRIM(1,1,IBLOCK(1),1),CPRIMU(1,1,IBLOCK(1),1),      &
     &             NBLOCK,ISGEN(IBLOCK(1)),                             &
     &             KAOVEC,KPRIM,mol%charge,                             &
     &             1,DSM,BASREF(1,basindex,1))                         
          NBLCK(basindex,1) = NBLOCK
        CASE ('WELLTEMP','EVENTEMP','GEOM')
          IBSFLAG=4
          CALL POS_FILE(infile, lupri, linelen,
     &         '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
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
          CALL POS_FILE(infile, lupri, linelen,
     &         '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
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

!     PRINT BASISSET
        NOORBTS=NOORBTS.AND.IQM(basindex,1).EQ.0
        IF (IQM(basindex,1).GT.0) THEN
          WRITE(lupri,'(3X,A)') 'Large component set:'
          WRITE(lupri,'(6X,A,I5,A,12I5)')                               &
     &    'Max.ang.quantum no.:',(IQM(basindex,1)-1),                   & 
     &    '  Blocks:',(JCO(i,basindex,1), i=1,  IQM(basindex,1))
          NHTYP  =   MAX(NHTYP,IQM(basindex,1))
        ENDIF

        DO i = 1, NBLCK(basindex,1)
           IF ((NUC(i,1).NE.NRC(i,1)).AND.(NRC(i,1).GT.0))              &
     &          CNTBAS = .TRUE.
        END DO
        IF ( CNTBAS .AND. (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
           DIRCON(1) = .TRUE.
        END IF
          
        IF(.NOT.TWOCOMP) THEN
          IF (ASSOCIATED(baspointer%next)) THEN
            IF (baspointer%next%name==baspointer%name)                    &
     &            baspointer=>baspointer%next
!     JHS060707 Move to small component.
!               TMP We can extend it to include a search for the small component
          END IF

          IF (.NOT.(baspointer%comp(1:5)=='SMALL')) THEN 
            WRITE(lupri,*) 'Small component not specified.'
            WRITE(lupri,*) 'Assuming kinetic balance.'
            IBSFLAG=4
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
!     POSITION the input file.
            CALL POS_FILE(infile, lupri, linelen,
     &           '**BASIS', '*BASIS', baspointer%name,baspointer%comp) 
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

          IF (IQM(basindex,2).GT.0) THEN
            WRITE(lupri,'(3X,A)') 'Small component set:'
            WRITE(lupri,'(6X,A,I5,A,12I5)')                               &
     &        'Max.ang.quantum no.:',(IQM(basindex,2)-1),                 &
     &        '  Blocks:',(JCO(i,basindex,2), i=1,IQM(basindex,2))
          ENDIF

          NHTYP = MAX(NHTYP,IQM(basindex,2))
          
          IF (IQM(basindex,2).GT.0) THEN
            IF(INPTST) THEN
               WRITE(lupri,'(3X,A)') 'Small component set:'
               WRITE(lupri,'(6X,A,I5,A,12I5)')                            &
     &              'Max.ang.quantum no.:',(IQM(basindex,2)-1),           &
     &              '  Blocks:',(JCO(i,basindex,2), i=1,IQM(basindex,2))
            ENDIF
            NHTYP     = MAX(NHTYP,IQM(basindex,2))
!            IF(NHTYP.GT.MXQN) GOTO 5000
          ENDIF
          
          DO i = 1, NBLCK(basindex,2)
            IF ((NUC(i,2).NE.NRC(i,2)).AND.(NRC(i,2).GT.0))               &
     &             CNTBAS = .TRUE.
          END DO
          IF ( CNTBAS .AND. IBSFLAG .NE. 4 ) DIRCON(2) = .TRUE.

        END IF

        IBLOCK(1:2)=IBLOCK(1:2)+NBLCK(basindex,1:2)
        NOORBT(basindex+1:NONT(basindex)+basindex) = NOORBTS
          
        DO i=3,4
!     JHS060707 TMP Include also the reading of fitsets.
          IBSFLAG=2
          IBLOCK(i)=IBLOCK(i)+NBLCK(basindex,i)
            
        END DO
        prevname=baspointer%name
        DO WHILE (associated(baspointer))
          IF (.NOT.(baspointer%name==prevname)) EXIT
          baspointer=>baspointer%next
        END DO
        
      END DO

      LINEAR = .FALSE.
      CLASS = 'N/A'
      IF (ADDSYM) THEN
         CALL SYMADD(WORK,LWORK,basindex,NSYMOP,NONT,KATOM,KASYM,
     &               CLASS,.false.)
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
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)

!
!     ***************************************************
!     ***** Process orbital and geometry input data *****
!     ****************************************************
!
      CALL BASPRO(WORK,LWORK,basindex,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,NRC,&
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,KBLOCK,KPRIM, &
     &            DOOWN, blk)
!
!     *************************************************
!     ***** Print orbital and geometry input data *****
!     *************************************************
!

      CALL BASOUT(WORK,LWORK,basindex,NONT,IQM,NBLCK,JCO,BASREF,        &
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,      &
     &            KBLOCK,KPRIM,HERMIT, blk)
!
!     TMP copy input to .out file.

!     Default threshold for distinguishing classical and
!     non-classical contributions to Coulomb integrals

      THRCLS = THRS
!
!     ****************************************
!     ***** Output on LUONEL *****************
!     ****************************************
!

      IF (HERMIT) CALL WRONEL(TTITLE,basindex,NONT,IQM,IFXYZ,KATOM,JCO2,&
     &                        KANG)

          DEALLOCATE(NONT)
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

        END SUBROUTINE WRITE_BASIS

        SUBROUTINE PRINT_ATOMS(firstatom, lupri)
          !TYPE(atom), POINTER, INTENT(IN) :: firstatom
          TYPE(atom), POINTER :: firstatom
          INTEGER, INTENT(IN) :: lupri
          TYPE(atom), POINTER :: mol

          CALL TITLER('XYZ coordinates read','*',118)
          WRITE(lupri,'(a,t13,a,t25,a,t37,a,t49,a,t61,a)') 'Label',     &
     &         'x','y','z','charge','nuc exp'
          mol=>firstatom
          DO
            IF (associated(mol)) THEN
              WRITE(lupri,'(a,t9, 5f12.6)')                             & 
     &          TRIM(mol%name), mol%x, mol%charge, mol%nucexp
              mol=>mol%next
            ELSE
              EXIT
            ENDIF
          END DO
        WRITE(lupri,'(a)') ' '

        END SUBROUTINE PRINT_ATOMS

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

        SUBROUTINE DEL_ATOMS(firstatom, lupri)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          !TYPE(atom),POINTER,INTENT(INOUT) :: firstatom
          TYPE(atom),POINTER :: firstatom
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
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
              exit
            ENDIF
          END DO

        END SUBROUTINE DEL_ATOMS

        SUBROUTINE DEL_BASIS(firstbas, lupri)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          !TYPE(basistype),POINTER,INTENT(INOUT) :: firstbas
          TYPE(basistype),POINTER :: firstbas
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
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

        SUBROUTINE CLEAR_XYZ(lupri)
!     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          CALL DEL_ATOMS(firstatom, lupri)
          CALL DEL_BASIS(firstbas, lupri)

        END SUBROUTINE CLEAR_XYZ


        SUBROUTINE READ_XYZ(molfile,inpfile, lupri)
#include "infpar.h"        
!     THE MAIN SUBROUTINE THAT READS AN XYZ FILE

          INTEGER, INTENT(IN) :: molfile,inpfile,lupri

          CALL READ_ATOMS(molfile,linelen,firstatom,lupri)
          CALL READ_BASIS(inpfile, linelen,lupri)

          IF(.NOT.(PARHER)) CALL PRINT_ATOMS(firstatom,lupri)
          IF(.NOT.(PARHER)) CALL PRINT_BASIS(firstbas,lupri)
            


        END SUBROUTINE READ_XYZ


      END MODULE READFILE

