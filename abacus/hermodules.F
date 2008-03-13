
C     !!! NEW MODULE FOR CONVERTING BETWEEN TYPES !!!
C     OTHER SOLLUTIONS (C++) ARE MORE THAN WELCOME
C     JHS, 11-04-2007
      MODULE CONVERSION
        PRIVATE
        PUBLIC INTTOCHAR
C     Maybe include priunit.h instead of passing output variable?
        INTEGER, SAVE      :: zerochar=-1
        INTEGER, PARAMETER :: maxint=11
        
      CONTAINS
        CHARACTER(LEN=maxint) FUNCTION INTTOCHAR(n)
#include "priunit.h"
          INTEGER, INTENT(IN) :: n
          INTEGER             :: i,j,n_temp

CMI        write(lupri,*) 'INTTOCHAR: n=',n
CMI ...  NASTY BUG with gfortran/g95 - spoiled integer !!!
CMI  ... needed help_variable
           n_temp = n

          IF (zerochar<0) zerochar=IACHAR('0')
          inttochar=' '
          DO i=0,maxint-1
            IF (10**i>n) EXIT
          END DO 
          IF (i==maxint) i=i-1
          DO j=i,0, -1
CMI         write(lupri,*) 'n,j,zerochar=',n,j,zerochar
            inttochar(i-j:i-j)=ACHAR(MOD((n_temp/(10**j)), 10)+zerochar)
CMI         write(lupri,*) 'n,j,zerochar=',n,j,zerochar
CMI         write(lupri,*) 'inttochar(i-j:i-j)=',inttochar(i-j:i-j)
          END DO
        END FUNCTION INTTOCHAR

      END MODULE CONVERSION

C     !!! MODULE FOR HANDLING FILES !!!
C     Contains
C       - OPEN_FILE(filename) or OPEN_FILE(filename, filenum,formatted,output)
C         Opens a file and returns the id. Which is a new number if only the
C         name is given, otherwise it is equal to filenum.
C         A negative or zero indicates an error.
C       - GET_LENGTH(filename) or GET_LENGTH(filenum)
C         Returns the maximum linelength of a file.
C       - GET_EXTENSION(filename)
C         Returns the first extension in the array filetypes for which
C         filename.extension is an existing filename.
C       - READ_VAR(filenum,linewidth,var,index,fixedwidth,skip)
C         Reads a variable of character, integer or real type starting
C         from position index.
C         By default text is read from index tot the next unquoted space
C         or one of the unqouted commentsigns: '#'
C         if fixedwidth is given all characters between index and index+
C         fixedwidth are read into var.
C         The type of variable var determines how the text will be interpreted.
C         If skip='TRUE', the filepointer is positioned to the next line.
C         Index is set to the first non read character.
C         If it is 0 the end of the file is reached.
C         If it smaler than 0 an error has occured.
C     For the errorcodes see the comments in the corresponding routines.
C    
C     JHS 05-07-07

      MODULE FILE
      IMPLICIT NONE

      PRIVATE
      PUBLIC OPEN_FILE, GET_LENGTH, GET_EXTENSION, READ_VAR

        CHARACTER(LEN=8), DIMENSION(2), PARAMETER :: 
     &       filetypes=(/'MOL', 'XYZ'/)
        INTEGER :: linewidth=0

        INTERFACE OPEN_FILE
          MODULE PROCEDURE OPEN_FILE_NAME, OPEN_FILE_NUM
        END INTERFACE
        
        INTERFACE READ_VAR
#if defined (SYS_SX64) 
CSK       SK: 16-11-2007 on NEC's 64-Bit SX-6/SX-8 machines
CSK       we do not have REAL! REAL*8 is the minimum...
          MODULE PROCEDURE READ_CHAR,READ_INT,READ_REAL8
#else
          MODULE PROCEDURE READ_CHAR,READ_INT,READ_REAL,READ_REAL8
#endif
        END INTERFACE  
        
        INTERFACE GET_LENGTH
          MODULE PROCEDURE GET_LENGTH_NUM, GET_LENGTH_NAME
        END INTERFACE


      CONTAINS

        INTEGER FUNCTION OPEN_FILE_NAME(filename)
C     Opens a file and returns it's identifying number (>0) if the operation was succesful.
C     If the opening was insuccesful a negative or 0 is returned. 
C      0  -   File does not exists.
C     -1  -   To many files opened.
C     -99  -   Error while executing the open statement.
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

        INTEGER FUNCTION OPEN_FILE_NUM(
     &     filename, filenum, formatted, lupri)
C     Opens a file and returns it's identifying number (>0) if the operation was succesful.
C     If the opening was insuccesful a negative or 0 is returned. 
C      0   -   File does not exists.
C     -1   -   To many files opened.
C     -2   -   Number already assigned
C     -3   -   File already opened with a different number.
C     -99  -   Error while executing the open statement.
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

        INTEGER FUNCTION GET_LENGTH_NUM(filenum, lupri)
C     Returns the maximum linelength of an already opened file.
C
          INTEGER, INTENT(IN) :: filenum
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          INTEGER             :: ios, length, newlength
          CHARACTER(LEN=1)   :: tmpchr

          REWIND(filenum)
          length=0
          DO
            newlength=0
            DO
              READ(filenum, '(a)', iostat=ios, ADVANCE='no') tmpchr
              IF (ios<0) THEN
                length=MAX(length, newlength)
                EXIT
              ELSE
                newlength=newlength+1
              END IF
            END DO
            IF (ios==-1) THEN
              EXIT
            ELSE IF (.NOT.(ios==-2).and.PRESENT(lupri)) THEN
              WRITE (lupri,'(A)') "WARNING problems encountered 
     &                       while reading file."
            END IF
          END DO
          REWIND(filenum)
        get_length_num=length
        END FUNCTION GET_LENGTH_NUM
      
        INTEGER FUNCTION GET_LENGTH_NAME(filename, lupri)
C     Opens a file, returns the maximum line length and closes the file again
C
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          INTEGER :: filenum
          
          filenum=open_file(filename)
          IF (filenum>0) THEN
            get_length_name=get_length_num(filenum)
            close(filenum)
          ELSE
            WRITE(lupri,'(a)') 'Could not open file. name=', filename
            get_length_name=-1
          END IF
        END FUNCTION GET_LENGTH_NAME

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
               
        SUBROUTINE READ_CHAR(filenum,linewidth,var,index,fixedwidth,
     &                       skip)
C     TMP adjust for an optional linewith.

          INTEGER, INTENT(IN)           :: filenum
          INTEGER, INTENT(IN)           :: linewidth
          INTEGER, INTENT(INOUT)        :: index
          CHARACTER(LEN=*), INTENT(OUT) :: var
          INTEGER, INTENT(IN), OPTIONAL :: fixedwidth
          LOGICAL, INTENT(IN), OPTIONAL :: skip

          LOGICAL                       :: quoted=.FALSE.
          LOGICAL                       :: lskip=.FALSE.
          INTEGER :: iprev, ios, i
          CHARACTER(LEN=linewidth) :: str
          
          if (present(skip)) lskip = skip
          var=' '; iprev=index

          READ(filenum, '(a)', ADVANCE='yes', iostat=ios) str

          IF (ios==0) THEN
            IF (present(fixedwidth)) THEN
              DO i=index,MIN(len(str), index+fixedwidth)
                SELECT CASE (str(i:i))
                CASE ('#')
                  index=1
                  EXIT                  
                END SELECT 
              END DO
              IF (i-2>index.OR.fixedwidth==1) 
     &          var=str(iprev:i-1)
              IF (lskip) THEN
                index=1
              ELSE IF (i-1>linewidth) THEN
                index=1
              ELSE
                BACKSPACE(filenum)
                index=i
              END IF
            ELSE  
              DO i=index, len(str)
                SELECT CASE (str(i:i))
                CASE ('#')
                  IF (.NOT.(quoted)) THEN
                    IF (iprev<i-1) 
     &                var=str(iprev:i)
                    index=1
                    RETURN
                  END IF
                CASE (' ')
                  IF (.not.(quoted)) THEN
                    IF (iprev==1.AND.iprev<i) THEN
                      var=str(iprev:i)
                      IF (lskip) THEN
                        index=1
                      ELSE
                        BACKSPACE(filenum)
                        index=i
                      END IF
                      RETURN
                    ELSE IF (iprev<i-1) THEN
                      var=str(iprev+1:i)
                      IF (lskip) THEN
                        index=1
                      ELSE
                        BACKSPACE(filenum)
                        index=i
                      END IF
                      RETURN
                    END IF
                    iprev=i
                  END IF
                CASE ('"')
                  quoted=.not.(quoted)
                  if (.not.(quoted)) THEN
                    var=str(iprev+1:i-1)
                    IF (lskip) THEN
                      index=1
                    ELSE
                      BACKSPACE(filenum)
                      index=i+1
                    END IF
                    RETURN
                  END IF                    
                END SELECT 
              END DO             

              var=str(iprev:i)
              index=1
            END IF
          ELSE
            index=-1
          END IF
        END SUBROUTINE READ_CHAR

        REAL(KIND=8) FUNCTION READ_NUM8(filenum,iprev,i,ios,idot)
          USE CONVERSION
#include "priunit.h"
          INTEGER, INTENT(IN)    :: filenum, i, iprev
          INTEGER, INTENT(IN), OPTIONAL :: idot
          INTEGER, INTENT(OUT) :: ios
          CHARACTER(LEN=64) :: chr
         
          CALL QENTER('READ_NUM8') 
   
CMI       write(lupri,*) 'idot,iprev,i=',idot,iprev,i

          IF (present(idot) .AND. idot>0) THEN
            chr='(t'//trim(inttochar(iprev+1))//', f'//
     &           trim(inttochar(i-iprev-1))//'.'//
     &           trim(inttochar(i-idot-1))//')'
          ELSE
            chr='(t'//trim(inttochar(iprev+1))//', f'//
     &           trim(inttochar(i-iprev-1))//')'
          END IF
          BACKSPACE(filenum)

CMI        write(lupri,*) 'chr=',chr

          READ(filenum, chr,ADVANCE='yes', iostat=ios) READ_NUM8

CMI        write(lupri,*) 'READ_NUM8=',READ_NUM8

          CALL QEXIT('READ_NUM8') 

        END FUNCTION READ_NUM8

        REAL FUNCTION READ_NUM(filenum,iprev,i,ios,idot)
          USE CONVERSION
          INTEGER, INTENT(IN)    :: filenum, i, iprev
          INTEGER, INTENT(IN), OPTIONAL :: idot
          INTEGER, INTENT(OUT) :: ios
          CHARACTER(LEN=64) :: chr
         
          IF (present(idot) .AND. idot>0) THEN
            chr='(t'//trim(inttochar(iprev+1))//', f'//
     &           trim(inttochar(i-iprev-1))//'.'//
     &           trim(inttochar(i-idot-1))//')'
          ELSE
            chr='(t'//trim(inttochar(iprev+1))//', f'//
     &           trim(inttochar(i-iprev-1))//')'
          END IF
          BACKSPACE(filenum)
          READ(filenum, chr,ADVANCE='yes', iostat=ios) READ_NUM
        END FUNCTION READ_NUM

      INTEGER FUNCTION READ_INTEGER(filenum,iprev,i,ios)
          USE CONVERSION
          INTEGER, INTENT(IN)    :: filenum, i, iprev
          INTEGER, INTENT(OUT) :: ios
          CHARACTER(LEN=64) :: chr
         
          chr='(t'//trim(inttochar(iprev+1))//', i'//
     &         trim(inttochar(i-iprev-1))//')'

          BACKSPACE(filenum)
          READ(filenum, chr,ADVANCE='yes', iostat=ios) READ_INTEGER
        END FUNCTION READ_INTEGER

        SUBROUTINE READ_REAL8(filenum,linewidth,var,index,fixedwidth,
     &                        skip)
C     TMP adjust for an optional linewith.
#include "priunit.h"
          INTEGER, INTENT(IN)           :: filenum
          INTEGER, INTENT(IN)           :: linewidth
          INTEGER, INTENT(INOUT)        :: index
          REAL(KIND=8), INTENT(OUT) :: var
          INTEGER, INTENT(IN), OPTIONAL :: fixedwidth
          LOGICAL, INTENT(IN), OPTIONAL :: skip
          
          LOGICAL                       :: lskip=.FALSE.
          INTEGER :: iprev, ios, idot=0,i
          CHARACTER(LEN=linewidth) :: str
          
          if (present(skip)) lskip = skip
          var=0; iprev=index; idot=0

          READ(filenum, '(a)', ADVANCE='yes', iostat=ios) str

          IF (ios==0) THEN
            IF (present(fixedwidth)) THEN
              DO i=index+1,MIN(len(str), index+fixedwidth)
                SELECT CASE (str(i:i))
                CASE ('#')
                  EXIT
                CASE ('.')
                  idot=i
                END SELECT               
              END DO
              IF (i-2>index.OR.fixedwidth==1) 
     &          var=READ_NUM8(filenum,index,i-1,ios,idot)
              IF (lskip) THEN 
                index=1
              ELSE IF (i>linewidth) THEN 
                index=1
              ELSE
                BACKSPACE(filenum)
                index=i
              END IF
            ELSE
              DO i=index, len(str)
                SELECT CASE (str(i:i))
                CASE ('#')
                  IF (iprev<i-1) 
     &              var=READ_NUM8(filenum,iprev,i,ios,idot)
                  index=1
                  RETURN
                CASE (' ')
                  IF (iprev<i-1.OR.(i==2.AND..NOT.(index==2))) THEN
                    var=READ_NUM8(filenum,iprev,i,ios,idot)
                    IF (lskip) THEN 
                      index=1
                    ELSE
                      BACKSPACE(filenum)
                    index=i
                    END IF
                    RETURN
                  END IF
                  iprev=i
                CASE ('.')
                  idot=i
                END SELECT 
              END DO
              var=READ_NUM8(filenum,iprev,i,ios,idot)
              index=1
            END IF
          ELSE
            index=-1
           END IF
 
CMI       write(LUPRI,*) 'var=',var

        END SUBROUTINE READ_REAL8

#if !defined (SYS_SX64)
        SUBROUTINE READ_REAL(filenum,linewidth,var,index,fixedwidth,
     &                       skip)
C     TMP adjust for an optional linewith.
#include "priunit.h"
          INTEGER, INTENT(IN)           :: filenum
          INTEGER, INTENT(IN)           :: linewidth
          INTEGER, INTENT(INOUT)        :: index
          REAL, INTENT(OUT) :: var
          INTEGER, INTENT(IN), OPTIONAL :: fixedwidth
          LOGICAL, INTENT(IN), OPTIONAL :: skip
          
          LOGICAL                       :: lskip=.FALSE.
          INTEGER :: iprev, ios, idot=0,i
          CHARACTER(LEN=linewidth) :: str
          if (present(skip)) lskip = skip
          var=0.0; iprev=index; idot=0
          READ(filenum, '(a)', ADVANCE='yes', iostat=ios) str


          IF (ios==0) THEN
            IF (present(fixedwidth)) THEN
              DO i=index+1,MIN(len(str), index+fixedwidth)
                SELECT CASE (str(i:i))
                CASE ('#')
                  EXIT
                CASE ('.')
                  idot=i
                END SELECT  
              END DO
              IF (i-2>index.OR.fixedwidth==1)
     &          var=READ_NUM(filenum,index,i-1,ios,idot)
              IF (lskip) THEN 
                index=1
              ELSE IF (i>linewidth) THEN 
                index=1
              ELSE
                BACKSPACE(filenum)
                index=i
              END IF
            ELSE
              DO i=index+1, len(str)
                SELECT CASE (str(i:i))
                CASE ('#')
                  IF (iprev<i-1.OR.i==2) 
     &              var=READ_NUM(filenum,iprev,i,ios,idot)
                  index=1
                  RETURN
                CASE (' ')
                  IF (iprev<i-1.OR.(i==2.AND..NOT.(index==2))) THEN
                    var=READ_NUM(filenum,iprev,i,ios,idot)
                    IF (lskip) THEN 
                      index=1
                    ELSE
                      BACKSPACE(filenum)
                    index=i
                    END IF
                    RETURN
                  END IF
                  iprev=i
                CASE ('.')
                  idot=i
                END SELECT 
              END DO
              var=READ_NUM(filenum,iprev,i,ios,idot)
              index=1
            END IF
          ELSE
            index=-1
           END IF

CMI       write(LUPRI,*) 'var=',var

        END SUBROUTINE READ_REAL

#endif /* !defined (SYS_SX64) */

        SUBROUTINE READ_INT(filenum,linewidth,var,index,fixedwidth,
     &                      skip)
C     TMP adjust for an optional linewith.
          INTEGER, INTENT(IN)           :: filenum
          INTEGER, INTENT(IN)           :: linewidth
          INTEGER, INTENT(INOUT)        :: index
          INTEGER, INTENT(OUT) :: var
          INTEGER, INTENT(IN), OPTIONAL :: fixedwidth
          LOGICAL, INTENT(IN), OPTIONAL :: skip
          
          LOGICAL                       :: lskip=.FALSE.
          INTEGER :: iprev, ios,i
          CHARACTER(LEN=linewidth) :: str
          
          if (present(skip)) lskip = skip
          var=0; iprev=index
          READ(filenum, '(a)', ADVANCE='yes', iostat=ios) str

          IF (ios==0) THEN
            IF (present(fixedwidth)) THEN
              DO i=index+1,MIN(len(str), index+fixedwidth)
                SELECT CASE (str(i:i))
                CASE ('#')
                  EXIT
                END SELECT  
              END DO
              IF (i-2>index.OR.fixedwidth==1)
     &          var=READ_INTEGER(filenum,index,i-1,ios)
              IF (lskip) THEN 
                index=1
              ELSE IF (i>linewidth) THEN 
                index=1
              ELSE
                BACKSPACE(filenum)
                index=i
              END IF
            ELSE
              DO i=index, len(str)
                SELECT CASE (str(i:i))
                CASE ('#')
                  IF (iprev<i-1) 
     &              var=READ_INTEGER(filenum,iprev,i,ios)
                  index=1
                  RETURN
                CASE (' ')
                  IF (iprev<i-1.OR.(i==2.AND..NOT.(index==2))) THEN
                    var=READ_INTEGER(filenum,iprev,i,ios)
                    IF (lskip) THEN 
                      index=1
                    ELSE
                      BACKSPACE(filenum)
                    index=i
                    END IF
                    RETURN
                  END IF
                  iprev=i
                END SELECT 
              END DO
              var=READ_INTEGER(filenum,iprev,i,ios)
              index=1
            END IF
          ELSE
            index=-1
           END IF
        END SUBROUTINE READ_INT

      END MODULE FILE


C     !!! MODULE FOR READING XYZ FILES !!!
C     Module for reading xyz files. 
C     The coordinates are given in the .xyz file.
C     The basisset and possible charges and nuclear exponents are 
C     specified in the .inp file
C
C     In the long run the readingsection for the .inp file move to dirac/dirrdn.F
C     JHS 05-07-07

      MODULE READFILE
        USE FILE
        IMPLICIT NONE

        PRIVATE
        PUBLIC READ_XYZ,WRITE_ATOMS,WRITE_BASIS,CLEAR_XYZ
C     DETERMINES HOW DIFFERENT MOLECULE FILES SHOULD BE READ.
C     !!! PAM CAPITALIZES THE EXTENSIONS !!!

        INTEGER, PARAMETER :: namelen=8, basislen=32
        REAL(KIND=8), PARAMETER :: D0 = 0.0D0
        INTEGER            :: i

        TYPE atom
C     Defines a list for storing atomic data, like coordinates, charges and nuclear exponents.
          CHARACTER(LEN=namelen)  :: name=' '
          REAL(KIND=8)            :: x(3)=(/(0.0,i=1, 3)/)
          REAL(KIND=8)            :: charge=0.0
          REAL(KIND=8)            :: nucexp=0.0
          TYPE(atom), POINTER     :: next=>NULL()
        END TYPE atom
        TYPE(atom), POINTER :: firstatom

        TYPE basistype
C     Defines a list for storing the information about the basis sets.
          CHARACTER(LEN=namelen)     :: name=' '
          CHARACTER(LEN=basislen)    :: type=' '
          CHARACTER(LEN=basislen)    :: comp=' '
          CHARACTER(LEN=basislen)    :: file=' '
C         The following can be assigned for explicitely given basissets or Family type.
          INTEGER, POINTER           :: JCO(:)
          REAL, POINTER              :: FAM(:)
          TYPE(basistype), POINTER   :: next=>NULL()
        END TYPE basistype
        TYPE(basistype), POINTER :: firstbas

C     Global variables
C     Symmetry of the molecule
        LOGICAL :: ADDSYM=.TRUE.
        INTEGER :: NSYMOP=0
        CHARACTER(LEN=1) :: KASYM(3,3)=' '



      CONTAINS


        SUBROUTINE ADD_ATOM(prevmol, fatom, name,x)
C     JHS050707 Adds a new node to the atom list.
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
C     JHS050707 Adds a new node to the basislist
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
C     JHS050707 Assigns an array with JCO components.
          TYPE(basistype)    :: bas
          INTEGER            :: JCO(:)
          INTEGER,TARGET,ALLOCATABLE,SAVE :: JCO2(:)
          
          ALLOCATE(JCO2(SIZE(JCO)))
          JCO2=JCO
          bas%jco=>JCO2
        END SUBROUTINE

        SUBROUTINE ADD_FAM(filenum,linelen,bas,index)
C     JHS050707 Assigns for different types of family basissets an array FAM. 
C               And reads the needed coefficients from file.
          INTEGER,INTENT(IN) :: filenum,linelen
          TYPE(basistype)    :: bas
          REAL,TARGET,ALLOCATABLE,SAVE :: FAM(:)
          INTEGER,INTENT(INOUT) :: index
          INTEGER               :: i

          SELECT CASE(bas%type)
          CASE('WELLTEMP')
            ALLOCATE(FAM(5))
            DO i=1,5
              CALL READ_VAR(filenum,linelen,FAM(i),index)
            END DO
            bas%fam=FAM
          CASE('EVENTEMP')
            ALLOCATE(FAM(5))
            DO i=1,3
              CALL READ_VAR(filenum,linelen,FAM(i),index)
            END DO
            bas%fam=FAM
          CASE('GEOM')
            ALLOCATE(FAM(5))
            CALL READ_VAR(filenum,linelen,FAM(i),index)
            CALL READ_VAR(filenum,linelen,FAM(i),index)
            FAM(3)=D0
            FAM(4)=D0
            CALL READ_VAR(filenum,linelen,FAM(i),index)

            bas%fam=FAM
          CASE('FAMILY')
            ALLOCATE(FAM(1))
            CALL READ_VAR(filenum,linelen,FAM(1),index)
            bas%fam=FAM
          CASE('DUALFAMILY')
            ALLOCATE(FAM(2))
            DO i=1,2
              CALL READ_VAR(filenum,linelen,FAM(i),index)
            END DO
            bas%fam=FAM
          END SELECT
        END SUBROUTINE

        SUBROUTINE READ_ATOMS(filenum, linelen,fatom,lupri)
C     JHS050707 Reads atom name and coordinates from the .XYZ file 
C               and writes it to the atoms list.
#include "infpar.h"
#include "cbirea.h"
          TYPE(atom), POINTER :: fatom
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
          TYPE(atom), POINTER                :: mol=>NULL()
          CHARACTER(LEN=namelen)             :: name
          REAL(KIND=8)                       :: x(3)
          INTEGER                            :: i,index
 
          CALL QENTER('READ_ATOMS')
          
          IF (.NOT.(PARHER).AND.IPREAD>1) 
     &         WRITE(LUPRI,*) 'Reading coordinates from xyz file...'

          REWIND(filenum)
          index=1
          DO
            name=' '
            x=0
            IF (index<0) EXIT
            CALL READ_VAR(filenum, linelen,name, index) 
            IF (index<=1) CYCLE
            CALL READ_VAR(filenum, linelen,x(1), index)
            IF (index<=1) CYCLE
            CALL READ_VAR(filenum, linelen,x(2), index)
            IF (index<=1) CYCLE
            CALL READ_VAR(filenum, linelen,x(3), index)

            IF (IPREAD >= 4) THEN
              write(lupri,*) 'READ_ATOMS: read coordinates:',
     &        x(1),x(2),x(3)
            ENDIF
C     JHS050707 If coordinates are missing the atom is not added to the list.
            CALL ADD_ATOM(mol,fatom,name,x)            
          END DO

          IF(.NOT.(PARHER).AND. IPREAD>1) WRITE(LUPRI,*) 'Done'

          CALL QEXIT('READ_ATOMS')

        END SUBROUTINE READ_ATOMS

      
        FUNCTION find_atom(fatom, name)
C     JHS050707 Returns a pointer to the next atom in the list with name 
C               or NULL() if is no such atom. 
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
C     JHS050707 Returns a pointer to the next basisset with name in the list with basissets
C               or NULL() if there is no such basis.
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
C     JHS050707 Counts the length of the atomlist.
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

        INTEGER FUNCTION find_char(key, chararray)
C     JHS050707 Returns the first index of key in chararray, or else len(chararray)+1
          CHARACTER(LEN=*), INTENT(IN) :: key
          CHARACTER(LEN=*), INTENT(IN) :: chararray(:)
        
          DO i=1,len(chararray)
            IF (chararray(i)==key) find_char=i
            EXIT
          END DO
        END FUNCTION find_char

        SUBROUTINE SORT_ATOMS(fatom,fbas,lupri)
C     JHS050707 Sorts the atomlist in the same order as the basissets 
C               delets atoms that don't have a associated basisset
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
C     Look wether we have already ordered this species of atoms
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
C     JHS050707 Sorts the basisinput in the same order as the atoms 
C               and issues warnings if no basis can be found
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
          
C     FIND DEFAULT BASIS.
          default=>find_basis(fbas,'DEFAULT')
          IF (PRESENT(usecomp).AND.usecomp) THEN
            DO WHILE(ASSOCIATED(mol))
              first=.TRUE.
              bas=>fbas
              DO WHILE(ASSOCIATED(bas))
                bas=>find_basis(bas, mol%name)
                IF (first.AND..NOT.(ASSOCIATED(bas))) THEN
                  IF (.NOT.(ASSOCIATED(default)))
     &              CALL QUIT('*** ERROR SORT_BASIS, no defaultbasis 
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
                  CALL QUIT('*** ERROR SORT_BASIS, basis 
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
C     JHS050707 For a given atomname the corresponding charge is searched for in nucdata.h
          CHARACTER(LEN=*),INTENT(IN) :: atomname
          INTEGER,OPTIONAL,INTENT(IN) :: LUPRI
#include "nucdata.h"
          INTEGER :: i
          
          find_charge=0.0
          DO i=1,SIZE(NUCLABEL)
            IF (atomname==NUCLABEL(i)) THEN
              find_charge=REAL(i)
              RETURN
            END IF
          END DO
          IF (PRESENT(LUPRI)) WRITE(LUPRI,*) 'WARING: 
     &         No charge found for element', atomname

        END FUNCTION


        SUBROUTINE READ_BASIS(filenum,linelen,lupri)
#include "codata.h"
#include "infpar.h"
#include "cbirea.h"
C     JHS050707 Reads basis names and possibly the coefficients typed in the DIRAC.INP file.
C               Reads also Charges, nuclear exponents and the units of the XYZ file.
C     Stores it in the the basis list.
C     !!! TMP: TESTED FOR: 
C                         - 'LARGE BASIS'
C     (Implemented also for the rest)
          INTEGER, INTENT(IN)                :: filenum, linelen,lupri
          TYPE(atom), POINTER                :: mol=>NULL()
          TYPE(atom)                         :: tmpmol
          TYPE(basistype),POINTER            :: bas=>NULL()
          TYPE(basistype)                    :: tmpbas

C     FOR UNIT CONVERSION
          REAL(KIND=8)      :: UNITS=1.0
          CHARACTER(LEN=32) :: UNITNAME='a.u.', tmpunit=' '
        
          CHARACTER(LEN=linelen)             :: str,section(3)
          INTEGER                            :: i,j,index,iqm
          CHARACTER(LEN=namelen)             :: name, prevname
          REAL(KIND=8)                       :: charge, nucexp
          LOGICAL                            :: CONT=.FALSE.
          CHARACTER(LEN=9) :: symtxt

C     FOR explicitely typed basissets.
          INTEGER,TARGET,ALLOCATABLE         :: JCO(:)

          IF (.NOT.(PARHER).AND.IPREAD>1) 
     &      WRITE(LUPRI,*) 'Reading basisset information from file...'

C     initialization for security
          REWIND(filenum)

          index=1
          section=' '

          NSYMOP=0
          ADDSYM=.TRUE.
          KASYM(:,:)=' '
          WRITE(LUPRI,*) 'INIT',NSYMOP,ADDSYM,KASYM

C     Standard for xyz is angstrom.
          UNITS=1/XTANG


          DO
            CALL READ_VAR(filenum,linelen,str,index, skip=.TRUE.)          

            IF (index<0) THEN
              EXIT
            ELSE IF (str(1:2)=='**') THEN
              section(1)=str              
            ELSE IF (str(1:1)=='*') THEN
              IF (str=='*END') section=' '
              section(2)=str
            ELSE IF (str(1:1)=='.') THEN
               section(3)=str
               SELECT CASE(section(2))
               CASE('*SYMMETRY')
                  
                  IF (section(3)=='.NOSYM') ADDSYM=.FALSE.
               END SELECT

            ELSE IF(section(1)=='**MOLECULE') THEN
              IF (section(2)=='*BASIS'.AND.section(3)=='.DEFAULT') THEN
C     JHS050707 In the subsection .DEFAULT we give only give the name of a 'LARGE BASIS' type
C               like cc-pVDZ etc. This will be the basisset for all atoms that are not
C               specified in the .SPECIAL section.
                BACKSPACE(filenum)
                CALL READ_VAR(filenum,linelen,str,index, skip=.TRUE.)
                IF (str>' ') THEN
                  tmpbas%name='DEFAULT'
                  tmpbas%comp='LARGE'
                  tmpbas%type='BASIS'
                  tmpbas%file=str
                  tmpbas%next=>NULL()
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  CYCLE
                END IF
              ELSE IF(section(2)=='*BASIS'.AND.
     &                section(3)=='.SPECIAL') THEN
C     JHS050707 For some or all atoms we can give an other basisset than the .DEFAULT
                BACKSPACE(filenum)
                tmpbas%name=' '
                tmpbas%comp=' '
                tmpbas%type=' '
                tmpbas%file=' '
                tmpbas%next=>NULL()
                
                CALL READ_VAR(filenum,linelen,tmpbas%name,index)
                IF (index<=1) CYCLE

                CALL READ_VAR(filenum,linelen,tmpbas%comp,index)
                IF (index<=1) CYCLE
                CALL READ_VAR(filenum,linelen,tmpbas%type,index)                
                IF (index<=1) THEN
                  tmpbas%file=tmpbas%comp
                  tmpbas%comp='LARGE'
                  tmpbas%type='BASIS'
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  CYCLE
                END IF

                SELECT CASE(tmpbas%type)
                CASE ('BASIS','MOLFBAS')                                
                  IF (index<=1) CYCLE              
                  CALL READ_VAR(filenum,linelen,tmpbas%file,index, 
     &                 skip=.TRUE.)
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
  
                CASE('INTGRL','EXPLICIT')
C     FIRST read rest of the line and then read the explicit coefficients.
C     TMP We read the coefficients in in the old way in write commonblock
                  CALL READ_VAR(filenum,linelen,iqm,index)
                  ALLOCATE(JCO(iqm))
                  DO i=1,iqm-1
                    CALL READ_VAR(filenum,linelen,JCO(i),index)
                  END DO                 
                  CALL READ_VAR(filenum,linelen,JCO(i),index,
     &                 skip=.TRUE.)
                  CALL ADD_JCO(tmpbas,JCO)
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  DEALLOCATE(JCO)
  
                CASE('WELLTEMP','EVENTEMP','GEOM','FAMILY','DUALFAMILY')
                  CALL ADD_FAM(filenum,linelen,tmpbas,index)
               
                  CALL READ_VAR(filenum,linelen,iqm,index)
                  ALLOCATE(JCO(iqm))
                  DO i=1,iqm-1
                    CALL READ_VAR(filenum,linelen,JCO(i),index)
                  END DO                 
                  CALL READ_VAR(filenum,linelen,JCO(i),index,
     &                 skip=.TRUE.)
                  CALL ADD_JCO(tmpbas,JCO)
                  CALL ADD_BASIS(bas,firstbas,tmpbas)
                  DEALLOCATE(JCO)
  
                CASE DEFAULT
                  READ(filenum,*)
                  CYCLE
                END SELECT
C              ELSE IF (section(2)=='*BASIS') THEN
C     JHS050707  TO BE IMPLEMENTED Specify an other coordinate system than cartesian.
C                SELECT CASE (section(3))
C                CASE ('.SPHERICAL'); coordintesys='S'
C                CASE ('.OWN'); 
C                 
      
             ELSE IF(section(2)=='*SYMMETRY') THEN
                BACKSPACE(filenum)
                index=1
                CALL READ_VAR(filenum,linelen,NSYMOP,index,skip=.FALSE.)                
                do i=1,3
                  do j=1,3
                    if (index<=0) then
                      write(LUPRI,*) 
     &                   'ERROR:hermodules.F reading symmetry.'
            CALL QUIT('*** ERROR READ_BASIS, symmetry not
     &                     correctly specified')
                    CALL READ_VAR(filenum,linelen,kasym(i,j),
     &           index,1,.FALSE.)
                    end if
                  end do
                end do
                  
              ELSE IF (section(2)=='*CENTERS') THEN
C     JHS050707 We specify the charges, the others are found in the table of nucdata.h
                BACKSPACE(filenum)
                name=' '
                charge=0.0
                nucexp=0.0
  
                CALL READ_VAR(filenum,linelen,name,index)
                IF (.NOT.(index>1)) CYCLE
                CALL READ_VAR(filenum,linelen,charge,index)
                IF (index>1) CALL READ_VAR(filenum,linelen,nucexp,index)
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
              ELSE IF(section(2)=='*UNITS') THEN
C     JHS050707 If a unit is specified, the coordinates in the .XYZ file are given in these 
C               units and we have to transform them to a.u.
C               The conversion factor may be given in Angstrom typed as 'ANGSTROM' or 'A'
C               Note: Angstroms are known and can be specified by 'A' or 'ANGSTROM'
                IF (str==' ') CYCLE
                SELECT CASE(str)
                CASE('ANGSTROM', 'A')
                  UNITS=1/XTANG
                  unitname='angstrom'
                CASE('AU')
                  UNITS=XTANG
                  unitname='au'
                CASE DEFAULT
                  BACKSPACE(filenum)
                  CALL READ_VAR(filenum, linelen,unitname,index)
                  CALL READ_VAR(filenum, linelen,UNITS,index)
                  tmpunit=' '
                  IF(index>1) CALL READ_VAR(filenum,
     &                 linelen,tmpunit,index,skip=.TRUE.)
                  SELECT CASE(tmpunit)
                  CASE ('ANGSTROM', 'A')
                    UNITS=UNITS/XTANG
                  END SELECT
                END SELECT
              END IF
            END IF
          END DO

          IF(.NOT.(ASSOCIATED(firstbas))) THEN
            WRITE(LUPRI,*) 'No basisset given in input file. 
     &            Specify a basisset in  **MOLECULE *BASIS .DEFAULT'
            CALL QUIT('*** ERROR SORT_BASIS, no defaultbasis given in 
     &           input file.')
          END IF
          
          IF(.NOT.(PARHER).AND.IPREAD>1) WRITE(LUPRI,*) 'Done'

C     JHS050707 We now have add all data that was explicitly added in the .INP file.
C     Add default basis.
          IF(.NOT.(PARHER).AND.IPREAD>1)
     &         WRITE(LUPRI,*) 'Processing basissets...'
          CALL SORT_BASIS(firstatom,firstbas,.TRUE.,lupri)
C     Sort the atoms according to the bases.
          CALL SORT_ATOMS(firstatom,firstbas,lupri)
C     Add charges if they are missing
          mol=>firstatom
          prevname=' '
          DO WHILE(ASSOCIATED(mol))
            IF (mol%charge==0) THEN
              IF(.NOT.(mol%name)==prevname) 
     &              charge=FIND_CHARGE(mol%name,lupri)                
              mol%charge=charge
            END IF
  
            prevname=mol%name
            mol=>mol%next
          END DO

          IF(.NOT.(PARHER).AND.IPREAD>1) WRITE(LUPRI,*) 'Done'
          
C     Convert units
          IF(.NOT.(UNITS==1.0) .AND. .NOT.(PARHER)) THEN
            WRITE(LUPRI,*) 'Coordinates entered in ', unitname
            WRITE(LUPRI,*) 'Converting ',unitname,'to atomic units'
            WRITE(LUPRI,*) 'Conversion factor=', UNITS 
            mol=>firstatom
            DO WHILE (ASSOCIATED(mol))
              mol%x=mol%x*UNITS
              mol=>mol%next
            END DO
          END IF

        END SUBROUTINE READ_BASIS


        SUBROUTINE WRITE_ATOMS(lupri)
C     JHS050707 Writes the data contained in the list of atoms to the common blocks.

      IMPLICIT NONE

C     CONSTANTS
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"

C     COMMONBLOCKS
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
         WRITE(lupri,*) 'Number of atoms specified greater max number 
     &        of atoms'
         WRITE(lupri,*) 'Adjust maxorb.h'
         CALL QUIT('*** ERROR (READIN) natom>mxcent')
      END IF

C     nuclei.h
      NUCIND=0
      mol=>firstatom
      natom=0

      DO

        IF (.NOT.(associated(mol))) THEN
          EXIT
        ELSE
          natom=natom+1
C     write the atomic data.
C     If we didn't read the nuclear model somewhere read the nuclear radius.

C    WRITE NUCLEAR SHAPE. TMP To be changed, now always D0 is taken
          IF (GAUNUC.AND.NINT(mol%charge).NE.0) THEN
               
            IF (mol%nucexp.EQ.D0) THEN
              CALL NUCSIZ(mol%charge,mol%nucexp)
            END IF
            IF (mol%nucexp .LE. 1.0D8) THEN
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'WARNING!!!'
              WRITE(lupri,*) 'Nuclear gaussian exponent small: ',
     &                         mol%nucexp
            END IF
            WRITE(lupri,'(/A,F8.3,A,1P,D20.10)')
     &       '  Nuclear Gaussian exponent for atom of charge',
     &       mol%charge,' :',mol%nucexp
          ELSE
            mol%nucexp=D0
          END IF

C     REALLY WRITE COMMONBLOCKS
          NAMN(natom)=mol%name
          CORD(:,natom)=mol%x
          NAMEX(3*natom)=TRIM(mol%name)//' z'
          NAMEX(3*natom-1)=TRIM(mol%name)//' y'
          NAMEX(3*natom-2)=TRIM(mol%name)//' x'
          CHARGE(natom)=mol%charge
          GNUEXP(natom)=mol%nucexp

CMI     ... control printout
          IF (IPREAD.GE.4) THEN
            write(LUPRI,*) 'NAMN(natom)=',NAMN(natom)
            write(LUPRI,*) 'mol%x:',mol%x(1),mol%x(2),mol%x(3)
            write(LUPRI,*) 'CORD(:,natom)',
     &      CORD(1,natom),CORD(2,natom),CORD(3,natom)
          ENDIF
          
          mol=>mol%next
        END IF
      END DO

C      IF (ADDSTO) CALL STOADD()
C
C     ***************************************
C     ***** Process symmetry input data *****
C     ***************************************
C     TMP DETECT SYM?

C      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)

C

      NUCIND=natom

      CALL QEXIT('WRITE_ATOMS')

      END SUBROUTINE WRITE_ATOMS

        SUBROUTINE POS_FILE(infile,lupri,linelen,
     &     section,subsection,name,comp)
C     JHS050707 A function for later use, positions a file to a specific keyword.
        
          INTEGER, INTENT(IN)   :: infile,lupri,linelen
          CHARACTER(LEN=*), INTENT(IN) ::section,subsection,name,comp
          CHARACTER(LEN=linelen) :: str,pressect,pressub
          INTEGER                :: index

          REWIND(infile)
          DO            
            CALL READ_VAR(infile,linelen,str,index, skip=.TRUE.)          
            IF (index<0) THEN
              WRITE(lupri,*) '*** ERROR KEY NOT FOUND. ***'
              CALL QUIT('*** ERROR KEY NOT FOUND. ***')
            ELSE IF (str(1:2)=='**') THEN
              pressect=str              
            ELSE IF (str(1:1)=='*') THEN
              pressub=str
            ELSE IF (pressect==section.AND.pressub==subsection) THEN
              IF (str==name) THEN
                BACKSPACE(infile)
                CALL READ_VAR(infile,linelen,str,index)
                CALL READ_VAR(infile,linelen,str,index,skip=.TRUE.)
                IF (str==comp) EXIT
              END IF
            END IF
          END DO
        END SUBROUTINE POS_FILE

        SUBROUTINE WRITE_BASIS(infile,lupri,WORK,LWORK,
     &                  KANG,KSETS,KBLOCK,KPRIM,HERMIT)
C     JHS050707 Writes the contained in the basisset list to the commonblocks.
C               Note the routine is higly similar tot the subroutine READI1 for
C               .MOL files, the main difference is the split between reading the files
C               and writing them to the commonblocks
      IMPLICIT NONE

C     CONSTANTS
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "codata.h"

C     COMMONBLOCKS
#include "clsfmm.h"
#include "cbirea.h"
      REAL(KIND=8), PARAMETER :: THRMIN = 1.D-15,DSM=1.0D-30

      INTEGER,INTENT(IN)   :: lupri,infile,KANG,KSETS,KBLOCK,KPRIM
      INTEGER,INTENT(INOUT)::LWORK

      INTEGER      :: KLAST,i,j,IBSFLAG, KAOVEC, NBLOCK,
     &                katom,basindex,index

      REAL(KIND=8) :: QEXP=D0
      REAL         :: MXSPD
      LOGICAL      :: HERMIT,DOOWN,ANG,NOORBTS,NOFITFS, 
     &                CNTBAS=.FALSE.,DIRCON(2)=.FALSE.,IFXYZ(3)
      CHARACTER(LEN=72) :: TTITLE(2)=' '
      CHARACTER(LEN=11) :: SYMTXT
      CHARACTER(LEN=15) :: CLASS
      CHARACTER(LEN=namelen) :: prevname

      TYPE(basistype),POINTER       :: baspointer
      TYPE(atom),POINTER            :: mol

      REAL(KIND=8):: WORK(LWORK) 
      REAL(KIND=8),ALLOCATABLE:: 
     &               ALPHA(:,:,:), 
     &               CPRIM(:,:,:,:),
     &               CPRIMU(:,:,:,:)
      INTEGER,ALLOCATABLE :: NONT(:),JCO(:,:,:),JCO2(:,:),IBLOCK(:),
     &                       NUC(:,:),IQM(:,:),NRC(:,:),NBLCK(:,:),
     &                       ISGEN(:)
      LOGICAL,ALLOCATABLE :: SEG(:,:)
      CHARACTER(LEN=80),ALLOCATABLE :: BASREF(:,:,:)

C     JHS050707 For type conversion purposes.
      INTEGER :: int8
      integer blk

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
C#include "primit.h"
C#include "shells.h"
C#include "symmet.h"
C#include "aosotr.h"
C#include "infpar.h"
C#include "dcbrkb.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dgroup.h"
#else
#include "gnrinf.h"
      LOGICAL LINEAR
#endif
C#include "ibtfun.h"
      
      
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
      
C     JHS050707 TMP Initialization (NEEDED?)

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

C     ccom.h
      THRS=THRMIN
C     JHS050707 TMP change this for other coordinte systems
      DOCART=.TRUE.
      DOSPHE=.FALSE.
      DOOWN=.FALSE.
      NHTYP=0

C     nuclei.h
      CALL SPHINP(infile, WORK, LWORK, DOOWN, MXSPD)

      baspointer=>firstbas
      mol=>firstatom
      basindex=0

      DO
C     Loop over all basissets
        NOORBTS=.TRUE.
        NOFITFS=.TRUE.

        IF (.NOT.(associated(baspointer))) THEN
          EXIT
        ELSE
          basindex=basindex+1
          NONT(basindex)=count_atoms(firstatom,baspointer%name)
          BASREF(1,basindex,1) = ' Reference not found in input'
          BASREF(2:10,basindex,1:3) = ' Not initialized'

C     JHS060707 TMP We could introduce some additional error checks?
          SELECT CASE (baspointer%type)
          CASE ('INTGRL','EXPLICIT')
            IBSFLAG=1
            IQM(basindex,1)=SIZE(baspointer%jco)
            DO i=1,IQM(basindex,1)
              JCO(i,basindex,1)=baspointer%jco(i)
CMI           write(lupri,*) 'JCO(i,basindex,1)=',JCO(i,basindex,1)
CMI           call flshfo(lupri)
            END DO           
            BASNAM=baspointer%file
            KAOVEC = KBLOCK + 1 - IBLOCK(1)

C     POSITION the input file.
            CALL POS_FILE(infile, lupri, GET_LENGTH(infile),
     &           '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
            CALL GTOINP(infile,IQM(basindex,1),JCO(1,basindex,1),
     &           NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &           SEG(IBLOCK(1),1),
     &           ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &           CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &           NBLOCK,KAOVEC,KPRIM)
            NBLCK(basindex,1) = NBLOCK
            BASREF(1,basindex,1) = 
     &            "Basis set typed explicitly in input file "
          CASE ('MOLFBAS')
            IBSFLAG=2
            BASNAM = baspointer%file
            KLAST = 1 + (1+KPRIM)*KPRIM
            IF (KLAST.GT.LWORK)
     &        CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
            DIRCON(1) = .TRUE.
            CALL BFGINP(basindex,IBLOCK(1),1,WORK(1),WORK(1+KPRIM),
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,
     &           baspointer%file,BASREF(1,1,1))
          CASE ('BASIS')
            IBSFLAG=3
            BASNAM = baspointer%file
            KAOVEC = KBLOCK + 1 - IBLOCK(1)

            NBLOCK=basindex+1
            mol=>find_atom(mol,baspointer%name)
            CALL BASLIB(IQM(basindex,1),JCO(1,basindex,1),
     &             NUC(IBLOCK(1),1),
     &             NRC(IBLOCK(1),1),SEG(IBLOCK(1),1),
     &             ALPHA(1,IBLOCK(1),1),
     &             CPRIM(1,1,IBLOCK(1),1),CPRIMU(1,1,IBLOCK(1),1),
     &             NBLOCK,ISGEN(IBLOCK(1)),
     &             KAOVEC,KPRIM,mol%charge,
     &             1,DSM,BASREF(1,basindex,1))
            NBLCK(basindex,1) = NBLOCK
          CASE ('WELLTEMP','EVENTEMP','GEOM')
            IBSFLAG=4
            CALL POS_FILE(infile, lupri, GET_LENGTH(infile),
     &           '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
            CALL GENFAMEXP(infile,0)
            i=SIZE(baspointer%fam)-1            
            FAMPAR(1:i)=baspointer%fam(1:i)
            NFAMEXP(1)=TRANSFER(baspointer%fam(i+1),int8)
            KAOVEC = KBLOCK + 1 - IBLOCK(1)
            CALL FAMBAS(infile,IQM(basindex,1),JCO(1,basindex,1),
     &           NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &           SEG(IBLOCK(1),1),
     &           ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &           CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &           NBLOCK,KAOVEC,KPRIM)
            NBLCK(basindex,1) = NBLOCK
            BASREF(1,basindex,1) = 
     &   "Well-tempered basis set typed explicitly in input file "
          CASE ('FAMILY','DUALFAMILY')
            IBSFLAG=4
            CALL POS_FILE(infile, lupri, GET_LENGTH(infile),
     &           '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
            IF (baspointer%type=='FAMILY') i=1
            IF (baspointer%type=='DUALFAMILY') i=2
            CALL GENFAMEXP(infile,i)
            NFAMEXP(1:i)=TRANSFER(baspointer%fam(1:i),int8)
            KAOVEC = KBLOCK + 1 - IBLOCK(1)
            CALL FAMBAS(infile,IQM(basindex,1),JCO(1,basindex,1),
     &           NUC(IBLOCK(1),1),NRC(IBLOCK(1),1),
     &           SEG(IBLOCK(1),1),
     &           ALPHA(1,IBLOCK(1),1),CPRIM(1,1,IBLOCK(1),1),
     &           CPRIMU(1,1,IBLOCK(1),1),ISGEN(IBLOCK(1)),
     &           NBLOCK,KAOVEC,KPRIM)
            NBLCK(basindex,1) = NBLOCK
            BASREF(1,basindex,1) = 
     &   "Well-tempered basis set typed explicitly in input file "
          CASE DEFAULT
            WRITE(lupri,*) 'BASISSET specification not implemented'
            CALL QUIT('ERROR: READING BASISSET')
          END SELECT

C     PRINT BASISSET
          NOORBTS=NOORBTS.AND.IQM(basindex,1).EQ.0
          IF (IQM(basindex,1).GT.0) THEN
            WRITE(lupri,'(3X,A)') 'Large component set:'
            WRITE(lupri,'(6X,A,I5,A,12I5)')
     &      'Max.ang.quantum no.:',(IQM(basindex,1)-1),
     &      '  Blocks:',(JCO(i,basindex,1), i=1,  IQM(basindex,1))
            NHTYP  =   MAX(NHTYP,IQM(basindex,1))
          ENDIF

          DO i = 1, NBLCK(basindex,1)
             IF ((NUC(i,1).NE.NRC(i,1)).AND.(NRC(i,1).GT.0)) 
     &            CNTBAS = .TRUE.
          END DO
          IF ( CNTBAS .AND. (IBSFLAG .EQ. 3 .OR. IBSFLAG .EQ. 1 )) THEN
             DIRCON(1) = .TRUE.
          END IF
          
          
          IF (ASSOCIATED(baspointer%next)) THEN
            IF (baspointer%next%name==baspointer%name) 
     &            baspointer=>baspointer%next
C     JHS060707 Move to small component.
C               TMP We can extend it to include a search for the small component
          END IF

          IF (.NOT.(baspointer%comp(1:5)=='SMALL')) THEN 
            WRITE(lupri,*) 'Small component not specified.'
            WRITE(lupri,*) 'Assuming kinetic balance.'
            IBSFLAG=4
            CALL KINBAL(basindex,IBLOCK(1),IBLOCK(2),
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,
     &           INPTST,CNTBAS)
            BASREF(1,basindex,2) = ' Derived from large component'

          ELSE IF (baspointer%comp=='INTGRL'.OR.
     &           baspointer%comp=='EXPLICIT') THEN
            IBSFLAG=1
            IQM(basindex,2)=SIZE(baspointer%jco)
            DO i=1,IQM(basindex,2)
              JCO(i,basindex,2)=baspointer%jco(i)
            END DO            
            BASNAM=baspointer%file
            KAOVEC = KBLOCK + 1 - IBLOCK(2)
C     POSITION the input file.
            CALL POS_FILE(infile, lupri, GET_LENGTH(infile),
     &           '**BASIS', '*BASIS', baspointer%name,baspointer%comp)
            CALL GTOINP(infile,IQM(basindex,2),JCO(1,basindex,2),
     &           NUC(IBLOCK(2),2),NRC(IBLOCK(2),2),
     &           SEG(IBLOCK(2),2),
     &           ALPHA(1,IBLOCK(2),2),CPRIM(1,1,IBLOCK(2),2),
     &           CPRIMU(1,1,IBLOCK(2),2),ISGEN(IBLOCK(2)),
     &           NBLOCK,KAOVEC,KPRIM)
            NBLCK(basindex,2) = NBLOCK
            BASREF(1,basindex,2) = 
     &            "Basis set typed explicitly in input file "
          ELSE IF (baspointer%type(1:7)=='MOLFBAS') THEN
            IBSFLAG=2
            BASNAM=baspointer%file
            KLAST = 1 + (1+KPRIM)*KPRIM
            IF (KLAST.GT.LWORK)
     &        CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
            DIRCON(1) = .TRUE.
            CALL BFGINP(basindex,IBLOCK(1),1,WORK(1),WORK(1+KPRIM),
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,
     &           baspointer%file,BASREF(1,1,1))

          ELSE IF (.NOT.(baspointer%type(1:6)=='KINBAL')) THEN
            IBSFLAG=4
            CALL KINBAL(basindex,
     &           IBLOCK(1),IBLOCK(2),
     &           IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,
     &           INPTST,CNTBAS)
            BASREF(1,basindex,2) = ' Derived from large component'
          END IF          

          NOORBTS = NOORBTS.AND.IQM(basindex,2).EQ.0

          IF (IQM(basindex,2).GT.0) THEN
            WRITE(lupri,'(3X,A)') 'Small component set:'
            WRITE(lupri,'(6X,A,I5,A,12I5)')
     &        'Max.ang.quantum no.:',(IQM(basindex,2)-1),
     &        '  Blocks:',(JCO(i,basindex,2), i=1,IQM(basindex,2))
          ENDIF

          NHTYP = MAX(NHTYP,IQM(basindex,2))
          
          IF (IQM(basindex,2).GT.0) THEN
            IF(INPTST) THEN
               WRITE(lupri,'(3X,A)') 'Small component set:'
               WRITE(lupri,'(6X,A,I5,A,12I5)')
     &              'Max.ang.quantum no.:',(IQM(basindex,2)-1),
     &              '  Blocks:',(JCO(i,basindex,2), i=1,IQM(basindex,2))
            ENDIF
            NHTYP     = MAX(NHTYP,IQM(basindex,2))
C            IF(NHTYP.GT.MXQN) GOTO 5000
          ENDIF
          
          DO i = 1, NBLCK(basindex,2)
            IF ((NUC(i,2).NE.NRC(i,2)).AND.(NRC(i,2).GT.0)) 
     &            CNTBAS = .TRUE.
          END DO
          IF ( CNTBAS .AND. IBSFLAG .NE. 4 ) DIRCON(2) = .TRUE.


          IBLOCK(1:2)=IBLOCK(1:2)+NBLCK(basindex,1:2)
          NOORBT(basindex+1:NONT(basindex)+basindex) = NOORBTS
          
          DO i=3,4
C     JHS060707 TMP Include also the reading of fitsets.
            IBSFLAG=2
            IBLOCK(i)=IBLOCK(i)+NBLCK(basindex,i)
            
          END DO
          prevname=baspointer%name
          DO WHILE (associated(baspointer))
            IF (.NOT.(baspointer%name==prevname)) EXIT
            baspointer=>baspointer%next
          END DO
        END IF
      END DO

      LINEAR = .FALSE.
      CLASS = 'N/A'
      IF (ADDSYM) THEN
         CALL SYMADD(WORK,LWORK,basindex,NSYMOP,NONT,KATOM,KASYM,CLASS)
         IF (CLASS(3:4).EQ.'oo') LINEAR = .TRUE.
      END IF
C
C     ****************************************************
C     ***** Adds an STO-3G basis if Extended Hueckel *****
C     *****        wanted as starting orbitals       *****
C     ****************************************************
C
C      IF (ADDSTO) CALL STOADD()
C
C     ***************************************
C     ***** Process symmetry input data *****
C     ***************************************
C
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, CLASS)

C
C     ***************************************************
C     ***** Process orbital and geometry input data *****
C     ****************************************************
C
      CALL BASPRO(WORK,LWORK,basindex,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,NRC,
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,KBLOCK,KPRIM,
     &            DOOWN, blk)
C
C     *************************************************
C     ***** Print orbital and geometry input data *****
C     *************************************************
C

      CALL BASOUT(WORK,LWORK,basindex,NONT,IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &            KBLOCK,KPRIM,HERMIT, blk)
C
C     TMP copy input to .out file.

C     Default threshold for distinguishing classical and
C     non-classical contributions to Coulomb integrals

      THRCLS = THRS
C
C     ****************************************
C     ***** Output on LUONEL *****************
C     ****************************************
C

      IF (HERMIT) CALL WRONEL(TTITLE,basindex,NONT,IQM,IFXYZ,KATOM,JCO2,
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
          WRITE(lupri,'(a,t13,a,t25,a,t37,a,t49,a,t61,a)') 'Label',
     &         'x','y','z','charge','nuc exp'
          mol=>firstatom
          DO
            IF (associated(mol)) THEN
              WRITE(lupri,'(a,t9, 5f12.6)') 
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
          WRITE(lupri,'(a,t9,a,t41,a,t73,a)')
     &         'Label','component', 'input type','file'
          bas=>firstbas
          DO
            IF (associated(bas)) THEN
              WRITE(lupri,'(a,t9,a,a,a)') 
     &          TRIM(bas%name),bas%comp,bas%type,bas%file
              IF(associated(bas%jco).OR.associated(bas%fam))
     &           WRITE(lupri,*) 'Exponents:', bas%jco, bas%fam
              bas=>bas%next


            ELSE
              EXIT
            ENDIF
          END DO
        WRITE(lupri,'(a)') ' '

        END SUBROUTINE PRINT_BASIS

        SUBROUTINE DEL_ATOMS(firstatom, lupri)
C     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          !TYPE(atom),POINTER,INTENT(INOUT) :: firstatom
          TYPE(atom),POINTER :: firstatom
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          TYPE(atom), POINTER :: mol

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
C     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          !TYPE(basistype),POINTER,INTENT(INOUT) :: firstbas
          TYPE(basistype),POINTER :: firstbas
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          TYPE(basistype), POINTER :: bas

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
C     FREES THE MEMORY RESERVED FOR THE LIST OF ATOMS.
        
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          CALL DEL_ATOMS(firstatom, lupri)
          CALL DEL_BASIS(firstbas, lupri)

        END SUBROUTINE CLEAR_XYZ


        SUBROUTINE READ_XYZ(molfile,inpfile, lupri)
#include "infpar.h"        
C     THE MAIN SUBROUTINE THAT READS AN XYZ FILE

          INTEGER, INTENT(IN) :: molfile,inpfile,lupri

          CALL READ_ATOMS(molfile,GET_LENGTH(molfile),firstatom,lupri)
          CALL READ_BASIS(inpfile, GET_LENGTH(inpfile),lupri)

          IF(.NOT.(PARHER)) CALL PRINT_ATOMS(firstatom,lupri)
          IF(.NOT.(PARHER)) CALL PRINT_BASIS(firstbas,lupri)
            


        END SUBROUTINE READ_XYZ


      END MODULE READFILE







