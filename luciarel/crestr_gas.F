      SUBROUTINE CRESTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NACOB,
     &                  IPRNT)
*
* A type of strings containing NEL electrons are given
* set up all possible ways of adding an electron to this type of strings
*
*------------------------------------------
*  Jeppe Olsen
*
*   Fixed machine dependent problem with IF blocks.
*   They are now split up so that array dimensions are
*   not violated.
*  Timo Fleig, January 2001
*------------------------------------------
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL + 1 electrons
* NEWORD : Reordering array for N+1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL+1 strings
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
C     INTEGER STRING,TI,TTO,OCC,UNOCC,STRIN2,Z
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(NORB,NSTINI),TTO(NORB,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =  000
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' CRESTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
C     WRITE(6,*) ' Reorder array NEWORD '
C     CALL IWRTMA(NEWORD,1,NSTINO,1,NSTINO)
*
      DO 1000 ISTRIN = 1,NSTINI
C?    write(6,*) ' Input string ',istrin,(string(i,istrin),i=1,nel)
        DO 100 IORB = IORBOF, IORBOF-1+NORB
C?      write(6,*) ' orbital ',iorb
           IPLACE = 0
           IF(NEL.EQ.0) THEN
             IPLACE = 1
             GOTO 11
           ELSE IF ( NEL .NE. 0 ) THEN
            DO 10 IEL = 1, NEL
              IF (IEL.eq.1.AND.IORB.lt.STRING(1,ISTRIN)) THEN
                IPLACE = 1
                GOTO 11
              ELSE IF (IEL.EQ.NEL.AND.IORB.GT.STRING(IEL,ISTRIN)) then
                IPLACE = IEL + 1
                goto 11
              else if (IORB.eq.STRING(IEL,ISTRIN)) THEN
                IPLACE = 0
                GOTO 11
              end if
              if (IEL.LT.NEL) then
                if (IORB.GT.STRING(IEL,ISTRIN).AND.
     &              IORB.LT.STRING(IEL+1,ISTRIN)) THEN
                  IPLACE = IEL+1
                  GOTO 11
                end if
              end if
   10       CONTINUE
           END IF
   11     CONTINUE
*
C?        write(6,*) ' iplace = ', iplace
          IF(IPLACE.NE.0) THEN
*. Generate next string
            DO 30 I = 1, IPLACE-1
   30       STRIN2(I) = STRING(I,ISTRIN)
            STRIN2(IPLACE) = IORB
            DO 40 I = IPLACE,NEL
   40       STRIN2(I+1) = STRING(I,ISTRIN)
C?          write(6,*) ' updated string (STRIN2) '
C?          call iwrtma(STRIN2,1,NEL+1,1,NEL+1)
            JSTRIN = ISTRNM(STRIN2,NACOB,NEL+1,Z,NEWORD,1)
C?          write(6,*) ' corresponding number ', JSTRIN
*
            TTO(IORB-IORBOF+1,ISTRIN) = JSTRIN
            IIISGN = (-1)**(IPLACE-1)
            IF(LSGSTR.NE.0)
     &      IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
            IF(IIISGN .EQ. -1 )
     &      TTO(IORB-IORBOF+1,ISTRIN) = - TTO(IORB-IORBOF+1,ISTRIN)
            TI(IORB-IORBOF+1,ISTRIN ) = IORB
          END IF
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(6,*) ' Output from CRESTR : '
        WRITE(6,*) '==================='
*
        WRITE(6,*)
        WRITE(6,*) ' Strings with an electron added  '
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,NORB)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,NORB)
        END DO
      END IF
*
      RETURN
      END
