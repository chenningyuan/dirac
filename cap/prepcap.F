C  /* Deck capinp */
      SUBROUTINE CAPINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C  Input section for CAP-module
C
C  Called from HAMINP (dirrdn.F)
C
C  Written by  Miro ILIAS, Febr. 2007, Tel Aviv
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "dcbham.h"
#include "dcbcap.h"
#include "capint.h"
      PARAMETER (NTABLE = 20)
      PARAMETER (D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF,RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.INPQD ','.INPCAP','.PRINT ','.OVLCAP',
     &            '.NOOVLP','.NONREL','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
#include "memint.h"
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF

      CALL TITLER('CAP operator defined','*',127)

C
C        *** Initialize ****
C
      CAPOVLP = .FALSE.
      NOOVLP = .FALSE.
      DONRCAP = .FALSE.
CMI   DONRCAP = .TRUE.
      IF (LEVYLE) DONRCAP = .TRUE.
C
C
C     Process input for CAPINP
C     ========================
C
      NEWDEF = (WORD .EQ. '*CAP   ')

      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GOTO ( 1, 2, 3, 4,
     &                      5, 6, 7, 8) , I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in CAPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in CAPINP.')
 1             CONTINUE
C&&& INPQD: read the PW number 'theta'
               READ(LUCMD,*) PW
               GO TO 100
 2             CONTINUE
C&&& INPCAP: read integration parameters
              DO I = 1, 3
                READ(LUCMD,*) X0(I),LAMDA(I)
                READ(LUCMD,*) NGR(I)
                IF (NGR(I).GT.MAXSPL) THEN
                 WRITE(LUPRI,'(2X,A,I4)')
     & 'CAPINT: Increase the parameter MAXSPL in the dcbcap. file !'//
     & ' Current value is ',MAXSPL
                 CALL QUIT(
     & 'CAPINT: Increase dimension parameter MAXSPL !')
                ENDIF
                READ(LUCMD,*)(GINTPT(I,J),J=1,NGR(I))
                READ(LUCMD,*)(NGRS(I,J),J=1,NGR(I)-1)
              ENDDO
              GO TO 100
 3               CONTINUE
C&&& Read print level
                 READ(LUCMD,*) IPRCAP
                 GO TO 100
 4               CONTINUE
C&&& Do only the comparison between DIRAC and CAP overlap integrals
                 CAPOVLP = .TRUE.
                 GO TO 100
 5               CONTINUE
C&&& NO comparison between DIRAC and CAP overlap integrals before calculation of CAP
                 NOOVLP  = .TRUE.
                 GO TO 100
 6               CONTINUE
C&&& Set the NONRELATIVISTIC version of the CAP operator
                 DONRCAP = .TRUE.
                 GO TO 100
 7               CONTINUE
C&&& XXXXXX
                 GO TO 100
 8               CONTINUE
                 GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in CAPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in CAPINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF(.NOT.CAP) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(4X,A)')
     & 'CAPINP: Set-up for CAP calculations'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)') ' * Print level: ',IPRCAP
C
      IF (.NOT.CAPOVLP) THEN
C
        IF (DONRCAP) THEN
         WRITE(LUPRI,'(2X,A)')
     &  '* Using nonrelativistic form of the CAP operator.'
         WRITE(LUPRI,'(A)')
     &  ' * CAP integrals  according to',
     &  '     Sajeev ....         '
        ELSE
         WRITE(LUPRI,'(2X,A)')
     &  '* Using four-component relativistic form of the CAP operator.'
        ENDIF

        WRITE(LUPRI,'(2X,A,D25.10)')
     &  '* read INPQD (theta,PW) = ',PW
        WRITE(LUPRI,'(2X,A)') 
     & 'read integration parameters:'
CMI    ...  write out integration parameters
        DO I = 1, 3
          WRITE(LUPRI,'(2X,A,2F26.8)')
     &    '* XO, LAMDA:',X0(I),LAMDA(I)
          WRITE(LUPRI,*) NGR(I)
          WRITE(LUPRI,'(2X,A)')
     &    '* integration subintervals:'
          WRITE(LUPRI,*) (GINTPT(I,J),J=1,NGR(I))
          WRITE(LUPRI,'(2X,A)')
     &    '* number of integration points per interval:'
          WRITE(LUPRI,*) (NGRS(I,J),J=1,NGR(I)-1)
        ENDDO
        IF (NOOVLP) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '* no comparison between CAP-DIRAC overlap integrals !'
        ENDIF
       ELSE
         WRITE(LUPRI,'(2X,A)') 
     &   '* only comparison between CAP and DIRAC overlap integrals.'
       ENDIF

  999 CONTINUE
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dosajeev */
      SUBROUTINE DOSAJEEV(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Prepare arrays as parameters for the core CAP routine
C and run the CAP routine which prepares PEMAT binary file
C
C Called from F90 routine CAP_STUFF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "dummy.h"
      DIMENSION WORK(LWORK)

      CALL QENTER('DOSAJEEV')

#include "memint.h"
      ISIZE1=200
      CALL MEMGET('INTE',KBASIS_N,ISIZE1,WORK,KFREE,LFREE)
      ISIZE2=250*3
      CALL MEMGET('INTE',KANG_POW,ISIZE2,WORK,KFREE,LFREE)

      ISIZE3=250
      CALL MEMGET('REAL',KORDEXPON,ISIZE3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORDCOE,ISIZE3,WORK,KFREE,LFREE)

      ISIZE4=200*3
      CALL MEMGET('REAL',KRMS,ISIZE4,WORK,KFREE,LFREE)

      IF (CAPOVLP) THEN
       CALL MEMGET('REAL',KPE_AOr,NNBBASX,WORK,KFREE,LFREE)

       CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr),DUMMY,        
     &            DUMMY,DUMMY, DUMMY,DUMMY, DUMMY,DUMMY,
     &            WORK(KFREE),LFREE)
      ELSE
       IF (DONRCAP) THEN
CMI    ISIZE5=25000
CMI    CALL MEMGET('REAL',KPE_AOr,ISIZE5,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPE_AOr,NNBBASX,WORK,KFREE,LFREE)
CMI    CALL MEMGET('REAL',KPE_AOi,ISIZE5,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPE_AOi,NNBBASX,WORK,KFREE,LFREE)

        CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr),WORK(KPE_AOi),
     &            DUMMY,DUMMY, DUMMY,DUMMY, DUMMY,DUMMY,
     &            WORK(KFREE),LFREE)

      ELSE
C       ... four-component CAP operator ...
        CALL MEMGET('REAL',KPEX_AOr,NNBBASX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPEX_AOi,NNBBASX,WORK,KFREE,LFREE)
C       ...
        CALL MEMGET('REAL',KPEY_AOr,NNBBASX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPEY_AOi,NNBBASX,WORK,KFREE,LFREE)
C       ...
        CALL MEMGET('REAL',KPEZ_AOr,NNBBASX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPEZ_AOi,NNBBASX,WORK,KFREE,LFREE)

        CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            DUMMY, DUMMY,     
     &          WORK(KPEX_AOr),WORK(KPEX_AOi),
     &          WORK(KPEY_AOr),WORK(KPEY_AOi),
     &          WORK(KPEZ_AOr),WORK(KPEZ_AOi),
     &            WORK(KFREE),LFREE)

      ENDIF
      ENDIF

C     ... release all the allocated memmory
      CALL MEMREL('DOSAJEEV',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('DOSAJEEV')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      LOGICAL FUNCTION CHPRES(IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"

      LOGICAL FILES, FILES2, IS_PEMAT1, IS_PEMAT2

C  ... check presence of files CAPLBL.inp PEMAT CAPDIRSYM.inp CAPDIR.inp
      FILES = .FALSE.
      INQUIRE(FILE='CAPLBL.inp',EXIST=FILES)
      FILES2 = FILES

C   ... treat PEMAT / PEMAT_ovlp
C     FILES = .FALSE.
      IS_PEMAT1 = .FALSE.
      INQUIRE(FILE='PEMAT',     EXIST=IS_PEMAT1)
C     INQUIRE(FILE='PEMAT_ovlp',EXIST=IS_PEMAT2)
C     IF (CAPOVLP) THEN
C       FILES2 = FILES2.AND.(IS_PEMAT2)
C     ELSE
C       FILES2 = FILES2.AND.(IS_PEMAT1)
C     ENDIF

C     FILES = .FALSE.
C     INQUIRE(FILE='CAPDIRSYM1.inp',EXIST=FILES)
C     FILES2 = FILES2.AND.FILES

C     FILES = .FALSE.
C     INQUIRE(FILE='CAPDIRSYM2.inp',EXIST=FILES)
C     FILES2 = FILES2.AND.FILES

      FILES = .FALSE.
      INQUIRE(FILE='CAPDIR.inp',EXIST=FILES)
      FILES2 = FILES2.AND.FILES

      CHPRES = FILES2

      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capadapt */
      SUBROUTINE CAPADAPT(ICAPBFO,IOBF1_2,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Adapt the CAP integrals from the previously obtained UNformatted file
C   PEMAT for the DIRAC calculations.
C
C   For this it is using generated text files
C
C   Called from  prepcap90.F90/CAP_STUFF
C
C   Written by Miro ILIAS, Tel Aviv, 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "dgroup.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"

      DIMENSION WORK(LWORK),ICAPBFO(1:NBASIS),ICAPBFO1_2(1:NBASIS)
      CHARACTER*8 LBL

      CALL QENTER('CAPADAPT')
  
C ======== first real part (CAPRE or CAP_OVLP)

#include "memint.h"
      CALL MEMGET('REAL',KPMRSAT,NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMRSA,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMR,N2BBASX,WORK,KFREE,LFREE)

      IF (.NOT.CAPOVLP) THEN
       IF (DONRCAP) THEN
C        ... non-relativistic CAP integrals
         CALL CAPREAD('CAP_RE  ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &          ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

         CALL CAPREAD('CAP_IM  ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &          ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

       ELSE
C        ...  relativistic CAP integrals
          CALL CAPREAD('XCAPD1R ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('XCAPD1I ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('YCAPD1R ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('YCAPD1I ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('ZCAPD1R ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('ZCAPD1I ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

         ENDIF

      ELSE
C      ... overlap integrals
          CALL CAPREAD('CAP_OVLP',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &                 ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

      ENDIF

C     ... release all the memory
      CALL MEMREL('CAPADAPT',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('CAPADAPT')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dircapao */
      SUBROUTINE DIRCAPAO(ICAPBFO,IOBF1_2,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Find out matching between CAP and DIRAC basis functions
C
C  Called from prepcap.F/CAP_STUFF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION WORK(LWORK),ICAPBFO(NBASIS),IOBF1_2(NBASIS)

      CALL QENTER('DIRCAPAO')

#include "memint.h"
C     CALL MEMGET('INTE',KICAPBFO,NBASIS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIOBF1,3*NBASIS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIOBF2,3*NBASIS,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KILC1,NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KILC2,NUCDEP,WORK,KFREE,LFREE)

C     CALL IZERO(WORK(KICAPBFO),NBASIS)
      CALL IZERO(WORK(KIOBF1),3*NBASIS)
      CALL IZERO(WORK(KIOBF2),3*NBASIS)
      CALL IZERO(WORK(KILC1),NUCDEP)
      CALL IZERO(WORK(KILC2),NUCDEP)

      CALL GET_DIR_CAP_AO_MATCH(ICAPBFO,IOBF1_2,WORK(KIOBF1),
     &       WORK(KIOBF2),WORK(KILC1),WORK(KILC2))

      CALL MEMREL('DIRCAPAO',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('DIRCAPAO')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck get_dir_cap_ao_match */
      SUBROUTINE GET_DIR_CAP_AO_MATCH(ICAPBFO,IOBF1_2,
     &             IOBF1,IOBF2,ILC1,ILC2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Called from  prepcap90.F90/CAP_DIRAC_AO
C
C   Needs CAPLBL.inp and CAPDIR.inp files, where CAPLBL.inp is generated in SAJEEV's
C  routine !
C
C   Returns filled ICAPBFO,IOBF1_2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "dgroup.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"

      DIMENSION IOBF1(1:3,1:NBASIS),IOBF2(1:3,1:NBASIS),
     &          ICAPBFO(1:NBASIS), IOBF1_2(1:NBASIS)
      DIMENSION ILC1(1:NUCDEP),ILC2(1:NUCDEP)
      INTEGER AONUM,CMNINDX
      CHARACTER*6 AO1,AO2
      LOGICAL CHPRES

      CALL QENTER('GET_DIR_CAP_AO_MATCH')

      IF (.NOT.CHPRES(IPRCAP)) THEN
         CALL QUIT('GET_DIR_CAP_AO_MATCH:...not all files for CAP !')
      ENDIF

      LUCD = 0
      CALL GPOPEN(LUCD,'CAPDIR.inp','OLD',' ','FORMATTED',
     &       IDUMMY,.FALSE.)

      LUCL = 0
      CALL GPOPEN(LUCL,'CAPLBL.inp','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      IF (IPRCAP.GE.3) THEN
       CALL HEADER('Output from CAPADAPT',-1)
       write(LUPRI,'(2X,A,I4/)')
     & 'Total number of atomic basis functions, NBASIS = ',NBASIS
      ENDIF

C ...  read DIRAC/CAP basis functions ...
C       .. skip text lines
      READ(LUCD,*)
      READ(LUCL,*)
      DO I = 1, NBASIS
C      ... fill in CAP abfu characteristics
        read(LUCL,*) II1,INC1,AO1
        IOBF1(1,I) = INC1 
        IOBF1(2,I) = AONUM(AO1)
        CALL IOBFSET(IOBF1,I,IPRCAP)

C      ... fill in DIRAC abfu characteristics
        read(LUCD,*) II2,INC2,AO2
        IOBF2(1,I) = INC2 
        IOBF2(2,I) = AONUM(AO2)
        CALL IOBFSET(IOBF2,I,IPRCAP)

        IF (IPRCAP.GE.3) THEN
          write(LUPRI,'(I4,A,I3,A,A4,A,I2,A,I3,A,D14.6)')
     &    II1,'> nucleus:',IOBF1(1,I),' funct: ',AO1,'/',IOBF1(2,I),
     &     ' >  CAP unique fu #:',
     &        IOBF1(3,I),
     &        ' factor:  ',GET_CAPAF_FAC(IOBF1(2,I))
          write(LUPRI,'(I4,A,I3,A,A4,A,I2,A,I3)')
     &    II2,'> nucleus:',IOBF2(1,I),' funct: ',AO2,'/',IOBF2(2,I),
     &     ' >dirac unique fu #:',
     &        IOBF2(3,I)
          write(LUPRI,*)
        ENDIF

      ENDDO

      CALL GPCLOSE(LUCD,'KEEP')
      CALL GPCLOSE(LUCL,'KEEP')

C ... according to the DIRAC basis function order fill ICAPBFO
      DO I = 1, NBASIS
        IX = 0
        DO J = 1 , NBASIS
C         IF (IOBF2(1,I).EQ.IOBF1(1,J).AND.
C    &        IOBF2(2,I).EQ.IOBF1(2,J).AND.
C    &        IOBF2(3,I).EQ.IOBF1(3,J)) IX = J        
          IF (IOBF1(1,I).EQ.IOBF2(1,J).AND.
     &        IOBF1(2,I).EQ.IOBF2(2,J).AND.
     &        IOBF1(3,I).EQ.IOBF2(3,J)) IX = J        
        ENDDO

        IF (IX.EQ.0) THEN
          WRITE(LUPRI,*) '  '
          WRITE(LUPRI,*) '  '
          CALL QUIT('CAPADAPT: order does not fit !!')
        ENDIF
 
C       ... save function
        ICAPBFO(I) = IX

      ENDDO

      IF (IPRCAP.GE.2) THEN
       WRITE(LUPRI,'(2X,A,I4)')
     & 'Matching between CAP and DIRAC atomic basis'//
     & ' functions - basic order, NBASIS fu=',NBASIS
       DO I = 1, NBASIS
         WRITE(LUPRI,'(1X,I3,A,I3,A,$)') I,'->',ICAPBFO(I),'|'
         IF (MOD(I,8).EQ.0) WRITE(LUPRI,*)
       ENDDO
       WRITE(LUPRI,*)
      ENDIF

C   ... fill IOBF1_2 array !
      DO I=1,NBASIS
        IOBF1_2(I) = IOBF1(2,I)
      ENDDO

      CALL QEXIT('GET_DIR_CAP_AO_MATCH')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck delete_ls */
      SUBROUTINE DELETE_LS(V,N,IPRINT,DEL_LS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  delete LS,SL blocks of the matrix ! basic functions are in basis DIRAC order
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      REAL*8 V(N,N)
      LOGICAL DEL_LS
      CALL QENTER('DELETE_LS')

      DO I = 1 ,NBASIS
      DO J = 1, NBASIS
       IF (DEL_LS) THEN
C    ... delete LS,SL blocks
         if (i.le.NLARGE.and.j.gt.NLARGE) then
           V(i,j) = 0.0D0
         endif
         if (i.gt.NLARGE.and.j.le.NLARGE) then
           V(i,j) = 0.0D0
         endif
       ELSE
C    ... delete LL,SS blocks
         if (i.le.NLARGE.and.j.le.NLARGE) then
           V(i,j) = 0.0D0
         endif
         if (i.gt.NLARGE.and.j.gt.NLARGE) then
           V(i,j) = 0.0D0
         endif
       ENDIF
      ENDDO
      ENDDO

      CALL QEXIT('DELETE_LS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck getcapim */
      SUBROUTINE GETCAPMO(IMAGPART,CAPI,CAPI1,WORK,KFREE,LFREE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Purpose: get CAP operator in MO basis into CAPI(*) 
C
C   First get it (using PRPMAO) as the real operator
C
C   IMAGPART - if true, read imaginary part of Ham. otherwise real
C
C   Called from ...  (prepcap.F)
C
C   Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      PARAMETER(D0=0.0D0, D2=2.0D0, D1 = 1.0D0, DM1 = -1.0D0)
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprp.h"
#include "dcbxpr.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbcap.h"
#include "dummy.h"

      DIMENSION CAPI(NORBT,NORBT,NZ),CAPI1(NORBT,NORBT,NZ),WORK(*)
      LOGICAL REARR,IMAGPART

      CALL QENTER('GETCAPMO')

      KFRSAV = KFREE

      CALL DZERO(CAPI,NORBT*NORBT*NZ)

C  ...array for MO coefficients
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C boson irreps identification (read with MO coeff.)
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)

      CALL MEMCHK('getcapmo 1',WORK,1)

C     open file for MO coefficients
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GETCAPIM')
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,WORK(KBEIG),
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUMMY,IDUMMY,
     &            TOTERG,IOPT,WORK,KFREE,LFREE)
      ENDIF

C     ... Open property file
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','GETCAPMO')

      IF (IMAGPART) THEN
       IF (DONRCAP) THEN
        INDXPR = IPCAPI
       ELSE
        INDXPR = IPCAP_DI
       ENDIF
      ELSE
       IF (DONRCAP) THEN
        INDXPR = IPCAPR
       ELSE
        INDXPR = IPCAP_DR
       ENDIF
      ENDIF

      ISYM  = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      IOPSY = JBTOF(ISYM-1,1)
      CALL DZERO(CAPI,N2ORBXQ)
      CALL PRPMAT(INDXPR,IOPSY,CAPI,WORK(KCMO),WORK(KBEIG),
     &               ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

      IF (IPRCAP.GE.5) THEN
       WRITE(LUPRI,'(2X,2A)')
     & 'GETCAPMOM: Get MO matrix of property: ',PRPNAM(INDXPR)
       CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

      CALL MEMCHK('getcapmo 2',WORK,1)

CMI   ...MOVE ee-part to the upper left part of the MO-matrix
C     REARR = .FALSE.
      REARR = .true.
C
      IF (REARR.AND..NOT.LEVYLE) THEN

        CALL DZERO(CAPI1,NORBT*NORBT*NZ)

        IOFF_E = 0
        IOFF_EP = 0
        DO IFSYM =1 , NFSYM
         IF ( NORB(IFSYM).GT.NESH(IFSYM) ) THEN
           DO IZ = 1, NZ
             DO I = 1, NESH(IFSYM)
             DO J = 1, NESH(IFSYM)
              CAPI1(I+IOFF_EP,J+IOFF_EP,IZ) = 
     &        CAPI(I+NPSH(IFSYM)+IOFF_EP,J+NPSH(IFSYM)+IOFF_EP,IZ)
             ENDDO
             ENDDO
           ENDDO
           IOFF_EP = NORB(IFSYM)
           IOFF_E  = NESH(IFSYM)
         ENDIF
        ENDDO

       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(2X,2A)')
     &  'GETCAPMO: Rearranged MO matrix of property: ',PRPNAM(INDXPR)
        CALL PRQMAT(CAPI1,NORBT,NORBT,NORBT,
     &           NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
       ENDIF

      ENDIF

      CALL MEMCHK('getcapmo 3',WORK,1)

      CALL MEMREL('GETCAPMO',WORK,KFRSAV,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETCAPMO')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck addcapmom */
      SUBROUTINE ADDCAPMOM(NSTRT,CAPMOM,FMOM,MULTIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: add MOLFDIR like CAPMOM to the FMOM
C
C  On input: CAPMOM - to be added to FMOM
C            MULTIM - flag to decide the multiplication by the imag.unit
C
C  Called from the MRCON1 routine (traone.F).
C
C  Written by Miro ILIAS, Tel Aviv, febr. 2007 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbcap.h"

      DIMENSION CAPMOM(NSTRT*2,NSTRT*2,2), FMOM(NSTRT*2,NSTRT*2,2)
      LOGICAL MULTIM

      CALL QENTER('ADDCAPMOM')

C ... note: FMOM and CAPIM are in MOLFDIR format !
      IF (IPRCAP.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPMOM:  *** Entering core Fock matrix in MOLFDIR'//
     &    ' format, FMOM ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)

          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPMOM:  *** Entering (MOLFDIR like) CAPMOM ***'
          CALL PRQMAT(CAPMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF

C     ... multiply CAPMOM with "i"
      IF (MULTIM) THEN
        DO I = 1, NSTRT*2
        DO J = 1, NSTRT*2
C         i.(a+ib) = -b + i.a
          TEMP = CAPMOM(I,J,2)
          CAPMOM(I,J,2) = CAPMOM(I,J,1)
          CAPMOM(I,J,1) = - TEMP
        ENDDO
        ENDDO
        IF (IPRCAP.GE.5) THEN
         WRITE(LUPRI,'(/3X,A)')
     &   'ADDCAPMOM:  *** Modified i.CAPMOM - mult.by "i"-unit ***'
         CALL PRQMAT(CAPMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &               IPQTOQ(1,0),LUPRI)
        ENDIF
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'ADDCAPMOM: Warning - no multiplication with imaginary unit !'
      ENDIF

C      ...Add CAPMOM into FMOM
      DO I = 1, NSTRT*2
      DO J = 1, NSTRT*2
       TEMP1 = FMOM(I,J,1)
       FMOM(I,J,1) = TEMP1 + CAPMOM(I,J,1)
       TEMP2 =  FMOM(I,J,2) 
       FMOM(I,J,2) = TEMP2 + CAPMOM(I,J,2)
      ENDDO
      ENDDO

      IF (IPRCAP.GE.5) THEN
       WRITE(LUPRI,'(/3X,A)')
     & 'ADDCAPMOM:  *** Core Fock matrix (molfdir format)'//
     & ' with CAPMOM added (NZ=2) ***'
       CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

      CALL QEXIT('ADDCAPMOM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck capmfdadapt */
      SUBROUTINE CAPMFDADAPT(CAP,NSTRT,NSTR,FMO,FMOM,WORK,LWORK) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Do the adaptation of the Fock matrix for CAP-FSCCSD calculations. Add CAP real
C and imaginary contributions.
C
C Input: CAP - logical
C          
C
C Output: FMOM with CAP (re,im) parts added
C
C Called from MRCON1/traone.F
C
C Written by Miro ILIAS, Mau 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbcap.h"
      DIMENSION WORK(LWORK)
      DIMENSION FMO(NORBT,NORBT,NZ),FMOM(NSTRT*2,NSTRT*2,2)
      DIMENSION NSTR(2)
      LOGICAL CAP
      CALL QENTER('CAPMFDADAPT')

      IF (CAP) THEN

       IF (IPRCAP.GE.5) THEN
        WRITE(LUPRI,'(/3X,A,I2/)')
     &   '*** Pure Core Fock matrix (DIRAC format) ***'
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &                      IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/3X,A,I2/)')
     &   '*** Pure Core Fock matrix (MOLFDIR format) ***'
        CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &               IPQTOQ(1,0),LUPRI)
       ENDIF
#include "memint.h"
C         ... allocate space for CAP_IM
       CALL MEMGET('REAL',KCAPI, NORBT*NORBT*NZ,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KCAPI1,NORBT*NORBT*NZ,WORK,KFREE,LFREE)

C      .... get CAP  operator - real part (nonrel/rel) in the MO basis
        CALL GETCAPMO(.FALSE.,WORK(KCAPI),WORK(KCAPI1),WORK,KFREE,LFREE)
C       ... allocate space for MOLFDIR like CAP
        CALL MEMGET('REAL',KCAPIM,NSTRT*NSTRT*8,WORK,KFREE,LFREE)
C       ... turn it in the MOLFDIR like form
        CALL MOLTODIR(WORK(KCAPI1),WORK(KCAPIM),NSTR,NSTRT,.TRUE.)
C      ... add MOLFDIR like CAP_IM, WORK(KCAPIM), to the imaginary part of FMOM (if .true.)
        CALL ADDCAPMOM(NSTRT,WORK(KCAPIM),FMOM,.FALSE.)

C      .... get CAP nonre/rl (as a "real" operator) in the MO basis
        CALL GETCAPMO(.TRUE.,WORK(KCAPI),WORK(KCAPI1),WORK,KFREE,LFREE)
C       ... turn it in the MOLFDIR like form
        CALL MOLTODIR(WORK(KCAPI1),WORK(KCAPIM),NSTR,NSTRT,.TRUE.)
C      ... add MOLFDIR like CAP_IM, WORK(KCAPIM), to the imaginary part of FMOM (if .true.)
        CALL ADDCAPMOM(NSTRT,WORK(KCAPIM),FMOM,.TRUE.)

C       ... release all the memory
        CALL MEMREL('CAPMFDADAPT',WORK,KWORK,KWORK,KFREE,LFREE)

C       ... control printout
        IF (IPRCAP.GE.5) THEN
           IF (NZ.EQ.1) THEN
            NZX=2
           ELSE
            NZX=NZ
           ENDIF
           WRITE(LUPRI,'(/3X,A,I2/)')
     &     'CAPMFDADAPT: *** Core Fock matrix with all CAP(re,im) '//
     &     'added (MOLFDIR format) ***'
           CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZX,2),
     &               IPQTOQ(1,0),LUPRI)

        ENDIF
      ENDIF
      CALL QEXIT('CAPMFDADAPT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck savecap */
      SUBROUTINE SAVECAP(LBL,CAPINTS,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Save CAP integrals (CAPINTS) in DIRAC symmetry adapted order 
C  into the AOPROPER file and add corresponding operator to the
C  the list according to the LBL
C 
C  List of permitted labels: 
C 
C             CAP_OVLP - overlap integrals for both relativistic and nonrelativistic mode
C
C             CAP_RE  - real and imaginary integrals for non-relativistic calculations
C             CAP_IM
C
C             CAPD1XR - real and imaginary integrals for relativistic calculations
C             CAPD1XI
C             CAPD1YR
C             CAPD1YI
C             CAPD1ZR
C             CAPD1ZI
C
C  Called from prepcap.F/CAPREAD
C
C  Written by Miro ILIAS, Tel Aviv, Febr./May 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0)
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbxpr.h"
#include "dcbcap.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
      DIMENSION CAPINTS(NNBBASX),WORK(LWORK)
      LOGICAL ISFILE, DOXPRIND1,DOXPRIND2
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      CHARACTER PNAME*16, PLABEL(3)*8
      CHARACTER*8 LBL
      DIMENSION PFAC(3)
      
      CALL QENTER('SAVECAP')
C
C        Write integrals to file
C
      INQUIRE(FILE='AOPROPER',EXIST=ISFILE)

      IF (.NOT.ISFILE) THEN
        WRITE(LUPRI,*) 'SAVECAP: AOPROPER NOT PRESENT !'
        CALL QUIT('SAVECAP: AOPROPER file is NOT PRESENT !')
      ENDIF

#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif

      CALL GETDAT(RTNLBL(1),RTNLBL(2))
C
C     Generate integral labels and write integrals into AOPROPER
C
      IF (LBL.EQ.'CAP_RE  '.OR.
     &    LBL.EQ.'CAP_IM  '.OR.
     &    LBL.EQ.'CAP_OVLP') THEN
        RTNLBL(2)(1:2) = 'SY'
        IF (LEVYLE) THEN
         COMB = 'TFFF'
        ELSE
         COMB = 'TFFT'
        ENDIF
        INTREP = 1
        DOXPRIND1=.TRUE.
        DOXPRIND2=.FALSE.
      ELSE 
        DOXPRIND1=.FALSE.
        DOXPRIND2=.FALSE.

        RTNLBL(2)(1:2) = 'AN'
        COMB = 'FTTF'
             IF (LBL.EQ.'XCAPD1R'.OR. LBL.EQ.'XCAPD1I') THEN
          INTREP = ISYMAX(1,1) + 1
        ELSE IF (LBL.EQ.'YCAPD1R'.OR. LBL.EQ.'YCAPD1I') THEN
          INTREP = ISYMAX(2,1) + 1
        ELSE IF (LBL.EQ.'ZCAPD1R'.OR. LBL.EQ.'ZCAPD1I') THEN
          INTREP = ISYMAX(3,1) + 1
        ELSE
          WRITE(LUPRI,*) 'unidentified SAVECAP label, LBL=',LBL
          CALL QUIT('SAVECAP: Unidentified label !')
        ENDIF
        IF (LBL.EQ.'ZCAPD1I') DOXPRIND2=.TRUE.
      ENDIF

      LABEL = LBL

      WRITE(RTNLBL(2)(3:4),'(I2)') INTREP

      RTNLBL(2)(5:8) = COMB
      CALL WRTPRO(CAPINTS,NNBBASX,LABEL,RTNLBL,IPRCAP)

C     ... control printout ...
      IF (IPRCAP.GE.2) THEN
        WRITE(LUPRI,'(/2X,A,A8,A,A8,A,A)')
     &  'SAVECAP: CAP integrals written to the AOPROPER'//
     &  ' file under the label: ',
     &  RTNLBL(1),' ',RTNLBL(2),' ',LABEL
      ENDIF

      IF (DOXPRIND1) THEN
C
C       CAP (nr_re,nr_im,overlap) integrals 
C     ======================================
C
      PLABEL(1) = LBL
      PNAME = LBL
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &      INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRCAP)

C             ... assign indexes
C ---------------------------------------------------
      IF      (LBL.EQ.'CAP_RE  ') THEN
        IPCAPR    = INDXPR
      ELSE IF (LBL.EQ.'CAP_IM  ') THEN
        IPCAPI    = INDXPR
      ELSE IF (LBL.EQ.'CAP_OVLP') THEN
        IPCAP_S   = INDXPR
      ELSE   
        write(LUPRI,*) 'SAVECAP: LBL=',LBL
        CALL QUIT('SAVECAP: Can not identify INDXPR label !')
      ENDIF
      IF (IPRCAP.GE.2) THEN
        WRITE(LUPRI,'(2X,A,I3)')
     &  'operator added to the list, index=',INDXPR
        CALL WRIXPR(1,INDXPR)
      ENDIF

      ENDIF

C *************************************************
C        ... last two labels ...
C *************************************************
      IF (DOXPRIND2) THEN
C     ... last 2 XPRIND

      PNAME    = 'CAP_D_Re '
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = -CVAL
      PFAC(2)  = -CVAL
      PFAC(3)  = -CVAL
      PLABEL(1) = 'XCAPD1R'
      PLABEL(2) = 'YCAPD1R'
      PLABEL(3) = 'ZCAPD1R'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRCAP)
      IPCAP_DR = INDXPR

      PNAME    = 'CAP_D_Im  '
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = -CVAL
      PFAC(2)  = -CVAL
      PFAC(3)  = -CVAL
      PLABEL(1) = 'XCAPD1I'
      PLABEL(2) = 'YCAPD1I'
      PLABEL(3) = 'ZCAPD1I'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRCAP)
      IPCAP_DI = INDXPR

      IF (IPRCAP.GE.2) THEN
        WRITE(LUPRI,'(2X,A,I3)')
     &  '* CAP_D_Re operator added to the list, index=',IPCAP_DR
        CALL WRIXPR(1,IPCAP_DR)
        WRITE(LUPRI,'(2X,A,I3)')
     &  '* CAP_D_Im operator added to the list, index=',IPCAP_DI
        CALL WRIXPR(1,IPCAP_DI)
      ENDIF

      ENDIF

      CALL QEXIT('SAVECAP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck capready */
      SUBROUTINE CAPREAD(LBL,CAP,CAPSA,CAPSAT,ICAPBFO,IOBF1_2,
     &                   WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   read integrals from the PEMAT file under LBL label into the CAP array
C
C   do the symmetry adaptation into the CAPSA array
C   
C   make triangular matrix of it and save it to the AOPROPER file   
C   
C   Called from prepcap.F/CAPADAPT
C
C Written by Miro ILIAS, May 2007, Tel Aviv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"
#include "dummy.h"
      DIMENSION CAP(N2BBASX),CAPSA(N2BBASX),CAPSAT(NNBBASX),
     &       ICAPBFO(NBASIS),IOBF1_2(NBASIS),WORK(LWORK)
      INTEGER CMNINDX
      LOGICAL FNDLAB
      CHARACTER*8 LBL

      CALL QENTER('CAPREAD')

      CALL DZERO(CAP,N2BBASX)
      CALL DZERO(CAPSA,N2BBASX)
      CALL DZERO(CAPSAT,NNBBASX)

      REWIND(LUPEMAT)
      IF (FNDLAB(LBL,LUPEMAT)) THEN 
        CALL READT(LUPEMAT,NNBBASX,CAPSA)
      ELSE
        write(LUPRI,*) 'CAPREAD: LBL,LUPEMAT:',LBL,LUPEMAT
        CALL QUIT('CAPREAD: Label not found !')
      ENDIF

      IF (IPRCAP.GE.2) THEN
       WRITE(LUPRI,'(/2X,A,A,A)')
     &  'CAPADAPT: CAP integrals of the label ',
     &  LBL,' succesfully read from PEMAT file into the CAPSA array.'
      ENDIF

C==================================================================
C       ... in the cycle resort the matrix elements
C    from CAPSA PEMAT order to the DIRAC basic order into WORK arrays
C==================================================================
      IX = 0
      DO I = 1, NBASIS
      DO J = 1, I
        IX = IX + 1
        IX1 = CMNINDX(I,J)
        IF (IX.NE.IX1) THEN
         write(LUPRI,*) 'I,J: ',IX1,IX
         write(LUPRI,*) 'common indexes ',IX1,IX
         CALL QUIT('CAPADAPT: wrong IX,IX1 !')
        ENDIF

C       -------------------------------------------------
C       .. store  integrals in the DIRAC basic order ..
C       -------------------------------------------------

        indx_i =  NBASIS*(ICAPBFO(I)-1) + ICAPBFO(J)
        indx_j =  NBASIS*(ICAPBFO(J)-1) + ICAPBFO(I)

C       CAP_FACT = GET_CAPAF_FAC(IOBF1(2,I))*GET_CAPAF_FAC(IOBF1(2,J))
        CAP_FACT = GET_CAPAF_FAC(IOBF1_2(I))*GET_CAPAF_FAC(IOBF1_2(J))

        CAP(indx_i) = CAPSA(IX)*CAP_FACT
        CAP(indx_j) = CAPSA(IX)*CAP_FACT

      ENDDO 
      ENDDO 

C   ... delete SL/LS blocks integrals ....
      IF (LBL.EQ.'CAP_RE  '.OR.
     &    LBL.EQ.'CAP_IM  '.OR.
     &    LBL.EQ.'CAP_OVLP') THEN
        CALL DELETE_LS(CAP,NBASIS,IPRCAP,.TRUE.)
      ELSE
        CALL DELETE_LS(CAP,NBASIS,IPRCAP,.FALSE.)
      ENDIF

      IF (IPRCAP.GE.4) THEN
       CALL HEADER ('CAPxx integrals in basic DIRAC order'//
     & ' (cleared blocks) ',-1)
       WRITE(LUPRI,'(2X,A)') 'The label is:',LBL
       CALL PRQMAT(CAP,NBASIS,NBASIS,NBASIS,NBASIS,1,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

C   ...internal DIRAC symmetry adaptation into CAPSA
      IF (LBL.EQ.'CAP_RE  '.OR.
     &    LBL.EQ.'CAP_IM  '.OR.
     &    LBL.EQ.'CAP_OVLP') THEN
        IREP=0
C       IREPD = IRQMAT(1,IREP)
C       IREPD = IRQMAT(JM4POS(1),IREP)
        IREPD = 0
      ELSE IF (LBL.EQ.'XCAPD1R '.OR.LBL.EQ.'XCAPD1I ') THEN
        IREP = ISYMAX(1,1)
        IREPD = IRQMAT(JM4POS(4),IREP) 
      ELSE IF (LBL.EQ.'YCAPD1R '.OR.LBL.EQ.'YCAPD1I ') THEN
        IREP = ISYMAX(2,1)
        IREPD = IRQMAT(JM4POS(3),IREP)
      ELSE IF (LBL.EQ.'ZCAPD1R '.OR.LBL.EQ.'ZCAPD1I ') THEN
        IREP = ISYMAX(3,1)
        IREPD = IRQMAT(JM4POS(2),IREP)
      ELSE 
        WRITE(LUPRI,*) 'CAPREAD ... LBL=',LBL
        CALL QUIT('CAPREAD: wrong LBL !')
      ENDIF

      CALL DZERO(CAPSA,N2BBASX)

      CALL MTAOSO(CAP,CAPSA,NBASIS,IREPD,IPRCAP)

C   ... make triangular matrix of CAPSA into CAPSAT
      CALL DAMATR(NBASIS,CAPSA,NBASIS,CAPSAT)

C   ... save CAPSAT of LBL to the AOPROPER file
      CALL SAVECAP(LBL,CAPSAT,WORK,LWORK)

      CALL QEXIT('CAPREAD')
      RETURN
      END
 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck cmnindx */
      INTEGER FUNCTION CMNINDX(I,J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Common index for lower triangular element (I >= J)
C
C   Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      IM = MAX(I,J)
      JM = MIN(I,J)
      CMNINDX = (IM*(IM-1)/2) + JM
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck iobfset */
      SUBROUTINE IOBFSET(IOBF,II,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Find out unique indexes related to atomic basis functions and
C their centers.
C
C  On input: II - actual position
C
C  On output: IOBF - filled 3rd parameter
C
C Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION IOBF(1:3,*)
      CALL QENTER('IOBFSET')

      IOBF(3,II) = 0
      INUCL = IOBF(1,II) 
      IAOBF = IOBF(2,II)
      DO I = II, 1, -1
       IF (INUCL.EQ.IOBF(1,I).AND.IAOBF.EQ.IOBF(2,I)) THEN
          IOBF(3,II) = IOBF(3,II) + 1
       ENDIF 
      ENDDO

      IF (IPRINT.GE.17) THEN
       WRITE(LUPRI,*) 'IOBFSET: for ',II,
     & ' DIRAC function unique indx=',IOBF(3,II)
      ENDIF

      CALL QEXIT('IOBFSET')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck aonum */
      INTEGER FUNCTION AONUM(AO)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Return the unique number of the atomic basis 
C  function in variable AO.
C
C  Functions are only up to F.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      CHARACTER AO*(*)

      IF (AO.EQ.'S'.OR.AO.EQ.'s') THEN 
        AONUM = 1
      ELSE IF (AO.EQ.'PX'.OR.AO.EQ.'px') THEN
        AONUM = 2
      ELSE IF (AO.EQ.'PY'.OR.AO.EQ.'py') THEN
        AONUM = 3
      ELSE IF (AO.EQ.'PZ'.OR.AO.EQ.'pz') THEN
        AONUM = 4
      ELSE IF (AO.EQ.'XX'.OR.AO.EQ.'dxx') THEN
        AONUM = 5
      ELSE IF (AO.EQ.'XY'.OR.AO.EQ.'dxy') THEN
        AONUM = 6
      ELSE IF (AO.EQ.'XZ'.OR.AO.EQ.'dxz') THEN
        AONUM = 7
      ELSE IF (AO.EQ.'YY'.OR.AO.EQ.'dyy') THEN
        AONUM = 8
      ELSE IF (AO.EQ.'YZ'.OR.AO.EQ.'dyz') THEN
        AONUM = 9
      ELSE IF (AO.EQ.'ZZ'.OR.AO.EQ.'dzz') THEN
        AONUM = 10
      ELSE IF (AO.EQ.'XXX'.OR.AO.EQ.'fxxx') THEN
        AONUM = 11
      ELSE IF (AO.EQ.'XXY'.OR.AO.EQ.'fxxy') THEN
        AONUM = 12
      ELSE IF (AO.EQ.'XXZ'.OR.AO.EQ.'fxxz') THEN
        AONUM = 13
      ELSE IF (AO.EQ.'YYX'.OR.AO.EQ.'fxyy') THEN
        AONUM = 14
      ELSE IF (AO.EQ.'XYZ'.OR.AO.EQ.'fxyz') THEN
        AONUM = 15
      ELSE IF (AO.EQ.'ZZX'.OR.AO.EQ.'fxzz') THEN
        AONUM = 16
      ELSE IF (AO.EQ.'YYY'.OR.AO.EQ.'fyyy') THEN
        AONUM = 17
      ELSE IF (AO.EQ.'YYZ'.OR.AO.EQ.'fyyz') THEN
        AONUM = 18
      ELSE IF (AO.EQ.'ZZY'.OR.AO.EQ.'fyzz') THEN
        AONUM = 19
      ELSE IF (AO.EQ.'ZZZ'.OR.AO.EQ.'fzzz') THEN
        AONUM = 20
      ELSE
       write(LUPRI,*) '>>> fu AONUM: entering atomic orbital:',AO,
     & ' is not in the list (which is up only to F functions) !'       
       CALL QUIT('fu AONUM: entering AO does not fit; basis'//
     & ' set up to f functions only !')
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION GET_CAPAF_FAC(I)
#include "implicit.h"
#include "priunit.h"
         IF (I.EQ.5.OR.I.EQ.8.OR.I.EQ.10) THEN
C     ... d xx,yy,zz - functions
            GET_CAPAF_FAC =  DSQRT(3.0D0)
C     ... f - functions: XXX, YYY, ZZZ
         ELSE IF (I.EQ.11.OR.I.EQ.17.OR.I.EQ.20) THEN
C        ELSE IF (I.GE.11.OR.I.LE.20) THEN
C           GET_CAPAF_FAC =  DSQRT(3)
C           GET_CAPAF_FAC =  1.0D0
            GET_CAPAF_FAC = DSQRT(15.0D0)
C        XXY, XXZ, YYX, ZZX, YYZ, ZZY,              
         ELSE IF (I.EQ.12.OR.I.EQ.13.OR.I.EQ.14.OR.
     &            I.EQ.16.OR.I.EQ.18.OR.I.EQ.19) THEN
            GET_CAPAF_FAC =  DSQRT(3.0D0)
C        XYZ
         ELSE IF (I.EQ.15) THEN
            GET_CAPAF_FAC =  1.0D0
         ELSE
            GET_CAPAF_FAC = 1.0D0
         ENDIF

      RETURN
      END

