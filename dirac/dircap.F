C
C
C  File including Sajeev's program for generation
C of the CAP one-electron integrals, plus interface
C routines.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capinp */
      SUBROUTINE CAPINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C   Input section for CAP-module
C
C  Written by  Miro ILIAS, Febr. 2007, Tel Aviv
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 20)
      PARAMETER (D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF,RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)

C#include "mxcent.h"
C#include "dcbgen.h"
C#include "dcbpsi.h"
#include "dcbham.h"
#include "dcbcap.h"
C
      SAVE SET
      DATA TABLE /'.INPQD ','.INPCAP','.PRINT ','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
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
C&&& INPQD: read the PW number
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
                 CALL QUIT('CAPINT: Increase dimension parameter !')
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
C&&& XXXXXX
                 GO TO 100
 5               CONTINUE
C&&& XXXXXX
                 GO TO 100
 6               CONTINUE
C&&& XXXXXX
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
     &   'CAPINP: Set-up for CAP calculation'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)') ' * Print level: ',IPRCAP
C
 90   CONTINUE
      WRITE(LUPRI,'(A)')
     &     ' * CAP integrals  according to',
     &     '          Sajeev ....         '
C
       WRITE(LUPRI,'(2X,A,D25.10)') 'read INPQD = ',PW
       WRITE(LUPRI,'(2X,A)') 
     & 'read integration parameters'
CMI    ...  write out integration parameters
       DO I=1,3
         WRITE(LUPRI,*) X0(I),LAMDA(I)
         WRITE(LUPRI,*) NGR(I)
         WRITE(LUPRI,*) (GINTPT(I,J),J=1,NGR(I))
         WRITE(LUPRI,*) (NGRS(I,J),J=1,NGR(I)-1)
       ENDDO

  999 CONTINUE
C
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capints */
      SUBROUTINE CAPINTS(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Interface for the CAP integral program.
C
C   Called from PAMINP.
C
C   Written by Miro ILIAS, Tel Aviv, January 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
C#include "dcbham.h"

      DIMENSION WORK(LWORK)
#include "memint.h"

      CALL QENTER('CAPINTS')

C      ... allocate memmory
      ISIZE1=10
      ISIZE1=20
      CALL MEMGET('INTE',KI,ISIZE1,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KF,ISIZE2,WORK,KFREE,LFREE)

      CALL SAJEEV(WORK,LWORK)

      CALL MEMREL('CAPINTS',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL CAPADAPT(WORK,LWORK)
            
      CALL QEXIT('CAPINTS')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      LOGICAL FUNCTION CHPRES(IPRHAM)
#include "implicit.h"
#include "priunit.h"

      LOGICAL FILES,FILES2

C  ... check presence of files CAPLBL.inp PEMAT CAPDIRSYM.inp CAPDIR.inp
      FILES = .FALSE.
      INQUIRE(FILE='CAPLBL.inp',EXIST=FILES)
      FILES2 = FILES

      FILES = .FALSE.
      INQUIRE(FILE='PEMAT',EXIST=FILES)
      FILES2 = FILES2.AND.FILES

      FILES = .FALSE.
      INQUIRE(FILE='CAPDIRSYM1.inp',EXIST=FILES)
      FILES2 = FILES2.AND.FILES

      FILES = .FALSE.
      INQUIRE(FILE='CAPDIRSYM2.inp',EXIST=FILES)
      FILES2 = FILES2.AND.FILES

      FILES = .FALSE.
      INQUIRE(FILE='CAPDIR.inp',EXIST=FILES)
      FILES2 = FILES2.AND.FILES

      CHPRES = FILES2

      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capadapt */
      SUBROUTINE CAPADAPT(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Adapt the CAP integrals from the formatted file PEMAT for the
C  DIRAC calculations.
C
C   Written by Miro ILIAS, Tel Aviv, 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER (MXABF=2000)
#include "dcbbas.h"
#include "dcbcap.h"
#include "mxcent.h"
C#include "maxorb.h"
#include "nuclei.h"
#include "dummy.h"

      DIMENSION WORK(LWORK)
      DIMENSION IOBF1(1:3,1:MXABF),IOBF2(1:3,1:MXABF),
     &          ICAPBFO(1:MXABF)
      DIMENSION ILC1(1:MXCENT),ILC2(1:MXCENT)
      INTEGER AONUM,CMNINDX
      CHARACTER*6 AO1,AO2
      LOGICAL CHPRES
      COMPLEX CAP
#include "memint.h"

      CALL QENTER('CAPADAPT')
 
      IF (NBASIS.GT.MXABF) THEN
       WRITE(LUPRI,*)
     & '>>> NBASIS=',NBASIS,' parameter MXABF=',MXABF 
       CALL QUIT('CAPADAPT: increase MXABF parameter !')
      ENDIF

      IF (.NOT.CHPRES(IPRCAP)) THEN
         CALL QUIT('...not all files for CAP !')
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
C     .. skip text lines
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
          write(LUPRI,'(I4,A,I3,A,A4,A,I3)')
     &    II1,'> nucleus:',INC1,' funct: ',AO1,
     &     '>  CAP unique fu #:',
     &        IOBF1(3,I)
          write(LUPRI,'(I4,A,I3,A,A4,A,I3)')
     &    II2,'> nucleus:',INC2,' funct: ',AO2,
     &     '>dirac unique fu #:',
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
          IF (IOBF2(1,I).EQ.IOBF1(1,J).AND.
     &        IOBF2(2,I).EQ.IOBF1(2,J).AND.
     &        IOBF2(3,I).EQ.IOBF1(3,J)) IX = J        
        ENDDO

        IF (IX.EQ.0) THEN
          WRITE(LUPRI,*) '  '
          WRITE(LUPRI,*) '  '
          CALL QUIT('CAPADAPT: order does not fit !!')
        ENDIF
 
        ICAPBFO(I) = IX

C       IF (IPRCAP.GE.2) THEN
C        WRITE(LUPRI,'(I4,A,I4)')
C    &   I,'. CAP bf is DIRAC unique fu:',ICAPBFO(I)
C       ENDIF
      ENDDO

      IF (IPRCAP.GE.2) THEN
       WRITE(LUPRI,'(2X,A,I4)')
     & 'Matching between CAP and DIRAC atomic basis'//
     & ' functions - basic order. NBASIS=',NBASIS
       DO I = 1, NBASIS
         WRITE(LUPRI,'(1X,I3,A,I3,A,$)') I,'->',ICAPBFO(I),'|'
         IF (MOD(I,8).EQ.0) WRITE(LUPRI,*)
       ENDDO
       WRITE(LUPRI,*)
      ENDIF


      LUPM = 0
      CALL GPOPEN(LUPM,'PEMAT','OLD',' ','FORMATTED',
     &       IDUMMY,.FALSE.)

      ISIZE = NBASIS*(NBASIS+1)/2
      IF (ISIZE.NE.NNBBASX) THEN
       WRITE(LUPRI,*) "ISIZE=",ISIZE," NNBBASX=",NNBBASX
       CALL QUIT('CAPADAPT: ISIZE <> NNBBASX !')
      ENDIF
      
      CALL MEMGET('REAL',KPMR,ISIZE,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMI,ISIZE,WORK,KFREE,LFREE)
C==================================================================
C     ... read PEMAT in cycle and resort the matrix elements
C    from CAP order to the DIRAC basic order into WORK arrays
C==================================================================
      IX = 0
      DO I = 1, NBASIS
      DO J = 1, I
        IX = IX + 1
        READ(LUPM,*) II,JJ,CAP
        CAPRE = DREAL(CAP)
CC      CAPIM = DIMAG(CAP)
        CAPIM = IMAG(CAP)
        IF (I.NE.II.OR.J.NE.JJ) THEN
         CALL QUIT('CAPADAPT: wrong indexes in PEMAT !')
        ENDIF
C       IX1 = (I*(I-1)/2) + J
        IX1 = CMNINDX(I,J)
        IF (IX.NE.IX1) THEN
         CALL QUIT('CAPADAPT: wrong IX !')
        ENDIF
C       write(LUPRI,*) II,JJ,CAP,' ',CAPRE,CAPIM
C       write(LUPRI,*) II,JJ,CAPRE,CAPIM

C       -------------------------------------------------
C       .. store the integrals in the DIRAC 'C1' order ..
C       -------------------------------------------------
        WORK(KPMR+CMNINDX(ICAPBFO(I),ICAPBFO(J))-1)=CAPRE
        WORK(KPMI+CMNINDX(ICAPBFO(I),ICAPBFO(J))-1)=CAPIM

      ENDDO 
      ENDDO 
      CALL GPCLOSE(LUPM,'KEEP')

C===========================================================
C
C Do the (HERMIT) symmetry adaptation of unsorted integrals
C
C===========================================================
      CALL MEMGET('REAL',KPMRSA,ISIZE,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMISA,ISIZE,WORK,KFREE,LFREE)

      IF (IPRCAP.GE.4) THEN
       WRITE(LUPRI,'(/,5X,A/)') 'CAPADAPT: Symmetry adaptation:'
      ENDIF

      LUCDS1 = 0
      CALL GPOPEN(LUCDS1,'CAPDIRSYM1.inp','OLD',' ','FORMATTED',
     &             IDUMMY,.FALSE.)
      LUCDS2 = 0
      CALL GPOPEN(LUCDS2,'CAPDIRSYM2.inp','OLD',' ','FORMATTED',
     &             IDUMMY,.FALSE.)
      READ(LUCDS1,*)
      DO I = 1, NBASIS
        READ(LUCDS1,*) II1,NLC1
        READ(LUCDS1,*) (ILC1(JJ1),JJ1=1,NLC1)

        IF (IPRCAP.GE.4) THEN
          write(LUPRI,*) 'I-b.fu:',II1,' no lin.com:',
     &    NLC1,':',(ILC1(JJ1),JJ1=1,NLC1)
        ENDIF

      READ(LUCDS2,*)
      DO J = 1, I
        READ(LUCDS2,*) II2,NLC2
        READ(LUCDS2,*) (ILC2(JJ2),JJ2=1,NLC2)
 
        IF (IPRCAP.GE.4) THEN
          write(LUPRI,*) 'J-b.fu:',II2,' no of lin.com:',
     &    NLC2,':',(ILC2(JJ2),JJ2=1,NLC2)
        ENDIF

C      ... symmetry adapt CAP_(I,J):
        CAPSARE = 0.0D0
        CAPSAIM = 0.0D0
        DO K =1, NLC1
        DO L =1, NLC2

           COMFACT = 1/DFLOAT(NLC1*NLC2)

           IF (ILC1(K).GT.0) THEN
             SGNFAC1 = +1
           ELSE
             SGNFAC1 = -1
           ENDIF

           IF (ILC2(L).GT.0) THEN
             SGNFAC2 = +1
           ELSE
             SGNFAC2 = -1
           ENDIF

           TOTFAC = COMFACT*SGNFAC1*SGNFAC2

           CAPRE = WORK(KPMR+CMNINDX(ABS(ILC1(K)),ABS(ILC2(L)))-1)
           CAPIM = WORK(KPMI+CMNINDX(ABS(ILC1(K)),ABS(ILC2(L)))-1)

           CAPSARE = CAPSARE + TOTFAC*CAPRE
           CAPSAIM = CAPSAIM + TOTFAC*CAPIM

        ENDDO        
        ENDDO        
        WORK(KPMRSA+CMNINDX(I,J)-1)=CAPSARE
        WORK(KPMISA+CMNINDX(I,J)-1)=CAPSAIM
      ENDDO
      REWIND(LUCDS2)
      ENDDO
      CALL GPCLOSE(LUCDS1,'KEEP')
      CALL GPCLOSE(LUCDS2,'KEEP')

C =======================================================================
C     Save HERMIT symmetry adapted integrals into the AOPROPER file.
C =======================================================================
      CALL SAVECAP(WORK(KPMRSA),WORK(KPMISA),IPRCAP)

      CALL CAPREADY(WORK,LWORK,IPRCAP)

      CALL MEMREL('CAPADAPT',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('CAPADAPT')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck getcapim */
      SUBROUTINE GETCAPI(CAPI,WORK,KFREE,LFREE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Purpose: get CAP_IM operator in MO basis into CAPI(*) 
C
C   First get it (using PRPMAO) as the real operator
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
#include "dummy.h"

      DIMENSION CAPI(NORBT,NORBT,NZ),WORK(*)
      LOGICAL SWITCH

      CALL QENTER('GETCAPI')

      KFRSAV = KFREE

      IPRINT = IPRHAM

C  ...array for MO coefficients
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)

C boson irreps identification (read with MO coeff.)
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
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
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','GETCAPIM')

      INDXPR = IPCAPI
      ISYM  = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      IOPSY = JBTOF(ISYM-1,1)
      CALL DZERO(CAPI,N2ORBXQ)
      CALL PRPMAT(INDXPR,IOPSY,CAPI,WORK(KCMO),WORK(KBEIG),
     &                  ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(2A)')
     & 'GETCAPI: Get MO matrix of property: ',PRPNAM(INDXPR)
       CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

C     ... make the real operator in MO basis "i"-operator
      SWITCH=.FALSE.
      IF (SWITCH) THEN
       IF (NCMOTQ.LT.NORBT*NORBT*NZ) THEN
         CALL QUIT('GETCAPIM: fix the size !!!')
       ENDIF
CMI    switch in MO-basis is NOT theoretically established !!!!
      IF (NZ.GT.1) THEN
C      ... switch between IZ=1 and IZ=2
       DO I = 1, NORBT
       DO J = 1, NORBT
        TEMP = CAPI(I,J,1)
        CAPI(I,J,1) = -CAPI(I,J,2)    
        CAPI(I,J,2) = TEMP
       ENDDO
       ENDDO

C      ... switch between IZ=3 and IZ=4
       IF (NZ.EQ.4) THEN
        SF = D1
CMI     SF = DM1
        DO I = 1, NORBT
        DO J = 1, NORBT
          TEMP = CAPI(I,J,3)
          CAPI(I,J,3) = SF*CAPI(I,J,4)
          CAPI(I,J,4) = -SF*TEMP
        ENDDO
        ENDDO
       ENDIF

       WRITE(LUPRI,'(2A)')
     & 'GETCAPIM: "i"-changed MO matrix of property: ',PRPNAM(INDXPR)
       CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

      ENDIF

      CALL MEMREL('MRCON1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
CC
CC    call quit('GETCAPIM ...')

      CALL QEXIT('GETCAPI')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck addcapim */
      SUBROUTINE ADDCAPIM(NSTRT,CAPIM,FMOM,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Purpose: add MOLFDIR like i.CAP_IM to the FMOM
C
C  Written by Miro ILIAS, Tel Aviv, febr. 2007 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"

      DIMENSION CAPIM(NSTRT*2,NSTRT*2,2), FMOM(NSTRT*2,NSTRT*2,2)

      CALL QENTER('ADDCAPIM')

      IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Entering core Fock matrix in MOLFDIR'//
     &    ' format, FMOM ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)

          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Entering (Hamiltonian like) CAPIM ***'
          CALL PRQMAT(CAPIM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF

C     ... multiply CAPIM with i:

      DO I = 1, NSTRT*2
      DO J = 1, NSTRT*2
       TEMP = CAPIM(I,J,2)
       CAPIM(I,J,2) = CAPIM(I,J,1)
       CAPIM(I,J,1) = - TEMP
C      FMOM(I,J,2) = CAPIM(I,J,1)
      ENDDO
      ENDDO

      IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Modified i.CAPIM ***'
          CALL PRQMAT(CAPIM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF

C      ...Add CAPIM into FMOM
      DO I = 1, NSTRT*2
      DO J = 1, NSTRT*2
       FMOM(I,J,1) = FMOM(I,J,1) + CAPIM(I,J,1)
       FMOM(I,J,2) = FMOM(I,J,2) + CAPIM(I,J,2)
      ENDDO
      ENDDO


      IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Core Fock matrix with CAP_IM (NZ=2) ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &                IPQTOQ(1,0),LUPRI)
      ENDIF

      CALL QEXIT('ADDCAPIM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck savecap */
      SUBROUTINE SAVECAP(CAPR,CAPI,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Save complex - CAPR,CAPI - integrals into AOPROPER
C
C Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbgen.h"
      DIMENSION CAPR(*),CAPI(*)
      LOGICAL ISFILE
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      
      CALL QENTER('SAVECAP')

      IF (IPRINT.GE.3) THEN
        CALL HEADER('Output from SAVECAP',-1)
      ENDIF

C     write(LUPRI,*) "ISIZE=",ISIZE," NNBBASX=",NNBBASX

C        Write integrals to file
C
      INQUIRE(FILE='AOPROPER',EXIST=ISFILE)

      IF (.NOT.ISFILE) THEN
        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,*) 'SAVECAP: AOPROPER NOT PRESENT !'
        ENDIF
C       CALL QUIT('SAVECAP: AOPROPER NOT PRESENT !')
      ENDIF

#if defined (SYS_IBM)
C        OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
C        OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif

C
C  Generate integral labels and write integrals into AOPROPER
C
      COMB = 'TFFT'
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
      LABEL = 'CAP_RE  '
      CALL WRTPRO(CAPR,NNBBASX,LABEL,RTNLBL,IPRINT)

      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      COMB = 'TFFT'
C     COMB = 'TTTT'
      RTNLBL(2)(5:8) = COMB
      LABEL = 'CAP_IM  '
      CALL WRTPRO(CAPI,NNBBASX,LABEL,RTNLBL,IPRINT)

      IF (IPRINT.GE.1) THEN
       WRITE(LUPRI,'(2X,A)')
     & 'CAP_RE and CAP_IM integrals written into AOPROPER file.'
      ENDIF

      CALL QEXIT('SAVECAP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck capready */
      SUBROUTINE CAPREADY(WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Add cap operators into integral list.
C
C   Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.0D0,DM2 = -2.0D0,DHALF=0.5D0)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbxpr.h"
#include "dcbham.h"
      CHARACTER PNAME*16, PLABEL(3)*8
      DIMENSION PFAC(3),WORK(LWORK)

      CALL QENTER('CAPREADY')

      IF (IPRINT.GE.7) THEN
        CALL HEADER('Output from CAPREADY',-1)
      ENDIF
C
C       CAP (complex) integrals
C     ============================
C
      PNAME            = 'CAP - real part '
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'CAP_RE  '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPCAPR           = INDXPR


      PNAME            = 'CAP - imaginary part '
C     PNAME            = 'CAP_IM  '
C     IPTYP            = 20
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      PLABEL(1)        = 'CAP_IM  '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)
      IPCAPI           = INDXPR

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CMI     keep rather as real symmetric...
C      IPRPTIM(IPCAPI) = -1
C      write(LUPRI,*) '...CAP_IM IPRPTIM=',IPRPTIM(IPCAPI)
C      CALL PRPSYM(IPCAPI,IPRINT+10)
C      CALL WRIXPR(0,IPCAPI)
C      write(LUPRI,*) 'CAP_IM INDXPR=',IPCAPI
C      write(LUPRI,*) 'CAP_IM IPRPSYM=',IPRPSYM(IPCAPI)
C      write(LUPRI,*) 'CAP_IM IPRPTIM=',IPRPTIM(IPCAPI)
CMI    ... once again
C      CALL PRPSYM(IPCAPI,IPRINT+10)

      CALL QEXIT('CAPREADY')
      RETURN
      END
 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck cmnindx */
      INTEGER FUNCTION CMNINDX(I,J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Common index for lower triangular element (I >= J)
C
C Written by Miro ILIAS, Tel Aviv, febr. 2007
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
C  function in variable AO
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
       CALL QUIT('fu AONUM: entering AO does not fit !!')
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sajeev */
      SUBROUTINE SAJEEV(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Main program from Sajeev. Adapted into DIRAC enviroment
C  as the subroutine.
C
C  On input: IPRINT - print level
C
C  called from CAPINTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! this file  reads  all input datas for CAP and generates  
! CAP 1-e integrals in MO basis                                        
! Sajeev.Y.       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CMI   IMPLICIT REAL*8 (A-H,O-Z)
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
      
      DIMENSION WORK(LWORK)

CMI    ... make allocations in the calling routine !!
      integer *4 basis_n(200)
      integer *4 ang_pow(250,3)

      real *8 ord_expon(250),ord_coe(250)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(200,3)

C     real *8 ZEL_AO(25000)
C     real *8 KE_AO(25000)

CMI   ... this variable is beeing sent into AOINTEGRAL
CMI   ... split it into real and imaginary parts !!!
      COMPLEX *16 PE_AO(25000)
C     real *8 AOMO(200,200)
      real *8 BETA0,PE0
      
      LOGICAL ISFILE

CMI   ... are all these allocations needed ???
CMI   complex *16 EL1_AO(25000),EL1_MO(25000),EL1_CI(10000000)
C     Real *8 TWOE(10000000)
C     complex *16 HAM_CI(10000000)
      COMPLEX *16 THETA,THETA2,CI,BETA
      data CI/(0.d0,-1.d0)/

      CALL QENTER('SAJEEV')

CMI    ... no need of MO coefficients from GAMESS now, they are read later in moltra
C     open(unit=3,file='AOMOCOE',STATUS='UNKNOWN',FORM='FORMATTED')

      IF (IPRCAP.GE.2) THEN
        CALL HEADER('*** Output from SAJEEV ***',-1)
      ENDIF

CMI   ... this is for the "theta" parameter
C     open(unit=8,file='INPQD',STATUS='OLD',FORM='FORMATTED')
C     READ(8,*)PW
C     close(8,status="keep")
      IF (IPRCAP.GE.3) THEN
        write(LUPRI,*) ' read INPQD (theta) =',PW
      ENDIF

      CALL INPREAD(N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,IPRCAP)

C      read(3,*)N_B
C      do i=1,N_B
C         DO j=1,N_B
C            READ(3,*)k,l,AOMO(i,j)
C         END DO
C      END DO

CMI ... check the presence of PEMAT in order to decide whether to run 
CMI   or not the CAP integral evaluation
      INQUIRE(FILE='PEMAT',EXIST=ISFILE)
      IF (.NOT.ISFILE) THEN
CMI   ... evaluate CAP complex atomic untegrals
        IF (IPRCAP.GE.0) THEN
         WRITE(LUPRI,'(2X,A)')
     &   'PEMAT files in not present, going to evaluate CAP integrals.'
        ENDIF
        CALL AOINTEGRAL(N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     *                 PW,PE_AO,IPRCAP)
        
CMI   ... write out atomic integrals into the file
C     CALL MOINTEGRAL(N_B,AOMO,PE_AO,EL1_MO)
         open(unit=14,file='PEMAT',STATUS='UNKNOWN',FORM='FORMATTED')
         do ii=1,n_b
            do jj=1,n_b
               IPU=IPOINT(ii,jj)
               if (jj.gt.ii)go to 301
               write(14,*)ii,jj,PE_AO(IPU)
            end do
 301        continue
         end do
         close(14)
         IF (IPRCAP.GE.0) THEN
          WRITE(LUPRI,'(2X,A)')
     &    'PEMAT file with CAP integrals generated.'
        ENDIF
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'SAJEEV: PEMAT file is present, not calculating CAP integrals.'
      ENDIF
C     call quit(' end in sajeev...')

      CALL QEXIT('SAJEEV')
      RETURN
      END
      
C***************************
C ... read the input ...
C***************************      
      SUBROUTINE INPREAD(
     & N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,IPRINT)
C ...
CMI   IMPLICIT REAL*8 (A-H,O-Z)
#include "implicit.h"
#include "priunit.h"
C#include "dummy.h"
      
      CHARACTER *1 shell_type,BLNK(5)
      CHARACTER *2 shell_p(3),shell_d(6),Shell_arr(200)   
      CHARACTER *3 shell_f(10)   
      CHARACTER *8 IA

CMI ... NEW ...basis up to "F"
      character *3 Shell_basis(200)

      integer *4 Shell_atom(200)
      
      integer *4 n_con_guass,nj,g_n,sh,n_c_g
      integer *4 angu_pow_S(1,3),angu_pow_P(3,3)
      integer *4 angu_pow_D(6,3),angu_pow_F(10,3)
      integer *4 n_b,ia_p
      integer *4 ang_pow(250,3)
      integer *4 n_pg_atom(50)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(200)
      
      real *8 re(50,3),expon,coe
      real *8 ord_expon(250),ord_coe(250)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(200,3)

      DATA BLNK /'E','S','P','D','F'/
      DATA shell_p /'PX','PY','PZ'/
      DATA shell_d /'XX','YY','ZZ','XY','XZ','YZ'/
      DATA shell_f /'XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ',
     *     'ZZX','ZZY','XYZ'/
      DATA angu_pow_S /0,0,0/
      DATA angu_pow_P /1,0,0,0,1,0,0,0,1/
      DATA angu_pow_D /2,0,0,1,1,0,
     *     0,2,0,1,0,1,
     *     0,0,2,0,1,1/
      DATA angu_pow_F /3,0,0,2,2,1,0,1,0,1,
     *     0,3,0,1,0,2,2,0,1,1,
     *     0,0,3,0,1,0,1,2,2,1/

#include "dummy.h"
      CALL QENTER('INPREAD')

C     CALL HEADER('*** Output from SAJEEV:INPREAD ***',-1)

      do i=1,200
        Basis_N(i)=0
      End do

      LUinp = 57
      CLOSE(LUinp)

C     OPEN(unit=5,file='CAPINT.inp',STATUS='OLD',FORM='FORMATTED')
C     OPEN(unit=LUinp,file='CAPINT.inp',STATUS='OLD',FORM='FORMATTED')

C     LUinp = 0
C     CALL GPOPEN(LUinp,'CAPINT.inp','OLD',' ','FORMATTED',
C    &            IDUMMY,.FALSE.)

      OPEN(unit=LUinp,file='CAPINT.inp',STATUS='OLD',FORM='FORMATTED')
CC    OPEN(unit=LUinp,file='CAPINT.inp',STATUS='UNKNOWN',
CC   &      FORM='FORMATTED')

C     READ(5,*) n_atoms
      READ(LUinp,*) n_atoms

      IF (IPRINT.GE.3) THEN
       write(lupri,*) 'INPREAD: n_atoms= ',n_atoms
      ENDIF

      g_n=1
      sh=0
      iiii=1

      DO  k=1,n_atoms

         READ(LUinp,*) IA

      IF (IPRINT.GE.3) THEN
         write(lupri,*) 'label  IA=',IA
      ENDIF

CMI      write(19,*)IA

         READ(LUinp,*)(re(k,i),i=1,3)

      IF (IPRINT.GE.3) THEN
         write(lupri,*) '> coordinate: ',(re(k,i),i=1,3)
      ENDIF
          
         Do ISNEW=1,1000  

           READ(LUinp,*) shell_type, n_con_guass 
           
           IF (IPRINT.GE.3) THEN
            write(lupri,*) 'shell_type=',shell_type
            write(lupri,*) 'n_con_guass=',n_con_guass
           ENDIF

CMI       'E', means end
            if(shell_type.eq.BLNK(1))go to 201
CMI       'S' shell
            if(shell_type.eq.blnk(2)) go to 200
CMI       'P' shell
            if(shell_type.eq.blnk(3)) go to 300
CMI       'D' shell
            if(shell_type.eq.blnk(4)) go to 400
CMI       'F' shell
            if(shell_type.eq.blnk(5)) go to 500
              
 200        continue
          
            DO n_c_g=1,n_con_guass
               sh=sh+1

               READ(LUinp,*)nj,expon,coe

          IF (IPRINT.GE.3) THEN
            write(lupri,*)'S: nj=',nj,' expon=',expon,' contr.coef=',coe
          ENDIF

               ord_expon(sh)=expon
               ord_coe(sh)=coe
               shell_arr(sh)=shell_type
CMI NEW 
               Shell_basis(g_n)=shell_type
               Shell_atom(g_n)=k

               do ia_p=1,3
                  ang_pow(sh,ia_p)=angu_pow_S(1,ia_p)
               end do
               basis_n(g_n)=n_con_guass

            END DO
              g_n=g_n+1
              go to 600
             
 300         continue
             
             DO n_c_g=1,n_con_guass
                sh=sh+1
                READ(LUinp,*)nj,expon,coe

           IF (IPRINT.GE.3) THEN
          write(lupri,*)'P: nj=',nj,' expon=',expon,' contr.coef=',coe
           ENDIF

                do i=1,3
                   shell_arr(sh+(i-1)*n_con_guass)=shell_p(i)
CMI NEW
                   Shell_basis(g_n+i-1)=shell_p(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe
                   do ia_p=1,3
                   ang_pow(sh+(i-1)*n_con_guass,ia_p)=angu_pow_p(i,ia_p)
                   end do
                end do                              
             END DO
          
          do  iji=1,3
             basis_n(g_n)=basis_n(g_n)+n_con_guass
             g_n=g_n+1
             end do
             sh=sh+2*n_con_guass
             go to 600
             
 400         continue
             
             DO   n_c_g=1,n_con_guass
                sh=sh+1
                READ(LUinp,*)nj,expon,coe

           IF (IPRINT.GE.3) THEN
          write(lupri,*)'D: nj=',nj,' expon=',expon,' contr.coef=',coe
           ENDIF

                do i=1,6
                   shell_arr(sh+(i-1)*n_con_guass)=shell_d(i)
CMI NEW
                   Shell_basis(g_n+i-1)=shell_d(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe
                   do ia_p=1,3
                 ang_pow(sh+(i-1)*n_con_guass,ia_p)= angu_pow_d(i,ia_p)
                   end do
                end do
             END DO
       
             do  iji=1,6
                basis_n(g_n)=basis_n(g_n)+n_con_guass
                g_n=g_n+1
             end do
             sh=sh+5*n_con_guass
             go to 600
             
 500         continue 
             
             DO  n_c_g=1,n_con_guass
                sh=sh+1
                READ(LUinp,*)nj,expon,coe

           IF (IPRINT.GE.3) THEN
          write(lupri,*)'F: nj=',nj,' expon=',expon,' contr.coef=',coe
           ENDIF

                do i=1,10
                   shell_arr(sh+(i-1)*n_con_guass)=shell_f(i)
CMI NEW
                   Shell_basis(g_n+i-1)=shell_f(i)
                   Shell_atom(g_n+i-1)=k

                   ord_expon(sh+(i-1)*n_con_guass)=expon
                   ord_coe(sh+(i-1)*n_con_guass)=coe
                   do ia_p=1,3
                   ang_pow(sh+(i-1)*n_con_guass,ia_p)=angu_pow_f(i,ia_p)
                   end do
                end do
             END DO
             
             do  iji=1,10
                basis_n(g_n)=basis_n(g_n)+n_con_guass
                g_n=g_n+1
             end do
             sh=sh+9*n_con_guass
             
 600         continue
             
          End  Do
 201       continue
          n_pg_atom(1)=1
          n_pg_atom(k+1)=sh
          do ii= n_pg_atom(k),sh
             do isw=1,3
                R_ms(iiii,isw)=re(k,isw)
             end do
             iiii=iiii+1
          end do
       END DO
       n_b=g_n-1
       
C     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
       write(16,*)'number  of primitive guassians=',sh
       write(16,*)'number of bais',n_b 
       write(16,*)'no. of pri gua in each basis'
CMI NEW
        LUCAP=0
        CALL GPOPEN(LUCAP,
     &  'CAPLBL.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)
        WRITE(LUCAP,*) 'CAP unique order of atomic basis functions'

       do i=1,n_b
C        write(16,*)i, basis_n(i)
         write(16,'(2i4,A6)')i,SHELL_ATOM(i),Shell_basis(i)
         WRITE(LUCAP,'(i4,2x,i3,2x,A6)')i,SHELL_ATOM(i),Shell_basis(i)
       end do

       CALL GPCLOSE(LUCAP,'KEEP')
       IF (IPRINT.GE.1) THEN
        WRITE(LUPRI,'(2X,A)')
     &  'File "CAPLBL.inp" with CAP unique order of atomic basis'//
     &  ' functions ready.'
       ENDIF
       
       write(16,*)'sh_ty,Pgua_expo,nor_pg_coe,angu_pow,cartesian'
       
       do i=1,sh
          do isw=1,3
             R_m(isw)=R_ms(i,isw)
             R_n(isw)=R_ms(i,isw)
             k_m(isw)=(ang_pow(i,isw))
             k_n(isw)=(ang_pow(i,isw))
          end do
          
          R_m(4)=ord_expon(i)
          R_n(4)=ord_expon(i)
          
          CALL PEM(R_m,R_n,K_m,K_n,PW,PEINT)
          ord_coe(i)=ord_coe(i)/sqrt(PEINT)
          write (16,*)i,shell_arr(i),ord_expon(i),ord_coe(i)
          write(16,*) (ang_pow(i,j),j=1,3)
          write(16,*) (R_ms(i,j),j=1,3)
       end do

CMI ... I have to close it due to interference with her1pro.F:3403 
CMI    close(19,status="keep")
       close(LUinp,status="keep")
C      CALL GPCLOSE(LUinp,'KEEP')

C       call quit('CAP:INPREAD')
       CALL QEXIT('INPREAD')
       RETURN
       END

C*******************************************************************************
      SUBROUTINE AOINTEGRAL(N_B,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     *     PW,PE_AO, IPRINT )      
CMI
C
C
C   entering PW is from INPQD file
C
C
C  Entering big PE_AO arrary !!!
C
C
#include "implicit.h"
#include "priunit.h"
CMI   IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 n_b
      integer *4 ang_pow(250,3)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(200)
      integer *4 basis_ord(300)
     
      real *8 ord_expon(250),ord_coe(250)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(200,3)

      COMPLEX *16 PEINT,PE_AO(25000)
      COMPLEX *16 BETA

      CALL QENTER('AOINTEGRAL') 
      
      do i=1,n_b
         do j=1,n_b
            IPU=IPOINT(i,j)
            PE_AO(IPU)=(0.d0,0.d0)
         end do
      end do

      do i=1,n_b
         if(i.eq.1)then
            iasp=0
         else
            iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
      end do
      
      DO ii=1,n_b
C         write(*,*)ii
      DO jj=1,ii
         IPU=IPOINT(ii,jj)

         if (ii.eq.1)then
            iangii=0
         else
            iangii=basis_ord(ii-1)
         endif
         if (jj.eq.1)then
            iangjj=0
         else
            iangjj=basis_ord(jj-1)
         endif
         
         Do i=iangii+1,basis_ord(ii)
         do j=iangjj+1,basis_ord(jj)

            do isw=1,3
               R_m(isw)=R_ms(i,isw)
               R_n(isw)=R_ms(j,isw)
               k_m(isw)=(ang_pow(i,isw))
               k_n(isw)=(ang_pow(j,isw))
            end do
            R_m(4)=ord_expon(i)
            R_n(4)=ord_expon(j)

            COMUL=ord_coe(i)*ord_coe(j)
C            CALL PEMC(R_m,R_n,K_m,K_n,BETA,PEINT)
            CALL CAPM(R_m,R_n,K_m,K_n,PW,PEINT)
          
            PE_AO(IPU)= PE_AO(IPU)+COMUL*PEINT
         end do
      End Do
C      If(ii.NE.jj) PE_AO(jj,ii)= PE_AO(ii,jj)
      END DO
      END DO
      
      CALL QEXIT('AOINTEGRAL') 
      RETURN
      END
     
C************************************************************
CMI ... not used here ...

C     SUBROUTINE MOINTEGRAL(N_B,AOMO,AOINT,MOINT)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     real *8 AOMO(200,200)
C     COMPLEX *16 AOINT(25000),MOINT(25000)
      
C     Do i=1,n_b
C        do j=1,i
C           IPU=IPOINT(i,j)
C           MOINT(IPU)=(0.d0,0.d0)
C        end do
C     End do
      
C     DO ii=1,n_b
C        Do jj=1,ii
C            IPU=IPOINT(ii,jj)
C           do i=1,n_b
C              do j=1,n_b
C                 IPS=IPOINT(i,j)
C                 if(j.gt.i)IPS=IPOINT(j,i)
C                 COMULT=AOMO(jj,j)*AOMO(ii,i)
C                 MOINT(IPU)= MOINT(IPU)+COMULT*AOINT(IPS)
C              end do
C           end do
C        End  Do
C     END DO
      
C     RETURN
C     END

C***************************************************
CMI ... not used here....
C     SUBROUTINE CIHAMILTONIAN(N_B,NSIZE,MOINT,DETINT)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     COMPLEX *16  MOINT(25000),DETINT(10000000)
C     COMPLEX *16 ELEM
C     DIMENSION INDEXI((N_B*(N_B+1))/2,(N_B*(N_B+1))/2)
     

C     DO 7 I=1,(N_B*(N_B+1))/2                                        
C        DO 8 J=1,I                                                     
C           INDEXI(I,J) = I*(I-1)/2 + J                                 
C           INDEXI(J,I) = INDEXI(I,J)                 
C8       CONTINUE                  
C7    CONTINUE       

C     OPEN (12,FILE='CIINDEX',FORM='FORMATTED',ACCESS='SEQUENTIAL')
C     REWIND 12
C     READ (12,*)NSIZE
C     DO ICI=1,NSIZE
C     DO JCI=1,ICI
C        Read(12,*)IDET,JDET,IACON1,IACON2,IBCON1,IBCON2
C        CALL HELEM(N_B,MOINT,IACON1,IBCON1,IACON2,IBCON2,IDET,
C    *        JDET,INDEXI,ELEM)
C        IPU=IPOINT(IDET,JDET)
C        DETINT(IPU)=ELEM
C     END DO
C     END DO

C         RETURN
C         END

C     --------------------------------------------------------          
CMI ... not used here....
C     SUBROUTINE HELEM(N_B,SINT1,IACON1,IBCON1,      
C    *     IACON2,IBCON2,IJ,JI,INDEXI,ELEM)               
C     --------------------------------------------------------          

C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     INTEGER DIFF1(2),DIFF2(2),IPOS1(2),IPOS2(2)      
C     COMPLEX *16 SINT1(25000),ELEM
C     DIMENSION IACON1(1),IBCON1(1)                             
C     DIMENSION IACON2(1),IBCON2(1)                             
C     DIMENSION INDEXI((N_B*(N_B+1))/2,(N_B*(N_B+1))/2)              
                 
C                                                                       
C     Returns the matrix element < K | H | L > where K, L are           
C     determinants.  Alpha and Beta occupations are stored in           
C     iacon1, ibcon1 for K and iacon2, ibcon2 for L.                    
C                                                                       
C     ELEM = (0.0D0,0.d0)   

C     NAELE=1
C     NBELE=1
C     NCOR=0
C                                                                       
C                                                                       
C    If Determinants are same                                           
C                                                                       
C     IF (IJ.EQ.JI) THEN                                                
C        DO 200 I=1,NAELE                                               
C           IA = IACON1(I)                                              
C           I1 = INDEXI(IA,IA)     
C            if(ij.eq.10)write(*,*)I1,   SINT1(I1)                        
C           ELEM = ELEM + SINT1(I1)

C200        END DO

C        DO 210 I=1,NBELE                                               
C           IB = IBCON1(I)                                              
C           I1 = INDEXI(IB,IB)                                           

C           ELEM = ELEM + SINT1(I1)
C210         END DO
C     RETURN                                                            
C     ENDIF                                                             
C                                                                       
C     IDEA=0                                                            
C     IDEB=0                                                            
C                                                                       
C     Different orbitals in first deterinant                            
C                                                                       
C     DO 20 I=NCOR+1,NAELE                                              
C        DO 15 J=NCOR+1,NAELE                                           
C           IF (IACON1(I).EQ.IACON2(J)) GOTO 20                         
C  15    CONTINUE                                                       
C        IDEA = IDEA + 1                                                
C        IF (IDEA.GT.2) RETURN                                          
C        DIFF1(IDEA) = IACON1(I)                                        
C        IPOS1(IDEA) = I                                                
C  20 CONTINUE                                                          
C                                                                       
C     DO 30 I=NCOR+1,NBELE                                              
C        DO 25 J=NCOR+1,NBELE                                           
C           IF (IBCON1(I).EQ.IBCON2(J)) GOTO 30                         
C  25    CONTINUE                                                       
C        IDEB = IDEB + 1                                                
C        IF (IDEA+IDEB.GT.2) RETURN                                     
C        DIFF1(IDEA+IDEB) = IBCON1(I)                                   
C        IPOS1(IDEA+IDEB) = I                                           
C  30 CONTINUE                                                          
C                                                                       
C                                                                       
C    To find the different orbitals in second determinant               
C                                                                       
C     IST = NCOR+1                                                      
C     DO 63 II=1,IDEA                                                   
C           DO 50 I=IST,NAELE                                           
C              DO 45 J=NCOR+1,NAELE                                     
C                 IF (IACON2(I).EQ.IACON1(J)) GOTO 50                   
C  45          CONTINUE                                                 
C              GOTO 60                                                  
C  50       CONTINUE                                                    
C                                                                       
C  60       DIFF2(II) = IACON2(I)                                       
C           IPOS2(II) = I                                               
C           IST = I+1                                                   
C  63 CONTINUE                                                          
C                                                                       
C     IST = NCOR+1                                                      
C     DO 163 II=1,IDEB                                                  
C           DO 150 I=IST,NBELE                                          
C              DO 145 J=NCOR+1,NBELE                                    
C                 IF (IBCON2(I).EQ.IBCON1(J)) GOTO 150                  
C 145          CONTINUE                                                 
C              GOTO 160                                                 
C 150       CONTINUE                                                    
C 160       DIFF2(II+IDEA) = IBCON2(I)                                  
C           IPOS2(II+IDEA) = I                                          
C           IST = I+1                                                   
C 163 CONTINUE                                                          

C    If determinants differ by one orbital                              

C      IF (IDEA+IDEB.EQ.1) THEN                                         

C   One particle density contribution                                   

C         IND1 = INDEXI(DIFF1(1),DIFF2(1))                               
C         IPERM = (-1)**(IPOS1(1)-IPOS2(1))
                             
C         ELEM = ELEM + IPERM*SINT1(IND1)  

C         ENDIF                                                         
C         RETURN                                                        

C     RETURN                                                            
C     END                  
   
      SUBROUTINE PEMC(R_m,R_n,K_m,K_n,BETA,PEINT)
      IMPLICIT REAL*8 (A-H,O-Z)

      
      real *8 R_CM(2),R_CN(2),R_KM(2),R_KN(2),R_LM(2),R_LN(2)
      real *8 R_m(4),R_n(4)
      integer *4 K_m(3),K_n(3)
      COMPLEX *16 BETA,OVRLPI,OVRLPJ,OVRLPK,PEINT
      
C      write(*,*)BETA
         PEINT=(0.d0,0.d0)
         i=1
         k=2
         L=3
         R_CM(1)=R_m(i)
         R_KM(1)=R_m(k)
         R_LM(1)=R_m(L)

         R_CN(1)=R_n(i)
         R_KN(1)=R_n(k)
         R_LN(1)=R_n(L)

         R_CM(2)=R_m(4)
         R_KM(2)=R_m(4)
         R_LM(2)=R_m(4)

         R_CN(2)=R_n(4) 
         R_KN(2)=R_n(4)
         R_LN(2)=R_n(4)
  
         K_CM=K_m(i)
         K_KM=K_m(k)
         K_LM=K_m(L)

         K_CN=K_n(i)
         k_KN=K_n(k)     
         K_LN=k_n(L)
         CALL SMATC(R_KM,R_KN,K_KM,K_KN,BETA,OVRLPI)
         CALL SMATC(R_LM,R_LN,K_LM,K_LN,BETA,OVRLPJ)
         CALL SMATC(R_CM,R_CN,K_CM,K_CN,BETA,OVRLPK)
         PEINT=PEINT+OVRLPI*OVRLPJ*OVRLPK
       RETURN
       END

      subroutine SMATC(R_m,R_n,K_m,K_n,BETA,Int_GR)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8  R_m(2),R_n(2)
      COMPLEX *16 BETA,INT_GR,PIG,PREMULT,PA,PB
      COMPLEX *16 GAM,P,GAMP,RMULT,RMULTI,FINTM
      data PI/3.14/
      EXTERNAL  IFACT
   
      GAM=R_m(2)+(R_n(2)+BETA)
      P=(R_m(1)*R_m(2)+R_n(1)*(R_n(2)+BETA))/GAM
      PA=P-R_m(1)
      PB=P-R_n(1)
      AB=R_m(1)-R_n(1) 

      DO i=0,(K_m+k_n)/2
         INT_GR=(0.d0,0.d0)
         if (i.eq.0)go to 150
     
 150     continue
         IFI=IFACT(2*i)
         GAMP=(2*gam)**(-i)
         PIG=(PI/GAM)**0.5
         RMULT=IFI*PIG*gamp
         CALL F_INTC(2*i,PA,PB,k_m,K_n,FINTM)
         RMULTI=RMULT*FINTM 
         INT_GR=INT_GR+RMULTI 
 100     continue
      END DO
      
      PREMULT=exp((-1.d0*R_m(2)*(R_n(2)+BETA)*((AB)**2))/gam)
C      if (DABS(premult).lt.0.000000000d0)PREMULT=0.000000000d0
      INT_GR=INT_GR*PREMULT
      RETURN
      END

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION IFACT(l)
      IMPLICIT REAL*8 (A-H,O-Z)
      ifact=1
      if((l-1).lt.2)go to 200
      do j=1,(l)-1
         IF ((J/2)*2.EQ.j)GO TO 100
         IFact=ifact*j
 100     CONTINUE
      end do
      go to 300
 200  ifact=1
 300  continue
      END 
      
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE F_INTC(k,PA,PB,k_m,K_n,FINTM) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8 Int
      EXTERNAL JFACT
      COMPLEX *16 PA,PB,FINTM
      COMPLEX *16 PAMUL,PBMUL,RADD

      FINTM=(0.d0,0.d0)
      do i=0,(K_m)
         do j=0,k_n
            if ((i+j).eq.k)then
               PAMUL= PA**(K_m-i)
               PBMUL=PB**(K_n-j)
               JFKM=JFACT(k_m)
               JFKN=JFACT(k_n)
               JFI=JFACT(i)
               JFJ=JFACT(j)
               JFIM=JFACT(K_m-i)
               JFJN=JFACT(K_n-j)
               RMFACT=JFKM/(JFI*JFIM)
               RNFACT=JFKN/(JFJ*JFJN)
               RADD=PAMUL*PBMUL*RMFACT*RNFACT
               FINTM=FINTM+RADD
            endif
         end do 
      end do
      RETURN
      END
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      FUNCTION JFACT(l)
      IMPLICIT REAL*8 (A-H,O-Z)
      jfact=1
      if((l).lt.1)go to 200
      do j=1,l
         jFact=jfact*j
      end do
      go to 300
 200  jfact=1
 300  continue
      END  
      
      SUBROUTINE AOMOCOE(COE)
      IMPLICIT REAL*8 (A-H,O-Z)
      real *8  coe(200,200)
      
      OPEN(4,FILE='COINP',FORM='FORMATTED',ACCESS='SEQUENTIAL',
     *     STATUS='UNKNOWN')
      OPEN(3,FILE='AOMOCOE',FORM='FORMATTED',ACCESS='SEQUENTIAL',
     *     STATUS='UNKNOWN')
      

      irb=5
      kr=1
      
      read(4,*)nbasis
      kl=nbasis/5
      if(kl*5.NE.nbasis)then
         kls=nbasis-(kl*5)
         kl=kl+1
      endif
      kr=0
      Do kkkk=1,nbasis
         Do kk=1,kl
            
            if (kk.eq.kl.and.kl*5.NE.NBASIS)then
               krb=kls
            else
               krb=5
            endif
            read(4,*)kcol,icol,(coe(kcol,is+(kk-1)*5),is=1,krb)
   
         kr=kr+5
      End do
      End do

 110  format(i3,i2,5E15.12)

      do i=1,nbasis
         write(3,*)(coe(i,j),j=1,nbasis)
      end do
      CLOSE(3)
      CLOSE(4)
      
      RETURN
      END

      FUNCTION IPOINT(ii,jj)
      IMPLICIT REAL*8 (A-H,O-Z)
      IPOINT=(ii*(ii-1)/2)+jj
      END 

C**********************************************
 
  
   
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE PEM(R_m,R_n,K_m,K_n,BETA,PEINT)
      IMPLICIT REAL*8 (A-H,O-Z)


      
      real *8 R_CM(2),R_CN(2),R_KM(2),R_KN(2),R_LM(2),R_LN(2)
      real *8 R_m(4),R_n(4)
      integer *4 K_m(3),K_n(3)
      real *8 BETA,OVRLPI,OVRLPJ,OVRLPK,PEINT
C      write(*,*)BETA
         PEINT=0.d0
         i=1
         k=2
         L=3
         R_CM(1)=R_m(i)
         R_KM(1)=R_m(k)
         R_LM(1)=R_m(L)

         R_CN(1)=R_n(i)
         R_KN(1)=R_n(k)
         R_LN(1)=R_n(L)

         R_CM(2)=R_m(4)
         R_KM(2)=R_m(4)
         R_LM(2)=R_m(4)

         R_CN(2)=R_n(4) 
         R_KN(2)=R_n(4)
         R_LN(2)=R_n(4)
  
         K_CM=K_m(i)
         K_KM=K_m(k)
         K_LM=K_m(L)

         K_CN=K_n(i)
         k_KN=K_n(k)     
         K_LN=k_n(L)
C         write(*,*)'SMAT'
         CALL SMAT(R_KM,R_KN,K_KM,K_KN,0.d0,OVRLPI)
C         write(*,*)OVERLPI
         CALL SMAT(R_LM,R_LN,K_LM,K_LN,0.d0,OVRLPJ)
         CALL SMAT(R_CM,R_CN,K_CM,K_CN,0.d0,OVRLPK)
         PEINT=PEINT+OVRLPI*OVRLPJ*OVRLPK
C         write(*,*)PEINT
       RETURN
       END
 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SMAT(R_m,R_n,K_m,K_n,BETA,Int_GR)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8  R_m(2),R_n(2)
      real *8 BETA, INT_GR
      real *8 GAM,P,PA,PB,GAMP,RMULT,RMULTI,FINTM
      data PI/3.14/
      EXTERNAL  IFACT
   
      GAM=R_m(2)+(R_n(2)+BETA)
      P=(R_m(1)*R_m(2)+R_n(1)*(R_n(2)+BETA))/GAM
      PA=P-R_m(1)
      PB=P-R_n(1)
      AB=R_m(1)-R_n(1) 

      DO i=0,(K_m+k_n)/2
C         if (i.eq.1)write(*,*)i
         INT_GR=0.d0
         if (i.eq.0)go to 150
C         if ((i/2)*2.ne.i)go to 100
     
 150     continue
         IFI=IFACT(2*i)
C         if (i.eq.4)write(*,*)IFACT(i),'IFA'
         GAMP=(2*gam)**(-i)
         PIG=(PI/GAM)**0.5
         RMULT=IFI*PIG*gamp
         CALL F_INT(2*i,PA,PB,k_m,K_n,FINTM)
         RMULTI=RMULT*FINTM 
         INT_GR=INT_GR+RMULTI 
 100     continue
      END DO
      
      PREMULT=exp((-1.d0*R_m(2)*(R_n(2)+BETA)*((AB)**2))/gam)
      if (dABS(premult).lt.0.000000000d0)PREMULT=0.000000000d0
      INT_GR=INT_GR*PREMULT
      RETURN
      END

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE F_INT(k,PA,PB,k_m,K_n,FINTM) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 k_m,k_n
      real *8 Int
C     REAL *8 PA,PB
      EXTERNAL JFACT
      real *8 PA,PB,FINTM
      real *8 PAMUL,PBMUL,RADD

      FINTM=(0.d0,0.d0)
      do i=0,(K_m)
         do j=0,k_n
            if ((i+j).eq.k)then
               PAMUL= PA**(K_m-i)
               PBMUL=PB**(K_n-j)
               JFKM=JFACT(k_m)
               JFKN=JFACT(k_n)
               JFI=JFACT(i)
               JFJ=JFACT(j)
               JFIM=JFACT(K_m-i)
               JFJN=JFACT(K_n-j)
               RMFACT=JFKM/(JFI*JFIM)
               RNFACT=JFKN/(JFJ*JFJN)
               RADD=PAMUL*PBMUL*RMFACT*RNFACT
               FINTM=FINTM+RADD
            endif
         end do 
      end do
      RETURN
      END



C*******************************************************
CMI entering theta is from the INPQD file  
      SUBROUTINE CAPM(R_m,R_n,K_m,K_n,theta,VCAP)
C      IMPLICIT REAL*8 (A-H,O-Z)
#include "implicit.h"
#include "priunit.h"

      integer *4 NSTEP
      integer *4 NGR
      
      real *8 R_CM(2),R_CN(2),R_KM(2),R_KN(2),R_LM(2),R_LN(2)
      real *8 R_m(4),R_n(4)
      integer *4 K_m(3),K_n(3)
      real *8 CAP_PAR(3)
CMI   REAL *8 PGRT(100),GINTPT(100)
      REAL *8 PGRT(100)
CMI   real *8 X0,theta,LAMDA,GRSTP(2)
      real *8 theta,GRSTP(2)
CMI   integer *4 NGRS(100)
      complex  *16 PE,PE1
      complex *16 VCAP1,VCAP,VCAPS(3)
      COMMON PGRT
#include "dcbcap.h"

      CALL QENTER('CAPM')

CMI   .... no need for reading this file !!!
C     open(unit=8,file='INPCAP',STATUS='UNKNOWN',FORM='FORMATTED')
C     REWIND (8)
     
      PE=(0.d0,0.d0)
   
      DO i=1,3
       
CMI      ... replaced by X0(1:3), LAMDA(1:3) in dcbcap.h
CMI      READ(8,*)X0,LAMDA
CMI      .... replaced by NGR(1:3) in dcbcap.h
CMI      READ(8,*)NGR
CMI      .... replaced by GINTPT(3,MAXSPL) in dcbcap.h
CMI      READ(8,*)(GINTPT(j),j=1,NGR)
CMI       .... replaced by NGRS(3,MAXSPL) in dcbcap.h
CMI      READ(8,*)(NGRS(j),j=1,NGR-1)
   
      CAP_PAR(3)=THETA
CMI   CAP_PAR(2)=LAMDA
      CAP_PAR(2)=LAMDA(I)
CMI   CAP_PAR(1)=X0
      CAP_PAR(1)=X0(I)
      
         if(i.eq.1)then
            k=2
            L=3
         endif
         if(i.eq.2)then
            k=1
            L=3
         endif
         if(i.eq.3)then
            k=1
            L=2
         endif


         R_CM(1)=R_m(i)
         R_KM(1)=R_m(k)
         R_LM(1)=R_m(L)

         R_CN(1)=R_n(i)
         R_KN(1)=R_n(k)
         R_LN(1)=R_n(L)

         R_CM(2)=R_m(4)
         R_KM(2)=R_m(4)
         R_LM(2)=R_m(4)

         R_CN(2)=R_n(4) 
         R_KN(2)=R_n(4)
         R_LN(2)=R_n(4)

         K_CM=K_m(i)
         K_KM=K_m(k)
         K_LM=K_m(L)

         K_CN=K_n(i)
         k_KN=K_n(k)     
         K_LN=k_n(L)
       
       
         VCAPS(i)=(0.d0,0.d0)

CMI      do j=1,NGR-1
         do j=1,NGR(i)-1
CMI         GRSTP(1)=GINTPT(j)
            GRSTP(1)=GINTPT(i,j)
CMI         GRSTP(2)=GINTPT(j+1)
            GRSTP(2)=GINTPT(i,j+1)
CMI         NSTEP=NGRS(j)
            NSTEP=NGRS(i,j)
            CALL CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_CM,R_CN,K_CM,K_CN,VCAP1)
            VCAPS(i)=VCAPS(i)+VCAP1
         end do

         CALL SMAT(R_KM,R_KN,K_KM,K_KN,0.d0,OVRLPK)
         CALL SMAT(R_LM,R_LN,K_LM,K_LN,0.d0,OVRLPL)
         VCAPS(i)=VCAPS(i)*OVRLPK*OVRLPL
C         write(23,*)OVRLPK,OVRLPL

      END DO ! of I=1,3
             
       VCAP=VCAPS(1)+VCAPS(2)+VCAPS(3)

CMI    CLOSE(8) 

      CALL QEXIT('CAPM')
      RETURN
      END

C*************************************************************

      subroutine CAPMAT1(CAP_PAR,GRSTP,NSTEP,R_m,R_n,K_m,K_n,VCAPT)
      IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 k_m,k_n
      integer *4 NSTEP
      real *8 R_m(2),R_n(2)
      real *8 CAP_PAR(3)
      real *8 M
      real *8 X0,theta,LAMDA,GRSTP(2)
      complex *16 FFx,fx,Dfx,D2fx
      complex *16 V0x,V1x,V2x
      
      complex *16 CI,ROT
      complex *16 PEM,XS,PEN

      COMPLEX *16 VCAPT,VCAPXT,VCAPYT,VCAPZT
      COMPLEX *16 VCAPX(4),VCAPY(2),VCAPZ(2)
      COMPLEX *16 VCAP,VADD,VF
      EXTERNAL V0,V1,V2,FF
      double complex v0,v1,v2,FF
      
      data CI/(0.0d0,1.0d0)/
      data M/1.0d00/
     

      STEP=(GRSTP(2)-GRSTP(1))/(NSTEP)
      X0=CAP_PAR(1)
      LAMDA=CAP_PAR(2)
      THETA=CAP_PAR(3)
    
      VCAPXT=(0.d0,0.d0)
      XI=GRSTP(1)
      itx=0
 
      DO ix=1,NSTEP
         X=(XI+(ix-1)*STEP)
         itx=itx+1

         XR=X-R_n(1)
         XRM=X-R_m(1)
         DGr=K_n*(XR**-1)
         DG=DGr-(2*R_n(2)*(XR))
       
         D2G1X=K_n*(K_n-1)*(XR**-2)
         D2G2X=2*R_n(2)*(2*K_n+1)
         D2G3X=4*(R_n(2)**2)*(XR**2)
         D2GX=D2G1X-D2G2X+D2G3X
         
         D2G=D2GX

         GG1M=(XRM**2)
         GG1N=(XR**2)
         GG1=R_m(2)*GG1M+R_n(2)*GG1N

         GGEX=(XR**K_n)*(XRM**K_m)
         GG=GGEX*DEXP(-GG1)
       
         v0x=V0(M,THETA,LAMDA,X0,X)
         V1x=v1(M,THETA,LAMDA,X0,X)*DG
         v2x=v2(M,THETA,LAMDA,X0,X)*D2G
C         VADD=1/x-1/FF(THETA,LAMDA,X0,X)
         
C         XDC=5.0
C         IF(X.ge.XDC)then
C         VF=X-XDC
C         VADD=-0.05*VF
C         else
         VADD=0.0d0
C         endif

         VCAP=(V0x+V1x+V2x+VADD)*GG

         VCAPX(itx)=VCAP
         if(ix.eq.1) go to 321 
         
         VCAPXT=VCAPXT+0.5*(VCAPX(2)+VCAPX(1))*STEP
         if(ix.ne.1)VCAPX(2)=VCAPX(1)
         itx=0
C         if((zabs(vcap)).lt.(0.00000000001))exit

 321     CONTINUE

      END DO
      VCAPT=VCAPXT*2
      RETURN
      END

C************************************************************************
C     Implements the basic functions from the paper.
C************************************************************************
      DOUBLE COMPLEX FUNCTION FF(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = DCOSH(LAMBDA*(X-X0))/DCOSH(LAMBDA*(X+X0))
      AUX = X + DLOG(AUX)/(2.0D0*LAMBDA)
      FF = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      FF = DCMPLX(X) + FF*DCMPLX(AUX)
      END


      DOUBLE COMPLEX FUNCTION F0(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = DTANH(LAMBDA*(X-X0)) - DTANH(LAMBDA*(X+X0))
      AUX = 1.0D0 + AUX/2.0D0
      F0 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F0 = (1.0D0,0.0D0) + F0*DCMPLX(AUX)
      END

      DOUBLE COMPLEX FUNCTION F1(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8  THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX
      AUX = (DTANH(LAMBDA*(X+X0))**2) - (DTANH(LAMBDA*(X-X0))**2)
      F1 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F1 = DCMPLX(LAMBDA)*F1*DCMPLX(AUX)/(2.0D0,0.0D0)
      END


      DOUBLE COMPLEX FUNCTION F2(THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 THETA,LAMBDA,X0,X
      DOUBLE PRECISION AUX1,AUX2
      AUX1 = DTANH(LAMBDA*(X+X0))*(1.0D0-(DTANH(LAMBDA*(X+X0))**2))
      AUX2 = DTANH(LAMBDA*(X-X0))*(1.0D0-(DTANH(LAMBDA*(X-X0))**2))
      F2 = CDEXP(+(0.0D0,1.0D0)*DCMPLX(THETA))-(1.0D0,0.0D0)
      F2 = DCMPLX(LAMBDA**2)*F2*(DCMPLX(AUX1)-DCMPLX(AUX2)) 
      END


      DOUBLE COMPLEX FUNCTION V0(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8  M,THETA,LAMBDA,X0,X
      DOUBLE COMPLEX CA1,CA2
      EXTERNAL F0,F1,F2
      DOUBLE COMPLEX F0,F1,F2
      CA1 = F2(THETA,LAMBDA,X0,X)/(F0(THETA,LAMBDA,X0,X)**3)
      CA1 = CA1/(DCMPLX(4.0D0*M))
      CA2 = (F1(THETA,LAMBDA,X0,X)**2)/(F0(THETA,LAMBDA,X0,X)**4)
      CA2 = CA2/DCMPLX(5.0D0*M/8.0D0)
      V0 = CA1 - CA2
      END

      DOUBLE COMPLEX FUNCTION V1(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 M,THETA,LAMBDA,X0,X
      EXTERNAL F0,F1
      DOUBLE COMPLEX F0,F1
      V1 = F1(THETA,LAMBDA,X0,X)/(F0(THETA,LAMBDA,X0,X)**3)
      V1 = V1/DCMPLX(M)
      END

      DOUBLE COMPLEX FUNCTION V2(M,THETA,LAMBDA,X0,X)
      IMPLICIT NONE
      real *8 M,THETA,LAMBDA,X0,X
      EXTERNAL F0
      DOUBLE COMPLEX F0
      V2 = (1.0D0,0.0D0)/(F0(THETA,LAMBDA,X0,X)**2)
      V2 = ((1.0D0,0.0D0)-V2)/DCMPLX(2.0D0*M)
      END
