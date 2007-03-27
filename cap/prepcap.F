C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
      PARAMETER (NTABLE = 20)
      PARAMETER (D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF,RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
#include "dcbham.h"
#include "dcbcap.h"
C
      SAVE SET
      DATA TABLE /'.INPQD ','.INPCAP','.PRINT ','.OVLCAP',
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

      CALL TITLER('CAP operator defined','*',127)

C
C  *** Initialize ****
C
      CAPOVLP = .FALSE.

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
C&&& Comparison between DIRAC and CAP overlap integrals
                 CAPOVLP = .TRUE.
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
     & 'CAPINP: Set-up for CAP calculations'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)') ' * Print level: ',IPRCAP
C
      IF (.NOT.CAPOVLP) THEN
        WRITE(LUPRI,'(A)')
     &  ' * CAP integrals  according to',
     &  '          Sajeev ....         '
C
        WRITE(LUPRI,'(2X,A,D25.10)') 'read INPQD = ',PW
        WRITE(LUPRI,'(2X,A)') 
     & 'read integration parameters:'
CMI    ...  write out integration parameters
        DO I=1,3
          WRITE(LUPRI,*) X0(I),LAMDA(I)
          WRITE(LUPRI,*) NGR(I)
          WRITE(LUPRI,*) (GINTPT(I,J),J=1,NGR(I))
          WRITE(LUPRI,*) (NGRS(I,J),J=1,NGR(I)-1)
        ENDDO
       ELSE
         WRITE(LUPRI,'(2X,A)') 
     &   'Only comparison between CAP and DIRAC overlap integrals.'
       ENDIF

  999 CONTINUE
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capints */
      SUBROUTINE CAPINTS(WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Interface for the CAP integral program.
C
C   Called from PAMINP.
C
C   Written by Miro ILIAS, Tel Aviv, January 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       USE CAP
       USE CAP, only: cap_init

#include "implicit.h"
#include "priunit.h"
C#include "dcbham.h"

      DIMENSION WORK(LWORK)
#include "memint.h"
      CALL QENTER('CAPINTS')

      CALL MEMCHK('CAPINTS',WORK,1)

C    ...  call F90 routine from the module cap
       call cap_init

      ISIZE1=200
      CALL MEMGET('INTE',KBASIS_N,ISIZE1,WORK,KFREE,LFREE)
      ISIZE2=250*3
      CALL MEMGET('INTE',KANG_POW,ISIZE2,WORK,KFREE,LFREE)

      ISIZE3=250
      CALL MEMGET('REAL',KORDEXPON,ISIZE3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORDCOE,ISIZE3,WORK,KFREE,LFREE)

      ISIZE4=200*3
      CALL MEMGET('REAL',KRMS,ISIZE4,WORK,KFREE,LFREE)

      ISIZE5=25000
      CALL MEMGET('REAL',KPE_AOr,ISIZE5,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPE_AOi,ISIZE5,WORK,KFREE,LFREE)

CMI:    allocated arrays enters the SAJEEV routine !
      CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr),WORK(KPE_AOi),
     &            WORK(KFREE),LFREE)

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
#include "dgroup.h"
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
      LOGICAL CHPRES,OLDADAPT
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
          write(LUPRI,'(I4,A,I3,A,A4,A,I2,A,I3)')
     &    II1,'> nucleus:',IOBF1(1,I),' funct: ',AO1,'/',IOBF1(2,I),
     &     ' >  CAP unique fu #:',
     &        IOBF1(3,I)
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

C       IF (IPRCAP.GE.2) THEN
C        WRITE(LUPRI,'(I4,A,I4)')
C    &   I,'. CAP bf is DIRAC unique fu:',ICAPBFO(I)
C       ENDIF
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

      LUPM = 0
      CALL GPOPEN(LUPM,'PEMAT','OLD',' ','FORMATTED',
     &       IDUMMY,.FALSE.)

      ISIZE1 = NBASIS*(NBASIS+1)/2
      ISIZE = ( NBASIS*NBASIS )
      IF (ISIZE1.NE.NNBBASX) THEN
       WRITE(LUPRI,*) "ISIZE=",ISIZE1," NNBBASX=",NNBBASX
       CALL QUIT('CAPADAPT: ISIZE1 <> NNBBASX !')
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
        CAPRE = REAL(CAP)
        CAPIM = IMAG(CAP)
        IF (I.NE.II.OR.J.NE.JJ) THEN
         CALL QUIT('CAPADAPT: wrong indexes in PEMAT !')
        ENDIF
        IX1 = CMNINDX(I,J)
        IF (IX.NE.IX1) THEN
         CALL QUIT('CAPADAPT: wrong IX !')
        ENDIF

C       -------------------------------------------------
C       .. store  integrals in the DIRAC basic order ..
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

      OLDADAPT = .FALSE.
C ======================================================================

      IF (OLDADAPT) THEN

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

      ELSE

       CALL NTAOSO(WORK(KPMR),WORK(KPMRSA),NBASIS,0,IPRCAP)

       IF (.NOT.CAPOVLP) THEN
         CALL NTAOSO(WORK(KPMI),WORK(KPMISA),NBASIS,0,IPRCAP)
       ELSE
         CALL DZERO(WORK(KPMISA),ISIZE)
       ENDIF

      ENDIF

C =======================================================================
C     Save HERMIT symmetry adapted integrals into the AOPROPER file.
C =======================================================================
      CALL SAVECAP(WORK(KPMRSA),WORK(KPMISA),IPRCAP)

C     CALL CAPREADY(WORK,LWORK,IPRCAP)
      CALL CAPREADY(WORK(KFREE),LFREE,IPRCAP)

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
C   Called from MRCON1 (traone.F)
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

      DIMENSION CAPI(NORBT,NORBT,NZ),WORK(*)

      CALL QENTER('GETCAPI')

      KFRSAV = KFREE

      IPRINT = IPRCAP

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
     &               ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(2X,2A)')
     & 'GETCAPI: Get MO matrix of property: ',PRPNAM(INDXPR)
       CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

      CALL MEMREL('GETCAPI',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GETCAPI')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck addcapim */
      SUBROUTINE ADDCAPIM(NSTRT,CAPIM,FMOM,MULTIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: add MOLFDIR like i.CAP_IM to the FMOM
C
C  On input: CAPIM - to be added to FMOM
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

      DIMENSION CAPIM(NSTRT*2,NSTRT*2,2), FMOM(NSTRT*2,NSTRT*2,2)
      LOGICAL MULTIM

      CALL QENTER('ADDCAPIM')

      IF (IPRCAP.GE.5) THEN
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

C     ... multiply CAPIM with i
      IF (MULTIM) THEN
        DO I = 1, NSTRT*2
        DO J = 1, NSTRT*2
C         i.(a+ib) = -b + i.a
          TEMP = CAPIM(I,J,2)
          CAPIM(I,J,2) = CAPIM(I,J,1)
          CAPIM(I,J,1) = - TEMP
        ENDDO
        ENDDO
        IF (IPRCAP.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Modified i.CAPIM ***'
           CALL PRQMAT(CAPIM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &                IPQTOQ(1,0),LUPRI)
        ENDIF
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'ADDCAPIM: Warning - no multiplication with imaginary unit !'
      ENDIF

C      ...Add CAPIM into FMOM
      DO I = 1, NSTRT*2
      DO J = 1, NSTRT*2
       FMOM(I,J,1) = FMOM(I,J,1) + CAPIM(I,J,1)
       FMOM(I,J,2) = FMOM(I,J,2) + CAPIM(I,J,2)
      ENDDO
      ENDDO

      IF (IPRCAP.GE.5) THEN
       WRITE(LUPRI,'(/3X,A)')
     & 'ADDCAPIM:  *** Core Fock matrix (molfdir format)'//
     & ' with i.CAP_IM added (NZ=2) ***'
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
C  Save CAPR,CAPI integrals into the AOPROPER file.
C
C  Written by Miro ILIAS, Tel Aviv, febr. 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbham.h"
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
CC    IF (LEVYLE.OR.BSS) THEN
      IF (LEVYLE) THEN
        COMB = 'TFFF'
      ELSE
        COMB = 'TFFT'
      ENDIF

      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
      LABEL = 'CAP_RE  '
      CALL WRTPRO(CAPR,NNBBASX,LABEL,RTNLBL,IPRINT)

      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
C     COMB = 'TFFT'
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
C   Add cap operators into integral list.
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
       CALL QUIT('fu AONUM: entering AO does not fit !!')
      ENDIF

      RETURN
      END

