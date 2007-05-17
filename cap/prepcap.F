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
     &            '.NOOVLP','.XXXXXX','.XXXXXX','.XXXXXX',
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
      NOOVLP = .FALSE.

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
        WRITE(LUPRI,'(2X,A,D25.10)') 'read INPQD (theta,PW) = ',PW
        WRITE(LUPRI,'(2X,A)') 
     & 'read integration parameters:'
CMI    ...  write out integration parameters
        DO I=1,3
          WRITE(LUPRI,*) X0(I),LAMDA(I)
          WRITE(LUPRI,*) NGR(I)
          WRITE(LUPRI,*) (GINTPT(I,J),J=1,NGR(I))
          WRITE(LUPRI,*) (NGRS(I,J),J=1,NGR(I)-1)
        ENDDO
        IF (NOOVLP) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '*No comparison of CAP-DIRAC overlap integrals !'
        ENDIF
       ELSE
         WRITE(LUPRI,'(2X,A)') 
     &   'Only comparison between CAP and DIRAC overlap integrals.'
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
C and run the CAP routine which prepares PEMAT/PEMAT_ovlp formatted files
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION WORK(LWORK)

      CALL QENTER('DOSAJEEV')

#include "memint.h"
!     CALL MEMCHK('CAPINTS',WORK,1)

C    ...  call F90 routine from the module cap
C       call cap_init

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
      INQUIRE(FILE='PEMAT_ovlp',EXIST=IS_PEMAT2)
      IF (CAPOVLP) THEN
        FILES2 = FILES2.AND.(IS_PEMAT2)
      ELSE
        FILES2 = FILES2.AND.(IS_PEMAT1)
      ENDIF

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
C   Adapt the CAP integrals from previously obtained formatted files 
C   PEMAT/PEMAT_ovlp for the DIRAC calculations.
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
      CHARACTER*10 PEMAT_NAME

      CALL QENTER('CAPADAPT')
  
C ======== first real part (CAPRE or CAP_OVLP)
#include "memint.h"

      IF (.NOT.CAPOVLP) THEN
         PEMAT_NAME="PEMAT"

         CALL MEMGET('REAL',KPMRSAT,NNBBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMRSAT),NNBBASX)
         CALL MEMGET('REAL',KPMRSA,N2BBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMRSA),N2BBASX)
         CALL MEMGET('REAL',KPMR,N2BBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMR),N2BBASX)
         CALL CAPREAD(PEMAT_NAME,WORK(KPMR),WORK(KPMRSA),
     &                 ICAPBFO,IOBF1_2,.TRUE.,IPRCAP)
C     Save HERMIT symmetry adapted integrals into the AOPROPER file.
         CALL DAMATR( NBASIS, WORK(KPMRSA), NBASIS, WORK(KPMRSAT) )
         CALL SAVECAP(WORK(KPMRSAT),1,WORK(KFREE),LFREE,IPRCAP)

         CALL DZERO(WORK(KPMR),N2BBASX)
         CALL DZERO(WORK(KPMRSA),N2BBASX)
         CALL DZERO(WORK(KPMRSAT),NNBBASX)

         CALL CAPREAD(PEMAT_NAME,WORK(KPMR),WORK(KPMRSA),
     &                 ICAPBFO,IOBF1_2,.FALSE.,IPRCAP)

         CALL DAMATR( NBASIS, WORK(KPMRSA), NBASIS, WORK(KPMRSAT) )
         CALL SAVECAP(WORK(KPMRSAT),2,WORK(KFREE),LFREE,IPRCAP)

      ELSE
         PEMAT_NAME="PEMAT_ovlp"

         CALL MEMGET('REAL',KPMRSAT,NNBBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMRSAT),NNBBASX)
         CALL MEMGET('REAL',KPMRSA,N2BBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMRSA),N2BBASX)
         CALL MEMGET('REAL',KPMR,N2BBASX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KPMR),N2BBASX)

         CALL CAPREAD(PEMAT_NAME,WORK(KPMR),WORK(KPMRSA),
     &                 ICAPBFO,IOBF1_2,.TRUE.,IPRCAP)

C     Save HERMIT symmetry adapted integrals into the AOPROPER file.
         CALL DAMATR( NBASIS, WORK(KPMRSA), NBASIS, WORK(KPMRSAT) )
         CALL SAVECAP(WORK(KPMRSAT),3,WORK(KFREE),LFREE,IPRCAP)

      ENDIF

      CALL MEMREL('CAPADAPT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('CAPADAPT')
      RETURN
      END 

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

      SUBROUTINE GET_DIR_CAP_AO_MATCH(ICAPBFO,IOBF1_2,
     &           IOBF1,IOBF2,ILC1,ILC2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C             Called from  prepcap90.F90/CAP_DIRAC_AO
C
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

      subroutine delete_ls(V,N,IPRINT)
!
! delete LS,SL blocks of the matrix ! basic functions are in basis DIRAC order
!
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      real*8 V(N,N)
      call qenter('delete_ls')

!     write(LUPRI,*) 'NBASIS,NLARGE,NSMALL:',NBASIS,NLARGE,NSMALL

      do i = 1 ,NBASIS
      do j = 1, NBASIS
       if (i.le.nlarge.and.j.gt.nlarge) then
          V(i,j) = 0.0D0
       endif
!      ------------------------------------------------------
       if (i.gt.nlarge.and.j.le.nlarge) then
          V(i,j) = 0.0D0
       endif
      enddo
      enddo

      call qexit('delete_ls')
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck getcapim */
      SUBROUTINE GETCAPI(CAPI,CAPI1,WORK,KFREE,LFREE)
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

      DIMENSION CAPI(NORBT,NORBT,NZ),CAPI1(NORBT,NORBT,NZ),WORK(*)
      LOGICAL REARR

      CALL QENTER('GETCAPI')

      KFRSAV = KFREE

      IPRINT = IPRCAP

      CALL DZERO(CAPI,NORBT*NORBT*NZ)

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

CMI   ...MOVE ee-part to the upper left part of the MO-matrix
C     REARR = .FALSE.
      REARR = .true.
C     IF (REARR.AND..NOT.LEVYLE) THEN

        CALL DZERO(CAPI1,NORBT*NORBT*NZ)

        IOFF_E = 0
        IOFF_EP = 0
        DO IFSYM =1 , NFSYM
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
       ENDDO

       IF (IPRINT.GE.5) THEN
         WRITE(LUPRI,'(2X,2A)')
     & 'GETCAPI: Rearranged MO matrix of property: ',PRPNAM(INDXPR)
        CALL PRQMAT(CAPI1,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
       ENDIF

C     ENDIF

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

C ... note: FMOM and CAPIM are in MOLFDIR format !
      IF (IPRCAP.GE.5) THEN
          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Entering core Fock matrix in MOLFDIR'//
     &    ' format, FMOM ***'
          CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &                IPQTOQ(1,0),LUPRI)

          WRITE(LUPRI,'(/3X,A)')
     &    'ADDCAPIM:  *** Entering (MOLFDIR like) CAPIM ***'
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
       TEMP1 = FMOM(I,J,1)
       FMOM(I,J,1) = TEMP1 + CAPIM(I,J,1)
       TEMP2 =  FMOM(I,J,2) 
       FMOM(I,J,2) = TEMP2 + CAPIM(I,J,2)
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
      SUBROUTINE SAVECAP(CAPINTS,ILABEL,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Save CAP integrals into the AOPROPER file and add corresponding operator to the
C the list according to the ILABEL 
C
C  ILABEL = 1 CAP_RE
C  ILABEL = 2 CAP_IM
C  ILABEL = 3 CAP_OVLP
C
C  Called from prepcap.F/CAPADAPT
C
C  Written by Miro ILIAS, Tel Aviv, febr. 2007
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
      DIMENSION CAPINTS(*),WORK(LWORK)
      LOGICAL ISFILE
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      CHARACTER PNAME*16, PLABEL(3)*8
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

C
C  Generate integral labels and write integrals into AOPROPER
C
CC    IF (LEVYLE.OR.BSS) THEN
      IF (LEVYLE) THEN
        COMB = 'TFFF'
      ELSE
        COMB = 'TFFT'
      ENDIF

      IF (ILABEL.EQ.1) THEN
        LABEL = 'CAP_RE  '
        PNAME = 'CAP - real part '
      ELSE IF (ILABEL.EQ.2) THEN
        LABEL = 'CAP_IM  '
        PNAME = 'CAP - imag part '
      ELSE IF (ILABEL.EQ.3) THEN
        LABEL = 'CAP_OVL '
        PNAME = 'CAP - overlap   '
      ENDIF

      PLABEL(1) = LABEL

      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
      CALL WRTPRO(CAPINTS,NNBBASX,LABEL,RTNLBL,IPRINT)

C
C        CAP integrals - re,im or overlap
C     ======================================
C
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &               INDXPR,ISYXPR,ITRXPR,WORK,LWORK,IPRINT)

      IF (ILABEL.EQ.1) THEN
        IPCAPR    = INDXPR
      ELSE IF (ILABEL.EQ.2) THEN
        IPCAPI    = INDXPR
      ELSE IF (ILABEL.EQ.3) THEN
        IPCAP_S   = INDXPR
      ENDIF

      IF (IPRINT.GE.1) THEN
        IF (ILABEL.EQ.1) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'SAVECAP: CAP_re integrals written into AOPROPER file.'
        ELSE IF (ILABEL.EQ.2) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'SAVECAP: CAP_im integrals written into AOPROPER file.'
        ELSE IF (ILABEL.EQ.3) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'SAVECAP: CAP_ovlp integrals written into AOPROPER file.'
        ENDIF
        IF (IPRINT.GE.2) THEN
           WRITE(LUPRI,'(2X,A,A,A,A4,A,I3)')
     &     ' Operator label (PLABEL):',PLABEL(1),
     &     ' COMB=',COMB,' operator index(INDXPR)=',INDXPR
        ENDIF
      ENDIF

      CALL QEXIT('SAVECAP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck capready */
      SUBROUTINE CAPREAD(PEMAT_NAME,CAP,CAPSA,
     &               ICAPBFO,IOBF1_2,REAL_PART,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    read integrals from the file named PEMAT_NAME into CAP 
C and do the symmetry adaptation intp CAPSA
C 
C  return CAPSA - CAP integrals (re,im,Sovlp) in Dirac symmetry adapted order 
C   depending of the "REAL_PART"
C 
C  Called from CAPADAPT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
      DIMENSION CAP(N2BBASX),CAPSA(N2BBASX),
     &          ICAPBFO(NBASIS),IOBF1_2(NBASIS)
      INTEGER CMNINDX
      LOGICAL REAL_PART
      CHARACTER*10 PEMAT_NAME
      COMPLEX*16 CAPINT
      CALL QENTER('CAPREAD')

      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A,A)') 'CAPADAPT: going to read the file ',
     &  PEMAT_NAME
      ENDIF

      LUPM = 0
      CALL GPOPEN(LUPM,PEMAT_NAME,'OLD',' ','FORMATTED',
     &       IDUMMY,.FALSE.)

C==================================================================
C     ... read PEMAT in the cycle and resort the matrix elements
C    from CAP order to the DIRAC basic order into WORK arrays
C==================================================================
      IX = 0
      DO I = 1, NBASIS
      DO J = 1, I
        IX = IX + 1
        READ(LUPM,*) II,JJ,CAPINT
        CAPRE = REAL(CAPINT)
        CAPIM = IMAG(CAPINT)
        IF (I.NE.II.OR.J.NE.JJ) THEN
         CALL QUIT(' wrong indexes in PEMAT !')
        ENDIF
        IX1 = CMNINDX(I,J)
        IF (IX.NE.IX1) THEN
         write(LUPRI,*) 'read II,JJ,CAPINT:',II,JJ,CAPINT
         write(LUPRI,*) 'I,J, common index ',I,J,IX1
         CALL QUIT('CAPADAPT: wrong IX !')
        ENDIF

C       -------------------------------------------------
C       .. store  integrals in the DIRAC basic order ..
C       -------------------------------------------------

         indx_i =  NBASIS*(ICAPBFO(I)-1) + ICAPBFO(J)  
         indx_j =  NBASIS*(ICAPBFO(J)-1) + ICAPBFO(I) 

C       CAP_FACT = GET_CAPAF_FAC(IOBF1(2,I))*GET_CAPAF_FAC(IOBF1(2,J))
        CAP_FACT = GET_CAPAF_FAC(IOBF1_2(I))*GET_CAPAF_FAC(IOBF1_2(J))

        IF (IPRINT.GE.10) THEN
          write(LUPRI,*) 'indx_i, indx_j, CAPRE, CAPRE*CAP_FACT:',
     &    indx_i, indx_j, CAPRE, CAPRE*CAP_FACT
        ENDIF

        IF (REAL_PART) THEN
          CAP(indx_i) = CAPRE*CAP_FACT
          CAP(indx_j) = CAPRE*CAP_FACT
        ELSE
          CAP(indx_i) = CAPIM*CAP_FACT
          CAP(indx_j) = CAPIM*CAP_FACT
        ENDIF

      ENDDO 
      ENDDO 

      CALL GPCLOSE(LUPM,'KEEP')

C! delete SL/LS blocks integrals
      call delete_ls(CAP,NBASIS,IPRCAP)

      IF (IPRINT.GE.4) THEN
       CALL HEADER ('CAPxx in basic DIRAC order',-1)
       CALL PRQMAT(CAP,NBASIS,NBASIS,NBASIS,NBASIS,1,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

! ...internal DIRAC symmetry adaptation
      CALL MTAOSO(CAP,CAPSA,NBASIS,0,IPRINT)

C     WRITE(LUPRI,'(/,5X,A/)') 'CAPADAPT: Symmetry adaptation:'

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


