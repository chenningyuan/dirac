!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C  /* Deck capinp */
      SUBROUTINE CAPINP(WORD,RESET)
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
C
      SAVE SET
      DATA TABLE /'.INPQD ','.INPCAP','.PRINT ','.OVLCAP',
     &            '.NOOVLP','.NONREL','.NOREAL','.NOIMAG',
     &            '.IMASRE','.DIPVEL','.CAPNEW','.DC    ',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
!#include "memint.h"
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
      NOOVLP  = .FALSE.
      DONRCAP = .FALSE.
CMI   IF (LEVYLE) DONRCAP = .TRUE.
      NOCAPRE = .FALSE.
      NOCAPIM = .FALSE.
      CAPIMRE = .FALSE.
      CAPDIPVEL = .FALSE.    
      CAPNEW = .FALSE.
      DCADD = .FALSE.
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
     &                      5, 6, 7, 8, 9, 10, 11, 12), I
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
C&&& Do not use real part of CAP in FSCC
                 NOCAPRE = .TRUE.
                 GO TO 100
 8               CONTINUE
C&&& Do not use imaginary part of the CAP operator in FSCC
                 NOCAPIM = .TRUE.
                 GO TO 100
 9               CONTINUE
C&&& Use imaginary part of CAP as REAL (in FSCC part)
                 CAPIMRE = .TRUE.
                 GO TO 100
 10              CONTINUE
C&&& Do the comparison between CAP and DIRAC dipole velocity integrals
                 CAPDIPVEL = .TRUE.
                 GO TO 100
 11              CONTINUE
C&&& New and more compact implementation of the CAP
                 CAPNEW = .TRUE.
                 GO TO 100
 12              CONTINUE
C&&& The DC field is on, read it's components
                  DCADD = .TRUE.
C              ... dc field strength components
                  READ(LUCMD,*) EPSDC(1), EPSDC(2), EPSDC(3)
C              ... dc field active coordinates
                  READ(LUCMD,*)  X0DC(1),  X0DC(2),  X0DC(3)
                 GO TO 100
 13              CONTINUE
C&&& XXXXXX
                 GO TO 100
 14              CONTINUE
C&&& XXXXXX
                 GO TO 100
 15              CONTINUE
C&&& XXXXXX
                 GO TO 100
 16              CONTINUE
C&&& XXXXXX
                 GO TO 100
 17              CONTINUE
C&&& XXXXXX
                 GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized in CAPINP.'
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
         WRITE(LUPRI,'(A)')
     &  ' * four-component CAP integrals  according to',
     &  '    Miro, Sajeev ....         '
        ENDIF

        IF (DONRCAP.AND.CAPDIPVEL) THEN
         WRITE(LUPRI,'(2X,A)')
     &  '* Can not compare CAP/DIRAC dipole velocity integrals in'//
     &  ' the nonrelativistic mode.'
          CAPDIPVEL = .FALSE.
        ENDIF

        IF (CAPDIPVEL) THEN
         WRITE(LUPRI,'(2X,A)')
     &  '* Comparison between CAP/DIRAC dipole velocity integrals in'//
     &  ' the relativistic mode (then calculation).'
          IF (DCADD) THEN
             WRITE(LUPRI,'(2X,A)')
     &  '...no electric field in this case'
          ENDIF
          DCADD = .FALSE.
          DO I = 1, 3
           IF (GINTPT(I,1).GT.0.0D0) THEN
             WRITE(LUPRI,'(2X,A,D12.5)')
     &       '...integration interval must start from ZERO !'//
     &       ' current value ',GINTPT(I,1)
             CALL QUIT('CAPINP: integration interval'//
     &                 ' must start from ZERO !')
           ENDIF
          ENDDO
        ENDIF

        IF (NOCAPRE) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '* NOT using REAL part of the CAP in the correlation step.' 
        ENDIF

        IF (NOCAPIM) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '* NOT using IMAGINARY part of the CAP'//
     &   ' in the correlation step.' 
        ENDIF

        IF (CAPIMRE.AND..NOT.NOCAPIM) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '* Using IMAGINARY part of the CAP'//
     &   ' as REAL in the correlation step.' 
        ENDIF

        IF (CAPIMRE.AND..NOT.NOCAPIM.AND..NOT.NOCAPRE) THEN
         WRITE(LUPRI,'(2X,A)')
     &   '* Using IMAGINARY part of the CAP as REAL together with'//
     &   ' the ordinary REAL part in the correlation step.' 
        ENDIF

        IF (CAPNEW) THEN
         WRITE(LUPRI,'(2X,A)')
     &'* New implementation of the CAP in progress...'
        ENDIF

        WRITE(LUPRI,'(2X,A,F14.5)')
     &  '* read theta (PW) parameter (former INPQD file) = ',PW
        WRITE(LUPRI,'(2X,A)') 
     & '---   read integration parameters:   ---'
CMI    ...  write out integration parameters
        DO I = 1, 3
          WRITE(LUPRI,'(2X,A,F14.5,A,F14.5)')
     &    '* XO=',X0(I),' LAMDA=',LAMDA(I)
          WRITE(LUPRI,'(2X,A,I3)')
     &    '* integration subintervals, #',NGR(I)
          WRITE(LUPRI,*) (GINTPT(I,J),J=1,NGR(I))
          IF (NGR(I).GT.1) THEN
           WRITE(LUPRI,'(2X,A)')
     &     '* number of integration points per interval:'
           WRITE(LUPRI,*) (NGRS(I,J),J=1,NGR(I)-1)
          ELSE
           CALL QUIT('must have least two subintervals...')
          ENDIF
        ENDDO
        IF (DCADD) THEN
C       ... print out DC field parameters
          WRITE(LUPRI,'(2X,A)')
     &   '* DC field is added to CAP with these parameters:'
          WRITE(LUPRI,'(2X,A,3F12.6)')
     &    '...   field strenghts (x,y,z), in a.u.:',
     &    EPSDC(1),EPSDC(2),EPSDC(3)
          WRITE(LUPRI,'(2X,A,3F12.6)')
     &    '...active coordinates (x,y,z), in a.u.:',
     &    X0DC(1),X0DC(2),X0DC(3)
          DO I=1,3
            IF (X0DC(I).LT.GINTPT(I,1).AND.EPSDC(I).GT.0.0D0) THEN
              WRITE(LUPRI,*)
     &       'Nonzero DC field must be active IN the CAP space'
              WRITE(LUPRI,'(2X,2(A,F12.6))')
     &         'DC field is starting from:',X0DC(I),
     &         ' while CAP is staring from:',GINTPT(I,1)
              CALL QUIT('increase the cube for Vdc !')
            ENDIF
          ENDDO
          WRITE(LUPRI,'(/2X,A)')
     &   '* WARNING: x,y,z components of the DC field break'//
     &   ' the overall symmetry as CAP is totally symmetric operator ! '
          WRITE(LUPRI,'(2X,A/)')
     &   '...for this you have to use the lowest C1 symmetry...'//
     &   ' '
        ENDIF
        IF (NOOVLP) THEN
         WRITE(LUPRI,'(2X,A)') 
     &   '* no comparison between CAP-DIRAC overlap integrals  !'
        ENDIF
       ELSE
         WRITE(LUPRI,'(2X,A)') 
     &   '* only comparison between CAP and DIRAC overlap integrals.'
         DCADD = .FALSE.
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
C
C  Prepare arrays as parameters for the core CAP routine
C and run the CAP routine which prepares PEMAT binary file.
C
C
C  Called from F90 routine CAP_STUFF
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbcap.h"
#include "dummy.h"

#include "mxcent.h"
#include "nuclei.h"

      DIMENSION WORK(LWORK)

      CALL QENTER('DOSAJEEV')

      write(LUPRI,*) 'DOSAJEEV: LWORK=',LWORK

      write(LUPRI,*) 'dosajeev NPBAS,NBASIS=',NPBAS,NBASIS

#include "memint.h"

      NPBAS1 = NPBAS + NPBAS
C     ... basis_n
      CALL MEMGET('INTE',KBASIS_N,NBASIS,WORK,KFREE,LFREE)
C     ... basis_ord
      CALL MEMGET('INTE',KBASIS_ORD,NBASIS,WORK,KFREE,LFREE)
C     ... ang_pow
C     CALL MEMGET('INTE',KANG_POW,NPBAS*3,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KANG_POW,NPBAS1*3,WORK,KFREE,LFREE)
C     ... or_expon
C     CALL MEMGET('REAL',KORDEXPON,NPBAS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORDEXPON,NPBAS1,WORK,KFREE,LFREE)
C     ... ord_coe
C     CALL MEMGET('REAL',KORDCOE,NPBAS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORDCOE,NPBAS1,WORK,KFREE,LFREE)
C     .. R_ms
C     CALL MEMGET('REAL',KRMS,NPBAS*3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRMS,NPBAS1*3,WORK,KFREE,LFREE)

      IF (CAPOVLP) THEN

       CALL MEMGET('REAL',KPE_AOr,NNBBASX,WORK,KFREE,LFREE)

       CALL MEMCHK('dosajeev 0.1',WORK,1)
       write(LUPRI,*) 'dosajeev KFREE,LFREE=',KFREE,LFREE

       CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),WORK(KBASIS_ORD),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr),DUMMY,        
     &            DUMMY,DUMMY, DUMMY,DUMMY,DUMMY,DUMMY,
CMI  &            NBASIS2,NANGPOW)
     &            WORK(KFREE),LFREE,NBASIS,NPBAS1)

       write(LUPRI,*) 'dosajeev 1 KFREE,LFREE=',KFREE,LFREE

       CALL MEMCHK('dosajeev 1',WORK,1)

      ELSE
       IF (DONRCAP) THEN
CMI    ISIZE5=25000
CMI    CALL MEMGET('REAL',KPE_AOr,ISIZE5,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPE_AOr,NNBBASX,WORK,KFREE,LFREE)
CMI    CALL MEMGET('REAL',KPE_AOi,ISIZE5,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPE_AOi,NNBBASX,WORK,KFREE,LFREE)

        CALL MEMCHK('dosajeev 2',WORK,1)

        CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),WORK(KBASIS_ORD),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr),WORK(KPE_AOi),
     &            DUMMY,DUMMY, DUMMY,DUMMY, DUMMY,DUMMY,
C    &            NBASIS2,NANGPOW)
     &            WORK(KFREE),LFREE,NBASIS,NPBAS1)

        CALL MEMCHK('dosajeev 3',WORK,1)

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

        CALL MEMCHK('dosajeev 4',WORK,1)

        IF (DCADD) THEN

          CALL MEMGET('REAL',KPE_AOr,NNBBASX,WORK,KFREE,LFREE)

          CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),WORK(KBASIS_ORD),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            WORK(KPE_AOr), DUMMY, 
     &            WORK(KPEX_AOr),WORK(KPEX_AOi),
     &            WORK(KPEY_AOr),WORK(KPEY_AOi),
     &            WORK(KPEZ_AOr),WORK(KPEZ_AOi),
C    &           NBASIS2,NANGPOW)
     &            WORK(KFREE),LFREE,NBASIS,NPBAS1)

        ELSE
          CALL SAJEEV(WORK(KBASIS_N),WORK(KANG_POW),WORK(KBASIS_ORD),
     &            WORK(KORDEXPON),WORK(KORDCOE),WORK(KRMS),
     &            DUMMY, DUMMY,     
     &            WORK(KPEX_AOr),WORK(KPEX_AOi),
     &            WORK(KPEY_AOr),WORK(KPEY_AOi),
     &            WORK(KPEZ_AOr),WORK(KPEZ_AOi),
C    &           NBASIS2,NANGPOW)
     &            WORK(KFREE),LFREE,NBASIS,NPBAS1)
        ENDIF

        CALL MEMCHK('dosajeev 5',WORK,1)

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

      write(LUPRI,*) 'CAPADAPT: LWORK=',LWORK

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
        IF (.NOT.CAPDIPVEL) THEN
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

          IF (DCADD) THEN ! Vdc contribution to the relativistic CAP_d
          CALL CAPREAD('CAPD_VDC',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &                  ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)
          ENDIF
        ELSE
C         ... CAP relativistic dipole velocity integrals
          CALL CAPREAD('XCAPDVE ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('YCAPDVE ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

          CALL CAPREAD('ZCAPDVE ',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &         ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)

        ENDIF

       ENDIF

      ELSE
C      ... overlap integrals
        CALL CAPREAD('CAP_OVLP',WORK(KPMR),WORK(KPMRSA),WORK(KPMRSAT),
     &                ICAPBFO,IOBF1_2,WORK(KFREE),LFREE)
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

      LUCD = -1
      CALL GPOPEN(LUCD,'CAPDIR.inp','OLD',' ','FORMATTED',
     &       IDUMMY,.FALSE.)

      LUCL = -1
      CALL GPOPEN(LUCL,'CAPLBL.inp','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      IF (IPRCAP.GE.3) THEN
       CALL HEADER('Output from CAPADAPT',-1)
       write(LUPRI,'(2X,A,I4/)')
     & 'Total number of atomic basis functions, NBASIS = ',NBASIS
      ENDIF

      !return

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
C         ...    printout of IOBF parameters 
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
          CALL FLSHFO(LUPRI)
        ENDIF

      ENDDO

      !return

      CALL GPCLOSE(LUCD,'KEEP')
      CALL GPCLOSE(LUCL,'KEEP')

C=======================================================================
C     ... according to the DIRAC basis function order fill ICAPBFO
C=======================================================================
      DO I = 1, NBASIS
        IX = 0
        DO J = 1 , NBASIS
          IF (IOBF1(1,I).EQ.IOBF2(1,J).AND.
     &        IOBF1(2,I).EQ.IOBF2(2,J).AND.
     &        IOBF1(3,I).EQ.IOBF2(3,J)) IX = J        
        ENDDO

        IF (IX.EQ.0) THEN
          WRITE(LUPRI,*) 'CAPADAPT: I,J =',I,J
          WRITE(LUPRI,*) 'IOBF1(1,I),IOBF2(1,J):',IOBF1(1,I),IOBF2(1,J)
          WRITE(LUPRI,*) 'IOBF1(2,I),IOBF2(2,J):',IOBF1(2,I),IOBF2(2,J)
          WRITE(LUPRI,*) 'IOBF1(3,I),IOBF2(3,J):',IOBF1(3,I),IOBF2(3,J)
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
       CALL FLSHFO(LUPRI)
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
C  delete LS,SL blocks of the matrix ! 
C
C  basic functions are in basis DIRAC order with L functions first, than S functions
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
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
         IF (LEVYLE.and.i.gt.NLARGE.and.j.gt.NLARGE) then
           V(i,j) = 0.0D0
         ENDIF
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
C   Purpose: get CAP operator in MO basis into CAPI(*)  - output CAPI1 array 
C
C   First get it (using PRPMAO) as the real operator
C
C   IMAGPART - if true, read imaginary part of Ham. otherwise real
C
C   Called from  CAPMFDADAPT (prepcap.F)
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
     &               ICMOQ,NORB,WORK,KFREE,LFREE,IPRCAP)

      IF (IPRCAP.GE.5) THEN
        CALL FLSHFO(LUPRI)
        WRITE(LUPRI,'(/2X,2A)')
     & 'GETCAPMO: Get MO matrix of property (to be readjusted): ',
     &  PRPNAM(INDXPR)
        CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

C .... for the relativistic case add Vdc operator (if wished)
      IF (.NOT.IMAGPART.AND..NOT.DONRCAP.AND.
     &                      .NOT.CAPDIPVEL.AND.DCADD) THEN

       INDXPR_OLD = INDXPR

       INDXPR = IPCAP_VDC
       ISYM  = IPRPSYM(INDXPR)
       IREP = ISYM - 1
       IOPSY = JBTOF(ISYM-1,1)
       CALL DZERO(CAPI1,N2ORBXQ)
       CALL PRPMAT(INDXPR,IOPSY,CAPI1,WORK(KCMO),WORK(KBEIG),
     &               ICMOQ,NORB,WORK,KFREE,LFREE,IPRCAP)

      IF (IPRCAP.GE.5) THEN
        CALL FLSHFO(LUPRI)
        WRITE(LUPRI,'(/2X,2A)')
     & 'GETCAPMO: Get MO matrix of scalar property Vdc'//
     & ' (to be readjusted): ',
     &  PRPNAM(INDXPR)
        CALL PRQMAT(CAPI1,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

C      ... add CAPI1 to CAPI
       CALL DAXPY(NORBT*NORBT*NZ,1.0D0,CAPI1,1,CAPI,1)

       IF (IPRCAP.GE.5) THEN
         CALL FLSHFO(LUPRI)
         WRITE(LUPRI,'(/2X,3A)')
     &  'GETCAPMO: Get MO matrix of scalar property Vdc+CAP_d'//
     &  ' (to be readjusted): ',
     &  PRPNAM(INDXPR_OLD),PRPNAM(INDXPR)
        CALL PRQMAT(CAPI,NORBT,NORBT,NORBT,
     &         NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
       ENDIF

C      ... all we muest have in CAPI
C      CALL DCOPY(NORBT*NORBT*NZ,CAPI1,1,CAPI,1)

      ENDIF

      CALL MEMCHK('getcapmo 2',WORK,1)

CMI   ...MOVE ee-part to the upper left part of the MO-matrix
C     REARR = .FALSE.
      REARR = .true.
C
C     IF (REARR.AND..NOT.LEVYLE) THEN
      IF (REARR.AND..NOT.LEVYLE.AND..NOT.TWOCOMP) THEN

        CALL DZERO(CAPI1,NORBT*NORBT*NZ)
C
C       CAPI -> CAPI1
C
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

      ELSE
C   ... no rearrangement, simply copy
       IF (IPRCAP.GE.5) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GETCAPMO: No rearrangement of the matrix.'
        CALL FLSHFO(LUPRI)
       ENDIF
       CALL DCOPY(NORBT*NORBT*NZ,CAPI,1,CAPI1,1)
      ENDIF

      IF (IPRCAP.GE.5) THEN
         IF (.NOT.IMAGPART.AND..NOT.DONRCAP.AND.
     &                      .NOT.CAPDIPVEL.AND.DCADD) THEN
        WRITE(LUPRI,'(/2X,2A,A,A)')
     &  'GETCAPMO: Final MO matrix of property (output):',
     &    PRPNAM(INDXPR_OLD),'+',PRPNAM(INDXPR)
         ELSE
        WRITE(LUPRI,'(/2X,2A)')
     &  'GETCAPMO: Final MO matrix of property (output):',
     &  PRPNAM(INDXPR)
         ENDIF

        CALL PRQMAT(CAPI1,NORBT,NORBT,NORBT,
     &        NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
      ENDIF

C     CALL MEMCHK('getcapmo 3',WORK,1)

      CALL MEMREL('GETCAPMO',WORK,KFRSAV,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETCAPMO')
      RETURN
      END 

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck addcapmom */
      SUBROUTINE ADDCAPMOM(NSTRT,CAPMOM,FMOM,MULTIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: add MOLFDIR like operator CAPMOM to the FMOM
C
C  On input: CAPMOM - to be added to FMOM
C            MULTIM - flag to decide the multiplication by the imag.unit
C
C  Called from CAPMFADAPT (prepcap.F).
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
        WRITE(LUPRI,'(2X,A)')
     &  'ADDCAPMOM: CAP_MO with imaginary unit to get CAP_IM.'
        IF (IPRCAP.GE.5) THEN
         WRITE(LUPRI,'(/3X,A)')
     &   'ADDCAPMOM:  *** Modified i.CAPMOM - mult.by "i"-unit ***'
         CALL PRQMAT(CAPMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(2,2),
     &               IPQTOQ(1,0),LUPRI)
        ENDIF
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'ADDCAPMOM: No multiplication of CAP_MO with imaginary unit.'
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
C Adapt of the Fock matrix for CAP-FSCCSD calculations. Add both CAP real
C and imaginary contributions.
C
C Input: CAP - logical
C
C Output: FMOM with CAP (re,im) parts added
C
C ===  Called from MRCON1/traone.F   ===
C
C Written by Miro ILIAS, May 2007, Tel Aviv
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
     &   '*** Core Fock matrix (in DIRAC format) ***'
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &                      IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/3X,A,I2/)')
     &   '*** Core Fock matrix (in MOLFDIR format) ***'
        CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MIN0(NZ,2),
     &               IPQTOQ(1,0),LUPRI)
       ENDIF

#include "memint.h"
C         ... allocate space for CAP_IM
       CALL MEMGET('REAL',KCAPI, NORBT*NORBT*NZ,WORK,KFREE,LFREE)
       CALL MEMGET('REAL',KCAPI1,NORBT*NORBT*NZ,WORK,KFREE,LFREE)

C       ... allocate space for MOLFDIR like CAP
       CALL MEMGET('REAL',KCAPIM,NSTRT*NSTRT*8,WORK,KFREE,LFREE)

       IF (.NOT.NOCAPRE) THEN
C       ... get CAP  operator - real part (nonrel/rel) in the MO basis
        CALL GETCAPMO(.FALSE.,WORK(KCAPI),WORK(KCAPI1),WORK,KFREE,LFREE)
C       ... turn it in the MOLFDIR like form
        CALL MOLTODIR(WORK(KCAPI1),WORK(KCAPIM),
     &                NSTR,NSTRT,.TRUE.,IPRCAP)
C       ... add MOLFDIR like CAP_IM, WORK(KCAPIM), to the imaginary part of FMOM (if .true.)
        CALL ADDCAPMOM(NSTRT,WORK(KCAPIM),FMOM,.FALSE.)
       ENDIF

       IF (.NOT.NOCAPIM) THEN
C      .... get CAP nonre/rel (as a "real" operator) - IMAGINARY part - in the MO basis
        CALL GETCAPMO(.TRUE.,WORK(KCAPI),WORK(KCAPI1),WORK,KFREE,LFREE)
C       ... turn it in the MOLFDIR like form
        CALL MOLTODIR(WORK(KCAPI1),WORK(KCAPIM),NSTR,NSTRT,
     &                .NOT.CAPIMRE,IPRCAP)
C      ... add MOLFDIR like CAP_IM, WORK(KCAPIM), to the imaginary part of FMOM (if .true.)
        CALL ADDCAPMOM(NSTRT,WORK(KCAPIM),FMOM,.TRUE.)
       ENDIF

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
     &  'CAPMFDADAPT: *** Core Fock matrix with complete CAP(re,im) '//
     &  'added (MOLFDIR format) ***'
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
C             XCAPD1R - real and imaginary integrals for relativistic calculations
C             XCAPD1I
C             YCAPD1R
C             YCAPD1I
C             ZCAPD1R
C             ZCAPD1I
C
C             CAPD_VDC - real scalar Vdc contribution to the relativistic CAP operator
C
C           ... dipole velocity integrals calculated by CAP
C             XCAPDVE
C             YCAPDVE
C             ZCAPDVE
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
      LOGICAL ISFILE, DOXPRIND1,DOXPRIND2, DOXPRIND3
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
     &    LBL.EQ.'CAPD_VDC'.OR.
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
             IF (LBL.EQ.'XCAPD1R'.OR. LBL.EQ.'XCAPD1I'.OR.
     &           LBL.EQ.'XCAPDVE') THEN
          INTREP = ISYMAX(1,1) + 1
        ELSE IF (LBL.EQ.'YCAPD1R'.OR. LBL.EQ.'YCAPD1I'.OR.
     &           LBL.EQ.'YCAPDVE') THEN
          INTREP = ISYMAX(2,1) + 1
        ELSE IF (LBL.EQ.'ZCAPD1R'.OR. LBL.EQ.'ZCAPD1I'.OR.
     &            LBL.EQ.'ZCAPDVE') THEN
          INTREP = ISYMAX(3,1) + 1
        ELSE
          WRITE(LUPRI,*) 'unidentified SAVECAP label, LBL=',LBL
          CALL QUIT('SAVECAP: Unidentified label !')
        ENDIF
        IF (LBL.EQ.'ZCAPD1I') DOXPRIND2=.TRUE.
      ENDIF

C     ... label for AOPROPER ...
      LABEL = LBL

      WRITE(RTNLBL(2)(3:4),'(I2)') INTREP

      RTNLBL(2)(5:8) = COMB

C     ... write integrals according to each individual label
      CALL WRTPRO(CAPINTS,NNBBASX,LABEL,RTNLBL,IPRCAP)

C     ... control printout ...
      IF (IPRCAP.GE.0) THEN
        WRITE(LUPRI,'(2X,A,A8,A,A8,A,A)')
     &  'SAVECAP: CAP integrals written to the AOPROPER'//
     &  ' file under the label: ',
     &  RTNLBL(1),' ',RTNLBL(2),' ',LABEL
      ENDIF

      IF (DOXPRIND1) THEN
C
C       CAP (nr_re,nr_im,overlap,Vdc) integrals of IPTYP=1
C     =============================================================
C
      PLABEL(1) = LBL
      PNAME = LBL
      IPTYP            = 1
      NPCOMP           = 1
      PFAC(1)          = D1
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &      INDXPR,ISYXPR,ITRXPR,IPRCAP)

C             ... assign indexes
C ---------------------------------------------------
      IF      (LBL.EQ.'CAP_RE  ') THEN
        IPCAPR    = INDXPR
      ELSE IF (LBL.EQ.'CAP_IM  ') THEN
        IPCAPI    = INDXPR
      ELSE IF (LBL.EQ.'CAP_OVLP') THEN
        IPCAP_S   = INDXPR
      ELSE IF (LBL.EQ.'CAPD_VDC') THEN
        IPCAP_VDC   = INDXPR
      ELSE   
        write(LUPRI,*) 'SAVECAP: LBL=',LBL
        CALL QUIT('SAVECAP: Can not identify INDXPR label !')
      ENDIF
      IF (IPRCAP.GE.2) THEN
        WRITE(LUPRI,'(2X,A,A,A,I3)')
     &  '> operator ',LBL,' added to the DIRAC operator'//
     &  ' list, index=',
     &  INDXPR
        CALL WRIXPR(1,INDXPR)
      ENDIF

      ENDIF

C *************************************************
C        ... last two labels ...
C *************************************************
      IF (DOXPRIND2) THEN
C     ... last 2 XPRIND

C     CVALX = 0.50D0*CVAL ! inly for checking ...
      CVALX = CVAL

      PNAME    = 'CAP_D_Re '
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = -CVALX
      PFAC(2)  = -CVALX
      PFAC(3)  = -CVALX
      PLABEL(1) = 'XCAPD1R'
      PLABEL(2) = 'YCAPD1R'
      PLABEL(3) = 'ZCAPD1R'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,IPRCAP)
      IPCAP_DR = INDXPR

      PNAME    = 'CAP_D_Im  '
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = -CVALX
      PFAC(2)  = -CVALX
      PFAC(3)  = -CVALX
      PLABEL(1) = 'XCAPD1I'
      PLABEL(2) = 'YCAPD1I'
      PLABEL(3) = 'ZCAPD1I'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,IPRCAP)
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

      IF (LBL.EQ.'ZCAPDVE') THEN

      PNAME     = 'CAP_dip_ve'
      IPTYP     = 8
      NPCOMP    = 3
      PFAC(1)   = -CVAL
      PFAC(2)   = -CVAL
      PFAC(3)   = -CVAL
      PLABEL(1) = 'XCAPDVE'
      PLABEL(2) = 'YCAPDVE'
      PLABEL(3) = 'ZCAPDVE'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,IPRCAP)
      IPCAP_DV = INDXPR

      IF (IPRCAP.GE.2) THEN
        WRITE(LUPRI,'(2X,A,I3)')
     &  '* CAP_dip_ve operator added to the list, index=',IPCAP_DV
        CALL WRIXPR(1,IPCAP_DV)
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
C  Purpose:
C=============
C   read integrals from the PEMAT file under the LBL label into the CAP array
C
C   do the symmetry adaptation into the CAPSA array
C   
C   make triangular matrix of it and SAVE integrals on the AOPROPER file   
C   
C  On iput:
C =============   LBL
C                 CAP, CAPSA, CAPSAT
C                 ICAPBFO
C                 IOBF1_2  
C   
C   Called from prepcap.F/CAPADAPT
C
C   Written by Miro ILIAS, May 2007, Tel Aviv
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
      LOGICAL FNDLAB,SYMMETRIC, LSDEL
      CHARACTER*8 LBL

      CALL QENTER('CAPREAD')

      CALL DZERO(CAP,N2BBASX)
      CALL DZERO(CAPSA,N2BBASX)
      CALL DZERO(CAPSAT,NNBBASX)

C    ... read integrals with label LBL from PEMAT into CAPSA
      REWIND(LUPEMAT)
      IF (FNDLAB(LBL,LUPEMAT)) THEN 
        CALL READT(LUPEMAT,NNBBASX,CAPSA)
      ELSE
        write(LUPRI,*) 'CAPREAD: LBL,LUPEMAT:',LBL,LUPEMAT
        CALL QUIT('CAPREAD: Label not found !')
      ENDIF

      IF (IPRCAP.GE.2) THEN
       WRITE(LUPRI,'(/2X,A,A,A)')
     & 'CAPADAPT: CAP integrals of the label ',
     & LBL,' succesfully read from PEMAT file into the CAPSA array.'
      ENDIF

C     ... evaluate symmetry of operator
      SYMMETRIC =  (LBL.EQ.'CAP_RE  '.OR. ! real part of CAPnr
     &              LBL.EQ.'CAP_IM  '.OR. ! imaginary part of CAPnr
!    &              LBL.EQ.'XCAPD1R '.OR. !  
!    &              LBL.EQ.'YCAPD1R '.OR. !  
!    &              LBL.EQ.'ZCAPD1R '.OR. !  
!    &              LBL.EQ.'XCAPD1I '.OR. !  
!    &              LBL.EQ.'YCAPD1I '.OR. !  
!    &              LBL.EQ.'ZCAPD1I '.OR. !  
     &              LBL.EQ.'CAPD_VDC'.OR. ! Vdc scalar field contribution
     &              LBL.EQ.'CAP_OVLP')    ! CAP overlap integrals


      IF (SYMMETRIC) THEN
        SGNFACT = 1.0D0
      ELSE
        SGNFACT = -1.0D0
      ENDIF

C      ... tell to delete LS blocks (for beeing carried with alpha operator)
      LSDEL =    (
     &              LBL.EQ.'XCAPD1R '.OR. !  
     &              LBL.EQ.'YCAPD1R '.OR. !  
     &              LBL.EQ.'ZCAPD1R '.OR. !  
     &              LBL.EQ.'XCAPD1I '.OR. !  
     &              LBL.EQ.'YCAPD1I '.OR. !  
     &              LBL.EQ.'ZCAPD1I '.OR. !  
     &              LBL.EQ.'XCAPDVE '.OR. !  
     &              LBL.EQ.'YCAPDVE '.OR. !  
     &              LBL.EQ.'ZCAPDVE ' )

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

C      ... consider changing sign only for CAP/DIRAC non-diagonal elements
        IF (I.EQ.J) THEN
          SGNFACT_SAVE = 1.0D00
        ELSE
          SGNFACT_SAVE = SGNFACT
        ENDIF

C       -------------------------------------------------
C        .. store  integrals in the DIRAC basic order ..
C       -------------------------------------------------

        indx_i =  NBASIS*(ICAPBFO(I)-1) + ICAPBFO(J) ! common ij index
        indx_j =  NBASIS*(ICAPBFO(J)-1) + ICAPBFO(I) ! common ji index

        CAP_FACT = GET_CAPAF_FAC(IOBF1_2(I))*GET_CAPAF_FAC(IOBF1_2(J))

        CAP(indx_i) = CAPSA(IX)*CAP_FACT
        CAP(indx_j) = CAPSA(IX)*CAP_FACT*SGNFACT_SAVE ! only for (j,i)
 
        IF (IPRCAP.GE.4) THEN
         WRITE(LUPRI,'(2X,A,1X,3(A,2I3,1X,D14.8))')
     &   LBL,'CAP:',I,J,CAPSA(IX),
     &   ' DIRACij:',ICAPBFO(I),ICAPBFO(J),CAP(indx_i),
     &   ' DIRACji:',ICAPBFO(J),ICAPBFO(I),CAP(indx_j)
        ENDIF

      ENDDO 
      ENDDO 

      IF (IPRCAP.GE.5) THEN
       CALL HEADER ('Full CAP integrals in basic DIRAC order'//
     & ' (NO cleared blocks) ',-1)
       WRITE(LUPRI,'(2X,A,A)') 'The label of integrals is:',LBL
       CALL PRQMAT(CAP,NBASIS,NBASIS,NBASIS,NBASIS,1,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

C   ... delete SL/LS blocks integrals ....
      IF (.NOT.LSDEL) THEN
        CALL DELETE_LS(CAP,NBASIS,IPRCAP,.TRUE.)
      ELSE
CMI  ... delete LL,SS blocks of integrals
        CALL DELETE_LS(CAP,NBASIS,IPRCAP,.FALSE.)
      ENDIF

      IF (IPRCAP.GE.4) THEN
       CALL HEADER ('CAPxx integrals in basic DIRAC order'//
     & ' (cleared blocks) ',-1)
       WRITE(LUPRI,'(2X,A,A)') 'The label of integrals is:',LBL
       WRITE(LUPRI,'(2X,A,L1)') 'Were off-diagonal LS blocks cleared ?',
     &                          LSDEL
       CALL PRQMAT(CAP,NBASIS,NBASIS,NBASIS,NBASIS,1,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

C   ...internal DIRAC symmetry adaptation into CAPSA
      IF (LBL.EQ.'CAP_RE  '.OR.
     &    LBL.EQ.'CAP_IM  '.OR.
     &    LBL.EQ.'CAPD_VDC'.OR.
!    &          LBL.EQ.'XCAPD1R '.OR. !  
!    &          LBL.EQ.'YCAPD1R '.OR. !  
!    &          LBL.EQ.'ZCAPD1R '.OR. !  
!    &          LBL.EQ.'XCAPD1I '.OR. !  
!    &          LBL.EQ.'YCAPD1I '.OR. !  
!    &          LBL.EQ.'ZCAPD1I '.OR. !  
     &    LBL.EQ.'CAP_OVLP') THEN
C       IREP=0
C       IREPD = IRQMAT(1,IREP)
C       IREPD = IRQMAT(JM4POS(1),IREP)
        IREPD = 0
      ELSE IF (LBL.EQ.'XCAPD1R '.OR.LBL.EQ.'XCAPD1I '.OR.
     &         LBL.EQ.'XCAPDVE ') THEN
        IREPD = ISYMAX(1,2)
C       IREPD = IRQMAT(JM4POS(4),IREP) 
C       IREPD = ISYMAX(1,2)
C       IREPD = 6
      ELSE IF (LBL.EQ.'YCAPD1R '.OR.LBL.EQ.'YCAPD1I '.OR.
     &         LBL.EQ.'YCAPDVE ') THEN
        IREPD = ISYMAX(2,2)
C       IREPD = IRQMAT(JM4POS(3),IREP)
C       IREPD = 5
      ELSE IF (LBL.EQ.'ZCAPD1R '.OR.LBL.EQ.'ZCAPD1I '.OR.
     &         LBL.EQ.'ZCAPDVE ') THEN
        IREPD = ISYMAX(3,2)
C       IREPD = IRQMAT(JM4POS(2),IREP)
C       IREPD = 3
      ELSE 
        WRITE(LUPRI,*) 'CAPREAD ... LBL=',LBL
        CALL QUIT('CAPREAD: wrong LBL !')
      ENDIF

      CALL DZERO(CAPSA,N2BBASX)

C   ... only full matrixes must be passed to the subroutine
      CALL MTAOSO(CAP,CAPSA,NBASIS,IREPD,IPRCAP)

      IF (IPRCAP.GE.4) THEN
       CALL HEADER ('CAPxx integrals in SO DIRAC order'//
     & ' (after MTAOSO) ',-1)
       WRITE(LUPRI,'(2X,A,A)')
     & 'The label of these integrals is:',LBL
       CALL PRQMAT(CAPSA,NBASIS,NBASIS,NBASIS,NBASIS,1,
     &             IPQTOQ(1,0),LUPRI)
      ENDIF

C   ... make triangular matrix of CAPSA into CAPSAT
      CALL DAMATR(NBASIS,CAPSA,NBASIS,CAPSAT)

C   ... save triangular CAPSAT of LBL to the AOPROPER file
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
C  Functions are only up to G
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
C   ... d shells ...
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
C   .... f shells ...
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
C    .... G shells
      ELSE IF (AO.EQ.'XXXX'.OR.AO.EQ.'g400') THEN
        AONUM = 21
      ELSE IF (AO.EQ.'XXXY'.OR.AO.EQ.'g310') THEN
        AONUM = 22
      ELSE IF (AO.EQ.'XXXZ'.OR.AO.EQ.'g301') THEN
        AONUM = 23
      ELSE IF (AO.EQ.'XXYY'.OR.AO.EQ.'g220') THEN
        AONUM = 24
      ELSE IF (AO.EQ.'XXYZ'.OR.AO.EQ.'g211') THEN
        AONUM = 25
      ELSE IF (AO.EQ.'XXZZ'.OR.AO.EQ.'g202') THEN
        AONUM = 26
      ELSE IF (AO.EQ.'XYYY'.OR.AO.EQ.'g130') THEN
        AONUM = 27
      ELSE IF (AO.EQ.'XYYZ'.OR.AO.EQ.'g121') THEN
        AONUM = 28
      ELSE IF (AO.EQ.'XZZZ'.OR.AO.EQ.'g103') THEN
        AONUM = 29
      ELSE IF (AO.EQ.'YYYY'.OR.AO.EQ.'g040') THEN
        AONUM = 30
      ELSE IF (AO.EQ.'YYYZ'.OR.AO.EQ.'g031') THEN
        AONUM = 31
      ELSE IF (AO.EQ.'YYZZ'.OR.AO.EQ.'g022') THEN
        AONUM = 32
      ELSE IF (AO.EQ.'YZZZ'.OR.AO.EQ.'g013') THEN
        AONUM = 33
      ELSE IF (AO.EQ.'XYZZ'.OR.AO.EQ.'g112') THEN
        AONUM = 34
      ELSE IF (AO.EQ.'ZZZZ'.OR.AO.EQ.'g004') THEN
        AONUM = 35
      ELSE
       write(LUPRI,*) '>>> fu AONUM: entering atomic orbital:',AO,
     & ' is not in the list (which is up only to G functions) !'       
       CALL QUIT('fu AONUM: entering AO does not fit; basis'//
     & ' set up to G functions only !')
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
     &            I.EQ.16.OR.I.EQ.18.OR.I.EQ.19.OR.
     &            I.EQ.28.OR.I.EQ.25.OR.I.EQ.34) THEN
            GET_CAPAF_FAC =  DSQRT(3.0D0)
C        XYZ
         ELSE IF (I.EQ.35.OR.I.EQ.30.OR.I.EQ.21) THEN
            GET_CAPAF_FAC =  DSQRT(105.0D0)
         ELSE IF (I.EQ.33.OR.I.EQ.31.OR.I.EQ.29.OR.I.EQ.23.OR.
     &            I.EQ.27.OR.I.EQ.22) THEN
            GET_CAPAF_FAC =  DSQRT(15.0D0)
         ELSE IF (I.EQ.32.OR.I.EQ.24.OR.I.EQ.26) THEN
            GET_CAPAF_FAC =  DSQRT(9.0D0)
         ELSE
            GET_CAPAF_FAC = 1.0D0
         ENDIF

      RETURN
      END

