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

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mp2noinp */
      SUBROUTINE MP2NOINP(WORD,RESET,WORK,LWORK)
C***********************************************************************
C
C     Input section for generation of MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (NTABLE = 8)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "../moltra/dcbtra.h"
#include "dcbmp2no.h"
C
      LOGICAL SET, NEWDEF, RESET,LBIT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
C  .MAXVIRT - maximum number of virtuals in MP2
C             if energy criterium gives more then max energy will be
C             reduced
C  .SCHEME  - transformation scheme used for MOLTRA (default taken from moltra)
C  .SELECT  - select min NO for active orbitals afterwards (relccsd,
C             lucita, kr-mcscf, luciarel, ...)
C             Default: e.g. "NO-occ 1.99 0.001 0.1" in analogy with the
C             "energy -10.0 20.0 0.2"          |--> 10%
C  .MULPOP  - perform Mulliken population analysis for MP2-NOs
C
C  FIXME more ideas for options: (18-Apr-2008)
C
C  .FOCKNO  - transform selected NO's to canonical Fock orbitals in
C             subspace (?? for better convergence in CI?)
C  
      SAVE SET
      DATA TABLE /'.ACTIVE','.PRINT ','.INTFLG','.F2CONT',
     &            '.MAX VS','.SCHEME','.SEL NO','.MULPOP'/
      DATA SET/.FALSE./
C
#include "ibtfun.h"
C
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
C     Initialize /CBIMP2NO/
C     =====================
C
      DEF_EN   = .TRUE.
      IPRMP2NO = 0
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      IGTINT = IBTAND(INTGEN/8,1)
      IFCMP2NO = 1
      IFEMP2NO = 1
C
      MP2NATPOP = .FALSE.
      MP2NATPOP2= .FALSE.
C     initialize string of active orbitals for 4-index/MP2 calculation
C     initialize string of active NO orbitals for subsequent correlation
C     step
      DO I = 1,NFSYM
        VECMP2NO(I) = 'energy -10.0 20.0 0.2 '
        VECNOSEL(I) = 'NO-occ 1.99999999 0.0001 0.1 '
      ENDDO
      IMAXVSP     = JMAXVSP
Caspg, default changed from scheme 4 to scheme 6
      IAOMOSCHEME = 6 
      
C
C     Process input from CBIMP2NO
C     =========================
C
      NEWDEF = (WORD .EQ. '*MP2 NO')
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
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in *MP2 NO.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *MP2 NO.')
    1          CONTINUE
C&&&& ACTIVE: Spinors to be included in the MP2 calculation
                  READ(LUCMD,'(A72)') (VECMP2NO(I),I=1,NFSYM)
                  DEF_EN = .FALSE.
               GO TO 100
    2          CONTINUE
C&&&& PRINT: General print level
                  READ(LUCMD,*) IPRMP2NO
               GO TO 100
    3          CONTINUE
C&&& INTFLG: Specify what two-itegrals should be included in this run
                  IF (IGTINT.EQ.1) THEN
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT,IGTINT
                  ELSE
                     READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  END IF
               GO TO 100
    4          CONTINUE
C&&& F2CONT: Specify Coulomb/exchange contributions
                  READ(LUCMD,*) IFCMP2NO,IFEMP2NO
               GO TO 100
    5          CONTINUE
C&&& MAX VS: Specify max. number of virt. spinors in NO-MP2 calculation
                  READ(LUCMD,*) IMAXVSP
               GO TO 100
    6          CONTINUE
C&&& SCHEME: Specify integral transformation scheme for MP2 (RELCCSD
C            module) calculation 
                  READ(LUCMD,*) IAOMOSCHEME
               GO TO 100
    7          CONTINUE
C&&& SEL NO: Select min NO for subsequent correlation calculation
                  READ(LUCMD,*) (VECNOSEL(I),I=1,NFSYM)
               GO TO 100
    8          CONTINUE
C&&& MULPOP: do Mulliken population analysis for MP2 NOs ('all' or 'active')
                  MP2NATPOP2 = .TRUE.
                  READ(LUCMD,*) SELPOPNAT
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in *MP2 NO.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *MP2 NO.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C
      INTMP2NO = ILLINT + 2*ISLINT + 4*ISSINT + 8*IGTINT
      IF2MP2NO = IFCMP2NO + 2*IFEMP2NO
C
C     Print section
C     =============
C
      IF(.NOT.DOMP2NO) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' *MP2 NO : MP2 natural orbitals'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I3)') ' * General print level:',IPRMP2NO
      WRITE(LUPRI,'(A)') ' * allowed active spinors in MP2-NO run:'
      DO I = 1,NFSYM
        WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &  '- spinors in fermion ircop ',FREP(I),' :',VECMP2NO(I)
      ENDDO
      WRITE(LUPRI,'(A)') ' * allowed range of NO spinors in 
     & correlation run:'
      DO I = 1,NFSYM
        WRITE(LUPRI,'(4X,A,A3,A,A72)')
     &  '- spinors in fermion ircop ',FREP(I),' :',VECNOSEL(I)
      ENDDO
      WRITE(LUPRI,'(A,I3)') ' * Integral transformation scheme: ',
     &                          IAOMOSCHEME
      WRITE(LUPRI,'(A,I3)') ' * Maximum number of virtual spinors: ',
     &                          IMAXVSP    
      IF(ILLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No LL contributions'
      IF(ISLINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SL contributions'
      IF(ISSINT.EQ.0) WRITE(LUPRI,'(A)') ' * No SS contributions'
      IF(IGTINT.EQ.0) WRITE(LUPRI,'(A)') ' * No GT contributions'
      IF( MP2NATPOP ) THEN
        WRITE(LUPRI,'(A)') ' * Run Mulliken population analysis'
        IF( SELPOPNAT .eq. 'all' .or. SELPOPNAT .eq. 'active' )THEN 
          WRITE(LUPRI,'(A,A6,A)') 
     &                     '   - print analysis for ',SELPOPNAT,' NOs'
        ELSE
          CALL Abend2( 'MP2NOINP: wrong keyword for Mull. pop. ana' )
        END IF
      END IF
      IF(IF2MP2NO.NE.3) THEN
        WRITE(LUPRI,'(A)')' * Two-electron contributions restricted to:'
        IF(LBIT(IF2MP2NO,1)) WRITE(LUPRI,'(3X,A)')
     &     '+ direct contribution'
        IF(LBIT(IF2MP2NO,2)) WRITE(LUPRI,'(3X,A)')
     &     '+ exchange contribution'
      ENDIF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mp2nogen */
      SUBROUTINE MP2NOGEN()
C***********************************************************************
C
C     Generate MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "aovec.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
#include "dcbfir.h"
#include "dcbgen.h"
#include "dcbmp2no.h"
#include "blocks.h"
C
      LOGICAL SAVEFLAGS(4), TOBE, FNDLAB
C
#include "memint.h"
      CALL TITLER('MP2 natural orbitals','*',125)
C
        WRITE(LUPRI,'(//A/A/A/)')
     &'     * The following modules will be invoked:    ',
     &'       - MOLTRA (integral transformation)        ',
     &'       - RELCCSD-MP2 to obtain MP2 density matrix'
C
C     run RELCCSD to get MP2dens
C
      RUN_CCMOD   = .TRUE.
      SKIP_MP2DEN = .FALSE.
      SKIP_ITRAFO = .FALSE.
C
C     test for availabilty of integral/MP2 nat.orb. transformation 
C     matrix
C
      INQUIRE(FILE='KRMCOLD',EXIST=TOBE)
      IF (TOBE) THEN
C       ... search for coefficients label
        CALL OPNFIL(LUKRMC,'KRMCOLD','UNKNOWN','MP2_NA')
        IF( FNDLAB('NEWNATOB',LUKRMC)) THEN 
           SKIP_MP2DEN = .TRUE.
           WRITE(LUPRI,'(//A/A)') 
     &'     ***** MP2 natural orbitals on file KRMCOLD present *****',
     &'     *****             using existing orbitals          *****'
        END IF
        CLOSE(LUKRMC,STATUS='KEEP')
      END IF
      INQUIRE(FILE='MDCINT',EXIST=TOBE)
      IF (TOBE) SKIP_ITRAFO = .TRUE.
C
      CALL MP2NOGE1()
C
      RUN_CCMOD  = .FALSE.
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mp2noge1 */
      SUBROUTINE MP2NOGE1()
C***********************************************************************
C     Generate MP2 natural orbitals
C
C     Written by S.Knecht and H.J.Aa. Jensen - March 2008
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DM1=-1.0D0,DM2=-2.0D0,D2=2.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmp2no.h"
#include "dcbham.h"
#include "infpar.h"
#if defined (VAR_MPI)
#include "mpif.h"
#endif
      LOGICAL FND,FNDLAB,TOBE
      DIMENSION IC(2),IE(2),ID(2), KVEC(2)
      INTEGER NR(2)
      CHARACTER CTMP*1, LABEL*8
      real(8), allocatable :: WORK(:)
      CALL QENTER('MP2NOGE1')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in MVOGE1')
C
      KFRSAV=KFREE
C
C     Initialization
C
      IMO_NAT     = NORBT*NESHT*NZ
      IOCC_NAT    = NESHT
      LUNATORB    = 98
C
      LUBUF = 22
C     full Fock matrix
      IFCMP2NO = 1
C     only direct (Coulomb) contributions 
      IF(IF2MP2NO.EQ.1) IFCMP2NO = 2
C     only exchange contributions
      IF(IF2MP2NO.EQ.2) IFCMP2NO = 3
      ISYMP2NO = 1
      IHRMP2NO = 1
C
C     Generate MP2NOs 
C
      IF( .NOT. SKIP_MP2DEN ) CALL DENMP2NO(WORK(KFREE),LFREE)
C
C     transf. matrix and NO occ. numbers on disk 
C
C     read / write / print section 
C     ------------
C
      IF( MYTID .eq. MPARID )THEN
C
C       memory allocation
C       -----------------
        CALL MEMGET('REAL',KOCC_NAT ,NESHT,WORK,KFREE,LFREE)
C
C       initialize
C       ----------
        CALL DZERO(WORK(KOCC_NAT),NESHT)
C
C       read occ. numbers and transf. matrix
C       ------------------------------------
C
        CALL OPNFIL(LUKRMC,'KRMCOLD','OLD','MP2NOG')
        CALL REAKRMC(LUKRMC,'MP2NOOCC',WORK(KOCC_NAT),NESHT)
C
        IF( .NOT. SKIP_MP2DEN ) THEN
C
           CALL MEMGET('REAL',KCMO_NAT ,NORBT*NESHT*NZ,WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KCMO_TRA ,NCMOTQ,WORK,KFREE,LFREE)
C          ... initialize
           CALL DZERO(WORK(KCMO_NAT),NORBT*NESHT*NZ)
           CALL DZERO(WORK(KCMO_TRA),NCMOTQ)
C          ... MP2 NO coefficients
           CALL REAKRMC(LUKRMC,'MP2NATOB',WORK(KCMO_NAT),NORBT*NESHT*NZ)
C          ... MO coefficients
           CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MP2NOG')
           CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO_TRA),XDUMMY,XDUMMY,
     &                 XDUMMY,2,WORK,KFREE,LFREE)
           CLOSE(LUCOEF,STATUS='KEEP')
C
C          transformation 
C          --------------
           CALL TDMP2NATMO(WORK(KCMO_TRA),WORK(KCMO_NAT),
     &                     WORK(KFREE),LFREE)
C
C          save to disk
C          ------------
           CALL WRTKRMC(LUKRMC,'NEWNATOB',WORK(KCMO_TRA),NCMOTQ)
C
C          release memory
C          --------------
           CALL MEMREL('MP2NOG.1',WORK,KFRSAV,KCMO_NAT,KFREE,LFREE)
C
        END IF
C
        CLOSE(LUKRMC,STATUS='KEEP')
C
C       print section 
C       -------------
        WRITE(LUPRI,'(//A/A/A)')
     &'   ************************************************************',
     &'   *************** MP2 natural orbital occupancies ************',
     &'   ************************************************************'
C
        IOFF2 = 1
C
        DO I = 1, NFSYM
           NOCCI = NESH(I)
           IF (NOCCI .EQ. 0) THEN
              WRITE(LUPRI,9001) I
           ELSE
              OCCSUM = DSUM(NOCCI,WORK(KOCC_NAT+IOFF2-1),1)
              WRITE(LUPRI,9002) I,(WORK(KOCC_NAT+IOFF2-2+J),J=1,NESH(I))
              WRITE(LUPRI,9003) OCCSUM
           END IF
           IOFF2 = IOFF2 + NESH(I)
        END DO
C
C       run Mulliken population analysis
C       --------------------------------
C
        IF( MP2NATPOP2 )THEN
C         count 'active NOs'
          IOFF2 = 1
          DO I = 1, NFSYM
            NACTMP2NO(I) = 0
            ICOUNT = 0
            NOCCI  = NESH(I)
            DO J = 1, NOCCI
              IF(WORK(KOCC_NAT+IOFF2-2+J) .ne. D0) ICOUNT = ICOUNT + 1
            END DO
            NACTMP2NO(I) = ICOUNT
            IOFF2 = IOFF2 + NESH(I)
          END DO
CSK       WRITE(LUPRI,*) ' active spinors for MULPOP:'
CSK       DO I = 1, NFSYM
CSK         WRITE(LUPRI,*) ' NACTMP2NO is',NACTMP2NO(I)
CSK       END DO
C         enable MP2-NO part in MULPOP
          MP2NATPOP = .TRUE.
          CALL MULPOP(WORK(KFREE),LFREE)
C
C         disable - we may want to run it again after a
C         successful MCSCF...
          MP2NATPOP = .FALSE.
        END IF
      END IF
C
C     pick out 'active' NO spinors according to 
C     list in VECNOSEL
C
      WRITE(LUPRI,'(//A/A/A)')
     &'   ************************************************************',
     &'   *************** selection of active NO spinors *************',
     &'   ************************************************************'
      NTOT  = 0
      IOFF2 = 1
      DO I = 1,NFSYM
        NR(I) = 1
        CALL MEMGET('INTE',KVEC(I),NESH(I),WORK,KFREE,LFREE)
        CALL NUMLS4(VECNOSEL(I),WORK(KVEC(I)),NESH(I),NR(I),
     &              NESH(I),NR(I),WORK(KOCC_NAT+IOFF2-1))
        IF( NR(I) .eq. 0 ) THEN
          WRITE(LUPRI,'(3X,A,I3)')
     &    '- No active spinors in fermion ircop ',I
        ELSE
          NTOT = NTOT + NR(I)
          CALL NOSELPRI(WORK(KVEC(I)),I,NR(I))
        ENDIF
        IOFF2 = IOFF2 + NESH(I)
      ENDDO
C
C     FIXME: - automatic 'active space selection for GAS-CI?'
C            - proper implementation of .FOCKNO option
C
C
C     memory deallocation
C     -------------------
      CALL MEMREL('MP2NOGE1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      call dealloc(WORK)
C
      CALL QEXIT('MP2NOGE1')
      RETURN
C
 9001 FORMAT(/'   Symmetry',I3,//,'   No occupied orbitals')
 9002 FORMAT(/'   Natural orbital occupation numbers, symmetry',
     &            I3,//,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck denmp2no */
      SUBROUTINE DENMP2NO(WORK,LWORK)
C***********************************************************************
C
C     Generate density matrix for chosen MP2NO vectors
C
C     1-/2-e integral transformation part based on routine PAMTRA
C
C     call RELCCSD-MP2 to obtain MP2 density matrix
C
C     written by S. Knecht - June 2008, Odense
C
C***********************************************************************
      use moltra_labeling
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmp2no.h"
#include "../moltra/dcbtra.h"
#include "../moltra/dcbtr3.h"
C
#include "maxorb.h"
#include "dcbgen.h"
      CHARACTER*24 DAYTID
      LOGICAL TOBE,TOBEK,TRASAM,ACTIVE
      DIMENSION WORK(*),NSTR(2,0:2,4),KVEC(2,5),KQ(2,4),KIBE(2,4),
     &          IR(2,2),NR(2),NQ(2,4),NQT(4),KE(2,4),
     &          NSPC(2,0:2),KQC(2),NSTRT(4)
C
#include "memint.h"
      KFRSAV = KFREE
C
C     skip integral transformation if MDCINT/MRCONEE are present
C     ----------------------------------------------------------
      IF( SKIP_ITRAFO ) GOTO 20
C
C     set default values
C
      TRA_ANTIS = .FALSE.
      TRA_MSOUT = .TRUE.
      NO2IND    = .FALSE.
      NO4IND    = .FALSE.
      NOPAIR    = .TRUE.
      PRPTRA    = .FALSE.
      RCORBS    = .FALSE.
      NOMDCINT  = .FALSE.
      MDCSCAT   = .FALSE.
      THROUT    = 1.0D-14
      SCRTRA    = 1.0D-14
      IPAR4BS   = -1
      ITRA_INTFLG = INTGEN
      ITRA_INTFL4 = ITRA_INTFLG
      ITRA_INTFL2 = ITRA_INTFLG
      ISTRAT      = IAOMOSCHEME
      IPRTRA      = IPRMP2NO
C 
      DO I = 1,NFSYM
        WRITE(TRA_CORSTR(I),'(A)') 'not specified'
        TRA_INDSTR(I)      = VECMP2NO(I)
        DO J = 1,2
          TRA2_INDSTR(J,I) = VECMP2NO(I)
        ENDDO
        DO J = 1,4
          TRA4_INDSTR(J,I) = VECMP2NO(I)
        ENDDO
      ENDDO
C
      CALL GETTIM(CPU1,WALL1)
      CALL TRAHI(IPRTRA,3)
C
C     Define number of double quaternionic classes
C
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Set up index arrays for active orbitals
C     =======================================
C
      IF(NOPAIR) THEN
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
C     NSTR(I): Number of Spinors to be Transformed for index I
C
      DO I = 1,4
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = - 1
          CALL  NUMLST(TRA4_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),
     &                NPSH(IFRP),NESH(IFRP),
     &                NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('MP2DEN.1',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO
C
      KFRSAV = KFREE

      CALL Make_Kramer_to_SpinorIndex (NFSYM,NSTR)
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(/A)') '* Orbital ranges for 4-index transformation:'
      DO I = 1, NFSYM
         CALL TRAPRI(4,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),
     &                 WORK(KVEC(I,3)),WORK(KVEC(I,4)),NSTR)
      ENDDO
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 4 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 4
         NQT(I) = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQT(I) = NQT(I) + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM)
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
      ENDDO
C
C     Allocate the memory for the coefficients, check their relations
C
      K = 2
      IF (NFSYM.EQ.1) K=1
      DO I = 1, 4
C
         DO J = I, 1, -1
            IF (TRASAM(WORK(KVEC(1,I)),WORK(KVEC(K,I)),
     &                 WORK(KVEC(1,J)),WORK(KVEC(K,J)),
     &          NSTR(1,0,I),NSTR(1,0,J))) ISAME(I) = J
         ENDDO
C
         IF (ISAME(I).EQ.I) THEN
            IF (I.NE.1)CALL Abend2('MP2DEN: All ranges should be equal')
            CALL MEMGET('REAL',KQ(1,I),NQT(I),WORK,KFREE,LFREE)
         ELSE
            KQ(1,I) = KQ(1,ISAME(I))
         ENDIF
         KQ(2,I) = KQ(1,I) + NQ(1,I)
C
      ENDDO
C
C     Selected eigenvalues
C
      DO I = 1,4
        NDIM = 0
        DO IFRP = 1, NFSYM
          NDIM = NDIM + NSTR(IFRP,0,I)
        ENDDO
        CALL MEMGET('REAL',KE(1,I),NDIM,WORK,KFREE,LFREE)
        KE(2,I) = KE(1,I) + NSTR(1,0,I)
      ENDDO
C
C     Allocate memory for the integer array with information about the
C     spinors.
C
      DO I = 1, 4
       DO IFRP = 1, NFSYM
        CALL MEMGET('INTE',KIBE(IFRP,I),NSTR(IFRP,0,I),WORK,KFREE,LFREE)
       ENDDO
      ENDDO
C
C     Get all coefficients
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(/A)')
     &   'MP2DEN ERROR: MO coefficient file not present'
         CALL Abend2('MP2DEN: No MO coefficient file !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MP2DEN')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &               WORK(KIBEIG),TOTERG,14,WORK,KFREE,LFREE)
         CLOSE(LUCOEF,STATUS='KEEP')
      END IF
C
C     Select the set that we need
C
      DO I = 1, 4
        DO IFRP = 1, NFSYM
          IF(NSTR(IFRP,0,I).GT.0) THEN
             CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   NSTR(IFRP,0,I),WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
             CALL SELEIG(WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &                 WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
             CALL SELIBEIG(WORK(KIBEIG),IORB(IFRP),IFRP,
     &                     WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                     NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          ENDIF
        ENDDO
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('MP2DEN.2',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
C     If necessary recanonize (virtual) orbitals
C     Assuming that all sets are equal, we use the first index !!!
C     ============================================================
C
      IF (RCORBS) THEN
C        Overdimension KCMO because it is used as Focksize scratch space
         CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
         CALL VIRCAN(WORK,KFREE,LFREE,IPRTRA,KQ(1,1),KE(1,1),
     &               KIBE(1,1),NDMOQR(1,1,1),ICMOQR(1,1),
     &               WORK(KCMO))
         CALL MEMREL('MP2DEN.3',WORK,KFRSAV,KCMO,KFREE,LFREE)
      ENDIF
C
C     Set up index arrays for core orbitals
C     =====================================
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 1
        IR(2,I) = NESH(I)
      ENDDO
C
      DO IFRP = 1,NFSYM
         CALL MEMGET('INTE',KVEC(IFRP,5),NR(IFRP),WORK,KFREE,LFREE)
         NSPC(IFRP,0) = - 1
C
C        If the user did not specify the core index range explictly we
C        take the occupied orbitals that do not belong to the active
C        set.
C
         IF (INDEX(TRA_CORSTR(IFRP),'not specified').NE.0) THEN
            CALL NUMCOR (NOCC(IFRP),NSTR(IFRP,1,1),WORK(KVEC(IFRP,1)),
     &                   NSPC(IFRP,1),WORK(KVEC(IFRP,5)))
            NSPC(IFRP,2) = 0
         ELSE
            CALL NUMLST(TRA_CORSTR(IFRP),WORK(KVEC(IFRP,5)),
     &                  NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                  IFRP,NSPC(IFRP,0))
            CALL ORBCNT(WORK(KVEC(IFRP,5)),NSPC(IFRP,0),
     &                  NPSH(IFRP),NESH(IFRP),
     &                  NSPC(IFRP,2),NSPC(IFRP,1))
         ENDIF
         NSPC(IFRP,0) = NSPC(IFRP,1) + NSPC(IFRP,2)
         CALL MEMREL('MP2DEN.4',WORK,KFRSAV,KVEC(IFRP,5),KFREE,LFREE)
         CALL MEMGET('INTE',KVEC(IFRP,5),NSPC(IFRP,0),
     &               WORK,KFREE,LFREE)
         NDMOQC(1,IFRP,1) = NFBAS(IFRP,0)
         NDMOQC(2,IFRP,1) = NSPC(IFRP,0)
         NDMOQC(1,IFRP,2) = NFBAS(IFRP,0)
         NDMOQC(2,IFRP,2) = NSPC(IFRP,0)
      ENDDO
C####################################################
C
C     number of active occupied
C
      NAOCCT=0
      DO IFRP = 1, NFSYM
         NAOCC(IFRP) = NOCC(IFRP) - NSPC(IFRP,0)
         NAOCCT = NAOCCT + NAOCC(IFRP)
      ENDDO
C        
C     number of inactive occupied
C    
      NIOCCT=0
      DO IFRP = 1, NFSYM
         NIOCC(IFRP) = NOCC(IFRP) - NAOCC(IFRP)
         NIOCCT = NIOCCT + NIOCC(IFRP)
      ENDDO
C
C     number of active virtuals
C     
      NAVIRT=0
      DO IFRP = 1, NFSYM
         NAVIR(IFRP) = NSTR(IFRP,0,1) - NAOCC(IFRP)
         NAVIRT = NAVIRT + NAVIR(IFRP)
      ENDDO
C
C     if selection criterion is set check for max. virt. 
C
      IF( DEF_EN )THEN
        IF( NAVIRT .gt. IMAXVSP ) THEN
          WRITE(LUPRI,*)' NOTE: number of active virt.:',NAVIRT
          WRITE(LUPRI,*)' maximum number of act. virt.:',IMAXVSP
          WRITE(LUPRI,*)' you may reconsider your choice of 
     & active space'
        END IF
      END IF
C     
C     number of inactive virtuals
C 
      NIVIRT=0
      DO IFRP = 1, NFSYM 
         NIVIR(IFRP) = NORB(IFRP) - NOCC(IFRP) - NPSH(IFRP) -NAVIR(IFRP)
         NIVIRT = NIVIRT + NIVIR(IFRP)
      ENDDO
C####################################################
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(/A)')
     &     '* Core orbital ranges for 2-index transformation:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,5)),0,0,0,NSPC)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
C     ICMOQC : Index to start of representation in the coefficient array
C
      NQT(1) = 0
      ICMOQC(1,1) = 1
      DO IFRP = 1, NFSYM
         NQ(IFRP,1) = NFBAS(IFRP,0)*NSPC(IFRP,0)*NZ
         NQT(1) = NQT(1) + NQ(IFRP,1)
         IF (IFRP.LT.NFSYM) ICMOQC(IFRP+1,1) = ICMOQC(IFRP,1)+NQ(IFRP,1)
      ENDDO
C
C     Allocate the memory for the coefficients
C
      CALL MEMGET('REAL',KQC(1),NQT(1),WORK,KFREE,LFREE)
      KQC(2) = KQC(1) + NQ(1,1)
C
C     Get all coefficients
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(/A)') 'MO coefficient file not present'
         CALL QUIT('PAMTRA: No MO coefficient file !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','MP2DEN')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),
     &               WORK(KIBEIG),TOTERG,14,WORK,KFREE,LFREE)
      ENDIF
C
C     Select the set that we need
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(KQC(IFRP)),
     &                NSPC(IFRP,0),WORK(KVEC(IFRP,5)),
     &                NSPC(IFRP,2),NSPC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      CALL MEMREL('MP2DEN.5',WORK,KFRSAV,KCMO,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
C
      IF (.NOT.NO2IND)
     &   CALL TRAONE(WORK,KFREE,LFREE,NSTR,NSTRT,NSPC,KQ,KQC,KE,KIBE)
C
      CALL FLSHFO(LUPRI)
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      IF (.NOT.NO4IND)
     &   CALL PAMTR1(WORK,KFREE,LFREE,IPRTRA,KQ,KE,KIBE,
     &               NDMOQR,ICMOQR,NSTR,.FALSE.,DUMMY,
     &               TRA_MSOUT,TRA_ANTIS,ITRA_INTFL4)
C
C      Print timing information
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
C
      IMINS  = INT(WALL)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(WALL) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(//A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total wall time used for transformation :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CPU    = CPU2 - CPU1
      IMINS  = INT(CPU)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(CPU) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total CPU  time used for transformation :',
     &      IHOURS,':',IMINS,':',ISECS
C
      WRITE(LUPRI,'(/A,F20.6)')
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)') ' Transformation ended at : ',DAYTID
C
      CALL FLSHFO(LUPRI)
C
 20   CONTINUE
C
C     Transformation finished - start MP2 density run
C     -----------------------------------------------
C
C     prepare and write RELCCSD MP2-NO input
C
      LMP2_DEN_INP = 99
      MP2NOFN      = 'MP2_NO.INP.0'
      MP2NOFN_L    = 12
      INFO_STRING(1) =
     &'&RELCCSD TIMING=T, DOSORT=T, DOENER=T, DOFOPR=T, DOSOPR=F &END'
      INFO_STRING(2) =
     &'&CCENER DOMP2=T, DOCCSD=F, DOCCSDT=F  &END                    '
      INFO_STRING(3) =
     &'&CCFOPR DOMP2G=T, NOZG=T, DONATORB=T  &END                    '

      OPEN(LMP2_DEN_INP,FILE=MP2NOFN(1:MP2NOFN_L),STATUS='NEW',
     &     FORM='FORMATTED',ACCESS='SEQUENTIAL')
      REWIND(LMP2_DEN_INP)
      WRITE(LMP2_DEN_INP,'(1X,A62)') INFO_STRING(1) 
      WRITE(LMP2_DEN_INP,'(1X,A62)') INFO_STRING(2) 
      WRITE(LMP2_DEN_INP,'(1X,A62)') INFO_STRING(3) 
C
      CALL PAMCCM()
C
      CLOSE(LMP2_DEN_INP,STATUS='KEEP')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck READMP2NAT */
      SUBROUTINE READMP2NAT(CMO,EIG,EIG_NAT)
C***********************************************************************
C
C     Read MP2-NO coefficients and occupations numbers from file 
C     KRMCOLD (has to be opened in calling routine). 
C     Expand to full EIG
C
C     written by S.Knecht - July 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmp2no.h"
      DIMENSION CMO(*), EIG(*), EIG_NAT(*)
C
C     initialize
      CALL DZERO(EIG,NORBT)
C     read MP2-NO vecs and occupancies from file KRMCOLD
      CALL REAKRMC(LUKRMC,'NEWNATOB',CMO,NCMOTQ)
      CALL REAKRMC(LUKRMC,'MP2NOOCC',EIG_NAT,NESHT)
C
      IOFF2 = 1
C
      DO IFRP = 1, NFSYM
C        ... occupancies
         CALL DCOPY(NESH(IFRP),EIG_NAT(IOFF2),1,
     &              EIG(IORB(IFRP)+1+NPSH(IFRP)),1)
C
         IOFF2 = IOFF2 + NESH(IFRP)
      END DO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TDMP2NATMO */
      SUBROUTINE TDMP2NATMO(CMO,CMO_NAT,WORK,LWORK)
C***********************************************************************
C
C     transform to MP2 natural orbitals using transformation matrix
C     CMO_NAT
C
C     written by S.Knecht - July 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbmp2no.h"
      DIMENSION CMO(*), CMO_NAT(NORBT,NESHT,NZ)
      DIMENSION WORK(*)
#include "memint.h"
C
C     Transform orbitals
C     ------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KCMO),NCMOTQ)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)

      IPRMP2NO_SAVE = IPRMP2NO
CSK   IPRMP2NO = 5
      IF( IPRMP2NO .ge. 5 ) THEN 
         CALL HEADER('initial orbitals:',-1)
         DO I = 1,NFSYM
            CALL HEADER('Coefficients for irrep '//CHAR(I+48),-1)
            WRITE(LUPRI,*) ' printing offset ICMOQ(I)+1',ICMOQ(I)+1
            WRITE(LUPRI,*) ' NFBAS(I,0), NORB(I)',NFBAS(I,0),NORB(I)
            CALL PRQMAT(WORK(KCMO+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDDO
      END IF
C
      IOFF2 = 1
C
      DO IFRP = 1, NFSYM
C
C        CMO(new) = CMO(old) * CMO_NAT
C
         NBASI    = NFBAS(IFRP,0)
         NORBI    = NORB(IFRP)
         NESHI    = NESH(IFRP)
C
         ICMO_ESH = ICMOQ(IFRP) + NBASI*NPSH(IFRP)
         IOFF1    = IORB(IFRP)  + NPSH(IFRP) + 1
C         
csk      WRITE(LUPRI,*) ' IOFF2, IOFF1, ICMO_ESH', 
csk  &                    IOFF2, IOFF1, ICMO_ESH
csk      WRITE(LUPRI,*) ' IORB(IFRP), NPSH(IFRP), NORBI', 
csk  &                    IORB(IFRP), NPSH(IFRP), NORBI
C        WRITE(LUPRI,*) ' diagonalizing offset:1+ICMO_ESH',1+ICMO_ESH
C
!                    m      n     k
         CALL QGEMM(NBASI,NESHI,NESHI,D1,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KCMO+ICMO_ESH),NBASI,NORBI,NZ, ! m x k
     &        'N','N',IPQTOQ(1,0),
     &        CMO_NAT(IOFF1,IOFF2,1),NORBT,NESHT,NZ, ! k x n
     &        D0,IPQTOQ(1,0),
     &        CMO(1+ICMO_ESH),NBASI,NORBI,NZ)        ! m x n
C
         IOFF2 = IOFF2 + NESH(IFRP)
C
      END DO
C
      IF( IPRMP2NO .ge. 5 ) THEN 
         CALL HEADER('final MP2 natural orbitals:',-1)
         DO I = 1,NFSYM
            CALL HEADER('Coefficients for irrep '//CHAR(I+48),-1)
            WRITE(LUPRI,*) ' printing offset ICMOQ(I)+1',ICMOQ(I)+1
            WRITE(LUPRI,*) ' NFBAS(I,0), NORB(I)',NFBAS(I,0),NORB(I)
            CALL PRQMAT(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &                  NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
      IPRMP2NO = IPRMP2NO_SAVE
C
C     Memory deallocation
C     -------------------
      CALL MEMREL('TDMP2NO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      END
