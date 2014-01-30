!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
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
C  /* Deck prpexp */
      SUBROUTINE PRPEXP(WF,EXPVAL,WORK,LWORK)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*),WORK(LWORK)
      CHARACTER WF*4
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NORBT   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST ,2*NEXPP   ,WORK,KFREE,LFREE)

CHJ TODO : check if 2 DMAT or NOPEN+1 DMAT is the correct ...
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX*2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT  ,N2BBASX*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)


CMI ... Needed for SQUARE matrixes of some property operators!
C     CALL MEMGET('REAL',KONEIN ,   N2BBASX,WORK,KFREE,LFREE)
C

#ifdef ANALYZE_PROPERTY_GRADIENT
      LUNIT = 14
      CALL OPNFIL(LUNIT,'PMAT','UNKNOWN','PAMEXP')
#endif
      CALL PRPEX1(WF,EXPVAL,WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),
     &            WORK(KFIRST),WORK(KFREE),LFREE)
#ifdef ANALYZE_PROPERTY_GRADIENT
      CLOSE(LUNIT,STATUS = 'KEEP')
#endif
C
C     Memory deallocation
C
      CALL MEMREL('PRPEXP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRPEXP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1 */
      SUBROUTINE PRPEX1(WF,EXPVAL,LEXPST,CMO,EIG,DMAT,PMAT,
     &                   FIRST,WORK,LWORK)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C     Last modifications: M.Ilias, March2006, Strasbourg
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DTOL=1.0D-16)
C
      LOGICAL FIRST(NZ),LEXPST(2,*),FILEX,IS_SQUARE
      DIMENSION EXPVAL(*),CMO(*),DMAT(*),PMAT(*),
     &          EIG(*),WORK(LWORK),NP(2,0:2)
      CHARACTER MXFORM*6,FMT*6
      CHARACTER WF*4
      INTEGER   GETIOPEN
      LOGICAL FNDLAB
      
C
#include "cbihr1.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbdhf.h"
#include "dcbprj.h"
C
C
C     Initialize
C
#include "memint.h"
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Initialize
C     ==========
C
      CALL DZERO(EXPVAL,NEXPP)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Generate density matrix
C     =======================
C
      IF (WF .EQ. 'DHF' .OR. WF .EQ. 'DFT') THEN
        INQUIRE(FILE='DFCOEF',EXIST=FILEX)
        IF(.NOT.FILEX) GOTO 1000
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PRPEXP')
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6,
     &               WORK,KFREE,LFREE)
        CLOSE(LUCOEF,STATUS='KEEP')
        IF (.NOT. ORBANA) THEN
C.........accumulate all density into one matrix
          FILEX = AOC
          AOC = .FALSE.
          CALL DENMAT(DMAT,CMO,IPREXP)
          AOC = FILEX
        END IF
      ELSE IF (WF .EQ. 'MP2 ' .OR. WF .EQ. 'CC  ') THEN
C
         IF (WF.EQ.'CC  ') 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'CC properties not implemented')
C
         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for MP2')
C
C        Get relaxed density matrix computed by RELCCSD.
C
         CALL GET_CC_DENSITY (WF,'AO',DMAT,IERR)
C
      ELSE IF (WF .EQ. 'KRMC') THEN
C
         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for MCSCF')
C
         INQUIRE(FILE='KRMCSCF', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "KRMCSCF" not found!')
C
         call dzero(dmat,N2BBASXQ*2)

         CALL OPNFIL(LUCOEF,'KRMCSCF','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'NEWORB  ',CMO,NCMOTQ)
         CLOSE(LUCOEF,STATUS='KEEP')

         INQUIRE(FILE='MCDENS', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "MCDENS" not found!')
         CALL OPNFIL(LUCOEF,'MCDENS','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'MC  DENS',DMAT(1),nasht**2*nz)
!        total DV = DV_act + DV_core
         IF (NASHT .GT. 0) THEN ! DV_act
           CALL DVSCT(DMAT,DMAT(1+N2BBASXQ))
         END IF
         call dzero(pmat,N2BBASXQ)
         CALL DCSCT(pmat) ! DV_core
         call daxpy(N2BBASXQ,1.0d0,pmat,1,DMAT(1+N2BBASXQ),1)
         call dzero(pmat,N2BBASXQ)
         CLOSE(LUCOEF,STATUS='KEEP')
C
C        Get AO density matrix.
!        input: C == U+ A U
C
         DO I = 1, NFSYM
            CALL QTRANS('MOAO','S',D0,
     &           NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &           DMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(1+N2BBASXQ+I2ORBX(I,I)),
     &           NORBT,NORBT,NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK,LWORK,0)
         END DO
C
C        Result is scaled by 2 in PRPEX2...
C
         CALL DSCAL(N2BBASXQ,0.5d0,DMAT,1)

!        WRITE(LUPRI,'(/2X,A,I1)') 'PRPEX1: DMAT of IZ=',1
!        CALL OUTPUT(DMAT,1,NTBAS(0),1,NTBAS(0),
!    &               NTBAS(0),NTBAS(0),1,LUPRI)
C
      ELSE
         WRITE(LUPRI,'(1X,2A)')
     &        '*** ERROR in PRPEX1 *** Unknown wave function: ',WF
         CALL QUIT('*** ERROR in PRPEX1 *** Unknown wave function')
      END IF
C
C     Projection analysis
C     ===================
C
      IF(PRJANA) THEN
        INQUIRE(FILE='DFPRJC',EXIST=FILEX)
        IF(.NOT.FILEX) THEN
          WRITE(LUPRI,'(A)')
     &      '** WARNING ! No projection coefficients found !',
     &      '    - Skipping expectation value projection analysis !'
          PRJANA = .FALSE.
          KRVEC = KFREE
          KBVEC = KFREE
          KPRMT = KFREE
        ELSE
          CALL OPNFIL(LUCOEF,'DFPRJC','OLD','PRPEXP')
          READ(LUCOEF) NREFS
          REWIND LUCOEF
          CALL MEMGET('INTE',KNR,NFSYM*NREFS,WORK,KFREE,LFREE)
          CALL REAPRJ(LUCOEF,NP,WORK(KNR),KRVEC,KBVEC,
     &                WORK,KFREE,LFREE)
          CLOSE(LUCOEF,STATUS='KEEP')
          NBMAX = 0
          NPFMT = 0
          DO IFRP = 1,NFSYM
            NBMAX = MAX(NBMAX,NP(IFRP,0))
            NPFMT = NPFMT + NP(IFRP,0)*NP(IFRP,0)*NZ
          ENDDO
          NN = NBMAX*(NBMAX+1)/2
          CALL MEMGET('REAL',KPFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('REAL',KDFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('INTE',KIBUF,NN,WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KINDX,NBMAX*2,WORK,KFREE,LFREE)
          NTOT = NREFS + 1
          NFDIM = NTOT*NTOT
          CALL MEMGET('REAL',KFMAT,NFDIM,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KPM  ,NTOT,WORK,KFREE,LFREE)
C.........Calculate density matrix in fragment basis
          CALL PRJEX4(WORK(KDFMT),WORK(KBVEC),NP(1,0),NP(1,2),IPREXP)
        ENDIF
      ELSE
        KRVEC = KFREE
        KBVEC = KFREE
        KPFMT = KFREE
      ENDIF
C
C
C     Expectation values normal procedure
C     ===================================
C
      IF(.NOT.ORBANA) THEN

        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)

        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          IS_SQUARE = (IPRPTIM(INDXPR).EQ.0)
          IF (IS_SQUARE) THEN
            IDIMENS = N2BBASX
            CALL MEMGET('REAL',KPTRI  ,NZ*IDIMENS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KOP1INT,   IDIMENS,WORK,KFREE,LFREE)
          ELSE
            IDIMENS = NNBBASX
            CALL MEMGET('REAL',KPTRI  ,NZ*IDIMENS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KOP1INT,   IDIMENS,WORK,KFREE,LFREE)
          ENDIF
C.........generate property matrix in SA-AO-basis into PMAT ...
CMI ... extended parameter list by the OP1INT array for square matrixes
          CALL PRPEX2(INDXPR,EXPVAL(I),.true.,LEXPST(1,I),IDIMENS,
     &                WORK(KPTRI),PMAT,WORK(KOP1INT),FIRST,
     &                WORK,KFREE,LFREE,IPREXP)

          IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 30

C.............calculate expectation value
          CALL PRPEX3(EXPVAL(I),INDXPR,DMAT,PMAT,FIRST,
     &        WORK,KFREE,LFREE,IPREXP)
 30       CONTINUE

          IF(PRJANA) THEN
C.........Transform property matrix to fragment orbital basis
CMI        CALL PRJEX1(WORK(KPFMT),IDIMENS,WORK(KPTRI),PMAT,WORK(KRVEC),
Caspg 2006-07-07
C quick hack to get the projection analysis to work, before IREPB and ITIMB
C were undefined, hence the crash.
C  irepb: Molecular point group symmetry
C  itimb: Time reversal symmetry
            IREPB = IPRPSYM(INDXPR)-1
            ITIMB = IPRPTIM(INDXPR)
            CALL PRJEX1(WORK(KPFMT),PMAT,WORK(KRVEC),
     &                  INDXPR,IREPB,ITIMB,NP(1,0),FIRST,IPREXP,
     &                  WORK,KFREE,LFREE)
Caspg
          ENDIF
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &       PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &       LEXPST(1,I),LEXPST(2,I)
          CALL FLSHFO(LUPRI)
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),IDIMENS,WORK(KPTRI),WORK(KPFMT),
     &                 WORK(KDFMT),WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),IPREXP,PMAT,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',76,0)
          ENDIF
          IF(PRPCAN) THEN
            CALL EXPCAN(INDXPR,PMAT,CMO,IPREXP,WORK,KFREE,LFREE)
          ENDIF
CMI     ... release the allocated memory ...
          CALL MEMREL('PRPEX1.nooana',WORK,KWORK,KPTRI,KFREE,LFREE)
        ENDDO

        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C     Expectation values; individual orbital contributions
C     ====================================================
C
      ELSE
        CALL TITLER('Expectation values','*',127)
        WRITE(LUPRI,'(51X,2(A2,1X))') 's0','t0'
        CALL PRSYMB(LUPRI,'-',76,0)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          IS_SQUARE = (IPRPTIM(INDXPR).EQ.0)
          IF (IS_SQUARE) THEN
            IDIMENS = N2BBASX
            CALL MEMGET('REAL',KPTRI  ,NZ*IDIMENS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KOP1INT,   IDIMENS,WORK,KFREE,LFREE)
          ELSE
            IDIMENS = NNBBASX
            CALL MEMGET('REAL',KPTRI  ,NZ*IDIMENS,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KOP1INT,   IDIMENS,WORK,KFREE,LFREE)
          ENDIF

          CALL PRSYMB(LUPRI,'-',60,4)
          CALL HEADER('Operator '//PRPNAM(INDXPR)//':',0)            
          CALL HEADER('Expectation value for individual orbitals',-1)
          EXPVAL(I) = D0
          DO 10 IFRP = 1,NFSYM
            IF(NORB(IFRP).EQ.0) GOTO 10
C...........generate property matrix in SO-basis
            CALL PRPEX2(INDXPR,EXPBUF,.true.,LEXPST(1,I),IDIMENS,
     &                  WORK(KPTRI),PMAT,WORK(KOP1INT),FIRST,
     &                  WORK,KFREE,LFREE,IPREXP)

            IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 20

            DO J = 1,NOCC(IFRP)
              CALL DENORB(DMAT,J,1,IFRP,CMO,IPREXP)
              EXPBUF = D0
C.............calculate expectation value
              CALL PRPEX3(EXPBUF,INDXPR,DMAT,PMAT,FIRST,
     &          WORK,KFREE,LFREE,IPREXP)
              IOPEN = GETIOPEN(J,IFRP)
              EXPBUF = EXPBUF*DF(IOPEN)
              EXPVAL(I) = EXPVAL(I) + EXPBUF
              FMT = MXFORM(EXPBUF,15)
              WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',1X,A4)')
     &          FREP(IFRP),J,EXPBUF,'a.u.'
              CALL FLSHFO(LUPRI)
            ENDDO

            IF(PRJANA) THEN
C...........Transform property matrix to fragment orbital basis
Caspg
C quick hack to get the projection analysis to work, before IREPB and ITIMB
C were undefined, hence the crash.
C  irepb: Molecular point group symmetry
C  itimb: Time reversal symmetry
      IREPB = IPRPSYM(INDXPR)-1
      ITIMB = IPRPTIM(INDXPR)
C
           CALL PRJEX1(WORK(KPFMT),PMAT,WORK(KRVEC),
     &                  INDXPR,IREPB,ITIMB,NP(1,0),FIRST,IPREXP,
     &                  WORK,KFREE,LFREE)
Caspg
            ENDIF
 10       CONTINUE
 20       CONTINUE
          CALL PRSYMB(LUPRI,'-',60,4)
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A8,3X,'//FMT//',1X,A4,3X,2(3X,L1))')
     &     'Total : ',EXPVAL(I),'a.u.',
     &     LEXPST(1,I),LEXPST(2,I)
          CALL PRSYMB(LUPRI,'-',60,4)
          IF(PRPCAN) THEN
            CALL EXPCAN(INDXPR,PMAT,CMO,IPREXP,WORK,KFREE,LFREE)
          ENDIF
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),IDIMENS,WORK(KPTRI),
     &                  WORK(KPFMT),WORK(KDFMT),WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),IPREXP,PMAT,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',76,0)
          ENDIF
CMI        ... release the memory
          CALL MEMREL('PRPEX1.orbana',WORK,KWORK,KPTRI,KFREE,LFREE)
        ENDDO
C        CALL PRSYMB(LUPRI,'-',76,0)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
      ENDIF
C
      CLOSE(LU1INT,STATUS = 'KEEP')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PRPEXP: Coefficient file not found !'
      CALL QUIT('PRPEXP: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex2 */
      SUBROUTINE PRPEX2(INDXPR,EXPVAL,ONLY_TS,LEXPST,IDIMENS,
     &                  PTRI,PMAT,OP1INT,FIRST,WORK,KFREE,LFREE,
     &                  IPRINT)
C***********************************************************************
C
C     if ONLY_TS = .true.
C        for totally symmetric operators return property matrix
C        for other symmetries return .FALSE. in
C            LEXPST(1) - molecular point group symmetry
C        or  LEXPST(2) - time reversal symmetry
C
C     if ONLY_TS = .false.
C        return property matrix anyway (i need it in OpenRSP --radovan)
C
C     On output: PMAT - the complete propery matrix in SA-AO basis !
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C                   Strasbourg/March 2006 - Miro ILIAS
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcbham.h"
C
!     radovan:                    ONLY_TS: allow only totally symmetric
      LOGICAL FIRST(NZ),LEXPST(2),ONLY_TS
C     DIMENSION PTRI(NNBBASX,NZ),PMAT(N2BBASX,NZ),
      DIMENSION PTRI(IDIMENS,NZ),PMAT(N2BBASX,NZ),
     &          OP1INT(IDIMENS),WORK(*)

      CALL QENTER('PRPEX2')

      EXPVAL = D0
C
C     Check symmetry
C     ==============
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LEXPST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.EQ.-1) LEXPST(2) = .TRUE.

CMI ... control output...
      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/2X,A,/1X,A,I2,1X,A,A,2I2)')
     & '*** Output from PRPEX2  ***: ',
     & ' Operator index and name:',
     &  INDXPR,PRPNAM(INDXPR),
     & ' IREP,ITIM:',IREP,ITIM
       IF (LEXPST(1)) WRITE(LUPRI,'(1X,A)')
     & ' Operator not symmetric! Do not calculate exp.value!'
      ENDIF
C
C     If not totally symmetric, return with zero value
C
      IF (ONLY_TS .AND. LEXPST(1).AND.(.NOT.PRJANA)) GOTO 90 
  109 CONTINUE
!.q 
C
C     Totally symmetric operator ==> proceed
C
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Operator to be read from file
C     =============================
C
CMI ... check it ! In the pure two-component mode there should be no reading from AOPROPER !!
      IOPTYP = IPRPTYP(INDXPR)
      IF(IOPTYP.EQ.0.AND..NOT.(BSS.AND.TWOCOMP)) THEN
        KLBL  = IPRPLBL(1,INDXPR)
        INQUIRE(FILE=PRPLBL(KLBL),EXIST=FIRST(1))
        IF(FIRST(1)) THEN
          CALL OPNFIL(2,PRPLBL(KLBL),'OLD','PRPEX2')
          CALL REAFCK(2,PMAT,.FALSE.,1)
          CLOSE(2,STATUS = 'KEEP')
          PFAC = FACPRP(1,INDXPR)
          IF(PFAC.NE.D1) CALL DSCAL(N2BBASXQ,PFAC,PMAT,1)
          CALL LSET(NZ,.FALSE.,FIRST)
        ELSE
          WRITE(LUPRI,'(A,A8,A,A16,A)')
     &      'Property file ',PRPLBL(KLBL),' for operator ',
     &      PRPNAM(INDXPR),' not found !'
          CALL QUIT('PRPEX2: Property file not found !')
        ENDIF
      ELSEIF (ITIM.EQ.0) THEN
        IF (IDIMENS.NE.N2BBASX) CALL QUIT('PRPEX2: wrong IDIMENS !')
        CALL PRPMSAO(INDXPR,IDIMENS,DUMMY,PMAT,OP1INT,
     &              FIRST,WORK(KFREE),LFREE,IPRINT)
      ELSE
        IF (IDIMENS.NE.NNBBASX) CALL QUIT('PRPEX2: wrong IDIMENS !')
        CALL PRPMSAO(INDXPR,IDIMENS,PTRI,PMAT,OP1INT,
     &              FIRST,WORK(KFREE),LFREE,IPRINT)
      ENDIF
C     
C
 90   CONTINUE     
      CALL QEXIT('PRPEX2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex3 */
      SUBROUTINE PRPEX3(EXPVAL,INDXPR,DMAT,PMAT,FIRST,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C     MI entry for LAO-project, May 2003
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbprp.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      LOGICAL FIRST(NZ),LEXPST(2)
      DIMENSION DMAT(N2BBASX,NZ),WORK(*),
     &          PMAT(N2BBASX,NZ)

      CALL QENTER('PRPEX3')

      IREP = IPRPSYM(INDXPR)-1
      ITIM = IPRPTIM(INDXPR)

      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,'(/A,A16)')
     & 'PRPEX3: Entering total SA-AO matrix of property ',PRPNAM(INDXPR)
        write(lupri,*) '   proerty operator IREP=',IREP,' ITIM=',ITIM
        CALL PRQMAT(PMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &                 IPQTOQ(1,IREP),LUPRI)
      ENDIF        

C========================================================================
C     Get expectation value, using ONLY the matrices containing
C     integrals (ensured by logical FIRST(1..NZ) )
C========================================================================
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        EXPVAL = EXPVAL + DDOT(N2BBASX,DMAT(1,IZ),1,PMAT(1,IZ),1)
        IF (IPRINT.GE.5) THEN
          WRITE(lupri,*)
          WRITE(LUPRI,*)
     &       ' == PRPEX3: accumulated EXPVAL (ORBANA=.false.)'//
     &       ' after IZ=',IZ,
     &       '/',NZ,' exp.val=',2*EXPVAL
          WRITE(LUPRI,'(/2X,A,I1)') 'PRPEX3: PMAT of IZ=',IZ
          CALL OUTPUT(PMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),1,LUPRI)
          WRITE(LUPRI,'(/2X,A,I1)') 'PRPEX3: DMAT of IZ=',IZ
          CALL OUTPUT(DMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),1,LUPRI)
        ENDIF
      ENDIF
      ENDDO
      EXPVAL = EXPVAL + EXPVAL
C
      CALL QEXIT('PRPEX3')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex1 */
      SUBROUTINE PRJEX1(PFO,PAO,RVEC,INDXPR,IREP,ITIM,NCNT,
     &                  FIRST,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform property matrix to fragment orbital basis
C
CMI   On input:  PAO - property matrix in SA-AO basis
C
C     Written by Trond Saue Apr 27 2003 
C     Last modifications: M.Ilias, March 2006, Strasbourg
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION PFO(*),PAO(N2BBASX,NZ),NCNT(*),RVEC(*),WORK(*)

      CALL QENTER('PRJEX1')
C
C     Transform on the fly
C
      IOPTYP = IPRPTYP(INDXPR)
      NDIM = 0
      DO IFRP = 1,NFSYM
        NDIM = NDIM + NCNT(IFRP)*NCNT(IFRP)*NZ
      ENDDO
      CALL DZERO(PFO,NDIM)
      FADD = D0

      DO IZ = 1,NZ

      IF(.NOT.FIRST(IZ)) THEN
         II = IZ
CC         IF(IOPTYP.EQ.0) THEN
C this is a fock matrix...
CC           II = IZ
CC         ELSE
CC           II = 1
CC CMI      ... this is deactivated as property matrix PAO is ready
CC CMI       IQ = IPQTOQ(IZ,IREP)
CC C         IH = IHQMAT(IQ,ITIM)
CC C         IF(IH.EQ.1) THEN
CC C           CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PAO)
CC C         ELSEIF(IH.EQ.2) THEN
CC C           CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PAO)
CC C         ENDIF
CC C         CALL BUTOBS(PAO,1,WORK(KFREE),LFREE)
CC         ENDIF

        IPA = 1
        IPF = 1
        IPR = 1
Caspg
        DO IFRP = 1,NFSYM
Caspg
C       write (6,*) 'debug prjex1 qtrans',ifrp,irep,iz
C       write (6,*) 'debug prjex1 qtrans',nfbas(ifrp,0),ncnt(ifrp)
C       write (6,*) 'debug prjex1 qtrans',(ipqtoq(i,irep),i=1,nz)
Caspg
          CALL QTRANS('AOMO','S',FADD,
     &         NFBAS(IFRP,0),NFBAS(IFRP,0),NCNT(IFRP),NCNT(IFRP),
     &         PAO(IPA,II),NTBAS(0),NTBAS(0),1,IPQTOQ(IZ,IREP),
     &         PFO(IPF),NCNT(IFRP),NCNT(IFRP),NZ,IPQTOQ(1,IREP),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         WORK(KFREE),LFREE,IPRINT)
          IPA = IPA + NFBAS(IFRP,0)*(NTBAS(0)+1)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
          IPR = IPR + NFBAS(IFRP,0)*NCNT(IFRP)*NZ
        ENDDO
        FADD = D1
      ENDIF
      ENDDO
C
#ifdef ANALYZE_PROPERTY_GRADIENT
      LUNIT = 14
      WRITE(LUNIT) PRPNAM(INDXPR)
      IPF = 1
      DO IFRP = 1,NFSYM
        NDIM = NCNT(IFRP)*NCNT(IFRP)*NZ
        CALL WRITT(LUNIT,NDIM,PFO(IPF))
        IPF = IPF + NDIM
      ENDDO
#endif
C
      IF(IPRINT.GE.2) THEN
        IPF = 1
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A16,A,I3)')
     &    'Property matrix of ',PRPNAM(INDXPR),
     &    'in fragment basis..',IFRP
          CALL PRQMAT(PFO(IPF),NCNT(IFRP),NCNT(IFRP),
     &         NCNT(IFRP),NCNT(IFRP),NZ,
     &         IPQTOQ(1,IREP),LUPRI)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
CTROND
          LUBUF = 35
          CALL OPNFIL(LUBUF,'PRPMAT','UNKNOWN','PAMEXP')
          NOST = NCNT(IFRP)*NCNT(IFRP)
          CALL WRITT(LUBUF,NOST,PFO)
          CLOSE(LUBUF,STATUS='KEEP')
        ENDDO       
      ENDIF
C
      CALL QEXIT('PRJEX1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex2 */
      SUBROUTINE PRJEX2(EXPVAL,IDIMENS,PTRI,PFO,DFO,BVEC,IBUF,
     &                  INDXPR,NCNT,NMO,
     &                  NR,NTOT,FMAT,PM,INDX,
     &                  IPRINT,PBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Make generalized population matrix
C
C     Written by Trond Saue Apr 27 2003 
C
CMI    MI/March 2006 ... seems that PTRI is needed...
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprj.h"
C
      DIMENSION PFO(*),DFO(*),NCNT(*),NMO(*),
C    &          BVEC(*),PTRI(NNBBASX),IBUF(*),PBUF(*),
     &          BVEC(*),PTRI(IDIMENS,NZ),IBUF(*),PBUF(*),
     &          NR(NREFS,NFSYM),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*),
     &          WORK(*)
      CHARACTER MXFORM*6,FMT*6
C
      N2TOT = NTOT*NTOT
      CALL DZERO(FMAT,N2TOT)
      CALL DZERO(PM,NTOT)
C
C     Make population matrix in fragment orbitals
C
      IPD = 1
      DO 10 IFRP = 1,NFSYM
        IF(NMO(IFRP).EQ.0) GOTO 10
        II = 0
        DO I = 1,NREFS
          DO J = 1,NR(I,IFRP)
            II = II + 1
            INDX(1,II) = I
            INDX(2,II) = J
          ENDDO
        ENDDO
        I = NREFS + 1
        DO J = 1,NMO(IFRP)
          II = II + 1
          INDX(1,II) = I
          INDX(2,II) = J
        ENDDO
        CALL PRJEX3(EXPVL2,IDIMENS,PTRI,PFO(IPD),DFO(IPD),NCNT(IFRP),
     &              NR(1,IFRP),NMO(IFRP),IBUF,PBUF,FMAT,PM,INDX,
     &              NTOT,INDXPR,IFRP,IPRINT)
        IPD = IPD + NCNT(IFRP)*NCNT(IFRP)*NZ
 10   CONTINUE
      CALL DSCAL(NTOT,D2,PM,1)
C
C     Give fragment contributions
C
      CALL HEADER('* Total fragment contributions',-1)
      WRITE(LUPRI,'(A)') ' * Intraatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        TMP = D2*FMAT(J,J)
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        FMT = MXFORM(TMP,25)
        WRITE(LUPRI,'(I5,2X,A6,16X,'//FMT//',E16.5)') 
     &    J,REFFIL(J),TMP,FAC
        FMT = MXFORM(PM(J),25)
        WRITE(LUPRI,'(9X,A,'//FMT//',E16.5)')
     &    '- principal moment: ',PM(J),PM(J)/EXPVAL
        TMP = TMP - PM(J)
        FMT = MXFORM(TMP,25)
        WRITE(LUPRI,'(9X,A,'//FMT//',E16.5)')
     &    '- hybridization   : ',TMP,TMP/EXPVAL
      ENDDO
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &      'Total: ',TOT,TOT/EXPVAL 
      WRITE(LUPRI,'(A)') ' * Interatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        DO I = 1,(J-1)
          TMP = D2*(FMAT(I,J)+FMAT(J,I))
          TOT = TOT + TMP
          FAC = TMP/EXPVAL
          FMT=MXFORM(TMP,25)
          WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
        ENDDO
      ENDDO
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &      'Total: ',TOT,TOT/EXPVAL
      WRITE(LUPRI,'(A)') ' * Polarization contributions:'
      J = NREFS+1
      TOT = D0
      DO I = 1,NREFS
        TMP = D2*(FMAT(I,J)+FMAT(J,I))
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        FMT=MXFORM(TMP,25)
        WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
      ENDDO
      TMP = D2*FMAT(J,J)
      TOT = TOT + TMP
      FAC = TMP/EXPVAL
      FMT=MXFORM(TMP,25)
      WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &    J,REFFIL(J),J,REFFIL(J),TMP,FAC
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &     'Total: ',TOT,TOT/EXPVAL
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_efg */
      SUBROUTINE NUC_EFG(IPRINT,DOATOM)
C***********************************************************************
C     
C     Calculates nuclear contribution to electric field gradients
C
C     Written by T.Enevoldsen  - Sep 1996
C     Last revision by W.A. de Jong - April 1999
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D0=0.0D0,D3 = 3.0D0)
      DIMENSION  GEOM(3*8*MXCENT)
      LOGICAL    DOATOM(*)
#include "nuclei.h"
#include "symmet.h"
#include "nqcc.h"
#include "ibtfun.h"
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
C
      CALL DZERO(EFGNUC,9*MXCENT)
      KK = 1
      NATOM = 0
      IATOMD = 0
      DO 10 I1 = 1, NUCIND
      IF(DOATOM(I1)) THEN
        DO 20 IDEG1 = 1,NUCDEG(I1)
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          LL = 1
          DO 100 I2 = 1, NUCIND
          IF (CHARGE(I2).EQ.D0) GOTO 100
            DO 110 IDEG2=1,NUCDEG(I2)
              IF (LL .NE. KK) THEN
                XCOOR = GEOM(LL) - GEOM(KK)
                YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                R5 = R2*R2*SQRT(R2)
                EFGNUC(1,1,NATOM) = EFGNUC(1,1,NATOM)
     &                            + CHARGE(I2)*(D3*XCOOR*XCOOR - R2)/R5
                EFGNUC(2,2,NATOM) = EFGNUC(2,2,NATOM)
     &                            + CHARGE(I2)*(D3*YCOOR*YCOOR - R2)/R5
                EFGNUC(3,3,NATOM) = EFGNUC(3,3,NATOM)
     &                            + CHARGE(I2)*(D3*ZCOOR*ZCOOR - R2)/R5
                EFGNUC(1,2,NATOM) = EFGNUC(1,2,NATOM)
     &                            + CHARGE(I2)*D3*XCOOR*YCOOR/R5
                EFGNUC(1,3,NATOM) = EFGNUC(1,3,NATOM)
     &                            + CHARGE(I2)*D3*XCOOR*ZCOOR/R5
                EFGNUC(2,3,NATOM) = EFGNUC(2,3,NATOM)
     &                            + CHARGE(I2)*D3*YCOOR*ZCOOR/R5
              END IF
              LL = LL + 3
 110        CONTINUE
 100      CONTINUE
          EFGNUC(2,1,NATOM) = EFGNUC(1,2,NATOM)
          EFGNUC(3,1,NATOM) = EFGNUC(1,3,NATOM)
          EFGNUC(3,2,NATOM) = EFGNUC(2,3,NATOM)
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,2A,/)')
     &             'Nuclear EFG Tensor of nucleus: ',NAMDEP(IATOMD)
            CALL OUTPUT(EFGNUC(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
          END IF
          KK = KK + 3
 20     CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(I1)
        KK = KK + 3*NUCDEG(I1)
      ENDIF
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck nuc_nef */
      SUBROUTINE NUC_EFN(IPRINT,DOATOM,EFNNUC,GEOM)
C***********************************************************************
C
C     Calculates nuclear contribution to the electric field at 
C     individual nuclei
C
C     Based on NUC_EFG
C     Written by T.Saue  - Feb 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D3 = 3.0D0)
#include "nuclei.h"
#include "symmet.h"
      LOGICAL    DOATOM(*)
      DIMENSION  GEOM(3*NUCDEP),EFNNUC(3,NUCDEP)
#include "ibtfun.h"
C
C
C     Initialization of GEOM
C     (copied from abacus/hergp.F:CMMASS)
C
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
                  DO KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
                  END DO
                  JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
      CALL DZERO(EFNNUC,3*NUCDEP)
      KK = 1
      NATOM  = 0
      IATOMD = 0
      DO 10 I1 = 1, NUCIND
      IF(DOATOM(I1)) THEN
        DO 20 IDEG1 = 1,NUCDEG(I1)
          NATOM = NATOM + 1
          IATOMD = IATOMD + 1
          LL = 1
          DO 100 I2 = 1, NUCIND
            IF (CHARGE(I2) .EQ. 0.0D0) GOTO 100
            DO 110 IDEG2 = 1,NUCDEG(I2)
              IF (LL .NE. KK) THEN
                XCOOR = GEOM(LL) - GEOM(KK)
                YCOOR = GEOM(LL + 1) - GEOM(KK + 1)
                ZCOOR = GEOM(LL + 2) - GEOM(KK + 2)
                R2 = XCOOR*XCOOR + YCOOR*YCOOR + ZCOOR*ZCOOR
                R3 = R2*SQRT(R2)
                EFNNUC(1,NATOM) = EFNNUC(1,NATOM)- CHARGE(I2)*XCOOR/R3
                EFNNUC(2,NATOM) = EFNNUC(2,NATOM)- CHARGE(I2)*YCOOR/R3
                EFNNUC(3,NATOM) = EFNNUC(3,NATOM)- CHARGE(I2)*ZCOOR/R3
              END IF
              LL = LL + 3
 110        CONTINUE
 100      CONTINUE
          IF (IPRINT .GE. 2) THEN
            WRITE (LUPRI ,'(/2X,2A,/)')
     &        'Nuclear electric field at nucleus: ',NAMDEP(IATOMD)
            IOFF = ICHAR('w')
            DO K = 1,3
               WRITE(LUPRI,'(4X,A1,A,E20.10)') CHAR(IOFF+K),
     &              '-component: ',EFNNUC(K,NATOM)
            ENDDO
          END IF
          KK = KK + 3
 20     CONTINUE
      ELSE
        IATOMD = IATOMD + NUCDEG(I1)      
        KK = KK + 3*NUCDEG(I1)
      ENDIF
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjex3 */
      SUBROUTINE PRJEX3(EXPVAL,IDIMENS,PTRI,PFO,DFO,NCNT,
     &                  NR,NMO,IBUF,PBUF,FMAT,PM,INDX,NTOT,
     &                  INDXPR,IFRP,IPRINT)
C***********************************************************************
C
C     Calculate generalized population matrix
C
C     Written by T. Saue Apr 27 2003
C     Last modifications: MI, March2006/Strasbourg
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0,D0=0.0D0,DTOL=1.0D-8)
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprj.h"
      DIMENSION PTRI(IDIMENS),PFO(NCNT,NCNT,NZ),DFO(NCNT,NCNT,NZ),
     &          NR(*),PBUF(*),IBUF(*),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*)
C
C
C     Make population matrix
C
      NTRI = NCNT*(NCNT+1)/2
      CALL DZERO(PTRI,NTRI)
      DO IZ = 1,NZ
        IJ = 0
        DO J = 1,NCNT
          JR = INDX(1,J)
          DO I = 1,(J-1)
            IR = INDX(1,I)
            IJ = IJ + 1
            TMP = D2*PFO(I,J,IZ)*DFO(I,J,IZ)
            FMAT(IR,JR) = FMAT(IR,JR) + TMP
            PTRI(IJ) = PTRI(IJ) + TMP
          ENDDO
          IJ = IJ + 1
          TMP = PFO(J,J,IZ)*DFO(J,J,IZ)
          PTRI(IJ) = PTRI(IJ) + TMP
          PM(JR) = PM(JR) + TMP
          FMAT(JR,JR) = FMAT(JR,JR) + TMP
        ENDDO
      ENDDO      
      CALL DSCAL(NTRI,D2,PTRI,1)
C      CALL DSCAL(NTOT,D2,PM,1)
      EXPVAL = DSUM(NTRI,PTRI,1)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(A,A3/A,A16,E16.10)') 
     &    '* Contribution from fermion ircop ',FREP(IFRP),
     &    '  to expectation value of property ',PRPNAM(INDXPR),EXPVAL
        WRITE(LUPRI,'(A)') '* Corresponding population matrix:'
        CALL PRMUTC(PTRI,NCNT,1,LUPRI)      
      ENDIF
C
C     Make ordered list of contributions
C
!     IF(IPRINT.GE.1) THEN
! We always want this list /hjaaj Sep 2012
        CALL HEADER('Ordered list of contributions',-1)            
        DO IJ = 1,NTRI
          PBUF(IJ) = -ABS(PTRI(IJ))
        ENDDO
        CALL INDEXX (NTRI,PBUF,IBUF)
        TMP = D0      
        WRITE(LUPRI,'(2(A4,10X),2A16,A10)')
     &    'Ind1','Ind2','Contribution','Accumulated','Fraction'
        DO M = 1,NTRI
          IJ = IBUF(M)
          IF(ABS(PTRI(IJ)/EXPVAL).LT.DTOL) GOTO 10
          J  = INT(SQRT(dble(2*IJ) + 0.25D0) + 0.4999D0)
          I  = IJ - J*(J-1)/2
          JR = INDX(1,J)
          JJ = INDX(2,J)
          IR = INDX(1,I)
          II = INDX(2,I)
          TMP = TMP + PTRI(IJ)
          FAC = TMP/EXPVAL
          WRITE(LUPRI,'(2(A6,I5,3X),1P,2E16.8,E16.5)')
     &          REFFIL(JR),JJ,REFFIL(IR),II,PTRI(IJ),TMP,FAC
        ENDDO
 10     CONTINUE
!     ENDIF
C      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjex4 */
      SUBROUTINE PRJEX4(DFO,BVEC,NFO,NMO,IPRINT)
C***********************************************************************
C
C     Make density matrix in fragment basis
C
C     NFO(IFRP) - number of fragment orbitals in fermion ircop IFRP
C     NMO(IFRP)   - number of molecular orbitals in fermion ircop IFRP
C
C     Written by Trond Saue Oct 21 2003 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION DFO(*),NFO(*),NMO(*),BVEC(*)
C
C     Make density matrix in fragment orbitals
C
      IPD = 1
      IPB = 1
      DO 10 IFRP = 1,NFSYM
        IF(NMO(IFRP).EQ.0) GOTO 10
        CALL DENST1(DFO(IPD),NFO(IFRP),NFO(IFRP),NZ,D1,D0,
     &              BVEC(IPB),NFO(IFRP),NMO(IFRP),
     &              1,NISH(IFRP),NFO(IFRP))  
        ISTART = NISH(IFRP)+1
        DO IOPEN = 1,NOPEN
          NVEC = NACSH(IFRP,IOPEN)
          CALL DENST1(DFO(IPD),NFO(IFRP),NFO(IFRP),NZ,DF(IOPEN),D1,
     &                BVEC(IPB),NFO(IFRP),NMO(IFRP),
     &                ISTART,NVEC,NFO(IFRP))  
          ISTART = ISTART + NVEC
        END DO

        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(A,I3)')
     &    'Density matrix in fragment basis..',IFRP
          CALL PRQMAT(DFO(IPD),NFO(IFRP),NFO(IFRP),
     &         NFO(IFRP),NFO(IFRP),NZ,
     &         IPQTOQ(1,0),LUPRI)
        ENDIF
        IPB = IPB + NFO(IFRP)*NMO(IFRP)*NZ
        IPD = IPD + NFO(IFRP)*NFO(IFRP)*NZ
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Exploc */
      SUBROUTINE EXPCAN(INDXPR,PAO,CMO,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate property canonical occupied orbitals,
C     that is orbitals that diagonalize the matrix of 
C     property INDXPR
C
C     Input:
C       PAO - property matrix in AO-basis
C       CMO - molecular orbital coefficients
C
C     Written by T. Saue Jan 20 2005 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dgroup.h"
      DIMENSION PAO(*),CMO(*),WORK(*)
C
      KFRSAV = KFREE
C
C     Read eigenvalues and boson irrep information
C
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBP,NORBT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','EXPCAN')
      CALL REACMO(LUCOEF,'DFCOEF',DUM,WORK(KEIG),WORK(KIBP),TOTERG,12,
     &               WORK,KFREE,LFREE)
      DO IFRP = 1,NFSYM
C       Memory allocation
        CALL MEMGET('REAL',KPMO,NISH(IFRP)*NISH(IFRP)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSCO,NISH(IFRP)*NFBAS(IFRP,0)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUCO,NISH(IFRP)*NISH(IFRP)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPEI,NISH(IFRP),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KSEL,NISH(IFRP),WORK,KFREE,LFREE)
C       Select occupied molecular orbitals
        CALL OCCUPY(WORK(KSEL),NISH(IFRP))
        CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KSCO),NISH(IFRP),
     &              WORK(KSEL),D0,NISH(IFRP),NFBAS(IFRP,0),NORB(IFRP))
        CALL MEMREL('EXPCAN.sel',WORK,KSEL,KSEL,KFREE,LFREE)
C       Generate canonical property orbitals
        CALL EXPCA1(INDXPR,IFRP,PAO(I2BASX(IFRP,IFRP)+1),
     &              WORK(KPMO),CMO(ICMOQ(IFRP)+1),WORK(KEIG+IORB(IFRP)),
     &              WORK(KSCO),WORK(KUCO),WORK(KPEI),
     &              IPRINT,WORK,KFREE,LFREE)
C      
        CALL MEMREL('EXPCAN.1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      ENDDO
C     Write new coefficients
      IUNIT = 2
      CALL OPNFIL(IUNIT,PRPNAM(INDXPR),'NEW','EXPCAN')  
      CALL WRICMO(IUNIT,CMO,WORK(KEIG),WORK(KIBP),TOTERG)
      CLOSE(IUNIT,STATUS='KEEP')
C     Read old coefficients
      REWIND LUCOEF
      CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IDUM,TOTERG,2,
     &               WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS='KEEP')
      CALL MEMREL('EXPCAN.2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Explo1 */
      SUBROUTINE EXPCA1(INDXPR,IFRP,PAO,PMO,CMO,EIG,SCO,UCO,PEIG,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate property canonical occupied orbitals,
C     that is orbitals that diagonalize the matrix of 
C     property INDXPR
C
C     Input:
C       PAO - property matrix in AO-basis
C       CMO - molecular orbital coefficients
C
C     Written by T. Saue Jan 20 2005 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      CHARACTER MXFORM*6,FMT*6
      DIMENSION PAO(*),PMO(*),SCO(*),UCO(*),CMO(*),EIG(*),
     &          PEIG(*),WORK(*)
C
      KFRSAV = KFREE
C
C     Generate property matrix
C
      CALL QTRANS('AOMO','S',D0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &            NISH(IFRP),NISH(IFRP),
     &            PAO,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            PMO,NISH(IFRP),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,*) '* Property matrix:'
        CALL PRQMAT(PMO,NISH(IFRP),NISH(IFRP),NISH(IFRP),NISH(IFRP),
     &            NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Diagonalize
C
      CALL QDIAG(NZ,NISH(IFRP),PMO,NISH(IFRP),NISH(IFRP),
     &           PEIG,1,UCO,NISH(IFRP),NISH(IFRP),
     &           WORK(KFREE),LFREE,IERR)
C
C     Write eigenvalues
C
      WRITE(LUPRI,*) '* Eigenvalues:'
      TMP = D0
      DO I = 1,NISH(IFRP)
        FMT = MXFORM(PEIG(I),15)
        WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',1X,A4)')
     &        FREP(1),I,PEIG(I),'a.u.'
      ENDDO
      TMP = DSUM(NISH(IFRP),PEIG,1)
      FMT = MXFORM(TMP,15)
      WRITE(LUPRI,'(4X,A8,3X,'//FMT//',1X,A4)')
     &     'Total : ',TMP,'a.u.'
      IOFF = NPSH(IFRP)+1
      CALL DCOPY(NISH(IFRP),PEIG,1,EIG(IOFF),1)
C
C     Write eigenvectors
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* Eigenvectors:'
        CALL PRQMAT(UCO,NISH(IFRP),NISH(IFRP),NISH(IFRP),NISH(IFRP),NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Backtransform
C
      IOFF = NFBAS(IFRP,0)*NPSH(IFRP)+1
      CALL QGEMM(NFBAS(IFRP,0),NISH(IFRP),NISH(IFRP),D1,
     &           'N','N',IPQTOQ(1,0),SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,
     &           'N','N',IPQTOQ(1,0),UCO,NISH(IFRP),NISH(IFRP),NZ,
     &           D0,IPQTOQ(1,0),CMO(IOFF),NFBAS(IFRP,0),NORB(IFRP),NZ)
      RETURN
      END
