C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dirgrd.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psigrd */
      SUBROUTINE PSIGRD(WORK,LWORK)
C*****************************************************************************
C
C     calculate gradient
C
C     Written by J. Thyssen - 1997/05
C
C     Last revision: jth 1997/07/09
C
C*****************************************************************************
#include <implicit.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <mxorb.h>
#include <mxcent.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL   LSAVE(3)
#include <infpar.h>
#include <ccom.h>
#include <dcbgen.h>
#include <dcbdhf.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <blocks.h>
#include <dcbgrd.h>
C
      CALL QENTER('PSIGRD')
#include <memint.h>
C
C     Initialize
C     ==========
C
C     If calculation of gradient requested without
C     calculating the wavefunction SETDHF must be called -
C     so we call it anyway.
C
      CALL SETDHF(IPRGRD)
C
C
C
C     Save DO2GRD
C
      DO I = 1,3
         LSAVE(I) = DO2GRD(I)
      END DO
C
C
C     Memory allocation
C
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,       WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),       WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT,N2BBASX*4*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,N2BBASX*4*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2TMT          ,WORK,KFREE,LFREE)
C
      CALL TITLER('Calculation of molecular gradient','*',110)
C
C     Get hold of the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = 1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
C     Calculate gradient
C
      CALL DHFGRD(WORK(KCMO),WORK(KDMAT),WORK(KFOCK),
     &            WORK(KTMAT),WORK(KFREE),LFREE,IPRGRD)
C
C     Release the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = -1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
C     Output of gradient
C
      CALL GRDOUT(IPRGRD)
C
C
C     Restore DO2GRD
C
      DO I = 1,3
         DO2GRD(I) = LSAVE(I)
      END DO
C
C
C     Memory deallocation
      CALL MEMREL('PSIGRD.DHFOUT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('PSIGRD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dhfgrd */
      SUBROUTINE DHFGRD(CMO,DMAT,FOCK,TMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Calculate molecular gradient
C
C     J. Thyssen, 1997/18/04
C    
C     Revision: 1997/11/11 - Parallel calculations is now possible
C     Remember: Thresholds
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER(D0 = 0.0D0)
C
      DIMENSION FOCK(*),CMO(*),TMAT(*),
     &          DMAT(N2BBASX,NZ,NFMAT),WORK(LWORK)
      DIMENSION DNTSKP(8)
      DIMENSION IQBF(4)
      CHARACTER*1 CTMP
      CHARACTER*3 STR
      LOGICAL   DOINT(2,2)
      LOGICAL   DOLS,DOSS,LBIT
      LOGICAL   NODV,NOPV,FIRST,NOCONT,ICEDIF
C     Timing variables
      DOUBLE PRECISION DPRTIM(4,2)
      SAVE DPRTIM
      DATA DPRTIM /D0,D0,D0,D0,D0,D0,D0,D0/
C
#include <dcbbas.h>
#include <dcbqal.h>
#include <nuclei.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbisol.h>
#include <dcbgen.h>
#include <dcbpsi.h>
#include <dcbdhf.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dorps.h>
C two el. part
#include <ccom.h>
#include <cbihr2.h>
#include <blocks.h>
#include <dcbgrd.h>
C
#include <exeinf.h>
#include <symmet.h>
C
      CALL QENTER('DHFGRD')
C
#include <memint.h>
C
      CALL GETTIM(CPUTO1,WLLTO1)
      MAXDER = 1
C  
C     Initialisation of timings
C
      CPUONE = D0
      CPUTLL = D0
      CPUTLS = D0
      CPUTSS = D0
      CPUTOT = D0
      WLLONE = D0
      WLLTLL = D0
      WLLTLS = D0
      WLLTSS = D0
      WLLTOT = D0
C
C    
      DO I = 1,2
         DO J = 1,2
            DOINT(I,J) = .TRUE.
         END DO
      END DO
C
      IF ( DOLVC .AND. DO2GRD(3) ) THEN
         DO2GRD(3) = .FALSE.
         WRITE(LUPRI,'(/XA)') 'INFO: LVCORR used in wavefunction, but'//
     &                        ' not in gradient =>'
         WRITE(LUPRI,'(XA/)') 'LVCORR gradient used (SS-gradient is'//
     &                        ' not calculated)'
      END IF
C
C     Do totally symmetric contributions
      DOREPS(0) = .TRUE.
      CALL SETSYM
C
      DO I = 1, 3*MXCENT
         DOPERT(I,1) = .TRUE.
         DOPERT(I,2) = .TRUE.
      END DO
      FTRONV = .TRUE.
      RELGRD = .TRUE.
      CALL TROINI
      CALL TROINV(WORK(KFREE),LWORK)
      RELGRD = .FALSE.
C
C     Nuclear repulsion
C     =================
C
      RELGRD = .TRUE.
      CALL NUCINI
      CALL NUCREP
      RELGRD = .FALSE.
C
C     "SS-integral contribution by classical repulsion of 
C      small component atomic charges"-gradient
C
      CALL LVCGRD(IPRINT,WORK(KFREE),LWORK)
C    
C     One electron part of gradient
C     =============================
C
      CALL HEADER('DHFGRD: One-electron part',-1)
      IF ( DPRTIM(1,1) .GT. D0 ) THEN
         CALL TIMDAT
         WRITE(LUPRI,'(/1X,A,F10.4,A)') 
     &               'Estimated CPU time:  ',DPRTIM(1,1),' seconds.'
         WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &               'Estimated WALL time: ',DPRTIM(1,2),' seconds.'
      END IF
C
C
      CALL GETTIM(CPU1,WALL1)
C
C     Get contravariant Fock matrix and density matrix
C     in symmetry distinct AO-basis
C
      CALL GTCNVF(FOCK,DMAT,CMO,IPRINT,WORK(KFREE),LWORK)
C
C
C     Call ONEDRV for construction of gradient
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
C
C     PROPTY = .TRUE.
C     MAXDIF = 1
C     DIFDIP = .FALSE.
C     SECDER = .FALSE.
C     NOPV   = .FALSE.
C     NOCV   = .FALSE.
C
C
C     Make one-electron part of gradient
C
      CALL ONEDR1(WORK(KSTHMA),DMAT,FOCK,WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),DUMMY,WORK(KFREE),LWRK,
     &            IPRINT,.TRUE.,1,.TRUE.,.TRUE.,.FALSE.,NTBAS(0),
     &            NNBBASX,NNBBAST,DOINT,.FALSE.)
C
C     Dispose temp. memory
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUONE = CPU2 - CPU1
      WLLONE = WALL2 - WALL1
      WRITE(LUPRI,'(/A,F10.4,A/)') 
     &      'Time for one-electron part of gradient    :',
     &      CPUONE,' seconds.'
      CALL FLSHFO(LUPRI)
C
C
C     Two electron part of gradient
C     =============================
C
      CALL HEADER('DHFGRD: Two-electron part',-1)
      CALL GETTIM(CPU11,WALL11)
C
C
      CALL DZERO(GRADEE,MXCOOR)
      CALL DZERO(GRADER(1,1),MXCOOR)
      CALL DZERO(GRADER(1,2),MXCOOR)
      CALL DZERO(GRADER(1,3),MXCOOR)
C
C     Get contravariant density matrix, no folding
C
      CALL GTCNVD(DMAT,CMO,IPRINT,WORK(KFREE),LWORK)
C
C
C     Prepare screening
C
C     AO-integrals
CC      IJOB   = 0
CC      IGTYP  = 1
C     Differentiated integrals
CC      ITYPE  = 2
CC      MAXDIF = 1
C     Misc.
CC      JATOM  = 0
CC      NOCONT = .FALSE.
CC      CALL GABGEN(IJOB,ITYPE,IGTYP,MAXDIF,JATOM,NOCONT,
CC     &            100,WORK(KFREE),LFREE)
CCC     &            IPRINT,WORK(KFREE),LFREE)
C
      IF ( DO2GRD(1) .AND. .NOT. LBIT(INTFLG,1)) THEN
         WRITE(LUPRI,'(/1X,A/)') 'WARNING: LL-integrals not used in '//
     &                  ' SCF calculation => '//
     &                  ' LL-Two-electron gradient .NOT. calculated!'
         DO2GRD(1) = .FALSE.
      END IF
C
      CALL CYESNO(STR,DO2GRD(1))
      WRITE(LUPRI,'(1X,2A/)') 'Calculate LL-Two-electron gradient: ',
     &                        STR
C
C
      IF ( DO2GRD(1) ) THEN
         CALL HEADER('DHFGRD: LL-Two-electron part',-1)
         IF ( DPRTIM(2,1) .GT. D0 ) THEN
            CALL TIMDAT
            WRITE(LUPRI,'(/1X,A,F10.4,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(2,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(2,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('DHFGRD: Skipping LL-Two-el.',-1)
         GOTO 100
      END IF
C
C
      CALL GETTIM(CPU1,WALL1)
C
C     LL-gradient
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      I2TYP  = 1
      INTCLASS  = 0
      ICEDIF = .FALSE.
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
C
      NEWBAS = .TRUE. 
      RELGRD = .TRUE.
      NDMAT   = NFMAT*4
      NDMRAO  = NSYMBL*NSYMBL*NDMAT
      CALL MEMGET('REAL',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
      LWRK = LWORK - KFREE + 1
      CALL DZERO(DNTSKP,8)
#if defined (VAR_MPI) || defined (VAR_PVM)
      IF (PARCAL) THEN
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL PARDRV(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &               RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &               IFTHRS,DUMMY,WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &               NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,IFTHRS,
     &               DUMMY,WORK(KDMRAO),DUMMY,DNTSKP,RELCAL)
#if defined (VAR_MPI) || defined (VAR_PVM)
      END IF
#endif
      RELGRD = .FALSE.
C
C
C     Dispose memory
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,IPRINT)
C
C     Copy LL-GRADEE gradient to GRADER(1,1)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,1),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C     Memory deallocation
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTLL = CPU2 - CPU1
      WLLTLL = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTLL = MAX(CPUTLL,DMXCPT)
         WLLTLL = MAX(WLLTLL,DMXWLT)
      END IF
      WRITE(LUPRI,'(/1X,A,F10.4,A/)') 
     &      'Time for LL-two-electron part of gradient :',
     &      CPUTLL,' seconds.'
C
      CALL FLSHFO(LUPRI)
 100  CONTINUE
C
C
C     Calculate estimated LS- and SS- gradient and skip if
C     they seem unnecessary.
C
C
      IF ( DOTRCK ) THEN
         CALL ESTGLS(DOLS)
         DO2GRD(2) = DO2GRD(2) .AND. DOLS
      END IF
C
C     Some output...
      IF ( DO2GRD(2) .AND. .NOT. LBIT(INTFLG,2) ) THEN
         WRITE(LUPRI,'(/1X,A/)') 'WARNING: LS-integrals not used in '//
     &                  ' SCF calculation => '//
     &                  ' LS-Two-electron gradient .NOT. calculated!'
         DO2GRD(2) = .FALSE.
      END IF
      CALL CYESNO(STR,DO2GRD(2))
      WRITE(LUPRI,'(/1X,2A/)') 'Calculate LS-Two-electron gradient: ',
     &                        STR
C
      IF ( DO2GRD(2) ) THEN
         CALL HEADER('DHFGRD: LS-Two-electron part',-1)
         IF ( DPRTIM(3,1) .GT. D0 ) THEN
            CALL TIMDAT
            WRITE(LUPRI,'(/1X,A,F10.4,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(3,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(3,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('DHFGRD: Skipping LS-Two-el.',-1)
         GOTO 200
      END IF
      CALL GETTIM(CPU1,WALL1)
C
C
C     LS-gradient
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      I2TYP  = 2
      INTCLASS  = 1
      ICEDIF = .FALSE.
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
      NEWBAS = .TRUE.
      RELGRD = .TRUE.
      NDMAT   = NFMAT*4
      NDMRAO  = NSYMBL*NSYMBL*NDMAT
      CALL MEMGET('REAL',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
      LWRK = LWORK - KFREE + 1
      CALL DZERO(DNTSKP,8)
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI) || defined (VAR_PVM)
      IF (PARCAL) THEN
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL PARDRV(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &               RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &               IFTHRS,DUMMY,WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &               NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,IFTHRS,
     &               DUMMY,WORK(KDMRAO),DUMMY,DNTSKP,RELCAL)
#if defined (VAR_MPI) || defined (VAR_PVM)
      END IF
#endif
      RELGRD = .FALSE.
C
C     Dispose memory
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,IPRINT)
C     Copy LS-GRADEE gradient to GRADER(1,2)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,2),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C     Dispose memory
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTLS = CPU2 - CPU1
      WLLTLS = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTLS = MAX(CPUTLS,DMXCPT)
         WLLTLS = MAX(WLLTLS,DMXWLT)
      END IF
      WRITE(LUPRI,'(/1X,A,F10.4,A/)') 
     &      'Time for LS-two-electron part of gradient :',
     &      CPUTLS,' seconds.'
C
C
C
C
      CALL FLSHFO(LUPRI)
 200  CONTINUE
C
      IF ( DOTRCK .AND. .NOT. DOLVC) THEN
         CALL ESTGSS(DOSS)
         DO2GRD(3) = DOSS .AND. DO2GRD(3)
      END IF
C
      IF ( DO2GRD(3) .AND. .NOT. DO2GRD(2) ) THEN
         WRITE(LUPRI,'(/1X,A/)') 'WARNING: SS-gradient .NOT. '//
     &       'calculated since LS-gradient is not calculated.'
         DO2GRD(3) = .FALSE.
      END IF
      CALL CYESNO(STR,DO2GRD(3))
      WRITE(LUPRI,'(1X,2A/)') 'Calculate SS-Two-electron gradient: ',
     &                        STR
C
      IF ( DO2GRD(3) ) THEN
         CALL HEADER('DHFGRD: SS-Two-electron part',-1)
         IF ( DPRTIM(4,1) .GT. D0 ) THEN
            CALL TIMDAT
            WRITE(LUPRI,'(/1X,A,F10.4,A)') 
     &                  'Estimated CPU time:  ',DPRTIM(4,1),' seconds.'
            WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &                  'Estimated WALL time: ',DPRTIM(4,2),' seconds.'
         END IF
      ELSE
         CALL HEADER('DHFGRD: Skipping SS-Two-el.',-1)
         GOTO 300
      END IF
      CALL GETTIM(CPU1,WALL1)
C
C
C     SS-gradient
C
C     IREPDM = IDUMMY
      ITYPE  = 2
C     GMAT   = DUMMY
C     INDX   = IDUMMY
C     INDXAB = IDUMMY
      MAXDIF = 1
      JATOM  = 0
      NODV   = (NFMAT.EQ.1)
      NOPV   = .TRUE.
      NOCONT = .FALSE.
      FIRST  = .TRUE.
C     NPOS   = IDUMMY
      MTOTTK = NSMLBL*(NSMLBL+1)/2
      I2TYP  = 3
      INTCLASS  = 2
      ICEDIF = .FALSE.
C     GABRAO = DUMMY
C     DMRSO  = DUMMY
C
C
      RELGRD = .TRUE.
      NEWBAS = .TRUE.
      NDMAT   = NFMAT*4
      NDMRAO  = NSYMBL*NSYMBL*NDMAT
      CALL MEMGET('REAL',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
      LWRK = LWORK - KFREE + 1
      CALL DZERO(DNTSKP,8)
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI) || defined (VAR_PVM)
      IF (PARCAL) THEN
         CALL MEMGET('INTE',KNPOS,MTOTTK,WORK,KFREE,LFREE)
         CALL PARDRV(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               ITYPE,MAXDIF,JATOM,NODV,NOPV,NOCONT,TKTIME,
     &               RTNTWO,FIRST,WORK(KNPOS),MTOTTK,I2TYP,ICEDIF,
     &               IFTHRS,DUMMY,WORK(KDMRAO),DUMMY,DNTSKP)
      ELSE
#endif
         CALL TWOINT(WORK(KFREE),LWRK,FOCK,DMAT,NDMAT,IDUMMY,IFCTYP,
     &               DUMMY,IDUMMY,IDUMMY,ITYPE,MAXDIF,JATOM,NODV,
     &               NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,IDUMMY,I2TYP,ICEDIF,IFTHRS,
     &               DUMMY,WORK(KDMRAO),DUMMY,DNTSKP,RELCAL)
#if defined (VAR_MPI) || defined (VAR_PVM)
      END IF
#endif
      RELGRD = .FALSE.
C
C     Dispose memory
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Symmetrization of GRADEE
C
      CALL MEMGET('REAL',KSKLTN,MXCOOR*2           ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAC  ,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDEX,3*NUCIND*(MAXREP+1),WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPOINT,3*NUCDEP*(MAXREP+1),WORK,KFREE,LFREE)
      CALL SKLEXP(WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),WORK(KPOINT),
     &               MAXDER,IPRINT)
C     Copy SS-GRADEE gradient to GRADER(1,3)
C
      CALL DCOPY(MXCOOR,GRADEE,1,GRADER(1,3),1)
      CALL DZERO(GRADEE,MXCOOR)
C
C
C     End of gradient calculation
C     ===========================
C
C     Memory deallocation
C
      CALL MEMREL('DHFGRD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL GETTIM(CPU2,WALL2)
C
      CPUTSS = CPU2 - CPU1
      WLLTSS = WALL2 - WALL1
      IF ( PARCAL ) THEN
         CPUTSS = MAX(CPUTSS,DMXCPT)
         WLLTSS = MAX(WLLTSS,DMXWLT)
      END IF
      WRITE(LUPRI,'(/1X,A,F10.4,A/)') 
     &      'Time for SS-two-electron part of gradient :',
     &      CPUTSS,' seconds.'
C
      CALL FLSHFO(LUPRI)
 300  CONTINUE
C
C
C     Output of timings
C
      CALL HEADER('Timings for DHFGRD',-1)
      WRITE(LUPRI,'(/1X,A,F10.4,A)') 
     &     'One-electron part of gradient    CPU  :',CPUONE,' seconds.'
      WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &     '                                 WALL :',WLLONE,' seconds.'
      IF ( DO2GRD(1) ) THEN
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      'LL-Two-electron part of gradient CPU  :',CPUTLL,' seconds.'
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      '                                 WALL :',WLLTLL,' seconds.'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &      'LL-Two-electron part SKIPPED!'
      END IF
      IF ( DO2GRD(2) ) THEN
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      'LS-Two-electron part of gradient CPU  :',CPUTLS,' seconds.'
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      '                                 WALL :',WLLTLS,' seconds.'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &      'LS-Two-electron part SKIPPED!'
      END IF
      IF ( DO2GRD(3) ) THEN
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      'SS-Two-electron part of gradient CPU  :',CPUTSS,' seconds.'
         WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      '                                 WALL :',WLLTSS,' seconds.'
      ELSE
         WRITE(LUPRI,'(1X,A)') 
     &      'SS-Two-electron part SKIPPED!'
      END IF
      WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      'Total two-electron part of grad. CPU  :',
     &      CPUTSS+CPUTLS+CPUTLL,' seconds.'
      WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &      '                                 WALL :',
     &      WLLTSS+WLLTLS+WLLTLL,' seconds.'
      CALL GETTIM(CPUTO2,WLLTO2)
      CPUTOT = CPUTO2 - CPUTO1
      WLLTOT = WLLTO2 - WLLTO2
      WRITE(LUPRI,'(1X,A,F10.4,A)') 
     &      'Total time in DHFGRD                  :',
     &      CPUTOT,' seconds.'
      WRITE(LUPRI,'(1X,A,F10.4,A/)') 
     &      '                                 WALL :',
     &      WLLTOT,' seconds.'
C
C     Save timings for later use
C
      DPRTIM(1,1) = CPUONE
      DPRTIM(2,1) = CPUTLL
      DPRTIM(3,1) = CPUTLS
      DPRTIM(4,1) = CPUTSS
      DPRTIM(1,2) = WLLONE
      DPRTIM(2,2) = WLLTLL
      DPRTIM(3,2) = WLLTLS
      DPRTIM(4,2) = WLLTSS
C
      CALL QEXIT('DHFGRD')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setsym */
      SUBROUTINE SETSYM
C*****************************************************************************
C
C     Temporary subroutine: sets DOSYM-array from abainf.h
C
C     jth, 1997/24/04
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <abainf.h>
      DO I=1,8
         DOSYM(I) = .TRUE.
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcnvf */
      SUBROUTINE GTCNVF(FOCK,DMAT,CMO,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver routine to GTCNV1
C
C     jth 1997/05/06
C     Last revision: 1997/05/06
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION FOCK(*),WORK(*),DMAT(*),CMO(*)
      CHARACTER*1 CTMP
C
#include <dcbbas.h>
#include <dcbgen.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dummy.h>
C
#include <memint.h>
C
C
      CALL QENTER('GTCNVF')
      IF ( IPRINT .GE. 5 ) CALL HEADER('Output from GTCNVF',-1)
C
C     Allocate memory for temporary memory for use in GTCNV1
C
      CALL MEMGET('REAL',KTMP1,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP2,N2BBASXQ*(NFMAT+1),WORK,KFREE,LFREE)
C
C     Read coefficient matrix and eigenvalues
C     ===================================================
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','DHFGRD')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,DUMMY,TOTERG,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Print section
C
      IF ( IPRINT .GE. 10 ) THEN
         DO I=1,NFSYM
           IF (NORB(I) .NE. 0) THEN
           WRITE(CTMP,'(I1)') I
           CALL HEADER('Coefficient matrix, Ircop '//CTMP,-1)
           CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           END IF
         END DO
      END IF
C
      LWRK   = LWORK - KFREE + 1
C
      CALL GTCNV1(FOCK,DMAT,WORK(KTMP1),WORK(KTMP2),CMO,
     &            IPRINT,WORK(KFREE),LWRK)
C
C     Dispose memory
C
      CALL MEMREL('GTCNVF',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNVF')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcnv1 */
      SUBROUTINE GTCNV1(FOCK,DMAT,DTMP1,DTMP2,CMO,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get contravariant Fock matrix in AO-basis and density matrix
C
C     jth 1997/05/02, ideas based on DSOFS1 by PRT and TUH
C     Last revision: 1997/06/19 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <dcbbas.h>
#include <dcbdhf.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbgen.h>
      PARAMETER ( D0 = 0.0 D00, D1 = 1.0 D00, D2 = 2.0 D00 ,
     &            LUBUF = 22, DP5 = 0.50D00)
C
      DIMENSION FOCK(*),CMO(*),WORK(*)
      DIMENSION DTMP2(N2BBASX,NZ,*),DMAT(*),DTMP1(N2BBASXQ)
      DIMENSION IQBF(4)
      LOGICAL   TOBE
C
      IFULL(I,J)     = I + (J-1) * NTBAS(0)
      IQOFF(IZ,IMAT) = (IZ-1) * N2BBASX + (IMAT-1) * N2BBASXQ
C
#include <memint.h>
C
      IP = IPRINT
C      IPRINT = 20
      CALL QENTER('GTCNV1')
      IF ( IPRINT .GT. 5 ) CALL HEADER('OUTPUT FROM GTCNV1',-1)
C
C     FOCK MATRIX 
C     ===========
C
C
      CALL MEMGET('REAL',KFMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDTMP3,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMOAO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2TMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KFMO),N2BBASXQ)
C
C     Get AO Fock matrices
C
C     Two-electron Fock matrix/matrices
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)
      IF(TOBE) THEN
         CALL OPNFIL(LUBUF,'DFFCK2','OLD','GETFC1')
         CALL REAFCK(LUBUF,DTMP2,.TRUE.,NFMAT)
         CLOSE(LUBUF,STATUS='KEEP')
      ELSE
         WRITE(LUPRI,'(A)')
     &   '*** ERROR: Gradient calc. without DFFCK2 '//
     *   'is not implemented yet'
         CALL QUIT('***ERROR: no DFFCK2')
      END IF
C
C     One electron Fock Matrix
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF (TOBE) THEN
         CALL OPNFIL(LUBUF,'DFFCK1','OLD','GTCNV1')
         CALL REAFCK(LUBUF,DTMP1,.TRUE.,1)
         CLOSE(LUBUF,STATUS='KEEP')
      ELSE
         WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GETFCK: No 1-Fock found. Regenerating.'
         CALL ONEFCK(DTMP1,IPRINT,WORK(KFREE),LFREE)
      END IF
C
C     FC = H1 + QC
C
      CALL DAXPY(N2BBASXQ,D1,DTMP1,1,DTMP2,1)
      CALL DZERO(WORK(KFMO),N2ORBTQ)
C
      IF (NASHT.GT.0) THEN
C
C        Transform FC to MO basis
C
         DO 10 I = 1,NFSYM
            IF (NASH(I).EQ.0) GOTO 10
            CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                  NORB(I),NORB(I),
     &                  DTMP2(I2BASX(I,I)+1,1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,0),
     &                  DTMP1(1+I2ORBT(I)),NORB(I),NORB(I),
     &                  NZA,IQBF,
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,IPRINT)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'FC in MO-basis, corep ',I
               CALL PRQMAT(DTMP1(1+I2ORBT(I)),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IADR = 1 + I2ORBT(I)+(NISH(I)+NPSH(I))*NORB(I)
            CALL DCOPY(NASH(I)*NORB(I),DTMP1(IADR),1,
     &                                 WORK(KFMO+IADR-1),1)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'Active part of F, corep ',I
               CALL PRQMAT(WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
 10      CONTINUE
C
C        FD = H1 + QC + QV 
C           = FC + QV
C
         CALL DAXPY(N2BBASXQ,DP5,DTMP2(1,1,2),1,DTMP2,1)
      END IF
C
      IF (NISHT.GT.0) THEN
         DO 20 I = 1,NFSYM
            IF (NISH(I).EQ.0) GOTO 20
            CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                  NORB(I),NORB(I),
     &                  DTMP2(I2BASX(I,I)+1,1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,0),
     &                  DTMP1(1+I2ORBT(I)),NORB(I),NORB(I),
     &                  NZA,IQBF,
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,IPRINT)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'FD in MO-basis, corep ',I
               CALL PRQMAT(DTMP1(1+I2ORBT(I)),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
            IADR = 1 + I2ORBT(I)+NPSH(I)*NORB(I)
            CALL DSCAL(NISH(I)*NORB(I),D2,DTMP1(IADR),1)
            CALL DCOPY(NISH(I)*NORB(I),DTMP1(IADR),1,
     &                                 WORK(KFMO+IADR-1),1)
            IF (IPRINT.GE.10) THEN
               WRITE(LUPRI,*) 'Total Fock matrix, corep ',I
               CALL PRQMAT(WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &                    NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            END IF
 20      CONTINUE
      END IF
C
C     Construct contravariant Fock matrix in AO-basis
C     Use DTMP2 as temporarily Fock matrix
C
C     FOCK = CCtFCCt = DFDt = CF(MO)Ct
C
      CALL DZERO(DTMP2,N2BBASXQ)
C
      DO I= 1,NFSYM
         IF (NORB(I) .NE. 0) THEN
            NZA = NZ
            CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                  NORB(I),NORB(I),
     &                  DTMP2(I2BASX(I,I)+1,1,1),NTBAS(0),NTBAS(0),
     &                  NZA,IQBF,
     &                  WORK(KFMO+I2ORBT(I)),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,
     &                  IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,IPRINT)
         END IF
      END DO
C
C     Print section
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('GTCNV1: FOCK After QTRANS to AO-basis',-1)
         CALL FLSHFO(LUPRI)
         CALL PRQMAT(DTMP2,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Transform to sorted basis
C
      CALL BSTOBU(DTMP2,NZ,WORK(KFREE),LFREE)
      IF(IPRINT.GE.5) THEN
         CALL HEADER('GTCNV1: FOCK After BSTOBU to unsorted basis',-1)
         CALL PRQMAT(DTMP2,NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Fold Fock matrix
C
      DO IZ = 1,NZ
         IF ( IZ .EQ. 1 ) THEN
            DSCALE = 1.0
         ELSE
            DSCALE = -1.0
         END IF
         DO IMU = 2,NTBAS(0)
            DO INU = 1,IMU-1
               IADR1  = IFULL(IMU,INU)
               IADR2  = IFULL(INU,IMU)
               DTMP2(IADR1,IZ,1) = DTMP2(IADR1,IZ,1) 
     &                           + DSCALE * DTMP2(IADR2,IZ,1)
               DTMP2(IADR2,IZ,1) = D0
            END DO
         END DO
      END DO
C
C
      IF ( IPRINT .GT. 5 ) THEN
         CALL HEADER('GTCNV1: FOCK After folding',-1)
         CALL PRQMAT(DTMP2,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &            NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Insert half-phases
C     ==================
C
      IF(NZ.LT.4) THEN
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,0)
            CALL Q2PHASE('F',IQ,1,DTMP2(1,IZ,1))
         ENDDO
         IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted FOCK with half-phases inserted',-1)
            CALL PRQMAT(DTMP2,NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZ,
     &                  IPQTOQ(1,0),LUPRI)
         ENDIF
      ENDIF
C
C
C     Transform fock matrix from QO to AO
C     ======================================
C
      CALL DZERO(FOCK,4*N2BBASX)
      IOFF = 0
      IPAR = JBTOF(0,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,0)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,0)
         IF(IPQ.GT.NZ) CALL QUIT('GTCNV1.FOCK: IPQ.GT.NZ !')
         IF(IPQ.LE.0 ) CALL QUIT('GTCNV1.FOCK: IPQ.LE.0 !')
         CALL DSYM1R(DTMP2(1,IPQ,1),FOCK(IOFF+1),
     &               NTBAS(0),IREPD,IPRINT)
         IOFF = IOFF + N2BBASX
      ENDDO
C
C
C     Construction of the density matrix
C     ==================================
C
C     Use DTMP2 temporarily density matrix
C
C     Since Trond insists in D_ij = \sum_a C_ia C_ja^*
C     and .NOT. D_ij = 2 \sum C_ia C_ja^*, we multiply DTMP2 with 2.0
C
C
C     Get density matrix
C
      CALL GENDEN(DTMP2,CMO,1,IPRINT)
      IF (NASHT.GT.0) THEN
         CALL GENDEN(DTMP2(1,1,2),CMO,2,IPRINT)
C
C        DVao + DCao
C
         CALL DAXPY(N2BBASXQ,DP5,DTMP2(1,1,2),1,
     &                          DTMP2,1)
      END IF
      CALL DSCAL(N2BBASXQ,D2,DTMP2,1)
C
C     Transform to sorted basis
C
      CALL BSTOBU(DTMP2(1,1,1),NZ,WORK(KFREE),LFREE)
      IF(IPRINT.GE.5) THEN
         CALL HEADER('GTCNV1: DMAT After BSTOBU to unsorted basis',
     &               -1)
         CALL PRQMAT(DTMP2(1,1,1),NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),NZ,
     &               IPQTOQ(1,0),LUPRI)
      END IF
C
C     Fold density matrix
C
      DO IZ = 1,NZ
         IF ( IZ .EQ. 1 ) THEN
            DSCALE = 1.0
         ELSE
            DSCALE = -1.0
         END IF
         DO IMU = 2,NTBAS(0)
            DO INU = 1,IMU-1
               IADR1 = IFULL(IMU,INU)
               IADR2 = IFULL(INU,IMU)
               DTMP2(IADR1,IZ,1) = DTMP2(IADR1,IZ,1) +
     &                            DSCALE * DTMP2(IADR2,IZ,1)
               DTMP2(IADR2,IZ,1) = D0
            END DO
         END DO
      END DO
C
C
      IF ( IPRINT .GT. 5 ) THEN
         CALL HEADER('GTCNV1: DMAT after folding',-1)
         CALL PRQMAT(DTMP2(1,1,1),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Transform density matrix to AO basis
C
C
C     Insert half-phases
C     ==================
C
      IREP = ISYMOP(1)-1
      IF(NZ.LT.4) THEN
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2PHASE('D',IQ,1,DTMP2(1,IZ,1))
         ENDDO
         IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted DMAT with half-phases inserted',-1)
            CALL PRQMAT(DTMP2(1,1,1),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         ENDIF
      ENDIF
C
C
C       Transform density matrix from QO to AO
C       ======================================
C
      CALL DZERO(DMAT,4*N2BBASX)
      IOFF = 0
      IPAR = JBTOF(0,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,0)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,0)
         IF(IPQ.GT.NZ) CALL QUIT('GTCNV1.DENSTY: IPQ.GT.NZ !')
         IF(IPQ.LE.0 ) CALL QUIT('GTCNV1.DENSTY: IPQ.LE.0 !')
         CALL DSYM1R(DTMP2(1,IPQ,1),DMAT(IOFF+1),
     &              NTBAS(0),IREPD,IPRINT)
         IOFF = IOFF + N2BBASX
      ENDDO
C
C
C
C     Final print section
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('GTCNV1: Final matrices.',-1)
         CALL HEADER('Density matrix:',-1)
         CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &               NTBAS(0),4,IQMULT,LUPRI)
         CALL HEADER('Contravariant fock matrix',-1)
         CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),
     &               NTBAS(0),4,IQMULT,LUPRI)
      END IF
C
C     Dispose used memory
C
      CALL MEMREL('GTCNV1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNV1')
      CALL FLSHFO(LUPRI)
C
      IPRINT = IP
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dsym1r */
      SUBROUTINE DSYM1R(ASO,AAO,NBAST,IREPDM,IPRINT)
C*****************************************************************************
C
C     Matrix transform: SO-basis to AO-basis
C     Based on HERMIT routine DSOTAO (880418  PRT)
C
C     Slightly changed to DSYM1R by jth 1997/05/09
C     Written by T.Saue - Sept 10 1995
C     Last revision: Sep 10 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
      DIMENSION ASO(NBAST,NBAST), AAO(NBAST,NBAST)
      DIMENSION IR(2)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <pgroup.h>
#include <ibtfun.h>
      IF (IPRINT .GE. 10) CALL HEADER('Subroutine DSYM1R',-1)
C
C     Loop over all irreps in molecule
C
      IR(1) = IREPDM
      IFXYZ = IBTXOR(ISYMAX(3,1),ISYMAX(3,2))
      IR(2) = IBTXOR(IREPDM,IFXYZ)
      ISTRA = 1
C      CALL DZERO(AAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         DO 150 II = 1,2
         NBA   = NCOS(IREPA+1,II)
         DO 200 I = ISTRA,ISTRA + NBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
               DO 450 JJ = 1,2
                  ID    = MOD(II+JJ,2)+1
                  NBB   = NCOS(IREPB+1,JJ)
                  IF (IBTXOR(IREPA,IREPB).EQ.IR(ID)) THEN
                  DO 500 J = ISTRB,ISTRB + NBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
Cjth - insert
                     MAB    = IBTOR(MULA,MULB)
                     KAB    = IBTAND(MULA,MULB)
                     HKAB   = FMULT(KAB)
Cjth - endinsert
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        AAO(INDA,INDB) = AAO(INDA,INDB)
     &                                 + HKAB*FACA*FACB*ASO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NBB
  450             CONTINUE
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NBA
  150    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('DSYM1R.INPUT: Total matrix in SO basis',-1)
         CALL OUTPUT(ASO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('DSYM1R.OUTPUT:Total matrix in AO basis',-1)
         CALL OUTPUT(AAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdout */
      SUBROUTINE GRDOUT(IPRINT)
C*****************************************************************************
C
C     (c) 1997/05/28 by J. Thyssen
C
C
C     Last revision: 1997/06/11 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
#include <nuclei.h>
#include <taymol.h>
#include <dummy.h>
C
#include <dcbpsi.h>
#include <dcbgrd.h>
C
      DIMENSION CSTRA(9*MXCOOR*MXCOOR),SCTRA(9*MXCOOR*MXCOOR)
      CHARACTER*(20) LABEL(3)
C
      CALL QENTER('GRDOUT')
C
      CALL TITLER('Final output from PSIGRD','*',110)
C
C     Sum up gradient
C
      DO I = 1,3*NUCDEP
         GRADER(I,0) = GRADER(I,1) + GRADER(I,2) + GRADER(I,3)
         GRADRO(I,0) = GRADRO(I,1) + GRADRO(I,2)
         GRADNU(I,0) = GRADNU(I,1) + GRADNU(I,2)
         GRADKN(I,0) = GRADKN(I,1) + GRADKN(I,2)
      END DO
C
      IPRINT = 2
      IF ( IPRINT .GE. 2 ) THEN
C
C       Kinetic energy
C       ==============
C
        CALL HEADER('Kinetic energy integral gradient',-1)
        LABEL(1) = 'Momentum'
        LABEL(2) = 'Beta mtx.'
        CALL PRIGDR(GRADKN,2,CSTRA,SCTRA,LABEL)
C
C       Nuclear attraction
C       ==================
C
        CALL HEADER('Nuclear attraction integral gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'SS'
        CALL PRIGDR(GRADNU,2,CSTRA,SCTRA,LABEL)
C
C       Reorthonormalization gradient
C       =============================
C
        CALL HEADER('Reorthonormalization gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'SS'
        CALL PRIGDR(GRADRO,2,CSTRA,SCTRA,LABEL)
C
C       Two-electron gradient
C       =====================
C
        CALL HEADER('Two-electron integral gradient',-1)
        LABEL(1) = 'LL'
        LABEL(2) = 'LS'
        LABEL(3) = 'SS'
        CALL PRIGDR(GRADER,3,CSTRA,SCTRA,LABEL)
C
C       "SS-integral contribution by classical repulsion of
C        small component atomic charges"-gradient
C       ===================================================
C
        IF ( DO2GRD(2) ) THEN
           CALL HEADE3('"SS-integral contribution by classical ',
     &                 'repulsion of small component atomic ',
     &                 'charges"-gradient',-1)
           IF ( DO2GRD(3) ) THEN
              WRITE(LUPRI,'(/A)') 'SCC gradient is NOT used, since '//
     &                            'SS-gradient is calculated.'
           ELSE
              WRITE(LUPRI,'(/A)') 'SCC gradient is used, since '//
     &                            'SS-gradient is NOT calculated.'
           END IF 
           CALL PRIGDR(GRADLV,0,CSTRA,SCTRA,LABEL)
        END IF
C       Nuclear repulsion
C       =================
C
        CALL HEADER('Nuclear repulsion gradient',-1)
        CALL PRIGDR(GRADNN,0,CSTRA,SCTRA,LABEL)
C
C
      END IF
C
C
      CALL ZERGRD
      CALL ADDGRD(GRADKN(1,0))
      CALL ADDGRD(GRADNU(1,0))
      CALL ADDGRD(GRADRO(1,0))
      CALL ADDGRD(GRADER(1,0))
      CALL ADDGRD(GRADNN)
      IF ( DO2GRD(2) .AND. .NOT. DO2GRD(3) ) THEN
         CALL ADDGRD(GRADLV)
      END IF
      CALL HEADER('Total molecular gradient',-1)
      CALL PRIGDR(GRDMOL,0,CSTRA,SCTRA,LABEL)
C
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GRDOUT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtcnvd */
      SUBROUTINE GTCNVD(DMAT,CMO,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/06/09 by J. Thyssen
C
C     Get contravariant denstity matrix without folding
C     Copied from GTCNV1
C     
C     Last revision: 1997/06/10 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(*),DMAT(N2BBASX,4,NFMAT),CMO(*)
      LOGICAL FOLD
      PARAMETER ( D1 = 1.0 D00, D2 = 2.0 D00 )
#include <dcbbas.h>
#include <dcbgen.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbdhf.h>
#include <dcbqal.h>
#include <dummy.h>
C
#include <ibtfun.h>
      IQOFF(IZ,IMAT) = (IZ-1) * N2BBASX + (IMAT-1) * N2BBASXQ
C
#include <memint.h>
C
C
      CALL QENTER('GTCNVD')
      IF ( IPRINT .GE. 5) CALL HEADER('Output from GTCNVD',-1)
C 
      CALL MEMGET('REAL',KTMP,N2BBASXQ,WORK,KFREE,LFREE)
C
C
C     Read coefficient matrix and eigenvalues
C     ===================================================
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','DHFGRD')
      CALL REACMO(LUCOEF,'DFCOEF',CMO,DUMMY,TOTERG,3)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Contruct density matrix
C
      DO IMAT = 1,NFMAT
         IBIT = IBTSHL(1,IMAT-1)
         CALL GENDEN(WORK(KTMP),CMO,IBIT,IPRINT)
C
C        Scale inactive density matrix with D2 (= occupation)
C
         IF (IMAT.EQ.1)
     &      CALL DSCAL(N2BBASXQ,D2,WORK(KTMP),1)
C
C        Transform to unsorted basis
C
         IREP = ISYMOP(IMAT)-1
         CALL BSTOBU(WORK(KTMP),NZ,WORK(KFREE),LFREE)
         IF(IPRINT.GE.5) THEN
            CALL HEADER('GTCNVD: DMAT After BSTOBU to unsorted basis',
     &                  -1)
            WRITE(LUPRI,'(A,I5)') '*** Matrix no. ',IMAT
            CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         END IF
C
C
C        Insert half-phases
C        ==================
C
         IREP = ISYMOP(IMAT)-1
         IF(NZ.LT.4) THEN
            DO IZ = 1,NZ
               IQ = IPQTOQ(IZ,IREP)
               CALL Q2PHASE('D',IQ,1,WORK(KTMP))
            ENDDO
            IF(IPRINT.GE.5) THEN
               CALL HEADER('Unsorted DMAT with half-phases inserted',-1)
               CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
            ENDIF
         ENDIF
C
C
C        Transform density matrix from QO to AO
C        ======================================
C
         CALL DZERO(DMAT(1,1,IMAT),4*N2BBASX)
         IREP = ISYMOP(IMAT)-1
         IPAR = JBTOF(IREP,1)
         DO IZ = 1,4
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
            IPQ   = IQTOPQ(IQ,IREP)
            IF(IPQ.GT.NZ) CALL QUIT('GTCNV1.DENSITY: IPQ.GT.NZ !')
            IF(IPQ.LE.0 ) CALL QUIT('GTCNV1.DENSITY: IPQ.LE.0 !')
            IADR = (IPQ-1) * N2BBASX
            CALL DTSOAO(WORK(KTMP+IADR),DMAT(1,IZ,IMAT),
     &                 NTBAS(0),IREPD,IPRINT)
         ENDDO
      ENDDO
      IF (NASHT.GT.0) THEN
         CALL DAXPY(N2BBASX*4,D1,DMAT(1,1,2),1,DMAT(1,1,1),1)
      END IF
      IF ( IPRINT .GE. 10 ) THEN
         DO IMAT = 1,NFMAT
C
C        Final print section
C
            CALL HEADER('GTCNVD: Final matrices.',-1)
            CALL HEADER('Density matrix:',-1)
            WRITE(LUPRI,*) 'Matrix number ',IMAT
            CALL PRQMAT(DMAT(1,1,IMAT),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),4,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C     Dispose used memory
C
      CALL MEMREL('GTCNVD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNVD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck estgls */
      SUBROUTINE ESTGLS(DOLS)
C*****************************************************************************
C
C     (c) 1997/06/27 by J. Thyssen
C
C     On basis of one-electron gradient,g(1), and two-electron LL-gradient, 
C     g(2,LL), estimate if two-electron LS- gradient is necessary. 
C
C     If ||g(2,LS)|| > 0.1 ||g(1)+g(2,LL)|| then we calculate g(2,LS).
C     As an estimate of ||g(2,LS)|| we use 0.005 ||g(2,LL)|| based
C     on empirical facts.
C
C     Since the LS-gradient is based on (LL|SS) or (SS|LL) integrals
C     it scales with 1/(c^2). 
C
C     Last revision: 1997/07/02 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
      PARAMETER( D0 = 0.0 D00 , DM2 = -2.0 D00, CSTD = 137.0359895)
      DIMENSION TGRAD(MXCOOR)
      LOGICAL DOLS
      CHARACTER*3 STR
#include <nuclei.h>
#include <symmet.h>
#include <dcbgen.h>
#include <dcbgrd.h>
C
      CALL QENTER('ESTGLS')
C
      DOLS = .FALSE.
C
      DO I = 1, 3*NUCDEP
C
C        Calculate "total" gradient
C
         TGRAD(I)  = GRADNN(I)   + 
     &               GRADRO(I,1) + GRADRO(I,2) +
     &               GRADKN(I,1) + GRADKN(I,2) +
     &               GRADNU(I,1) + GRADNU(I,2) +
     &               GRADER(I,1)
C
      END DO
C
C     Calculate norms
C
      TOTNRM = DDOT(3*NUCDEP,TGRAD,      1,TGRAD,      1)
      TOTNRM = SQRT(TOTNRM)
      DLLNRM = DDOT(3*NUCDEP,GRADER(1,1),1,GRADER(1,1),1)
      DLLNRM = SQRT(DLLNRM)
C
C     If ||g(2,LL)|| / ||g(1)+g(2,LL)|| > 20 then calculate LS-gradient
C
      CFAC = (( CSTD / CVAL) ** ( DM2 ) )
      IF ( TOTNRM .GT. D0 ) THEN
         DL = DLLNRM / TOTNRM
         DR = 2.0 D01 * CFAC
         IF ( DL .GE. DR ) THEN
            DOLS = .TRUE.
         END IF
      ELSE
         DOLS = .TRUE.
      END IF
      WRITE(LUPRI,'(/1X,A)') 'Output from ESTGLS'
      WRITE(LUPRI,1000) 'Norm of g(1)+g(2,LL): ',TOTNRM
      WRITE(LUPRI,1000) 'Norm of g(2,LL)     : ',DLLNRM
      WRITE(LUPRI,1000) 'Est. norm of g(2,LS): ',DR
      CALL CYESNO(STR,DOLS)
      WRITE(LUPRI,'(1X,A,A3/)') 'Calculate g(2,LS)   : ',STR
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('ESTGLS')
      RETURN
 1000 FORMAT(1X,A,F16.9)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck estgss */
      SUBROUTINE ESTGSS(DOSS)
C*****************************************************************************
C
C     (c) 1997/07/02 by J. Thyssen
C
C     Do something similar as in ESTGLS
C
C     Last revision: 1997/07/02 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
      PARAMETER( D0 = 0.0 D00 , DM4 = -4.0 D00, CSTD = 137.0359895, 
     &           D2 = 2.0 D00 , DP1 = 1.0D-01)
      DIMENSION TGRAD(MXCOOR)
      LOGICAL DOSS
      CHARACTER*3 STR
#include <nuclei.h>
#include <symmet.h>
#include <dcbgen.h>
#include <dcbgrd.h>
#include <optinf.h>
C
      CALL QENTER('ESTGSS')
C
      DOSS = .FALSE.
C
      DO I = 1, 3*NUCDEP
C
C        Calculate "total" gradient
C
         TGRAD(I)  = GRADNN(I)   + 
     &               GRADRO(I,1) + GRADRO(I,2) +
     &               GRADKN(I,1) + GRADKN(I,2) +
     &               GRADNU(I,1) + GRADNU(I,2) +
     &               GRADER(I,1) + GRADER(I,2)
C
      END DO
C
C     Calculate norms
C
      TOTNRM = DDOT(3*NUCDEP,TGRAD,      1,TGRAD,      1)
      TOTNRM = SQRT(TOTNRM)
      DLLNRM = DDOT(3*NUCDEP,GRADER(1,1),1,GRADER(1,1),1)
      DLLNRM = SQRT(DLLNRM)
      DLSNRM = DDOT(3*NUCDEP,GRADER(1,2),1,GRADER(1,2),1)
      DLSNRM = SQRT(DLSNRM)
C
C     As estimate for ||g(2,SS)|| use ||g(2,LS)||^2 / ||g(2,LL)||
C
      CFAC = (( CSTD / CVAL) ** ( DM4 ) )
      IF ( TOTNRM .GT. D0 ) THEN
         DR = ( DLSNRM ** D2 ) / DLLNRM
         DL = CFAC * DP1 * MAX( TOTNRM, GRDTHR )
         IF ( DR .GE. DL ) THEN
            DOSS = .TRUE.
         END IF
      ELSE
         DOSS = .TRUE.
      END IF
C
      WRITE(LUPRI,'(/1X,A)') 'Output from ESTGSS'
      WRITE(LUPRI,1000) 'Norm of g(1)+g(2,LL)+g(2,LS): ',TOTNRM
      WRITE(LUPRI,1000) 'Norm of g(2,LL)             : ',DLLNRM
      WRITE(LUPRI,1000) 'Norm of g(2,LS)             : ',DLSNRM
      WRITE(LUPRI,1000) 'Estimated norm of g(2,SS)   : ',DR    
      CALL CYESNO(STR,DOSS)
      WRITE(LUPRI,'(1X,A,A3/)') 'Calculate g(2,SS)           : ',
     &                              STR
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('ESTGSS')
      RETURN
 1000 FORMAT(1X,A,F16.9)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cyesno */
      SUBROUTINE CYESNO(STR,BOOL)
C*****************************************************************************
C
C     (c) 1997/07/08 by J. Thyssen
C
C     Return yes or no depending on BOOL
C
C     Last revision: 1997/07/08 jth
C
C*****************************************************************************
#include <implicit.h>
      CHARACTER*(*) STR
      LOGICAL BOOL
C
      IF ( BOOL ) THEN
         STR = 'yes'
      ELSE
         STR = 'no'
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck exedir */
      SUBROUTINE EXEDIR(WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/07/09 by J. Thyssen
C
C     EXEcute DIRac; Calculate wavefunction
C
C     Last revision: 1997/07/08 jth
C
C*****************************************************************************
#include <implicit.h>
      DIMENSION WORK(LWORK)
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C 
#include <priunit.h>
#include <dcbgen.h>
#include <dcbdhf.h>
C
      CALL QENTER('EXEDIR')
C
C     Delete DFCYCL so we don't restart DIIS
C     
      CALL DELFILE(LUCYCL,'DFCYCL')
C
C     Delete AOPROPER so we don't use in old integrals
C
#if defined (SYS_IBM)
      CALL DELFILE(LU1INT,'APROPER')
#else
      CALL DELFILE(LU1INT,'AOPROPER')
#endif
C
C     To avoid readin twice in the 0. iteration skip
C     the first time
C
C     Input and initialization
C
      IF ( FIRST ) THEN
         FIRST = .FALSE.
      ELSE
         CALL PAMINP(WORK,LWORK)
         CALL PAMSET(WORK,LWORK)
      END IF
C     
C     Get wavefunction
C
      CALL PAMPSI(WORK,LWORK)
C
C     Calculate properties
C
      IF (DOPRP) CALL PAMPRP(WORK,LWORK)
C
C     Analyze wavefunction
C
      IF(DOANA) CALL PAMANA(WORK,LWORK)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('EXEDIR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck delfile */
      SUBROUTINE DELFILE(IUNIT,FILENAME)
C*****************************************************************************
C
C     (c) 1997/07/09 by J. Thyssen
C
C     Delete file FILENAME
C
C     Last revision: 1997/07/09 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*(*) FILENAME
      LOGICAL EX
C
      CALL QENTER('DELFILE')
C
      INQUIRE(FILE=FILENAME,EXIST=EX)
      IF ( EX ) THEN
         CALL OPNFIL(IUNIT,FILENAME,'UNKNOWN','DELFILE')
         CLOSE(IUNIT,STATUS='DELETE')
         INQUIRE(FILE=FILENAME,EXIST=EX)
         IF ( EX ) WRITE(LUPRI,'(//1X,3A//)') 'WARNING: ',
     &      FILENAME,' could .NOT. be deleted!' 
      END IF
C
C
      CALL QEXIT('DELFILE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcgrd */
      SUBROUTINE LVCGRD(IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/09/23 by J. Thyssen
C
C     "SS-integral contribution by classical repulsion of
C      small component atomic charges"-gradient
C
C     Last revision: 1997/09/23 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
C
#include <nuclei.h>
#include <frame.h>
      DIMENSION WORK(LWORK)
C
#include <memint.h>
C
      CALL QENTER('LVCGRD')
C
      CALL MEMGET('REAL',KSQC,NUCIND,WORK,KFREE,LFREE)
      CALL LVCGR1(WORK(KSQC),IPRINT)
      CALL MEMREL('LVCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('LVCGRD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcgr1 */
      SUBROUTINE LVCGR1(SQC,IPRINT)
C*****************************************************************************
C
C     (c) 1997/09/23 by J. Thyssen
C
C     "SS-integral contribution by classical repulsion of
C      small component atomic charges"-gradient
C
C     Last revision: 1997/09/23 jth
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0,
     *           THREE = 3.0D0)
#include <nuclei.h>
#include <frame.h>
#include <dorps.h>
#include <symmet.h>
C
      DIMENSION SQC(NUCIND)
#include <dcbgrd.h>
#include <ibtfun.h>
C
C
      CALL QENTER('LVCGR1')
C
C     Get small component charges
C
      CALL FLSHFO(LUPRI)
      CALL SCDENS(CHARGE,SQC,NUCIND,IPRINT)
C
C     This is copied from ABACUS/NUCREP
C
      NCDEP3 = 3*NUCDEP
      CALL DZERO(GRADLV,NCDEP3)
C
C
C     Run over symmetry-independent nuclei A
C
      DO 200 NCENTA = 1, NUCIND
         CHARGA = SQC(NCENTA)
         IF (ABS(CHARGA) .GT. ZERO) THEN
            NAX    = 3*NCENTA - 2
            NAY    = 3*NCENTA - 1
            NAZ    = 3*NCENTA
            CORDXA = CORD(1,NCENTA)
            CORDYA = CORD(2,NCENTA)
            CORDZA = CORD(3,NCENTA)
            MULA   = ISTBNU(NCENTA)
C
C           Run over symmetry-independent nuclei B
C
            DO 300 NCENTB =  1, NCENTA
               CHARGB = SQC(NCENTB)
               IF (ABS(CHARGB) .GT. ZERO) THEN
                  NBX    = 3*NCENTB - 2
                  NBY    = 3*NCENTB - 1
                  NBZ    = 3*NCENTB
                  CORBX0 = CORD(1,NCENTB)
                  CORBY0 = CORD(2,NCENTB)
                  CORBZ0 = CORD(3,NCENTB)
                  MULB   = ISTBNU(NCENTB)
C
                  MAB    = IBTOR (MULA,MULB)
                  KAB    = IBTAND(MULA,MULB)
                  HKAB   = FMULT(KAB)
                  CROSS = ONE
                  IF (NCENTA .EQ. NCENTB) THEN
                     HKAB = HALF*HKAB
                     CROSS = TWO
                  END IF
C
C                 Run over symmetry independent charge distributions
C
                  DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .EQ. 0) THEN
                     ICENTA = NUCNUM(NCENTA,1)
                     ICENTB = NUCNUM(NCENTB,ISYMOP+1)
                  IF (ICENTA .EQ. ICENTB) GO TO 400
C
                     SBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                     SBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                     SBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                     XAB = CORDXA - SBX*CORBX0
                     YAB = CORDYA - SBY*CORBY0
                     ZAB = CORDZA - SBZ*CORBZ0
C
                     XAB2   = XAB*XAB
                     YAB2   = YAB*YAB
                     ZAB2   = ZAB*ZAB
                     RAB2   = XAB2 + YAB2 + ZAB2
                     RAB1   = SQRT(RAB2)
                     ZZR3IN = - HKAB*CHARGA*CHARGB/(RAB1*RAB2)
                     VNUCX  = XAB*ZZR3IN
                     VNUCY  = YAB*ZZR3IN
                     VNUCZ  = ZAB*ZZR3IN
C
C                    ********************
C                    ***** Gradient *****
C                    ********************
C
C                    Totally symmetric contributions only
C
                     IF (DOREPS(0)) THEN
                        IAX  = IPTCNT(NAX,0,1)
                        IAY  = IPTCNT(NAY,0,1)
                        IAZ  = IPTCNT(NAZ,0,1)
                        IBX  = IPTCNT(NBX,0,1)
                        IBY  = IPTCNT(NBY,0,1)
                        IBZ  = IPTCNT(NBZ,0,1)
                        IF (IAX.NE.0) GRADLV(IAX) =GRADLV(IAX)+VNUCX
                        IF (IAY.NE.0) GRADLV(IAY) =GRADLV(IAY)+VNUCY
                        IF (IAZ.NE.0) GRADLV(IAZ) =GRADLV(IAZ)+VNUCZ
                        IF (IBX.NE.0) GRADLV(IBX) =GRADLV(IBX)-SBX*VNUCX
                        IF (IBY.NE.0) GRADLV(IBY) =GRADLV(IBY)-SBY*VNUCY
                        IF (IBZ.NE.0) GRADLV(IBZ) =GRADLV(IBZ)-SBZ*VNUCZ
                     END IF
                  END IF
  400             CONTINUE
               END IF
  300       CONTINUE
         END IF
  200 CONTINUE
C
C
      CALL QEXIT('LVCGR1')
C
      RETURN
      END
C
C
      SUBROUTINE HEADE3(HEAD1,HEAD2,HEAD3,IN)
#include <implicit.h>
      CHARACTER*(*) HEAD1,HEAD2,HEAD3
#include <priunit.h>
C
      IF (IN .GE. 0) THEN
         I1 = IN + 1
         I2 = I1
         I3 = I1
      ELSE
         L1 = LEN(HEAD1)
         L2 = LEN(HEAD2)
         L3 = LEN(HEAD3)
         LM = MAX(L1,L2,L3)
         I1 = (72 - L1) /2 + 1
         I2 = (72 - L2) /2 + 1
         I3 = (72 - L3) /2 + 1
         IM = (72 - LM) /2 + 1
      END IF
      WRITE(LUPRI, '(//,80A)') (' ',I=1,I1), HEAD1
      WRITE(LUPRI, '(80A)')    (' ',I=1,I2), HEAD2
      WRITE(LUPRI, '(80A)')    (' ',I=1,I3), HEAD3
      WRITE(LUPRI, '(80A)') (' ',I=1,IM), ('-',I=1,LM)
      WRITE(LUPRI, '()')
      RETURN
      END
C  /* Deck prigdr */
      SUBROUTINE PRIGDR(GRAD,ICOMPO,CSTRA,SCTRA,LABEL)
C
C     Written by Joern Thyssen 1997/06/20
C     Last revision: 1997/06/20 jth
C
C     Print relativistic gradient
C     Based on PRIGRD
C     GRAD is asummed to be an array GRAD(MXCOOR,ICOMPO), where
C     (1) GRAD(0) is total gradient
C     (2) GRAD(1) is LABEL(1) component of gradient
C     (3) GRAD(2) is LABEL(2) c... etc.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
C
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION GRAD(MXCOOR,0:ICOMPO), CGRAD(MXCOOR,0:3), 
     &          CSTRA(*), SCTRA(*)
      CHARACTER*(*) LABEL(ICOMPO)
      CHARACTER*10  TLABEL
C
      DIMENSION DNORM(0:3)
C
C
      IF (.NOT.DOSYM(1)) THEN
         WRITE (LUPRI,'(2X,A)')
     &  ' Gradient not calculated - '//
     &  ' totally symmetric distortions not requested.'
      ELSE IF (MAXREP .EQ. 0) THEN
         IOFF = 0
         DO 100 IATOM = 1, NUCDEP
            TLABEL = 'Total     '
            IF ( ICOMPO .EQ. 0 ) TLABEL = '          '
            WRITE (LUPRI,1000) NAMDEP(IATOM),TLABEL,
     &                         (GRAD(IOFF+J,0),J=1,3)
            DO K = 1, ICOMPO
               WRITE (LUPRI,1000) '      ',LABEL(K),
     &                            (GRAD(IOFF+J,K),J=1,3)
            END DO
            IF ( ICOMPO .GT. 0 ) WRITE(LUPRI,'()')
            IOFF = IOFF + 3
  100    CONTINUE
      ELSE
         TLABEL = 'Total     '
         IF ( ICOMPO .EQ. 0 ) TLABEL = '          '
         DO 200 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(16X,A6,3X,A10,F17.10)') 
     &        NAMEX(IPTCOR(I,1)),TLABEL,GRAD(I,0)
            DO J = 1, ICOMPO
               WRITE (LUPRI,'(16X,A6,3X,A10,F17.10)') 
     &           '      ',LABEL(J),GRAD(I,J)
            END DO
  200    CONTINUE
         WRITE (LUPRI,'(//)')
         DO K = 0, ICOMPO
            CALL TRAGRD(GRAD(1,K),CGRAD(1,K),CSTRA,SCTRA,NCRREP(0,1),
     &                  3*NUCDEP)
         END DO
         IOFF = 0
         DO 300 IATOM = 1, NUCDEP
            WRITE (LUPRI,1000) NAMDEP(IATOM),TLABEL,
     &                         (CGRAD(IOFF+J,0),J=1,3)
            DO K = 1, ICOMPO
               WRITE (LUPRI,1000) '      ',LABEL(K),
     &                            (CGRAD(IOFF+J,K),J=1,3)
            END DO
            IF ( ICOMPO .GT. 0 ) WRITE(LUPRI,'()')
            IOFF = IOFF + 3
  300    CONTINUE
C
      END IF
C 
C     Calculate weight statistics
C
      IF ( ICOMPO .GT. 0 ) THEN
         DNORM(0) = 0.0D00
         DO K = 1,ICOMPO
            DNORM(K) = DDOT(3*NUCDEP,GRAD(1,K),1,GRAD(1,K),1)
            DNORM(K) = SQRT(DNORM(K))
            DNORM(0) = DNORM(0) + DNORM(K)
         END DO
         WRITE(LUPRI,'(A)') 'Weights:'
         DO K = 1,ICOMPO
            WRITE(LUPRI,'(1X,A10,5X,F13.9,A)') LABEL(K),
     &            1.0D02 * DNORM(K) / DNORM(0),'%'
         END DO
      END IF
      WRITE (LUPRI,'(//)')
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT (1X,A6,1X,A10,F17.10,2F21.10)
      END
