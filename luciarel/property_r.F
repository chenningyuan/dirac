!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETINT_NTS(YINT,INDEX,IRI,IDIM,DMFC,IPRNT)
C***********************************************************************
C
C     fetch integral to YINT for non-totally symmetric operators 
C     from Molfdir-type array DMFC.  
C
C     order in DIRAC (Molfdir format) matrix DMFC:
C                                 _
C           p spinors followed by p
C                                                 _     _
C           if NFSYM = 2 (inversion symmetry): p, p, q, q
C
C    *********************************************************
C    * NOTE: the ordering is GA space independent, e.g.      *
C    *                                                       *
C    * .GASSH   yields the very same DMFC matrix as   .GASSH *
C    * 3                                              2      *
C    * 1 1                                            1 1    *
C    * 0 2                                            1 2    *
C    * 1 0                                                   *
C    *********************************************************
C
C     The present routine exhibits a new and more general offset 
C     calculation for accessing the integral matrix DMFC compared to 
C     GETIND_NZ3. The old offset calculation in GETIND_NZ3 is *NOT* valid 
C     for general one-elctron operators.
C     The new one looks a little more complicated though...
C
C     Written by S. Knecht - Dec 2009
C
C     Last revision :
C
C***********************************************************************
!     use memory_allocator
#include "implicit.h"
#include "dgroup.h"
#include "dcborb.h"
!
      DIMENSION DMFC(*), ICSAVE(4), INDEX(4)
      integer jc_fsym, ja_fsym, nasht2, n2ashx2
!
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRNT)
!
!     set dimensions for DMFC matrix. Note that we here access the
!                         _   _
!     matrix in spinor (p,p,q,q) and NOT quaternion-packed Kramers-paired 
!     basis.
!     E.g. IASH(2) refers in the usual Dirac notation to the 
!     offset in fermion symmetry 2 wrt to the Kramers-paired basis. 
!     This means that we must multiply the offset by 2 to get the 
!     correct spinor offset.
!
      nasht2  = NASHT  * 2
      n2ashx2 = nasht2 * nasht2
!
!     input print
c     IF(NTEST .ge. 5)THEN
c       WRITE(6,*) ' GETIND_NTS : '
c       WRITE(6,*) '  No. of active orbitals : ',NASHT
c       WRITE(6,*) '  Compound index dimension: ',NNASHX
c       WRITE(6,'(A,6X,4I3)') '  INDEX array: ',(INDEX(I),I=1,4,1)
c       WRITE(6,*) '1-el. integral array in GETIND_NTS:'
c       CALL WRTMAT(DMFC,IDIM*2,IDIM*2,IDIM*2,IDIM*2)
c     END IF
!
      NCT = 0
      DO IC=1,4,1
        IF( INDEX(IC).ne. 0)THEN
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        END IF
      END DO
!
!     NCT = 2: pick out one-electron integral identified by 
!     two indices IORI (creator index) and JORI (annihilator index)
      IORI = INDEX(ICSAVE(1))
      JORI = INDEX(ICSAVE(2))
!
!     get fermion symmetry of creator and annihilator 
      jc_fsym = 1
      ja_fsym = 1
      IF( INDEX(ICSAVE(1)) .gt. NASH(1)) jc_fsym = 2
      IF( INDEX(ICSAVE(2)) .gt. NASH(1)) ja_fsym = 2
!
!     compute input address on array DMFC:
!     compute a relative creator and annihilator spinor 
!     (not Kramers pair) index
!
!     1. creator index
      I_REL = IORI + ((ICSAVE(1)-1)*NASH(jc_fsym))+IASH(jc_fsym)
!
!     2. annihilator index
      J_REL = JORI + ((ICSAVE(2)-3)*NASH(ja_fsym))+IASH(ja_fsym)
!
!     combine I_REL and J_REL to compute the address IIN
!     --------------------------------------------------
!     first  line: offset wrt real/imaginary
!     second line: offset wrt annihilator index
!     third  line: offset wrt the creator index
      IIN =  ((IRI-1)  * n2ashx2) +
     &       ((J_REL-1)* nasht2)  +
     &         I_REL
!
!     print section
c     if (NTEST.ge.1) then
c       write(6,*) 'orig. index pair ',IORI,JORI,' address :',IIN
c       write(6,'(A,I5)') ' I_REL = ',I_REL
c       write(6,'(A,I5)') ' J_REL = ',J_REL
c       DO LLLL = 1, (n2ashx2 * NZ)
c         WRITE(6,*) ' element ',LLLL,'in matrix DMFC is = ',DMFC(LLLL)
c       END DO
c     end if
!
!     fetch integral
      YINT = DMFC(IIN)
!
!     sign for imaginary 1-electron integrals:
!     (according to h_(IJ) = -h_(I_J_), reference are the 
!     unbarred spinors.)
!
!     stefan: the comment above does not apply anymore to the new 
!     sorting scheme. here, we pick out the proper integral h_(I_J_) 
!     by default.
!     if (IRI.eq.2.and.INDEX(1).eq.0) YINT = - YINT

!     write(6,*) '1-electron integral.'
!     write(6,'(3I4,2X,1E14.8)') I_REL,J_REL,IIN,YINT
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_IAIB_DIRAC(IASP_KRCI_SYM,IBSP_KRCI_SYM,ELM_ABSMAX,
     &                          AMAT,NROW,NCOL,LRQ,LCQ,JZ,IQP,NTEST)
C***********************************************************************
C
C     get symmetry indices of spinors which are connected by the 
C     largest absolute element ELM_ABSMAX in matrix AMAT(LROW,LCOL,JZ).
C
C     order in DIRAC (Molfdir format) matrix:
C                                 _
C           p spinors followed by p
C                                                 _     _
C           if NFSYM = 2 (inversion symmetry): p, p, q, q
C
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMAT(NROW,NCOL,JZ), IQP(*)
      CHARACTER M(4)*13
      LOGICAL INDICES_FOUND
      INDICES_FOUND = .FALSE.
      M(1) = 'Matrix A real'
      M(2) = 'Matrix A imag'
      M(3) = 'Matrix B real'
      M(4) = 'Matrix B imag'
C
      DO 10 IZ = 1, JZ
        IQ = IQP(IZ)
        IF( NTEST .gt. 5 )THEN
           WRITE(LUPRI,'(/A,A13,A,2I4/)')
     &     '*** ',M(IQ),' part ***;    NROW,NCOL:',NROW,NCOL
        END IF
        CALL GET_KRCI_SPINOR_IND(AMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,1,
     &                           ELM_ABSMAX,IASP_KRCI_SYM,IBSP_KRCI_SYM,
     &                           INDICES_FOUND,NTEST)
 10   CONTINUE
C
C     error exit if no indices pair is found
C
      IF( .NOT. INDICES_FOUND ) THEN
        WRITE(LUPRI,'(/A)')'  *** ERROR in GET_IAIB_DIRAC ***'//
     &                     ' failed to find indices pair'//
     &                     ' of DIRAC spinors for largest element in'//
     &                     ' property matrix.'
        CALL QUIT('  *** ERROR in GET_IAIB_DIRAC *** failed to find 
     &             DIRAC spinor indices pair')

      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_KRCI_SPINOR_IND(AMATRX,ROWLOW,ROWHI,COLLOW,
     &                               COLHI,ROWDIM,COLDIM,NCTL,
     &                               ELM_ABSMAX,
     &                               IASP_KRCI_SYM,IBSP_KRCI_SYM,
     &                               INDICES_FOUND,NTEST)
C***********************************************************************
C
C     find spinor indices IASP_DC and IBSP_DC and determine the 
C     corresponding symmetry indices IASP_KRCI_SYM and IBSP_KRCI_SYM.
C
C     order in DIRAC (Molfdir format) matrix:
C                                 _
C           p spinors followed by p
C                                                 _     _
C           if NFSYM = 2 (inversion symmetry): p, p, q, q
C
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxpdim.inc"
#include "lucinp.inc"
#include "symm.inc"
#include "cgas.inc"
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      LOGICAL INDICES_FOUND
C
      CALL QENTER('GET_KRCI_SPINOR_IND')
C
C     ... check if we already have what we want...
C
      IF (INDICES_FOUND) GOTO 999
C
      IASP_KRCI_SYM = 0
      IBSP_KRCI_SYM = 0
      IASP_DC       = 0
      IBSP_DC       = 0
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       CONTINUE
            DO 7 I = BEGIN,LAST 
              IF( INDICES_FOUND ) GO TO 7
              IF( DABS(AMATRX(K,I)) .eq. ELM_ABSMAX )THEN
                 IASP_DC = K
                 IBSP_DC = I
                 INDICES_FOUND = .TRUE.
C                WRITE(LUPRI,*) ' indices K, I',K,I
C                WRITE(LUPRI,*) ' AMATRX(K,I)',AMATRX(K,I)
              END IF
    7       CONTINUE
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
C
C     row and column index determined - translate to LUICAREL symmetry
C     indices
C
      IF( .NOT. INDICES_FOUND ) GOTO 999
C
C     test print:    
C
      IF( NTEST .gt. 0 ) THEN
        WRITE(LUPRI,'(/A,I5,A,I5)') '  SPINOR indices: index a',
     &                                 IASP_DC,' and SPINOR '//
     &                                 'index b',IBSP_DC
      END IF
C
C     loop over all spinors (GAS1, GAS2, ...) in a given irrep 
C     until we find the correct index:
C
C     hardwired order in C2h:  Ag  Bg ||  Au  Bu
C
C     gerade   unbarred spinors: Ag
C     gerade     barred spinors: Bg
C
C     ungerade unbarred spinors: Au
C     ungerade   barred spinors: Bu
C
C     --> therefore the following loop structure is chosen:
C
C     loop over irreps
C         loop over GA spaces for unbarred spinors
C         end loop
C         loop over GA spaces for barred spinors
C         end loop
C     end loop
C
C     get "A spinor" symmetry for spinor index IASP_DC
C
      ILOOP_INDEX = 0
      IXSYM_REP   = 0
      INUMORB_GAS = 0
 10   CONTINUE
      IXSYM_REP = IXSYM_REP + 1
C
C     ... unbarred spinors p
      IGAS = 0
 20   CONTINUE
      IGAS = IGAS + 1
      INUMORB_GAS = NGSOB(IXSYM_REP,IGAS)
      IF( INUMORB_GAS .gt. 0) THEN
         ISNUMORBL = 0
 22      CONTINUE
         ILOOP_INDEX = ILOOP_INDEX + 1
         ISNUMORBL   = ISNUMORBL + 1
C        ... index matches
         IF( ILOOP_INDEX .eq. IASP_DC ) GOTO 30
C        ... next spinor in this GA space
         IF ( ISNUMORBL .lt. INUMORB_GAS ) GOTO 22
      END IF
C     ... next GA space     
      IF( IGAS .lt. NGAS ) GOTO 20
C                        _
C     ... barred spinors p
      IGAS = 0
 25   CONTINUE
      IGAS = IGAS + 1
      INUMORB_GAS = NGSOB2(IXSYM_REP,IGAS)
      IF( INUMORB_GAS .gt. 0) THEN
         ISNUMORBL = 0
 27      CONTINUE
         ILOOP_INDEX = ILOOP_INDEX + 1
         ISNUMORBL   = ISNUMORBL + 1
C        ... index matches
         IF( ILOOP_INDEX .eq. IASP_DC ) GOTO 30
C        ... next spinor in this GA space
         IF ( ISNUMORBL .lt. INUMORB_GAS ) GOTO 27
      END IF
C     ... next GA space     
      IF( IGAS .lt. NGAS ) GOTO 25
C     ... goto next irrep - start at GAS 1 for unbarred p
      IF( IXSYM_REP .lt. NSMOB ) GOTO 10
C
 30   CONTINUE
C
C     assign corresponding irrep.
C     =========================
      IASP_KRCI_SYM = IXSYM_REP
C     =========================
C
C     get "B spinor" symmetry for spinor index IBSP_DC
C
      ILOOP_INDEX = 0
      IXSYM_REP   = 0
C
 40   CONTINUE
      IXSYM_REP = IXSYM_REP + 1
C
C     ... unbarred spinors p
      IGAS = 0
 50   CONTINUE
      IGAS = IGAS + 1
      INUMORB_GAS = NGSOB(IXSYM_REP,IGAS)
      IF( INUMORB_GAS .gt. 0) THEN
         ISNUMORBL = 0
 52      CONTINUE
         ILOOP_INDEX = ILOOP_INDEX + 1
         ISNUMORBL   = ISNUMORBL + 1
C        ... index matches
         IF( ILOOP_INDEX .eq. IBSP_DC ) GOTO 60
C        ... next spinor in this GA space
         IF ( ISNUMORBL .lt. INUMORB_GAS ) GOTO 52
      END IF
C     ... next GA space     
      IF( IGAS .lt. NGAS ) GOTO 50
C                        _
C     ... barred spinors p
      IGAS = 0
 55   CONTINUE
      IGAS = IGAS + 1
      INUMORB_GAS = NGSOB2(IXSYM_REP,IGAS)
      IF( INUMORB_GAS .gt. 0) THEN
         ISNUMORBL = 0
 57      CONTINUE
         ILOOP_INDEX = ILOOP_INDEX + 1
         ISNUMORBL   = ISNUMORBL + 1
C        ... index matches
         IF( ILOOP_INDEX .eq. IBSP_DC ) GOTO 60
C        ... next spinor in this GA space
         IF ( ISNUMORBL .lt. INUMORB_GAS ) GOTO 57
      END IF
C     ... next GA space     
      IF( IGAS .lt. NGAS ) GOTO 55
C     ... goto next irrep - start at GAS 1 for unbarred p
      IF( IXSYM_REP .lt. NSMOB ) GOTO 40
C
 60   CONTINUE
C     assign corresponding irrep.
C     =========================
      IBSP_KRCI_SYM = IXSYM_REP
C     =========================
C
      IF( NTEST .gt. 0 )THEN
        WRITE(LUPRI,'(/A,/A)') 
     &  '  GET_KRCI_SPINOR_IND:        IASP_DC          IBSP_DC   ',
     &  '  (indices pair found)     _____________    _____________'
C                                                  
        WRITE(LUPRI,'(I35,I17)') IASP_DC,IBSP_DC
        WRITE(LUPRI,'(/A,/A)') 
     &  '                           IASP_KRCI_SYM    IBSP_KRCI_SYM',
     &  '  (irreps assigned)        _____________    _____________'
C                                                  
        WRITE(LUPRI,'(I35,I17)') IASP_KRCI_SYM,IBSP_KRCI_SYM
      END IF
 999  CALL QEXIT('GET_KRCI_SPINOR_IND')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_XPROP_REOD_KRCI(LUF1MAT,WRK,LWRK,CIRUN)
C***********************************************************************
C
C     property operator matrix in DIRAC (Molfdir format) order. 
C     Reorder to KR-CI format according to SIGDEN formalism and store 
C     on file KRMC_FOCK (file handle LUF1MAT).
C                                  
C     Written by S. Knecht - Sep 2008
C
C     Last revision :
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
* boson symmetry info
#include "maxorb.h"
#include "../krmc/dcbbos.h"
#include "dcborb.h"
#include "dgroup.h"
* memory
#include "wrkspc.inc"
#include "ipoist8.inc"
C print flags
#include "cprnt.inc"
C required common blocks
#include "mxpdim.inc"
#include "lucinp.inc"
#include "crun.inc"
#include "cintfo.inc"
#include "picasso.inc"
#include "ctcc.inc"
#include "glbbas.inc"
C information for property run
#include "krciprop.h"
#include "dcbxpr.h"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
      CHARACTER CIRUN*6, FILELAB*8, FILELAB_REOD*8, RODLAB*2
      LOGICAL MCSCF, RUN_SHMEM_T
      DIMENSION WRK(LWRK)
      real*8, allocatable :: XXXM(:)
C
      CALL QENTER('GET_XPROP_REOD_KRCI')
C
      RODLAB = 'RO'
      MZ = MIN(NZ,2)
      IMEMCK = 0
      MCSCF = .true.
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C     output file handle in KR-CI module
      LUWRT = LUPRI
C
C     Compute offset from WRK to WORK
      CALL COMPOFF(WRK,WORK,K_OFFSET)
C     *k_offset = (long) (wrk-work+1);
      MXPWRD = LWRK
      KADD = MXPWRD
C
C     initialize memory
      call memmar(K_OFFSET,KADD,'INI   ',IDUMMY,'DUMMY')
C
C     Internal string information (stored in WORK, bases in /STRBAS/)
      CALL STRINF_GAS_REL(WORK(1),IPRSTR)
C
C     Internal subspaces
      CALL LCISPC_REL(LCSBLK,IPRCIX)
C
C     IDGSTYP determines the integral classes and corresponding operator
C     classes depending on the double group in use.
C
      IDGSTYP = -1
C
      IF ( PNTGRP .le. 4 ) THEN
        IDGSTYP = 1
      ELSE IF ( PNTGRP .eq. 6 .or. PNTGRP .eq. 5 ) THEN
        IDGSTYP = 3
      ELSE IF ( PNTGRP .ge. 7 .and. PNTGRP .le. 9 ) THEN
        IDGSTYP = 2
      END IF
C
C     operator type
      IF (IDGSTYP .le. 2 ) THEN
         IHTYPE = 4
      ELSE IF (IDGSTYP.eq.3) THEN
         IHTYPE = 5
      ELSE
        WRITE(LUWRT,'(/A,I5)')'  *** ERROR in GET_XPROP_REOD_KRCI.'//
     &                        ' Unallowed value of IDGSTYP :',IDGSTYP
        CALL QUIT(' *** ERROR in GET_XPROP_REOD_KRCI: 
     &               Unallowed value of IDGSTYP.')
      END IF
C     transfer to common block in krciprop.h
      IHTYPE_X  = IHTYPE
C     set memory mark
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK ',IDUM,'PROPRO')
C
C     allocate for property operator matrix
C
      CALL MEMMAR(KXPROPM,4*N2ASHX*MZ,'ADDL  ',2,'IXPRPM')
C     do reordering - this flag controls behavior in subroutine ITRCTL
      DOPROPREOD  = .TRUE.
C
C     ... loop over property operators
      DO 100 II = 1, NPROP_KRCI
C
        ISYM_T = ISYMOPPRP_KRCI(II)
C       ... nothing to do for "inactive" operators
        IF( ISYM_T .le. 0 ) GOTO 100
C
C       ... no two-electron terms
        NO_TTWO         = 1
        IPRHAM_CI       = 0
        INDXPR          = LPROP_KRCI(II)
        ISYM_XXX        = IPRPSYM(INDXPR)
C
C
        CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
C
C       generate symmetry info for PICASSO
        CALL LUC2PIC(ISO1,IPRTRA)
C
C       read property matrix from file
        CALL DZERO(WORK(KXPROPM),4*N2ASHX*MZ)
        FILELAB = PRPNAM(INDXPR)(1:8)
        CALL KRCI_PRPFILE(LUF1MAT,FILELAB,WORK(KXPROPM),
     &                    4*N2ASHX*MZ,1)
C       ... reorder always at this point: RUN_SHMEM_T = .FALSE.
        RUN_SHMEM_T = .FALSE.
#if defined DEBUG_LUCIPROP
        WRITE(LUWRT,'(/A,A8/)') '  Reorder step for property matrix: ',
     &                             FILELAB
        CALL WRTMATMN(WORK(KXPROPM),1,4*N2ASHX*MZ,1,4*N2ASHX*MZ,LUWRT)
#endif
        CALL PICASSO(INTSPACE,ISO1,IPRTRA,MCSCF,NL2D,IUB2D,ISM2D,
     &               WORK(KXPROPM),DUMMY,DUMMY,DUMMY,IDUMMY,IBOSYM,
     &               CIRUN,RUN_SHMEM_T)
C
C       write sorted property matrix back to file
        WRITE(FILELAB_REOD,'(A6,A2)') PRPNAM(INDXPR)(1:6),RODLAB
#if defined DEBUG_LUCIPROP
!SK DEBUG
        WRITE(LUWRT,'(/A,A8/)') '  Final matrix on file with label: ',
     &                             FILELAB_REOD
        WRITE(LUWRT,'(A,I18/)') '  total number of integrals on KT_CC ',
     &                             NALLINT
        IF( FILELAB_REOD .eq. 'Y dipoRO')THEN
          WRITE(LUWRT,'(/A)') '  DEBUG READING DUMMY Y-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call alloc(XXXM,NALLINT)
          allocate(XXXM(NALLINT))
          open(unit=19,file='y_int_test', 
     &         status='UNKNOWN',form='FORMATTED')
          rewind(19)
          READ(19,'(1D22.14)') (XXXM(I),I=1,NALLINT)
          close(19)
          LINT_WRITE = NALLINT
          WRITE(LUWRT,'(/A)') '  DEBUG WRITING DUMMY Y-OPERATOR MATRIX'
          CALL WRTMATMN(XXXM,1,LINT_WRITE,1,LINT_WRITE,LUWRT)
          CALL DCOPY(LINT_WRITE,XXXM,1,WORK(KT_CC),1)
          WRITE(LUWRT,'(/A/)') '  USING DEBUG DUMMY Y-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call dealloc(XXXM)
          deallocate(XXXM)
        END IF
        IF( FILELAB_REOD .eq. 'X dipoRO')THEN
          WRITE(LUWRT,'(/A)') '  DEBUG READING DUMMY X-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call alloc(XXXM,NALLINT)
          allocate(XXXM(NALLINT))
          open(unit=19,file='x_int_test', 
     &         status='UNKNOWN',form='FORMATTED')
          rewind(19)
          READ(19,'(1D22.14)') (XXXM(I),I=1,NALLINT)
          close(19)
          LINT_WRITE = NALLINT
          WRITE(LUWRT,'(/A)') '  DEBUG WRITING DUMMY X-OPERATOR MATRIX'
          CALL WRTMATMN(XXXM,1,LINT_WRITE,1,LINT_WRITE,LUWRT)
          CALL DCOPY(LINT_WRITE,XXXM,1,WORK(KT_CC),1)
          WRITE(LUWRT,'(/A/)') '  USING DEBUG DUMMY X-OPERATOR MATRIX'
!radovan: presently alloc cannot handle integer(8) size with ILP32
!         call dealloc(XXXM)
          deallocate(XXXM)
        END IF
!SK DEBUG
#endif
C
C       store number of non-zero integrals for this operator
        IN1ELPRP_KRCI(II) = 0
        LINT_WRITE        = NALLINT
        IN1ELPRP_KRCI(II) = LINT_WRITE
#if defined DEBUG_LUCIPROP
        CALL WRTMATMN(WORK(KT_CC),1,LINT_WRITE,1,LINT_WRITE,LUWRT)
#endif
        CALL KRCI_PRPFILE(LUF1MAT,FILELAB_REOD,WORK(KT_CC),
     &                    LINT_WRITE,0)
        XT_TEST = 0.0D0
        XT_TEST = DDOT(LINT_WRITE,WORK(KT_CC),1,WORK(KT_CC),1)
C
C       reset property operator entries for operator with zero entries 
C       in reordered integral matrix
        IF( XT_TEST .eq. 0.0D0 ) THEN
           ISYMOPPRP_KRCI(II) = 0
           IN1ELPRP_KRCI(II)  = 0
        END IF
 100  CONTINUE
C     reordering done.
      DOPROPREOD  = .FALSE.
C
C     release memory
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM',IDUM,'PROPRO')
C
 999  CALL QEXIT('GET_XPROP_REOD_KRCI')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XPROP_KRCI(XPROPMAT,IHTYPE,IPROCLIST,IGROUPLIST,CIRUN)
C***********************************************************************
C
C     calculate property matrix. 
C     Reordered property integrals (KR-CI SIGDEN formalism) are stored 
C     on file KRMC_FOCK (file handle LUPROP_MAT).
C                                  
C     Written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
      use krci_properties
#include "implicit.h"
#include "priunit.h"
* boson symmetry info
#include "maxorb.h"
#include "../krmc/dcbbos.h"
#include "dcborb.h"
#include "dgroup.h"
* memory
#include "wrkspc.inc"
#include "ipoist8.inc"
C print flags
#include "cprnt.inc"
C required common blocks
#include "mxpdim.inc"
#include "symm.inc"
#include "lucinp.inc"
#include "crun.inc"
#include "cands.inc"
#include "cicisp.inc"
#include "cgas.inc"
#include "oper.inc"
#include "cintfo.inc"
#include "clunit.inc"
#include "picasso.inc"
#include "ctcc.inc"
#include "glbbas.inc"
#include "cstate.inc"
C information for property run
#include "krciprop.h"
#include "dcbxpr.h"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
#endif
#include "parluci.h"
      DIMENSION XPROPMAT(NPROP_ROOTS_KRCI,NPROP_ROOTS_KRCI,NPROP_KRCI,*)
      DIMENSION IPROCLIST(*), IGROUPLIST(*)
      CHARACTER CIRUN*6,FILELAB*8,FILELAB_REOD*8,RODLAB*2,x_NAME*800
      LOGICAL MCSCF, RUN_SHMEM_T, COPY_NEW_C, EX
      INTEGER*8 IL_KTCC8
      integer mz
C
      CALL QENTER('XPROP_KRCI')
      NTESTX = 10
C
C     set memory mark
      IDUM = 0
      CALL MEMMAR(KDUMMY,IDUM,'MARK  ',IDUM,'XPROP ')
C
      RODLAB = 'RO'
      MZ = MIN(NZ,2)
C     set required variables
      IMEMCK  = 0
      IPERTOP = 0
      INTIMP  = 7
C     ... no two-electron excitations
      NO_TTWO   = 1
      MCSCF     = .TRUE.
      DOSIGPROP = .TRUE.
C     hardwired to 1
      ICSPC  = NCMBSPC
      ISSPC  = NCMBSPC
C     real or complex eigenvectors
      IRC = 2
      IF( NZ .eq. 1 ) IRC = 1
      IRC_SAVE = IRC
C     scratch
      NBATXS = 0
C
C     property calculation setup section
C     __________________________________
C
      IF(MYPROC .eq. MASTER ) REWIND LUPROP_MAT
C
C     determine max. dimension of T and VECx( x = 1-3)
C     ------------------------------------------------
      IL_KTCC          = 1
      L_COMBI_MAX      = 0
      L_COMBI_MAX_SCR  = 0
      L_COMBI_MAX_SCRX = 0
C
      DO I = 1, NPROP_KRCI
        IL_KTCC  = MAX(IL_KTCC,IN1ELPRP_KRCI(I))
      END DO
      IL_KTCC8 = IL_KTCC
C
      CALL IZERO(IXSYMLIST,MXPROPKRCI_SYM)
C
C     IXSYMLIST will now be initialized and later adapted in 
C     Z_BLKFO_XPROP such that, e.g. IXSYMLIST(1) == '4', which is 
C     therefore the first "active" symmetry irrep, gives the 
C     right offset for NPARBLOCK(*,'1').
C
      DO IPROP_X = 1, NPROP_KRCI
         IF( ISYMOPPRP_KRCI(IPROP_X) .gt. 0 ) THEN
           DO I = 1, NPROP_ROOTS_KRCI
              ILC_C = XISPSM(ISYMEIG_KRCI(I),ICSPC)
              JSSM   = 
     &        IDBGMULT(ISYMOPPRP_KRCI(IPROP_X),ISYMEIG_KRCI(I))
              ILC_S  = XISPSM(JSSM,ISSPC)
C             mark "active" C/sigma symmetry irreps
              IXSYMLIST(ISYMEIG_KRCI(I)) = ISYMEIG_KRCI(I)
              IXSYMLIST(JSSM)  = JSSM
              L_COMBI_MAX_SCRX = MAX(ILC_C,ILC_S)
              L_COMBI_MAX_SCR  = MAX(L_COMBI_MAX_SCR,L_COMBI_MAX_SCRX)
           END DO
         END IF
      END DO
C
C     store symmetry pairs - active pairs are marked by 1.
      CALL MEMMAR(KXSYMPAIRS,MXPROPKRCI_SYM**2,'ADDL  ',1,'KXPAIR')
      CALL IZERO(WORK(KXSYMPAIRS),MXPROPKRCI_SYM**2)
C
      CALL GET_SC_PAIRS(WORK(KXSYMPAIRS),MXPROPKRCI_SYM,NPROP_KRCI,
     &                  ISYMOPPRP_KRCI,NPROP_ROOTS_KRCI,ISYMEIG_KRCI)
C
C     largest number of dets - used for MPI file offset calculation
      L_COMBI_MAX = MAX(L_COMBI_MAX_SCR,L_COMBI_MAX)
C     initialize L_COMBI to avoid crash in block distribution 
      L_COMBI     = L_COMBI_MAX
C
C     T block allocation
      CALL MEMMAR_I8(KT_CC,IL_KTCC8,'ADDL  ',2,'KT_CC ')
      CALL DZERO8(WORK(KT_CC),IL_KTCC8)
C     number of blocks per symmetry irrep
      CALL MEMMAR(KXSYM_BLK,MXPROPKRCI_SYM,'ADDL  ',1,'KXSYMB')
      CALL IZERO(WORK(KXSYM_BLK),MXPROPKRCI_SYM)
C
C     set block (parallel: distribution list) arrays and vector dimensions 
      CALL Z_BLKFO_XPROP(IDC,NMS2VAL,ICSPC,KNODELIST,KBLOCKLIST,
     &                   KRCCTOS,KXSYMDISTR,WORK(KXSYM_BLK),
     &                   WORK(KXSYMPAIRS),IRC,IGROUPLIST)
C
      CALL GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C
C     set file offset/view for MPI-I/O files
#if defined (VAR_MPI2)
C     set active blocks for C scratch file ILUC
      MY_ACT_BLK_ALL = NUM_BLOCKS
C     general offset
      MY_LUC_OFF = 0
C
      CALL CALC_OFF_MPI_FILE_XPRP(IGROUPLIST,WORK(KBLOCKLIST),
     &                            WORK(KNODELIST),NPROP_ROOTS_KRCI,
     &                            NPROP_ROOTS_KRCI)
C
C     length for allocation of file arrays
C     WRITE(LUWRT,*) ' NPROP_ROOTS_KRCI,MY_ACT_BLK2,NUM_BLOCKS,...',
C    &NPROP_ROOTS_KRCI,MY_ACT_BLK2,NUM_BLOCKS,NUM_BLOCKS2
      IALL_LU1 = NPROP_ROOTS_KRCI * MY_ACT_BLK2
      IALL_LU2 = NPROP_ROOTS_KRCI * MY_ACT_BLK2
      IALL_LUC = NUM_BLOCKS2
C     WRITE(LUWRT,*) ' RMEMCHK before list alloc',
C    &                 IALL_LU1,IALL_LU2,IALL_LUC
C     CALL RMEMCHK
C
C     allocate file arrays
      CALL MEMMAR(KILU1LIST,IALL_LU1,'ADDS  ',1,'LU1LST')
      CALL MEMMAR(KILU2LIST,IALL_LU2,'ADDS  ',1,'LU2LST')
      CALL MEMMAR(KILUCLIST,IALL_LUC,'ADDS  ',1,'LUCLST')
C     initialize ...
      CALL IZERO(WORK(KILU1LIST),IALL_LU1)
      CALL IZERO(WORK(KILU2LIST),IALL_LU2)
      CALL IZERO(WORK(KILUCLIST),IALL_LUC)
C
C     T block distribution - reordered arrays of property integrals
      CALL DIST_XPROP(WORK(KT_CC),LUPROP_MAT)
C     WRITE(LUWRT,*) ' RMEMCHK before DIST_XCVEC'
C     CALL RMEMCHK
C
C     distribute C vectors from MASTER to slaves for all active irreps
      CALL DIST_XCVEC(WORK(KVEC1),LUC,ILU1,WORK(KILU1LIST),
     &                WORK(KNODELIST),WORK(KBLOCKLIST),WORK(KXSYM_BLK),
     &                IXSYMLIST,IXSYMACT,IRC)
C     WRITE(LUWRT,*) ' RMEMCHK after DIST_XCVEC'
C     CALL RMEMCHK
C     scratch arrays
      CALL MEMMAR(KXR,NPROP_ROOTS_KRCI,'ADDL  ',2,'KXR   ')
      IF(MZ.eq.2) THEN 
        CALL MEMMAR(KXI,NPROP_ROOTS_KRCI,'ADDL  ',2,'KXI   ')
      ELSE
        CALL MEMMAR(KXI,               0,'ADDL  ',2,'KXI   ')
      END IF
      CALL DZERO(WORK(KXR),NPROP_ROOTS_KRCI)
      IF( MZ .eq. 2 ) CALL DZERO(WORK(KXI),NPROP_ROOTS_KRCI)
C
C     end of setup section
C     --------------------
#endif
C
C     computation section
C     ___________________
C
C     initialize start codes for ICSM_OLD and NVC_LUC (# of vec per irrep)
      ICSM_OLD  = -1
      NVC_LUC   = -3
#if defined (VAR_MPI2)
C     ... C batches of blocks 
      CALL MEMMAR(KPCLBXX,MXNTTS,  'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEXX,MXNTTS,  'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1XX,MXNTTS,  'ADDL  ',1,'CI1BTX')
      CALL MEMMAR(KPCIBXX,8*MXNTTS,'ADDL  ',1,'CIBTX ')
C     ... sigma batches of blocks
      CALL MEMMAR(KPCLBXXS,MXNTTS,  'ADDL  ',1,'CLBTXY')
      CALL MEMMAR(KPCLEXXS,MXNTTS,  'ADDL  ',1,'CLEBTY')
      CALL MEMMAR(KPCI1XXS,MXNTTS,  'ADDL  ',1,'CI1BTY')
      CALL MEMMAR(KPCIBXXS,8*MXNTTS,'ADDL  ',1,'CIBTXY')
#endif
      WRITE(LUWRT,'(/A/)') 
     & '  ********* start of property matrix calculation *********'
C
C     outer loop over all eigenstates
      DO 100 ISTATE_X = 1, NPROP_ROOTS_KRCI
C
         ICSM     = ISYMEIG_KRCI(ISTATE_X)
         L_COMBI  = XISPSM(ICSM,ICSPC)
         IF( ICSM .ne. ICSM_OLD )
     &     WRITE(LUWRT,'(/A,I4/)') 
     &     '  * looping now over all states of symmetry',ICSM
#if defined (VAR_MPI2)
C        erase T operator file offset
         IOFFSET_SCRATCH = 0
C        get block partitioning for new active irrep
         IF( ICSM .ne. ICSM_OLD )THEN
C          which "active" C irrep?
           IX_CT = IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ICSM)
C
           CALL Z_BLKFO_XPRP_PAR(IDC,NMS2VAL,ICSPC,ICSM,WORK(KPCLBXX),
     &                           WORK(KPCLEXX),WORK(KPCI1XX),
     &                           WORK(KPCIBXX),NBATX,WORK(KNODELIST),
     &                           IX_CT)
         END IF
#else
         NVC_LUC = NVC_LUC + 1
         IF( ICSM .ne. ICSM_OLD )THEN
C          close C vector file for previous irrep (-2 == start code +1)
           IF( NVC_LUC .ne. -2 ) CLOSE(UNIT=LUC,STATUS='KEEP')
C          reset NVC_LUC and open new C vector file
           NVC_LUC = 1
           CALL GET_LUC_FILE(LUC,ICSM,EX)
         END IF
         REWIND LUPROP_MAT
#endif
C        new ICSM_OLD
         ICSM_OLD = ICSM
C        copy new C vector blocks from ILU1 to ILUC
         COPY_NEW_C = .TRUE.
C
C        loop over all property operators
         DO 200 IPROP_X = 1, NPROP_KRCI
C
            ISYM_T = ISYMOPPRP_KRCI(IPROP_X)
C           only zero elements in property operator matrix
            IF( ISYM_T .le. 0 ) GOTO 200
C
C           property pointer
            INDXPR      = LPROP_KRCI(IPROP_X)
C
C           reset operator specific common blocks
            CALL SET_HOP_DBG(IHTYPE,NO_TTWO,ISYM_T,IPRHAM_CI)
C
            LENGTH_BLK = IN1ELPRP_KRCI(IPROP_X)
C           erase old property integrals
            CALL DZERO8(WORK(KT_CC),IL_KTCC8)
C
C           read new property integrals from disk
#if defined (VAR_MPI2)
            CALL MPI_FILE_READ_AT(ILPRP_X,IOFFSET_SCRATCH,WORK(KT_CC),
     &                            LENGTH_BLK,MPI_REAL8,ISTAT,IERR)
            IOFFSET_SCRATCH = IOFFSET_SCRATCH + LENGTH_BLK
C
#ifdef MOD_DEBUG
C
C           debug print
            WRITE(LUWRT,*) ' matrix of property:',PRPNAM(INDXPR)
            CALL WRTMATMN(WORK(KT_CC),1,LENGTH_BLK,1,LENGTH_BLK,LUWRT)
#endif
#else
            WRITE(FILELAB_REOD,'(A6,A2)') PRPNAM(INDXPR)(1:6),RODLAB
            CALL KRCI_PRPFILE(LUPROP_MAT,FILELAB_REOD,WORK(KT_CC),
     &                        LENGTH_BLK,1)
#endif
            ISSM      = IDBGMULT(ICSM,ISYM_T)
            L_COMBI_S =   XISPSM(ISSM,ISSPC)
C           which "active" S irrep?
            IX_ST     = IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ISSM)
C
C           compute new sigma vector = T_prop x C
C           -------------------------------------
            CALL XSIGDEN_CTRL(WORK(KVEC1),WORK(KVEC2),WORK(KT_CC),
     &                        COPY_NEW_C,NVC_LUC,MZ
#if defined (VAR_MPI2)
     &                        ,WORK(KILU1LIST),WORK(KILU2LIST),
     &                         WORK(KILUCLIST),WORK(KPCLBXX),
     &                         WORK(KPCLEXX),WORK(KPCI1XX),
     &                         WORK(KPCIBXX),NBATX,WORK(KBLOCKLIST),
     &                         WORK(KNODELIST),IGROUPLIST,IPROCLIST,
     &                         WORK(KRCCTOS),WORK(KXSYMPAIRS),ICSM,ISSM,
     &                         IX_CT,IX_ST,ISTATE_X-1
#endif
     &                        )

C
C           inner loop over all eigenstates
#if defined (VAR_MPI2)
C 
#ifdef MOD_DEBUG
            WRITE(LUWRT,*) '<S_ISSM|prop|C_ICSM> for C_# and IPROP_X',
     &                         ISSM,PRPNAM(INDXPR),ICSM,ISTATE_X,IPROP_X
#endif
           CALL Z_BLKFO_XPRP_PAR(IDC,NMS2VAL,ISSPC,ISSM,WORK(KPCLBXXS),
     &                           WORK(KPCLEXXS),WORK(KPCI1XXS),
     &                           WORK(KPCIBXXS),NBATXS,WORK(KNODELIST),
     &                           IX_ST)

            CALL INPRDD_REAL_CPLX_XPRP(ILU2,ILU1,WORK(KVEC1),WORK(KVEC2)
     &                                ,WORK(KXR),WORK(KXI),NBATXS,
     &                                 WORK(KPCLBXXS),WORK(KPCLEXXS),
     &                                 WORK(KPCI1XXS),WORK(KPCIBXXS),
     &                                 MY_LU2_OFF,MY_LU1_OFF,
     &                                 WORK(KILU2LIST),WORK(KILU1LIST),
     &                                 NPROP_ROOTS_KRCI,ISSM,
     &                                 ISYMEIG_KRCI,MZ,1,ISTATE_X)
C
C           debug printing
#ifdef MOD_DEBUG
C
            WRITE(LUWRT,*) ' WORK(KXR) '
            CALL WRTMATMN(WORK(KXR),1,NPROP_ROOTS_KRCI,1,
     &                    NPROP_ROOTS_KRCI,LUWRT)
            IF( MZ.eq.2)THEN
              WRITE(LUWRT,*) ' WORK(KXI) '
              CALL WRTMATMN(WORK(KXI),1,NPROP_ROOTS_KRCI,1,
     &                      NPROP_ROOTS_KRCI,LUWRT)
            END IF
#endif
C
C           lz is pure imaginary operator - scale with -1.0
C           FIXME: introduce proper scaling for each imaginary 
C                  operator
            DMONE = -1.0D0
            IF( PRPNAM(INDXPR) .eq. 'z-angular mom' )THEN
              CALL DSCAL(NPROP_ROOTS_KRCI,DMONE,WORK(KXR),1)
              IF(MZ.eq.2) CALL DSCAL(NPROP_ROOTS_KRCI,DMONE,WORK(KXI),1)
            END IF
C
            CALL DCOPY(NPROP_ROOTS_KRCI,WORK(KXR),1,
     &                 XPROPMAT(1,ISTATE_X,IPROP_X,1),1)
            IF( MZ .eq. 2 ) 
     &      CALL DCOPY(NPROP_ROOTS_KRCI,WORK(KXI),1,
     &                 XPROPMAT(1,ISTATE_X,IPROP_X,2),1)
C
C           DO 300 ISTATE_Y = 1, NPROP_ROOTS_KRCI
C             IF( ISSM .ne. ISYMEIG_KRCI(ISTATE_Y) ) GOTO 300
C             XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,1) =
C    &                                              WORK(KXR+ISTATE_Y-1)
C             WRITE(LUWRT,*) 'ISTATE_Y,ISTATE_X, PRPNAM(INDXPR) and'//
C    &                       ' XPROPMAT(Y,X,PROP,1)',ISTATE_Y,ISTATE_X,
C    &                        PRPNAM(INDXPR),
C    &        XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,1)
C             IF( MZ .eq. 2 ) XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,2) =
C    &                                              WORK(KXI+ISTATE_Y-1)
C             WRITE(LUWRT,*) 'ISTATE_Y,ISTATE_X, PRPNAM(INDXPR) and'//
C    &                       ' XPROPMAT(Y,X,PROP,2)',ISTATE_Y,ISTATE_X,
C    &                        PRPNAM(INDXPR),
C    &        XPROPMAT(ISTATE_Y,ISTATE_X,IPROP_X,2)
C300        CONTINUE
C           ^ end of inner loop over all eigenstates
#else
#ifdef MOD_DEBUG
            WRITE(LUWRT,*) '<S_ISSM|prop|C_ICSM> for C_# and IPROP_X',
     &                         ISSM,PRPNAM(INDXPR),ICSM,ISTATE_X,IPROP_X
#endif
!           close current c-vector file
            close(unit=luc,status="keep")
            is_correct_sym = 0
            
            do istate_y = 1, nprop_roots_krci
C
              call get_inner_product_sigmaprop_dot_c(work(kvec1),
     &                                               work(kvec2),
     &                                               xpropmat(
     &                                               istate_y,
     &                                               istate_x,
     &                                               iprop_x,1),
     &                                               xpropmat(
     &                                               istate_y,
     &                                               istate_x,
     &                                               iprop_x,mz),
     &                                               issm,            
     &                                               isymeig_krci(
     &                                               istate_y),
     &                                               luhc,           
     &                                               is_correct_sym,    
     &                                               mz,     
     &                                               xsymflab(
     &                                               isymeig_krci(
     &                                               istate_y)))

              if(prpnam(indxpr).eq.'z-angular mom')then
                            xpropmat(istate_y,istate_x,iprop_x,1) = 
     &                     -xpropmat(istate_y,istate_x,iprop_x,1)
                if(mz.eq.2) xpropmat(istate_y,istate_x,iprop_x,2) = 
     &                     -xpropmat(istate_y,istate_x,iprop_x,2)
              end if
C
            end do
C           ^ end of inner loop over all eigenstates

C           reopen current c-vector file
            call get_luc_file(luc,icsm,ex)
#endif
 200     CONTINUE
C        ^ end of loop over all properties
 100  CONTINUE
C     ^ end of outer loop over all eigenstates
C
      WRITE(LUWRT,'(/A/)') 
     & '  *********  end of property matrix calculation  *********'
C
#if defined (VAR_MPI2)
C
C     update of complete property matrix
C     ----------------------------------
      CALL MEMMAR(KXPROP_SCR,LXPRPKRCI,'ADDL  ',2,'XPRPSC')
      CALL DZERO(WORK(KXPROP_SCR),LXPRPKRCI)
C
      CALL REDVEC_REL(XPROPMAT,WORK(KXPROP_SCR),LXPRPKRCI,2,
     &                MPI_SUM,MPI_COMM_WORLD,MASTER)
      IF( MYPROC .eq. MASTER) 
     &    CALL DCOPY(LXPRPKRCI,WORK(KXPROP_SCR),1,XPROPMAT,1)
#endif
C
C     print section
C     _____________
C
      NTESTX = 00
      IF( NTESTX .gt. 5 .and. myproc.eq.0)THEN
        WRITE(LUWRT,*) ' calculation of XPROPMAT finished'
        WRITE(LUWRT,*) ' content of XPROPMAT'
        CALL WRTMATMN(XPROPMAT,1,LXPRPKRCI,1,LXPRPKRCI,LUWRT)
      END IF
      NTESTX = 00
C     release memory
 9999 call MEMMAR(KDUMMY,IDUMMY,'FLUSM ',IDUM,'XPROP ')
C
#if defined (VAR_MPI2)
C
C     ... close scratch file with reordered property integrals (opened
C     in DIST_XPROP)
      CALL MPI_FILE_CLOSE(ILPRP_X,IERR)
#endif
C
C     set back to default
      DOSIGPROP = .FALSE.
C      
      CALL QEXIT('XPROP_KRCI')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     end of file property_r.F
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
