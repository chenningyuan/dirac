C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mpglvc */
      SUBROUTINE MPGLVC(DMO,XLAGRAN,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: Driver routine for the construction of the MP2
C              Lagrangian.
C
C     INPUT:   DMO - second order density matrix, not sym. packed
C
C     OUTPUT:  XLAGRAN - The Lagrangian, not sym. packed
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION DMO(*)
      DIMENSION XLAGRAN(*)
      DIMENSION WORK(*)
C
      CALL QENTER('MPGLVC')
C
      CALL DZERO(XLAGRAN,NORBT*NORBT*NZ)
C
C     Calculate density matrix dependent part of Lagrangian
C     =====================================================
C
      CALL LVCDEN(DMO,XLAGRAN,WORK,KFREE,LFREE)
C
      CALL QEXIT('MPGLVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcden */
      SUBROUTINE LVCDEN(DMO,XLAGRAN,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: generate density matrix dependent part of Lagrangian
C
C     INPUT:   DMO - second order density matrix
C
C     OUTPUT:  XLAGRAN - The Lagrangian
C
C     REMARKS: 
C
C              L_(pq)=-2<kp||jq>DM_(jk) - 2<cp||bq>DM_(bc)
C
C              with: j and k active occupied orbitals
C                    b and c active virtual orbitals
C
C              pq corresponds the following type of pairs:
C
C                    active virtual - active occupied
C                    active virtual - inactive occupied
C
C                    inactive virtual - active occupied
C                    inactive virtual - inactive occupied
C
C                    positronic - active occupied
C                    positronic - inactive occupied
C              
C              using the definitions for the density matrices:
C
C              DA_(yx)(o-o)=C_(yj)(D_(jk))(C+)_(xk) for the 
C              occupied-occupied block
C                 and
C              DA_(yx)(v-v)=C_(yb)(D_(bc))(C+)_(xc) for the
C              virtual-virtual block
C
C              We can write
C
C              L_(pq)=-2[g_(pq,xy)-g_(py,xq)]DA_(yx)
C
C              with DA_(yx)=DA_(yx)(o-o) + DA_(yx)(v-v)
C
C              L_(pq) now has the same form as the G-matrix (i.e. 
C              a Fock type matrix) used in the sigma vector, which is
C              used in RPA (See eq 67 of T.Saue and H.J.Aa Jensen,
C              J. Chem. Phys. 118, 522 (2003)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "mxorb.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
#include "blocks.h"
#include "dcbfir.h"
C
      DIMENSION DMO(*)
      DIMENSION XLAGRAN(*)
      DIMENSION WORK(*)
C
      CALL QENTER('LVCDEN')
C
      KFRSAV = KFREE
C            
C     New NPOS(*) array for HER_PARDRV, reset dcbfir.h to true
C            
      NPOS = 0
      IF(PARCAL) THEN
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2 
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      ENDIF       
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
C
C     Memory allocation
C     ====================================================
      CALL MEMGET('REAL',KCMO   ,N2BBASXQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXLAGAO,NZ*N2BBASX,WORK,KFREE,LFREE)
C
      CALL LVCDEN1(DMO,XLAGRAN,WORK(KCMO),WORK(KDMAT),WORK(KXLAGAO),
     &             WORK,KFREE,LFREE)
C
C     Release memory
C     ==============
      CALL MEMREL('LVCDEN',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('LVCDEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcden1 */
      SUBROUTINE LVCDEN1(DMO,XLAGRAN,CMO,DMAT,XLAGAO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: generate density matrix dependent part of Lagrangian
C
C     INPUT:   DMO - second order density matrix
C
C     OUTPUT:  XLAGRAN - The Lagrangian in spinor (mo) basis
C
C     Other arguments:
C              CMO  - spinor (mo) coefficients
C              DMAT - AO density matrix (or better to say density 
C                     matrix transformed to 'orthonormal SO' basis 
C              XLAGAO - Lagrangian transformed to AO basis
C
C     REMARKS:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "cbihr2.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
C
      LOGICAL TOBE
      DIMENSION DMO(*)
      DIMENSION XLAGRAN(*)
      DIMENSION CMO(*), DMAT(*), XLAGAO(N2BBASXQ)
      DIMENSION WORK(*)
C
      CALL QENTER('LVCDEN1')
C
C     a) Get coefficients
C     -------------------
C
      LUBUF = 1
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','LVCDEN1')
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,DUM,TOTERG,2,WORK,KFREE,LFREE)
      CLOSE(LUBUF,STATUS='KEEP')
C
C     b) Generate AO density matrix
C     -----------------------------
C
c     CALL PRQMAT(DMO,NORBT,NORBT,NORBT,NORBT,NZ,
c    &            IPQTOQ(1,0),LUPRI)
      CALL DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C
C     c) Generate AO Lagrangian
C     -------------------------
C
C     jvs: put integral flag now by hand, later via commonblock
C
      ILL = 1
      ILS = 1
      ISS = 0
      INTFLG = ILL + 2 * ILS + 4 * ISS 
C
      ISYMLAG = 1 ! Lagrangian is totally symmetric
      IHRMOP = 1  ! Hermitian operator
      IFCKOP = 1  ! 'Fock' matrix type
      NRFMAT = 1  ! nr of Fock-type matrices
C
      CALL DZERO(XLAGAO,N2BBASXQ)
C
      CALL TWOFCK(ISYMLAG,IHRMOP,IFCKOP,
     &            XLAGAO,DMAT,NRFMAT,NPOS,INTFLG,
     &            IPRTWO,WORK(KFREE),LFREE)
C
C     d) Transform Lagrangian to MO basis
C     -----------------------------------
C
      CALL LTAOMO(XLAGAO,XLAGRAN,CMO,WORK,KFREE,LFREE)
C
      CALL QEXIT('LVCDEN1')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'LVCDEN1: Coefficient file not found !'
      CALL QUIT('LVCDEN1: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmp2ao */
      SUBROUTINE DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: Transform MO MP2 densitry matrix to quaternionic
C              symmetry-packed AO matrix
C
C     INPUT:   DMO  - MO density matrix
C              CMO  - spinor (MO) coefficients
C
C     OUTPUT:  DMAT - AO MP2 denstity matrix 
C
C     REMARKS:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER(D0 = 0.0D0, DP5 = 0.5D00)
C         
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),DMO(NORBT,NORBT,NZ),WORK(*)
C           
      EXTERNAL BDQALGEBRA
C
      CALL DZERO(DMAT,N2BBASXQ)
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for first index
C
          NS = NORB(IFSYM)
          IS = IORB(IFSYM) + 1
          JS = ICMOQ(IFSYM) + 1
C
C        Parameters for first index
C
          NI  = NORB(IFSYM)
          II  = IORB(IFSYM) + 1
          JI  = ICMOQ(IFSYM) + 1
C
          CALL QTRANS('MOAO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                DMAT(I2BASX(IFSYM,IFSYM)+1,1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                DMO(IS,II,1),NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1)
C
      ENDDO
c       CALL HEADER('DMP2AO: AO density matrix',-1)
c       CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
c    &              NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
c       CALL FLSHFO(LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ltaomo */
      SUBROUTINE LTAOMO(XLAGAO,XLAGRAN,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: Transform Lagrangian from AO-basis to MO-basis
C
C     INPUT:   XLAGAO - Lagrangian in AO-basis
C              CMO  - spinor (MO) coefficients
C
C     OUTPUT:  XLAGRAN - The Lagrangian in spinor (MO) basis, not
C                        symmetry packed
C
C     REMARKS:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1=1.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION XLAGRAN(NORBT,NORBT,NZ)
      DIMENSION XLAGAO(N2BBASXQ), CMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('LTAOMO')
C
C     +--------------+
C     |   e-e part   |
C     +--------------+
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for 1st index: virtuals
C
         NVIRTL= NORB(IFSYM) - NPSH(IFSYM) - NOCC(IFSYM)
         NS = NVIRTL      ! nr of virtuals for ifsym
         IS = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM) + 1
         JS = ICMOQ(IFSYM) + (NPSH(IFSYM)+NOCC(IFSYM))*NFBAS(IFSYM,0) +1
C
C        Parameters for 2nd index: occupied
C
         NI = NOCC(IFSYM) 
         II = IORB(IFSYM) + NPSH(IFSYM) + 1
         JI = ICMOQ(IFSYM) + NPSH(IFSYM)*NFBAS(IFSYM,0) + 1
C
       IF(NS.NE.0.AND.NI.NE.0) THEN 
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                XLAGAO(I2BASX(IFSYM,IFSYM)+1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                XLAGRAN(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1) ! last one is printlevel
       ENDIF
      ENDDO
C
C     +--------------+
C     |   e-p part   |
C     +--------------+
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for 1st index: positrons
C
         NS = NPSH(IFSYM)      ! nr of positrons for ifsym
         IS = IORB(IFSYM) + 1
         JS = ICMOQ(IFSYM) +1
C
C        Parameters for 2nd index: occupied
C
         NI = NOCC(IFSYM)
         II = IORB(IFSYM) + NPSH(IFSYM) + 1
         JI = ICMOQ(IFSYM) + NPSH(IFSYM)*NFBAS(IFSYM,0) + 1
C
       IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                XLAGAO(I2BASX(IFSYM,IFSYM)+1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                XLAGRAN(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1) ! last one is printlevel
       ENDIF
      ENDDO
C
C     Scale the Lagrangian by -1
C     ==========================
      CALL DSCAL(NORBT*NORBT*NZ,-D1,XLAGRAN,1)
c     do iz = 1, nz
c        do j = 1, norbt
c           do i= 1, norbt
c              print*, 'lagr=', iz, i, j, XLAGRAN(i,j,iz) 
c           enddo
c        enddo
c     enddo
c     print*,'lagrangian'
c     CALL PRQMAT(XLAGRAN,NORBT,NORBT,NORBT,NORBT,NZ,
c    &            IPQTOQ(1,0),LUPRI)
C
      CALL QEXIT('LTAOMO')
      RETURN
      END
