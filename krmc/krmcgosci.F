C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_startvec */
      SUBROUTINE GOSCI_STARTVEC(CMO,CREF,WORK,LWORK)
C*****************************************************************************
C
C     Generate CI start vectors from GOSCI by diagonalizing the
C     full CI matrix.
C
C     Input : CMO - MO coefficents
C
C     Output: CREF - CI vector
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbopt.h>
#include <dgroup.h>
C
      DIMENSION WORK(*), CREF(NZCONF,*)
C
#include <memint.h>
C
      CALL QENTER('GOSCI_STARTVEC')
C
C     Calculate core Hamiltonian:
C     ---------------------------
C
C     The core Hamiltonian is the active-active part of FC,
C     a.k.a. FCAC.
C
      CALL MEMGET('REAL',KFCAC,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Calculate FCmo:
C
      CALL RGETFC(CMO,WORK(KFC),ECORE,WORK(KFREE),LFREE)
C
C     Get active-active part of FCmo:
C
      CALL RGETAC(WORK(KFC),WORK(KFCAC),IPROPT)
C
C     Transform quaternion FCACmo to molfdir type FCAC.
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MIN0(NZ,2),
     $     WORK,KFREE,LFREE)
      CALL QFC2MFC(WORK(KFCAC),WORK(KMOLFFC),IPROPT)
      CALL GOSCI_SETONE(WORK(KMOLFFC),2*NASHT)
C
      CALL MEMREL('after RGETFC',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Read two-electron integrals:
C     ----------------------------
C
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
C
C
      CALL MEMGET('REAL',KVUUUU,NASHT**4 * 4, WORK,KFREE,LFREE)
      CALL REAINTS(WORK(KVUUUU),DUMMY,DUMMY,DUMMY,NASHT,NORBT,1,JTRLVL,
     &     IPRNT,WORK(KFREE),LFREE)
      CALL D2MUUUU(WORK(KVUUUU),WORK(KMUUUU),10)
      CALL MEMREL('GOSCI_STARTVEC.D2MUUU',WORK,KWORK,KVUUUU,KFREE,LFREE)
C
C     Allocate memory for matrices:
C
      NZCONF2 = NZCONF ** 2
      CALL MEMGET('REAL',KHR,NZCONF2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHI,NZCONF2,WORK,KFREE,LFREE)
C
C     Construct Hamiltonian:
C     ----------------------
C
      CALL DZERO(WORK(KHR),NZCONF2)
      CALL DZERO(WORK(KHI),NZCONF2)
      CALL GOSCI_MAKEH(WORK(KHR),WORK(KHI), WORK(KMUUUU))
C
C
C     Diagonalize Hamiltonian:
C     ------------------------
C
      IRP = IOPT_SYMMETRY
      IF ( MOD ( NAELEC, 2 ) .EQ. 0 ) THEN
C        ...boson symmetry
         IRRP = IGOSCI_NREP() + IRP
      ELSE
C        ...fermion symmetry
         IRRP = IRP
      END IF
C
      IF ( IPROPT .GE. 10 ) THEN
         CALL HEADER('GOSCI_STARTVEC: Full CI matrix',-1)
         CALL OUTPUT(WORK(KHR),1,NZCONF,1,NZCONF,NZCONF,NZCONF,1,LUPRI)
         CALL OUTPUT(WORK(KHI),1,NZCONF,1,NZCONF,NZCONF,NZCONF,
     $           1,LUPRI)
      END IF
      CALL FLSHFO(LUPRI)
C
C
      CALL DIAGH(IRP,WORK(KHR),WORK(KHI))
C
C     Get CI vector:
C     --------------
C
C     The CI vector is stored in colunms; the lowest eigenstate is in
c     the first column.
C     
      IOFF = ( IOPT_ISTATE - 1 ) * NZCONF
      IF ( NZ .EQ. 1 ) THEN
         DO I = 1,NZCONF
            CREF(I,1) = WORK(KHR + IOFF + I - 1 )
         END DO
      ELSE
         DO I = 1,NZCONF
            CREF(I,1) = WORK(KHR + IOFF + I - 1 )
            CREF(I,2) = WORK(KHI + IOFF + I - 1 )
         END DO
      END IF
C
C DEBUG DEBUG DEBUG
      CREF(1,1) = 1.0D00
      DO I = 2, 4
         CREF(I,1) = 0.0D00
      END DO
C
      IF ( IPROPT .GE. 5 ) THEN
         CALL OUTPUT(WORK(KHR),1,NZCONF,1,NZCONF,NZCONF,NZCONF,1,LUPRI)
         CALL HEADER('GOSCI_STARTVEC: CI vector',-1)
         CALL PRQMAT(CREF,NZCONF,1,NZCONF,1,
     $        MIN0(NZ,2),IPQTOQ(1,0),LUPRI)
      END IF
C
C     Output section:
C     ---------------
C
      CALL PRTOUT(IRP,IRRP,WORK(KHR),WORK(KHI))
C
 10   CALL MEMREL('GOSCI_STARTVEC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('GOSCI_STARTVEC')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmulta */
      SUBROUTINE GOSCI_IREPS (NZ,NFSYM,NSYMRPA,REPA)
C*****************************************************************************
C     
C     Rip-off from GMULTA from MOLTRA./jth-03072000
C
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C
C     Luuk Visscher, Mar 18 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*4 REPA(8)
C
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
         NSYMRPA = 1
         REPA(1) = '   A'
         REPA(2) = '   a'
      ELSEIF (NFSYM.EQ.1.AND.NZ.LE.2) THEN
         NSYMRPA = 2
         REPA(1) = '  1E'
         REPA(2) = '  2E'
         REPA(3) = '   a'
         REPA(4) = '   b'
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
         NSYMRPA = 2
         REPA(1) = '  AG'
         REPA(2) = '  AU'
         REPA(3) = '  ag'
         REPA(4) = '  au'
      ELSEIF (NFSYM.EQ.2.AND.NZ.LE.2) THEN
         NSYMRPA = 4
         REPA(1) = ' 1Eg'
         REPA(2) = ' 2Eg'
         REPA(3) = ' 1Eu'
         REPA(4) = ' 2Eu'
         REPA(5) = '  ag'
         REPA(6) = '  bg'
         REPA(7) = '  au'
         REPA(8) = '  bu'
      ELSE
         CALL QUIT('Multiplication table is not defined')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck QFC2MFC */
      SUBROUTINE QFC2MFC(DQFC,DMFC,IPRINT)
C*****************************************************************************
C
C
C     Input : DQFC  - FCAC matrix in quaternion basis
C
C     Output: DMFC  - FCAC matrix in molfdir basis
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION DQFC(NASHT,NASHT,NZ)
      DIMENSION DMFC(2*NASHT,2*NASHT,2)
C
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output from QFC2MFC',-1)
         CALL HEADER('Input matrix:',-1)
         CALL PRQMAT(DQFC,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      IF ( NZ .EQ. 1 ) THEN
         CALL DZERO(DMFC,4 * NASHT**2)
      ELSE
         CALL DZERO(DMFC,4 * NASHT**2)
      END IF
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
         IOFF3 = IOFF3 + NASH(IFSYM)
         DO J = 1, NASH(IFSYM)
            CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,1),1,
     &           DMFC(1+IOFF2,J+IOFF2,1),1)
            CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,1),1,
     &           DMFC(1+IOFF3,J+IOFF3,1),1)
         ENDDO
         IF (NZ.GE.2) THEN
            DO J = 1, NASH(IFSYM)
               CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,2),1,
     &              DMFC(1+IOFF2,J+IOFF2,2),1)
               CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,2),1,
     &              DMFC(1+IOFF3,J+IOFF3,2),1)
               CALL DSCAL(NASH(IFSYM),DM1,DMFC(1+IOFF3,J+IOFF3,2),1)
            ENDDO
            IF (NZ.GE.4) THEN
               DO J = 1, NASH(IFSYM)
                  CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,3),1,
     &                 DMFC(1+IOFF2,J+IOFF3,1),1)
                  CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,3),1,
     &                 DMFC(1+IOFF3,J+IOFF2,1),1)
                  CALL DSCAL(NASH(IFSYM),DM1,DMFC(1+IOFF3,J+IOFF2,1),1)
                  CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,4),1,
     &                 DMFC(1+IOFF2,J+IOFF3,2),1)
                  CALL DCOPY(NASH(IFSYM),DQFC (1+IOFF1,J+IOFF1,4),1,
     &                 DMFC(1+IOFF3,J+IOFF2,2),1)
               ENDDO
            ENDIF
         ENDIF
         IOFF1 = IOFF1 + NASH(IFSYM)
         IOFF2 = IOFF2 + 2*NASH(IFSYM)
         IOFF3 = IOFF3 + NASH(IFSYM)
      ENDDO
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output matrix:',-1)
         CALL PRQMAT(DMFC,2*NASHT,2*NASHT,2*NASHT,2*NASHT,MIN0(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mfc2qfc */
      SUBROUTINE MFC2QFC(DMFC,DQFC,IPRINT)
C*****************************************************************************
C
C
C     Input : DMFC  - FCAC matrix in molfdir basis
C
C     Output: DQFC  - FCAC matrix in quaternion basis
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION DQFC(NASHT,NASHT,NZ)
      DIMENSION DMFC(2*NASHT,2*NASHT,2)
C
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output from MFC2QFC',-1)
         CALL HEADER('Input matrix:',-1)
         CALL PRQMAT(DMFC,2*NASHT,2*NASHT,2*NASHT,2*NASHT,MIN0(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL DZERO(DQFC, NASHT**2 * NZ )
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
         IOFF3 = IOFF3 + NASH(IFSYM)
         DO J = 1, NASH(IFSYM)
            CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF2,1), 1,
     $           DQFC (1+IOFF1,J+IOFF1,1),1 )
         ENDDO
         IF (NZ.GE.2) THEN
            DO J = 1, NASH(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF2,2), 1,
     $              DQFC (1+IOFF1,J+IOFF1,2),1 )
            ENDDO
         END IF
         IF (NZ.GE.4) THEN
            DO J = 1, NASH(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF3,1), 1,
     $              DQFC (1+IOFF1,J+IOFF1,3), 1)
               CALL DCOPY(NASH(IFSYM), DMFC(1+IOFF2,J+IOFF3,2), 1,
     $              DQFC (1+IOFF1,J+IOFF1,4), 1)
            ENDDO
         ENDIF
         IOFF1 = IOFF1 + NASH(IFSYM)
         IOFF2 = IOFF2 + 2*NASH(IFSYM)
         IOFF3 = IOFF3 + NASH(IFSYM)
      ENDDO
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output matrix:',-1)
         CALL PRQMAT(DQFC,NASHT,NASHT,NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mfq2qfq */
      SUBROUTINE MFQ2QFQ(DMFQ,DQFQ,IPRINT)
C*****************************************************************************
C
C
C     Input : DMFQ  - FQ matrix in molfdir basis
C
C     Output: DQFQ  - FQ matrix in quaternion basis
C
C
C     Written by J. Thyssen - Aug 22 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION DQFQ(NASHT,NORBT,NZ)
      DIMENSION DMFQ(2*NASHT,2*NORBT,2)
C
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output from MFQ2QFQ',-1)
         CALL HEADER('Input matrix:',-1)
         CALL PRQMAT(DMFQ,2*NASHT,2*NORBT,2*NASHT,2*NORBT,MIN0(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL DZERO(DQFQ, NASHT * NORBT * NZ )
      IOFF11 = 0
      IOFF12 = 0
      IOFF21 = 0
      IOFF22 = 0
      IOFF31 = 0
      IOFF32 = 0
      DO IFSYM = 1, NFSYM
         IOFF31 = IOFF31 + NASH(IFSYM)
         IOFF32 = IOFF32 + NORB(IFSYM)
         DO J = 1, NORB(IFSYM)
            CALL DCOPY(NASH(IFSYM), DMFQ(1+IOFF21,J+IOFF22,1), 1,
     $           DQFQ (1+IOFF11,J+IOFF12,1),1 )
         ENDDO
         IF (NZ.GE.2) THEN
            DO J = 1, NORB(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFQ(1+IOFF21,J+IOFF22,2), 1,
     $              DQFQ (1+IOFF11,J+IOFF12,2),1 )
            ENDDO
         END IF
         IF (NZ.GE.4) THEN
            DO J = 1, NORB(IFSYM)
               CALL DCOPY(NASH(IFSYM), DMFQ(1+IOFF21,J+IOFF32,1), 1,
     $              DQFQ (1+IOFF11,J+IOFF12,3), 1)
               CALL DCOPY(NASH(IFSYM), DMFQ(1+IOFF21,J+IOFF32,2), 1,
     $              DQFQ (1+IOFF11,J+IOFF12,4), 1)
            ENDDO
         ENDIF
         IOFF11 = IOFF11 + NASH(IFSYM)
         IOFF12 = IOFF12 + NORB(IFSYM)
         IOFF21 = IOFF21 + 2*NASH(IFSYM)
         IOFF22 = IOFF22 + 2*NORB(IFSYM)
         IOFF31 = IOFF31 + NASH(IFSYM)
         IOFF32 = IOFF32 + NORB(IFSYM)
      ENDDO
C
      IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('Output matrix:',-1)
         CALL PRQMAT(DQFQ,NASHT,NORBT,NASHT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_gendet */
      SUBROUTINE GOSCI_GENDET(NDET)
C*****************************************************************************
C
C     Driver routine.
C     Generate determinants for GOSCI.
C
C     Input : none
C
C     Output: NDET - number of determinants
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbopt.h>
#include <dgroup.h>
C
C     Set GOSCI common blocks
C
C     IRP  : fermion irrep
C     IRRP : symmetry of wave function 
C            1 .. NREP           : fermion irrep
C            NREP + 1 .. 2 * NREP: boson irrep
C
C             fermion               boson
C     C1   :  A                     a 
C     Ci   :  AG, AU                ag, au
C
C     C2   :  1E, 2E                a, b
C     Cs   :  1E, 2E                a, b
C     C2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
C     C2v  :  1E, 2E                a, b
C     D2h  :  1Eg, 2Eg, 1Eu, 2Eu    ag, bg, au, bu
C
      IF ( MOD ( NAELEC, 2 ) .EQ. 0 ) THEN
C        ...boson symmetry
         IRRP = IGOSCI_NREP() + IOPT_SYMMETRY
      ELSE
C        ...fermion symmetry
         IRRP = IOPT_SYMMETRY
      END IF
      write(6,*) 'gosci_gendet: nrep',IGOSCI_NREP()
      write(6,*) 'gosci_gendet: SYMM',IOPT_SYMMETRY
      write(6,*) 'gosci_gendet: IRRP',IRRP
C
      CALL GENDET(IDUMMY,IRRP)
C
      NDET = IGOSCI_NDET()
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_setone */
      SUBROUTINE GOSCI_SETONE(ONE_IN,NDIM)
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D00 )
C
#include "../dirrci/goscom.inc"
#include <dgroup.h>
C
      DIMENSION ONE_IN(NDIM,NDIM,*)
C
      DO J = 1,NDIM
         DO I = 1,NDIM
            ONER(I,J) = ONE_IN(I,J,1)
            IF ( NZ .GT. 1 ) THEN
               ONEI(I,J) = ONE_IN(I,J,2)
            ELSE
               ONEI(I,J) = D0
            END IF
         END DO
      END DO
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igosci_nrep */
      INTEGER FUNCTION IGOSCI_NREP()
#include <implicit.h>
#include "../dirrci/goscom.inc"
      IGOSCI_NREP = NREP
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igosci_ndet */
      INTEGER FUNCTION IGOSCI_NDET()
#include <implicit.h>
#include "../dirrci/goscom.inc"
      IGOSCI_NDET = NDET
      RETURN 
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_setcb */
      SUBROUTINE GOSCI_SETCB()
C*****************************************************************************
C
C     Driver routine.
C     Setup GOSCI common blocks.
C
C
C     Input : none
C
C     Output: none
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dcbopt.h>
#include <dgroup.h>
C
C     Set GOSCI common blocks
C
      write(6,*) 'ipropt ',ipropt
      CALL GOSCI_SETCB1(2*NASHT,NAELEC,IPROPT,NASH)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_setcb1 */
      SUBROUTINE GOSCI_SETCB1(IN_NORB,IN_NELEC,IN_IPRNT,NASH)
C*****************************************************************************
C
C     Setup GOSCI common blocks.
C
C     Input : IN_NORB, IN_NELEC, IN_PRNT
C             transferred from outside, since dcbopt.h is incompatible with
C             with goscom.inc
C
C     Output: 
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( D0 = 0.0D00 )
C
#include "../dirrci/goscom.inc"
C
#include <dgroup.h>
C
      DIMENSION MULTB_TMP(8,8), NASH(2)
C
C
C     /GODATIN/
C
C     - NORB  : number of active spinors (i.e., equal to 2 * NASHT)
C     - NELEC : number of active electrons
C     - IREF  : bitpattern for reference determinant 
C     - IREFE : frozen electrons
C               (used in conjuction with bitpattern in IREF)
C     - NREP  : number of fermion irreps
C     - IPRNT : print level
C
C
      NORB  = IN_NORB
      NELEC = IN_NELEC
      IPRNT = IN_IPRNT
      IREF  = 0
C     ...NREP is set below
      IREFE = 0
C
C     /GOBRINT/
C
C     - BREIT : Breit integrals
C
      BREIT = .FALSE.
C
C     /GODETS/
C
C     - IDET  : determinant bit strings
C               (set up by call to GENDET)
C
C     - NDET  : number of determinants
C               (set up by call to GENDET)
C
C     /GODPTAB/
C
C     - MULTB : symmetry multiplication table
C     - IRPAMO: irrep for each MO 
C
C     /GOCHR/
C
C     - REPNA : string with all symmetries (boson + fermion)
C
C     Set MULTB and REPNA
C
      CALL GMULTA (NZ,NFSYM,NREP,REPNA,MULTB_TMP)
      DO I = 1, 2 * NREP
         MULTB(I,0) = I
      END DO
      DO I = 1,8
         DO J = 1,8
            MULTB(I,J) = MULTB_TMP(I,J)
         END DO
      END DO
C
C     Set IRPAMO
C
      II    = 0
      IREPA = 0
      DO I = 1, NFSYM
         IREPA = IREPA + 1
         DO ISTR = 1, NASH(I)
            II = II + 1
            IRPAMO(II) = IREPA
         ENDDO
         IF (NZ.NE.4) IREPA = IREPA + 1
         DO ISTR = NASH(I) + 1, 2 * NASH(I)
            II = II + 1
            IRPAMO(II) = IREPA
         ENDDO
      ENDDO
C
C     /GOSBAS/
C
C     - EPS   : orbital energies (unused)
C     - OCC   : (unused)
C     - COPCOF: (unused)
C     - ECORE : core energy
C     - EAVER : average energy
C
      ECORE = D0
      EAVER = D0
C     
C     /GOSCIN/
C
C     - ONER  : core Hamiltonian (real part)
C     - ONEI  : core Hamiltonian (imaginary part)
C     
C     ONER and ONEI are setup by GOSCI_STARTVEC call
C
C     /GOSRESULT/
C
C     - D     :
C     - NDETI :
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rmakdm */
      SUBROUTINE GOSCI_RMAKDM(CL,CR,DV,PV,IKR)
C*****************************************************************************
C
C     Construct one- and two-electron active density matrix
C
C     D_{pq} = < CL | X^+_{pq} | CR >
C
C     P_{pqrs} = < CL | x^+_{pq,rs} | CR >
C
C     Input:
C        CL, CR  - CI vectors
C
C     Output:
C        DV      - one electron active density matrix
C        PV      - two electron active density matrix
C
C     WARNING: this code is not for the fainthearted CI expert, the code
C     is as un-optimized as possible... :-)
C
C     DV has dimension (2*NASHT)**2 * MIN0(NZ,2).
C     I haven't yet found the quaternion PV, so it has dimensions
C     (2*NASHT)**4 * MIN0(NZ,2). Triangular packing would be better,
C     but whatever... The code is now very straightforward.
C
C     Written by J. Thyssen - Oct 22 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00, DP5 = 0.5D00 )
C
#include "../dirrci/goscom.inc"
#include <dgroup.h>
C
      DIMENSION CL(NDET,NZ),CR(NDET,NZ),DV(NORB,NORB,2)
      DIMENSION PV(NORB,NORB,NORB,NORB,2)
      DIMENSION IKR(2,*)
      LOGICAL LBIT
C
C
#include <sdpre.h>
#include <ibtfun.h>
C
      CALL QENTER('GOSCI_RMAKDM')
C
C     Construct the IKR pointer
C     -------------------------
C
      CALL GOSCI_MKIKR(IKR,IPRNT)
C
C
C     Construct one-particle density matrix
C     -------------------------------------
C
C
      call header('real(CL)',-1)
      call output(CL(1,1),1,NDET,1,1,NDET,1,1,LUPRI)
      if ( nz .ge. 2 ) then
         call header('imag(CL)',-1)
         call output(CL(1,2),1,NDET,1,1,NDET,1,1,LUPRI)
      end if
      call header('real(CR)',-1)
      call output(CR(1,1),1,NDET,1,1,NDET,1,1,LUPRI)
      if ( nz .ge. 2 ) then
         call header('imag(CR)',-1)
         call output(CR(1,2),1,NDET,1,1,NDET,1,1,LUPRI)
      end if
      CALL DZERO(DV, NORB*NORB*MIN0(NZ,2) )
C
      DO 10 IU = 1, NORB
         IKU = IKR(1,IU)
         IKUSIGN = IKR(2,IU)
         DO 11 IV = 1, NORB
            IKV = IKR(1,IV)
            IKSIGN = IKUSIGN * IKR(2,IV)
c            write(6,*) '  V = ',IV
            DO 20 I = 1, NDET
               ID = IDET(I)
c               write(6,*) '    I = ',I
c               write(6,*) '    ID = ',ID
c               write(6,*) '    LBIT(ID,IU) = ',LBIT(ID,IU)
C
C              Act with u^+ on ID, i.e. remove orbital u in ID.
C              If orbital u not present in ID goto next determinant.
C
               IF ( .NOT. LBIT(ID,IU) ) GO TO 20
               IDD = IBTXOR ( ID, IBTSHL(1,IU-1) )
c               write(6,*) '    IDD = ',IDD
C
C              Find sign:
C              Find number of bits lower than IU.
C
               IDNBITS = NBITS ( IBTAND( ID, 2** ( IU - 1 ) - 1 ) )
c               write(6,*) '    IDNBITS = ',
c    $              IDNBITS, IBTAND(ID, 2**(IU-1) - 1 )
C
               IF ( MOD(IDNBITS,2) .EQ. 0) THEN
                  IDSIGN = -1
               ELSE
                  IDSIGN = 1
               END IF
c               write(6,*) '    IDSIGN = ',IDSIGN
c               write(6,*) '    CREF(I)',
c     $              (CREF(I,K), K = 1,MIN0(NZ,2))
C
               DO 21 J = 1, NDET
                  JD = IDET(J)
c                  write(6,*) '      J = ',J
c                  write(6,*) '      JD = ',JD
C     
                  IF ( .NOT. LBIT(JD,IV) ) GO TO 21
                  JDD = IBTXOR (JD, IBTSHL(1,IV-1) )
c                  write(6,*) '      bla = ',2 ** NORB - 1,
c     $                 IBTSHL(1,IV-1), 2**NORB - 1 - IBTSHL(1,IV-1)
c                  write(6,*) '      JDD = ',JDD
                  IOVLP = IBTXOR(IDD,JDD)
c                  write(6,*) '      XOR(IDD,JDD) = ',IOVLP
                  IF ( IOVLP .NE. 0 ) GO TO 21
                  JDNBITS = NBITS ( IBTAND( JD, 2** ( IV - 1 ) - 1 ) )
c                  write(6,*) '      JDNBITS = ',
c     $              JDNBITS, IBTAND(JD, 2**(IV-1) - 1 )
C
                  IF ( MOD(JDNBITS,2) .EQ. 0) THEN
                     JDSIGN = -1
                  ELSE
                     JDSIGN = 1
                  END IF
c                  write(6,*) '      JDSIGN = ',JDSIGN
c                  write(6,*) '      CREF(J)',
c     $                 (CREF(J,K), K = 1,MIN0(NZ,2))
C
                  IF ( NZ .EQ. 1 ) THEN
                     DV(IU,IV,1) = DV(IU,IV,1) +
     &                    CL(I,1) * CR(J,1) *
     $                    SDPRE(IDSIGN * JDSIGN)
                     DV(IKV,IKU,1) = DV(IKV,IKU,1) +
     &                    CL(I,1) * CR(J,1) *
     $                    SDPRE(IKSIGN * IDSIGN * JDSIGN)
                  ELSE
                     DV(IU,IV,1) = DV(IU,IV,1) + 
     $                    ( CL(I,1) * CR(J,1) +
     $                      CL(I,2) * CR(J,2) ) *
     $                    SDPRE(IDSIGN * JDSIGN)
                     DV(IU,IV,2) = DV(IU,IV,1) + 
     $                    ( - CL(I,2) * CR(J,1) 
     $                      + CL(I,1) * CR(J,2) ) *
     $                    SDPRE(IDSIGN * JDSIGN)
                     DV(IKV,IKU,1) = DV(IKV,IKU,1) + 
     $                    ( CL(I,1) * CR(J,1) +
     $                      CL(I,2) * CR(J,2) ) *
     $                    SDPRE(IKSIGN * IDSIGN * JDSIGN)
                     DV(IKV,IKU,2) = DV(IKV,IKU,1) + 
     $                    ( - CL(I,2) * CR(J,1) 
     $                      + CL(I,1) * CR(J,2) ) *
     $                    SDPRE(IKSIGN * IDSIGN * JDSIGN)
                  END IF
c                 write(6,*) '      DV(U,V) = ',
c    $                 (DV(IU,IV,K),K=1,MIN0(NZ,2))
 21            CONTINUE
 20         CONTINUE
 11      CONTINUE
 10   CONTINUE
      call header('real(DV) in Molfdir-basis',-1)
      call output(DV(1,1,1),1,NORB,1,NORB,NORB,NORB,1,LUPRI)
      if ( nz .ge. 2 ) then
         call header('imag(DV) in Molfdir-basis',-1)
         call output(DV(1,1,2),1,NORB,1,NORB,NORB,NORB,1,LUPRI)
      end if
C
C     Calculate PV
C
      CALL DZERO(PV, NORB ** 4 * MIN0(NZ,2) )
C
      DO 100 IU = 1, NORB
         IKU = IKR(1,IU)
         IKUSIGN = IKR(2,IU)
         DO 101 IV = 1, NORB
            IKV = IKR(1,IV)
            IKUVSIGN = IKUSIGN * IKR(2,IV)
            DO 102 IX = 1, NORB
               IKX = IKR(1,IX)
               IKXSIGN = IKR(2,IX)
               DO 103 IY = 1, NORB
                  IKY = IKR(1,IY)
                  IKXYSIGN = IKXSIGN * IKR(2,IY)
                  DO 200 I = 1, NDET
                     ID = IDET(I)
C
C     Act with u^+x^+ on <ID|, i.e. remove orbital u and x in |ID>. 
C     If orbital u or v not present in ID goto next determinant. 
C
                     IF ( .NOT. LBIT(ID,IU) ) GO TO 200
                     IDD = IBTXOR (ID, IBTSHL(1,IU-1) )
                     IF ( .NOT. LBIT(IDD,IX) ) GO TO 200
C
C     Find sign:
C     Find number of bits lower than IU.
C
                     IDNBITS = NBITS
     $                    ( IBTAND( ID, 2** ( IU - 1 ) - 1 ) )
C
                     IDSIGN = 1
                     IF ( MOD(IDNBITS,2) .EQ. 0) THEN
                        IDSIGN = -1
                     END IF
C
C     Find number of bits lower than IX.
C
                     IDNBITS = NBITS
     $                    ( IBTAND( IDD, 2** ( IX - 1 ) - 1 ) )
                     IDD = IBTXOR (IDD, IBTSHL(1,IX-1) )
C
                     IF ( MOD(IDNBITS,2) .EQ. 0) THEN
                        IDSIGN = -IDSIGN
                     END IF
C
                     DO 201 J = 1, NDET
                        JD = IDET(J)
C     
                        IF ( .NOT. LBIT(JD,IV) ) GO TO 201
                        JDD = IBTXOR (JD, IBTSHL(1,IV-1))
                        IF ( .NOT. LBIT(JDD,IY) ) GO TO 201
C
                        JDNBITS = NBITS
     $                       ( IBTAND( JD, 2** ( IV - 1 ) - 1) )
C
                        JDSIGN = 1
                        IF ( MOD(JDNBITS,2) .EQ. 0) THEN
                           JDSIGN = -1
                        END IF
C
                        JDNBITS = NBITS
     $                       ( IBTAND( JDD, 2** ( IY - 1 ) - 1) )
                        JDD = IBTXOR (JDD, IBTSHL(1,IY-1))
C
                        IF ( MOD(JDNBITS,2) .EQ. 0) THEN
                           JDSIGN = -JDSIGN
                        END IF
C
                        IOVLP = IBTXOR(IDD,JDD)
                        IF ( IOVLP .NE. 0 ) GO TO 201
C
                        IF ( NZ .EQ. 1 ) THEN
C
C                          Real CI expansion:
C
C                          C^*(I) * C(J) * signs =
C                          C(I) * C(J) * signs
C
                           PVELEM =  ( CL(I,1) * CR(J,1) ) *
     $                          SDPRE(IDSIGN * JDSIGN)
C
                           PV(IU,IV,IX,IY,1) = PV(IU,IV,IX,IY,1) +
     $                          PVELEM
C
                           PV(IKV,IKU,IX,IY,1) = PV(IKV,IKU,IX,IY,1) +
     $                          SDPRE(IKUVSIGN) * PVELEM
C
                           PV(IU,IV,IKY,IKX,1) = PV(IU,IV,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN) * PVELEM
C
                           PV(IKV,IKU,IKY,IKX,1) =
     $                          PV(IKV,IKU,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN * IKUVSIGN) * PVELEM
C
                        ELSE
C
C                          Complex CI expansion:
C     
C                          C^*(I) * C(J) * signs =
C                          (a-bi) * (c+di) * signs =
C                          ( (ac+bd) + (-bc+ad)i ) * signs
C
                           PVELER =
     $                          ( CL(I,1) * CR(J,1) +
     $                            CL(I,2) * CR(J,2) ) *
     $                          SDPRE(IDSIGN * JDSIGN)
                           PVELEI =
     $                          ( - CL(I,2) * CR(J,1)
     $                            + CL(I,1) * CR(J,2) ) *
     $                          SDPRE(IDSIGN * JDSIGN)
C
                           PV(IU,IV,IX,IY,1) = PV(IU,IV,IX,IY,1) +
     $                          PVELER
C
                           PV(IKV,IKU,IX,IY,1) = PV(IKV,IKU,IX,IY,1) +
     $                          SDPRE(IKUVSIGN) * PVELER
C
                           PV(IU,IV,IKY,IKX,1) = PV(IU,IV,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN) * PVELER
C
                           PV(IKV,IKU,IKY,IKX,1) =
     $                          PV(IKV,IKU,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN * IKUVSIGN) * PVELER
C
                           PV(IU,IV,IX,IY,1) = PV(IU,IV,IX,IY,1) +
     $                          PVELEI
C
                           PV(IKV,IKU,IX,IY,1) = PV(IKV,IKU,IX,IY,1) +
     $                          SDPRE(IKUVSIGN) * PVELEI
C
                           PV(IU,IV,IKY,IKX,1) = PV(IU,IV,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN) * PVELEI
C
                           PV(IKV,IKU,IKY,IKX,1) =
     $                          PV(IKV,IKU,IKY,IKX,1) +
     $                          SDPRE(IKXYSIGN * IKUVSIGN) * PVELEI
C
                        END IF
C                        
 201                 CONTINUE
 200              CONTINUE
 103           CONTINUE
 102        CONTINUE
 101     CONTINUE
 100  CONTINUE
      call header('PV Molfdir-basis',-1)
      WRITE(LUPRI,'(4X,A,14X,A,10X,A)')
     $     'Index','Real(PV)','Imag(PV)'
      do iu = 1,norb
         do iv = 1,norb
            do ix = 1,norb
               do iy = 1,norb
                  IF ( NZ .EQ. 1 ) THEN
                     IF ( ABS ( PV(IU,IV,IX,IY,1) ) .GT. D0 ) THEN
                        WRITE(LUPRI,'(4I3,2X,1P,D17.10)')
     $                       IU,IV,IX,IY,PV(IU,IV,IX,IY,1)
                     END IF
                  ELSE
                     if ( abs ( pv(iu,iv,ix,iy,1) + pv(iu,iv,ix,iy,2) )
     $                    .gt. d0 ) then
                        WRITE(LUPRI,'(4I3,2X,1P,D17.10,4X,D17.10)')
     $                       IU,IV,IX,IY,
     $                       PV(IU,IV,IX,IY,1),PV(IU,IV,IX,IY,2)
                     END IF
                  end if
               end do
            end do
         end do
      end do
C
C
C     Check that P{uv,xx} = (N-1) D_{uv}
C     
      DO IU = 1, NORB
         DO IV = 1, NORB
            TRACER = D0
            TRACEI = D0
            DO IX = 1, NORB
               TRACER = TRACER + PV(IU,IV,IX,IX,1)
               IF ( NZ .GE. 1 ) THEN
                  TRACEI = TRACEI + PV(IU,IV,IX,IX,2)
               END IF
            END DO
            WRITE(6,'(A,I2,A,I2,A,1P,D13.6,A,D13.6,A)')
     $           'real P(',IU,',',IV,',*,*) = ',DP5 * TRACER,
     $           ' (expected ',(NELEC - 1) * DV(IU,IV,1),')'
            IF ( NZ .GT. 1 ) THEN
               WRITE(6,'(A,I2,A,I2,A,1P,D13.6,A,D13.6,A)')
     $              'imag P(',IU,',',IV,',*,*) = ',DP5 * TRACEI,
     $              ' (expected ',(NELEC - 1) * DV(IU,IV,2),')'
            END IF
         END DO
      END DO
C
      CALL QEXIT('GOSCI_RMAKDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_mkikr */
      SUBROUTINE GOSCI_MKIKR(IKR,IPRINT)
C*****************************************************************************
C
C     For each orbital in molfdir basis find it's kramers partner
C     and determine whether it is a p or p-bar type orbital.
C
C     Input:
C
C     Output:
C        IKR - the Kramers pair pointer
C
C     Written by J. Thyssen - Aug 9 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION IKR(2,*)
C
#include <dgroup.h>
#include <dcborb.h>
C
C     The orbitals are listed as
C
C     Gerade orbitals                   Ungerade orbitals     
C     1 2 3 ... 1-bar 2-bar 3-bar ...   4 5 6 ... 4-bar 5-bar 6-bar
C
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
C
            IKR( 1, I + 2 * IASH(IFSYM) ) =
     $           I + 2 * IASH(IFSYM) + NASH(IFSYM)
            IKR( 2, I + 2 * IASH(IFSYM) ) = 1
C           the kramers partner
            IKR( 1, IKR(1, I + 2 * IASH(IFSYM) ) ) = I + 2 * IASH(IFSYM)
            IKR( 2, IKR(1, I + 2 * IASH(IFSYM) ) ) = -1
         END DO
      END DO
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('Output from GOSCI_MKIKR',-1)
         WRITE(LUPRI,*) 'Orbital   Kramers partner    Type'
         DO I = 1, 2*NASHT
            WRITE(LUPRI,9001)
     $           I,IKR(1,I),IKR(2,I)
         END DO
         CALL FLSHFO(LUPRI)
      END IF
 9001 FORMAT(2X,I3,13X,I3,13X,I2)
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rgeth2 */
      SUBROUTINE GOSCI_RGETH2(FQ,NDIM,H2AC,IDXAC,IPKR,PV)
C*****************************************************************************
C
C     Calculate FQ and get H2AC (two-electron integrals
C     with active indeces).
C
C     The integral reading code is a rip from RDTWO in cigosc.F
C
C     Input:
C
C     Output:
C        FQ   -   aux. FQ fock matrix
C        H2AC -   two-el. int. with act. ind.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D00, D0 = 0.0D00, DP5 = 0.5D00 )
C
#include "../dirrci/goscom.inc"
#include <dcbopt.h>
#include <dgroup.h>
C
      DIMENSION H2AC(NORB,NORB,NORB,NORB,2)
      DIMENSION PV(NORB,NORB,NORB,NORB,2)
      DIMENSION FQ(NORB,NDIM,NZ)
      DIMENSION INDK(N3),INDL(N3),BUF(N3,2)
      DIMENSION KR(-N1:N1)
      DIMENSION IDXAC(4,*), IPKR(2,*)
      CHARACTER DATEX*10, TIMEX*8
C
#include <sdpre.h>
C
      CALL QENTER('GOSCI_RGETH2')
C
      write(6,*) 'jtrlvl',jtrlvl
      IF ( JTRLVL .EQ. 0 ) THEN
C
         WRITE(LUPRI,'(/,A,/,A)')
     $        '*** WARNING from GOSCI_RGETH2 ***',
     $        'Only calculating H2AC.'
C
C        Read and distribute integrals.
C
         CALL DZERO(H2AC,NORB**4 * MIN0(NZ,2) )
         MDINT = 1
         OPEN(MDINT,FILE='MDCINT',FORM='UNFORMATTED')
         Read (mdint,err=10000,end=10000) datex,timex,nkr,
     &        (KR(I),KR(-I),I=1,NKR) 
C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nkrcoul = nkr
         nuniq = 0
 10      read (mdint,ERR=10010,END=10010) ikr,jkr,nint,
     &        (indk(inz),indl(inz),inz=1,nint),
     &        (buf(inz,1),buf(inz,2),inz=1,nint)
         If (ikr.EQ.0) GoTo 20  ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
         i = kr(ikr)
         itr = kr(-ikr)
         j = kr(jkr)
         jtr = kr(-jkr)
         SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
         do inz = 1, nint
            nuniq = nuniq + 1
            kkr = indk(inz)
            k = kr(kkr)
            ktr = kr(-kkr)
            lkr = indl(inz)
            l = kr(lkr)
            ltr = kr(-lkr)
            SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
            write(6,'(A,8I3)') 'i,j,k,l',I,J,K,L,ITR,JTR,KTR,LTR
            write(6,*) '   --> ',BUF(INZ,1)
            DO IZ = 1, MIN0(NZ,2) 
               H2AC(I,J,K,L,IZ) = BUF(INZ,IZ)
               H2AC(JTR,ITR,K,L,IZ) = BUF(INZ,IZ) * SDPRE( SignIJ )
               H2AC(I,J,LTR,KTR,IZ) = BUF(INZ,IZ) * SDPRE( SignKL )
               H2AC(JTR,ITR,LTR,KTR,IZ) =
     $              BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL )
C--> Complex conjugate plus time-reversed partners
               IZSIGN = 1
               IF ( IZ .EQ. 2 ) IZSIGN = -1
               H2AC(J,I,L,K,IZ) =
     $              BUF(INZ,IZ) * SDPRE( IZSIGN )
               H2AC(ITR,JTR,L,K,IZ) =
     $              BUF(INZ,IZ) * SDPRE( SignIJ * IZSIGN )
               H2AC(J,I,KTR,LTR,IZ) =
     $              BUF(INZ,IZ) * SDPRE( SignKL * IZSIGN )
               H2AC(ITR,JTR,KTR,LTR,IZ) =
     $              BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL * IZSIGN )
C--> Particle interchanged plus time-reversed partners
               H2AC(K,L,I,J,IZ) = 
     $              BUF(INZ,IZ)
               H2AC(LTR,KTR,I,J,IZ) =
     $              BUF(INZ,IZ) * SDPRE( SignKL )
               H2AC(K,L,JTR,ITR,IZ) =
     $              BUF(INZ,IZ) * SDPRE( SignIJ )
               H2AC(LTR,KTR,I,J,IZ) = 
     $              BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL )
C--   > Particle interchanged and complex conjugated plus time-reversed
c     partners
               H2AC(L,K,J,I,IZ) =
     $              BUF(INZ,IZ)
               H2AC(KTR,LTR,J,I,IZ) =
     $              BUF(INZ,IZ) * SDPRE ( SignKL * IZSIGN )
               H2AC(L,K,ITR,JTR,IZ) =
     $              BUF(INZ,IZ) * SDPRE ( SignIJ * IZSIGN )
               H2AC(KTR,LTR,ITR,JTR,IZ) =
     $              BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL * IZSIGN )
            Enddo
         Enddo
         Goto 10
 20      Continue
         CLOSE(MDINT)
         Write (6,1000) "Coulomb",datex,timex,nuniq
         call header('H2AC Molfdir-basis',-1)
         WRITE(LUPRI,'(4X,A,14X,A,10X,A)')
     $        'Index','Real(H2AC)','Imag(H2AC)'
         do iu = 1,norb
            do iv = 1,norb
               do ix = 1,norb
                  do iy = 1,norb
                     IF ( NZ .EQ. 1 ) THEN
                        IF ( ABS ( H2AC(IU,IV,IX,IY,1) ) .GT. D0 ) THEN
                           WRITE(LUPRI,'(4I3,2X,1P,D17.10)')
     $                          IU,IV,IX,IY,H2AC(IU,IV,IX,IY,1)
                        END IF
                     ELSE
                        if ( abs ( H2AC(iu,iv,ix,iy,1) +
     $                       H2AC(iu,iv,ix,iy,2) ).gt. d0 ) then
                           WRITE(LUPRI,'(4I3,2X,1P,D17.10,4X,D17.10)')
     $                          IU,IV,IX,IY,
     $                          H2AC(IU,IV,IX,IY,1),H2AC(IU,IV,IX,IY,2)
                        END IF
                     end if
                  end do
               end do
            end do
         end do
      ELSE
C
C     Construct the IKR pointer
C     -------------------------
C
         CALL GOSCI_MKIKR(IPKR,IPRNT)
C
C
C        The indeces read from the file are ordered as
C        Gerade orbitals                 Ungerade orbitals 
C        1 2 3 ... 1-bar 2-bar 3-bar ... 4 5 6 ... 6-bar 7-bar
C
C        From the level of the transformation we know the offset
C        for each index. For example, for level 1 transformations
C
C        Index 1: (general index)
C
C        1 .. NISH(1)  is inactive orbitals,
C        NISH(1)+1 .. NISH(1)+NASH(1)  is the active orbitals
C        NISH(1)+NASH(1)+1 .. NESH(1)  is (electronic) virtual orbitals
C
C        NESH(1)+1 .. NESH(1)+NISH(1)  is inactive orbitals (Kramers partners)
C        ...etc...
C
C        2 * NESH(1) + 1 .. 2 * NESH(1) + NISH(2)  inactive orbitals (irrep 2)
C        ...etc...
C
C        Index 2-4: (active indeces)
C
C        1..NASH(1)    is the active orbitals
C
C
         IF ( ITRLVL .GT. 0 ) THEN
            NGENERAL = ITRLVL
         ELSE
            NGENERAL = 4
         END IF
         CALL GOSCI_GTIDX(NGENERAL,IDXAC,IPROPT)
C
         CALL DZERO(H2AC,NORB**4 * MIN0(NZ,2) )
         CALL DZERO(FQ,  NORB * NDIM * MIN0(NZ,2) )
         MDINT = 1
         OPEN(MDINT,FILE='MDCINT',FORM='UNFORMATTED')
         Read (mdint,err=10000,end=10000) datex,timex,nkr,
     &        (KR(I),KR(-I),I=1,NKR) 
         write(6,*) 'nkr = ',nkr
         write(6,*) 'kr  + ',(kr(i),i=1,nkr)
         write(6,*) 'kr  - ',(kr(-i),i=1,nkr)
C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nkrcoul = nkr
         nuniq = 0
 100     read (mdint,ERR=10010,END=10010) ikr,jkr,nint,
     &        (indk(inz),indl(inz),inz=1,nint),
     &        (buf(inz,1),buf(inz,2),inz=1,nint)
         If (ikr.EQ.0) GoTo 200  ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
         i = kr(ikr)
         itr = kr(-ikr)
         j = kr(jkr)
         jtr = kr(-jkr)
         SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
         do inz = 1, nint
            nuniq = nuniq + 1
            kkr = indk(inz)
            k = kr(kkr)
            ktr = kr(-kkr)
            lkr = indl(inz)
            l = kr(lkr)
            ltr = kr(-lkr)
            SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C
C           If all 4-indeces are active integrals go into
C           H2AC.
C
            IAC   = IDXAC(1,I)
            ITRAC = IDXAC(1,ITR)
            JAC   = IDXAC(2,J)
            JTRAC = IDXAC(2,JTR)
            KAC   = IDXAC(3,K)
            KTRAC = IDXAC(3,KTR)
            LAC   = IDXAC(4,L)
            LTRAC = IDXAC(4,LTR)
C
            write(6,'(A,8I3)') 'i,j,k,l',I,J,K,L,ITR,JTR,KTR,LTR
            write(6,'(A,8I3)') '   --> i,j,k,l',
     $           IAC,JAC,KAC,LAC,ITRAC,JTRAC,KTRAC,LTRAC
            write(6,*) '   --> ',BUF(INZ,1)
#ifdef UNDEF
            if ( itrac .lt. 0 .or. itrac .gt. 16 ) then
               do iii = 1,4
                  do jjj = 1,2*nkr
                     write(6,*) 'idxac =',iii,jjj,idxac(iii,jjj)
                  end do
               end do
               call quit('yrk!')
            end if
#endif
            IF ( IAC .NE. 0 .AND. JAC .NE. 0 .AND.
     $           KAC .NE. 0 .AND. LAC .NE. 0 ) THEN
C
C--> Original integral plus time-reversed partners
               DO IZ = 1, MIN0(NZ,2) 
                  H2AC(IAC,JAC,KAC,LAC,IZ) =
     $                 BUF(INZ,IZ)
                  H2AC(JTRAC,ITRAC,KAC,LAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignIJ )
                  H2AC(IAC,JAC,LTRAC,KTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignKL )
                  H2AC(JTRAC,ITRAC,LTRAC,KTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL )
C--> Complex conjugate plus time-reversed partners
                  IZSIGN = 1
                  IF ( IZ .EQ. 2 ) IZSIGN = -1
                  H2AC(JAC,IAC,LAC,KAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( IZSIGN )
                  H2AC(ITRAC,JTRAC,LAC,KAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignIJ * IZSIGN )
                  H2AC(JAC,IAC,KTRAC,LTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignKL * IZSIGN )
                  H2AC(ITRAC,JTRAC,KTRAC,LTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL * IZSIGN )
C--> Particle interchanged plus time-reversed partners
                  H2AC(KAC,LAC,IAC,JAC,IZ) = 
     $                 BUF(INZ,IZ)
                  H2AC(LTRAC,KTRAC,IAC,JAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignKL )
                  H2AC(KAC,LAC,JTRAC,ITRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE( SignIJ )
                  H2AC(LTRAC,KTRAC,IAC,JAC,IZ) = 
     $                 BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL )
C--   > Particle interchanged and complex conjugated plus time-reversed
c     partners
                  H2AC(LAC,KAC,JAC,IAC,IZ) =
     $                 BUF(INZ,IZ)
                  H2AC(KTRAC,LTRAC,JAC,IAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE ( SignKL * IZSIGN )
                  H2AC(LAC,KAC,ITRAC,JTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE ( SignIJ * IZSIGN )
                  H2AC(KTRAC,LTRAC,ITRAC,JTRAC,IZ) =
     $                 BUF(INZ,IZ) * SDPRE ( SignIJ * SignKL * IZSIGN )
               Enddo
            END IF
C
C           Construct FQ in Molfdir basis:
C           At least one index must be general, and
C           J, K, L must be active indeces.
C 
            IF ( NGENERAL .GT. 0 .AND.
     $           JAC .NE. 0 .AND. KAC .NE. 0 .AND. LAC .NE. 0) THEN
C
C              FQ(U,P) = sum_{VXY} P_{UV,XY} (PV|XY)
C              Ie. some integral (PV|XY) contributes all
C              FQ(U,P) where U is any active index.
C
C
               DO IU = 1, NORB
                  IF ( NZ .EQ. 1 ) THEN
                     write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                    'before fq(',IU,',',I,') = ',FQ(IU,I,1)
                     if ( iac .ne. 0 ) then
                        write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                       'before fq(',IU,',',JTR,') = ',FQ(IU,JTR,1)
                        write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                       'before fq(',IU,',',J,') = ',FQ(IU,J,1)
                     end if
                     write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                    'before fq(',IU,',',ITR,') = ',FQ(IU,ITR,1)
C
C                    Kramers conjungation and/or complex conj.
C
C                    the integral: (pv|xy)
C                    FQ(u,p) = FQ(u,p) + P(uv,xy) (pv|xy)
                     FQ(IU,I,1) = FQ(IU,I,1) +
     $                    PV(IU,JAC,KAC,LAC,1) * BUF(INZ,1)
                     IF ( IAC .NE. 0 ) THEN
C                       the integral: (vbar pbar|xy)
C                       FQ(u,vbar) = FQ(u,vbar) + P(u pbar,xy) (vbar pbar|xy)
                        FQ(IU,JTR,1) = FQ(IU,JTR,1) +
     $                       PV(IU,ITRAC,KAC,LAC,1) * BUF(INZ,1)
                     END IF
                     IF ( I .NE. J .OR. K .NE. L ) THEN
C                       the integral: (vp|yx)^*
C                       FQ(u,v) = FQ(u,v) + P(up,yx) (vp|yx)^*
                        IF ( IAC .NE. 0 ) THEN
                           FQ(IU,J,1) = FQ(IU,J,1) +
     $                          PV(IU,IAC,LAC,KAC,1) * BUF(INZ,1)
                        END IF
C                       the integral: (pbar vbar|yx)
C                       FQ(u,pbar) = FQ(u,pbar) + P(u vbar,yx) (pbar vbar|yx)^*
                        FQ(IU,ITR,1) = FQ(IU,ITR,1) +
     $                       PV(IU,JTRAC,LAC,KAC,1) * BUF(INZ,1)
                     END IF
#ifdef UNDEF
                     write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                    'after fq(',IU,',',I,') = ',FQ(IU,I,1)
                     write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                    'after fq(',IU,',',JTR,') = ',FQ(IU,JTR,1)
                     if ( jac .ne. 0 ) then
                        write(6,'(A,I1,A,I1,A,1P,D17.10)')
     $                       'after fq(',IU,',',J,') = ',FQ(IU,J,1)
                     end if
                     write(6,'(A,I1,A,I1,A,1P,D17.10,/)')
     $                    'after fq(',IU,',',ITR,') = ',FQ(IU,ITR,1)
#endif
                  ELSE
                     CALL QUIT('NZ > 1 FQ not implemented...')
                     FQ(IU,I,1) = FQ(IU,I,1) +
     $                    PV(IU,JAC,KAC,LAC,1) * BUF(INZ,1) -
     $                    PV(IU,JAC,KAC,LAC,2) * BUF(INZ,2)
                     FQ(IU,I,1) = FQ(IU,I,1) +
     $                    PV(JTRAC,IKU,KAC,LAC,1) * BUF(INZ,1) -
     $                    PV(JTRAC,IKU,KAC,LAC,2) * BUF(INZ,2)
                     FQ(IU,I,1) = FQ(IU,I,1) +
     $                    PV(IU,JAC,LTRAC,KTRAC,1) * BUF(INZ,1) -
     $                    PV(IU,JAC,LTRAC,KTRAC,2) * BUF(INZ,2)
                     FQ(IU,I,1) = FQ(IU,I,1) + 
     $                    PV(JTRAC,IKU,LTRAC,KTRAC,1) * BUF(INZ,1) -
     $                    PV(JTRAC,IKU,LTRAC,KTRAC,2) * BUF(INZ,2)
                     FQ(IU,I,2) = FQ(IU,I,2) +
     $                    PV(IU,JAC,KAC,LAC,1) * BUF(INZ,2) +
     $                    PV(IU,JAC,KAC,LAC,2) * BUF(INZ,1)
                     FQ(IU,I,2) = FQ(IU,I,2) +
     $                    PV(JTRAC,IKU,KAC,LAC,1) * BUF(INZ,2) +
     $                    PV(JTRAC,IKU,KAC,LAC,2) * BUF(INZ,1)
                     FQ(IU,I,2) = FQ(IU,I,2) +
     $                    PV(IU,JAC,LTRAC,KTRAC,1) * BUF(INZ,2) +
     $                    PV(IU,JAC,LTRAC,KTRAC,2) * BUF(INZ,1)
                     FQ(IU,I,2) = FQ(IU,I,2) + 
     $                    PV(JTRAC,IKU,LTRAC,KTRAC,1) * BUF(INZ,2) +
     $                    PV(JTRAC,IKU,LTRAC,KTRAC,2) * BUF(INZ,1)
                     
                  END IF
               END DO
C     
            END IF
C     

         Enddo
         Goto 100
 200     Continue
         close (mdcint)
         Write (6,1000) "Coulomb",datex,timex,nuniq
         call header('H2AC Molfdir-basis',-1)
         WRITE(LUPRI,'(4X,A,14X,A,10X,A)')
     $        'Index','Real(H2AC)','Imag(H2AC)'
         do iu = 1,norb
            do iv = 1,norb
               do ix = 1,norb
                  do iy = 1,norb
                     IF ( NZ .EQ. 1 ) THEN
                        IF ( ABS ( H2AC(IU,IV,IX,IY,1) ) .GT. D0 ) THEN
                           WRITE(LUPRI,'(4I3,2X,1P,D17.10)')
     $                          IU,IV,IX,IY,H2AC(IU,IV,IX,IY,1)
                        END IF
                     ELSE
                        if ( abs ( H2AC(iu,iv,ix,iy,1) +
     $                       H2AC(iu,iv,ix,iy,2) ).gt. d0 ) then
                           WRITE(LUPRI,'(4I3,2X,1P,D17.10,4X,D17.10)')
     $                          IU,IV,IX,IY,
     $                          H2AC(IU,IV,IX,IY,1),H2AC(IU,IV,IX,IY,2)
                        END IF
                     end if
                  end do
               end do
            end do
         end do
         call header('fq in molfdir basis',-1)
         call prqmat(FQ,NORB,NKR*2,NORB,NKR*2,MIN0(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL QEXIT('GOSCI_RGETH2')
      RETURN
 1000 Format (/1X,A7,' integral file was generated at ',A10,1x,A8,
     &     /' Read',I8,' unique integrals')
10000 STOP 'END READING MDCINT'
10010 STOP 'ERROR READING MDCINT'
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_gtidx */
      SUBROUTINE GOSCI_GTIDX(NGENERAL,IOFF,IPRINT)
C*****************************************************************************
C
C     Calculate index vector from general index onto active index.
C
C     Input:
C        NGENERAL - number of general indeces
C
C     Output:
C        IOFF - offsets
C
C     Written by J. Thyssen - Aug 9 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION IOFF(4,3,*)
      DIMENSION IOFFP(2), IOFFT(2)
      CHARACTER CBARRED(2)*8
C
      DATA CBARRED /'unbarred','barred'/
C
#include <dgroup.h>
#include <dcborb.h>
#include <dcbopt.h>
C
C     The orbitals are listed as
C
C     Gerade orbitals                   Ungerade orbitals     
C     1 2 3 ... 1-bar 2-bar 3-bar ...   4 5 6 ... 4-bar 5-bar 6-bar
C
C
C     General indeces
C
      IF ( OPT_NOPFQ ) THEN
         IOFFP (1) = 0
         IOFFP (2) = 0
         IOFFT (1) = NESH(1)
         IOFFT (2) = NESH(2)
         IOFFTT = NESHT
      ELSE
         IOFFP (1) = NPSH(1)
         IOFFP (2) = NPSH(2)
         IOFFT (1) = NORB(1)
         IOFFT (2) = NORB(2)
         IOFFTT = NORBT
      END IF
      IOFF1 = 0
      IOFFJ = 0
C
      CALL IZERO ( IOFF, 16 * IOFFTT )
C
      DO IFSYM = 1,NFSYM
C        set p:
         DO J = 1, IOFFT(IFSYM)
            DO I = 1, NGENERAL
               IOFF(I,2,J+IOFFJ) = 1
               IOFF(I,3,J+IOFFJ) = J
            END DO
         END DO
C        set pbar:
         DO J = IOFFT(IFSYM) + 1, IOFFT(IFSYM) * 2
            DO I = 1, NGENERAL
               IOFF(I,2,J+IOFFJ) = -1
               IOFF(I,3,J+IOFFJ) = J - IOFFT(IFSYM)
            END DO
         END DO
C        active: u
         DO J = IOFFP(IFSYM) + NISH(IFSYM)+1, IOFFP(IFSYM) + NOCC(IFSYM)
            IOFF1 = IOFF1 + 1
            DO I = 1, NGENERAL
               IOFF(I,1,J + IOFFJ) = IOFF1
               IOFF(I,2,J + IOFFJ) = 1
            END DO
         END DO
C        active: ubar
         DO J = IOFFT(IFSYM) + IOFFP(IFSYM) + NISH(IFSYM) + 1,
     $        IOFFT(IFSYM) + IOFFP(IFSYM) + NOCC(IFSYM)
            IOFF1 = IOFF1 + 1
            DO I = 1, NGENERAL
               IOFF(I,1,J + IOFFJ) = IOFF1
               IOFF(I,2,J+ IOFFJ) = -1
            END DO
         END DO
         IOFFJ = IOFFJ + 2 * IOFFT(IFSYM)
      END DO
C
C     Active indeces:
C     - the offset is trivially 0
C
      IOFFJ = 0
      DO IFSYM = 1, NFSYM
         DO J = 1, NASH(IFSYM)
            DO I = NGENERAL + 1, 4
               IOFF(I,1,J + IOFFJ) = J + IOFFJ
               IOFF(I,2,J + IOFFJ) = 1
            END DO
         END DO
         DO J = NASH(IFSYM) + 1, NASH(IFSYM) * 2
            DO I = NGENERAL + 1, 4
               IOFF(I,1,J + IOFFJ) = J + IOFFJ
               IOFF(I,2,J + IOFFJ) = -1
            END DO
         END DO
         IOFFJ = IOFFJ + 2 * NASH(IFSYM)
      END DO
C
      IF ( IPRINT .GE. 10 .AND. NGENERAL .GT. 0) THEN
         CALL HEADER('Output from GOSCI_GTOFF',-1)
         DO J = 1 , 4
            WRITE(LUPRI,'(A,I3,A)')
     $           'Orbital   corresponding active orbital (index ',J,')'
            DO I = 1, 2 * IOFFTT
               IF ( IOFF(J,1,I) .GT. 0 ) THEN
                  WRITE(LUPRI,9001) I,IOFF(J,1,I),
     &                 CBARRED((3-IOFF(J,2,I))/2), IOFF(J,3,I)
               ELSE
                  WRITE(LUPRI,9002) I,'none',
     &                 CBARRED((3-IOFF(J,2,I))/2), IOFF(J,3,I)
               END IF
            END DO
         END DO
         CALL FLSHFO(LUPRI)
      END IF
 9001 FORMAT(2X,I3,13X,I3,10X,A10,5X,I3)
 9002 FORMAT(2X,I3,13X,A4,9X,A10,5X,I3)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_makeh */
      SUBROUTINE GOSCI_MAKEH(HR,HI,VMUUUU)
C*****************************************************************************
C
C     Construct Hamiltonian for use in gosci.
C
C     Input:
C        VMUUUU - two electron integrals.
C
C     Output:
C        HR - real part of Hamiltonian
C        HI - imag part of Hamiltonian
C
C     Written by J. Thyssen - Nov 6 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION HR(NDET,NDET),HI(NDET,NDET)
      DIMENSION VMUUUU(NORB,NORB,NORB,NORB,*)
#include "../dirrci/goscom.inc"
C
      DO I = 1, NDET
        DO J = 1, I
          CALL GOSCI_HMAT (IDET(I),IDET(J),HR(I,J),HI(I,J),VMUUUU)
          HR (J,I) =  HR (I, J)
          HI (J,I) = -HI (I, J)
        ENDDO
      ENDDO
      PRINT*, ' ** HR **'
      DO 100 I=1,NDET
        WRITE(6,1010) (HR(I,J),J=1,NDET)
100     CONTINUE
      PRINT*, ' ** HI **'
      DO 110 I=1,NDET
        WRITE(6,1010) (HI(I,J),J=1,NDET)
110   CONTINUE
CTROND
 1010 FORMAT(15F8.4)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_hmat */
      SUBROUTINE GOSCI_HMAT(ID,JD,HMR,HMI,VMUUUU)
C*****************************************************************************
C
C     Construct Hamiltonian element I,J.
C
C     Input:
C        ID, JD: bit strings for determinant I and J.
C        HMR, HMI: Hamiltonian matrix elements
C        VMUUUU - two electron integrals.
C
C     Output:
C        HMR, HMI: Hamiltonian matrix elements
C
C     Written by J. Thyssen - Nov 6 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION VMUUUU(NORB,NORB,NORB,NORB,*)
#include "../dirrci/goscom.inc"
C
#include <ibtfun.h>
C
      HMR = 0.0D0
      HMI = 0.0D0
C
      IDV = IBTXOR (ID, JD)
      IT = NBITS (IDV)
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
         DO 30 I = 1, NORB
            IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
               HMR = HMR + ONER (I, I)
               DO 20 J = I + 1, NORB
                  IF (IBTAND (ID, IBTSHL (1, J - 1)) .NE. 0) THEN
C---------------------
C                    H = H + (I, I, J, J)
C---------------------
                     HMR = HMR + VMUUUU(I,I,J,J,1)
                     HMI = HMI + VMUUUU(I,I,J,J,2)
C---------------------
C                    H = H - (I, J, J, I)
C---------------------
                     HMR = HMR - VMUUUU(I,J,J,I,1)
                     HMI = HMI - VMUUUU(I,J,J,I,2)
                  ENDIF
 20            CONTINUE
            END IF
 30      CONTINUE
         RETURN
      END IF
C
      IF(IT.EQ.2) THEN
C
C        ===================================
C        DETERMINANTS DIFFER BY ONE FUNCTION
C        ===================================
C
         IA = IBTAND (IDV, ID)
         IB = IBTAND (IDV, JD)
         DO 35 I = 1, NORB
            IF (IBTAND (IA, IBTSHL (1, I - 1)) .NE. 0) K = I
            IF (IBTAND (IB, IBTSHL (1, I - 1)) .NE. 0) L = I
 35      CONTINUE
C        
         HMR = HMR + ONER (K, L)
         HMI = HMI + ONEI (K, L)
C
         DO 40 I = 1, NORB
            IF ((I .EQ. K) .OR. (I .EQ. L)) GO TO 40
            IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
C---------------------
C              H = H + (I, I, K, L)
C---------------------
               HMR = HMR + VMUUUU(I,I,K,L,1)
               HMI = HMI + VMUUUU(I,I,K,L,2)
C---------------------
C              H = H - (I, L, K, I)
C---------------------
               HMR = HMR - VMUUUU(I,L,K,I,1)
               HMI = HMI - VMUUUU(I,L,K,I,2)
            END IF
 40      CONTINUE
         PH = PHASE (K, ID) * PHASE (L, JD)
         HMR = PH * HMR
         HMI = PH * HMI
         RETURN
      END IF
C
      IF(IT.EQ.4) THEN
C
C        ====================================
C        DETERMINANTS DIFFER BY TWO FUNCTIONS
C        ====================================
C        
         IA = IBTAND (IDV, ID)
         IB = IBTAND (IDV, JD)
         L = 0
         N = 0
         DO 50 I = 1, NORB
            IF (IBTAND (IA, IBTSHL (1, I - 1)) .NE. 0) THEN
               K = L
               L = I
            END IF
            IF (IBTAND (IB, IBTSHL (1, I - 1)) .NE. 0) THEN
               M = N
               N = I
            END IF
 50      CONTINUE
C
C---------------------
C        H = H + (K, M, L, N)
C---------------------
         HMR = HMR + VMUUUU(K,M,L,N,1)
         HMI = HMI + VMUUUU(K,M,L,N,2)
C---------------------
C        H = H - (K, N, L, M)
C---------------------
         HMR = HMR - VMUUUU(K,N,L,M,1)
         HMI = HMI - VMUUUU(K,N,L,M,2)
         PH = PHASE (K, ID) * PHASE (L, ID) *
     +        PHASE (M, JD) * PHASE (N, JD)
         HMR = HMR * PH
         HMI = HMI * PH
         RETURN
      END IF
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C     
      HMR = 0.0D0
      HMI = 0.0D0
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_sigma */
      SUBROUTINE GOSCI_SIGMA(BCVEC,SCVEC,IBCVC,FCAC,H2AC,
     &     IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Calculate conf. contribution to conf. sigma vector.
C     (Equation 3.29 in the KRMC paper (JCP 104, 4083 (1996)).
C     
C     Input:
C       BCVEC: the CI trial vectors.
C       IBCVC: boson symmetry of the CI trial vectors.
C       GPCI : the CI gradient (only needed for MC wave func. optimization).
C       FCAC : active-active part of FC fock matrix 
C       H2AC : two-electron integrals.
C
C     Output:
C       SCVEC: the CI sigma vectors.
C
C     Written by J. Thyssen - Nov 23 2000
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( D1 = 1.0D00, DM1 = -D1, D0 = 0.0D00, D2 = 2.0D00 )
C
      DIMENSION BCVEC(NZCONFQ,*)
      DIMENSION SCVEC(NZCONFQ,*)
      DIMENSION WORK(*),IBCVC(*)
      DIMENSION FCAC(*),H2AC(*)
C
#include <dcbxrs.h>
#include <dcborb.h>
#include <dgroup.h>
C
#include <memint.h>
C
      CALL QENTER('GOSCI_SIGMA')
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('Output from GOSCI_SIGMA:',-1)
      END IF
C
C
C     *********************************
C     *** Calculate CI sigma vector ***
C     *********************************
C
C     In princinple we would do this direct, i.e. calculating a
C     standard CI sigma vector. Instead we do it the ``easy'' way:
C
C     - calculate CI matrix
C     - do a matrix x vector multiplication
C
C     Calculate CI matrix
C     ===================
C
C
C     Transform quarternion FCAC (Dirac MO-basis) to Molfdir MO-basis
C     ---------------------------------------------------------------
C
      CALL MEMGET('REAL',KMOLFFC,(2*NASHT)*(2*NASHT)*MIN0(NZ,2),
     $     WORK,KFREE,LFREE)
      CALL QFC2MFC(FCAC,WORK(KMOLFFC),IPRINT)
      CALL GOSCI_SETONE(WORK(KMOLFFC),2*NASHT)
      CALL MEMREL('after RGETFC',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Transform integrals from Dirac (NZ,3) basis to Molfdir basis
C     ------------------------------------------------------------
C
      CALL MEMGET('REAL',KMUUUU,(2*NASHT)**4 * 2, WORK,KFREE,LFREE)
      CALL DZERO(WORK(KMUUUU),(2*NASHT)**4 * 2)
C
      CALL DNZ32M(H2AC,WORK(KMUUUU),IPRINT)
C
C     Construct Hamiltonian:
C     ----------------------
C
      NZCONF2 = NZCONF ** 2
      CALL MEMGET('REAL',KHR,NZCONF2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHI,NZCONF2,WORK,KFREE,LFREE)
C
      CALL DZERO(WORK(KHR),NZCONF2)
      CALL DZERO(WORK(KHI),NZCONF2)
      CALL GOSCI_MAKEH(WORK(KHR),WORK(KHI), WORK(KMUUUU))
C
C     
C     *******************************
C     *** Loop over trial vectors ***
C     *******************************
C
C
      DO I  = 1, NCSIM
C
C
C        **********************
C        *** <\mu|H|B> part ***
C        **********************
C
C        Real part: real(H) * real(BCVEC) - imag(H) * imag(BCVEC)
C
C        real(H) * real(BCVEC):
C
         CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &        WORK(KHR),NZCONF,
     &        BCVEC(1,I),NZCONF,
     &        D0,SCVEC(1,I),NZCONF)
C
         IF (NZ .GE. 2) THEN
C           
C           imag(H) * imag(BCVEC):
C
            CALL DGEMM('N','N',NZCONF,1,NZCONF,DM1,
     &           WORK(KHI),NZCONF,
     &           BCVEC(1+NZCONF,I),NZCONF,
     &           D1,SCVEC(1,I),NZCONF)
C
C
C           Imaginary part: real(H) * imag(BCVEC) + imag(H) * real(BCVEC)
C
C           real(H) * imag(BCVEC):
C
            CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &           WORK(KHR),NZCONF,
     &           BCVEC(1+NZCONF,I),NZCONF,
     &           D0,SCVEC(1+NZCONF,I),NZCONF)
C
C           imag(H) * real(BCVEC):
C
            CALL DGEMM('N','N',NZCONF,1,NZCONF,D1,
     &           WORK(KHI),NZCONF,
     &           BCVEC(1,1),NZCONF,
     &           D1,SCVEC(1+NZCONF,I),NZCONF)
C
         END IF
C
         IF ( IPRINT .GE. 0 ) THEN
            WRITE(LUPRI,'(1X,A,I2,A)') 'H * BCVEC(',I,'), real part'
            CALL OUTPUT(SCVEC(1,I),1,NZCONF,1,1,NZCONF,1,1,LUPRI)
            IF (NZ .GE. 2) THEN
               WRITE(LUPRI,'(1X,A,I2,A)') 'H * BCVEC(',I,'), imag part'
               CALL OUTPUT(SCVEC(1+NZCONF,I),1,NZCONF,1,1,
     &              NZCONF,1,1,LUPRI)
            END IF
         END IF
C
      END DO
C
      CALL MEMREL('GOSCI_SIGMA',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT ('GOSCI_SIGMA')
      RETURN
      END
#ifdef UNDEF
Description of the molfdir MRCONEE file

NORB, BREIT, ECORE

- NORB   : the number of active orbitals
- BREIT  : Breit interaction
- ECORE  : core energy

NSYMRP, (REPN(IRP), IRP = 1, NSYMRP )

- NSYMRP : number of fermion representations
- REPN   : string with name of fermion symmetry

NREP, (REPNA(IRP), IRP = 1, 2 * NREP )

- NREP   : number of irreps (number of boson + fermion irreps = 2 * NREP)
- REPNA  : string with name of symmetry
           (first is fermion irreps, then boson)

( ( MULTB(I,J), J = 1, 2 * NREP ), I = 1, 2 * NREP )

- MULTB  : multiplication table

( IRPMO(IMO), IRPAMO(IMO), EPS(IMO), IMO = 1,NORB)

- IRPMO  : fermion irreps for each MO
- IRPAMO : irrep (fermion or boson) irreps for each MO
- EPS    : spinor energy

( ( ONER(IMO,JMO), ONEI(IMO,JMO), JMO = 1, NORB ), IMO = 1, NORB )

- ONER   : real part of core Hamiltonian
- ONEI   : imaginary part of core Hamiltonian

#endif
