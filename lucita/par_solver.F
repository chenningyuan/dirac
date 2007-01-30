C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
#if defined (VAR_MPI)
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAS_PAR(NAEL,IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,H,XA,XB,SCR,RJ,RK,
     &           NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &           IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
* Parallel adaption in January 2007, Stefan Knecht
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK),NPARBLOCK(NBLOCK)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMAT(H,1,NORB,1,NORB)
        WRITE(6,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' I12 = ',I12
      END IF
*
*  Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
*              +         J(I,J) * NIA*NJB
*
*. K goes to J - K
      IF(I12.EQ.2)
     &CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
      IDET = 0
      ITDET = 0
      IF(LUDIA.NE.0) CALL REWINO(LUDIA)
*
      DO 100 IBLK = 1, NBLOCK
*
*. write only blocks to disc that are necessary for this node
        IF(NPARBLOCK(IBLK).NE.MYPROC) GOTO 100
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
        IOFF =  1
        DO IA = 1, NSSOA(IASM,IATP)
          EAA = 0.0D0
          DO IEL = 1, NAEL
            IAEL = IASTR(IEL,IA)
            EAA = EAA + H(IAEL)
            IF(I12.EQ.2) THEN
              DO JEL = 1, NAEL
                EAA =   EAA + 0.5D0*RK(IASTR(JEL,IA),IAEL )
              END DO
            END IF
          END DO
          RJKAA(IA-IOFF+1) = EAA
        END DO
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
          IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)
*
          DO IEL = 1, NBEL
            IBEL = IBSTR(IEL,IB)
            HB = HB + H(IBEL )
*
            IF(I12.EQ.2) THEN
              DO JEL = 1, NBEL
                RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
              END DO
*
              DO IORB = 1, NORB
                XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
              END DO
            END IF
          END DO
          EB = HB + 0.5D0*RJBB + ECORE
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
*
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA))
            END DO
            DIAG(IDET) = X
            IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
          END DO
*         ^ End of loop over alpha strings|
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            if(myproc.eq.master)then
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMAT(DIAG,1,IDET,1,IDET)
            endif
          END IF
          CALL ITODS(IDET,1,-1,LUDIA)
          CALL TODSC(DIAG,IDET,-1,LUDIA)
          IDET = 0
        END IF
 100  CONTINUE
*        ^ End of loop over blocks

      IF(NTEST.GE.5) WRITE(6,*)
     &' Number of diagonal elements generated ',ITDET
CSK      WRITE(6,*) ' Number of diagonal elements generated by node',ITDET,myproc
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        WRITE(6,*) ' CIDIAGONAL '
        CALL WRTMAT(DIAG(1),1,IDET,1,IDET)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAT_PAR(DIAG,LUDIA,ECORE,ICISTR,I12,
     &                      IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK),NPARBLOCK(NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('CIDIA')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' GASDIAT_PAR speaking'
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
      END IF
*
**. Local memory
*
      IDUM = 0
      CALL MEMMAN(IDUM,  IDUM,    'MARK  ',IDUM,'GASDIA')
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(WORK(KLH1D))
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
      IF( LUDIA .GT. 0 ) CALL REWINO(LUDIA)
      CALL GASDIAS_PAR(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &     NACOB,DIAG,NSMST,WORK(KLH1D),
     &     WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &     WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &     WORK(KLRJKA),I12,IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*.Flush local memory
      CALL MEMMAN(IDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('CIDIA')
*
C?    Call Abend2( ' Stefan forced me to stop after GASDIA ' )
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GATHER_LOW_PAR(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                          RECVARRAY,TESTARR,
     &                          IRECVARRAY,ITESTARRAY,IRANKARR,
     &                          IGATCOMM,NTESTG)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      DIMENSION SUBVAL(*),ISCAT(*),RECVARRAY(*),TESTARR(*)
      DIMENSION IRECVARRAY(*),ITESTARRAY(*),IRANKARR(NSUBMX)
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*)'SUBROUTINE GATHER_LOW_PAR entered'
      END IF        
      IZERO = 0
      VMAX = 1.0D99
      VMIN = -1.0D99
      IF(NTEST.GE.100) 
     & WRITE(6,*)'SET VMAX TO 1.0D99 and VMIN TO -1.0D99',VMAX,VMIN
      CALL SETVEC(TESTARR,VMAX,NSUBMX*NMPROC)
      CALL SETVEC(RECVARRAY,VMAX,NSUBMX*NMPROC)
      CALL ISETVC(IRECVARRAY,IZERO,NSUBMX*NMPROC)
      CALL ISETVC(ITESTARRAY,IZERO,NSUBMX*NMPROC)
      ITESTPURPOS = 0
      IF(MYPROC.NE.MASTER.AND.ITESTPURPOS.EQ.1)THEN
        CALL SETVEC(SUBVAL,VMAX,NSUBMX)
        CALL ISETVC(ISCAT,IZERO,NSUBMX)
      END IF
      ICOUNT = 1
      DO IL = 1, NSUBMX
        IRANKARR(IL) = ICOUNT
        ICOUNT = ICOUNT + 1
      ENDDO
*
*. gather lowest elements from nodes on all nodes - we need no master
      CALL MPI_ALLGATHER(SUBVAL,NSUBMX,MPI_REAL8,RECVARRAY,NSUBMX,
     &                  MPI_REAL8,IGATCOMM,IERR)
      CALL MPI_ALLGATHER(ISCAT,NSUBMX,my_MPI_INTEGER,IRECVARRAY,NSUBMX,
     &                  my_MPI_INTEGER,IGATCOMM,IERR)
      CALL MPI_ALLGATHER(IRANKARR,NSUBMX,my_MPI_INTEGER,ITESTARRAY,
     &                   NSUBMX,my_MPI_INTEGER,IGATCOMM,IERR)
*
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10)THEN
        WRITE(6,*)'before search: RECVARRAY and IRECVARRAY'
        CALL WRTMAT(RECVARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC)
        CALL IWRTMA(IRECVARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC)
        WRITE(6,*)'before search: ITESTARRAY'
        CALL IWRTMA(ITESTARRAY,1,NSUBMX*NMPROC,1,NSUBMX*NMPROC)
      END IF
                    
*. we need the NSUBMX lowest elements from all over the world - sorted
      ISEARCH = 0
      VMINLAST = VMIN
1000  CONTINUE
*
      VMINTMP = VMAX
      IMINRANK = 0
      IMINPLACE = NSUBMX + 1
      
      DO 300 II = 1,NSUBMX*NMPROC
*
         TEMPMIN   = RECVARRAY(II)
         ITEMPMIN  = IRECVARRAY(II)
         ITEMPRANK = ITESTARRAY(II) 
         ITEMPN    = II
*
         IF(ITEMPMIN.GT.0)THEN
           IF(TEMPMIN.LE.VMINTMP.AND.TEMPMIN.NE.VMAX)THEN
             IF(ITEMPRANK.EQ.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP)THEN
*               
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
               ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                 IF(ITEMPMIN.LT.IMINPLACE)THEN
*                 
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
                 END IF
               ELSE 
                 GOTO 300
               END IF
             ELSE IF(ITEMPRANK.EQ.1.AND.IMINRANK.NE.1)THEN
*                 
               IMINPLACE = ITEMPMIN
               VMINTMP   = TEMPMIN
               IMINNUMB  = ITEMPN
               IMINRANK  = ITEMPRANK
*                 
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
               END IF
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.NE.1)THEN
               IF(ITEMPRANK.EQ.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                   IF(ITEMPMIN.GT.IMINPLACE)THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                   END IF
                 ENDIF
               ELSE IF(ITEMPRANK.LT.IMINRANK) THEN
                 IF(TEMPMIN.LE.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 END IF
               ELSE IF(ITEMPRANK.GT.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 END IF                 
               END IF
             END IF
           END IF
         END IF
*
*
 300  CONTINUE
*
      ISEARCH = ISEARCH + 1
*.    VMINTMP should be the lowest value w.r.t. ISEARCH 
      VMINLAST = VMINTMP
*. test writing
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'VMINLAST',VMINLAST
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'PLACE FOR VMINLAST',IMINPLACE
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'NUMBER',IMINNUMB
      IF(MYPROC.EQ.MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'IMINRANK',IMINRANK
*. end of test writing
*. put the result to permanent storage in memory
      SUBVAL(ISEARCH) = VMINTMP
      ISCAT(ISEARCH) = IMINPLACE 
      RECVARRAY(IMINNUMB) = VMAX
      IRECVARRAY(IMINNUMB) = 0
      ITESTARRAY(IMINNUMB) = -1

      IF(ISEARCH.LT.NSUBMX) GOTO 1000
1001  CONTINUE
*
      END
                
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0LNSL_PAR(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM,SIGNTP)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMAT(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
C               INVMAT(A,B,MATDIM,NDIM)
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMAT(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0M1TD_PAR(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                     H0,SHIFT,WORK,XH0PSX,
     &                     NP1,NP2,NQ,VEC1,VEC2,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , September 1993
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar.inc"
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10)
     & write(6,*) ' H0M1TD_PAR , NPQDM = ', NPQDM
      IF(NPQDM.NE.0)THEN
       WRITE(6,*)'Problem in H0M1TD_PAR: NPQDM is not 0:',NPQDM
       Call Abend2('Problem in H0M1TD_PAR during parallel execution
     & detected')
      ENDIF
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
*. Obtain subspace components of input vector
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        IZERO = 0
        CALL GATVCD_PAR(LUIN,LBLK,NPQDM,IPNTR,WORK(KLV1),VEC1,
     &                 IZERO)
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL_PAR(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &                  WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &                  NTEST )
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
      CALL DMTVDS_PAR(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &               IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TD_PAR '
        WRITE(6,*) ' ============================= '
*. Note : works only if result vector is first file on LUOUT
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MICDV6_PAR(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,
     &                  NVAR,LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,
     &                  NINVEC,APROJ,AVEC,WORK,IPRTXX,NPRDIM,H0,IPNTR,
     &                  NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,E_CONV,IROOTHOMING,
     &                  LUWRTOUT,BSCRFAC,IBKINF,NBLOCKDN,IBLOCKI,SCRRED)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV6 
*
* parallel adaption, Stefan Knecht and Hans Joergen Aa. Jensen, 
* winter 2007 
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scratch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
       DIMENSION ISTAT(MPI_STATUS_SIZE)
       INTEGER STATUS_ARRAY(MPI_STATUS_SIZE,NUMNOD), IERR, IREQ(NUMNOD)
#include "lucipar.inc"
       DIMENSION VEC1(*),VEC2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*), SCRRED(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*      SCRRED : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
       REAL*8 INPRDD, INPROD 
       DOUBLE PRECISION tottime,endtime,starttime
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*      characters used for timing
       CHARACTER WALLTSTEP*12,  WPART22*12
*. Notice XJEP is also used for ROOTHOMING, should be allocated
* outside (for roothoming :dim = 3*MAXVEC )
       DIMENSION XJEP(10000)
       INTEGER   IXJEP(10000)
       DIMENSION BSCRFAC(*)
       DIMENSION IBKINF(*),NBLOCKDN(IBLOCKI)
       INTEGER NDATATYPE,NREDTYPE
       REAL*8 REDSCRVAR
*
       WRITE(LUWRTOUT,*) 
     & ' '
       WRITE(LUWRTOUT,*) 
     & ' '
       WRITE(LUWRTOUT,*) 
     & '                        :::::::::::::::::::::::'
       WRITE(LUWRTOUT,*)          
     & '                        : Entering MICDV6_PAR :'
       WRITE(LUWRTOUT,*) 
     & '                        :::::::::::::::::::::::'
       WRITE(LUWRTOUT,*) 
     & ' '
       WRITE(LUWRTOUT,*) 
     & ' '
*
       IPICO = 0
       IF(IPICO.EQ.0) THEN
C?       WRITE(6,*) ' Variational  solver '
       END IF
*
       IPRT = 000
c      IPRT = IPRTXX
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRTOUT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRTOUT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(6,*) ' Enforced stop on MICDV6 '
         Call Abend1( 20 )
       END IF
*
       IF(IROOTHOMING.EQ.1) THEN
         WRITE(LUWRTOUT,*) ' Root homing performed '
       END IF
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
CSK       CALL REWINO(LU1)
CSK       IF(MYPROC.EQ.MASTER) THEN
       IROUNDIT = 0
 4    CONTINUE
CSK       DO I = 1, NINVEC
CSK         IF(IROUNDIT.EQ.MYPROC)THEN
CSK           write(6,*)'initial vector',I,'of myproc',myproc
CSK           CALL WRTVCD(VEC1,LU1,0,-1)
CSK         END IF
CSK         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CSK       END DO
CSK       IROUNDIT = IROUNDIT + 1
CSK       IF(IROUNDIT.LT.NMPROC) GOTO 4
       
CSK       END IF
                    
                                                  
       CALL REWINO(LU1)
       CALL REWINO(LU2)
*
       DO IVEC = 1,NINVEC
         CALL REWINO(LU5)
         CALL REWINO(LU6)
         CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         CALL MV7_PAR(VEC1,VEC2,LU5,LU6,BSCRFAC,IBKINF,NBLOCKDN,IBLOCKI)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU6)
         CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
           DO JVEC = 1, IVEC
              CALL REWINO(LU5)
              IJ = IVEC*(IVEC-1)/2 + JVEC
              APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
           END DO
*          ^ End of loop over JVEC
       END DO
*      ^ End of loop over IVEC
*. end of initial loop
*
*. every node has its part of APROJ --> reduce to each node!
        NDATATYPE = 2
        NREDTYPE = -1
        IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
        CALL REDVEC(APROJ,SCRRED,IREDL,NDATATYPE,MPI_SUM,
     &              MPI_COMM_WORLD,NREDTYPE)
        CALL COPVEC(SCRRED,APROJ,IREDL)
*
CSK       IF(MYPROC.EQ.MASTER+1) THEN
CSK         WRITE(6,*)'MAXVEC = ',MAXVEC
CSK         CALL WRTMAT(APROJ,1,MAXVEC*(MAXVEC+1)/2,1,MAXVEC*(MAXVEC+1)/2)
CSK       ENDIF
        
           WRITE(LUWRTOUT,'(A)') 
     &'         ::::::::::::::::::::::::::::::::::::::::::::::
     &:::::::: '
           WRITE(LUWRTOUT,'(A,1X,I4,1X,A)')
     &'         ::     Parallel calculation will use',
     & NMPROC,'CPU''s     ::'
           WRITE(LUWRTOUT,'(A)') 
     &'         ::::::::::::::::::::::::::::::::::::::::::::::
     &:::::::: '
           WRITE(LUWRTOUT,'(A)')' '
*
       IF( IPRT .GE.5 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
CSK       IF(MYPROC.EQ.MASTER) THEN
         CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
         CALL EIGEN(WORK(KAPROJ),AVEC,NINVEC,0,1)
         DO IROOT = 1, NROOT
           EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
         END DO
CSK       END IF
*. Master shares all eigenvalues with the nodes
*. SK: disabled because of the fact that all nodes can also execute the same 
*. code above, saves communication.
CSK       DO IROOT = 1, NROOT
CSK         CALL MPI_BCAST(EIG(1,IROOT),1,MPI_REAL8,MASTER,
CSK     &                  MPI_COMM_WORLD,IERR)
CSK       END DO
CSK       CALL MPI_BCAST(AVEC,MAXVEC**2,MPI_REAL8,MASTER,
CSK     &                MPI_COMM_WORLD,IERR)
*
*. end of disabled code!
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,*) ' Eigenvalues of initial iteration (with
     &   shift)',myproc
         WRITE(6,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(6,*) ' Initial set of eigen values (no shift) ',myproc
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
C     TRAVCD(VEC1,VEC2,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
       CALL REWINO(LU3)
       CALL TRAVCD_PAR(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,
     &                 LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINO(LU3)
       CALL TRAVCD_PAR(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,
     &                 LBLK,LU4,LU5)
*
       IF(IPRT.GE.600) THEN
         WRITE(6,*) ' Initial set of eigenvectors '
         CALL REWINO(LU1)
         DO IROOT = 1, NROOT
           WRITE(6,1066) IROOT, LBLK, LBL
 1066 FORMAT(/'lets see IROOT =', 3I4/)
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(6,*) ' Initial set of sigma vectors '
         CALL REWINO(LU2)
         DO IROOT = 1, NROOT
           WRITE(6,1067) IROOT
 1067 FORMAT(/'lets see IROOT again=', I4/)
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       ZERO = 0.0D0
       CALL SETVEC(APROJ,ZERO,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT
       IF (MAXIT .EQ. 1 ) GOTO  1001
*
* ======================
*. Loop over iterations
* ======================
*
CSK         print*,'NP1,NP2,NQ ',NP1,NP2,NQ
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        WALLITR1 = MPI_WTIME()
        write(LUWRTOUT,*)
        write(LUWRTOUT,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRTOUT,*) '_______________________'
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINO(LU1)
       CALL REWINO(LU2)
*
       DO 100 IROOT = 1, NROOT
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINO(LU7)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINO(LU7)
         CALL REWINO(LU5)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(6,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         REDSCRVAR = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
         CALL REDVEC(REDSCRVAR,SCRRED,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
         CALL COPVEC(SCRRED,REDSCRVAR,1)
         RNORM = SQRT(REDSCRVAR)
         RNRM(ITER-1,IROOT) = RNORM
         WRITE(LUWRTOUT,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
*
         CALL FLSHFO(6)
*
         IF(RNORM.LT. TEST .OR.
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
CSK         print*,'NP1,NP2,NQ ',NP1,NP2,NQ
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(6,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU5)
           CALL REWINO(LU6)
           CALL H0M1TD_PAR(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C     H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
CSK        print*,'IOLSTM ',IOLSTM
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD_PAR(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             REDSCRVAR = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
             CALL REDVEC(REDSCRVAR,GAMMA,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
CSK             WRITE(LUWRTOUT,*) ' Gamma = ', Gamma
             IF(IPRT.GE.1000) WRITE(6,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINO(LU5)
             CALL REWINO(LU7)
             SCRRED(1) = VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK)
             CALL REDVEC(SCRRED(1),REDSCRVAR,1,2,MPI_SUM,
     &                   MPI_COMM_WORLD,-1)
              VNORM = SQRT(REDSCRVAR)
              IF(IPRT.GE.1000) write(6,*) ' VNORM = ', VNORM
CSK              write(LUWRTOUT,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                REDSCRVAR = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                CALL REDVEC(REDSCRVAR,DELTA,1,2,MPI_SUM,
     &                      MPI_COMM_WORLD,-1)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*) ' DELTA,GAMMA,
     &          FACTOR,MYPROC'
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*)   DELTA,GAMMA,
     &          FACTOR,MYPROC
CSK                WRITE(LUWRTOUT,*) ' DELTA,GAMMA,FACTOR,MYPROC'
CSK                WRITE(LUWRTOUT,*)   DELTA,GAMMA,FACTOR,MYPROC
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
*.. Vectors on LU1
*
           CALL REWINE( LU1 ,LBLK)
           DO IVEC = 1, NROOT
               CALL REWINO(LU6)
               WORK(IVEC) = -INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
           END DO
           CALL REDVEC(WORK,SCRRED,NROOT,2,MPI_SUM,MPI_COMM_WORLD,-1)
           CALL COPVEC(SCRRED,WORK,NROOT)
           IF(IPRT.GE.1000) THEN
             Write(LUWRTOUT,*) ' Overlap with vectors on LU1',myproc
             call wrtmat(work,1,nroot,1,nroot)
           END IF
           ONE = 1.0D0
           CALL REWINO(LU1)
           CALL MVCSMD(LU1,WORK,LU7,LU5,VEC1,VEC2,NROOT,1,LBLK)
           CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU6,LU5,1,LBLK)
           IF(IPRT.GE.1000) THEN
             write(6,*) ' orthogonalized to vectors on LU1'
             CALL WRTVCD(VEC1,LU5,1,LBLK)
           END IF
*.. Vectors on LU3
           IF(NVEC+IADD-1-NROOT.GT.0) THEN
             CALL REWINE( LU3 ,LBLK)
             DO IVEC = 1, NVEC+IADD-1-NROOT
                 CALL REWINO(LU6)
                 WORK(IVEC) = -INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
             END DO
             CALL REDVEC(WORK,SCRRED,NVEC+IADD-1-NROOT,2,MPI_SUM,
     &                   MPI_COMM_WORLD,-1)
             CALL COPVEC(SCRRED,WORK,NVEC+IADD-1-NROOT)
             ONE = 1.0D0
             CALL REWINO(LU3)
             CALL MVCSMD(LU3,WORK,LU7,LU6,VEC1,VEC2,NVEC+IADD-1-NROOT,
     &            1,LBLK)
             CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU5,LU6,1,LBLK)
           ELSE
             CALL REWINO(LU3)
             CALL COPVCD(LU5,LU6,VEC1,1,LBLK)
           END IF
*. 1.4 Normalize vector
           REDSCRVAR = INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
           CALL REDVEC(REDSCRVAR,SCALE,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU6,LBLK)
           CALL SCLVCD(LU6,LU3,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU6,LU7,FACTOR,VEC1,1,LBLK)
             WRITE(6,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU7,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) THEN
         GOTO  1001
      END IF
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
      IF(IPRT.GE.1000) THEN
       IF (MYPROC.EQ.MASTER) THEN
        WRITE(LUWRTOUT,*) ' Vectors on LU3'
        WRITE(LUWRTOUT,*) 'NVEC-NROOT+IADD =',NVEC-NROOT+IADD
        CALL REWINO(LU3)
        DO IVEC = 1, NVEC-NROOT+IADD
          CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
       END IF
       CALL REWINO(LU3)
      END IF
*
*. 
      starttime = MPI_WTIME()
*
      CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
      CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
      IMUSTRED = 0
      
 150  DO IVEC = 1, IADD
        CALL REWINE(LU5,LBLK)
        CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
*.      Check timings for sigma vector generation
***********************************************************************
        sigmatime = MPI_WTIME()
        CALL MV7_PAR(VEC1,VEC2,LU5,LU6,BSCRFAC,IBKINF,NBLOCKDN,IBLOCKI)
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        WRITE(LUWRTOUT,9400) WALLTID
***********************************************************************
        CALL REWINE(LU6,LBLK)
        CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*2.2 Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO JVEC = 1, NROOT
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
          IMUSTRED = IMUSTRED + 1
          IF(IVEC.EQ.1.AND.JVEC.EQ.1) ISTRED = IJ
        END DO

        CALL REWINE(LU3,LBLK)
        DO JVEC = NROOT+1, NVEC+IVEC
         CALL REWINE(LU6,LBLK)
         IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
         APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
         IMUSTRED = IMUSTRED + 1
        END DO
      END DO
*     /\ End do over new trial vectors
      CALL REDVEC(APROJ(ISTRED),SCRRED,IMUSTRED,2,MPI_SUM,
     &            MPI_COMM_WORLD,-1)
      CALL COPVEC(SCRRED,APROJ(ISTRED),IMUSTRED)
      
CSK      IF(MYPROC.EQ.MASTER) THEN
CSK        write(LUWRTOUT,*)'APROJ is'
CSK        CALL WRTMAT(APROJ,1,MAXVEC*(MAXVEC+1)/2,1,MAXVEC*(MAXVEC+1)/2)
CSK      ENDIF
      
C
      endtime = MPI_WTIME()
      tottime = 0.D0
      tottime = endtime - starttime
*. TIMING FOR PARTS 2.1 - 2.2 
      WPART22 = SECTID(tottime) 
      WRITE(LUWRTOUT,9350) WPART22
*
*. 2.3 Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN(WORK(KAPROJ),AVEC,NVEC,0,1)
*
CSK Jan07: since IPICO is set to 0, I comment this part out
C      IF(IPICO.NE.0) THEN
C        E0VAR = WORK(KAPROJ)
C        C0VAR = AVEC(1)
C        C1VAR = AVEC(2)
C        C1NRM = SQRT(C0VAR**2+C1VAR**2)
C*. overwrite with pert solution
C        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
C        AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
C        E0PERT = AVEC(1)**2*APROJ(1)
C     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
C     &         + AVEC(2)**2*APROJ(3)
C        WORK(KAPROJ) = E0PERT
C        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
C        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
C        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
C      END IF
CSK end-of 
*
      IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       CALL REDVEC(XJEP,SCRRED,NVEC,2,MPI_SUM,MPI_COMM_WORLD,-1)
       CALL COPVEC(SCRRED,XJEP,NVEC)
       WRITE(LUWRTOUT,*)
     & ' Norm of projections to previous vector space '
       CALL WRTMAT(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(XJEP,ONEM,NVEC)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(LUWRTOUT,*) ' Warning : Some elements lost in sorting '
         WRITE(LUWRTOUT,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
       END DO
       DO INEW = 1, NVEC
         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' Reordered WORK and AVEC arrays '
         CALL PRSYM(WORK,NVEC)
         CALL WRTMAT(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
*
      IF(IPRT .GE. 3 ) THEN
        WRITE(LUWRTOUT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
        WRITE(LUWRTOUT,'(A)') ' Norm of Residuals (Previous it) '
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
      END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRTOUT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRTOUT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,NVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*. 3.1 : Current wave function approximations, collect on LU7
*
*  check timing of part3 in loop iteration
      CALL GETTIM(CPUSTEP3,WALLSTEP3)
      WALLSTEP3 = MPI_WTIME()
*
      IF(IPRT.GE.1000) THEN
      write(LUWRTOUT,*) ' I am going to reset '
      write(LUWRTOUT,*) ' nroot, nvec ', nroot,nvec
      END IF
      IF(IPRT.GE.1000) THEN
        WRITE(LUWRTOUT,*) ' Initial vectors on LU1'
        CALL REWINO(LU1)
        DO IVEC = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
        END DO
        WRITE(LUWRTOUT,*) ' Initial vectors on LU3'
        CALL REWINO(LU3)
        DO IVEC = 1, NROOT
           CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
      END IF
*
      CALL REWINE( LU7,LBLK)
      DO IROOT = 1, NROOT
*. From LU1 to LU5
        CALL REWINO(LU1)
        CALL MVCSMD_PAR(LU1,AVEC((IROOT-1)*NVEC+1),
     &  LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 1'
*. and add LU3 to LU5
        CALL REWINO(LU3)
        ONE = 1.0D0
        CALL MVCSMD2_PAR(LU3,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &  LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 2'
*
        REDSCR = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
        CALL REDVEC(REDSCR,XNORM,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
        CALL REWINE(LU5,LBLK)
        SCALE  = 1.0D0/SQRT(XNORM)
CSK        WRITE(6,*) ' SCALE = ', SCALE
        WORK(IROOT) = SCALE
*. scale LU5 => LU7
        CALL REWINO(LU5)
        CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 3'
      END DO
*. Transfer C vectors from LU7 to LU1
      CALL REWINE( LU7,LBLK)
      CALL REWINE( LU1,LBLK)
      DO IVEC = 1,NROOT
        CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
      END DO
      IF(IPRT.GE.1000)
     &write(6,*) ' end of c reset, part 4'
      IF(IPRT.GE.1000) THEN
        WRITE(6,*) ' Reset C-vectors on LU1 (myproc)',myproc
        CALL REWINO(LU1)
        DO IVEC = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
        END DO
      END IF
*. 3.2 : corresponding sigma vectors
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU2 to LU5
          CALL REWINO(LU2)
          CALL MVCSMD_PAR(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 1'
*. and add LU4 to LU5
          CALL REWINO(LU4)
          CALL MVCSMD2_PAR(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 2'
          SCALE  = WORK(IROOT)
*. scale LU5 => LU7
          CALL REWINO(LU5)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(6,*) ' end of s reset, part 3'
        END DO
*. Transfer sigma  vectors from LU7 to LU2
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU2,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 4'
        NNVEC = NROOT
*. Overlap between first vectors on LU1 and LU2
        REDSCR = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
        CALL REDVEC(REDSCR,E11,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
CSK        WRITE(LUWRTOUT,*) ' MYPROC: <Lu1!Lu2> = ', MYPROC,E11
        IF(3*NROOT.LE.MAXVEC) THEN
*
*. Orthogonalize the
*. last set of correction vectors to the current
*. eigenvectors on LU1, and save on LU2
*. Overlap with root approximations
*. Start of last set of trial vectors
          ISTART = NVEC-NROOT-IADD+1
          CALL SKPVCD(LU3,ISTART-1,VEC1,1,LBLK)
*
          CALL REWINO(LU5)
          DO JVEC = 1, IADD
*. Orthogonalize to vectors on LU1
            CALL REWINO(LU7)
            CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
            IF(IPRT.GE.1000) THEN
            write(6,*) ' Initial vector on LU7'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
            CALL REWINO(LU1)
            IMUSTRED = 0
            DO IROOT = 1, NROOT
              CALL REWINO(LU7)
              WORK(IROOT+(JVEC-1)*2*NROOT) = 
     &       -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
              IF(IROOT.EQ.1) ISTRED = IROOT+(JVEC-1)*2*NROOT
              IMUSTRED = IMUSTRED + 1
            END DO
            CALL REDVEC(WORK(ISTRED),SCRRED,IMUSTRED,2,
     &                  MPI_SUM,MPI_COMM_WORLD,-1)
            CALL COPVEC(SCRRED,WORK(ISTRED),IMUSTRED)
                 
            IF(IPRT.GE.1000)
     &      write(6,*) ' c, part1 finito'
*. And to trial vectors on LU5
            CALL REWINO(LU5)
            IMUSTRED = 0
            DO KVEC = 1, JVEC-1
              CALL REWINO(LU7)
              WORK(NROOT+KVEC+(JVEC-1)*2*NROOT) = 
     &       -INPRDD(VEC1,VEC2,LU7,LU5,0,LBLK)
              IF(KVEC.EQ.1) ISTRED = NROOT+KVEC+(JVEC-1)*2*NROOT
              IMUSTRED = IMUSTRED + 1
            END DO
            IF(IMUSTRED.GT.0)THEN
              CALL REDVEC(WORK(ISTRED),SCRRED,IMUSTRED,2,
     &                    MPI_SUM,MPI_COMM_WORLD,-1)
              CALL COPVEC(SCRRED,WORK(ISTRED),IMUSTRED)
            ENDIF
            WORK(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
            IF(IPRT.GE.1000)
     &      write(6,*) ' c, part2 finito'
*
            ONE = 1.0D0
            CALL MVCSMD2_PAR(LU1,WORK(1+(JVEC-1)*2*NROOT),ONE ,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            ONE = 1.0D0
            CALL MVCSMD2_PAR(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000) THEN
              write(6,*) ' c, part4 finito'
              write(6,*) ' Vector after sec ort'
              CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
*
            REDSCR = INPRDD(VEC1,VEC1,LU7,LU7,1,LBLK)
            CALL REDVEC(REDSCR,FACTOR,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
            SCALE = 1.0D0/SQRT(FACTOR)
            CALL SCALVE(WORK((JVEC-1)*2*NROOT+1),SCALE,
     &           NROOT+JVEC)
            CALL REWINO(LU7)
            CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &       write(6,*) ' c, part6 finito'
         END DO
*        /\ End of loop over orthogonalized directions
*. Transfer modified directions to LU3
         CALL REWINO(LU3)
         CALL REWINO(LU5)
         DO JVEC =1, IADD
           CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000) THEN
           write(6,*) ' c, part7 finito'
           WRITE(6,*) ' Additional trial vectors on LU3 myproc',myproc
           CALL REWINO(LU3)
           DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
           END DO
         END IF
* Sigma vectors corresponding to orthogonalized directions
         CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
         CALL REWINO(LU5)
         DO JVEC = 1, IADD
            CALL REWINO(LU7)
            CALL COPVCD(LU4,LU7,VEC1,0,LBLK)
*
            FACT = WORK(NROOT+JVEC+(JVEC-1)*2*NROOT)
            CALL MVCSMD2_PAR(LU2,WORK(1+(JVEC-1)*2*NROOT),FACT,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 1 finito '
            ONE = 1.0D0
            CALL MVCSMD2_PAR(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 2 finito '
            CALL REWINO(LU7)
            CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
            IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 3 finito '
         END DO
*        /\ End of loop over orthogonalized directions
*. Copy back to LU4
         CALL REWINO(LU4)
         CALL REWINO(LU5)
         DO JVEC = 1, IADD
            CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000)
     &   write(LUWRTOUT,*) ' s, part 4 finito '
*
         NNVEC = NROOT + IADD
       END IF
*      /\ End if more than NROOT vectors to be reset
*
       IF(IPRT.GE.1000) THEN
         WRITE(LUWRTOUT,*)' Additional sigma vectors on LU4'
         CALL REWINO(LU4)
         DO JVEC = 1, IADD
           CALL WRTVCD(VEC1,LU4,0,LBLK)
         END DO
       END IF
CSK      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CSK      IF(IODENSEPAR.EQ.1) THEN
CSK         CALL Abend2('Stefan forced me to stop at s part4 finito... 
CSK     & DV6_PAR')
CSK      END IF
*
       NVEC = NNVEC
*. New subspace
*. Calculate subspace Hamiltonian from actual vectors
*. on disc
       IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian'
       CALL REWINO(LU1)
       CALL REWINO(LU3)
       DO IVEC = 1, NVEC
*
         CALL REWINO(LU5)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINO(LU2)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
           IF(IVEC.EQ.1.AND.JVEC.EQ.1) ISTRED = IJ
           IMUSTRED = IMUSTRED + 1
         END DO
         CALL REWINO(LU4)
         DO JVEC = NROOT+1,IVEC
           CALL REWINO(LU5)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
           IMUSTRED = IMUSTRED + 1
         END DO
       END DO
       CALL REDVEC(APROJ(ISTRED),SCRRED,IMUSTRED,2,MPI_SUM,
     &             MPI_COMM_WORLD,-1)
       CALL COPVEC(SCRRED,APROJ(ISTRED),IMUSTRED)
       if (IPRT.ge.2) then
         write(6,*) ' Reset hamiltonian'
         call prsym(aproj,nvec)
       end if
*. Test : Orthogonality of new vectors
CSK       CALL REWINO(LU1)
CSK       CALL REWINO(LU3)
CSK       DO IVEC = 1, NVEC
CSK*
CSK         CALL REWINO(LU5)
CSK         IF(IVEC.LE.NROOT) THEN
CSK           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
CSK         ELSE
CSK           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
CSK         END IF
CSK*
CSK         CALL REWINO(LU1)
CSK         DO JVEC = 1, MIN(IVEC,NROOT)
CSK           CALL REWINO(LU5)
CSK           IJ = IVEC*(IVEC-1)/2+JVEC
CSK           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU1,0,LBLK)
CSK         END DO
CSK         CALL REWINO(LU3)
CSK         DO JVEC = NROOT+1,IVEC
CSK           CALL REWINO(LU5)
CSK           IJ = IVEC*(IVEC-1)/2+JVEC
CSK           XJEP(IJ) = INPRDD(VEC1,VEC2,LU5,LU3,0,LBLK)
CSK         END DO
CSK       END DO
CSK       if (IPRT.ge.2) then
CSK         write(6,*) ' Overlap matrix    '
CSK         call prsym(xjep,nvec)
CSK       end if
*
*
*  check timing of part3 in loop iteration
       WALLSTEP4 = MPI_WTIME()
       WALLTSTEP = SECTID(WALLSTEP4-WALLSTEP3)
       WRITE(LUWRTOUT,9600) WALLTSTEP
*
*  Timing of this iteration
       CALL GETTIM(CPUITR2,WALLITR2)
       CPUTID = SECTID(CPUITR2-CPUITR1)
       WALLITR2 = MPI_WTIME()
       WALLTID = SECTID(WALLITR2-WALLITR1)
       WRITE(LUWRTOUT,9300) WALLTID
*. End of resetting business
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
* ( End of loop over iterations )
*
* let's synchronize all processes!
*
 1003  DO 1004 I = 1, 2
        CALL MPI_BARRIER(ICOMM, IERR)
        TIME1 = MPI_Wtime()
        WRITE(LUWRTOUT,*) 'I AM IN SYNC PROC., MYPROC = ', MYPROC
        CALL MPI_BARRIER(ICOMM, IERR)
        TOTAL_TIME = MPI_Wtime() - TIME1
 1004  CONTINUE
*
       WRITE(LUWRTOUT,*) 'SYNC. PROCESS HAS USED', TOTAL_TIME, 'SEC.'
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(LUWRTOUT,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRTOUT,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          write(LUWRTOUT,*)
          write(LUWRTOUT,*) '------------------------'
          write(LUWRTOUT,*) 'Root number  ',IROOT
          write(LUWRTOUT,*) '------------------------'
          DO I=1,ITER
            WRITE(LUWRTOUT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
          WRITE(LUWRTOUT,*)
          WRITE(LUWRTOUT,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRTOUT,*)
     &    '============================================'
          WRITE(LUWRTOUT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(LUWRTOUT,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(LUWRTOUT,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        write(LUWRTOUT,*)
 1310   FORMAT
     &  (1H0,' Iteration point        Eigenvalue         Residual ')
        DO IROOT = 1, NROOT
          write(LUWRTOUT,*) '------------------------'
          write(LUWRTOUT,*) 'Root number  ',IROOT
          write(LUWRTOUT,*) '------------------------'
          DO I=1,ITER
            WRITE(LUWRTOUT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
        end do
        write(LUWRTOUT,*)
        write(LUWRTOUT,*) '**********************************'//
     &             '**************************'
        write(LUWRTOUT,*) '   Iter  Root       Energy        '//
     &             'RESIDUAL     RESRATIO '
        write(LUWRTOUT,*) '**********************************'//
     &             '**************************'
        write(LUWRTOUT,*)
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           WRITE(LUWRTOUT,'(3X,I3,2X,I3,3X,F18.10,2X,E10.3,3X,E10.3)')
     &               ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &               RNRM(1,IROOT)/RNRM(ITER,IROOT)
           write(LUWRTOUT,'(A,F18.10)') 
     &               ' Final energy ',FINEIG(IROOT)
 1601   CONTINUE
      END IF
*
      CALL FLSHFO(6)
*
      Call Add_Info('E_CI',FINEIG,NROOT)
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION         : ',A)
 9400 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL         : ',A)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION: ',A)
 9350 FORMAT(' >>>  WALL TIME FOR PART 2.1 - 2.2            : ',A)
      END
***********************************************************************
#else
      SUBROUTINE PAR_DUMMY
      END 
#endif /* defined (VAR_MPI)*/
