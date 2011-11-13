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
C
C
C VERSION : $Revision: 1.26.2.2 $
C DATE    : $Date: 2004/09/01 08:16:07 $
C FILE    : pamgrd.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgrd */
      SUBROUTINE PAMGRD(INDXPR,GPCI,GPOE,GPOP,
     &                  JOPSY,JCONF,JXOPE,JXOPP,NCONF,NXOPE,NXOPP,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct gradient for property of index INDXPR in /XCBPRP/
C     GPC - configurational part
C     GPO - orbital part
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
C
      DIMENSION GPCI(*),GPOE(*),GPOP(*),JCONF(*),JXOPE(*),JXOPP(*),
     &          WORK(*)
      CALL QENTER('PAMGRD')
C
C     Consistency check
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
      IF(IOPSY.NE.JOPSY) THEN
        WRITE(LUPRI,'(A)') 
     &    '*** ERROR(PAMGRD): Incorrect fermion ircop of property ',
     &    PRPNAM(INDXPR)
        CALL QUIT('*** ERROR(PAMGRD) Incorrect fermion ircop')
      ENDIF
C
C     Configurational part
C     ====================
C
      IF(NCONF.GT.0) THEN
C        CALL GRDCI(INDXPR,GPC,WORK,KFREE,LFREE,IPRINT)
      ENDIF      
C
C     Orbital part
C     ============
C
      NXOPT = NXOPE + NXOPP
      IF(NXOPT.GT.0) THEN
        CALL GRDORB(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &              WORK,KFREE,LFREE,IPRINT)
      ENDIF
C
      CALL QEXIT('PAMGRD')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdorb */
      SUBROUTINE GRDORB(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | g  | 
C     E[1]_X = | g* | ; 
C
C     where
C
C     g_{pq} = <0|[-q(+)p,H_X]|0> = -h_{X;ps}D_{qs} + h_{X;rq}D_{rp}
C
C     Active-inactive block:
C
C     g_{ai} = -h_{X;ai}
C
C     Orbital classes:
C
C       Inactive orbitals  :  i,j,k,l
C       Active orbitals    :  v,w,x,y,z
C       Secondary orbitals :  a,b,c,d,e
C       General orbitals   :  m,n,o,p,q,r,s,t,u,v
C
C     GRDORB constructs g_{pq}.
C     It employs the orbital rotation index JXOP(p,q)
C     with the following ordering:
C
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C          AB : active-active
C
C     Indices are therefore permuted...
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
C
      DIMENSION GPOE(*),GPOP(*),JXOPE(*),JXOPP(*),WORK(*)
C
      CALL QENTER('GRDORB')
      KFRSAV = KFREE
C
      CALL MEMGET('REAL',KPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO  ,NCMOTQ ,WORK,KFREE,LFREE)
      IF(SPINFR.OR.LEVYLE) THEN
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
      ELSE
        KIBRP = KFREE
      ENDIF
C
      CALL GRDOR1(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &            WORK(KPRPMO),WORK(KCMO),WORK(KIBRP),
     &            WORK,KFREE,LFREE,IPRINT) 
C
      CALL MEMREL('GRDORB',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('GRDORB')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdor1 */
      SUBROUTINE GRDOR1(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,NXOPP,IOPSY,
     &                  PRPMO,CMO,IBEIG,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | g  | 
C     E[1]_X = | g* | ; 
C
C     where
C
C     g_{pq} = <0|[-q(+)p,H_X]|0> = -h_{X;ps}D_{qs} + h_{X;rq}D_{rp}
C
C     Active-inactive block:
C
C     g_{ai} = -h_{X;ai}
C
C     Orbital classes:
C
C       Inactive orbitals  :  i,j,k,l
C       Active orbitals    :  v,w,x,y,z
C       Secondary orbitals :  a,b,c,d,e
C       General orbitals   :  m,n,o,p,q,r,s,t,u,v
C
C     GRDORB constructs g_{pq}.
C     It employs the orbital rotation index JXOP(p,q)
C     with the following ordering:
C
C          IA : inactive-active   
C          IS : inactive-secondary
C          AS : active-secondary
C          AB : active-active
C
C     Indices are therefore permuted...
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 (hjaaj, new parameter list for PRPMAT)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbham.h"
C
      DIMENSION GPOE(*),GPOP(*),JXOPE(*),JXOPP(*),
     &          PRPMO(NORBT,NORBT,NZ),CMO(*),IBEIG(*),WORK(*)
C
      CALL QENTER('GRDOR1')
C
C     Read coefficients
C     =================
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GRDOR1')
      IF(SPINFR.OR.LEVYLE) THEN
        CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IBEIG,TOTERG,10,
     &              WORK,KFREE,LFREE)
      ELSE
        CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,DUM,TOTERG,2,
     &              WORK,KFREE,LFREE)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Get property matrix in MO-basis and scale by factor -1
C     ======================================================
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL PRPMAT(INDXPR,IOPSY,PRPMO,CMO,IBEIG,ICMOQ,NORB,
     &            WORK,KFREE,LFREE,IPRINT)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Hermicity of operator
C     =====================
C
      ITIM = IPRPTIM(INDXPR)
      IREP = IPRPSYM(INDXPR) - 1
      IQ   = IPQTOQ(1,IREP)
      IH   = IHQMAT(IQ,ITIM)
C
C
C
C     Construct gradient - electronic part
C     ====================================
C
      KFRSAV = KFREE
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRINT,.FALSE.,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
      IF(NXOPE.GT.0) THEN
C
C
        CALL DZERO(GPOE,NZ*NXOPE)
        CALL PRPORB(GPOE,JXOPE,NXOPE,PRPMO,NORBT,WORK(KDV),NZ,IH)

c       no scaling
        CALL DSCAL(NZ*NXOPE,DM1,GPOE,1)
      ENDIF
C
C     Construct gradient - positronic part
C     ====================================
C
      IF(NXOPP.GT.0) THEN
        CALL DZERO(GPOP,NZ*NXOPP)
        CALL PRPORB(GPOP,JXOPP,NXOPP,PRPMO,NORBT,WORK(KDV),NZ,IH)
c       no scaling
        CALL DSCAL(NZ*NXOPP,DM1,GPOP,1)
      ENDIF
      CALL MEMREL('GRDOR1',WORK,1,KFRSAV,KFREE,LFREE)
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from GRDORB',-1)
        IREP = IPRPSYM(INDXPR)-1
        IF(NXOPE.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-e) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOE,NXOPE,1,NXOPE,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
        IF(NXOPP.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-p) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOP,NXOPP,1,NXOPP,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
      ENDIF
      CALL QEXIT ('GRDOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xopgat */
      SUBROUTINE XOPGAT(QVEC,JXOP,NXOP,QMAT,NDIM,NZ)
C***********************************************************************
C
C     Gather matrix elements into vector using JXOP
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ)
C
      DO IZ = 1,NZ
#include "ivdep.h"
        DO IXOP = 1,NXOP
          IA = JXOP(1,IXOP)
          IM = JXOP(2,IXOP)
          QVEC(IXOP,IZ) = QMAT(IM,IA,IZ)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xopsct */
      SUBROUTINE XOPSCT(QVEC,JXOP,NXOP,QMAT,NDIM,NZ)
C***********************************************************************
C
C     Scatter vector elements into matrix using JXOP
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ)
C
      DO IZ = 1,NZ
#include "ivdep.h"
        DO IXOP = 1,NXOP
          IA = JXOP(1,IXOP)
          IM = JXOP(2,IXOP)
          QMAT(IM,IA,IZ) = QVEC(IXOP,IZ)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpmat */
      SUBROUTINE PRPMAT(INDXPR,IOPSY,PRPMO,CMO,IBEIG,IPCMO,NPMO,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct MO-matrix for property INDXPR listed in /XCBXPR/
C     -- NPMO(i) is dimension of output property matrix 
C                in fermion symmetry i
C     -- IPCMO(i) is off-set to first MO-coeff in fermion symmetry i
C        NPMO = NORB and IPCMO = ICMOQ gives full PRPMO matrix
C        NPMO = NASH and IPCMO = xxxxx gives      PRPMO matrix over 
C                                active orbitals
C
C     Written by T.Saue July 1996
C     Revisions:
C       Feb 17 2000 (hjaaj, IPCMO(2), NPMO(2) parameters for 
C                    partial MO transf.)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION PRPMO(*),CMO(*),IBEIG(*),IPCMO(*),NPMO(*),WORK(*)
C
      CALL QENTER('PRPMAT')
      KFRSAV = KFREE
C
C     Memory allocation
C
      CALL MEMGET('REAL',KPRPAO ,N2BBASXQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTRI  ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN ,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST ,NZ        ,WORK,KFREE,LFREE)
C
      NPMOT = NPMO(1)
      IF (NFSYM.EQ.2) NPMOT = NPMOT + NPMO(2)
      N2PMOX = NPMOT*NPMOT
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,'(/A,4I10)')
     &      ' Input to PRPMAT: IPCMO, NPMO:',
     &      (IPCMO(I),NPMO(I),I = 1,NFSYM)
D        WRITE(LUPRI,*)' NPMOT,N2PMOX:',NPMOT,N2PMOX
      END IF
      CALL PRPMA1(INDXPR,PRPMO,CMO,IBEIG,IPCMO,NPMO,N2PMOX,
     &            WORK(KPRPAO),WORK(KPTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)

      CALL MEMREL('PRPMAT',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPMAT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpma1 */
      SUBROUTINE PRPMA1(INDXPR,PRPMO,CMO,IBEIG,IPCMO,NPMO,N2PMOX,PRPAO,
     &                  PTRI,OP1INT,FIRST,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct MO-matrix for property INDXPR listed in /XCBXPR/
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 (hjaaj, see PRPMAT)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbxpr.h"
C
      LOGICAL FIRST(*)
      DIMENSION PRPMO(N2PMOX,NZ),CMO(*),IBEIG(*),IPCMO(2),NPMO(2),
     &          PRPAO(N2BBASX,NZ),PTRI(NNBBASX,NZ),
     &          OP1INT(NNBBASX),WORK(*)
      DIMENSION IPMO(2), I2PMOX(2,2)
C
      CALL QENTER('PRPMA1')
C
C     Set up index information for PRPMO(N2PMOX,NZ)
C
      IPMO(1) = 0
      IPMO(2) = NPMO(1)
      NPMOT = NPMO(1)
      IF (NFSYM.EQ.2) NPMOT = NPMOT + NPMO(2)
      IOFF = 0
      DO J = 1,NFSYM
        DO I = 1,NFSYM
          I2PMOX(I,J) = IOFF + IPMO(I)
        END DO
        IOFF = IOFF + NPMO(J)*NPMOT
      END DO
C
C     Extract symmetry information
C     =============================
C
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)

C ... check for anti-/symmetry of the property matrix/MI

      IF (.NOT.(ITIM.EQ.1.OR.ITIM.EQ.-1)) THEN
       WRITE(LUPRI,'(/A,I3,1X,2A)')
     & '  PRPMA1: Property operator no.',INDXPR,PRPNAM(INDXPR),
     & ' is NOT (anti-)hermitian!'
       CALL QUIT('PRPMA1: Non-hermitian property!')
      ENDIF

C
C
C     Get property matrix in AO-basis
C     ===============================
C     
      CALL LSET(NZ,.TRUE.,FIRST)
      CALL PRPMAO(LU1INT,INDXPR,PTRI,OP1INT,FIRST,WORK(KFREE),
     &           LFREE,IPRINT)   
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(PRPAO(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,ITIM)
          IF    (IH.EQ.1) THEN
            CALL DSPTSI(NTBAS(0),PTRI(1,IZ),PRPAO(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),PTRI(1,IZ),PRPAO(1,IZ))
          ENDIF
        ENDIF
      ENDDO
      CALL BUTOBS(PRPAO,NZ,WORK(KFREE),LFREE)
C
C     Transform to MO-basis only the blocks of correct fermion symmetry
C     =================================================================
C         
      CALL DZERO(PRPMO,N2PMOX*NZ)
      IOPSY = JBTOF(IREP,1)
      DO I1 = 1,NFSYM
         I2  = MOD(I1+IOPSY,2) + 1
           IF(SPINFR.AND.((ITIM.EQ.1).OR.NOSPIN)) THEN
           CALL QBTRANS(IREP,'AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NPMO(I1),NPMO(I2),
     &           PRPAO(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
     &           PRPMO(I2PMOX(I1,I2)+1,1),NPMOT,NPMOT,NZ,IPQTOQ(1,IREP),
     &           CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           IBEIG(IORB(I1)+1),
     &           CMO(IPCMO(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           IBEIG(IORB(I2)+1),
     &           WORK(KFREE),LFREE,IPRINT)
         ELSE
           CALL QTRANS('AOMO','S',D0,
     &           NFBAS(I1,0),NFBAS(I2,0),NPMO(I1),NPMO(I2),
     &           PRPAO(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),NZ,
     &           IPQTOQ(1,IREP),
     &           PRPMO(I2PMOX(I1,I2)+1,1),NPMOT,NPMOT,NZ,IPQTOQ(1,IREP),
     &           CMO(IPCMO(I1)+1),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &           CMO(IPCMO(I2)+1),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRINT)
         ENDIF
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Output from PRPMA1',-1)
        WRITE(LUPRI,'(A,A16)') 
     &    '*MO-matrix of property ',PRPNAM(INDXPR)
        CALL PRQMAT(PRPMO,NPMOT,NPMOT,NPMOT,NPMOT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
        IF(IPRINT.GE.10) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '*AO-matrix of property ',PRPNAM(INDXPR)
        CALL PRQMAT(PRPAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDIF
C
      CALL QEXIT('PRPMA1')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gpget */
      SUBROUTINE GPGET(INDXPR,GP,JXPAR,NXPAR,CMO,IBEIG,ITYP,
     &                 WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C
C              | (q+|Q) |   | -<0|[q ,Q]|Q> |   | e[1]_Q* |
C     E[1]_Q = | (q |Q) | = |  <0|[q+,Q]|Q> | = | e[1]_Q  |
C
C     In RPA q = a(+)m, where a is occupied and m secondary.
C     This routine will construct e[1]_Q*
C
C     e[1]_Q*(ma) = -<0|[a(+)m,Q]|Q> 
C                 = -<0|a(+)[m,Q]|0> - <0|[a(+),Q]m|0>
C                 = -<0|a(+)[m,Q]|0> - <0|,m(+)[a,Q]m|0>(*) = -Q_(ma)
C
C     GRDORB employes the orbital rotation index JXOP(a,m).
C     Indices are therefore permuted...
C
C     Written by T.Saue - July 1996
C     Last revision: Jul 10 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbibn.h"
      DIMENSION GP(*),CMO(*),IBEIG(*),JXPAR(*),WORK(*)
C
      CALL QENTER('GPGET')
      KFRSAV = KFREE
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
      IF(ITYP.EQ.JBENDX) THEN
        IF(NXPAR.GE.0) THEN
          CALL MEMGET('REAL',KPRPMO,N2BBASXQ,WORK,KFREE,LFREE)
          CALL GPOGET(INDXPR,GP,JXPAR,NXPAR,WORK(KPRPMO),CMO,
     &                IBEIG,WORK,KFREE,LFREE,IPRINT)
          CALL MEMREL('GPGET.e-e',WORK,1,KFRSAV,KFREE,LFREE)
          IF(IPRINT.GE.4) THEN
            CALL HEADER('Output from GPGET',-1)
            IREP = ISYM - 1
            WRITE(LUPRI,'(A,A16)') 
     &      '* Gradient (e-e) of property ',PRPNAM(INDXPR)
            CALL PRQMAT(GP,NXPAR,1,NXPAR,1,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
      ELSEIF(ITYP.EQ.JBPNDX) THEN
        IF(NXPAR.GE.0) THEN
          CALL MEMGET('REAL',KPRPMO,N2BBASXQ,WORK,KFREE,LFREE)
          CALL GPOGET(INDXPR,GP,JXPAR,NXPAR,WORK(KPRPMO),CMO,
     &                IBEIG,WORK,KFREE,LFREE,IPRINT)
          CALL MEMREL('GPGET.e-p',WORK,1,KFRSAV,KFREE,LFREE)
          IF(IPRINT.GE.4) THEN
            CALL HEADER('Output from GPGET',-1)
            IREP = ISYM-1
            WRITE(LUPRI,'(A,A16)') 
     &      '* Gradient (e-p) of property ',PRPNAM(INDXPR)
              CALL PRQMAT(GP,NXPAR,1,NXPAR,1,NZ,
     &                    IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
      ELSEIF(ITYP.EQ.JBCNDX) THEN
        WRITE(LUPRI,'(A)') 
     &    'GPGET ERROR: Configurational part not yet written !'
        CALL QUIT('GPGET: Configurational part not yet written !')
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'GPGET ERROR: ',
     &    'unknown ITYP = ',ITYP
        CALL QUIT('GPGET: Unknown ITYP !')
      ENDIF
C
C
      CALL QEXIT('GPGET')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GPOGET */
      SUBROUTINE GPOGET(INDXPR,GPO,JXOP,NXOP,PRPMO,CMO,
     &                  IBEIG,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Construct orbital part of gradient for property INDXPR 
C     listed in /XCBXPR/
C
C     The property gradient is defined by
C
C              | (q+|Q) |   | -<0|[q ,Q]|Q> |   | e[1]_Q* |
C     E[1]_Q = | (q |Q) | = |  <0|[q+,Q]|Q> | = | e[1]_Q  |
C
C     In RPA q = a(+)m, where a is occupied and m secondary.
C     This routine will construct e[1]_Q*
C
C     e[1]_Q*(ma) = -<0|[a(+)m,Q]|Q> 
C                 = -<0|a(+)[m,Q]|0> - <0|[a(+),Q]m|0>
C                 = -<0|a(+)[m,Q]|0> - <0|,m(+)[a,Q]m|0>(*) = -Q_(ma)
C
C     GRDORB employes the orbital rotation index JXOP(a,m).
C     Indices are therefore permuted...
C
C     Written by T.Saue July 1996
C     Last revision: Feb 17 2000 hjaaj (new parameter list for PRPMAT)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
C
      DIMENSION GPO(*),PRPMO(NORBT,NORBT,NZ),CMO(*),
     &          IBEIG(*),JXOP(NXOP),WORK(*)
C
      CALL QENTER('GPOGET')
C
C     Get property matrix and scale by factor -1
C     ==========================================
C
      ISYM  = IPRPSYM(INDXPR)
      IOPSY = JBTOF(ISYM-1,1)
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL PRPMAT(INDXPR,IOPSY,PRPMO,CMO,IBEIG,ICMOQ,NORB,
     &            WORK,KFREE,LFREE,IPRINT)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Hermicity of operator
C     =====================
C
      ITIM = IPRPTIM(INDXPR)
      IREP = IPRPSYM(INDXPR) - 1
      IQ   = IPQTOQ(1,IREP)
      IH   = IHQMAT(IQ,ITIM)
C
C
C
C     Construct gradient - electronic part
C     ====================================
C
      KFRSAV = KFREE
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHXQ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRINT,.FALSE.,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
      CALL DZERO(GPO,NZ*NXOP)
      CALL PRPORB(GPO,JXOP,NXOP,PRPMO,NORBT,WORK(KDV),NZ,IH)
C
      CALL MEMREL('GPOGET',WORK,1,KFRSAV,KFREE,LFREE)
c       no scaling
      CALL DSCAL(NZ*NXOP,DM1,GPO,1)
      CALL QEXIT ('GPOGET')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prporb */
      SUBROUTINE PRPORB(QVEC,JXOP,NXOP,QMAT,NDIM,DV,NZ,IH)
C***********************************************************************
C
C     Construct 1-eletron property gradient
C
C     Input:
C        JXOP    - list of orbital rotations
C        QMAT    - one-electron property matrix
C        UDV     - active density matrix (symmetry packed)
C        IH      - hermicity
C           IH = 1: g(k,l) = sum_m Dlq Qkq - ( sum_m Dmk Qml )^*
C           IH = 2: g(k,l) = sum_m Dlq Qkq + ( sum_m Dmk Qml )^*
C
C     Written by J. Thyssen May 7 1999
C     Last revision 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
      PARAMETER ( D2 = 2.0D00, D1 = 1.0D00, DM1 = -1.0D00 )
C
#include "dcborb.h"
#include "mx2orb.h"
#include "dcbidx.h"
C
      DIMENSION QVEC(NXOP,NZ),JXOP(2,NXOP),QMAT(NDIM,NDIM,NZ),
     &          DV(*)
C
C
      IF (IH .EQ. 1) THEN
         DS = D1
      ELSE IF (IH. EQ. 2) THEN
         DS = DM1
      END IF
C
      DO IG = 1,NXOP
         K = JXOP(1,IG)
         L = JXOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         ISYMK = IFSMO(K)
         ISYML = IFSMO(L)
C
C        First index:
C
         IF (ITYPK .EQ. JTINAC) THEN
C
C        ...inactive
C
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) + D2 * QMAT(L,K,IZ)
            END DO
C
         ELSE IF (ITYPK .EQ. JTACT) THEN
C
C        ...active
C           TODO : this code is written for MCSCF, it is not tested
C                  for open shell RPA /Trond 2001, hjaaj July 2002
C
            KDV = K - (NPSH(ISYMK) + NISH(ISYMK)) - IORB(ISYMK)
            JDV = I2ASHT(ISYMK) + (KDV-1)*NASH(ISYMK) + 1
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) 
     &            + DDOT(NASH(ISYMK),DV(JDV),1,QMAT(L,K,IZ),1)
               JDV = JDV + N2ASH(ISYMK)
            END DO
C
         END IF
C
C        Second index:
C
         IF (ITYPL .EQ. JTACT) THEN
C
C           ...active
C           TODO : this code is written for MCSCF, it is not tested
C                  for open shell RPA /Trond 2001, hjaaj July 2002
C
            KDV = L - (NPSH(ISYML) + NISH(ISYML)) - IORB(ISYML)
            JDV = I2ASHT(ISYML) + (KDV-1)*NASH(ISYML) + 1
            DO IZ = 1,NZ
               QVEC(IG,IZ) = QVEC(IG,IZ) 
     &            + DS * DDOT(NASH(ISYML),DV(JDV),1,QMAT(K,L,IZ),1)
               JDV = JDV + N2ASH(ISYML)
            END DO
C
         END IF
C

      END DO
      RETURN
      END
