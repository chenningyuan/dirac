C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dirac, a relativistic ab initio electronic structure program,
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cl2nz3 */
      SUBROUTINE CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C*****************************************************************************
C
C     Transform from 4INDXXXX0-integral-file format to (NZ,3) format.
C
C     Input:
C        VCL   : elements in (*,NCLASS) format
C        NELEM : # of elements
C        IREPIJ: irep of VCL
C       
C     Output:
C        VNZ3  : elements in (*,NZ,3) format
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"     
#include "dcbham.h"     
C
      DIMENSION VCL(NELEM,NZ,NZ,*)
      DIMENSION VNZ3(NELEM,NZ,3)
      DIMENSION IREPIJB(2),IREPKLB(2,*),IKLB(2)
C
#include "sdpre.h"
C
      CALL DZERO(VNZ3,NELEM*NZ*3)
C
      IF (SPINFR.OR.LEVYLE) THEN
C
C     Special code that takes out the quaternion phase factors
C     The resulting integrals have full 8-fold permutational
C     symmetry and can be fed into non-relativistic codes.
C     Only one class of integrals (all-unbarred) suffices for
C     this purpose.
C
      IREPABR = 0
      DO IREPAB = 0, NBSYM-1
         IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1
C
C           Zero out the spin-orbit coupling for electron 1
C
            IF (IREPAB.NE.IREPIJB(1)) GOTO 100
            DO I = 1, NELEM
               IF (IREPAB.NE.IREPKLB(1,IKLB(IREPIJ)+I)) THEN
C
C                  Zero out the spin-orbit coupling for electron 2
C
                   VCL(I,1,1,IREPABR) = 0.D0
               ELSE
C
C                  Take out quaternion phase factor by multiplying
C                  with the compound phase factor that was introduced
C                  in the consecutive multiplications
C
                   VCL(I,1,1,IREPABR) =
     &             VCL(I,1,1,IREPABR) *
     &             IREPIJB(2)*IREPKLB(2,IKLB(IREPIJ)+I)
               ENDIF
            ENDDO
            IFAC = 1
            CALL DAXPY (NELEM,SDPRE(IFAC),
     &                  VCL(1,1,1,IREPABR),1,
     &                  VNZ3(1,1,1),1)
 100        CONTINUE
         ENDIF
      ENDDO
C
      ELSE ! spinfree
      DO ICLASS = 1, 3
         IREPABR = 0
         DO IREPAB = 0, NBSYM - 1
            IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
               IREPABR = IREPABR + 1
               DO IZ1 = 1, NZ
                  IQ1 = IPQTOQ(IZ1,IREPAB)
                  DO IZ2 = 1,NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
                     IFAC = IFACNZ3(ICLASS,1,IQ1,IQ2)
                     IF (IFAC .NE. 0) THEN
                        CALL DAXPY(NELEM,SDPRE(IFAC),
     &                       VCL(1,IZ2,IZ1,IREPABR),1,
     &                       VNZ3(1,1,ICLASS),1)
                     END IF
                     IF (NZ .GT. 1) THEN
                        IFAC = IFACNZ3(ICLASS,2,IQ1,IQ2)
                        IF (IFAC .NE. 0) THEN
                           CALL DAXPY(NELEM,SDPRE(IFAC),
     &                          VCL(1,IZ2,IZ1,IREPABR),1,
     &                          VNZ3(1,2,ICLASS),1)
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
      IF (NZ .EQ. 4) THEN
         DO ICLASS = 1, 3
            IREPABR = 0
            DO IREPAB = 0, NBSYM - 1
               IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ1 = 1, NZ
                     IQ1 = IPQTOQ(IZ1,IREPAB)
                     DO IZ2 = 1,NZ
                        IQ2 = IPQTOQ(IZ2,IREPAB)
                        IFAC = IFACNZ3(ICLASS,3,IQ1,IQ2)
                        IF (IFAC .NE. 0) THEN
                           CALL DAXPY(NELEM,SDPRE(IFAC),
     &                          VCL(1,IZ2,IZ1,IREPABR),1,
     &                          VNZ3(1,3,ICLASS),1)
                        END IF
                        IFAC = IFACNZ3(ICLASS,4,IQ1,IQ2)
                        IF (IFAC .NE. 0) THEN
                           CALL DAXPY(NELEM,SDPRE(IFAC),
     &                          VCL(1,IZ2,IZ1,IREPABR),1,
     &                          VNZ3(1,4,ICLASS),1)
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
C
      ENDIF ! spinfree
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck x2uu */
      SUBROUTINE X2UU(HX,H2XY)
C*****************************************************************************
C
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "dgroup.h"     
#include "dcbopt.h"     
#include "dcborb.h"     
#include "pgroup.h"     
#include "dcbidx.h"
C
C
      DIMENSION HX(NORBT,NORBT,NZ,3)
      DIMENSION H2XY(NASHT,NASHT,NZ,3)
C
      DO I3 = 1, 3
         DO IZ = 1, NZ
            DO IU = 1, NASHT
               IUG = IDXU2G(IU)
               DO IV = 1, NASHT
                  IVG = IDXU2G(IV)
                  H2XY(IU,IV,IZ,I3) = HX(IUG,IVG,IZ,I3)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck x2gu */
      SUBROUTINE X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C*****************************************************************************
C
C     Make (ug) block integrals from (gu) block.
C
C     Input:
C        GUINT: (gu) integrals
C        NORBI,NASHI: dimensions
C       
C     Output:
C        UGINT: (ug) integrals.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "dgroup.h"     
#include "dcbopt.h"     
#include "dcborb.h"     
#include "pgroup.h"     
#include "dcbidx.h"
C
C
      DIMENSION HX(NORBT,NORBT,NZ,3)
      DIMENSION H2XY(NORBT,NASHT,NZ,3)
      DIMENSION H2YX(NORBT,NASHT,NZ,3)
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3)
      LOGICAL   ONLYXY, DOH2AC, DOFQ
C
      IF (DOH2AC) THEN
C
         DO I3 = 1, 3
            DO IZ = 1, NZ
               DO IU = 1, NASHT
                  IUG = IDXU2G(IU)
                  DO IV = 1, NASHT
                     IVG = IDXU2G(IV)
                     H2AC(IU,IV,IUU,IZ,I3) = HX(IUG,IVG,IZ,I3)
                  END DO
               END DO
            END DO
         END DO
C
C        Print section
C
         IF (IPROPT .GE. 30 .AND. DOH2AC) THEN
            WRITE(LUPRI,'(1X,A,I3,A)') '(X2GU) H2AC(*,*,',IUU,'): '
            DO I3 = 1, 3
               DO IZ = 1, NZ
                  WRITE(LUPRI,'(1X,A,I3,A,I3,A)')
     &                 'Matrix no. (',IZ,',',I3,'): '
                  CALL OUTPUT(H2AC(1,1,IUU,IZ,I3),1,NASHT,1,NASHT,
     &                 NASHT,NASHT,1,LUPRI)
               END DO
            END DO
         END IF
      END IF
C
C     Collect (gu) and (ug) elements.
C
      IF (DOFQ) THEN
         CALL DZERO(H2XY, NORBT * NASHT * NZ * 3)
         IF (.NOT. ONLYXY) CALL DZERO(H2YX, NORBT * NASHT * NZ * 3)
C
         DO I3 = 1, 3
            DO IZ = 1, NZ
               DO IG = 1, NORBT
                  DO IU = 1, NASHT
                     IUG = IDXU2G(IU)
                     H2XY(IG,IU,IZ,I3) = HX(IG,IUG,IZ,I3)
                     IF (.NOT. ONLYXY)
     &                    H2YX(IG,IU,IZ,I3) = HX(IUG,IG,IZ,I3)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gu2ug */
      SUBROUTINE GU2UGT(GUINT,UGTINT,IREPIJ,NORBI,NASHI)
C*****************************************************************************
C
C     Make (ug) block integrals from (gu) block.
C
C     Input:
C        GUINT: (gu) integrals
C        NORBI,NASHI: dimensions
C       
C     Output:
C        UGINT: (ug) integrals.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dgroup.h"     
#include "pgroup.h"     
C
      DIMENSION GUINT(NORBI,NASHI,NZ,NZ,*)
      DIMENSION UGTINT(NORBI,NASHI,NZ,NZ,*)
C
#include "sdpre.h"
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DCOPY(NORBI*NASHI*NCLASS,GUINT,1,UGTINT,1)
C
C
C     All matrices are either symmetric or anti-symmetric.
C
C     p = (a b)^T, q = (c d)^T
C
C     (pq) = p^{\dagger}q = (a^* b^*) (c d)^T = a^* c  + b^* d
C     (qp) = q^{\dagger}p = (c^* d^*) (a b)^T = c^* a  + d^* d = (pq)^*
C
C     (pQ) = p^{\dagger}Q = (a^* b^*) (-d^* c^*)^T 
C                         = - a^* d^* + b^* c^*
C     (qP) = q^{\dagger}P = (c^* d^*) (-b^* a^*)^T
C                         = - c^* b^* + d^* a^* = - (pQ)
C
C     Hence,
C
C     Re(pq) is symmetric,
C     Im(pq), Re(pQ), Im(pQ) are all anti-symmetric.
C
      IREPABR = 0
      DO IREPAB = 0, NBSYM - 1
         IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1
C           
            DO IZ2 = 1,NZ
C              
               IQ2 = IPQTOQ(IZ2,IREPAB)
C              
               IF (IQ2 .EQ. 1) THEN
C                 ...symmetric (scale block with 1)
C                 no-op
               ELSE
C                 ...anti symmetric (scale block with -1)
                  CALL DSCAL(NZ*NASHI*NORBI,DM1,
     &                 UGTINT(1,1,1,IZ2,IREPABR),1)
               END IF
            END DO
         END IF
      END DO
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ifacnz3 */
      FUNCTION IFACNZ3(I3,IZ,IQ1,IQ2)
C*****************************************************************************
C
C     Factors for integrals.
C
C     Written by Joern Thyssen - Nov 15 2000
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
C
      COMMON /SMOELF/ IFAC(4,3,4,4)
C
      LOGICAL FIRST
      DATA FIRST / .TRUE. /
      SAVE FIRST
C
      IF (FIRST) THEN
C        Initialize IFAC array
         FIRST = .FALSE.
         CALL IZERO(IFAC,3*4*4*4)
C
C        IZ = 1
C
C        Re(pq|rs) = Re(rs|pq) = (Re(rs)|Re(pq)) - (Im(pq)|Im(pq))
C                                    +(1,1)            -(2,2)
C
         IFAC(1,1,1,1) = +1
         IFAC(1,1,2,2) = -1
C
C        Re(pQ|rS) = Re(rS|pQ) = (Re(rS)|Re(pQ)) - (Im(rS)|Im(pQ))
C                                    +(3,3)            -(4,4)
C
         IFAC(1,2,3,3) = +1
         IFAC(1,2,4,4) = -1
C
C        Re(Pq|rS) = Re(rS|Pq) = - (Re(rS)|Re(pQ)) - (Im(rS)|Im(pQ))
C                                      -(3,3)            -(4,4)
C
         IFAC(1,3,3,3) = -1
         IFAC(1,3,4,4) = -1
C
         IF (NZ .GE. 2) THEN
C
C           Im(pq|rs) = Im(rs|pq) = (Re(rs)|Im(pq)) + (Im(rs)|Re(pq)) + 
C                                       +(1,2)            +(2,1)
C
            IFAC(2,1,1,2) = +1
            IFAC(2,1,2,1) = +1
C
C           Im(pQ|rS) = Im(rS|pQ) = (Re(rS)|Im(pQ)) + (Im(rS)|Re(pQ))
C                                       +(3,4)            +(4,3)
C
            IFAC(2,2,3,4) = +1
            IFAC(2,2,4,3) = +1
C
C           Im(Pq|rS) = Im(rS|Pq) = - (Re(rS)|Im(pQ)) + (Im(rS)|Im(pQ))
C                                         +(3,4)            +(1,2)
C
            IFAC(2,3,3,4) = +1
            IFAC(2,3,4,3) = -1
         END IF
         IF (NZ .EQ. 4) THEN
C
C           Re(pq|rS) = Re(rS|pq) = (Re(rS)|Re(pq)) - (Im(rS)|Im(pq))
C                                       +(3,1)            -(4,2)
C
            IFAC(3,1,3,1) = +1
            IFAC(3,1,4,2) = -1
C
C           Im(pq|rS) = Im(rS|pq) = (Re(rS)|Im(pq)) + (Im(rS)|Re(pq))
C                                       +(3,2)            +(4,1)
C
            IFAC(4,1,3,2) = +1
            IFAC(4,1,4,1) = +1
C
C           Re(pQ|rs) = Re(rs|pQ) = (Re(rs)|Re(pQ)) - (Im(rs)|Im(pQ))
C                                       +(1,3)            -(2,4)
C
            IFAC(3,2,1,3) = +1
            IFAC(3,2,2,4) = -1
C
C           Im(pQ|rs) = Im(rs|pQ) = (Re(rs)|Im(pQ)) + (Im(rs)|Re(pQ))
C                                       +(1,4)            +(2,3)
C
            IFAC(4,2,1,4) = +1
            IFAC(4,2,2,3) = +1
C
C           Re(Pq|rs) = Re(rs|Pq) = - (Re(rs)|Re(pQ)) - (Im(rs)|Im(pQ))
C                                         -(1,3)            -(2,4)
C
            IFAC(3,3,1,3) = -1
            IFAC(3,3,2,4) = -1
C
C           Im(Pq|rs) = Im(rs|Pq) = (Re(rs)|Im(pQ)) - (Im(rs)|Re(pQ))
C                                       +(1,4)            -(2,3)
C
            IFAC(4,3,1,4) = +1
            IFAC(4,3,2,3) = -1
         END IF
      END IF
C
      IFACNZ3 = IFAC(IZ,I3,IQ1,IQ2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck btch2f */
      SUBROUTINE BTCH2F(HT,IZ2,IREPABR,HX,IREPIJ,NKL)
C*****************************************************************************
C
C     Unpack packed integrals from HB to HF.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C     
      DIMENSION HT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ)
C
      CALL DZERO(HX,N2ORBXQ)
C
      IF (.NOT. OPT_NOPFQ) THEN
C        
         DO IZ1 = 1, NZ
            IKL = 1
            DO IREPL = 1, NFSYM
               IREPK = MOD (IREPL + IREPIJ, 2) + 1
               DO IROW = 1, NORB(IREPL)
                  CALL DCOPY(NORB(IREPK),HT(IKL,IZ1,IZ2,IREPABR),1,
     &                 HX(1+IORB(IREPK),IROW+IORB(IREPL),IZ1),1)
                  IKL = IKL + NORB(IREPK)
               END DO
            END DO
         END DO
C
      ELSE
C        
         DO IZ1 = 1, NZ
            IKL = 1
            DO IREPL = 1, NFSYM
               IREPK = MOD (IREPL + IREPIJ, 2) + 1
               DO IROW = 1, NESH(IREPL)
                  CALL DCOPY(NESH(IREPK),HT(IKL,IZ1,IZ2,IREPABR),1,
     &                 HX(1+IORB(IREPK)+NPSH(IREPK),
     &                    IROW+IORB(IREPL)+NPSH(IREPL),IZ1),1)
                  IKL = IKL + NESH(IREPK)
               END DO
            END DO
         END DO
C
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intb2x */
      SUBROUTINE INTB2X(HT,IZ2,IREPABR,IREPAB,HX,IREPIJ,NKL)
C*****************************************************************************
C
C     Unpack packed integrals from HT to HX.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C     
      DIMENSION HT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ,NZ,*)
C
      DIMENSION JCOL(2), NCOL(2)
      DIMENSION JROW(2), NROW(2)
      LOGICAL DOFULL
C
      IF (JTRLVL .EQ. 0) THEN
C
C        HT is symmetry-packed (uu) integrals
C
C        (uu) is a ``symmetric'' index so don't symmetrize matrix.
C
         DOFULL = .FALSE.
C
      ELSE IF (JTRLVL .EQ. 1) THEN
C
C        HT is symmetry-packed (gu) or (eu) integrals.
C
C        Do symmetrization:
C
         DOFULL = .TRUE.
C
      ELSE 
C
C        HT is symmetry-packed (gg), (ge), or (ee) integrals.
C
         IF (NZXOPP .EQ. 0) THEN
C
C           (ee) integrals.
C
C           (ee) is a ``symmetric'' index so don't symmetrize matrix.
C           
            DOFULL = .FALSE.
C
         ELSE IF (.NOT. OPT_NOPFQ) THEN
C
C           (gg) integrals.
C
C           (gg) is a ``symmetric'' index so don't symmetrize matrix.
C           
            DOFULL = .FALSE.
C
C
         ELSE
C
C           (ge) integrals.
C
            DOFULL = .TRUE.
C
         END IF
C
      END IF
C
C
      DO IZ1 = 1, NZ
         IKL = 0
         IQ1 = IPQTOQ(IZ1,IREPAB)
         DO IREPL = 1, NFSYM
            IREPK = MOD (IREPL + IREPIJ, 2) + 1
            DO ICOL = 1, NIDX4(IREPL)
               JJ = IDXT2G(ICOL+IIDX4(IREPL),4)
               DO IROW = 1, NIDX3(IREPK)
                  II = IDXT2G(IROW+IIDX3(IREPK),3)
                  IKL = IKL + 1
                  HX(II,JJ,IZ1,IZ2,IREPABR) =
     &                 HT(IKL,IZ1,IZ2,IREPABR)
               END DO
            END DO
C
            IF (DOFULL) THEN
C
C              Calculate full matrix
C
               IF (IQ1 .EQ. 1) THEN
                  DO ICOL = 1, NIDX3(IREPL)
                     JJ = IDXT2G(ICOL+IIDX3(IREPL),3)
                     DO IROW = 1, NIDX4(IREPK)
                        II = IDXT2G(IROW+IIDX4(IREPK),4)
                        HX(II,JJ,IZ1,IZ2,IREPABR) =
     &                       HX(JJ,II,IZ1,IZ2,IREPABR)
                     END DO
                  END DO
               ELSE
                  DO ICOL = 1, NIDX3(IREPL)
                     JJ = IDXT2G(ICOL+IIDX3(IREPL),3)
                     DO IROW = 1, NIDX4(IREPK)
                        II = IDXT2G(IROW+IIDX4(IREPK),4)
                        HX(II,JJ,IZ1,IZ2,IREPABR) =
     &                       - HX(JJ,II,IZ1,IZ2,IREPABR)
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GU2F(GUINT,IZ2,IREPAB,IREPABR,HX,IREPIJ,NKL)
C*****************************************************************************
C
C     Unpack packed integrals from HB to HF.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C     
      DIMENSION GUINT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ)
C
      CALL DZERO(HX,N2ORBXQ)
C
      DO IZ1 = 1, NZ
         IQ1 = IPQTOQ(IZ1,IREPAB)
         IF (IQ1.EQ.1) THEN
            IFAC = 1
         ELSE
            IFAC = -1
         END IF
C
         IKL = 1
         DO IREPL = 1, NFSYM
            IREPK = MOD (IREPL + IREPIJ, 2) + 1
C
C           Copy (gu) integrals into matrix:
C
            DO ICOL = 1, NASH(IREPL)
               CALL DCOPY(NORB(IREPK),GUINT(IKL,IZ1,IZ2,IREPABR),1,
     &              HX(1+IORB(IREPK),
     &              ICOL+IORB(IREPL)+NPSH(IREPL)+NISH(IREPL),IZ1),1)
               IKL = IKL + NORB(IREPK)
            END DO
C
C           Full matrix
C
            DO I = 1, NASH(IREPK)
               DO J = 1, NORB(IREPL)
                  HX(I+IORB(IREPK)+NPSH(IREPL)+NISH(IREPL),
     &                 J+IORB(IREPL),IZ1) =
     &                 IFAC *
     &                 HX(J+IORB(IREPL),
     &                 I+IORB(IREPK)+NPSH(IREPL)+NISH(IREPL),IZ1)
               END DO
            END DO
C
         END DO
C
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck f2btch */
      SUBROUTINE F2BTCH(HX,HT,IZ2,IREPABR,IREPIJ,NKL)
C*****************************************************************************
C
C     Pack integrals from HT into HX.
C     I.e. go from dimensioning N2ORBTQ to N2ORBXQ.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C     
      DIMENSION HT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ)
C
C
      IF (.NOT. OPT_NOPFQ) THEN
C        
         DO IZ1 = 1, NZ
            CALL DZERO(HT(1,IZ1,IZ2,IREPABR),NKL)
            IKL = 1
            DO IREPL = 1, NFSYM
               IREPK = MOD (IREPL + IREPIJ, 2) + 1
               DO IROW = 1, NORB(IREPL)
                  CALL DCOPY(NORB(IREPK),
     &                 HX(1+IORB(IREPK),IROW+IORB(IREPL),IZ1),1,
     &                 HT(IKL,IZ1,IZ2,IREPABR),1)
                  IKL = IKL + NORB(IREPK)
               END DO
            END DO
         END DO
C
      ELSE
C        
         DO IZ1 = 1, NZ
            CALL DZERO(HT(1,IZ1,IZ2,IREPABR),NKL)
            IKL = 1
            DO IREPL = 1, NFSYM
               IREPK = MOD (IREPL + IREPIJ, 2) + 1
               DO IROW = 1, NESH(IREPL)
                  CALL DCOPY(NESH(IREPK),
     &                 HX(1+IORB(IREPK)+NPSH(IREPK),
     &                    IROW+IORB(IREPL)+NPSH(IREPL),IZ1),1,
     &                 HT(IKL,IZ1,IZ2,IREPABR),1)
                  IKL = IKL + NESH(IREPK)
               END DO
            END DO
         END DO
C
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck f2gu */
      SUBROUTINE F2GU(HX,HT,IZ2,IREPABR,IREPIJ,NKL)
C*****************************************************************************
C
C     Pack full (gg) matrix into symmetry packed (gu) matrix.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C     
      DIMENSION HT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ)
C
C
      DO IZ1 = 1, NZ
         CALL DZERO(HT(1,IZ1,IZ2,IREPABR),NKL)
         IKL = 1
         DO IREPL = 1, NFSYM
            IREPK = MOD (IREPL + IREPIJ, 2) + 1
            DO IROW = 1, NASH(IREPL)
               CALL DCOPY(NORB(IREPK),
     &              HX(1+IORB(IREPK),
     &              IROW+IORB(IREPL)+NPSH(IREPL)+NISH(IREPL),IZ1),1,
     &              HT(IKL,IZ1,IZ2,IREPABR),1)
               IKL = IKL + NORB(IREPK)
            END DO
         END DO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck blkzero */
      SUBROUTINE BLKZERO(J1MO,N1MO,J2MO,N2MO,AMAT,NRDIM,NCDIM,
     &              NZ,NFSYM,IREPIJ)
C*****************************************************************************
C
C
C     Zero out block defined by the offsets J{1,2}MO, index lengths
C     N{1,2}MO. Only zero out the quaternion blocks with symmetry IREPIJ.
C
C
C*****************************************************************************
#include "implicit.h"
C
#include "consts.h"
C
      DIMENSION AMAT(NRDIM,NCDIM,NZ)
      DIMENSION J1MO(NFSYM), N1MO(NFSYM)
      DIMENSION J2MO(NFSYM), N2MO(NFSYM)
C
      DO IZ = 1, NZ
         DO JFSYM = 1, NFSYM
            IFSYM = MOD ( JFSYM + IREPIJ, 2 ) + 1
            DO J = J2MO(JFSYM), J2MO(JFSYM) + N2MO(JFSYM) - 1
               DO I = J1MO(IFSYM), J1MO(IFSYM) + N1MO(IFSYM) - 1
                  AMAT(I,J,IZ) = D0
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck blkunit */
      SUBROUTINE BLKUNIT(JMO,NMO,AMAT,NRDIM,NCDIM,NZ,NFSYM)
C*****************************************************************************
C
C
C     Zero out block defined by the offsets J{1,2}MO, index lengths
C     N{1,2}MO. Only zero out the quaternion blocks with symmetry IREPIJ.
C
C
C*****************************************************************************
#include "implicit.h"
C
#include "consts.h"
C
      DIMENSION AMAT(NRDIM,NCDIM,NZ)
      DIMENSION JMO(NFSYM), NMO(NFSYM)
C
C     First zero matrix
C
      CALL BLKZERO(JMO,NMO,JMO,NMO,AMAT,NRDIM,NCDIM,NZ,NFSYM,1)
C
C     Insert 1's in diagonal
C
      DO JFSYM = 1, NFSYM
         DO J = JMO(JFSYM), JMO(JFSYM) + NMO(JFSYM) - 1
            AMAT(J,J,1) = D1
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rprci */
      SUBROUTINE RPRCI(CIVEC,NZCONF,NZ,CIPRG,THRPCI,LUPRI)
C*****************************************************************************
C
C     Print ci vector.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Nov 23 2000
C     July 2002 hjaaj: nothing CIPRG specific here!
C
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION CIVEC(NZCONF,*)
      DIMENSION IP(4)
      CHARACTER CIPRG*(*)
C
chj   IF (CIPRG .EQ. 'GOSCI') THEN
         IF (NZ .EQ. 1) THEN
            WRITE(LUPRI,9000)
            DO I = 1, NZCONF
               CIWGHT = CIVEC(I,1)**2
               IF (CIWGHT .GT. THRPCI) WRITE(LUPRI,9001)
     &              I,CIVEC(I,1),CIWGHT
            END DO
         ELSE
            WRITE(LUPRI,9002)
            DO I = 1, NZCONF
               CIWGHT = CIVEC(I,1)**2 + CIVEC(I,2)**2
               IF (CIWGHT .GT. THRPCI) WRITE(LUPRI,9003)
     &              I,(CIVEC(I,IZ),IZ=1,2),CIWGHT
            END DO
         END IF
 9000    FORMAT(/' Real CI vector:'/
     &        T5,'Det. no.',T20,'Real part of C_i',
     &        T45,'Weight'/
     &        1X,52('-'))
 9001    FORMAT(1P,T5,I5,T20,D20.12,T45,0P,F8.6)
 9002    FORMAT(/' Complex CI vector:'/
     &        T5,'Det. no.',T20,'Real part of C_i',
     &        T45,'Imag part of C_i',T70,'Weight'/
     &        1X,78('-'))
 9003    FORMAT(1P,T5,I5,T20,D20.12,T45,D20.12,T70,0P,F8.6)
         WRITE(LUPRI,'(/)')
chj   END IF
C
      RETURN
      END
C  /* Deck prvec */
      SUBROUTINE RPRVEC(NDIM,VEC,INCVEC,NZ,IPQ,THRESH,MAXLIN,LUOUT)
C
C 19-Aug-1989 hjaaj
C   print VEC(1:NDIM*INCVEC:INCVEC)
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector to be printed
C   INCVEC    : Increment between each element in vector VEC
C   NZ        : real, complex, or quaternion
C   IPQ       : index packed quaternion to quaternion
C   THRESH    : Print threshold for vector with unit norm
C               (if THRESH .lt. 0 then -THRESH is used without
C                renormalization).
C   MAXLIN    : max. lines of output with vector elements
C   LUOUT     : output unit
C
C
#include "implicit.h"
      DIMENSION VEC(NDIM*INCVEC,NZ)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (D1LOW = D1 - 1.0D-10, D1HGH = D1 + 1.0D-10)
      LOGICAL   VSCALE
      DIMENSION IPQ(4)
#include "dcbqal.h"
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &     'No print from PRVEC because NDIM .le. 0'
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &      'No print from PRVEC because INCVEC .le. 0'
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) RETURN
C
C
C
      C2NRM = DNORM2(NDIM*NZ,VEC,INCVEC)
      IF (THRESH .LE. D0) THEN
         THRES2 = -THRESH
         VSCALE = .FALSE.
      ELSE
         THRES2 =  THRESH * C2NRM
         VSCALE = (C2NRM .LT. D1LOW .OR. C2NRM .GT. D1HGH)
      END IF
      IF (VSCALE) THEN
         WRITE (LUOUT,'(/2A,1P,D12.4,/)')
     &      ' Print of vector elements (vector scaled to unit norm)',
     &      ' larger than',ABS(THRESH)
         SCALE = D1 / C2NRM
      ELSE
         WRITE (LUOUT,'(/A,1P,D12.4,/)')
     &      ' Print of vector elements larger than',ABS(THRESH)
         SCALE = D1
      END IF
C
      IPR   = 0
      IZER  = 0
      C2SUM = D0
C
      IF (VSCALE) THEN
         WRITE(LUOUT,9000) 
      ELSE
         WRITE(LUOUT,9001) 
      END IF
 9000 FORMAT(5X,'Element',3X,'Coefficient (scaled to unit norm)',
     &     /,1X,78('-'))
 9001 FORMAT(5X,'Element',3X,'Coefficient',/,1X,78('-'))
C
      DO 300 I = 1, NDIM
         NA = 1 + (I-1)*INCVEC
         DNNA = D0
         DO IZ = 1, NZ
            DNNA = DNNA + VEC(NA,IZ)**2
         END DO
         IF (SQRT(DNNA).LE.THRES2 .OR. IPR .GE. MAXLIN) THEN
            C2SUM = C2SUM + DNNA
            IF (DNNA .EQ. D0) IZER = IZER + 1
         ELSE
            IF (MOD(IPR,5) .EQ. 0 .AND. IPR .NE. 0) WRITE (LUOUT,'()')
            IPR = IPR + 1
            IF (VSCALE) THEN
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,50)
     &                 I,VEC(NA,1),
     &                 SCALE*VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,51)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,52)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            ELSE
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,60)I,VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,61)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,62)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            END IF
 50         FORMAT(1X,I10,3X,1P,D10.2,0P,'(',F10.6,')')
 51         FORMAT(1X,I10,3X,1P,D10.2,2X,D10.2,1X,A1,2X,
     &           0P,'(',F10.6,2X,F10.6,1X,A1,')')
 52         FORMAT(1X,I10,3X,1P,D10.2,2X,3(D10.2,1X,A1,2X),
     &           0P,'(',F10.6,2X,3(F10.6,1X,A1,2X),')')
 60         FORMAT(1X,I10,3X,F20.10)
 61         FORMAT(1X,I10,3X,F20.10,2X,F20.10,1X,A1,2X)
 62         FORMAT(1X,I10,3X,F20.10,2X,3(F20.10,1X,A1,2X))
         END IF
  300 CONTINUE
      IF (IPR .GE. MAXLIN) THEN
C
C     *** We have reached the print limit
C
         WRITE (LUOUT,910) IPR
      END IF
  910 FORMAT(/' Print limit of',I6,' elements has been reached.')
      C2SUM = SQRT(C2SUM)
      IF (IZER .EQ. NDIM) THEN
         WRITE (LUOUT,920) NDIM
      ELSE
         WRITE (LUOUT,930) NDIM,NDIM-IPR,IZER,C2SUM,C2NRM
      END IF
  920 FORMAT(/' Length of vector                      :',I10,
     *       /' All elements are zero.',/)
  930 FORMAT(/' Length of vector                      :',I10,
     *       /' Number of elements not printed        :',I10,
     *       /' Number of zero elements               :',I10,
     *       /' Total norm of coefficients not printed:',F10.6,
     *       /' (the coefficients are normalized to    ',F10.6,')',/)
      RETURN
C
C End of PRVEC
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rappbos */
      SUBROUTINE RAPPBOS(CMO,IBOSYM,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Calculates approximate boson symmetry of orbitals.
C
C     The "approximate" boson symmetry of an spinor is the max over all
C     symmetries,S, of:
C
C     <i|P(S)|i> = \sum_{\mu\nu} c_{\mu i}^* c_{\nu i} <\mu|P(S)|\nu>
C                = \sum_{\mu\nu\in S} c_{\mu i}^* c_{\nu i} S_{\mu\nu}
C
C     We use Luuk and Timo's APPBOS for the dirty work.
C
C     Input:
C        CMO : MO coefficients
C
C     Output:
C        IBOSYM: "approximate" boson symmetry of orbital.
C
C     Written by J. Thyssen - Jan 13 2001
C
C*****************************************************************************
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
#include "dcbbas.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
C
      DIMENSION CMO(*)
      DIMENSION IBOSYM(*)
      DIMENSION WORK(*)
C
      DIMENSION WEIGHT(8)
C
#include "sdpre.h"
#include "ibtfun.h"
#include "memint.h"
C
C     Matrix to hold orbitals.
C
      CALL MEMGET('REAL',KCMO,N2ORBTQ,WORK,KFREE,LFREE)
C
C     Construct T^{\dagger} S
C     -----------------------
C
C     Fetch T
C
      CALL MEMGET('REAL',KTMAT,N2TMT,WORK,KFREE,LFREE)
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','RAPPBOS')
      CALL READT(LUTMAT,N2TMT,WORK(KTMAT))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Calculate T^{\dagger} S
C
      CALL MEMGET('REAL',KSMOAO,N2TMT,WORK,KFREE,LFREE)
      CALL MKSMOAO(WORK(KSMOAO),WORK(KTMAT),IPRINT,WORK(KFREE),LFREE)
C
C     Transform orbitals to orthonormal basis
C     ---------------------------------------
C
      DO I = 1, NFSYM
c        write(6,*) 'cmo'
c        CALL PRQMAT(CMO(1+I2ORBT(I)),
c    &        NFBAS(I,0),NORB(I),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
c    &        LUPRI)
         CALL QGEMM(NORB(I),NORB(I),NFBAS(I,0),D1,
     &        'N','N',IPQTOQ(1,0),WORK(KSMOAO+I2TMT(I)),
     &        NORB(I),NFBAS(I,0),NZT,
     &        'N','N',IPQTOQ(1,0),CMO(1+ICMOQ(I)),
     &        NFBAS(I,0),NORB(I),NZ,
     &        D0,IPQTOQ(1,0),WORK(KCMO+I2ORBT(I)),
     &        NORB(I),NORB(I),NZ)
c        write(6,*) 'cmo in orthonormal basis'
c        CALL PRQMAT(WORK(KCMO+I2ORBT(I)),
c    &        NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
c    &        LUPRI)
      end do
C
C     Identify approximate boson symmetry
C     -----------------------------------
C
      CALL MEMREL('RAPPBOS.trans',WORK,KWORK,KTMAT,KFREE,LFREE)
C
      CALL MEMGET('INTE',KINDX,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRDUM,4*NORBT,WORK,KFREE,LFREE)
C
      DO I = 1, NFSYM
         CALL APPBOS(I,NORB(I),WORK(KINDX),IBOSYM(1+IORB(I)),
     &        WORK(KCMO+I2ORBT(I)),WORK(KRDUM))
      END DO
C
c     IF (IPRINT .GE. 3) THEN
         WRITE(LUPRI,'(/A)')
     &        ' (RAPPBOS) Approximate boson symmetries: '
         WRITE(LUPRI,9000) (REP(IBOSYM(I)),I=1,NORBT)
c     END IF
 9000 FORMAT(10X,15(A4))
C
      CALL MEMREL('RAPPBOS',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgmcmo */
      SUBROUTINE RGMCMO(CMO,WORK,LWORK)
C*****************************************************************************
C
C     Gram-Schmidt orthogonalize coefficients.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C
C*****************************************************************************
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"
C
      DIMENSION CMO(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('RGMCMO')
C
C     Matrix to hold orbitals.
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSAO),D1,WORK(KFREE),LFREE)
C
C     Gram-Schmidt orthogonalize coefficients
C     
      DO I = 1, NFSYM
         CALL MEMGET('REAL',KSCR1,NFBAS(I,0)*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KSCR2,NORB(I),WORK,KFREE,LFREE)
         WRITE(LUPRI,'(1X,A,I3)')
     &        'RGMCMO: Coefficients for symmetry ',I
         CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
C
         CALL QNORM(WORK(KSAO),CMO(1+ICMOQ(I)),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &        IRETUR,WORK(KSCR1),WORK(KSCR2))
C
         CALL MEMREL('RGMCMO.QNORM',WORK,KWORK,KSCR1,KFREE,LFREE)
C
         WRITE(LUPRI,'(1X,A,I3)')
     &        'RGMCMO: Coefficients for symmetry after GM',I
         CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
      END DO
      CALL MEMREL('RGMCMO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RGMCMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnorm */
      SUBROUTINE QNORM(S,VC,NBASI,NORBI,NZ,IP,IRETUR,SCR1,SCR2)
C*****************************************************************************
C
C     Gram-Schmidt orthogonalize coefficients.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C
C*****************************************************************************
#include "priunit.h"
#include "implicit.h"
C
#include "thrzer.h"
#include "consts.h"
      PARAMETER (THRRND = 0.9D0, THNORM = 1.0D-6)
C
      DIMENSION S(NBASI,NBASI)
      DIMENSION VC(NBASI,NORBI,NZ)
      DIMENSION IP(*)
      DIMENSION SCR1(NBASI,NZ)
      DIMENSION SCR2(NORBI)
C
      DO 20 I = 1, NORBI
         ITURN = 0
 1       ITURN = ITURN + 1
         IROUND = 0
C
C        Calculate S * VC(I)
C
         CALL QGEMM(NBASI,1,NBASI,D1,
     &        'N','N',IP,S,NBASI,NBASI,1,
     &        'N','N',IP,VC(1,I,1),NBASI,NORBI,NZ,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
         write(6,*) 'tnorm = ',tnorm
         IF (TNORM .LT. THRZER) THEN
            IRETUR = -1
            RETURN
         END IF
C
C        Normalize vector
C
         TNORM = D1 / SQRT(TNORM)
         DO IZ = 1, NZ
            CALL DSCAL(NBASI,TNORM,VC(1,I,IZ),1)
         END DO
C
C        Calculate S * VC(I)
C
         CALL QGEMM(NBASI,1,NBASI,D1,
     &        'N','N',IP,S,NBASI,NBASI,1,
     &        'N','N',IP,VC(1,I,1),NBASI,NORBI,NZ,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
C
C        Begin normalization loop
C
         DO J = 1, I - 1
            T = D0
            DO IZ = 1, NZ
               T = T +  DDOT(NBASI,VC(1,J,IZ),1,SCR1(1,IZ),1)
            END DO
            TNORM = TNORM - T * T
            SCR2(J) = -T
         END DO
C
         IF (TNORM .LT. THRZER) THEN
C           ... zero vector after orthogonalization
            IRETUR = -1
            RETURN
         END IF
C
         IF (TNORM .LT. THNORM) IROUND = IROUND + 1
C
         TNORM = D1 / SQRT (TNORM)
         DO J = 1, I - 1
            SCR2(J) = SCR2(J) * TNORM
         END DO
         SCR2(I) = TNORM
C
C        Replace VC(*,I)
C
C        SCR1 = \sum_{k=1,i} SCR2(k) * VC(*,k)
C
         CALL QGEMM(NBASI,1,I,D1,
     &        'N','N',IP,VC,NBASI,NORBI,NZ,
     &        'N','N',IP,SCR2,NORBI,1,1,
     &        D0,IP,SCR1,NBASI,1,NZ)
C
         DO IZ = 1, NZ
            CALL DCOPY(NBASI,SCR1(1,IZ),1,VC(1,I,IZ),1)
         END DO
C
         IF (ITURN .EQ. 1 .AND. IROUND .GT. 0) THEN
            write(6,*) 'INFO: second round in QNORM'
            GOTO 1
         END IF
         IF (IROUND .GT. 0) CALL QUIT('QNORM internal error')
 20   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck appfq */
      SUBROUTINE APPFQ(CTL,FQ,DV,FV,WORK,LWORK)
C*****************************************************************************
C
C     Approximate positronic part of FQ with (DV*FV).
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
C     
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION DV(NASHT,NASHT,NZ)
      DIMENSION FV(NORBT,NORBT,NZ)
      DIMENSION WORK(*)
      CHARACTER CTL*1
C
#include "memint.h"
C
      CALL QENTER('APPFQ')
C
      IF (IPROPT .GE. 20) THEN
         WRITE(LUPRI,'(1X,A)') '(APPFQ) FV'
         CALL PRQMAT(FV,NORBT,NORBT,NORBT,NORBT,
     &        NZ,IPQTOQ(1,0),LUPRI)
         WRITE(LUPRI,'(1X,A)') '        DV'
         CALL PRQMAT(DV,NASHT,NASHT,NASHT,NASHT,
     &        NZ,IPQTOQ(1,0),LUPRI)
         WRITE(LUPRI,'(1X,A)') '(APPFQ) FQ before add'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,
     &        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      IF (CTL .EQ. 'P') THEN
C        
C        We are called to generate the positronic part of FQ.
C
         DO I = 1, NFSYM
            CALL QGEMM(NPSH(I),NASH(I),NASH(I),D1,
     &           'N','N',IPQTOQ(1,0),
     &           FV(1+IORB(I),1+IORB(I)+NPSH(I)+NISH(I),1),
     &           NORBT,NORBT,NZ,
     &           'T','N',IPQTOQ(1,0),
     &           DV(1+IASH(I),1+IASH(I),1),NASHT,NASHT,NZ,
     &           D0,IPQTOQ(1,0),
     &           FQ(1+IORB(I),1+IASH(I),1),NORBT,NASHT,NZ)
         END DO
C
      ELSE IF (CTL .EQ. 'G') THEN
C
C        We are called to generate the whole FQ.
C     
         DO I = 1, NFSYM
            CALL QGEMM(NORB(I),NASH(I),NASH(I),D1,
     &           'N','N',IPQTOQ(1,0),
     &           FV(1+IORB(I),1+IORB(I)+NPSH(I)+NISH(I),1),
     &           NORBT,NORBT,NZ,
     &           'T','N',IPQTOQ(1,0),
     &           DV(1+IASH(I),1+IASH(I),1),NASHT,NASHT,NZ,
     &           D0,IPQTOQ(1,0),
     &           FQ(1+IORB(I),1+IASH(I),1),NORBT,NASHT,NZ)
         END DO
C
      ELSE
C
         WRITE(LUPRI,'(1X,A,A)')
     &        '*** ERROR in APPFQ *** Invalid control parameter: ',CTL
         CALL QUIT('*** ERROR in APPFQ *** Invalid control parameter!')
      END IF
#ifdef UNDEF
C
C     old code that find || FQ - DV * FV ||
C
C        We are called to check the approximation:
C
         CALL MEMGET('REAL',KFQ,NORBT*NASHT*NZ,WORK,KFREE,LFREE)
         CALL DCOPY(NORBT*NASHT*NZ,FQ,1,WORK(KFQ),1)
C
         DO I = 1, NFSYM
            CALL QGEMM(NPSH(I),NASH(I),NASH(I),D1,
     &           'N','N',IPQTOQ(1,0),
     &           FV(1+IORB(I),1+IORB(I)+NPSH(I)+NISH(I),1),
     &           NORBT,NORBT,NZ,
     &           'T','N',IPQTOQ(1,0),
     &           DV(1+IASH(I),1+IASH(I),1),NASHT,NASHT,NZ,
     &           D0,IPQTOQ(1,0),
     &           WORK(KFQ+IORB(I)+IASH(I)*NORBT),
     &           NORBT,NASHT,NZ)
         END DO
         write(6,*) 'approximate pfq'
         call prqmat(work(kfq),norbt,nasht,norbt,nasht,nz,
     &        ipqtoq(1,0),lupri)
         CALL DAXPY(NORBT*NASHT*NZ,DM1,FQ,1,WORK(KFQ),1)
         write(6,*) 'PFQ - FQ'
         call prqmat(work(kfq),norbt,nasht,norbt,nasht,nz,
     &        ipqtoq(1,0),lupri)
         DN = DNORM2(NORBT*NASHT*NZ,WORK(KFQ),1)
         write(6,*) 'norm of diff',DN
C
         CALL MEMREL('APPFQ',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
C
#endif
      IF (IPROPT .GE. 20) THEN
         WRITE(LUPRI,'(1X,A)') '(APPFQ) FQ after add'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,
     &        NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL QEXIT('APPFQ')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dcsct */
      SUBROUTINE DCSCT(DC)
C*****************************************************************************
C
C     Return N2ORBXQ DC.
C
C     Input:
C
C     Output:
C        DV - in-active density matrix with dimenion N2ORBX
C
C     Written by J. Thyssen - Nov 24 2000
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION DC(NORBT,NORBT,NZ)
C
      CALL DZERO(DC,N2ORBXQ)
      DO IFSYM = 1, NFSYM
         DO I = IORB(IFSYM)+NPSH(IFSYM)+1,
     &        IORB(IFSYM)+NPSH(IFSYM)+NISH(IFSYM)
            DC(I,I,1) = D2
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reomat */
      SUBROUTINE MATGAT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C*****************************************************************************
C
C     Construct BMAT from AMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ)
      DIMENSION BMAT(NRB,NCB,NZ)
      DIMENSION IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               BMAT(K,J,IZ) = AMAT(KK,JJ,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reomat */
      SUBROUTINE MATSCT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C*****************************************************************************
C
C     Construct BMAT from AMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ)
      DIMENSION BMAT(NRB,NCB,NZ)
      DIMENSION IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               AMAT(KK,JJ,IZ) = BMAT(K,J,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotfc */
      SUBROUTINE ROTFC(C,U,FC,WORK,LWORK)
C*****************************************************************************
C
C     Transform FC with U.
C     (C controls which block are transformed)
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER C*1
      DIMENSION FC(*)
      DIMENSION WORK(*)
      DIMENSION U(*)
C
      DIMENSION NFC(2), IFC(2), IU(2)
C
#include "memint.h"
C
      CALL QENTER('ROTFC')
C
      CALL MEMGET('REAL',KFC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DCOPY(N2ORBXQ,FC,1,WORK(KFC),1)
C
      CALL MEMGET('REAL',KU,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DUNIT(WORK(KU),NORBT)
C
      IF (C .EQ. 'E') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NESHT,NESHT,IDXE2G,NZ)
      ELSE IF (C .EQ. 'P') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NPSHT,NPSHT,IDXP2G,NZ)
      ELSE IF (C .EQ. 'I') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NISHT,NISHT,IDXI2G,NZ)
      ELSE IF (C .EQ. 'A') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NASHT,NASHT,IDXU2G,NZ)
      ELSE IF (C .EQ. 'V') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NSSHT,NSSHT,IDXV2G,NZ)
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in ROTFC *** Unknown control character: ',C
         CALL QUIT('*** ERROR in ROTFC *** Unknown control character')
      END IF
C
      CALL ROTFC1(C,WORK(KU),FC,WORK(KFC),WORK(KFREE),LFREE)
C
      CALL MEMREL('ROTFC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ROTFC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotfc1 */
      SUBROUTINE ROTFC1(C,U,FC,FCOLD,WORK,LWORK)
C*****************************************************************************
C
C     Transform FC with U.
C     (C controls which block are transformed)
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "dcbidx.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CHARACTER C*1
      DIMENSION FC(*), FCOLD(*)
      DIMENSION U(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
c     write(6,*) 'rotfc1 trans. mat'
c     call prqmat(u,norbt,norbt,norbt,norbt,nz,ipqtoq(1,0),lupri)
C
      CALL QENTER('ROTFC1')
C
      DO I = 1, NFSYM
         CALL QTRANS('AOMO','S',D0,NORB(I),NORB(I),NORB(I),NORB(I),
     &        FCOLD(1+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &        FC(1+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &        U(1+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &        U(1+I2ORBX(I,I)),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,-1)
      END DO
C
      CALL QEXIT('ROTFC1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotfc */
      SUBROUTINE ROTFQ(C,U,FQ,WORK,LWORK)
C*****************************************************************************
C
C     Transform FQ with U.
C     (C controls which block are transformed)
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CHARACTER C*1
      DIMENSION FQ(NORBT,NASHT,*)
      DIMENSION WORK(*)
      DIMENSION U(*)
C
      DIMENSION NFC(2), IFC(2), IU(2)
C
#include "memint.h"
C
      IF (MCTYPE .LT. JMCMIN) RETURN
C     ... return if second-order HF
C
      CALL QENTER('ROTFQ')
C
      CALL MEMGET('REAL',KFQ,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KFQ),N2ORBXQ)
      DO I = 1, NFSYM
         CALL DCOPY(NASH(I)*NORBT,
     &        FQ(1,1+IASH(I),1),1,
     &        WORK(KFQ+(IORB(I)+NPSH(I)+NISH(I))*NORBT),1)
      END DO
C
c     write(6,*) 'fq unpacked'
c     call prqmat(work(kfq),norbt,norbt,norbt,norbt,nz,
c    &     ipqtoq(1,0),lupri)
C
      CALL ROTFC(C,U,WORK(KFQ),WORK(KFREE),LFREE)
C
c     write(6,*) 'fq unpacked after rot'
c     call prqmat(work(kfq),norbt,norbt,norbt,norbt,nz,
c    &     ipqtoq(1,0),lupri)
C
      DO I = 1, NFSYM
         CALL DCOPY(NASH(I)*NORBT,
     &        WORK(KFQ+(IORB(I)+NPSH(I)+NISH(I))*NORBT),1,
     &        FQ(1,1+IASH(I),1),1)
      END DO
C
      CALL MEMREL('ROTFQ',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ROTFQ')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotorb */
      SUBROUTINE ROTORB(C,U,CMO,WORK,LWORK)
C*****************************************************************************
C
C     Transform FC with U.
C     (C controls which block are transformed)
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "dcbbas.h"
#include "dcbidx.h"
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER C*1
      DIMENSION CMO(*)
      DIMENSION WORK(*)
      DIMENSION U(*)
C
      DIMENSION NFC(2), IFC(2), IU(2)
C
#include "memint.h"
C
      CALL QENTER('ROTFC')
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
      CALL MEMGET('REAL',KU,N2ORBXQ,WORK,KFREE,LFREE)
      CALL DUNIT(WORK(KU),NORBT)
C
      IF (C .EQ. 'E') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NESHT,NESHT,IDXE2G,NZ)
      ELSE IF (C .EQ. 'P') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NPSHT,NPSHT,IDXP2G,NZ)
      ELSE IF (C .EQ. 'I') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NISHT,NISHT,IDXI2G,NZ)
      ELSE IF (C .EQ. 'A') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NASHT,NASHT,IDXU2G,NZ)
      ELSE IF (C .EQ. 'V') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NSSHT,NSSHT,IDXV2G,NZ)
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in ROTFC *** Unknown control character: ',C
         CALL QUIT('*** ERROR in ROTFC *** Unknown control character')
      END IF
C
C     Calculate CMO(new) = CMO(old) U
C
      DO I = 1, NFSYM
         NBASI = NFBAS(I,0)
         NORBI = NORB(I)
         CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KCMO+ICMOQ(I)),NBASI,NORBI,NZ,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KU+I2ORBX(I,I)),NORBT,NORBT,NZ,
     &        D0,IPQTOQ(1,0),
     &        CMO(1+ICMOQ(I)),NBASI,NORBI,NZ)
      END DO
C
      CALL MEMREL('ROTFC',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ROTFC')
      RETURN
      END
