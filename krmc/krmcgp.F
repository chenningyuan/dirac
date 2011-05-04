!dirac_copyright_start
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
!dirac_copyright_end

C VERSION : $Revision$
C DATE    : $Date$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,
     &                  set_qphaseFAC)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbham.h"
C
      DIMENSION VCL(NELEM,NZ,NZ,*)
      DIMENSION VNZ3(NELEM,NZ,3)
      DIMENSION IREPIJB(2),IREPKLB(2,*)
      logical   set_qphaseFAC
      real(8)   qphase_FAC
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
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(6,'(2x,a,1i4)') ' (CL2NZ3) IREPIJB(1) = ',IREPIJB(1)
!     write(6,'(2x,a,1i4)') ' (CL2NZ3) IREPIJB(2) = ',IREPIJB(2)
#endif
      IREPABR = 0
      DO IREPAB = 0, NBSYM-1
         IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1

C
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(6,'(/2x,a,2i4)') ' (CL2NZ3) IREPAB and IREPIJB(1) = ',
!    &                         IREPAB, IREPIJB(1)
#endif
C
C           Zero out the spin-orbit coupling for electron 1
            IF (IREPAB.NE.IREPIJB(1)) GOTO 100

            DO I = 1, NELEM
!     write(6,'(/2x,a,3i4)') '(CL2NZ3) element I and irepab vs kl(1,i)',
!    &                        i,irepab,IREPKLB(1,I)
               IF (IREPAB.NE.IREPKLB(1,I)) THEN
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
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(6,'(/2x,a,2i4,D15.6)') ' (CL2NZ3) element I,1,1,IREPABR = ',
!    &                         I,IREPABR,VCL(I,1,1,IREPABR)
!     write(6,'(/2x,a,2i4)') ' (CL2NZ3) multiplication factors are'//
!    &                       ' IREPIJB(2) * IREPKLB(2,I)', 
!    &                         IREPIJB(2), IREPKLB(2,I)
#endif
                   if(set_qphaseFAC)then
                     qphase_FAC         =     IREPIJB(2)*IREPKLB(2,I)
                   else
                     qphase_FAC         = abs(IREPIJB(2)*IREPKLB(2,I))
                   end if
                   VCL(I,1,1,IREPABR) = VCL(I,1,1,IREPABR) * qphase_FAC
               ENDIF
            ENDDO
            CALL DAXPY(NELEM,1.0D0,
     &                 VCL(1,1,1,IREPABR),1,
     &                 VNZ3(1,1,1),1)
 100        CONTINUE
         ENDIF
      ENDDO
C
      ELSE ! not spinfree

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
!        ---------------------------------------------------------
!        zero out integrals that should be zero in linear symmetry
!        but aren't due to numerical noise
!        ---------------------------------------------------------
!        FIXME!!! - SK - August 2010. (linzero is in moltra/traout.F)
!        the routine below does not work yet for classes 2 and 3 but i
!        am working on it...
!
         if(linear)then
!          call linzero_nz3(iclass,irepijb,irepklb(1,1),nelem,vnz3)
         end if

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
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE X2UU(HX,H2XY)
C***********************************************************************
C
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "pgroup.h"
#include "maxorb.h"
#include "dcbidx.h"
C
C
      DIMENSION HX(NORBT,NORBT,NZ,3)
      DIMENSION H2XY(NASHT,NASHT,NZ,3)
C
      if(.not.spinfr_krmc)then
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
      else
        DO IU = 1, NASHT
           IUG = IDXU2G(IU)
           DO IV = 1, NASHT
              IVG = IDXU2G(IV)
              H2XY(IU,IV,1,1) = HX(IUG,IVG,1,1)
           END DO
        END DO
      end if

      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C***********************************************************************
C
C     Make (ug) block integrals from (gu) block.
C
C     Input:
C        GUINT: (gu) integrals
C
C     Output:
C        UGINT: (ug) integrals.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "dgroup.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "pgroup.h"
#include "maxorb.h"
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
         if(.not.spinfr_krmc)then
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
         else 
           DO IU = 1, NASHT
              IUG = IDXU2G(IU)
              DO IV = 1, NASHT
                 IVG = IDXU2G(IV)
                 H2AC(IU,IV,IUU,1,1) = HX(IUG,IVG,1,1)
              END DO
           END DO
         end if
C
C        Print section
C
         IF(IPROPT .GE. 30 .AND. DOH2AC) THEN
#if defined MCSCF_DEBUG_SPINFREE
!          WRITE(LUPRI,'(1X,A)') '(X2GU) input HX(*,*,*,*): '
!          DO I3 = 1, 3
!             DO IZ = 1, NZ
!                WRITE(LUPRI,'(1X,A,I3,A,I3,A)')
!    &                'Matrix no. (',IZ,',',I3,'): '
!                CALL OUTPUT(HX(1,1,IZ,I3),1,norbt,1,norbt,
!    &                norbt,norbt,1,LUPRI)
!             END DO
!          END DO
#endif
           WRITE(LUPRI,'(1X,A,I3,A)') '(X2GU) output H2AC(*,*,',IUU,'):'
           DO I3 = 1, 3
              DO IZ = 1, NZ
                 WRITE(LUPRI,'(1X,A,I3,A,I3,A)')
     &                'Matrix no. (',IZ,',',I3,'): '
                 CALL OUTPUT(H2AC(1,1,IUU,IZ,I3),1,NASHT,1,NASHT,
     &                       NASHT,NASHT,1,LUPRI)
              END DO
           END DO
         END IF
      END IF

!     get rid of noise
      CALL DZERO(H2XY, n2orbx * nz * nz * nbsym / nfsym)
      IF (.NOT. ONLYXY) 
     &CALL DZERO(H2YX, n2orbx * nz * nz * nbsym / nfsym)
C
C     Collect (gu) and (ug) elements.
C
      IF (DOFQ) THEN
!        CALL DZERO(H2XY, NORBT * NASHT * NZ * 3)
!        IF (.NOT. ONLYXY) CALL DZERO(H2YX, NORBT * NASHT * NZ * 3)
C
C
         if(.not.spinfr_krmc)then
           DO I3 = 1, 3
              DO IZ = 1, NZ
                 DO IG = 1, NORBT
                    DO IU = 1, NASHT
                       IUG = IDXU2G(IU)
                       H2XY(IG,IU,IZ,I3) = HX(IG,IUG,IZ,I3)
                       IF (.NOT. ONLYXY)
     &                 H2YX(IG,IU,IZ,I3) = HX(IUG,IG,IZ,I3)
                    END DO
                 END DO
              END DO
           END DO
         else  ! I3 = 1 and NZ = 1 for spin-free
           DO IG = 1, NORBT
             DO IU = 1, NASHT
                IUG = IDXU2G(IU)
                H2XY(IG,IU,1,1) = HX(IG,IUG,1,1)
                IF (.NOT. ONLYXY)    
     &          H2YX(IG,IU,1,1) = HX(IUG,IG,1,1)
#if defined MCSCF_DEBUG_SPINFREE
!               if(HX(IUG,IG,1,1).ne.0.0D0) 
!    &          write(lupri,'(2x,a,1p,2d15.6)') ' (X2GU) H2X(act,sec)'//
!    &          'is nonzero... compare ig,iug and iug,ig',
!    &          HX(IG,IUG,1,1), HX(IUG,IG,1,1)
#endif
!                                     ^  ^ this is NOT a typo...
!                                 (HX matrix is symmetric and
!                                 H2X(act,sec) is not calculated for
!                                 spin-free)
             END DO
           END DO
         end if ! spinfree KRMC
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GU2UGT(GUINT,UGTINT,IREPIJ,NORBI,NASHI)
C***********************************************************************
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
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION IFACNZ3(I3,IZ,IQ1,IQ2)
C***********************************************************************
C
C     Factors for integrals.
C
C     Written by Joern Thyssen - Nov 15 2000
C
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BTCH2F(HT,IZ2,IREPABR,HX,IREPIJ,NKL)
C***********************************************************************
C
C     Unpack packed integrals from HB to HF.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INTB2X(HT,IZ2,IREPABR,IREPAB,HX,IREPIJ,NKL)
C***********************************************************************
C
C     Unpack packed integrals from HT to HX.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION HT(NKL,NZ,NZ,*)
      DIMENSION HX(NORBT,NORBT,NZ,NZ,*)
C
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
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(lupri,'(/2x,a,l)') ' (INTB2X) symmetrize matrix == ',DOFULL
!     write(lupri,'(/2x,a,i4)') ' (INTB2X) IREPABR == ',IREPABR
#endif
C
      DO IZ1 = 1, NZ
         IKL = 0
         DO IREPL = 1, NFSYM
            IREPK = MOD (IREPL + IREPIJ, 2) + 1
            IF (DOFULL) THEN
               IQ1 = IPQTOQ(IZ1,IREPAB)
               IF (IQ1 .EQ. 1) THEN
                  QFAC = 1.0D0
               ELSE
                  QFAC = -1.0D0
               END IF
               DO ICOL = 1, NIDX4(IREPL)
                  JJ = IDXT2G(ICOL+IIDX4(IREPL),4)
                  DO IROW = 1, NIDX3(IREPK)
                     II = IDXT2G(IROW+IIDX3(IREPK),3)
                     IKL = IKL + 1
                     HTVALUE = HT(IKL,IZ1,IZ2,IREPABR)
                     HX(II,JJ,IZ1,IZ2,IREPABR) = HTVALUE
                     HX(JJ,II,IZ1,IZ2,IREPABR) = QFAC*HTVALUE
                     END DO
               END DO
            ELSE
!           stefan: we end up here
               DO ICOL = 1, NIDX4(IREPL)
                  JJ = IDXT2G(ICOL+IIDX4(IREPL),4)
!#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(lupri,'(2x,a,i4,a,i6)') 
!    &   ' (INTB2X) ICOL ==',ICOL, ': index JJ =',JJ
!#endif
                  DO IROW = 1, NIDX3(IREPK)
                     II = IDXT2G(IROW+IIDX3(IREPK),3)
                     IKL = IKL + 1
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(lupri,'(2x,a,i4,a,i6,a,i6)') 
!    &   ' (INTB2X) IROW ==',IROW, ': index II =',II, ' index IKL',IKL
!     write(lupri,'(2x,a,i4)') 
!    &   ' (INTB2X) its absolute address is: (jj-1)*norbt + ii',
!    &     (jj-1)*norbt + ii
#endif
                     HX(II,JJ,IZ1,IZ2,IREPABR) =
     &                    HT(IKL,IZ1,IZ2,IREPABR)
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!     write(lupri,'(2x,a,1D15.6)') 
!    &   ' (INTB2X) new integral is =',HX(II,JJ,IZ1,IZ2,IREPABR)
#endif
                  END DO
               END DO
            END IF
#ifdef HJ_DISABLED_18jun2004
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
#endif
         END DO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GU2F(GUINT,IZ2,IREPAB,IREPABR,HX,IREPIJ,NKL)
C***********************************************************************
C
C     Unpack packed integrals from HB to HF.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE F2BTCH(HX,HT,IZ2,IREPABR,IREPIJ,NKL)
C***********************************************************************
C
C     Pack integrals from HT into HX.
C     I.e. go from dimensioning N2ORBTQ to N2ORBXQ.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE F2GU(HX,HT,IZ2,IREPABR,IREPIJ,NKL)
C***********************************************************************
C
C     Pack full (gg) matrix into symmetry packed (gu) matrix.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLKZERO(J1MO,N1MO,J2MO,N2MO,AMAT,NRDIM,NCDIM,
     &              NZ,NFSYM,IREPIJ)
C***********************************************************************
C
C
C     Zero out block defined by the offsets J{1,2}MO, index lengths
C     N{1,2}MO. Only zero out the quaternion blocks with symmetry IREPIJ.
C
C
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLKUNIT(JMO,NMO,AMAT,NRDIM,NCDIM,NZ,NFSYM)
C***********************************************************************
C
C
C     Zero out block defined by the offsets J{1,2}MO, index lengths
C     N{1,2}MO. Only zero out the quaternion blocks with symmetry IREPIJ.
C
C
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RPRCI(CIVEC,NZCONF,NZ,CIPRG,THRPCI,LUPRI)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
C
      DIMENSION CIVEC(NZCONF,*)
      DIMENSION IP(4)
      CHARACTER CIPRG*(*)
C
chj   IF (CIPRG .EQ. 'GASCIP') THEN
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RAPPBOS(CMO,IBOSYM,IBEIG,IPRINT,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
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
      DIMENSION IBOSYM(*), IBEIG(*)
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
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','RAPPBOS')
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

#if defined MCSCF_DEBUG_SPINFREE
!     additional check: ibosym <--> ibeig
      if(spinfr.or.levyle)then
        
        do  i = 1, norbt
          if(ibosym(i) .ne. ibeig(i))then
            WRITE(LUPRI,'(/A)')
     &        ' (RAPPBOS) ibosym: Approximate boson symmetries: '
         WRITE(LUPRI,9000) (REP(IBOSYM(j)),j=1,NORBT)
            WRITE(LUPRI,'(/A)')
     &        ' (RAPPBOS) ibeig: Approximate boson symmetries: '
         WRITE(LUPRI,9000) (REP(ibeig(j)),j=1,NORBT)
            call quit('*** error in RAPPBOS: 
     &      wrong approximate boson symmetries. ***')
          else
            WRITE(LUPRI,'(A,i6)')
     &' (RAPPBOS) boson symmetry from IBOSYM and IBEIG for orbital # =',
     &  i         
              WRITE(LUPRI,'(10x,2a4/)') REP(IBOSYM(i)), REP(IBEIG(i))
          end if
        end do
      end if
#endif

 9000 FORMAT(10X,15(A4))
C
      CALL MEMREL('RAPPBOS',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RGMCMO(CMO,WORK,LWORK)
C***********************************************************************
C
C     Gram-Schmidt orthogonalize coefficients.
C
C     Input: CMO
C
C     Output: CMO re-orthonomalized using the Gram-Schmidt algorithm
C
C     Written by J. Thyssen - Jan 23 2001
C     Revised Aug 2008, H.J.Aa.Jensen (control output with IPROPT)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"
C
      DIMENSION CMO(*), WORK(*)
C
#include "memint.h"
C
      CALL QENTER('RGMCMO')
C
C     Matrix to hold orbitals.
C
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
      CALL GTOVLX(WORK(KSAO),D1)
      WRITE (LUPRI,'(/A)')
     &   ' -> Gram-Schmidt orthonormalizing MO coefficients.'
C
C     Gram-Schmidt orthogonalize coefficients
C
      JSAO = KSAO
      DO I = 1, NFSYM
         NORBI = NORB(I)
         NFBASI = NFBAS(I,0)
         CALL MEMGET('REAL',KSCR1,NFBASI*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KSCR2,NORBI,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KSAO1,NFBASI*NFBASI,WORK,KFREE,LFREE)
         JSAO1 = KSAO1
         DO J = 1, NFBASI
            CALL DCOPY(NFBASI,WORK(JSAO),1,WORK(JSAO1),1)
            JSAO  = JSAO  + NTBAS(0)
            JSAO1 = JSAO1 + NFBASI
         END DO
         JSAO = JSAO + NFBASI
         CALL MEMCHK('after copy of overlap mat',WORK,1)
         IF (IPROPT .GE. 10) THEN
            WRITE(LUPRI,'(A,I3)')
     &        ' RGMCMO: Coefficients for symmetry ',I
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBASI,NORBI,
     &        NFBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IF ( NPSH(I) .EQ. 0) THEN
            CALL QNORM(WORK(KSAO1),CMO(1+ICMOQ(I)),NFBASI,NORBI,
     &                 NZ,IPQTOQ(1,0),IRETUR,WORK(KSCR1),WORK(KSCR2))
         CALL MEMCHK('after first QNORM',WORK,1)
         ELSE
C           --- change CMO order from positron,electron to
C               electron,positron for the Gram-Schmidt
C               orthonormalization,
C               because we want to start with the occupied orbitals
            CALL MEMGET('REAL',KCMO1,NFBASI*NORBI*NZ,WORK,KFREE,LFREE)
            NPSHI = NPSH(I)
            NESHI = NESH(I)
            LPCMO = NPSHI * NFBASI
            LECMO = NESHI * NFBASI
            JCMO  = 1+ICMOQ(I)
            JCMO1 = KCMO1
            DO IZ = 1,NZ
               CALL DCOPY(LPCMO,CMO(JCMO),1,WORK(JCMO1+LECMO),1)
               CALL DCOPY(LECMO,CMO(JCMO+LPCMO),1,WORK(JCMO1),1)
               JCMO  = JCMO  + NORBI*NFBASI
               JCMO1 = JCMO1 + NORBI*NFBASI
            END DO
            CALL QNORM(WORK(KSAO1),WORK(KCMO1),NFBASI,NORBI,
     &                 NZ,IPQTOQ(1,0),IRETUR,WORK(KSCR1),WORK(KSCR2))
            JCMO  = 1+ICMOQ(I)
            JCMO1 = KCMO1
            DO IZ = 1,NZ
               CALL DCOPY(LPCMO,WORK(JCMO1+LECMO),1,CMO(JCMO),1)
               CALL DCOPY(LECMO,WORK(JCMO1),1,CMO(JCMO+LPCMO),1)
               JCMO  = JCMO  + NORBI*NFBASI
               JCMO1 = JCMO1 + NORBI*NFBASI
            END DO
         END IF
C
         CALL MEMREL('RGMCMO.QNORM',WORK,KWORK,KSCR1,KFREE,LFREE)
C
         IF (IPROPT .GE. 10) THEN
            WRITE(LUPRI,'(A,I3)')
     &        ' RGMCMO: Coefficients for symmetry after GM',I
            CALL PRQMAT(CMO(1+ICMOQ(I)),NFBASI,NORBI,
     &        NFBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END DO
      CALL MEMREL('RGMCMO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RGMCMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE APPFQ(CTL,FQ,DV,FV,WORK,LWORK)
C***********************************************************************
C
C     Approximate positronic part of FQ with (DV*FV).
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision :
C
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DCSCT(DC)
C***********************************************************************
C
C     Return N2ORBXQ DC.
C
C     Input:
C
C     Output:
C        DC - inactive density matrix with dimension N2ORBX
C
C     Written by J. Thyssen - Nov 24 2000
C
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ROTFC(C,U,FC,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
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
      ELSE IF (C .EQ. 'S') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NSSHT,NSSHT,IDXS2G,NZ)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ROTFC1(C,U,FC,FCOLD,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "maxorb.h"
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ROTFQ(C,U,FQ,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ROTORB(C,U,CMO,WORK,LWORK)
C***********************************************************************
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
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "dcbbas.h"
#include "maxorb.h"
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
      CALL QENTER('ROTORB')
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
      ELSE IF (C .EQ. 'S') THEN
         CALL MATSCT(WORK(KU),NORBT,NORBT,U,NSSHT,NSSHT,IDXS2G,NZ)
      ELSE
         WRITE(LUPRI,'(2A)')
     &        '*** ERROR in ROTORB *** Unknown control character: ',C
         CALL QUIT('*** ERROR in ROTORB *** Unknown control character')
      END IF
C
C     Calculate CMO(new) = CMO(old) * U
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
      CALL MEMREL('ROTORB',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ROTORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIAGC(C,AM,EIG,EIGVEC,NDIM,WORK,LWORK)
C***********************************************************************
C
C     diagonalize square matrix AM, return eigenvectors and eigenvalues.
C
C     C controls which type of matrix is diagonalized - dimension NDIM2
C
C     Written by S. Knecht - July 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "maxorb.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
C
      CHARACTER C*1
      DIMENSION EIG(*), EIGVEC(NDIM,NDIM,NZ), AM(NDIM,NDIM,NZ)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('DIAGC')
C
      IOFF2 = 1
C
C     special matrix X
      IF (C .EQ. 'X') GOTO 100
C
      DO IFSYM = 1, NFSYM
C
        IF (C .EQ. 'E') THEN
           NDIM2 = NESH(IFSYM)
        ELSE IF (C .EQ. 'P') THEN
           NDIM2 = NPSH(IFSYM)
        ELSE IF (C .EQ. 'I') THEN
           NDIM2 = NISH(IFSYM)
        ELSE IF (C .EQ. 'A') THEN
           NDIM2 = NASH(IFSYM)
        ELSE IF (C .EQ. 'V') THEN
           NDIM2 = NVSH(IFSYM)
        ELSE
           WRITE(LUPRI,'(2A)')
     &          '*** ERROR in DIAGC *** Unknown control character: ',C
           CALL QUIT('*** ERROR in DIAGC *** Unknown control character')
        END IF
C
C       diagonalize AM
C
C
        IF( NZ .eq. 1 ) THEN
C          ... calculate eigenvalues and eigenvectors
           IJOB   = 1
C          ... ascending order
           IORDER = 1
           IPACK  = 0
           CALL RSJACO(NDIM,NDIM2,NDIM2,AM(IOFF2,IOFF2,1),
     &                 EIG(IOFF2),IJOB,IORDER,IPACK,
     &                 EIGVEC(IOFF2,IOFF2,1),WORK(KFREE),LFREE)
        ELSE
C          ... calculate eigenvalues and eigenvectors
           MATZ = 1
           CALL QDIAG(NZ,NDIM2,AM(IOFF2,IOFF2,1),NDIM,NDIM,
     &                EIG(IOFF2),MATZ,EIGVEC(IOFF2,IOFF2,1),NDIM,NDIM2,
     &                WORK(KFREE),LFREE,IERR)
           IF( IERR .ne. 0 ) THEN
                  WRITE(LUPRI,'(/,1X,4A,I4)')
     &                 '*** ERROR in DIAGC ***: ',
     &                 '*** ERROR for control parameter: ',C,
     &                 'QDIAG failed with error code ',IERR
                  CALL QUIT('*** ERROR in DIAGC: block diag ***')
           END IF
        ENDIF
C       ... new offset
        IOFF2 = IOFF2 + NDIM2
      ENDDO
C
C     special matrix X
 100  IF (C .EQ. 'X') THEN
        NDIM2 = NDIM
C
C       ... calculate eigenvalues and eigenvectors
        IF (NZ .EQ. 1) THEN
           IJOB   = 1
C          ... no sorting
           IORDER = 0
           IPACK  = 0
           CALL RSJACO(NDIM,NDIM2,NDIM2,AM,EIG,IJOB,IORDER,IPACK,
     &                 EIGVEC,WORK(KFREE),LFREE)
        ELSE
           MATZ = 1
           CALL   QDIAG(NZ,NDIM2,AM,NDIM,NDIM,EIG,MATZ,EIGVEC,
     &                  NDIM,NDIM2,WORK(KFREE),LFREE,IERR)
           IF (IERR .NE. 0) THEN
              WRITE(LUPRI,'(/,1X,4A,I4)')
     &        '*** ERROR in DIAGC ***: ',
     &        '*** ERROR for control parameter: ',C,
     &        'QDIAG failed with error code ',IERR
              CALL QUIT('*** ERROR in DIAGC: block diag ***')
           END IF
        END IF
      END IF
C
      CALL QEXIT('DIAGC')
      END
!***********************************************************************
!     SUBROUTINE ONEDENS_DIAG_BLOCK(WORK,LWORK)
!***********************************************************************
!
!     diagonalize matrix in terms of mj-blocks, 
!     return eigenvectors and eigenvalues.
!
!
!     written by S. Knecht - April 2010
!
!***********************************************************************
!#include "implicit.h"
!#include "priunit.h"
!
!     DIMENSION EIG(*), EIGVEC(NDIM,NDIM,NZ), AM(NDIM,NDIM,NZ)
!     DIMENSION WORK(*)
!
!#include "memint.h"
!     IOFF2 = 1
!
!     ... calculate eigenvalues and eigenvectors
!     IJOB   = 1
!     ... no ordering
!     IORDER = 0
!     IPACK  = 0
!     CALL RSJACO(NDIM,NDIM2,NDIM2,AM(IOFF2,IOFF2,1),
!    &            EIG(IOFF2),IJOB,IORDER,IPACK,
!    &            EIGVEC(IOFF2,IOFF2,1),WORK(KFREE),LFREE)
!
!     END
!***********************************************************************
      subroutine block_one_dens(mjvector,amat)
!***********************************************************************
!
!     zero numerical noise in mj-blocked matrix amat
!
!
!     written by S. Knecht - April 2010
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"
      dimension mjvector(norbt), amat(nasht,nasht,nz)
      parameter(xtol_dens = 1.0D-3, d0=0.0D0)
      logical call_quit

      call_quit = .false.
      ioff1 = 0
      ioff2 = 0
      dmj_offmax = 1.0D-15
!
!     loop over active orbitals in each fermion corep
      do ifsym = 1, nfsym
        do mjorb = 1, nash(ifsym)
          mjvec_mjorb = mjvector(ioff2+npsh(ifsym)+nish(ifsym)+mjorb)

!         write(lupri,'(A,2I5)') 'mjvec_mjorb,mjorb:',mjvec_mjorb,mjorb
!         write(lupri,'(A,I5)') 'off-mj #1:',
!    &    ioff2+npsh(ifsym)+nish(ifsym)+mjorb

          do njorb = mjorb+1,nash(ifsym)

!         write(lupri,'(A,2I5)') 'mjvector(njorb),njorb:',
!    &    mjvector(ioff2+npsh(ifsym)+nish(ifsym)+njorb),njorb
!         write(lupri,'(A,I5)') 'off-mj #2:',
!    &    ioff2+npsh(ifsym)+nish(ifsym)+njorb

            if( mjvector(ioff2+npsh(ifsym)+nish(ifsym)+njorb) 
     &          .ne. mjvec_mjorb )then
!
!             test element
              do iz = 1, nz
                dmj_offmax =
     &          max(dmj_offmax,dabs(amat(ioff1+mjorb,ioff1+njorb,iz)))
              end do
!
!             stop if symmetry-breaking element is too large and print
!             diagnosis
!               
              if(dmj_offmax .gt. xtol_dens) then
                write(lupri,'(/A,I2,A)') '*** warning: mj-symmetry '//
     &        'breaking off-diagonal element larger than tolerance'//
     &        ' in fermion irrep',ifsym,'. ***'
                write(lupri,'(/A,1F14.9,A,1F14.9)') 
     &        'tolerance:',xtol_dens,' off-diagonal element:',dmj_offmax
                write(lupri,'(/A,1I4,A,1I4)') 
     &        'adress #1:',mjorb,' adress #2:',njorb
                write(lupri,'(/A,1I4,A,1I4)') 
     &        'mj-value #1:',mjvec_mjorb,' mj-value #2:',
     &         mjvector(ioff2+npsh(ifsym)+nish(ifsym)+njorb)
                dmj_offmax = 1.0D-15
                if( call_quit )then
                  write(lupri,'(/A)') ' *** full mj-vector *** '
		      call iwrtmamn(mjvector,1,norbt,1,norbt,lupri)
                   do ix = 1, nfsym
                     write(lupri,'(/1X,A,I3)')
     &               '*** full density matrix for fermion irrep ***',ix
           call prqmat(amat(1+iash(ix),1+iash(ix),1),nash(ix),nash(ix),
     &                 nasht,nasht,nz,ipqtoq(1,0),lupri)
                   end do
                  call quit(
     & ' *** error in block_one_dens: symmetry-breaking element 
     &   too large. Cannot continue. ***')
                end if
              end if

              do iz = 1, nz
                amat(ioff1+mjorb,ioff1+njorb,iz) = d0
                amat(ioff1+njorb,ioff1+mjorb,iz) = d0
              end do

            end if

          end do
        end do
        ioff1 = ioff1 + nash(ifsym)
        ioff2 = ioff2 + norb(ifsym)
      end do
#if defined MCSCF_DEBUG	
      write(lupri,'(1X,A,1F14.9)') 
     &	' * block_one_dens info: largest mj-symmetry breaking'//
     &      ' off-diagonal element in the DV matrix is',dmj_offmax
#endif

      end subroutine block_one_dens
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TRKRLUCI(CIRUN,ISYNC_NODES)
C***********************************************************************
C
C
C     Transfer information from KRMC/KRCI common blocks to
C     LUCIAREL common blocks.
C     Replaces old LUCI_TRACOM and LUCI_TRAVA.
C
C     ISYNC_NODES controls communication between MASTER and slaves.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision:
C
C***********************************************************************
      use symmetry_setup_krci
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "pgroup.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dcbbos.h"
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/crun.inc"
#include "../luciarel/cprnt.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/oper.inc"
#include "../luciarel/cands.inc"
#include "../luciarel/noccn_inf.inc"
#include "krmc_shmem.h"
#include "krmcluci_inf.h"
#include "krciprop.h"
C
      CHARACTER*6 CIRUN
      logical isopen, ex_mcdf
!     scratch
      integer, allocatable  :: orb_sym_vec(:)
      integer               :: orb_sym_fh
      character (len=11)    :: orb_sym_fn
      character (len= 4)    :: mynum_str
 
!     basic symmetry information
#ifdef VAR_MPI
      if(isync_nodes .eq. 1)then
       CALL MPI_BCAST(GROUP, 3,MPI_CHARACTER,MASTER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(spinfr,1,MPI_LOGICAL,  MASTER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(linear,1,MPI_LOGICAL,  MASTER,MPI_COMM_WORLD,IERR)
      end if
#endif

!     distribute orbital symmetry information (if applicable)
      if(linear.or.spinfr)then

        call alloc(orb_sym_vec,norbt)
        orb_sym_vec = 0

!       read data from master file krmcscf or dfcoef
        if(myproc.eq.master)then
          inquire(file="KRMCSCF",opened=isopen,exist=ex_mcdf)
          if(ex_mcdf)then
            if(.not.isopen) open(lukrmc,status="old",
     &                           form="unformatted")
            rewind(lukrmc)
            if(linear)then
              call ireakrmc(lukrmc,'MJVEC   ',orb_sym_vec,norbt)
            else if(spinfr)then
              call ireakrmc(lukrmc,'IBEIG   ',orb_sym_vec,norbt)
            end if
            if(.not.isopen) close(lukrmc,status="keep")
          else
            call opnfil(lucoef,'DFCOEF','OLD','mjlist ')
            call reacmo_no_work(lucoef,'DFCOEF',dummy,dummy,orb_sym_vec,
     &                          dummy,8)
            close(lucoef,status='KEEP')
          end if
        end if
!       communicate data
#ifdef VAR_MPI
        if(isync_nodes .eq. 1)then
          call mpi_bcast(orb_sym_vec,norbt,my_MPI_INTEGER,master,
     &                   MPI_COMM_WORLD,ierr)
        end if
#endif
!       store orbital symmetry information on file 'orbsym.xxxx' 
!       where xxxx is the process id.
        orb_sym_fh = 99
        call num2str(myproc,mynum_str)
        write(orb_sym_fn,'(a7,a4)') 'orbsym.',mynum_str
        open(orb_sym_fh,file=orb_sym_fn,status='replace',
     &       form='unformatted',action='readwrite',position='rewind')
        write(orb_sym_fh) (orb_sym_vec(i), i=1,norbt)
        close(orb_sym_fh)

        call dealloc(orb_sym_vec)
      end if

!     set double group and point group (for historical reasons and it is
!     still widely used in the code)

!     initialize PNTGRP and DOUGRP
      pntgrp = 0
      dougrp = 0

      IF (GROUP.eq.'C1 ') THEN
         PNTGRP      = 5
         DOUGRP      = 8
         NIRR_DG     = 2
         NIRR_PN     = 1
         IRIQ        = 2
      ELSE IF (GROUP.eq.'Ci ') THEN
         PNTGRP      = 6
         DOUGRP      = 7
         NIRR_DG     = 4
         NIRR_PN     = 2
         IRIQ        = 2
      ELSE IF (GROUP.eq.'Cs ') THEN
         PNTGRP      = 7
         DOUGRP      = 6
         NIRR_DG     = 4
         NIRR_PN     = 2
         IRIQ        = 2
      ELSE IF (GROUP.eq.'C2 ') THEN
         PNTGRP      = 8
         DOUGRP      = 5
         NIRR_DG     = 4
         NIRR_PN     = 2
         IRIQ        = 2
      ELSE IF (GROUP.eq.'D2 ') THEN
!        D2 treated as C2 inside KR-CI
         PNTGRP      = 8
         DOUGRP      = 5
         NIRR_DG     = 4
         NIRR_PN     = 2
         IRIQ        = 1
      ELSE IF (GROUP.eq.'C2h') THEN
         PNTGRP      = 9
         DOUGRP      = 4
         NIRR_DG     = 8
         NIRR_PN     = 4
         IRIQ        = 2
      ELSE IF (GROUP.eq.'C2v') THEN
         if(linear)then
!          Cinfv (known as C2v + linear) treated as C32 inside KR-CI
           PNTGRP      = 10
           DOUGRP      = 10
           NIRR_DG     = 128
           NIRR_PN     = NIRR_DG/2
           IRIQ        = 1
         else
!          C2v treated as C2 inside KR-CI
           PNTGRP      = 8
           DOUGRP      = 5
           NIRR_DG     = 4
           NIRR_PN     = 2
           IRIQ        = 1
         end if
      ELSE IF (GROUP.eq.'D2h') THEN
         if(linear)then
!          Dinfh (known as D2h + linear) treated as C16h inside KR-CI
           PNTGRP      = 11
           DOUGRP      = 11
           NIRR_DG     = 128
           NIRR_PN     = NIRR_DG/2
           IRIQ        = 1
         else if(spinfr.and..not.linear)then
!          D2h spinfree
           PNTGRP      = 9
           DOUGRP      = 9
           NIRR_DG     = 8
           NIRR_PN     = 8
           IRIQ        = 1
           print *,
     & '*** error in trkrluci: no spinfree formalism yet in KR-CI. ***'
           call quit(
     & '*** error in trkrluci: no spinfree formalism in KR-CI! ***')
         else
!          D2h treated as C2h inside KR-CI
           PNTGRP      = 9
           DOUGRP      = 4
           NIRR_DG     = 8
           NIRR_PN     = 4
           IRIQ        = 1
         end if
      ELSE
         call quit(
     & '*** error in trkrluci: unknown double group in KR-CI 
     &  specified! ***')
      END IF

!     consistency check
      if(max(2,NIRR_DG).gt.mxndgirr)then
        write(6,*) '   *** error in trkrluci: number of double'//
     &  ' group irreps exceeds current limit of 128 (mxpdim.inc).***'
        call quit('*** error in trkrluci. max number of double group
     &  irreps exceeded.***')
      end if

!     set symmetry relations for excitation, etc.
      NSMOB       = NIRR_DG
      NSMSX       = NIRR_DG
      NSMDX       = NIRR_DG
      NSMST       = NIRR_DG
      NSMCI       = NIRR_DG
      NSMXT       = NIRR_DG
      itssx       = 1
      itsdx       = 1
      itsxt       = 1

C     Reference symmetry
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) 
     &   CALL MPI_BCAST(IOPT_SYMMETRY,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      IREFSM = IOPT_SYMMETRY
C     max reference symmetry
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) 
     &   CALL MPI_BCAST(NMAX_SYM,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      NMAX_SYM_CI = NMAX_SYM
C
C     KRCI_CVECS.x vector file extensions
      CALL DETXFLAB(XSYMFLAB,NMAX_SYM_CI)
C     Roots to be obtained
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) 
     &   CALL MPI_BCAST(NCIROOT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      NROOT = NCIROOT
      DO IRT = 1,NROOT
         IROOT(IRT) = IRT
      END DO
C     Restart option - special treatment for some CIRUN flags
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) 
     &   CALL MPI_BCAST(IRESTRK,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      IRESTR = IRESTRK
      IF (CIRUN .EQ. 'SIGMA ' .OR. CIRUN .EQ. 'SIGMAD'
     &    .OR. (CIRUN .EQ. 'CIINII' .AND. IRESTR .EQ. 1)
     &    .OR. (CIRUN .EQ. 'KR-CI ' .AND. IRESTR .EQ. 1) )THEN
         IRESTR = 1
      ELSE
         IRESTR = 0
      END IF
C     Integral environment
      INTIMP = 7
      ENVIRO(1:6) = 'DIRAC '
C     Integrals all held in memory
      INCORE = 1
C     print flags
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) 
     &   CALL MPI_BCAST(IPROPT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      IF (IPROPT.EQ.1) THEN
         IPRSTR    = 0 * IPROPT
         IPRCIX    = 0 * IPROPT
         IPRORB    = 0 * IPROPT
         IPRDIA    = 0 * IPROPT
         IPRXT     = 0 * IPROPT
         IPROCC    = 0 * IPROPT
         IPRDEN    = 0 * IPROPT
         IPRSIG    = 0 * IPROPT
         IPRTRA    = 0 * IPROPT
         IPRHAM_CI = 0 * IPROPT
      ELSE
         IPRSTR    = 1   * IPROPT
         IPRCIX    = 0   * IPROPT
         IPRORB    = 1   * IPROPT
         IPRDIA    = 5   * IPROPT
         IPRXT     = 0   * IPROPT
         IPROCC    = 2   * IPROPT
         IPRDEN    = 2   * IPROPT
         IPRSIG    = 10  * IPROPT
         IPRTRA    = 1   * IPROPT
         IPRHAM_CI = 100 * IPROPT
      END IF
C     Maximum number of CI vectors
#if defined (VAR_MPI)
       IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(MXCIV,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      MXCIV_CI = MXCIV
      IF( MXCIV_CI .eq. 0 ) THEN
          MXCIV_CI = 3 * NCIROOT
      END IF
C     CI storage mode
      ICISTR = 3
C     No integrals (we use this to calculate just #determinants)
      IF( CIRUN .eq. 'NDET  ') THEN
         NOINT = 1
      ELSE
         NOINT = 0
      END IF
C     Dimension of resolution matrices (used in SIGDEN_REL3)
!     MXINKA = 50 
      MXINKA = 64 ! testing a slightly larger value...
C     Core energy from DIRAC
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(EMY,1,MPI_REAL8,MASTER,MPI_COMM_WORLD,IERR)
#endif
      ECORE = EMY
C     GAS shells and constraints.
      IDOGAS = 1
      NCISPC = 1
      NCMBSPC = NCISPC
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(NGAS_DC,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
      NGAS = NGAS_DC
      IF( MXPNGAS .ne. MXMCGAS) THEN
        WRITE(LUPRI,*) '*** ERROR on common blocks *** ' //
     &  'MXPNGAS (LUCIAREL) .ne. MXMCGAS (DIRAC):',
     &                 MXPNGAS,MXMCGAS
        CALL QUIT('*** ERROR on common blocks *** ' //
     &  'MXPNGAS (LUCIAREL) .ne. MXMCGAS (DIRAC)')
      END IF
#if defined (VAR_MPI)
C     ... just in case these variable have not 
C     been initialized on the outside.
      IF( ISYNC_NODES .eq. 1 ) THEN
        CALL MPI_BCAST(OPT_UCIBOS,1,MPI_LOGICAL,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NBSYM,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NGSH,2*MXPNGAS,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NGSHT,MXPNGAS,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(IORB,NFSYM,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NPSH,NFSYM,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NISH,NFSYM,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(NASH,NFSYM,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(IBOSYM,MXCORB,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
      END IF
#endif
            
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(NGSP,2*MXPNGAS,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      DO I = 1, NGAS
         CALL ADD_TO_GSS(IGSOCCX,I,NGSP(1,I),2)
         CALL ADD_TO_GSS(IGSOCCX,I+MXPNGAS,NGSP(2,I),2)
CSK      IGSOCCX(I          ) = NGSP(1,I)
CSK      IGSOCCX(I + MXPNGAS) = NGSP(2,I)
      END DO
C     Energy convergence threshold
      IF (NZ.EQ.1) THEN
        THRES_E = 1.0D-10
      ELSE
        THRES_E = 0.5D-5
      END IF
C     Sequence of CI calculations (here one!)
      NSEQCI(1) = 1
C     Maximum number of CI iterations
      ISEQCI(1,1) = 1
      IF( CIRUN .eq. 'CIINII' .or. CIRUN .eq. 'KR-CI ' )THEN
#if defined (VAR_MPI)
        IF( ISYNC_NODES .eq. 1 )
     &    CALL MPI_BCAST(MAXCIT,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
#endif
        ISEQCI(1,1) = MAXCIT
      END IF
C     No combinations of degenerate initial vectors
      INIDEG = 0
C     Density matrices requested
      IDENSI = 0
      IF( CIRUN .eq. 'DENS1 ' ) THEN
        IDENSI = 3
      ELSE IF ( CIRUN .eq. 'DENS2 ' ) THEN 
        IDENSI = 4
      END IF
C     Type of operator wrt. relativity
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) THEN
        CALL MPI_BCAST(SPINFR,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(LEVYLE,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      END IF
#endif
C     Relativistic CI
      ISPINFREE = 0
      IF( SPINFR .or. LEVYLE ) THEN
        ISPINFREE = 1
      END IF
C     Reference value for two times M_K projection
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(MK2REF,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
      MK2REF_CI = MK2REF
C     Reference value for two times Delta M_K projection (coupling range)
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(MK2DEL,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
      MK2DEL_CI = MK2DEL
      IF ((.not.SPINFR).and.(.not.LEVYLE).and.MK2DEL_CI.eq.0) THEN
           WRITE(LUWRT,*) 'Spin-dependent calculation and MK2DEL = 0.'
           WRITE(LUWRT,*) 'This is not allowed.'
           CALL QUIT( ' *** ERROR *** Spin-dependent calculation'//
     &                ' and MK2DEL = 0.' ) 
      END IF
C     COMHAM option
      IF (ISEQCI(1,1) .LT. 0 ) THEN
C        ISEQCI(1,1) = MAXCIT from KRMC/KRCI
C        MAXCIT .lt. 0 is a code for doing CI
C        with explicit CI Hamiltonian (not iterative);
C        however, ISEQCI must be .gt. 0 for LUCIAREL to work
         IDCOMH = 1
         ISEQCI(1,1) = 1
      ELSE
C        Use iterative Davidson-Olsen type algorithm
C        with MAXCIT iterations
         IDCOMH = 0
      END IF
C     Always two-particle (full) Hamiltonian 
      IHAM12 = 2
C     calculate property operator symmetry in double group subgroup
C     notation
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(DOSYMOPPRP,1,MPI_LOGICAL,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
C     reorder property operator matrix (already done by the MASTER...)
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(DOPROPREOD,1,MPI_LOGICAL,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
C     Wave function analysis
      IPRWFC = 0
      IF (CIRUN .EQ. 'ANALYZ') IPRWFC = 1
C     Number of active electrons and 2-el. routine calls (sigma)
      NACTEL = IGSOCCX(NGAS,2,1)
      IELCI = 2       ! default
C     Construction of diagonal is controlled by IRESTR flag (above)
      IDIAG = 1
C     No explicit Hamiltonian
      MXP1 = 0
      MXP2 = 0
      MXQ = 0
C     S and L are no good quantum numbers. Don't use combinations.
      PSSIGN = 0.D0
      PLSIGN = 0.D0
      IDC = 1
C     C(Ka,Jb,j) matrices will be used in alpha-beta loop
      ICJKAIB = 1
C     Only 1 GAS space combination in use.
      LCMBSPC(1) = 1
      ICMBSPC(1,1) = 1
C     No use of F + Lambda(H-F) instead of H as operator
      XLAMBDA = 1.0D0
C     No iterative natural spinors
      NATITER = 1
C     No class selection performed.
      ICLSSEL = 0
C     Never skip evaluation of energy (only used in pico3, anyway)
      ISKIPEI = 0
C     No property calculations here (yet)
      NPROP = 0
C     SIGDEN route is default (excitation class formalism)
      SIGDEN_ROUTE = 'SIGDEN'
      if(SIGDEN_ROUTE .eq.'TRADIT') 
     & call quit(' *** error in trkrluci: TRADIT route
     & no longer supported.***')
C     integral handling in parallel runs
      IIOMOD_REL = 0
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(IDOPARIO,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
      IIOMOD_REL = IDOPARIO
C     multiplier for max memory allocation guess
      ISMEMFAC = 0
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &  CALL MPI_BCAST(IMEMFAC,1,my_MPI_INTEGER,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
      ISMEMFAC = IMEMFAC
C     default maximum C/sigma block size
      IMAXLBLKSZ = 100 000 000
C     running a nat. orb. occ. num. run?
      NOOCCN_DENS_RUN = .FALSE.
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(NATOLCR,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      NOOCCN_DENS_RUN = NATOLCR
C     state number for density calculation
      IDENSI_STATE = - 1
      IF( CIRUN .eq. 'DENS1 ' .or. CIRUN .eq. 'DENS2 ') THEN
#if defined (VAR_MPI)
        IF( ISYNC_NODES .eq. 1 )
     &    CALL MPI_BCAST(IDENSLR_STATE,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
#endif
        IDENSI_STATE = IDENSLR_STATE
      END IF
C     running in "shared memory" mode? 
      SHARED_M    =   .FALSE.
      LEVEL_SM    = - 1
      IT_SHL      = - 1
      IC_SHL      = - 1
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CSHMEMO,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      SHARED_M = CSHMEMO
      IF( SHARED_M )THEN
#if defined (VAR_MPI)
        IF( ISYNC_NODES .eq. 1 )
     &    CALL MPI_BCAST(ISHMEM_TYPE,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
#endif
        LEVEL_SM = ISHMEM_TYPE
      END IF
      IF( LEVEL_SM .eq. 0 ) THEN
C       disable shared memory usage
        IT_SHL = - 1
        IC_SHL = - 1
        SHARED_M = .FALSE. 
        LEVEL_SM = - 1
      ELSE IF( LEVEL_SM .eq. 1 ) THEN
        IT_SHL =   0
        IC_SHL = - 1
      ELSE IF( LEVEL_SM .eq. 2 ) THEN
        IT_SHL = - 1
        IC_SHL =   0
      ELSE IF( LEVEL_SM .eq. 3 ) THEN
        IT_SHL =   0
        IC_SHL =   0
      ELSE IF( LEVEL_SM .eq. 4 ) THEN
        IT_SHL =   1
        IC_SHL =   0
      ELSE IF( LEVEL_SM .eq. 5 ) THEN
        IT_SHL =   0
        IC_SHL =   1
      ELSE IF( LEVEL_SM .eq. 6 ) THEN
        IT_SHL =   1
        IC_SHL =   1
      END IF
      IF( CIRUN .eq. 'IJKLRO' )THEN
        IF( IT_SHL .ge. 0 )THEN
          IT_SHL = - 2
        END IF
      END IF
C     write check point file during large-scale CI optimization 
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CHCKPT_WRT,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      CHECKPOINT_LUCIX = CHCKPT_WRT
C     read from / write to block distribution file KRCI_BLOCKDIST.x
C     (x=sym)
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CPBLCK_FILE,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      CPBLCK_FILE_LUCIX = CPBLCK_FILE
C     using GAS-scheme specific integral file resp. writing it
      REORD_IJKL = .FALSE.
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CINT_REORD,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      REORD_IJKL = CINT_REORD
C     using revised integral allocation scheme                
      SPLIT_IJKL = .FALSE.
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CINT_SPLIT,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      SPLIT_IJKL = CINT_SPLIT
C     using low-memory scheme for integral resorting
      LOWSRT_IJKL = .FALSE.
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(CINT_LOWSRT,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#endif
      LOWSRT_IJKL = CINT_LOWSRT
C
C     if (CIRUN .eq. PROP1) disable all 'shared memory' stuff
C                           --> implement C-sharing (IC_SHL) !?!
      IF( CIRUN .eq. 'PROP1 ')THEN
C        disable shared memory usage for now
         IT_SHL      = - 1
         IC_SHL      = - 1
         SHARED_M    = .FALSE.
         LEVEL_SM    = - 1
         REORD_IJKL  = .FALSE.
         SPLIT_IJKL  = .FALSE.
         LOWSRT_IJKL = .FALSE.
      END IF
C     sigma calculation in property run 
      DOSIGPROP = .FALSE.
C
C     do not write file CIDIA_REL.x
#if defined (VAR_MPI)
      NO_CDIAF = .TRUE.
      IF( MYPROC .eq. MASTER )THEN
         IF( MOPT_MXMACRO .ge. 0 ) NO_CDIAF = .FALSE.
      END IF
      IF( ISYNC_NODES .eq. 1 )
     &   CALL MPI_BCAST(NO_CDIAF,1,MPI_LOGICAL,MASTER,
     &                  MPI_COMM_WORLD,IERR)
#else
      NO_CDIAF = .FALSE.
#endif
C
C     flag for pure KR-CI run - enables reuse of already existing
C     DIAPAR.1
      CKRCIONLY = .FALSE.
      IF( CIRUN .eq. 'KR-CI ') CKRCIONLY = .TRUE.
#if defined (VAR_MPI)
      IF( ISYNC_NODES .eq. 1 ) THEN
        IF( PNTGRP .ge. 2 .and. PNTGRP .le. 4 ) THEN
          CALL MPI_BCAST(IREFML,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF( PNTGRP .eq. 3 .or. PNTGRP .eq. 4 ) THEN
          CALL MPI_BCAST(IREFPA,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF( PNTGRP .eq. 4 ) THEN
          CALL MPI_BCAST(IREFL,1,my_MPI_INTEGER,MASTER,
     &                   MPI_COMM_WORLD,IERR)
        END IF
        IF(IDIAG.EQ.2) THEN
          CALL MPI_BCAST(E_THRE,1,MPI_REAL8,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(C_THRE,1,MPI_REAL8,MASTER,
     &                   MPI_COMM_WORLD,IERR)
          IF(ICLSSEL.EQ.1) THEN
            CALL MPI_BCAST(E_CONV,1,MPI_REAL8,MASTER,
     &                     MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(C_CONV,1,MPI_REAL8,MASTER,
     &                     MPI_COMM_WORLD,IERR)
          END IF
        END IF
      END IF 
#endif

      END 
!**********************************************************************

      subroutine print_ci_info()
!**********************************************************************
      use symmetry_setup_krci
#include "implicit.h"
#include "priunit.h"
#include "parluci.h"
#include "pgroup.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbbos.h"
#include "../luciarel/mxpdim.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/crun.inc"
#include "../luciarel/cprnt.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/oper.inc"
#include "../luciarel/cands.inc"
#include "../luciarel/noccn_inf.inc"
#ifdef VAR_MPI
#include "mpif.h"
#include "krmc_shmem.h"
#endif
#include "krmcluci_inf.h"
#include "krciprop.h"
!**********************************************************************
*
*. Type of calculation
*
      WRITE(LUPRI,*) '************************'
      WRITE(LUPRI,*) '*  Type of calculation *'
      WRITE(LUPRI,*) '************************'

      WRITE(LUPRI,'(A)')
     &'     Relativistic calculation in double group '
*
*. Type of reference state
      WRITE(LUPRI,*) '*************************************'
      WRITE(LUPRI,*) '*  Symmetry and spin of CI vectors  *'
      WRITE(LUPRI,*) '*************************************'
*. Point group
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ D2H'
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ C inf v'
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ D inf h'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUPRI,'(A)')
     &  '     Point group ............ O3'
      ELSE IF(PNTGRP.EQ.5) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C1'
      ELSE IF(PNTGRP.EQ.6) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Ci'
      ELSE IF(PNTGRP.EQ.7) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Cs'
      ELSE IF(PNTGRP.EQ.8) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C2'
      ELSE IF(PNTGRP.EQ.9) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... C2H'
      ELSE IF(PNTGRP.EQ.10) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Cinf'
      ELSE IF(PNTGRP.EQ.11) THEN
        WRITE(LUPRI,'(A)')
     &  '     Double group ........... Cinfh'
      END IF
*.Spatial symmetry
      IF(PNTGRP.EQ.1) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     Spatial symmetry ....... ', IREFSM
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Ungerade'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(LUPRI,'(A,I1)')
     &  '     L  value ............... ', IREFL
        WRITE(LUPRI,'(A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(LUPRI,'(A)')
     &  '     Parity   ..............  Ungerade'
       ELSE IF(PNTGRP.GE.5)THEN
        print *,'@timo+stefan: check this print (krmcgp.F) here...'
        WRITE(LUPRI,'(A,I1)')
     &  '     Double group symmetry... ', IREFSM
      END IF
*.Spin
        WRITE(LUPRI,'(A,I2)')
     &    '     2*reference value of MK ', MK2REF_CI
        WRITE(LUPRI,'(A,I2,A,I2)')
     &    '     Allowed interval of MK2 ',
     &          MK2REF_CI-MK2DEL_CI,' to ', MK2REF_CI+MK2DEL_CI
*.Number of active electrons
      WRITE(LUPRI,'(A,I2)')
     &  '     Active electrons .....  ', NACTEL
*
*. GAS space
*
      WRITE(LUPRI,*) ' *************************'
      WRITE(LUPRI,*) ' Generalized active space '
      WRITE(LUPRI,*) ' *************************'
        WRITE(LUPRI,'(A/A)')
     &  ' Orbital subspaces in molecular pointgroup',
     &  ' ========================================='
      WRITE(LUPRI,'(A,10I4,A)')
     &  '                Irrep ',(I,I = 1,nirr_dg/2)
      WRITE(LUPRI,'(A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,nirr_dg/2 )
      DO IGAS = 1, NGAS
        WRITE(LUPRI,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, nirr_dg/2 )
      END DO
*
      WRITE(LUPRI,'(A,I3)')
     & ' Number of Occupation spaces : ',NCISPC
      WRITE(LUPRI,*)
      DO ICISPC = 1, NCISPC
      WRITE(LUPRI,'(A,I3)')
     &' Bounds on accumulated occupations for subspace : ',ICISPC
      WRITE(LUPRI,'(A)')
     & ' ====================================================== '
      WRITE(LUPRI,'(A)')
      WRITE(LUPRI,'(A)') '         Min. occ    Max. occ '
      WRITE(LUPRI,'(A)') '         ========    ======== '
      DO IGAS = 1, NGAS
        WRITE(LUPRI,'(A,I2,3X,I3,9X,I3)')
     &  '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
      END DO
      END DO
*
      WRITE(6,'(/A,I3)')
     &' Number of CI spaces included : ', NCMBSPC
      WRITE(6,*)
      DO JCMBSPC = 1, NCMBSPC
        WRITE(6,*)
        WRITE(6,'(A,I3)') ' Information about CI space ', JCMBSPC
        WRITE(6,'(A)')    ' =================================='
        WRITE(6,'(A,I3)')
     &  ' Number of GAS spaces included ', LCMBSPC(JCMBSPC)
        WRITE(6,'(A,10I3)')
     &  ' GAS spaces included           ',
     &    (ICMBSPC(II,JCMBSPC),II=1,LCMBSPC(JCMBSPC))
*
        WRITE(6,'(A,I3)') ' Number of calculations in this CI space ',
     &  NSEQCI(JCMBSPC)
        WRITE(6,'(A)')   '  Calculations in this subspace '
        WRITE(6,'(A)')   '  =============================='
        DO JSEQ = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQ,JCMBSPC).LT.0.AND.ISEQCI(JSEQ,JCMBSPC).GT.-5)
     &    THEN
            WRITE(6,'(A,I3)')
     &      '       Vector free calculation at level ',
     &      -ISEQCI(JSEQ,JCMBSPC)
          ELSE IF( ISEQCI(JSEQ,JCMBSPC).GT.0 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Normal CI with max. iterations = ',
     &      ISEQCI(JSEQ,JCMBSPC)
          ELSE IF (ISEQCI(JSEQ,JCMBSPC).EQ.-5 ) THEN
           WRITE(6,'(A,I3)')
     &      '       Perturbation calculation '
          END IF
        END DO
      END DO
*
      WRITE(LUPRI,*) '***********'
      WRITE(LUPRI,*) '*  Roots  *'
      WRITE(LUPRI,*) '*********** '
      WRITE(LUPRI,'(A,I3)')
     &  '     Number of roots to be obtained ', NROOT
      WRITE(LUPRI,'(A,(20I3))')
     &  '     Roots to be obtained ', (IROOT(I),I=1, NROOT )
*
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '**************************'
      WRITE(LUPRI,*) '*  Run time definitions  *'
      WRITE(LUPRI,*) '************************** '
      WRITE(LUPRI,*)
*. Program environment
      WRITE(6,'(A,A6)') '      Program environment... ', ENVIRO
*. No. of particle Hamiltonian
      if (IHAM12.eq.1) then
        WRITE(LUPRI,'(A)')
     &  '     Only 1-particle Hamiltonian under consideration '
      else if (IHAM12.eq.2) then
        WRITE(LUPRI,'(A)')
     &  '     Full 1- and 2-particle Hamiltonian active       '
      else if (IHAM12.eq.3) then
        WRITE(LUPRI,'(A)')
     &  '     Only 2-particle Hamiltonian under consideration '
      end if
*. Integral import
      IF(NOINT.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     No integrals will be read in       '
      ELSE IF(NOINT.EQ.0) THEN
      IF(INTIMP.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in MOLCAS format '
      ELSE IF(INTIMP.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in LUCAS format '
      ELSE IF(INTIMP.EQ.3) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in formatted form (E22.15) ',
     &  '      From unit 13'
        WRITE(LUPRI,'(A)')
     &  '     All integrals of correct symmetry combination read in'
      ELSE IF(INTIMP.EQ.5) THEN
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in DALTON/SIRIUS format '
      else if (INTIMP.eq.6) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals in Kramers reduced spinor list format '
        WRITE(LUPRI,'(A)')
     &  '       Dirac-Coulomb Hamiltonian '
      else if (INTIMP.eq.7) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals in (NZ,3) format '
        WRITE(LUPRI,'(A)')
     &  '       Dirac-Coulomb Hamiltonian '
      else if (INTIMP.eq.8) then
        WRITE(LUPRI,'(A)')
     &  '     Integrals read in in MOLECULE-SWEDEN format '
        WRITE(LUPRI,'(A)')
     &  '     one-electron ints from unit 2, transformation '
        WRITE(LUPRI,'(A)')
     &  '     spin-orbit ints from SOMFINT_*, trafo '
        WRITE(LUPRI,'(A)')
     &  '     two-electron ints from unit 13, transformed '
      END IF
      WRITE(LUPRI,'(/A/)')
     &  '     CI optimization performed with SD''s '
*. Initial approximation to vectors
      IF(IRESTR.EQ.1)THEN
         WRITE(LUPRI,'(A)')
     &  '     Restarted calculation '
      ELSE
         IF(MXP1.NE.0) THEN
           WRITE(LUPRI,'(A)')
     &  '     Initial vectors obtained from explicit Hamiltonian'
         ELSE IF(MXP1.EQ.0) THEN
           WRITE(LUPRI,'(A)')
     &  '     Initial vectors obtained from diagonal'
         END IF
      END IF
*. Handling of degeneracies of initial vectors
      IF(INIDEG.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Symmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.-1) THEN
        WRITE(LUPRI,'(A)')
     &  '     Antisymmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     No combination of degenerate initial vectors'
      END IF
*. CI storage mode
        WRITE(LUPRI,'(/A)')
     &  '     SIGDEN route for sigma vectors and densities invoked   '
      IF(ICISTR.EQ.1) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry blocks and two vectors will be held in core '
      ELSE IF( ICISTR.EQ.2) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry blocks will be held in core '
      ELSE IF( ICISTR.EQ.3) THEN
        WRITE(LUPRI,*)
     &  '     3 symmetry-occ-occ blocks will be held in core '
      END IF
      IF(LCSBLK.NE.0) WRITE(LUPRI,'(A,I10)')
     &  '      Smallest allowed size of sigma- and C-batch ',LCSBLK
      WRITE(LUPRI,'(A,I4)')
     &  '     Dimension of block of resolution strings ', MXINKA
      IF(ICJKAIB.EQ.1) THEN
        WRITE(LUPRI,*)
     &  '     C(Ka,Jb,j) matrices will be used in alpha-beta loop '
      END IF
*
      IF(IDENSI.EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '     No calculation of density matrices  '
      ELSE IF(IDENSI.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &  '     One-body density matrix calculated '
      ELSE IF(IDENSI.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &  '     One- and two-body density matrices  calculated '
      else if (IDENSI.eq.3) then
        WRITE(LUPRI,'(A)')
     &  '     Spin-dependent one-body density matrix calculated '
        WRITE(LUPRI,'(A)')
     &  '     Spin-dependent one- and two-body densities
     &        calculated '
      END IF
*
*
*. Diagonalization information
      WRITE(LUPRI,'(A/A)')
     &  '     CI diagonalization : ',
     &  '     ==================== '
*. Subspace Hamiltonian
      IF (MXP1+MXP2+MXQ .EQ.0) THEN
        WRITE(LUPRI,'(A)')
     &  '        No subspace Hamiltonian '
      ELSE
        WRITE(LUPRI,'(A,3I4)')
     &  '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
      END IF
*. Diagonalizer
      IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
        WRITE(LUPRI,'(A)')
     &    '        Diagonalizer : MINDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.EQ.2) THEN
        WRITE(LUPRI,'(A)')
     &    '        Diagonalizer : MICDV4 '
      ELSE IF(IDIAG.EQ.2) THEN
      WRITE(LUPRI,'(A)')
     &  '        Diagonalizer : PICO2  '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.3) THEN
        IF (NZ.EQ.1) THEN
          WRITE(LUPRI,'(A)')
     &      '        Diagonalizer : MICDV6 '
        ELSE
          WRITE(LUPRI,'(A)')
     &      '        Diagonalizer : CMICDV '
        END IF
      END IF
*
      IF(ISKIPEI.EQ.1) THEN
        WRITE(LUPRI,'(A/A)')
     &  '        Initial energy evaluations skipped after first calc',
     &  '        (Only active in connection with TERACI )'
      END IF
*
*. Number of iterations
      WRITE(LUPRI,'(A,I2)')
     &  '        Allowed number of iterations    ',MAXIT
*. Number of CI vectors in subspace
      WRITE(LUPRI,'(A,I2)')
     &  '        Allowed Dimension of CI subspace ',MXCIV_CI
*
      WRITE(LUPRI,'(A,1P,E11.5/)')
     &  '        Convergence threshold for energy ',THRES_E
      IF(IDIAG.EQ.2) THEN
        WRITE(LUPRI,'(A,1P,E11.5)')
     &   '        Individual second order energy threshold ',E_THRE
        WRITE(LUPRI,'(A,1P,E11.5)')
     &   '        Individual first order wavefunction threshold ',C_THRE
        IF(ICLSSEL.EQ.1) THEN
         WRITE(LUPRI,'(/A/A)')
     &   '         Class selection will be performed : ',
     &   '         =================================== '
         WRITE(LUPRI,'(A,1P,E11.5)')
     &    '          Total second order energy threshold ',E_CONV
         WRITE(LUPRI,'(A,1P,E11.5)')
     &    '          Total first order wavefunction threshold ',C_CONV
        ELSE
         WRITE(LUPRI,'(A)')
     &'            No class selection in iterative procedure '
        END IF
      END IF
*
      END IF
C     ^ end if IPROPT for input printing.
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADD_TO_GSS(IMAT,IADR,IVALUE,IJOB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IMAT(*)
      IF( IJOB .eq. 1 )THEN
        IMAT(IADR) = IMAT(IADR) + 1
      ELSE 
        IMAT(IADR) = IVALUE
      END IF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Function SUM_DIAG_ELM */
      REAL*8 FUNCTION SUM_DIAG_ELM(C,AMAT,NRA,NCA)
C***********************************************************************
C
C     Sum diagonal elements of matrix AMAT elements. 
C
C     Character 'C' controls which part of the matrix is summed up.
C
C     Written by S. Knecht - Aug 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "maxorb.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
      PARAMETER (  D0 = 0.0D0 )
      CHARACTER C*1
      DIMENSION AMAT(NRA,NCA,*)
C
      SUMDIAG = D0
      DO IFSYM = 1, NFSYM
C
        IF (C .EQ. 'E') THEN
           NDIM2 = NESH(IFSYM)
           IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
        ELSE IF (C .EQ. 'P') THEN
           NDIM2 = NPSH(IFSYM)
           IOFF1 = IORB(IFSYM) + 1
        ELSE IF (C .EQ. 'I') THEN
           NDIM2 = NISH(IFSYM)
           IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + 1
        ELSE IF (C .EQ. 'A') THEN
           NDIM2 = NASH(IFSYM)
           IOFF1 = IORB(IFSYM) + NPSH(IFSYM) + NISH(IFSYM) + 1
        ELSE IF (C .EQ. 'V') THEN
           NDIM2 = NVSH(IFSYM)
           IOFF1 = IORB(IFSYM)+NPSH(IFSYM)+NISH(IFSYM)+NASH(IFSYM)+1
        ELSE
           WRITE(LUPRI,'(2A)')
     &    '*** ERROR in SUM_DIAG_ELM *** Unknown control character: ',C
          CALL Abend2('* ERROR in SUM_DIAG_ELM * Unknown control 
     &                   character')
        END IF
C
        INUM_S = 0
        DO I = 1, NDIM2
           SUMDIAG = SUMDIAG + AMAT(IOFF1+INUM_S,IOFF1+INUM_S,1)
           INUM_S = INUM_S + 1
        END DO
      END DO
      SUM_DIAG_ELM = SUMDIAG
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_XOMEGALIST(OMLISTREOD,XEIGVC,NXSTATES,
     &                          STATEREP,XEIGVL)
C***********************************************************************
C
C     determine largest element of a given eigenvector to assign the 
C     corresponding eigenstate.
C
C     Written by S. Knecht - Dec 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER STATEREP*4
      DIMENSION OMLISTREOD(*), XEIGVL(*), XEIGVC(NXSTATES,NXSTATES,2)
C
      XLARGE  =  0.0D0
      XMDUM   = -1.0D98
      XXMDUM  = -1.0D-98
      DO IXSTATE = 1, NXSTATES
        OMLISTREOD(IXSTATE) = XMDUM
      END DO
      DO IXSTATE = 1, NXSTATES
        XMAXVL = 0.0D0
C
C       find largest element for NZ == 1 and NZ == 2 of a given state
        IMAXNONE = IDAMAX(NXSTATES,XEIGVC(1,IXSTATE,1),1)
        IMAXNTWO = IDAMAX(NXSTATES,XEIGVC(1,IXSTATE,2),1)
        IMAXIND  = IMAXNONE
        IF( ABS(XEIGVC(IMAXNTWO,IXSTATE,2)) .gt.
     &      ABS(XEIGVC(IMAXNONE,IXSTATE,1))) IMAXIND = IMAXNTWO
        IMAXNTWO = 0
C
C       check if this index is already assigned (which may happen in
C       some cases - H_2 molecule with 1s1p function) - if yes - goto 
C       second largest element.
        IF( OMLISTREOD(IMAXIND) .eq. XMDUM )THEN
          OMLISTREOD(IMAXIND) = XEIGVL(IXSTATE)
        ELSE
          XLARGE = ABS(XEIGVC(IMAXIND,IXSTATE,1))
          IF( IMAXIND .eq. IMAXNTWO ) XLARGE = 
     &                                ABS(XEIGVC(IMAXIND,IXSTATE,2))
          IMAXNONE2 = IFIND_NELM(XEIGVC(1,IXSTATE,1),NXSTATES,
     &                           ABS(XXMDUM),IMAXNONE)
          IMAXNTWO2 = IFIND_NELM(XEIGVC(1,IXSTATE,2),NXSTATES,
     &                           ABS(XXMDUM),IMAXNTWO)
          IMAXIND2  = IMAXNONE2
          IF( ABS(XEIGVC(IMAXNTWO2,IXSTATE,2)) .gt.
     &        ABS(XEIGVC(IMAXNONE2,IXSTATE,1))) IMAXIND2 = IMAXNTWO2
          IF( OMLISTREOD(IMAXIND2) .eq. XMDUM )THEN
            OMLISTREOD(IMAXIND2) = XEIGVL(IXSTATE)
          ELSE
            WRITE(LUPRI,'(/A,I5)')
     &    ' *** ERROR in GET_XOMEGALIST *** failed to assign '//
     &    ' eigenvector == ',IXSTATE
            WRITE(LUPRI,'(/A,2I5)')
     &    ' *** IMAXIND and IMAXIND2 are:',IMAXIND,IMAXIND2
            WRITE(LUPRI,'(/A)') ' OMLISTREOD so far'
            CALL WRTMATMN(OMLISTREOD,1,NXSTATES,1,NXSTATES,LUPRI)
            CALL QUIT('* ERROR in GET_XOMEGALIST * assignment of 
     &                an Omega quantum number to a given state failed')
          END IF
        END IF
      END DO
C
C     list of eigenstates is filled - print
      DO I = 1, NXSTATES
        WRITE(LUPRI,'(I11,8X,A4,6X,1F16.10)')
     &                I,STATEREP,OMLISTREOD(I)
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Function IFIND_NELM */
      INTEGER FUNCTION IFIND_NELM(AVEC,NX,XMAX,MXINDEX)
C***********************************************************************
C
C     find 2nd largest element in AVEC.
C
C     Written by S. Knecht - Dec 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AVEC(NX)
C
      IFIND_NELM = 0
      INDEX_SCR  = 0
      XVAL       = 0.0D0
      DO I = 1, NX
        IF( I .ne. MXINDEX ) THEN
          XVAL = AVEC(I)
          IF( ABS(XVAL) .ge. XMAX ) INDEX_SCR = I
        END IF
      END DO
      IFIND_NELM = INDEX_SCR
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZERO_KAPPA_MJ(XKAPPA,IMJ2,IPRINT)
C***********************************************************************
C
C     Written by S. Knecht - Jan 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
      DIMENSION XKAPPA(NORBT,NORBT,NZ), IMJ2(NORBT)
C
      DO IFRP = 1,NFSYM
        DO IZ = 1, NZ
          DO J = 1, NORB(IFRP)
            DO I = 1, NORB(IFRP)
        IF( IMJ2(IORB(IFRP)+I) .ne. IMJ2(IORB(IFRP)+J)) 
     &  XKAPPA(IORB(IFRP)+I,IORB(IFRP)+J,IZ) = 0.0D0
            END DO
          END DO
        END DO
      END DO
      END
C=========================================================================
C end of file krmcgp.F
C=========================================================================
