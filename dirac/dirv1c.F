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
C...   group of the licence holder may not take any part of Dirac,
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
C FILE    : dirv1c.F
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type1drv */
      SUBROUTINE TYPE1DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,
     &                     QMAT,IPRINT) 
C*****************************************************************************
C
C    TYPE1DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the smlv1c type1 model.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2)
#include <dcbbas.h>
#include <dcbham.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbisol.h>
C
      CALL QENTER('TYPE1DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
      IF(DOINT(1,2) .OR. DOINT(2,1) ) THEN
         CALL QUIT ('TYPE1DRV called with LS or SL type integrals.')
      END IF
#include <memint.h>
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ... set charges. 
c         (We use the nuclear exponents)
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      IF(DOINT(1,1)) THEN
         IF (ICTLV1C(2) .EQ. 1) THEN 
            CALL SCDEN2(CHARGE,WORK(KCHRG),NUCIND,0)
         ELSE
            IOFF = 1
            DO I = 1,NUCIND
               IF (IOFF .GT. NUCDEP) 
     &             CALL QUIT('TYPE1DRV : IOFF .GT. NUCDEP')
               DO J = 1,NUCDEP
                  WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,2)
               END DO
               IOFF = IOFF + FMULT(ISTBNU(I))
            ENDDO
         END IF
      ELSE IF(DOINT(2,2)) THEN
         DO I = 1,NUCIND
            WORK(KCHRG-1+I) = -CHARGE(I)
         ENDDO
      END IF
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from TYPE1DRV',-1)
        WRITE(LUPRI,'(/1X,A)')
     &       '* Charges used for correction  :'
              WRITE(LUPRI,'(4X,A4,F12.8)') 
     &       (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
      ENDIF
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            GNUEXP,WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL MEMCHK('Before_symupk',WORK,KWORK,KFREE,LFREE)
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE1DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND() + 1.0D0
      TIME   = TIMEND - TIMSTR
      TIMSTR = TIMEND
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE1DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE1DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type2drv */
      SUBROUTINE TYPE2DRV(OP1INT,WORK,LWORK,DOINT,OMITVNUC,IPRINT) 
C*****************************************************************************
C
C    TYPE2DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the smlv1c type 2 model.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include <dcbbas.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbisol.h>
C
      CALL QENTER('TYPE2DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
#include <memint.h>
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KVCMAT),NNBBASX)
c   
c     ... set charges and gaussian exponents.
c
      CALL DZERO(WORK(KCHRG),NUCIND)
      DO I = 1,NUCIND
         WORK(KXGEXP-1+I) = GNUEXP(I)
         WORK(KCHRG-1+I) = -CHARGE(I)
      ENDDO
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from TYPE2DRV',-1)
        WRITE(LUPRI,'(/1X,A)')
     &       '* Charges used for correction  :'
              WRITE(LUPRI,'(4X,A4,F12.8)') 
     &       (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
      ENDIF
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals
c
      CALL CATDR1(WORK(KVCMAT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL MEMCHK('Before_symupk',WORK,KWORK,KFREE,LFREE)
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE2DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND() + 1.0D0
      TIME   = TIMEND - TIMSTR
      TIMSTR = TIMEND
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE2DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE2DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3drv */
      SUBROUTINE TYPE3DRV(OP1INT,WORK,LWORK,DOINT,ADDQIC,
     &                    OMITVNUC,IPRINT,QMAT)
C*****************************************************************************
C
C    TYPE3DRV : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the smlv1c type3 model.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
C      PARAMETER (D0 = 0.0D0, QTHRS = 1.0D-10)
      PARAMETER (D0 = 0.0D0, QTHRS = 1.0D-8)
      LOGICAL DOINT(2,2), ADDQIC(2)
      CHARACTER*4 OMITVNUC
      DIMENSION OP1INT(NNBBASX),WORK(LWORK),QMAT(NUCDEP,NUCDEP,2,*)
#include <dcbbas.h>
#include <dcbham.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbisol.h>
C
      CALL QENTER('TYPE3DRV')
c
c     Start timer
c
      TIMSTR = SECOND()
#include <memint.h>
C
      CALL MEMGET('REAL',KVCMAT,  NNBBASX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KVCMAT),NNBBASX)
C
C     Initialization of charges and coordinates
C
      NNDEP = NUCDEP*(NUCDEP+1)/2
      IF(ICTLV1C(2).EQ. 2) THEN
C
C     ... Atomic centered mulliken charges
C
         MXCENTC = NUCDEP
         KCENTC  = NUCDEP
      ELSE 
C
C     ... Distributed charges
C
         CALL TYPE3MXC(MXCENTC,NNDEP,QMAT,ADDQIC,QTHRS)
         KCENTC  = 0
      END IF
      IF (MXCENTC.LE.0) CALL QUIT('TYPE3DRV: MXCENTC.le.0')
      LFCBA = MXCENTC*NUCDEP**2
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  MXCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*MXCENTC,WORK,KFREE,LFREE)
C     ... these are not used now, should be removed.
C         TODO: (write new AVECA1 etc. instead of AVENA1 etc.)
      CALL MEMGET('REAL',KSIGNC,3*KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  KCENTC,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  KCENTC,WORK,KFREE,LFREE)
C
      IF(ICTLV1C(2).EQ. 2) THEN
C         ... Get charges and exponents
         CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KCHRG),NUCIND)
         IOFF = 1
         DO I = 1,NUCIND
            IF (IOFF .GT. NUCDEP) 
     &          CALL QUIT('TYPE1DRV : IOFF .GT. NUCDEP')
            DO J = 1,NUCDEP
               IF (ADDQIC(1))
     &             WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,1,1)
               IF (ADDQIC(2))
     &             WORK(KCHRG-1+I) = WORK(KCHRG-1+I) + QMAT(J,IOFF,2,1)
            END DO
            IOFF = IOFF + FMULT(ISTBNU(I))
            WORK(KXGEXP-1+I) = GNUEXP(I)
Ctest        ... finite size exponents for q_L (for Iodine,Carbon,Hydrogen)
C            IF(DOINT(2,2)) THEN
C               IF (NINT(CHARGE(I)) .EQ. 53)   WORK(KXGEXP-1+I) = 0.2D0
C               IF (NINT(CHARGE(I)) .EQ.  6)   WORK(KXGEXP-1+I) = 1.0D0
C               IF (NINT(CHARGE(I)) .EQ.  1)   WORK(KXGEXP-1+I) = 1.5D0
C            END IF
         ENDDO
         IF(IPRINT.GE.3) THEN
            CALL HEADER('Output from TYPE1DRV',-1)
            WRITE(LUPRI,'(/1X,A)')
     &         '* Charges used for correction  :'
                WRITE(LUPRI,'(4X,A4,F12.8)')
     &        (NAMN(I),WORK(KCHRG-1+I),I = 1,NUCIND)
         ENDIF
         CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &               WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),
     &               WORK(KJSYMC),WORK(KJCENT),OMITVNUC,IPRINT)
         MXCENTC = NCENTC
      ELSE
         CALL MEMGET('REAL',KCOORA,3*NUCDEP,WORK,KFREE,LFREE)
         CALL GETACORD(WORK(KCOORA))
         CALL TYPE3DR1(MXCENTC,WORK(KCOORC),WORK(KCOORA),WORK(KGEXP),
     &                QMAT,ADDQIC,QTHRS,WORK(KFCBA),OMITVNUC,
     &                NNDEP,IPRINT)
         IF(IPRINT .GE. 2) THEN
            CALL MEMGET('REAL',KACHRG,2*NUCIND,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KACHRG),2*NUCIND)
C
C           ... Calculate atomic charges
C
            IOFF = 1
            QLT  = D0
            QST  = D0
            DO I = 1,2*NUCIND,2
               IF (IOFF .GT. NUCDEP) 
     &            CALL QUIT('TYPE3DRV : IOFF .GT. NUCDEP')
                  DO J = 1,NUCDEP
                     WORK(KACHRG-1+I)=WORK(KACHRG-1+I)+QMAT(J,IOFF,1,1)
                     WORK(KACHRG  +I)=WORK(KACHRG  +I)+QMAT(J,IOFF,2,1)
                  END DO
               IOFF = IOFF + FMULT(ISTBNU((I+1)/2))
               QLT = QLT + WORK(KACHRG-1+I)*FMULT(ISTBNU((I+1)/2))
               QST = QST + WORK(KACHRG  +I)*FMULT(ISTBNU((I+1)/2))
            ENDDO
            WRITE(LUPRI,'(/A/)') 
     &      '  ** Total atomic charges (Mulliken) **'
            WRITE(LUPRI,'(A/)') 
     &      '  Atom                 Q_l                     Q_s'
            WRITE(LUPRI,'(4X,A4,1P,2D25.12)') 
     &      (NAMN((I+1)/2),WORK(KACHRG-1+I),WORK(KACHRG+I),
     &      I = 1,2*NUCIND,2)
            WRITE(LUPRI,'(2X,A)') 
     &      '--------------------------------------------------------'
            WRITE(LUPRI,'(A,1P,2D25.12)') 
     &      '  Sum   ',QLT,QST
         END IF
         CALL MEMREL('TYPE3DRV.TYPE3DR1',WORK,KWORK,KCOORA,KFREE,LFREE)
C
      END IF
C
C    ... Calculate integrals
C
      CALL CATDR1(WORK(KVCMAT),DUMMY,MXCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LFREE,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack nuclear attraction integrals
C
      CALL SYMUPK(WORK(KVCMAT),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('TYPE3DRV.CATDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
c     Stop timer
c
      TIMEND = SECOND() + 1.0D0
      TIME   = TIMEND - TIMSTR
      TIMSTR = TIMEND
      WRITE (LUPRI,'()')
      CALL TIMTXT('>>> Time used in TYPE3DRV is',TIME,LUPRI)
      CALL QEXIT('TYPE3DRV')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3mxc */
      SUBROUTINE TYPE3MXC(MXCENTC,NNDEP,QMAT,ADDQIC,QTHRS)
C*****************************************************************************
C
C    TYPE3MXC : Find MXCENTC, the number of sites needed to
C               calculate smlv1c type3 correction.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER (D0 = 0.0D0 ,  D2 = 2.0D0)
      LOGICAL ADDQIC(2)
      DIMENSION QMAT(NUCDEP,NUCDEP,2)
      CALL QENTER('TYPE3MXC')
C
      MXCENTC = 0
      IJ = 0
      DO J = 1, NUCDEP
         DO I = 1,J
            IJ = IJ + 1
            TMP = D0
            IF(ADDQIC(1)) TMP = QMAT(I,J,1) + TMP
            IF(ADDQIC(2)) TMP = QMAT(I,J,2) + TMP
            IF(I.NE.J) TMP = D2*TMP
            IF (ABS(TMP) .GT. QTHRS) MXCENTC = MXCENTC + 1
         ENDDO
      ENDDO
      CALL QEXIT('TYPE3MXC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getacord */
      SUBROUTINE GETACORD(COORA)
C*****************************************************************************
C
C    GETACORD : Make list atomic coordinates
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <symmet.h>
#include <pgroup.h>
      PARAMETER (D0 = 0.0D0 , DP5 = 0.5D0 , D2 = 2.0D0)
      DIMENSION COORA(3,*)
      CALL QENTER('GETACORD')
C
C     Make the full matrix of cartisian coordinates from CORD(NUCIND)
C
      JATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         IF (MULT(MULCNT) .EQ. 1) THEN
            JATOM = JATOM + 1
            COORA(1,JATOM) = CORD(1,ICENT)
            COORA(2,JATOM) = CORD(2,ICENT)
            COORA(3,JATOM) = CORD(3,ICENT)
         ELSE
            DO 200 ISYMOP = 0, MAXOPR
               IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  COORA(1,JATOM) = 
     &                 PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                  COORA(2,JATOM) = 
     &                 PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                  COORA(3,JATOM) = 
     &                 PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
      CALL QEXIT('GETACORD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck type3dr1 */
      SUBROUTINE TYPE3DR1(MXCENTC,COORC,COORA,XGNUEXP,QMAT,
     &                   ADDQIC,QTHRS,FACCBA,OMITVNUC,NNDEP,IPRINT)
C*****************************************************************************
C
C    TYPE3DR1 : Driver for coulomb attraction integrals,
C               used to mimic the multicenter integrals
C               neglected in the smlv1c type3 model.
C               Centers with charges smaller than QTHRS
C               are neglected.
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <dcbham.h>
      PARAMETER (D0 = 0.0D0 , D2 = 2.0D0 , DP5 = 0.5D0)
      DIMENSION FACCBA(MXCENTC,NUCDEP,NUCDEP),QMAT(NUCDEP,NUCDEP,2,*),
     &          XGNUEXP(MXCENTC),COORC(3,MXCENTC),COORA(3,NUCDEP)
      LOGICAL NO_CCC, NO_CAA,NO_CCA,NO_CBA, ADDQIC(2)
      CHARACTER*4 OMITVNUC
      CALL QENTER('TYPE3DR1')
C     <A | C | B> = CAB
      NO_CCC = (OMITVNUC(1:1) .EQ. 'T')
      NO_CAA = (OMITVNUC(2:2) .EQ. 'T')
      NO_CCA = (OMITVNUC(3:3) .EQ. 'T')
      NO_CBA = (OMITVNUC(4:4) .EQ. 'T')
      DO I = 1,MXCENTC
C     Initialize to point charge distributions
C     (change later if disired )
         XGNUEXP(I) = D0
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('Output from TYPE3DR1',-1)
         WRITE(LUPRI,'(/A)') 
     &        '  ** Distribution of charges  **'
         IF(ICTLV1C(2).EQ.0) THEN
            WRITE(LUPRI,'(A/)') 
     &      '     (Charges placed where dipole moment is zero)'
         ELSE IF(ICTLV1C(2).EQ.3) THEN
            WRITE(LUPRI,'(A/)')
     &      '     (Charges placed at midpoints between nuclei)'
         ENDIF
         WRITE(LUPRI,'(/,4X,A,10X,A,17X,A,24X,A,24X,A,/)') 
     &        'NCENTC','-CHARGE(NCENTC)','X','Y','Z'
      END IF
C
      N2DEP = NUCDEP*NUCDEP
      QSKIP = D0
      NCENTC = 0
      IJ = 0
      DO J = 1, NUCDEP
C
C        off-diagonal terms
C
         DO I = 1,J-1
            IJ = IJ + 1
            QTMP = D0
            IF(ADDQIC(1)) QTMP = QTMP - QMAT(I,J,1,1)
            IF(ADDQIC(2)) QTMP = QTMP - QMAT(I,J,2,1)
            QTMP = D2*QTMP
            IF (ABS(QTMP) .GT. QTHRS) THEN
               NCENTC = NCENTC + 1
               IF (ICTLV1C(2).EQ.0) THEN
                  DO K = 1,3
                     CTMP = D0
                     IF(ADDQIC(1)) CTMP = CTMP - QMAT(I,J,1,1+K)
                     IF(ADDQIC(2)) CTMP = CTMP - QMAT(I,J,2,1+K)
                     COORC(K,NCENTC) = D2 * CTMP / QTMP
                  END DO
C                  IF (ADDQIC(2)) XGNUEXP(NCENTC) = 100.0D0
C                  IF (ADDQIC(1)) XGNUEXP(NCENTC) = 2.0D0
               ELSE IF (ICTLV1C(2).EQ.3) THEN
                  COORC(1,NCENTC) = DP5*(COORA(1,I) + COORA(1,J))
                  COORC(2,NCENTC) = DP5*(COORA(2,I) + COORA(2,J))
                  COORC(3,NCENTC) = DP5*(COORA(3,I) + COORA(3,J))
               END IF
               IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,5X,1P,4D25.12)')
     &         IJ,QTMP,COORC(1,NCENTC),COORC(2,NCENTC),COORC(3,NCENTC)
C
               CALL DCOPY(N2DEP,QTMP,0, FACCBA(NCENTC,1,1),MXCENTC)
C
C              Zero elements of FACCBA which should not be included
C              -- all off-diagonal included in this version
C
               IF (IPRINT.GT.5) THEN
                   WRITE(LUPRI,'(/4X,A)') 'FACCBA for this center:'
                   CALL OUTPUT(FACCBA(NCENTC,1,1),1,1,1,N2DEP,
     &                 MXCENTC,N2DEP,1,LUPRI)
                   WRITE(LUPRI,'()')
               END IF
            ELSE
              IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,1P,D30.12,A)')
     &              IJ,QTMP, '         skipped!'
              QSKIP = QSKIP + ABS(QTMP)
            END IF
C
         END DO
C
C        diagonal terms
C
         IJ = IJ + 1
         QTMP = D0
         IF(ADDQIC(1)) QTMP = QTMP - QMAT(J,J,1,1)
         IF(ADDQIC(2)) QTMP = QTMP - QMAT(J,J,2,1)
         IF (ABS(QTMP) .GT. QTHRS) THEN
            NCENTC = NCENTC + 1
            IF (ICTLV1C(2).EQ.0) THEN
               DO K = 1,3
                  CTMP = D0
                  IF(ADDQIC(1)) CTMP = CTMP - QMAT(J,J,1,1+K)
                  IF(ADDQIC(2)) CTMP = CTMP - QMAT(J,J,2,1+K)
                  COORC(K,NCENTC) = CTMP / QTMP
               END DO
C               IF (ADDQIC(2)) XGNUEXP(NCENTC) = 100.0D0
C               IF (ADDQIC(1)) XGNUEXP(NCENTC) = 2.0D0
            ELSE IF (ICTLV1C(2).EQ.3) THEN
               COORC(1,NCENTC) = COORA(1,J)
               COORC(2,NCENTC) = COORA(2,J)
               COORC(3,NCENTC) = COORA(3,J)
            END IF
            IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,5X,4D25.12)')
     &          IJ,QTMP,COORC(1,NCENTC),COORC(2,NCENTC),COORC(3,NCENTC)
C
            CALL DCOPY(N2DEP,QTMP,0, FACCBA(NCENTC,1,1),MXCENTC)
C
C        Zero elements which should not be included
C
            IF (NO_CCC) FACCBA(NCENTC,J,J) = D0
            IF (NO_CAA) THEN
               DO I = 1, NUCDEP
                  IF(I.NE.J) FACCBA(NCENTC,I,I) = D0
               END DO
            END IF
            IF (NO_CCA) THEN
               DO I = 1, NUCDEP
                  IF (I.NE.J) THEN
                     FACCBA(NCENTC,J,I) = D0
                     FACCBA(NCENTC,I,J) = D0
                  END IF
               END DO
            END IF
            IF (NO_CBA) THEN
               DO K = 1, NUCDEP
                  IF (K.NE.J) THEN
                     DO I = 1, NUCDEP
                           IF (I.NE.J) FACCBA(NCENTC,I,K) = D0
                     END DO
                  END IF
               END DO
            END IF
            IF (IPRINT.GT.5) THEN
                WRITE(LUPRI,'(/4X,A)') 'FACCBA for this center:'
                CALL OUTPUT(FACCBA(NCENTC,1,1),1,1,1,N2DEP,
     &                      MXCENTC,N2DEP,1,LUPRI)
                WRITE(LUPRI,'()')
            END IF
         ELSE
           IF (IPRINT.GE.2) WRITE(LUPRI,'(4X,I3,1P,D30.12,A)')
     &         IJ,QTMP, '         skipped!'
           QSKIP = QSKIP + ABS(QTMP)
         END IF
      END DO
C
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(/A,I17)') 
     &               '     Number of charges skipped : ', NNDEP-NCENTC
         WRITE(LUPRI,'(/A,1P,D17.9)') 
     &               '     Amount of charge skipped  : ', QSKIP
         WRITE(LUPRI,'(A/)') 
     &               '     ( sum_i |Q_skip(i)| )'
      ENDIF
      IF (NCENTC.NE.MXCENTC) CALL QUIT('TYPE3DR1 : NCENTC .NE. MXCENTC')
      CALL QEXIT('TYPE3DR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck catcorr */
      SUBROUTINE CATCORR(F1MAT,F2MAT,DMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     CATCORR : Generate SS nuclear attraction integrals to mimic
C     two-electron integrals for the smlv1c models.
C
C     Written by hjaaj & jkp   -   oct. 2001
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL DOINT(2,2),ADDQIC(2)
      CHARACTER*4 OMITVNUC
      DIMENSION WORK(LWORK),F1MAT(*),F2MAT(*),DMAT(*)
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbham.h>
#include <cbihr1.h>
C
      CALL QENTER('CATCORR')
#include <memint.h>
C
C     Initialization
C
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .FALSE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
C
C     Memory allocation
C
      IF(INTV1C.EQ.3 .AND. ICTLV1C(2) .EQ. 0) THEN
          NMPOLOP = 4
      ELSE IF (INTV1C.EQ.2) THEN
          NMPOLOP = 0
      ELSE
         NMPOLOP = 1
      END IF
      CALL MEMGET('REAL',KQMAT,NUCDEP*NUCDEP*2*NMPOLOP,
     &             WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEMT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
c     ... Get matrix of small- and large-component charges
C
      KREL = KFREE
      IF(INTV1C.NE.2) THEN
        IF(NOPEN.GT.1) THEN
          CALL QUIT(
     &   'open shell not implemented correctly in CATCORR yet')
          CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
          CALL DCOPY(N2BBASX,DMAT,1,WORK(KBUF),1)
          DO ISHELL = 1, NOPEN
            JDMAT = 1+N2BBASXQ*ISHELL
            CALL DAXPY(N2BBASX,DF(ISHELL),DMAT(JDMAT),1,
     &                 WORK(KBUF),1)
         ENDDO
          CALL LSQMAT(WORK(KQMAT),WORK(KBUF),NMPOLOP,WORK(KFREE),
     &               LFREE,IPRINT)
        ELSE
          CALL LSQMAT(WORK(KQMAT),DMAT,NMPOLOP,WORK(KFREE),LFREE,
     &              IPRINT)
        ENDIF
      ENDIF
c
C     Omit selected terms depending on OMITVNUC
c     (see CATNUC-routine (ababus/her1cat.F) for detail.)
c
      OMITVNUC    = 'FFFF'
c
c     determine model
c
      TIMSTR = SECOND()
      IF (SMLV1C) THEN
         IF (INTV1C .EQ. 1) THEN
c           ... SMLV1C type 1
c           ... we need to correct all but < C | V_C | C > type integrals
c
c           ... correct <L|Vnuc|L> LL-type integrals
c
            DOINT(1,1)  = .TRUE.
            OMITVNUC    = 'TFFF'
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE1DRV FOR LL-TYPE INTEGRALS  >>>'
            CALL TYPE1DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate LL energy correction (Charge = q_s)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL ESMLV1C(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,1)
c
c            correct <S|Vnuc|S> SS-type integrals
c
            DOINT(1,1)  = .FALSE.
            DOINT(2,2)  = .TRUE.
            OMITVNUC    = 'TFFF'
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE1DRV FOR SS-TYPE INTEGRALS  >>>'
            CALL MEMGET('REAL',KONEIN2,NNBBASX,WORK,KFREE,LFREE)
            CALL TYPE1DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,WORK(KQMAT),IPRINT)
c
c           ... calculate SS energy correction (Charge = -Z_nuc)
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL ESMLV1C(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
            CALL DAXPY(NNBBASX,D1,WORK(KONEIN2),1,WORK(KONEIN),1)
            CALL MEMREL('CATCORR.SMLV1C1',WORK,KWORK,KONEIN2,KFREE,
     &                   LFREE)
         ELSE IF (INTV1C .EQ. 3) THEN
c           ... SMLV1C type 3
c           ... we need to correct all but < C | V_C | C > type integrals
c           ... correct <L|Vnuc|L> LL-type integrals
            OMITVNUC = 'TFFF'
            DOINT(1,1)  = .TRUE.
            DOINT(2,2)  = .FALSE.
            ADDQIC(1) = .FALSE.
            ADDQIC(2) = .TRUE.
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE3DRV FOR LL-TYPE INTEGRALS  >>>'
            CALL TYPE3DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  ADDQIC,OMITVNUC,IPRINT,WORK(KQMAT))
c
c           ... calculate energy correction due to q_s
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL ESMLV1C(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,1)
c           ... correct <S|Vnuc|S> SS-type integrals
            OMITVNUC = 'TFFF'
            DOINT(1,1)  = .FALSE.
            DOINT(2,2)  = .TRUE.
            ADDQIC(1) = .TRUE.
            ADDQIC(2) = .TRUE.
            CALL MEMGET('REAL',KONEIN2,NNBBASX,WORK,KFREE,LFREE)
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE3DRV FOR SS-TYPE INTEGRALS  >>>'
            CALL TYPE3DRV(WORK(KONEIN2),WORK(KFREE),LFREE,DOINT,
     &                  ADDQIC,OMITVNUC,IPRINT,WORK(KQMAT))
c
c           ... calculate energy correction 
c
            CALL DSPTSI(NTBAS(0),WORK(KONEIN2),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL ESMLV1C(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
            CALL DAXPY(NNBBASX,D1,WORK(KONEIN2),1,WORK(KONEIN),1)
            CALL MEMREL('CATCORR.TYPE3',WORK,KWORK,KONEIN2,KFREE,
     &                   LFREE)
         ELSE
c           ... SMLV1C type 2
c           ... we need to correct < X | V_C | Y > type integrals, X.ne.Y
c           ... correct SS-type integrals
            DOINT(2,2) = .TRUE.
            OMITVNUC    = 'TTFF'
c
c           ... calculate energy correction due to q_s
c
            WRITE(LUPRI,'(//1X,A)') 
     &      '<<<  ENTERING TYPE2DRV FOR SS-TYPE INTEGRALS  >>>'
            CALL TYPE2DRV(WORK(KONEIN),WORK(KFREE),LFREE,DOINT,
     &                  OMITVNUC,IPRINT)
            CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
            CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
            CALL ESMLV1C(WORK(KONEMT),DMAT,WORK(KFREE),LFREE,2)
         END IF
      ELSE
         CALL QUIT('CATCORR not implemented for .not. SMLV1C')
      END IF
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
      IF(IPRINT.GE.6) THEN
         IF (INTV1C .EQ. 3) THEN
             CALL HEADER('CATCORR: SMLV1C type 3 correction to nuclear 
     &                    attraction integrals',-1)
         ELSE IF (INTV1C .EQ. 1) THEN
             CALL HEADER('CATCORR: SMLV1C type 1 correction to nuclear 
     &                    attraction integrals',-1)
         ELSE
             CALL HEADER('CATCORR: SMLV1C type 2 correction to nuclear 
     &                    attraction integrals',-1)
         END IF
         CALL PRQMAT(WORK(KONEMT),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      ENDIF

C     Add correction to Fock-matrix
C
      IF (INTV1C .EQ. 3) THEN
c        ... SMLV1C type 3 : correct two-electron Fock-matrix
         CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F2MAT,1)
         IF(IPRINT.GE.6) THEN
            CALL HEADER('CATCORR: SMLV1C corrected F2MAT',-1)
            CALL PRQMAT(F2MAT,NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
         END IF
      ELSE
c        ... SMLV1C type 1 or 2 : correct one-electron Fock-matrix
         CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,F1MAT,1)
         IF(IPRINT.GE.6) THEN
            CALL HEADER('CATCORR: SMLV1C corrected F1MAT',-1)
            CALL PRQMAT(F1MAT,NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
CC
CC     Memory deallocation
      CALL MEMREL('CATCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('CATCORR')
      RETURN
      END
C*****************************************************************************
