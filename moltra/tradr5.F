C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDRV5 */
      SUBROUTINE TRDRV5(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                  NSTR,INDX,TRIAN,KQ,KE,KIBE,
     &                  IC,I2TYP,IBSTRT,IBEND,DRIJ,GABRAO,DINTSKP)
C
C     This is a modified version of TRDRV5 that makes the
C     Lagrangian in a mixed AO,MO form.
C
C     Formulas :
C
C     xxxx
C
C     Written by Joost van Stralen & Luuk Visscher, march 2003
C     Screening by T.Saue Sep 8 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "maxaqn.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxorb.h"
#include "dcbgen.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "twosta.h"
      LOGICAL NOPV, NODV, TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION IJP12(0:7)
      DIMENSION IHTM(0:7),MHTM(2),NHTM(2), IOFFB(0:7)
      DIMENSION IRECC(2),NBFSZ(2),NBUF(2),NIJBUF(2),NHTIM(2),IBOFF(2)
      DIMENSION WORK(*),GABRAO(*),DRIJ(*),DINTSKP(*)
C
      CALL QENTER('TRDRV5')
      KFRSAV = KFREE
C
C     We transform all shells in the first halftransformation
C
      CALL ICOPY ( 32,NSPCK,1,NSPCK12,1)
      CALL ICOPY (196,ISPCK,1,ISPCK12,1)
      CALL ICOPY ( 24,NBBAS,1,NBBAS1, 1)
      CALL ICOPY ( 24,NBBAS,1,NBBAS2, 1)
      CALL ICOPY ( 24,IBBAS,1,IBBAS1, 1)
      CALL ICOPY ( 24,IBBAS,1,IBBAS2, 1)
      CALL ICOPY (  2,IBAS, 1,IBAS1,  1)
      CALL ICOPY (  2,IBAS, 1,IBAS2,  1)
C
C     Copy the first half of the NDMOQ and ICMOQ arrays
C
      CALL ICOPY (8,NDMOQR,1,NDMOQS,1)
      CALL ICOPY (8,NDMOQR,1,NDMOQT,1)
      CALL ICOPY (4,ICMOQR,1,ICMOQS,1)
      CALL ICOPY (4,ICMOQR,1,ICMOQT,1)
C
C     Set up symmetry packing for reduced shell space in the second
C     half transformation. We have one index that has a range of 
C     active shells and one index that should belong to a particular
C     shell (ISHLA).
C
      CALL GISPCK(NSPCK34,ISPCK34,ISHLA,IBSTRT,IBEND,
     &     NBBAS3,NBBAS4,IBBAS3,IBBAS4,IBAS3,IBAS4)
C
C     Generate pointers for the reduced coefficient matrices
C     As we leave the last index untransformed we need only
C     one coefficient matrix and not two as in TRDRV4.
C     Also we need only the square form as a quater transformation
C     cannot be done on triangular arrays.
C
      NQ4S = 0
      DO IFSYM = 1, NFSYM
         NROW4 = 0
         DO IREPR = 0, NBSYM-1 
            IF (IFSYM.EQ.JBTOF(IREPR,IC)) THEN
               NROW4 = NROW4 + NBBAS4(IREPR,IC)
            ENDIF
         ENDDO
C
         NDMOQS(1,IFSYM,4) = NROW4
         NDMOQS(2,IFSYM,4) = NDMOQR(2,IFSYM,4)
C
         ICMOQS(IFSYM,4) = NQ4S + 1
         NQ4S = NQ4S + NROW4 * NDMOQR(2,IFSYM,4) * NZ
      ENDDO
C
C     Get the reduced coefficient matrices (only active shells)
C
      CALL MEMGET('REAL',KQ4S,NQ4S,WORK,KFREE,LFREE)
      DO IFSYM = 1, NFSYM
        IF(NDMOQR(2,IFSYM,4)*NDMOQS(2,IFSYM,4).NE.0) THEN
          CALL QREDUC (IBSTRT,IBEND,IFSYM,WORK(KQ(IFSYM,4)),
     &                NDMOQR(1,IFSYM,4),NDMOQR(2,IFSYM,4),
     &                WORK(KQ4S+ICMOQS(IFSYM,4)-1),
     &                NDMOQS(1,IFSYM,4),NDMOQS(2,IFSYM,4))
        ENDIF
      ENDDO
C
C     Make index array for half-transformed integrals
C
      NHTMT = 0
      DO IREPIJ = 1, NFSYM
         MHTM(IREPIJ) = 0
         NHTM(IREPIJ) = 0
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IHTM(IREPAB) = MHTM(IREPIJ) + NHTMT
               MHTM(IREPIJ) = MHTM(IREPIJ) + NSPCK34(IREPAB,0)*NZ
            ENDIF
            NHTM(IREPIJ) = MHTM(IREPIJ)*NFPCK12(IREPIJ)
         ENDDO
         NHTMT = NHTMT + NHTM(IREPIJ)
      ENDDO
C
      IF (I2TYP.LT.0) THEN
         GOTO 999
      ELSEIF (I2TYP.EQ.1) THEN
         ICS12 = 1
         ICF12 = 1
      ELSEIF (I2TYP.EQ.2) THEN
         ICS12 = 2
         ICF12 = 2
      ELSE
         ICS12 = 1
         ICF12 = 2
      ENDIF
C
C     Memory requirements for HERMIT:
C
      MWHER = MXMEMHER(IC,I2TYP)
C
C     Loop over the shells
C
      CALL IZERO(IOFFB,8)
      DO ISHLB = IBSTRT, IBEND
C
         KFRSAV3 = KFREE
         NINSHA  = NINSH(ISHLA,-1,INDX,IPRINT)
         NINSHB  = NINSH(ISHLB,-1,INDX,IPRINT)
         IF (IPRINT .GE. 3) WRITE(LUPRI,'(/A,2I3)')
     &      '    ---> TRDRV5 : integrals for shell A/B',ISHLA,ISHLB
C
         LDXAB = 5*NINSHA*NINSHB
         CALL MEMGET('INTE',KDXAB  ,LDXAB  ,WORK,KFREE,LFREE)
C
C           Set up gather array for distributions, 
C           calculate memory requirements for this batch.
C
C           NOTE : The actual dimension of INDXAB is kept in the
C           common block dcbtra.h. When we want to parallelize it may
C           be better to allocate it always as MXINSH*MXINSH*4,
C           or (even better) to send it also down to HERMIT.
C           Luuk. TODO
C
C           First a dummy run to calculate the memory requirements
C           for GMAT. We assume that HMAT is not more than twice GMAT,
C           and reserve the same memory also for use in HERMIT.
C
C        The size of HMAT is approximately: 
C           NUMDIS * NFPCK12(IREPIJ) * NZ / NPASS
C
C        (1) overhead: Subtract 10% or 1MW from available memory 
C
         MFREE = LFREE - MIN(LFREE/10, 1 000 000)
C
C        (2) We must be able to hold
C            (a) GMAT + HERMIT
C            (b) GMAT + HMAT
C            simultaneously in memory.
C 
C
         IF (MWHER .GT. MFREE) CALL QUIT('Not enough memory for HERMIT')
C            
         MWHMAT = 0
         IT = 0
 2000    CONTINUE
            IT = IT + 1
            NSIZG = MFREE - MAX(MWHMAT,MWHER)
C
C           Calculate NUMDIS and NSIZG
C           --------------------------
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &                       .FALSE.,INDX,WORK(KFREE),IJP12,NSIZG,IDUM)
C
C           Setup INDXAB & IJPASS
C           ---------------------
C
            CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.FALSE.,
     &                       INDX,WORK(KDXAB),IJP12,NSIZG,WORK(KIJPASS))
C
C           Calculate MWHMAT:
C           -----------------
C
            MWHMAT = MWHMATCALC(WORK(KIJPASS),NPASS,NFPCK12)
C
            CALL MEMREL('TRDRV5.mwhmat calc',WORK,1,KIJPASS,KFREE,LFREE)
C
C           Check that we fit in memory
C           ---------------------------
C
            IF (IPRINT .GE. 3) THEN
               WRITE(LUPRI,'(A,I3,A/4(A,I9))')
     &         '       - TRDRV5 memory estimate iteration no. ',IT,' :',
     &         '         Free: ',MFREE,', Hermit: ',MWHER,
     &         ', HMAT ',MWHMAT,', GMAT ',NSIZG
            END IF
            CALL FLSHFO(LUPRI)
            IF ( (NSIZG + MAX(MWHMAT,MWHER) ) .GT. MFREE ) THEN
               IF ( MWHMAT .GT. MFREE) THEN
                  MWHMAT = MFREE * 9 / 10
               END IF
               GOTO 2000
            END IF
C
C        Final calculation of INDXAB & IJPASS
C
         CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
         NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.FALSE.,
     &                    INDX,WORK(KDXAB),IJP12,NSIZG,WORK(KIJPASS))
C
         IF (IPRINT .GE. 3) THEN
            WRITE (LUPRI,'(//,4(1X,A,2I4,/))')
     &      ' Calculation of integral distribution set:',ISHLA,ISHLB
            WRITE (LUPRI,'(//,4(1X,A,I8,/))')
     &      ' Number of distributions in this set:     ',NUMDIS,
     &      ' Maximum size of scalar integral batch:   ',NSIZG,
     &      ' Number of passes :                       ',NPASS
         END IF
C
         IF (NPASS.GT.1) THEN
C
C           Initialize scalar integral buffer file and arrays
C           Use the space reserved for GMAT now as buffer storage
C
            LGFIL = LUTRA2
            NGBFSZ = (IRAT*NSIZG)/((1+IRAT)*NPASS)
C
C           Allocate the buffer arrays in WORK and set KGMAT to the
C           start of this allocation section. We now only need to
C           send down GMAT into HERMIT to find this section.
C
            CALL MEMGET('INTE',KALGREC,NPASS,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KALGBUF,NPASS,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KAIGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KARGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
C
C           As we refer to the places inside the section we need
C           to subtract KGMAT.
C
            KGMAT  = KALGREC
            KLGREC = KALGREC - KGMAT + 1
            KLGBUF = KALGBUF - KGMAT + 1
            KIGBUF = KAIGBUF - KGMAT + 1
            KRGBUF = KARGBUF - KGMAT + 1
         ELSE
C
C           Write directly into the array GMAT
C
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
            KLGBUF = 1
            KALGREC = KFREE
         ENDIF
C
C        Get distributions
C
         NODV = NASHT.EQ.0
         NOPV = NASHT.LT.2
         CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &               NODV,NOPV,GABRAO,DRIJ,DINTSKP,SCRTRA,
     &               WORK(KFREE),LFREE,IPRINT)
C
C        Do the transformation of the first index pair and continue
C        by contracting these integrals with the T-amplitudes to obtain
C        the Lagrangian contributions.
C        half-transformed integrals to file.
C
C        De-allocate buffers and allocate GMAT if we had multiple passes.
C
         IF (NPASS.GT.1) THEN
            CALL MEMREL('TRDRV5.1',WORK,1,KALGBUF,KFREE,LFREE)
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
         ENDIF
C
         CALL MS3RLAG(WORK,KFREE,LFREE,IPRINT,.TRUE.,.TRUE.,IOFFB,
     &                WORK(KALGREC),ICS12,ICF12,IC,
     &                TRIAN(1),NSTR(1,0,1),NSTR(1,0,2),
     &                WORK(KIJPASS),WORK(KQ(1,1)),WORK(KQ(1,2)),
     &                WORK(KGMAT),NSIZG,
     &                MHTM,IHTM,INDX,WORK(KDXAB),IRECC,WORK(KJBOFF),
     &                IBOFF,NBUF,NIJBUF,NBFSZ,WORK(KRBUF),
     &                WORK(KLBUF1))
C
         CALL MEMREL('TRDRV5.2',WORK,1,KFRSAV3,KFREE,LFREE)
C
C        Delete scalar integral buffer file
C
         IF (NPASS.GT.1) THEN
            CALL DELGBUF(LGFIL)
         ENDIF
C
         DO IREPB = 0,NBSYM-1
            IOFFB(IREPB) = IOFFB(IREPB) + NINSH(ISHLB,IREPB,INDX,IPRINT)
         END DO
      ENDDO
C
C     Release all memory and exit
C
 999  CONTINUE
      CALL MEMREL('TRDRV5.5',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDRV5')
C
      RETURN
      END
