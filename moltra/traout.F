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
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : traout.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTI */
      SUBROUTINE MDINTI(NSTR,IKRA,IKRB,IPRINT)
C
C     Initialize MDINT file
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION NSTR(2),IKRA(*),IKRB(*)
      CHARACTER*8 CDATE,TIMEX,DATEX*10
C
C     The date is kept in a 10-character variable in RELCCSD
C
      CALL GETDAT (CDATE,TIMEX)
      DATEX = ' '//CDATE//' '
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NSTR(IFSYM)
            II = II + 1
            JJ = JJ + 1
            IKRA(II) = JJ
            IKRB(II) = JJ + NSTR(IFSYM)
         ENDDO
         JJ = JJ + NSTR(IFSYM)
      ENDDO
C           
      OPEN (LUMLF2,FILE='MDCINT',FORM='UNFORMATTED')
      WRITE (LUMLF2) DATEX,TIMEX,II,(IKRA(I),IKRB(I),I=1,II)
      NTKL = 0
      N2EFIL = 0
C
      IF (IPRINT.GE.5) WRITE (LUPRI,'(//3X,A)') 'Integrals written :'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTF */
      SUBROUTINE MDINTF
C
C     Finish MDINT file
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
      WRITE (LUMLF2) 0,0,0
      CLOSE (LUMLF2,STATUS='KEEP')
      WRITE (LUPRI,1000)
C
 1000 FORMAT (//' MOLFDIR file MDCINT is written')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTN */
      SUBROUTINE MDINTN (NZKL)
C
C     Split MDINT file if necessary
C     Written by Luuk Visscher, Spring 2000
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
      CHARACTER*8 NEXTFIL
C
      NTKL = NTKL + NZKL
C     Every integral takes at most 32 bytes on 4 byte integer machines
C     Make a new file if we have written more than 60 million integrals.
C     This solves the 2Gb problem on Linux.
      IF (NTKL.GT.60 000 000) THEN
         N2EFIL = N2EFIL + 1
         IF (N2EFIL.LT.10) THEN
            WRITE (NEXTFIL,'(A7,I1)') 'MDCINT0',N2EFIL
         ELSEIF (N2EFIL.LT.100) THEN
            WRITE (NEXTFIL,'(A6,I2)') 'MDCINT',N2EFIL
         ELSE
            STOP "100 integral files is enough !"
         ENDIF
         WRITE (LUMLF2) 0,N2EFIL,0
         WRITE (LUMLF2) NEXTFIL
         CLOSE (LUMLF2,STATUS='KEEP')
         OPEN (LUMLF2,FILE=NEXTFIL,FORM='UNFORMATTED')
         NTKL = 0
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTB */
      SUBROUTINE MDINTB(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,FMAT,
     &                  NFM,INDK,INDL,RKL)
C
C     Write block of integrals to MDINT file. This is the old MOLFDIR
C     format and should work with both RELCCSD and DIRRCI.
C
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION NSTR3(2),NSTR4(2)
      DIMENSION FMAT(NFM,2),INDK(*),INDL(*)
      COMPLEX*16 RKL(*)
C
      NZKL = 0
      IKR1 = IKR
      JKR1 = JKR
      KL = 0
      LL = 0
      DO IREPL = 1, NFSYM
         IREPK = MOD(IREPL+IREPIJ,2) + 1
         DO L = 1, NSTR4(IREPL)
            LL = LL + 1
            DO K = 1, NSTR3(IREPK)
               IF (IREPK.EQ.1) THEN
                  KK = K
               ELSE
                  KK = K + NSTR3(1)
               ENDIF
               KL = KL + 1
               IF (DABS(FMAT(KL,1)).GT.THROUT.OR.
     &             DABS(FMAT(KL,2)).GT.THROUT) THEN
                  NZKL = NZKL + 1
                  INDK(NZKL) = KK
                  INDL(NZKL) = LL
                  RKL(NZKL) = DCMPLX(FMAT(KL,1),FMAT(KL,2))
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      IF (ICLASS.EQ.1) THEN
         CONTINUE
      ELSEIF (ICLASS.EQ.2) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.3) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.4) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.9) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.10) THEN
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.11) THEN
         IKR1 = - IKR1
      ELSEIF (ICLASS.EQ.12) THEN
         JKR1 = - JKR1
      ELSE
         CALL QUIT ('Class not defined in MDINTB')
      ENDIF
C
C     We may  exceed the maximum file length. Check !
C
      CALL MDINTN (NZKL)
C
      WRITE (LUMLF2) IKR1,JKR1,NZKL,
     &               (INDK(INZ),INDL(INZ),INZ=1,NZKL),
     &               (RKL(INZ),INZ=1,NZKL)
C
C     We should send in the print flag here, and make this an option.
C
      IF (IPRINT.GE.5) THEN
         DO INZ = 1, NZKL
            WRITE (LUPRI,'(4I4,2f30.16)') IKR1,JKR1,INDK(INZ),INDL(INZ),
     &                                    RKL(INZ)
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SYMFINT */
      SUBROUTINE SYMFINT(IPRINT,IREPIJ,IREPIJB,IREPKLB,IJ,IKR,JKR,
     &                   NSTR3,NSTR4,NKL,FMAT,WORK,KFREE,LFREE)
C
C     Symmetrize the final integrals and write them to file
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dcbham.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      INTEGER NSTR3(2),NSTR4(2)
      DIMENSION FMAT(NKL,NZ,NZ,*)
      DIMENSION IREPIJB(2),IREPKLB(2,NKL)
C
      CALL MEMGET('REAL',KFTMAT,NKL*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDK,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDL,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRKL,2*NKL,WORK,KFREE,LFREE)
C
      IF (SPINFR.OR.LEVYLE) THEN
C
C     Special code that takes out the quaternion phase factors
C     The resulting integrals have full 8-fold permutational
C     symmetry and can be fed into non-relativistic codes.
C     Only one class of integrals (all-unbarred) suffices for
C     this purpose.
C
      CALL DZERO(WORK(KFTMAT),2*NKL)
      IREPABR = 0
      DO IREPAB = 0, NBSYM-1
         IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1
C
C           Zero out the spin-orbit coupling for electron 1
C
            IF (IREPAB.NE.IREPIJB(1)) GOTO 100
            DO I = 1, NKL
               IF (IREPAB.NE.IREPKLB(1,I)) THEN
C
C                  Zero out the spin-orbit coupling for electron 2
C
                   FMAT(I,1,1,IREPABR) = 0.D0
               ELSE
C
C                  Take out quaternion phase factor by multiplying
C                  with the compound phase factor that was introduced
C                  in the consecutive multiplications
C 
                   FMAT(I,1,1,IREPABR) = 
     &             FMAT(I,1,1,IREPABR) *
     &             IREPIJB(2)*IREPKLB(2,I)
               ENDIF
            ENDDO
            FACTR = 1.D0
            CALL DAXPY (NKL,FACTR,
     &                  FMAT(1,1,1,IREPABR),1,
     &                  WORK(KFTMAT),1)
 100        CONTINUE
         ENDIF
      ENDDO
C
C     Only one class remains in the spinfree case
C
      ICLASS = 1
      CALL MDINTB(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &            WORK(KFTMAT),NKL,
     &            WORK(KINDK),WORK(KINDL),
     &            WORK(KRKL))
C
      ELSE ! spinfree
C
C     Real and complex classes
C     They are always written as complex for MOLFDIR
C
      DO ICLASS = 1, 4
         CALL DZERO(WORK(KFTMAT),2*NKL)
         IREPABR = 0
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IREPABR = IREPABR + 1
               DO IZ1 = 1, NZ
                  IQ1 = IPQTOQ(IZ1,IREPAB)
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
                     IFAC = INTCL(ICLASS,IQ1,IQ2)
                     IF (IFAC.NE.0) THEN
                        FACTR = IFAC
                        CALL DAXPY (NKL,FACTR,
     &                              FMAT(1,IZ2,IZ1,IREPABR),1,
     &                              WORK(KFTMAT),1)
                     ENDIF
                     IF (NZ.GT.1) THEN
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                                 FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                 WORK(KFTMAT+NKL),1)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         CALL MDINTB(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &                           WORK(KFTMAT),NKL,
     &                           WORK(KINDK),WORK(KINDL),
     &                           WORK(KRKL))
      ENDDO
C
C     Quaternion classes
C
      IF (NZ.EQ.4) THEN
         DO ICLASS = 9, 12
            CALL DZERO(WORK(KFTMAT),2*NKL)
            IREPABR = 0
            DO IREPAB = 0, NBSYM-1
               IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ1 = 1, NZ
                     IQ1 = IPQTOQ(IZ1,IREPAB)
                     DO IZ2 = 1, NZ
                        IQ2 = IPQTOQ(IZ2,IREPAB)
                        IFAC = INTCL(ICLASS,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                            FACTR = IFAC
                            CALL DAXPY (NKL,FACTR,
     &                                  FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                  WORK(KFTMAT),1)
                        ENDIF
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                                 FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                 WORK(KFTMAT+NKL),1)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            CALL MDINTB(IPRINT,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &                  WORK(KFTMAT),NKL,
     &                  WORK(KINDK),WORK(KINDL),
     &                  WORK(KRKL))
         ENDDO
      ENDIF
C
      ENDIF ! spinfree
C
      CALL MEMREL('SYMFINT',WORK,1,KFTMAT,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MKINDXKR*/
      SUBROUTINE MKINDXKR (NSTR,TRIAN,INDXKR)
C
C     Make an array that gives the Kramers pairs from the ij-value
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      LOGICAL TRIAN
      INTEGER NSTR(2,0:2,2)
      DIMENSION INDXKR(2,*)
C
      IJ = 0
      DO IREPIJ =1, NFSYM
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = 1, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    INDXKR(1,IJ) = I
                    INDXKR(2,IJ) = J
c                   write(6,*) '!trian, IJ = ',IJ,I,J
                    IF (IREPI.NE.1) INDXKR(1,IJ) = I + NSTR(IREPI-1,0,1)
                    IF (IREPJ.NE.1) INDXKR(2,IJ) = J + NSTR(IREPJ-1,0,1)
                 ENDDO
              ENDDO
            ELSEIF (IREPI.EQ.IREPJ) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = J, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    INDXKR(1,IJ) = I
                    INDXKR(2,IJ) = J
c                   write(6,*) 'trian, IJ = ',IJ,I,J
                    IF (IREPI.NE.1) INDXKR(1,IJ) = I + NSTR(IREPI-1,0,1)
                    IF (IREPJ.NE.1) INDXKR(2,IJ) = J + NSTR(IREPJ-1,0,1)
                 ENDDO
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRI */
      SUBROUTINE MDSCRI(WORK,KFREE,LFREE,NODE)
C
C     Initialize MDINT scratch files
C     Written by Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      CHARACTER*10 FNODE,LFNAM,FNAM*5
C
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open direct access scratch file for 4-indexed transformed integrals
C
      DO IREPIJ = 1, NFSYM
         IRECL = IRECLEN(NQQCLASS*NFPCK34(IREPIJ),0,0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUMLF3+IREPIJ,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
         CALL MEMGET('REAL',KFTMAT,NKL*NQQCLASS,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFTMAT),NKL*NQQCLASS)
         DO IJ = 1, NIJ
            WRITE (LUMLF3+IREPIJ,REC=IJ) 
     &            (WORK(IX),IX=KFTMAT,KFTMAT+NKL*NQQCLASS-1)
         ENDDO
         CALL MEMREL('MDSCRI',WORK,1,KFTMAT,KFREE,LFREE)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRR */
      SUBROUTINE MDSCRR(IPRINT,IREPIJ,IJ,NKL,FMAT,WORK,KFREE,LFREE)
C
C     Read the final integrals from the scratch file
C
C     Written by Luuk Visscher, July 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      DIMENSION FMAT(*)
C
C     Read integrals.
C
      READ (LUMLF3+IREPIJ,REC=IJ) (FMAT(I),I=1,NKL*NQQCLASS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRU */
      SUBROUTINE MDSCRU(IPRINT,IREPIJ,IJ,NKL,FMAT,WORK,KFREE,LFREE)
C
C     Update the scratch file with final integrals
C
C     Written by Luuk Visscher, July 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      DIMENSION FMAT(*)
C
C     Read integrals from previous batch and add up.
C
      CALL MEMGET('REAL',KFTMAT,NKL*NQQCLASS,WORK,KFREE,LFREE)
      READ (LUMLF3+IREPIJ,REC=IJ) (WORK(IX),IX=KFTMAT,
     &                              KFTMAT+NKL*NQQCLASS-1)
      FACTR = 1
      CALL DAXPY (NQQCLASS*NKL,FACTR,WORK(KFTMAT),1,FMAT,1)
      WRITE (LUMLF3+IREPIJ,REC=IJ) (FMAT(I),I=1,NKL*NQQCLASS)
      CALL MEMREL('MDSCRU',WORK,1,KFTMAT,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRW */
      SUBROUTINE MDSCRW(IPRINT,IREPIJ,IJ,NKL,FMAT,WORK,KFREE,LFREE)
C
C     Write the final integrals to scratch file.
C
C     Written by Luuk Visscher, July 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      DIMENSION FMAT(*)
C
      WRITE (LUMLF3+IREPIJ,REC=IJ) (FMAT(I),I=1,NKL*NQQCLASS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRF */
      SUBROUTINE MDSCRF
C
C     Delete MDSCR file
C     Written by Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
C
      DO IREPIJ = 1, NFSYM
         IF (.NOT. NOMDCINT) THEN
C           delete the files
            CLOSE (LUMLF3+IREPIJ,STATUS='DELETE')
         ELSE
C           close the files
            CLOSE (LUMLF3+IREPIJ,STATUS='KEEP')
         END IF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MKINDXB*/
      SUBROUTINE MKINDXB (NSTR,IBE11,IBE21,IBE12,IBE22,TRIAN,INDXB)
C
C    Make an array that gives the compound boson irrep from the ij-value
C    It also gives the sign of the quaternion phase factor that this
C    density obtains in the multiplication of the MO coefficients. This
C    sign is required when we want to backtransform (in spinfree calcs)
C    the integrals such that only class 1 integrals with full 8-fold
C    permutation symmetry remain.
C
C    Written by Luuk Visscher, April 1997
C    Last modification, Februari 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
#include <dcbqal.h>
      LOGICAL TRIAN
      INTEGER NSTR(2,0:2,2)
      DIMENSION IBE11(*),IBE21(*),IBE12(*),IBE22(*)
      DIMENSION INDXB(2,*)
#include <ibtfun.h>
C
      IJ = 0
      DO IREPIJ =1, NFSYM
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = 1, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    IF (IREPI.EQ.1) IREPIB = IBE11(I)
                    IF (IREPI.EQ.2) IREPIB = IBE21(I)
                    IF (IREPJ.EQ.1) IREPJB = IBE12(J)
                    IF (IREPJ.EQ.2) IREPJB = IBE22(J)
C
C                   The compound boson irrep can be obtained
C                   by an xor operation.
C
                    IREPIJB = IBTXOR(IREPIB,IREPJB)
                    INDXB(1,IJ) = IREPIJB
C
C                   We first determine which quaternion phase
C                   factors the individual functions receive and
C                   than obtain the compound phase factor by
C                   quaternion multiplication.
C
                    IQI = IPQTOQ(1,IREPIB)
                    IQJ = IPQTOQ(1,IREPJB)
C
C                   The sign of the quaternion multiplication is given by iqphase,
C                   but since IQI is conjugated we also have to 
C                   multiply by IQSIGN (-1 for q-imaginary, 1 for real)
C
                    INDXB(2,IJ) = IQPHASE(IQI,IQJ,1)*IQSIGN(IQI,2,1)
                 ENDDO
              ENDDO
            ELSEIF (IREPI.EQ.IREPJ) THEN
C
C             Same procedure for lower trangular case..
C
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = J, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    IF (IREPI.EQ.1) IREPIB = IBE11(I)
                    IF (IREPI.EQ.2) IREPIB = IBE21(I)
                    IF (IREPJ.EQ.1) IREPJB = IBE12(J)
                    IF (IREPJ.EQ.2) IREPJB = IBE22(J)
                    IREPIJB = IBTXOR(IREPIB,IREPJB)
                    INDXB(1,IJ) = IREPIJB
                    IQI = IPQTOQ(1,IREPIB)
                    IQJ = IPQTOQ(1,IREPJB)
                    INDXB(2,IJ) = IQPHASE(IQI,IQJ,1)*IQSIGN(IQI,2,1)
                 ENDDO
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INIMDINFO */
      SUBROUTINE INIMDINFO
C
C     Initialize MDINFO file
C     Authored by Joern Thyssen, 1/11/2000.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <dcbtra.h>
#include <dgroup.h>
C
      CHARACTER*8 DATEX,TIMEX
C
C     Open file
C
      OPEN (LUMLF2,FILE='MDINFO',FORM='UNFORMATTED')
C
C     Write date and time
C
      CALL GETDAT (DATEX,TIMEX)
      WRITE (LUMLF2) DATEX,TIMEX
C
C     The rest of the file is written in various places. The contents is:
C
C     Number of IJ indices for each fermion symmetry
C     Number of KL indices for each fermion symmetry
C     I and J indeces
C     Boson irrep for IJ indeces
C     Boson irrep for KL indeces
C
C
      RETURN
      END
