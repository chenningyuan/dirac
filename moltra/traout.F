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
      SUBROUTINE SYMFINT(IPRINT,IREPIJ,IJ,IKR,JKR,NSTR3,NSTR4,NKL,FMAT,
     &                   IQMO,WORK,KFREE,LFREE)
C
C     Symmetrize the final integrals and write them to file
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION WORK(*)
      INTEGER NSTR3(2),NSTR4(2)
      DIMENSION IQMO(4,0:7)
      DIMENSION FMAT(NKL,NZ,NZ,*)
C
      CALL MEMGET('REAL',KFTMAT,NKL*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDK,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDL,NKL,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRKL,2*NKL,WORK,KFREE,LFREE)
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
                  IQ1 = IQMO(IZ1,IREPAB)
C
C                 Fix for the complex groups. We need to redefine the 
C                 spin-quantization axis. This should work...
C
                  IF (NZ.EQ.2) IQ1 = IQCMPL(IQ1)
C
                  DO IZ2 = 1, NZ
                     IQ2 = IQMO(IZ2,IREPAB)
                     IF (NZ.EQ.2) IQ2 = IQCMPL(IQ2)
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
            DO IREPAB = 0, NBSYM-1
               IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
                  DO IZ1 = 1, NZ
                     IQ1 = IQMO(IZ1,IREPAB)
                     DO IZ2 = 1, NZ
                        IQ2 = IQMO(IZ2,IREPAB)
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
      SUBROUTINE MDSCRI(WORK,KFREE,LFREE,NFPCK12,NFPCK34)
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
      DIMENSION NFPCK12(2),NFPCK34(2)
      CHARACTER*10 FNODE,LFNAM,FNAM*5
C
C     Open direct access scratch file for 4-indexed transformed integrals
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
      DO IREPIJ = 1, NFSYM
         IRECL = IRECLEN(NCLASS*NFPCK34(IREPIJ),0,0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         call quit('the statement below cannot be executed')
c        FNODE = LFNAM(FNAM)
         OPEN (LUMLF3+IREPIJ,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
         CALL MEMGET('REAL',KFTMAT,NKL*NCLASS,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KFTMAT),NKL*NCLASS)
         DO IJ = 1, NIJ
            WRITE (LUMLF3+IREPIJ,REC=IJ) 
     &            (WORK(IX),IX=KFTMAT,KFTMAT+NKL*NCLASS-1)
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
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Read integrals.
C
      READ (LUMLF3+IREPIJ,REC=IJ) (FMAT(I),I=1,NKL*NCLASS)
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
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Read integrals from previous batch and add up.
C
      CALL MEMGET('REAL',KFTMAT,NKL*NCLASS,WORK,KFREE,LFREE)
      READ (LUMLF3+IREPIJ,REC=IJ) (WORK(IX),IX=KFTMAT,
     &                              KFTMAT+NKL*NCLASS-1)
      FACTR = 1
      CALL DAXPY (NCLASS*NKL,FACTR,WORK(KFTMAT),1,FMAT,1)
      WRITE (LUMLF3+IREPIJ,REC=IJ) (FMAT(I),I=1,NKL*NCLASS)
      CALL MEMREL('MDSCRU',WORK,1,KFTMAT,KFREE,LFREE)
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
         CLOSE (LUMLF3+IREPIJ,STATUS='DELETE')
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck IRECLEN*/
      FUNCTION IRECLEN (NREAL,NINTE,NLOGI)
C
C     Returns the record length for a record consisting of the specified
C     number of variables.
C
C     Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      INTEGER NREAL,NINTE,NLOGI,IRECLEN
C
#if defined (SYS_VAX)
      IRECL = 2 * NREAL + 1 * NINTE + 1 * NLOGI
#elif defined (SYS_CRAY) || defined (SYS_T3D)
      IRECL = 8 * NREAL + 8 * NINTE + 8 * NLOGI
#elif defined (SYS_IRIX)
C     SGI IRIX counts recl in *4 words
      IRECL = 2 * NREAL + 1 * NINTE + 1 * NLOGI
#elif defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_HPUX) || defined (SYS_LINUX)
      IRECL = 8 * NREAL + 4 * NINTE + 4 * NLOGI
#elif defined (SYS_IBM)
      IRECL = 8 * NREAL + 4 * NINTE + 4 * NLOGI
#else
     ===== insert definition for specification of record length
#endif
C
      IRECLEN = IRECL
C
      RETURN
      END
