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
C FILE    : tradr2.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR2T */
      SUBROUTINE TRDR2T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &                  NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2)
C
C     Written by Luuk Visscher, Februari 1997
C     Modifications: Jon K. Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),INDX(3,*),NSTRT(4)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDR2T')
      KFRSAV = KFREE
C
C     Check that no negative energy solutions are active.
C     Active negative energy solutions is not implemented
C     and is NOT a good idea unless one knows exactly
C     what one is doing...
C                                     Jon
C
      DO IND = 1,4
         DO IF = 1,NFSYM
            IF (NSTR(IF,2,IND).NE.0) 
     &           CALL QUIT
     &           ('TRDR2T: Negative energy solutions in active space!')
         ENDDO
      ENDDO
C
      DO I = 1,4
         NSTRT(I) = 0
         DO IFRP = 1, NFSYM
            NSTRT(I) = NSTRT(I) + NSTR(IFRP,0,I)
         ENDDO
      ENDDO
      IF (NSTRT(1).NE.NSTRT(3)) CALL QUIT
     &     ('TRD2T: Different number of active for index 1 and 3!')
      IF (NSTRT(2).NE.NSTRT(4)) CALL QUIT
     &     ('TRD2T: Different number of active for index 2 and 4!')
C     
C     Divide index transformation in tasks with different batches
C     of spinors. First find number of tasks.
C
      CALL TSKCHK(MAXTSK,NSTR,LFREE,LENTSK,NSTRT,IPRINT)
      CALL MEMGET('INTE',K1,MAXTSK*2*3*4,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',K2,MAXTSK*2*4,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',K3,MAXTSK*4,WORK,KFREE,LFREE)
C
      CALL TRDR22(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &     NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,
     &     WORK(K1),WORK(K2),MAXTSK,LENTSK,NSTRT,WORK(K3))
C     
      CALL MEMREL('TRDR2T',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDR2T')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR22 */
      SUBROUTINE TRDR22(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &     NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,
     &     NSTR1,NSTRS,MAXTSK,LENTSK,NSTRT,NOCC)
C
C     Written by Luuk Visscher, Februari 1997
C     and Jon K Laerdahl, June 1997 (A VERY warm day in Oslo..)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbmp2.h>
      LOGICAL ANTIS,LMP2
      LOGICAL IBEQJB
      DIMENSION NSTR(2,0:2,4),NSTR1(2,0:2,4,MAXTSK),NSTRS(2,4,MAXTSK)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),NSTRT(4),INDX(3,*)
      DIMENSION ICMOQR1(2,4),NOCC(MAXTSK,4)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION WORK(*),TIM1(7)
C
      CALL QENTER('TRDR22')
C     
      CALL IZERO(NSTR1,MAXTSK*2*3*4)
      CALL IZERO(NSTRS,MAXTSK*2*4)
      DO 100 IJTSK = 1,MAXTSK
         CALL UNPKIJ(IJTSK,ISPNR,JSPNR)
         NOCC(IJTSK,1) = (ISPNR-1)*LENTSK + 1
         NOCC(IJTSK,3) = (JSPNR-1)*LENTSK + 1
         NOCC(IJTSK,2) = ISPNR*LENTSK
         NOCC(IJTSK,4) = JSPNR*LENTSK
C     Last Batch of Kramers pairs?         
         IF (NOCC(IJTSK,2).GT.NSTRT(1)) NOCC(IJTSK,2)=NSTRT(1)
         IF (NOCC(IJTSK,4).GT.NSTRT(3)) NOCC(IJTSK,4)=NSTRT(3)
         IF(IPRINT.GE.1) THEN
            WRITE(LUPRI,'(A,3I5)') 'Task-info 1:',IJTSK,ISPNR,JSPNR
            WRITE(LUPRI,'(12X,4I5)') (NOCC(IJTSK,I),I=1,4)
         ENDIF
         IF(IPRINT.GE.1) WRITE(LUPRI,'(A)') 'Task-info 2:'
         DO IN = 0,1
            IND1=IN*2+1
            IND2=IN*2+2
            LSTIR1 = NSTR(1,0,IND1)
C     
C     Whole batch in irrep 1
C     
            IF((NOCC(IJTSK,IND1).LE.LSTIR1).AND.
     &           (NOCC(IJTSK,IND2).LE.LSTIR1)) THEN
               NSTRS(1,IND1,IJTSK)=NOCC(IJTSK,IND1)
               NSTR1(1,0,IND1,IJTSK)=
     &              NOCC(IJTSK,IND2)-NOCC(IJTSK,IND1)+1
               NSTRS(2,IND1,IJTSK)=-1
               NSTR1(2,0,IND1,IJTSK)=0
C     
C     Batch contributions from both irreps
C     
            ELSEIF((NOCC(IJTSK,IND1).LE.LSTIR1).AND.
     &              (NOCC(IJTSK,IND2).GT.LSTIR1)) THEN
               NSTRS(1,IND1,IJTSK)=NOCC(IJTSK,IND1)
               NSTR1(1,0,IND1,IJTSK)=LSTIR1-NOCC(IJTSK,IND1)+1
               NSTRS(2,IND1,IJTSK)=1
               NSTR1(2,0,IND1,IJTSK)=NOCC(IJTSK,IND2)-LSTIR1
C     
C     Whole batch in irrep 2
C         
            ELSE
               NSTRS(1,IND1,IJTSK)=-1
               NSTR1(1,0,IND1,IJTSK)=0
               NSTRS(2,IND1,IJTSK)=NOCC(IJTSK,IND1)-LSTIR1
               NSTR1(2,0,IND1,IJTSK)=
     &              NOCC(IJTSK,IND2)-NOCC(IJTSK,IND1)+1
            ENDIF
         ENDDO
C     
         CALL ICOPY (6,NSTR(1,0,2),1,NSTR1(1,0,2,IJTSK),1)
         CALL ICOPY (6,NSTR(1,0,4),1,NSTR1(1,0,4,IJTSK),1)
         NSTRS(1,2,IJTSK) = 1
         NSTRS(1,4,IJTSK) = 1
         NSTRS(2,2,IJTSK) = 1
         NSTRS(2,4,IJTSK) = 1
         DO IN = 0,1
            IND = IN*2+1
            DO IF = 1,NFSYM
               NSTR1(IF,1,IND,IJTSK) = NSTR1(IF,0,IND,IJTSK)
               NSTR1(IF,2,IND,IJTSK) = 0
               IF(IPRINT.GE.1) 
     &              WRITE(LUPRI,'(7I8)') IF,IND,IJTSK, 
     &              NSTRS(IF,IND,IJTSK),NSTR1(IF,0,IND,IJTSK),
     &              NSTR1(IF,1,IND,IJTSK),NSTR1(IF,2,IND,IJTSK)
            ENDDO
         ENDDO
 100  CONTINUE
C     
C
C     Zero MP2 energy and start looping over tasks
C
      IF (LMP2) THEN
         WRITE(LUPRI,'(/,15X,A,/)') '**** Calculating MP2 energy. ****'
         EMP2 = 0.0D00
         CALL DZERO(TIM1,7)
      ENDIF
      DO ITASK = 1, MAXTSK
C
         WRITE(LUPRI,'(A,I3,A,I3)') 'Doing IJ-batch ',ITASK,
     &        ' of ',MAXTSK

C     
C     Give relative offset for the coefficients
C     
         DO INDEX = 1, 4
            DO I = 1, NFSYM
               ICMOQR1(I,INDEX) = ICMOQR(I,INDEX)
     &              + (NSTRS(I,INDEX,ITASK)-1)*NFBAS(I,0)
            ENDDO
         ENDDO
         IOFF1 = NOCC(ITASK,1)
         IOFF2 = NOCC(ITASK,3)
C     
C     I and J same range -> IBEQJB = .TRUE.
C     
         IBEQJB = .FALSE.
         IBEQJB = IOFF1.EQ.IOFF2
         CALL TRDR23(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR1,
     &        NSTR1(1,0,1,ITASK),ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,
     &        E1(IOFF1),E2,E1(IOFF2),E2,IBEQJB,TIM1)
C
      ENDDO
      IF (LMP2) THEN
         WRITE(LUPRI,'(A)') ''
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A)') 
     &        'Final Relativistic Direct MP2 results:'
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(/15X,A,F20.10)') '@ MP2 energy       :', EMP2
         WRITE(LUPRI,'(15X,A,F20.10)')  '  DHF energy       :', DHFEMP
         WRITE(LUPRI,'(15X,A,F20.10,/)') '@ Total MP2 energy :', 
     &        DHFEMP+EMP2
         CALL TIMTXT('Time used in MP2CAL is: ',TIM1(5),LUPRI)
      ENDIF
      CALL TIMTXT('Time used in CALDIS is: ',TIM1(1),LUPRI)
      CALL TIMTXT('Time used in MS4IN1 is: ',TIM1(2),LUPRI)
      CALL TIMTXT('Time used in TQTRNS is: ',TIM1(3),LUPRI)
      CALL TIMTXT('Time used in LQTRNS is: ',TIM1(4),LUPRI)
C
      CALL QEXIT('TRDR22')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR23 */
      SUBROUTINE TRDR23(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &     NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,E3,E4,IBEQJB,TIM1)
C
C     Written by Luuk Visscher, December 1996
C     Modifications by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxorb.h>
#include <dcbgen.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbmpt.h>
#include <dgroup.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
      LOGICAL ANTIS, LMP2, NODV, NOPV,L2TYP(0:3),IBEQJB
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION NSTRT(4), ITQMT(0:7,0:7,2)
      DIMENSION IQMO(4,0:7),NZMO(0:7)
      DIMENSION IQP34(4,0:7),NZ34(0:7)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION IJP12(0:7)
      DIMENSION IJP34(2)
      DIMENSION NFPCK12(2),IFPCK12(2,2)
      DIMENSION NFPCK34(2),IFPCK34(2,2)
      DIMENSION IHM(0:7)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*),E3(*),E4(*)
      DIMENSION WORK(*),TIM1(7)
C
#include <ibtfun.h>
C
C     Notice: Range for index 1, 2, 3 and 4 are not
C     the same. E1(3) is NOT necessarily = E3(3)
C     and in MP2 E1(3).NE.E2(3).
C
C     To implement MSOUT: Get absolute position of
C     active occupied Kramers pairs from outside.
C     Apart from this Luuk has already the output
C     routine for MOLFIR integrals ready. Use this..
C
      CALL QENTER('TRDR23')
      KFRSAV = KFREE
C
C     Define the packing of the 2-index transformed integrals
C     We do not consider sorting of positronic and electronic
C     spinors at this point.
C
      DO IREPIJ =1, NFSYM
         NFPCK12(IREPIJ) = 0
         NFPCK34(IREPIJ) = 0
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
            IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
            NFPCK12(IREPIJ) = NFPCK12(IREPIJ) +
     &                        NSTR(IREPI,0,1)*NSTR(IREPJ,0,2)
            NFPCK34(IREPIJ) = NFPCK34(IREPIJ) +
     &                        NSTR(IREPI,0,3)*NSTR(IREPJ,0,4)
         ENDDO
      ENDDO
C
      IF(IPRINT.GE.1) THEN
         WRITE(LUPRI,'(A,2I8)') 'NFPCK12',(NFPCK12(I), I=1,NFSYM)
         WRITE(LUPRI,'(A,2I8)') 'NFPCK34',(NFPCK34(I), I=1,NFSYM)
         WRITE(LUPRI,'(A,4I8)') 'IFPCK12',
     &        ((IFPCK12(I,J), I=1,NFSYM), J=1,NFSYM)
         WRITE(LUPRI,'(A,4I8)') 'IFPCK34',
     &        ((IFPCK34(I,J), I=1,NFSYM), J=1,NFSYM)
      ENDIF
C
C     Make index array for 3/4-transformed integrals
C
      NTQMT = 0
      CALL IZERO(ITQMT,8*2*8)
      DO IREPRS = 0, NBSYM-1
         IREPIJ = JBTOF(IREPRS,1) ! Get parity of rs 
         NIJ = NFPCK12(IREPIJ)  ! Number of ij pairs 
         DO IREPS = 0, NBSYM-1
            DO 11 IC = 1,2
               IREPR = IBTXOR(IREPRS,IREPS)
               NR = NBBAS(IREPR,IC) ! Number of r functs
               ITQMT(IREPRS,IREPS,IC) = NTQMT ! Pointers
               IREPL = JBTOF(IREPS,IC) ! Parity that s contributes to
               NS = NBBAS(IREPS,IC) ! Number of s functs
               IF (IREPIJ.EQ.1) THEN
                  IREPK=IREPL
               ELSE
                  IREPK=3-IREPL
               ENDIF
c              WRITE(LUPRI,'(A,I10,10I6)') 
c    &              'NTQMT',NTQMT,NZ,NS,NIJ,NSTR(IREPK,0,3),IREPK,
c    &              IREPIJ,IREPS,IREPRS,NR,IREPR
               IF (NR.EQ.0) GOTO 11
               NTQMT = NTQMT + NIJ*NSTR(IREPK,0,3)*NS*NZ*NZ
 11         CONTINUE
         ENDDO
      ENDDO
C
C
C     Allocate memory for 3/4 Transformed integrals
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1000) NTQMT*8./(1024.*1024.)
C
      IF (LFREE.LT.NTQMT) THEN
            WRITE(LUPRI,'(/A,/A,/A)') 
     &           '### Not enough memory for 3/4 transf. matrix.   ###',
     &           '### Reduce memory requirements by setting       ###',
     &           '### .IJTSK to a lower value!                    ###'
C
      ENDIF
      CALL MEMGET('REAL',KTQMAT,NTQMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTQMAT),NTQMT)
C
C     Loop over (XX|LL) and (XX|SS) integrals
C
      IF (TESTLS) WRITE(LUPRI,'(/A)') 'YES: Doing (LL|SS) ONLY.'
      IF (TESTSL) WRITE(LUPRI,'(/A)') 'YES: Doing (SS|LL) ONLY.'
C
c     print*,'LFREE',LFREE
C
C     ----------------
      DO 100 IC = 1, 2 
C     ----------------
C
         IF ((IC.EQ.1).AND.TESTLS) GOTO 100
         IF ((IC.EQ.2).AND.TESTSL) GOTO 100
         CALL SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
         IF (I2TYP.LT.0) THEN
            GOTO 100
         ELSEIF (I2TYP.EQ.1) THEN
            ICS = 1
            ICF = 1
         ELSEIF (I2TYP.EQ.2) THEN
            ICS = 2
            ICF = 2
         ELSE
            ICS = 1
            ICF = 2
         ENDIF
C
C        Loop over the shells
C
         DO ISHLA = IASTRT, IASMAX
         DO ISHLB = IBSTRT, ISHLA
C
            KFRSAV3 = KFREE
            NINSHA = NINSH(ISHLA,-1,INDX,IPRINT)
            NINSHB = NINSH(ISHLB,-1,INDX,IPRINT)
C
C           Set up gather array for distributions, 
C           calculate memory requirements for this batch.
C
C           NOTE : The actual dimension of INDXAB is kept in the
C           common block dcbtra.h. When we want to parallelize it may
C           be better to allocate it always as MXINSH*MXINSH*4,
C           or (even better) to send it also down to HERMIT.
C           Luuk.
C
C           First a dummy run to calculate the memory requirements
C           for GMAT. Max size of GMAT is set in MAXSCL from
C           .SCLMEM in input or chosen as 1/3 of LFREE.
C     
            IF(MAXSCL.EQ.0) THEN
               NSIZG = LFREE/3
            ELSE
               NSIZG = MAXSCL
               
               IF (MAXSCL.GT.LFREE/100*80) WRITE(LUPRI,'(2(1X,A,I10))')
     &              '.GT. 80%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
               IF (MAXSCL.GT.LFREE/100*95) THEN
                  WRITE(LUPRI,'(2(1X,A,I10))')
     &                 '.GT. 95%.., LFREE:',LFREE,'MAXSCL:',MAXSCL
                  CALL QUIT('GMAT too large...')
               ENDIF
            ENDIF
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &           .TRUE.,INDX,IDUM,IJP12,NSIZG,IDUM)
C
            LDXAB = 5*NINSHA*NINSHB
            CALL MEMGET('INTE',KDXAB,LDXAB,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
C           Set INDXAB and INDPASS
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.TRUE.,
     &                       INDX,WORK(KDXAB),IJP12,NSIZG,WORK(KIJPASS))
C
            NSIZH = 0
            DO IREPAB = 0, NBSYM-1
               IHM(IREPAB) = NSIZH
               IREPIJ = JBTOF(IREPAB,1)
               NSIZH = NSIZH + IJP12(IREPAB)*NFPCK12(IREPIJ)
            ENDDO
C
            CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
            CALL DZERO (WORK(KHMAT),NSIZH*NZ)
C
            IF (IPRINT .GE. 3) THEN
               WRITE (LUPRI,
     &              '(//,1X,A,2I8/,2(1X,A,I8/),2(1X,A,I8,F10.3/))')
     &         ' Calculation of integral distribution set:',
     &              ISHLA,ISHLB,
     &         ' Number of distributions in this set:     ',NUMDIS,
     &         ' Number of passes :                       ',NPASS,
     &         ' Size of scalar integral batch:           ',NSIZG,
     &              NSIZG*8./(1024.*1024.),
     &         ' Size of half-transformed integral batch: ',NSIZH,
     &              NSIZH*8./(1024.*1024.)
            END IF
C
            IF (NPASS.GT.1) THEN
C
C              Initialize scalar integral buffer file and arrays
C              Use the space reserved for GMAT now as buffer storage
C
               LGFIL = LUTRA2
               NGBFSZ = (IRAT*NSIZG)/((1+IRAT)*NPASS)
C
C              Allocate the buffer arrays in WORK and set KGMAT to the
C              start of this allocation section. We now only need to
C              send down GMAT into HERMIT to find this section.
C
               CALL MEMGET('INTE',KALGREC,NPASS,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KALGBUF,NPASS,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KAIGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
               CALL MEMGET('REAL',KARGBUF,NPASS*NGBFSZ,WORK,KFREE,LFREE)
C
C              As we refer to the places inside the section we need
C              to subtract KGMAT.
C
               KGMAT  = KALGREC
               KLGREC = KALGREC - KGMAT + 1
               KLGBUF = KALGBUF - KGMAT + 1
               KIGBUF = KAIGBUF - KGMAT + 1
               KRGBUF = KARGBUF - KGMAT + 1
            ELSE
C
C              Write directly into the array GMAT
C
               CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
               CALL DZERO(WORK(KGMAT),NSIZG)
               KLGBUF = 1
            ENDIF
C
C           Get distributions
C
            NODV = NASHT.EQ.0
            NOPV = NASHT.LT.2
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &                  NODV,NOPV,WORK(KFREE),LFREE,IPRINT)
C
            CALL TIMEC('CALDIS',TIMSTR,TIMEND) 
            TIM1(1) = TIM1(1) + TIMEND - TIMSTR
C
C           Do the transformation of the first index pair and sort the
C           half-transformed integrals to memory
C
C           De-allocate buffers and allocate GMAT if necessary.
C
            IF (NPASS.GT.1) THEN
               CALL MEMREL('TRDRV3.1',WORK,1,KALGBUF,KFREE,LFREE)
               CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
               CALL DZERO(WORK(KGMAT),NSIZG)
            ENDIF
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL MS4IN1E(WORK,KFREE,LFREE,IPRINT,WORK(KALGREC),ICS,ICF,
     &                   NSTR(1,0,1),NSTR(1,0,2),
     &                   NFPCK12,IFPCK12,NDMOQR,ICMOQR,WORK(KIJPASS),
     &                   Q1,Q2,WORK(KGMAT),NSIZG,NZMO,IQMO,
     &                   WORK(KHMAT)) 
            CALL TIMEC('MS4IN1',TIMSTR,TIMEND) 
            TIM1(2) = TIM1(2) + TIMEND - TIMSTR
c           write (LUPRI,'(//A)') ' The 2/4 transformed integrals:'
c           do ix = 1, NSIZH
c              print*,ix,work(KHMAT+ix-1)
c           enddo
C
            CALL MEMREL('TRDRV3.2',WORK,1,KGMAT,KFREE,LFREE)
C
C           Delete scalar integral buffer file
C
            IF (NPASS.GT.1) THEN
               CALL DELGBUF(LGFIL)
            ENDIF
C
C
C     3/4 transformation
C
            CALL TIMEC('START ',TIMSTR,TIMEND)
            CALL TQTRNS (WORK,KFREE,LFREE,IPRINT,IC,IHM,
     &           INDX,WORK(KDXAB),NFPCK12,WORK(KHMAT),ITQMT,
     &           ICMOQR(1,3),NDMOQR(1,1,3),WORK(KTQMAT),NSTR(1,0,3),
     &           IQP34,NZ34,Q3)
C
            CALL TIMEC('TQTRNS',TIMSTR,TIMEND) 
            TIM1(3) = TIM1(3) + TIMEND - TIMSTR
            CALL MEMREL('TRDR23',WORK,1,KFRSAV3,KFREE,LFREE)
C
         ENDDO
         ENDDO
C
C----------------
 100     CONTINUE
C----------------
C
c      write (LUPRI,'(//A)') ' The 3/4 transformed integrals:'
c      do ix = 1, NTQMT
c         print*,ix,work(KTQMAT+ix-1)
c      enddo
C
      N4QMT = 0
      DO IREPIJ= 1,NFSYM
         DO IREPRS = 0, NBSYM-1
            IF(IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
               NIJ = NFPCK12(IREPIJ) ! Number of ij pairs 
               NKL = NFPCK34(IREPIJ) ! Number of kl pairs 
               N4QMT = N4QMT + NIJ*NKL*NZ**2
            ENDIF
         ENDDO
      ENDDO
C
C
C     Allocate memory for 4/4 Transformed integrals
C
      IF (IPRINT .GE. 1) WRITE (LUPRI,1001) N4QMT*8./(1024.*1024.)
      CALL MEMGET('REAL',K4QMAT,N4QMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(K4QMAT),N4QMT)
C     
C     4/4 transformation
C
      ICS2 = 1
      ICF2 = 2
      L2TYP(1) = MOD(INTFLG,2).EQ.1
      L2TYP(2) = MOD((INTFLG/2),2).EQ.1
      L2TYP(3) = MOD(((INTFLG/2)/2),2).EQ.1
      IF ((.NOT.L2TYP(1)).AND.(.NOT.L2TYP(2))) ICS2 = 2
      IF ((.NOT.L2TYP(2)).AND.(.NOT.L2TYP(3))) ICF2 = 1
      IF (TESTLS) ICS2 = 2
      IF (TESTSL) ICF2 = 1
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
      CALL LQTRNS (WORK,KFREE,LFREE,IPRINT,
     &     NFPCK12,NFPCK34,ITQMT,ICMOQR(1,4),NDMOQR(1,1,4),
     &     IQMO,NZMO,WORK(KTQMAT),NSTR(1,0,3),NSTR(1,0,4),
     &     IQP34,NZ34,Q4,WORK(K4QMAT),IFPCK34,ICS2,ICF2)
      CALL TIMEC('TQTRNS',TIMSTR,TIMEND) 
      TIM1(4) = TIM1(4) + TIMEND - TIMSTR
C     
c      write (LUPRI,'(//A)') ' The fully transformed integrals'
c      do ix = 1, n4qmt
c         print*,ix,work(k4qmat+ix-1)
c      enddo
C
      CALL TIMEC('START ',TIMSTR,TIMEND)
      DO I = 1,4
         NSTRT(I) = 0
         DO IFRP = 1, NFSYM
            NSTRT(I) = NSTRT(I) + NSTR(IFRP,0,I)
         ENDDO
      ENDDO
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(A)') 'Eigenvalues used in MP2:'
         CALL PRIEIN(E1,E2,E3,E4,
     &        NSTRT(1),NSTRT(2),NSTRT(3),NSTRT(4),NSTR)
      ENDIF
C
      IF (NZ.EQ.1) THEN
         NCLASS = 4
      ELSEIF (NZ.EQ.2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
C
      CALL DEFCL (NCLASS)
C
C     Allocate memory for naive antisymmetrize
C     -> Fully transformed integrals are
C     few compared to other arrays...
C     Storing final integrals in
C     TQMAT if there is enough space.
C
      LSIZ=16
      IF (NZ.LT.4) LSIZ=8
      DO I=1,4
         LSIZ=LSIZ*NSTRT(I)
      ENDDO
      IF(NTQMT.GE.LSIZ) THEN
         KANT = KTQMAT
         IF(IPRINT.GE.1)
     &        WRITE(LUPRI,'(A)') 
     &        '## Storing final integs. in 3/4-matrix. ##'
      ELSE
         IF(IPRINT.GE.1)
     &        WRITE(LUPRI,'(A)') 
     &        '## Final integs. .GT. 3/4-matrix. ##'
         CALL MEMGET('REAL',KANT,LSIZ,WORK,KFREE,LFREE)
      ENDIF
      IF (IPRINT .GE. 1) WRITE (LUPRI,1002) LSIZ*8./(1024.*1024.)
      CALL DZERO(WORK(KANT),LSIZ)
C
C     FQMAT stored as (ij,iz1,kl,iz2,Pk,Grs,Pkl) and the order
C     of ij (and kl) is gg,uu,ug,gu or as (i,j,iz1,k,l,iz2,Pk,Grs,Pkl)
C
C     F.ex. for Pkl=g, iz2 loops faster than loop over gg/uu for kl,
C     but slower than kl within gg/uu
C
      IOFF3 = 0
      DO IREPKL = 1,NFSYM
c       print*,' Fermion symmetry :',irepkl
         CALL MEMGET('REAL',KFTMAT,2*NFPCK12(IREPKL),WORK,KFREE,LFREE)
         LKR = 0
         DO IREPL = 1, NFSYM
            IREPK = MOD(IREPL+IREPKL,2) + 1
C     Offset for k and l symmetry (i.e. for Pkl=g, offset for uu)
            IOFF2 = IFPCK34(IREPK,IREPL)*NFPCK12(IREPKL)*NZ*NZ
            DO L = 1, NSTR(IREPL,0,4)
               LKR = LKR + 1
               DO K = 1, NSTR(IREPK,0,3)
                  IF (IREPK.EQ.1) THEN
                     KKR = K
                  ELSE
                     KKR = NSTR(1,0,3) + K
                  ENDIF
c     print*,' Batch labels (KKR,LKR):',kkr,lkr
C     
C     Real and complex classes
C     
C     Offset for K and L
                  IOFF1 = NFPCK12(IREPKL)*NZ*
     &                 (NSTR(IREPK,0,3)*(L-1)+(K-1))
                  DO JCLASS = 1, 4
                     CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                     ICNT1 = 0
                     DO IREPRS = 0, NBSYM-1
                        IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                           IOFF6 = ICNT1
                           ICNT1=ICNT1
     &                          +NFPCK12(IREPKL)*NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                           IRPCN1 = 0
                           DO IZ2 = 1, NZ
C     Offset for iz2 
                              IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                             NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C
                              IRPCN1 = IRPCN1 + 1
                              IQ2 = IQMO(IZ2,IREPRS)
                              IRPCN2 = 0
                              DO IZ1 = 1, NZ
C     Offest for iz1
                                 IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                 IRPCN2 = IRPCN2 + 1
                                 IQ1 = IQMO(IZ1,IREPRS)
                                 IFAC=INTCL(JCLASS,IQ1,IQ2)
                                 IOFF=IOFF1+IOFF2+IOFF3
                                 IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                 IF (IFAC.NE.0) THEN
                                    FACTR=IFAC
                                    CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                   WORK(K4QMAT+IOFF),1,
     &                                   WORK(KFTMAT),1)
                                 ENDIF
C 
C     Non-real groups
C
                                 IF (NZ.GT.1) THEN
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT+
     &                                      NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                    WORK(KFTMAT),WORK(KANT),
     &                    NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                    NSTRT(3),NSTRT(4),LSIZ)
                  ENDDO
C     
C     Quaternion groups
C     
                  IF(NZ.EQ.4) THEN
                     DO JCLASS = 9, 12
                        CALL DZERO(WORK(KFTMAT),2*NFPCK12(IREPKL))
                        ICNT1 = 0
                        DO IREPRS = 0, NBSYM-1
                           IF (IREPKL.EQ.JBTOF(IREPRS,1)) THEN
C     Counter for size of IREPKL=1 contribution and offset for
C     IREPKL
                              IOFF6 = ICNT1
                              ICNT1=ICNT1
     &                             +NFPCK12(IREPKL)
     &                             *NFPCK34(IREPKL)*NZ*NZ
C     Offset for Grs
                              IRPCN1 = 0
                              DO IZ2 = 1, NZ
C     Offset for iz2 
                                 IOFF4 = NFPCK12(IREPKL)*NZ*IRPCN1*
     &                                NSTR(IREPK,0,3)*NSTR(IREPL,0,4)
C     
                                 IRPCN1 = IRPCN1 + 1
                                 IQ2 = IQMO(IZ2,IREPRS)
                                 IRPCN2 = 0
                                 DO IZ1 = 1, NZ
C     Offest for iz1
                                    IOFF5 = NFPCK12(IREPKL)*IRPCN2
                                    IRPCN2 = IRPCN2 + 1
                                    IQ1 = IQMO(IZ1,IREPRS)
                                    IFAC=INTCL(JCLASS,IQ1,IQ2)
                                    IOFF=IOFF1+IOFF2+IOFF3
                                    IOFF=IOFF+IOFF4+IOFF5+IOFF6
C
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT),1)
                                    ENDIF
                                    IFAC=INTCL(JCLASS+4,IQ1,IQ2)
                                    IF (IFAC.NE.0) THEN
                                       FACTR=IFAC
                                       CALL DAXPY(NFPCK12(IREPKL),FACTR,
     &                                      WORK(K4QMAT+IOFF),1,
     &                                      WORK(KFTMAT
     &                                      +NFPCK12(IREPKL)),1)
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDIF
                        ENDDO
                        CALL MP2FI1(JCLASS,IREPKL,NSTR,KKR,LKR,
     &                       WORK(KFTMAT),WORK(KANT),
     &                       NFPCK12(IREPKL),NSTRT(1),NSTRT(2),
     &                       NSTRT(3),NSTRT(4),LSIZ)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL MEMREL('TRDR23',WORK,1,KFTMAT,KFREE,LFREE)
         IOFF3=ICNT1
      ENDDO
C
      CALL CONMP2(WORK(KANT),E1,E2,E3,E4,NSTRT(1),NSTRT(2),
     &     NSTRT(3),NSTRT(4),NSTR,IBEQJB,IPRINT)
      CALL TIMEC('CONMP2',TIMSTR,TIMEND) 
      TIM1(5) = TIM1(5) + TIMEND - TIMSTR
C     
C     Release all memory and exit
C
      CALL MEMREL('TRDR23',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('TRDR23')
C
      RETURN
C
 1000 FORMAT (/' Storing all 34-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1002 FORMAT (/' Storing all final integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TQTRNS */
      SUBROUTINE TQTRNS (WORK,KFREE,LFREE,IPRINT,IC,IHM,
     &     INDX,INDXAB,NFPCK12,HMAT,ITQMT,ICMOQR,NDMOQR,
     &     TQMAT,NSTR3,IQP34,NZ34,Q3)
C
C     Written by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <symmet.h>
C
      PARAMETER (D0 = 0.0D00, D1=1.0D00)
      INTEGER   NSTR3(2)
      DIMENSION NDMOQR(2,2),ICMOQR(2)
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5)
      DIMENSION NFPCK12(2),IQP12(4),ITQMT(0:7,0:7,2)
      DIMENSION IQP34(4,0:7),NZ34(0:7)
      DIMENSION IHM(0:7)
      DIMENSION HMAT(*),TQMAT(*),Q3(*),QCOF(4)
      DIMENSION WORK(*)
      DIMENSION IND(2)
C
#include <ibtfun.h>
C     
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL QENTER('TQTRNS')
C
C     Define HMAT as real matrices.
C
      NZ12 = 1
      IQP12(1) = 1
      IQP12(2) = 0
      IQP12(3) = 0
      IQP12(4) = 0
C
      KFRSAV = KFREE
      NIJMAX = 0
C     
C     TQMAT is stored as (ij,iz1,k,s,iz2,ic,Gs,Grs).
C     For real groups this is identical to
C     (ij,iz1,k,iz2,s,ic,Gs,Grs) which is generated
C     directly in QAXPY. For non-real groups we first
C     generate (ij,iz1,k,iz2* which is accumulated in TQMAT.
C     The correct format for TQMAT is needed in LQTRNS
C     for QGEMM: (n,s,iz)*(s,l,iz) = (n,l,iz)
C
      DO IFRP=1,NFSYM
         IF(NFPCK12(IFRP)*NZ**2.GT.NIJMAX) 
     &        NIJMAX = NFPCK12(IFRP)*NZ**2
      ENDDO
      CALL MEMGET('REAL',KT3MAT,NIJMAX,WORK,KFREE,LFREE)
C
      DO20 INDAR = 1, NINSHA
         DO 10 INDBR = 1, NINSHB
            INDAB  = INDXAB(INDAR,INDBR,1)
            IF (INDAB.EQ.0) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
C     Pointer AP	
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            NIJNZ1 = NFPCK12(IREPIJ)*NZ
            NIJNZ2 = NFPCK12(IREPIJ)*NZ**2
C     IOFF for HMAT
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
C     This is actual index R and S in this irrep and IC.
            INDA1 = INDA - ICOS(IREPA+1,IC)
            INDB1 = INDB - ICOS(IREPB+1,IC)
C
C     Loop over K-transformations. 
C     Accumulate first C_rk(rs), then C_sk(rs).
C     Don't accumulate twice if functions A.EQ.B
C
            KLST = 2
            IF (INDA.EQ.INDB) KLST = 1
            INDR = INDA1
            INDS = INDB1
            IREPR = IREPA
            IREPS = IREPB
            DO KTYP = 1,KLST
               IREPK = JBTOF(IREPR,IC)
               NK = NSTR3(IREPK)
               NS = NBBAS(IREPS,IC)
               NIJKS = NK*NS*NIJNZ1
C     
C     Offset for TQMAT.
C     
               JOFF1 = ITQMT(IREPAB,IREPS,IC) 
     &              + (INDS-1)*NK*NIJNZ1+1
C     
C     Offset for Coeffs
C     
               KOFF1 = ICMOQR(IREPK) + IBBAS(IREPR,IC)
     &               - IBAS(IREPK) + INDR - 1
C
C     Dimensions of Q3
C
               NRQ1 = NDMOQR(1,IREPK)
               NCQ1 = NDMOQR(2,IREPK)
               IF (IC.EQ.1) THEN
                  IREPR2 = IREPR
               ELSE
                  IREPR2 = IBTXOR(IXYZ,IREPR)
               ENDIF
C
C     Making IQP34 and NZ34
C
               CALL IQPACK(IQP12,NZ12,IPQTOQ(1,IREPR2),
     &              NZ,IQP34(1,IREPR2),NZ34(IREPR2))
C     
C     Start looping over K index
C
               DO K = 1,NK
C     
                  JOFF2 = JOFF1 + (K-1)*NIJNZ1
                  KOFF = KOFF1 + (K-1)*NRQ1
                  DO I = 1,NZ
                     QCOF(I) = Q3(KOFF+(I-1)*NRQ1*NCQ1)
                  ENDDO
C     
C     Non-real groups
C
                  IF (NZ.GT.1) THEN
                     CALL DZERO(WORK(KT3MAT),NIJNZ2)
                     CALL QAXPY(NIJNZ1,'C','N',QCOF,IPQTOQ(1,IREPR2),
     &                    NZ,HMAT(IOFF),NIJNZ1,1,IQP12,NZ12,
     &                    WORK(KT3MAT),NIJNZ1,1,IQP34(1,IREPR2),
     &                    NZ34(IREPR2))
C     
                     DO IZ = 1,NZ
                        JOFF = JOFF2 + (IZ-1)*NIJKS
                        CALL DAXPY(NIJNZ1,D1,
     &                       WORK(KT3MAT+NIJNZ1*(IZ-1)),
     &                       1,TQMAT(JOFF),1)
                     ENDDO
C
C     Real groups
C
                  ELSE
                     CALL QAXPY(NIJNZ1,'C','N',QCOF,IPQTOQ(1,IREPR2),
     &                    NZ,HMAT(IOFF),NIJNZ1,1,IQP12,NZ12,
     &                    TQMAT(JOFF2),NIJNZ1,1,IQP34(1,IREPR2),
     &                    NZ34(IREPR2))
                  ENDIF
               ENDDO
C     
C     Accumulate C_sk(rs)
C     
               INDR = INDB1
               INDS = INDA1
               IREPR = IREPB
               IREPS = IREPA
            ENDDO
 10      CONTINUE
 20   CONTINUE
      CALL MEMREL('TQTRNS',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TQTRNS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LQTRNS */
      SUBROUTINE LQTRNS (WORK,KFREE,LFREE,IPRINT,
     &     NFPCK12,NFPCK34,ITQMT,ICMOQR,NDMOQR,IQMO,NZMO,
     &     TQMAT,NSTR3,NSTR4,IQP34,NZ34,Q4,FQMAT,IFPCK,
     &     ICS,ICF)
C
C     Written by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <dcbbas.h>
#include <symmet.h>
C
      PARAMETER (D0 = 0.0D00, D1=1.0D00)
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION NDMOQR(2,2),ICMOQR(2)
      DIMENSION NFPCK12(2),NFPCK34(2),ITQMT(0:7,0:7,2),IFPCK(2,2)
      DIMENSION IQP34(4,0:7),NZ34(0:7)
      DIMENSION IQMO(4,0:7),NZMO(0:7)
      DIMENSION TQMAT(*),Q4(*),FQMAT(*)
      DIMENSION WORK(*)
C
#include <ibtfun.h>
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL QENTER('LQTRNS')
C
C     NZ34 and IQP34 imported from 3/4 transformation.
C
      IOFF=1
      DO IREPIJ = 1,NFSYM
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
         DO IREPRS = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
               DO IREPS = 0, NBSYM-1
                  DO IC = ICS,ICF
                     IREPR = IBTXOR(IREPRS,IREPS)
                     IREPK = JBTOF(IREPR,IC)
                     IREPL = JBTOF(IREPS,IC)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     N1 = NIJ*NK*NZ**2
                     N1 = NIJ*NK*NZ
                     N2 = NBBAS(IREPS,IC)
                     N3 = NL
                     IF (NBBAS(IREPR,IC).EQ.0) GOTO 100
C     
C     Offset for TQMAT.
C     
                     JOFF = ITQMT(IREPRS,IREPS,IC) + 1
C     
C     Offset for FQMAT.
C     
                     IOFF1 = IOFF + 
     &                    IFPCK(IREPK,IREPL)*NFPCK12(IREPIJ)*NZ**2
C
C     Offset for Coeffs
C     
                     KOFF = ICMOQR(IREPL) + IBBAS(IREPS,IC) 
     &                    - IBAS(IREPL)
C     
C     Dimensions of Q4
C     
                     NRQ1 = NDMOQR(1,IREPL)
                     NCQ1 = NDMOQR(2,IREPL)
C
                     IF (IC.EQ.1) THEN
                        IREPR2 = IREPR
                     ELSE
                        IREPR2 = IBTXOR(IXYZ,IREPR)
                     ENDIF
C
                     IF (IC.EQ.1) THEN
                        IREPS2 = IREPS
                     ELSE
                        IREPS2 = IBTXOR(IXYZ,IREPS)
                     ENDIF
C     
                     IF (N1*N2.GT.0) THEN
                        CALL IQPACK(IQP34(1,IREPR2),NZ34(IREPR2),
     &                       IPQTOQ(1,IREPS2),NZ,IQMO(1,IREPRS),
     &                       NZMO(IREPRS))
C
                        CALL QGEMM(N1,N3,N2,D1,
     &                       'N','N',IQP34(1,IREPR2),TQMAT(JOFF),
     &                       N1,N2,NZ34(IREPR2),
     &                       'N','N',IPQTOQ(1,IREPS2),Q4(KOFF),
     &                       NRQ1,NCQ1,NZ,
     &                       D1,IQMO(1,IREPRS),FQMAT(IOFF1),
     &                       N1,N3,NZMO(IREPRS))
                     ENDIF
 100                 CONTINUE
                  ENDDO
               ENDDO
               IOFF = IOFF + NIJ*NKL*NZ*NZ
            ENDIF
         ENDDO
      ENDDO
      CALL QEXIT('LQTRNS')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MP2FI1 */
      SUBROUTINE MP2FI1(ICLASS,IREPKL,NSTR,KKR,LKR,FMAT,SPMAT,
     &     NFM,NI,NJ,NK,NL,LS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION NSTR(2,0:2,4)
      DIMENSION FMAT(NFM,2)
      DIMENSION SPMAT(NI,NJ,NK,NL,LS)
C
c       WRITE(LUPRI,'(A,4I4)') 'MP2FI1:', NI,NJ,NK,NL
      M = 0
      JJ = 0
      DO IREPJ = 1, NFSYM
         IREPI = MOD(IREPJ+IREPKL,2) + 1
         DO J = 1, NSTR(IREPJ,0,2)
            JJ = JJ + 1
            DO I = 1, NSTR(IREPI,0,1)
               IF (IREPI.EQ.1) THEN
                  II = I
               ELSE
                  II = I + NSTR(1,0,1)
               ENDIF
               M = M + 1
c               WRITE(LUPRI,'(I4,A,4I4,2X,F20.10)') 
c     &              ICLASS,'##',II,JJ,KKR,LKR,FMAT(M,1)
               SPMAT(II,JJ,KKR,LKR,ICLASS)  =FMAT(M,1)
               SPMAT(II,JJ,KKR,LKR,ICLASS+4)=FMAT(M,2)
c               print*,'MP2AA',FMAT(M,1),FMAT(M,2),II,JJ,KKR,LKR,ICLASS
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck timec */
      SUBROUTINE TIMEC(TEXT,TIMSTR,TIMEND)
#include <implicit.h>
#include <priunit.h>
      CHARACTER TEXT*6
C
      IF (TEXT .EQ. 'START ') THEN
         TIMSTR = SECOND()
      ELSE
         TIMEND = SECOND()
         TIME   = TIMEND - TIMSTR
C        CALL TIMTXT(' Time used in '//TEXT//' is',TIME, LUPRI)
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck tskchk */
      SUBROUTINE TSKCHK(MAXTSK,NSTR,LWORK,LENTSK,NSTRT,IPRINT)
C
C     Written by Jon Laerdahl, June 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <dcbmpt.h>
      DIMENSION NSTRT(4),NSTR(2,0:2,4)
C
C     Maximum number of active occupied spinors in each
C     task is either set in input (IJTSK) or
C     on the basis of avaliable memory.
C
      WRITE(LUPRI,'(A,I10/)') 
     &     'Memory available for Direct MP2 is (Words):',
     &     LWORK
      IF (IJTSK.GT.0) THEN
         LENTSK = IJTSK
         IF(IJTSK.GT.NSTRT(1)) LENTSK = NSTRT(1)
         IF (IPRINT.GE.0) THEN
            WRITE(LUPRI,'(A,/A,I10)') 
     &           'Maximum number of occupied spinors I and J in batch',
     &           'set in input:', IJTSK
         ENDIF
      ELSE
         LENTSK = NSTRT(1)
         WRITE(LUPRI,'(A,/A,I6,A)') 
     &        'Maximum number of occupied spinors in batch not set'//
     &        ' in input.', 
     &        'Testing ALL spinors :',NSTRT(1),' in batch.'
         IF (IPRINT.GE.1) THEN
            WRITE(LUPRI,'(/A,/A,/A)') 
     &           '### If not enough memory is avaliable:          ###',
     &           '### Reduce memory requirements by setting       ###',
     &           '### .IJTSK to a lower value                     ###'
         ENDIF
      ENDIF
      IOLD = LENTSK
C
C     Reduce LENTSK until we have enough space...
C     May be done automatically when memory requirements are known.
C
c      IF (LENTSK.NE.IOLD) THEN
c         IF (IPRINT.GE.1) THEN
c            WRITE(LUPRI,'(/A,I5)') 'Batch Size Shrinked To:', LENTSK
c            WRITE(LUPRI,'(A,I10)') 'Memory Needed In Transformation:',
c     &           INEED
c         ENDIF
c      ENDIF
      IF (IPRINT.EQ.1) THEN
         WRITE(LUPRI,'(/A,I5)') 
     &        'Total Number of Active Occupied Spinors:', NSTRT(1)
         WRITE(LUPRI,'(A,I5)') 
     &        'Number in Max Batch:', LENTSK
      ENDIF
      IF ((MOD(NSTRT(1),LENTSK)).EQ.0) THEN
         NUM = NSTRT(1)/LENTSK
         IF (IPRINT.GE.1) 
     &        WRITE(LUPRI,'(/A)') 'All Batches have same Size.'
      ELSE
         NUM = NSTRT(1)/LENTSK + 1
         IF (IPRINT.GE.1) 
     &        WRITE(LUPRI,'(/A)') 'Last Batch has reduced Size.'
      ENDIF
      MAXTSK = NUM*(NUM+1)/2
      IF (IPRINT.GE.1) THEN
         WRITE(LUPRI,'(/A,I5)')'Number of tasks:', MAXTSK
         WRITE(LUPRI,'(A,I5)') 'Number of Batches of I:', NUM
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN1E*/
      SUBROUTINE MS4IN1E(WORK,KFREE,LFREE,IPRINT,LGREC,ICS,ICF,
     &                   NSTR1,NSTR2,NFPCK12,IFPCK12,
     &                   NDMOQR,ICMOQR,IJPASS,Q1,Q2,
     &                   GMAT,NSIZG,NZMO,IQMO,HMAT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher May 1997.
C
C     PURPOSE : Driver for the first index transformation.
C               Loop over batches for a given shell combination and
C               write the half-transformed integrals to memory
C
C     Input :
C
C     - IPRINT       Print flag
C     - NPASS        Number of passes through the scalar integrals
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK12      Number of spinor pairs for each compound symmetry
C     - IFPCK12      Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPASS       Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C     - NSIZG        Size of GMAT array
C     - MTHM         Number of scalar integral pairs that contribute
C                    to a given fermion ircop (multiplied with NZ)
C     - ITHM         Start address of a boson representation on the
C                    file of half-transformed integrals.
C     - INDX         For each boson function the boson irrep and the
C                    position in the shell
C
C     Output :
C
C     - NZMO         Number of non-zero quaternion components of HMAT
C     - IQMO         Quaternion packing of HMAT
C     -              Updated values of buffer blocks
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(NSIZG)
      DIMENSION NFPCK12(2),IFPCK12(2,2)
      DIMENSION IJPASS(0:7,NPASS)
      DIMENSION IQAO(4),IQMO(4,0:7),NZMO(0:7)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
      DIMENSION HMAT(*)
C
C     For the buffered input
C
      DIMENSION LGREC(NPASS)
C
C     For the buffered output
C
      DIMENSION IHM(0:7)
C
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN1E')
      KFRSAV = KFREE
C
      IHOFF = 1
      DO IPASS = 1, NPASS
C
C        Get the scalar integrals that are processed in this pass
C        The integral are already in GMAT when NPASS = 1
C
         IF (NPASS.GT.1) THEN
             CALL DZERO(GMAT,NSIZG)
             IREC = LGREC(IPASS)
             CALL MEMGET('REAL',KRBUF,NGBFSZ,WORK,KFREE,LFREE)
             CALL MEMGET('INTE',KLBUF,NGBFSZ,WORK,KFREE,LFREE)
             CALL READGBF(LGFIL,IREC,WORK(KRBUF),WORK(KLBUF),NGBFSZ,
     &                    GMAT)
             CALL MEMREL('MS4IN1E',WORK,KFRSAV,KRBUF,KFREE,LFREE)
         ENDIF
C
C        Calculate the size of the H-matrix for this pass
C
         NSIZH = 0
         DO IREPAB = 0, NBSYM-1
            IHM(IREPAB) = NSIZH
            IREPIJ = JBTOF(IREPAB,1)
            NSIZH = NSIZH + IJPASS(IREPAB,IPASS)*NFPCK12(IREPIJ)
         ENDDO
C
C        Do first step of 4-index transformation :
C        Transform first pair of indices
C
         CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &                .FALSE.,NSTR1,NSTR2,NFPCK12,IFPCK12,NDMOQR,
     &                ICMOQR,IJPASS(0,IPASS),Q1,Q2,
     &                GMAT,HMAT(IHOFF),NZMO,IQMO)
         IHOFF = IHOFF + NSIZH*NZ
C
C        We can also rewrite code and call the 3/4 transformation here.
C        This reduces the memory requirements for HMAT.
C
      ENDDO
C
      CALL QEXIT('MS4IN1E')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRIEIN */
      SUBROUTINE PRIEIN(E1,E2,E3,E4,N1,N2,N3,N4,NSTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
      DIMENSION E1(N1),E2(N2),E3(N3),E4(N4)
      DIMENSION NSTR(2,0:2,4)
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 1, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E1(NX+IXX), IXX=1,NSTR(IFRP,0,1))
         NX = NX + NSTR(IFRP,0,1)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 2, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E2(NX+IXX), IXX=1,NSTR(IFRP,0,2))
         NX = NX + NSTR(IFRP,0,2)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 3, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E3(NX+IXX), IXX=1,NSTR(IFRP,0,3))
         NX = NX + NSTR(IFRP,0,3)
      ENDDO
      NX = 0
      DO IFRP = 1, NFSYM
         WRITE(LUPRI,'(A,I1)') 'Index 4, Irrep ',IFRP
         WRITE(LUPRI,'(5D15.8)') 
     &        (E4(NX+IXX), IXX=1,NSTR(IFRP,0,4))
         NX = NX + NSTR(IFRP,0,4)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CONMP2 */
      SUBROUTINE CONMP2(SPMAT,E1,E2,E3,E4,N1,N2,N3,N4,NSTR,
     &     IBEQJB,IPRINT)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <dcbmp2.h>
#include <dcbmpt.h>
      LOGICAL IBEQJB
      DIMENSION NSTR(2,0:2,4)
      DIMENSION E1(N1),E2(N2),E3(N3),E4(N4)
      DIMENSION SPMAT(N1,N2,N3,N4,16)
      DIMENSION ICLDAT(16)
      DATA (ICLDAT(I),I=1,16) /1,1,2,4,3,3,4,2,9,9,10,12,11,11,12,10/
      IF (IPRINT.LT.10) GOTO 123
      DO ICL = 1,4
         WRITE(LUPRI,'(/A,I5/)') 'CLASS:', ICL
         ICL2 = ICL + 4
         DO II = 1,N1
            DO IA = 1,N2
               DO IJ = 1,N3
                  DO IB = 1,N4
                     WRITE(LUPRI,'(4I4,2F20.10)')
     &                    II,IA,IJ,IB,SPMAT(II,IA,IJ,IB,ICL),
     &                    SPMAT(II,IA,IJ,IB,ICL2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO ICL = 9,12
         WRITE(LUPRI,'(/A,I5/)') 'CLASS:', ICL
         ICL2 = ICL + 4
         DO II = 1,N1
            DO IA = 1,N2
               DO IJ = 1,N3
                  DO IB = 1,N4
                     WRITE(LUPRI,'(4I4,2F20.10)')
     &                    II,IA,IJ,IB,SPMAT(II,IA,IJ,IB,ICL),
     &                    SPMAT(II,IA,IJ,IB,ICL2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
 123  CONTINUE
      IF (TESTLS.OR.TESTSL) THEN
         WRITE(LUPRI,'(A)') 
     &        '(LL|SS)- or (SS|LL)-ONLY test does only work without'
         WRITE(LUPRI,'(A)') 
     &        'call to ENERC1. Recompile and remove these lines....'
         WRITE(LUPRI,'(A)') 
     &        'Do NOT loop over more than single batch of I and J!'
         CALL QUIT('CONMP2: Recompile to do TESTSL or TESTLS!')
      ENDIF
C     
C     Real and complex groups, Class 1..8
C     
      EC = 0d00
      DO ICL = 0,3
         ICL1= ICL*2+1
         ICL2= ICL1+1
c         CALL ENERC1(E,E1,E2,E3,E4,N1,N2,N3,N4,
         CALL ENERC2(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &        SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &        SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &        SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &        SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,
     &        IPRINT)
         EC=EC+E
      ENDDO
C
C     Quaternion groups, Class 9..16
C
      IF (NZ.EQ.4) THEN
         DO ICL = 4,7
            ICL1= ICL*2+1
            ICL2= ICL1+1
c            CALL ENERC1(E,E1,E2,E3,E4,N1,N2,N3,N4,
            CALL ENERC2(E,E1,E2,E3,E4,N1,N2,N3,N4,
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL1)+4),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)),
     &           SPMAT(1,1,1,1,ICLDAT(ICL2)+4),IBEQJB,IPRINT)
            EC=EC+E
         ENDDO
      ENDIF
      EMP2 = EMP2 + EC
      IF(IPRINT.GE.1) THEN
         WRITE(LUPRI,'(A,D20.10)') 'Contribution, MP2', EC
         WRITE(LUPRI,'(A,D20.10)') 'Accumulated, MP2 ', EMP2
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ENERC2 */
      SUBROUTINE ENERC2(E,EIG1,EIG2,EIG3,EIG4,N1,N2,N3,N4,
     & FMR1,FMI1,FMR2,FMI2,IBEQJB,IPRINT)
C
C     Written by Jon K. Laerdahl 1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      INTEGER I,A,J,B
      LOGICAL IBEQJB
      DIMENSION
     &     FMR1(N1,N2,N3,N4),FMI1(N1,N2,N3,N4),
     &     FMR2(N1,N2,N3,N4),FMI2(N1,N2,N3,N4)
      DIMENSION EIG1(N1),EIG2(N2),EIG3(N3),EIG4(N4)
C
C     IBEQJB is true if i and j are in same block
C     a and b should always be in same block!
C     -> NE2 = NE4 AND same range too! 
C
      E = 0.D0
      DO 100 I=1,N1
         DO 110 A=1,N2
            IF (IBEQJB) THEN
               JTOP = I
            ELSE
               JTOP = N3
            ENDIF
            DO 120 J=1,JTOP
               DO 130 B=1,A
c                  WRITE(LUPRI,'(/A,4I4)') 'aa', I,A,J,B
                  AR = FMR1(I,A,J,B)-FMR2(I,B,J,A)
                  AI = FMI1(I,A,J,B)-FMI2(I,B,J,A)
c                  WRITE(LUPRI,'(6F15.8)') 
c    &                 FMR1(I,A,J,B),FMR2(I,B,J,A),
c    &                 FMI1(I,A,J,B),FMI2(I,B,J,A),AR,AI
                  EN = AR**2 + AI**2
                  DIV = EIG1(I)+EIG3(J)-EIG2(A)
     &  -EIG4(B)
cc                  IF (DABS(EN).GT.1D-15)
c          WRITE(LUPRI,'(6F15.8)') 
c     & DIV,EIG1(I),EIG1(J),EIG2(A),EIG2(B),EN
                  DELIJ = 0.d00
                  IF ((I.EQ.J).AND.IBEQJB) DELIJ = 1.d00
                  IF (A.EQ.B) THEN
                     DELAB = 1.d00
                  ELSE
                     DELAB = 0.d00
                  ENDIF
                  E = E + (2.d00-DELIJ)*(2.d00-DELAB)*EN/(2.d00*DIV)
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF(IPRINT.GE.1) 
     &     WRITE(LUPRI,'(A,D25.16)')  'MP2-energy, 1/4 matrix    :',E
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck ENERC1 */
      SUBROUTINE ENERC1(E,EIG1,EIG2,EIG3,EIG4,N1,N2,N3,N4,
     &     FMR1,FMI1,FMR2,FMI2,IBEQJB,IPRINT)
C     
C     Written by Jon K. Laerdahl 1995
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      INTEGER I,A,J,B
      LOGICAL IBEQJB
      DIMENSION
     &     FMR1(N1,N2,N3,N4),FMI1(N1,N2,N3,N4),
     &     FMR2(N1,N2,N3,N4),FMI2(N1,N2,N3,N4)
      DIMENSION EIG1(N1),EIG2(N2),EIG3(N3),EIG4(N4)
C     
      E = 0.D0
      DO 100 I=1,N1
         DO 110 A=1,N2
            DO 120 J=1,N3
               DO 130 B=1,N4
                  AR = FMR1(I,A,J,B)-FMR2(I,B,J,A)
                  AI = FMI1(I,A,J,B)-FMI2(I,B,J,A)
                  EN = AR**2 + AI**2
                  DIV = EIG1(I)+EIG3(J)-EIG2(A)-EIG4(B)
                  E = E + EN/(2.d00*DIV)
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF(IPRINT.GE.1) 
     &     WRITE(LUPRI,'(A,D25.16)')  'MP2-energy, 1/4 matrix    :',E
      RETURN
      END




