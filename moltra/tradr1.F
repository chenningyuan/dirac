C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : tradr1.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR1T */
      SUBROUTINE TRDR1T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &                  NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2)
C
C     Written by Luuk Visscher, Februari 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <dcbgen.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dgroup.h>
      PARAMETER (MAXTSK=1)
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),NSTR1(2,0:2,4,MAXTSK),INDX(3,*)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION NSTRS(2,4,MAXTSK),ICMOQR1(2,4)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDR1T')
      KFRSAV = KFREE
C
C     Divide index transformation in task with different batches
C     of spinors
C
C     CALL DRV1T (WORK,KFREE,LFREE,MAXTSK,NSTR,NTASK,NSTR1,NSTRS)
C
C     NOTE 1 :
C     DRV1T remains to be written, number of tasks will depend on memory
C     availabe and required for the three-index transformed integrals
C     At the moment : assume that we need only one run through.
C     
      ntask = 1
      call icopy (24,NSTR(1,0,1),1,NSTR1(1,0,1,1),1)
      call icopy (8,1,0,NSTRS,1)
C
C     NOTE 2 :
C     Obviously the best place to parallelize this scheme is over this
C     loop, if the number of task is sufficiently large compared to the 
C     number of processors. More fine-grain parallelization is possible
C     by parallelizing over the calls to CALDIS and MS4IN1 in the actual
C     driver TRDRV1
C
      DO ITASK = 1, NTASK
C
C        Give relative offset for the coefficients
C
         DO INDEX = 1, 4
            DO I = 1, NFSYM
               ICMOQR1(I,INDEX) = ICMOQR(I,INDEX)
     &                          + (NSTRS(I,INDEX,ITASK)-1)*NFBAS(I,0)
            ENDDO
         ENDDO
C
         CALL TRDRV1(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR1,
     &               NSTR1(1,0,1,ITASK),ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,
     &               E1,E2)
C
      ENDDO
C
      CALL MEMREL('TRDR1T',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDR1T')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDRV1 */
      SUBROUTINE TRDRV1(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,
     &                  NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2)
C
C     Written by Luuk Visscher, december 1996
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
C
#include <mxorb.h>
#include <dcbgen.h>
#include <dcbtra.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAL ANTIS, LMP2, NOPV, NODV
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION IQMO(4,0:7),NZMO(0:7)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)
      DIMENSION IJP12(0:7)
      DIMENSION IJP34(2)
      DIMENSION NFPCK12(2),IFPCK12(2,2)
      DIMENSION NFPCK34(2),IFPCK34(2,2)
      DIMENSION IHTM(0:7),IHM(0:7),IFM(0:7),NHTM(2)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*),E1(*),E2(*)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDRV1')
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
C     Make index array for half-transformed integrals
C
      NHTM(1) = 0
      NHTM(2) = 0
      NHTMT = 0
      DO IREPIJ = 1, NFSYM
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IHTM(IREPAB) = NHTM(IREPIJ) + NHTMT
               NHTM(IREPIJ) = NHTM(IREPIJ) + NSPCK(IREPAB,0)*NZ
            ENDIF
         ENDDO
         NHTMT = NHTMT + NHTM(IREPIJ)*NFPCK12(IREPIJ)
      ENDDO
C
C     Allocate memory for half transformed integrals
C
      IF (IPRINT .GE. 0) WRITE (LUPRI,1000) NHTMT*8./(1024.*1024.)
      CALL MEMGET('REAL',KHTMAT,NHTMT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KHTMAT),NHTMT)
C
C     Loop over (XX|LL) and (XX|SS) integrals
C
C     ----------------
      DO 100 IC = 1, 2
C     ----------------
C
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
            NINSHA  = NINSH(ISHLA,-1,INDX,IPRINT)
            NINSHB  = NINSH(ISHLB,-1,INDX,IPRINT)
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
            LDXAB = 5*NINSHA*NINSHB
            CALL MEMGET('INTE',KDXAB,LDXAB,WORK,KFREE,LFREE)
            NSIZG = LFREE
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.FALSE.,.TRUE.,
     &                       INDX,WORK(KDXAB),IJP12,NSIZG,IDUM)
            NSIZH = 0
            DO IREPAB = 0, NBSYM-1
               IHM(IREPAB) = NSIZH
               IREPIJ = JBTOF(IREPAB,1)
               NSIZH = NSIZH + IJP12(IREPAB)*NFPCK12(IREPIJ)
            ENDDO
C
            IF (IPRINT .GE. 3) THEN
               WRITE (LUPRI,'(//,4(1X,A,2I4,/))')
     &         ' Calculation of integral distribution set:',ISHLA,ISHLB
               WRITE (LUPRI,'(//,4(1X,A,I8,/))')
     &         ' Number of distributions in this set:     ',NUMDIS,
     &         ' Size of scalar integral batch:           ',NSIZG,
     &         ' Size of half-transformed integral batch: ',NSIZH
            END IF
C
            CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KHMAT),NSIZH*NZ)
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
C
C           Get distributions
C
            NODV = NASHT.EQ.0
            NOPV = NASHT.LT.2
            CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &                  NODV,NOPV,WORK(KFREE),LFREE,IPRINT)
C
C           Do first step of 4-index transformation :
C           Transform first pair of indices
C 
            CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,.FALSE.,
     &                   NSTR(1,0,1),NSTR(1,0,2),
     &                   NFPCK12,IFPCK12,NDMOQR,ICMOQR,IJP12,Q1,Q2,
     &                   WORK(KGMAT),WORK(KHMAT),NZMO,IQMO)
C
C           Scatter half-transformed integrals to the right position
C
            CALL MS4INS (WORK,KFREE,LFREE,IPRINT,IC,IHM,NHTM,IHTM,INDX,
     &                   WORK(KDXAB),NFPCK12,WORK(KHMAT),WORK(KHTMAT))
C
            CALL MEMREL('TRDRV1',WORK,1,KFRSAV3,KFREE,LFREE)
C
         ENDDO
         ENDDO
C
C----------------
 100     CONTINUE
C----------------
C
      NSIZF = 0
      DO IREPAB = 0, NBSYM-1
         IREPIJ = JBTOF(IREPAB,1)
         NSIZF = NSIZF + NFPCK34(IREPIJ)*NZ*NZ*NFPCK12(IREPIJ)
      ENDDO
C
      IF (IPRINT .GE. 0) WRITE (LUPRI,1001) NSIZF*8./(1024.*1024.)
      CALL MEMGET('REAL',KFMAT,NSIZF,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KFMAT),NSIZF)
C
      icount1 = 0
      icount2 = 0
      do ix = 1, nhtmt
c        print*,ix,work(khtmat+ix-1)
         if (abs(work(khtmat+ix-1)).gt.1.0d-16) then
            icount1 = icount1 + 1
         else
            icount2 = icount2 + 1
         endif
      enddo
      print*,' Non-zero half-transformed integrals ',icount1
      print*,' Zero half-transformed integrals     ',icount2
      print*,' Total half-transformed integrals    ',icount1+icount2
      ICS = 1
      ICF = 2
      CALL MS4IN2 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &             NSTR(1,0,3),NSTR(1,0,4),
     &             NFPCK34,IFPCK34,NDMOQR(1,1,3),ICMOQR(1,3),
     &             NFPCK12,Q3,Q4,
     &             WORK(KHTMAT),WORK(KFMAT),NZMO,IQMO)
c     write (LUPRI,'(//A)') ' The fully transformed integrals'
c     do ix = 1, nsizf
c        print*,ix,work(kfmat+ix-1)
c     enddo
C
c     IF (ANTIS) WRITE(LUPRI,'(A)') 'Antisymmetrizing'
c     IF (LMP2) THEN
c        WRITE(LUPRI,'(A)') 'DOING MP2'
C         CALL MP2SEC(E1,E2,NSTR)
c     ENDIF
C
C     Initialize MOLFDIR integral file
C
      CALL MEMGET('INTE',KIKRA,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIKRB,NORBT,WORK,KFREE,LFREE)
      CALL MDINTI(NSTR,WORK(KIKRA),WORK(KIKRB),IPRINT)
      CALL MEMREL('TRDRV1',WORK,1,KIKRA,KFREE,LFREE)
C
C     Number of boson symmetry for each fermion irrep (parity)
C
      NBSYMP = NBSYM/NFSYM
C
C
C     Complete the integrals to the desired classes
C
      IOFF = 0
      DO IREPIJ = 1, NFSYM
         JKR = 0
         IJ = 0
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            DO J = 1, NSTR(IREPJ,0,2)
               JKR = JKR + 1
               DO I = 1, NSTR(IREPI,0,1)
                  IF (IREPI.EQ.1) THEN
                     IKR = I
                  ELSE
                     IKR = NSTR(1,0,1) + I
                  ENDIF
                  IJ = IJ + 1
                  CALL SYMFINT(IPRINT,IREPIJ,IJ,IKR,JKR,
     &                 NSTR(1,0,3),NSTR(1,0,4),
     &                 NFPCK34(IREPIJ),WORK(KFMAT+IOFF),IQMO,
     &                 WORK,KFREE,LFREE)
                  IOFF = IOFF + NFPCK34(IREPIJ)*NZ*NZ*NBSYMP
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     Close integral file
C
      CALL MDINTF
C
C     Release all memory and exit
C
      CALL MEMREL('TRDRV1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDRV1')
C
      RETURN
C
 1000 FORMAT (/' Storing all half-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MS4INS*/
      SUBROUTINE MS4INS (WORK,KFREE,LFREE,IPRINT,IC,IHM,NHTM,IHTM,
     &                   INDX,INDXAB,NFPCK12,HMAT,HTMAT)
C
C     Sort half-transformed integrals
C
C     Luuk Visscher
C
C
#include <implicit.h>
#include <priunit.h>
#include <dgroup.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbtra.h>
#include <dcbbas.h>
C
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,4)
      DIMENSION NFPCK12(2),IFPCK12(2,2)
      DIMENSION IHTM(0:7),IHM(0:7),NHTM(2)
      DIMENSION HMAT(*),HTMAT(*)
      DIMENSION WORK(*)
      DIMENSION IND(2)
C
#include <ibtfun.h>
C
C     Loop over the functions inside the shell-block.
C     We do it this way rather than just over the boson irreps of the
C     density and than over the pairs within that irrep, because we
C     need to access the information in INDXAB. The alternative gives
C     less jumping through memory, but requires introduction of a new 
C     n^2 type array in which indab is stored.
C
      DO 20 INDAR = 1, NINSHA
         DO 10 INDBR = 1, NINSHB
            INDAB  = INDXAB(INDAR,INDBR,1)
            IF (INDAB.EQ.0) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
            INDA1 = INDA - ICOS(IREPA+1,IC)
            INDB1 = INDB - ICOS(IREPB+1,IC)
C
C           The offset within the untransformed indexes. These were
C           the right-hand side (slowest varying) and are sorted to
C           become the left hand side. Symmetrization is also done
C           in this step (the batch indices were triangular).
C           
            JOFFAB = ISPCK(IREPA,IREPB,IC)
     &             + (INDB1-1)*NBBAS(IREPA,IC)+INDA1
            JOFFBA = ISPCK(IREPB,IREPA,IC)
     &             + (INDA1-1)*NBBAS(IREPB,IC)+INDB1
            JOFF1 = IHTM(IREPAB) + JOFFAB
            JOFF2 = IHTM(IREPAB) + JOFFBA
            DO IZ = 1, NZ
C
C              We copy all transformed indices that now become slowest
C              varying (right-hand). This gives a stride.
C
               CALL DCOPY(NFPCK12(IREPIJ),HMAT(IOFF),1,
     &                   HTMAT(JOFF1),NHTM(IREPIJ))
               IF (JOFF1.NE.JOFF2)
     &         CALL DCOPY(NFPCK12(IREPIJ),HMAT(IOFF),1,
     &                   HTMAT(JOFF2),NHTM(IREPIJ))
               IOFF  = IOFF  + NFPCK12(IREPIJ)
               JOFF1 = JOFF1 + NSPCK(IREPAB,0)
               JOFF2 = JOFF2 + NSPCK(IREPAB,0)
            ENDDO
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN2 */
      SUBROUTINE MS4IN2 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR1,NSTR2,
     & NFPCK,IFPCK,NDMOQR,ICMOQR,IJPAIRS,Q1,Q2,GMAT,HMAT,NZMO,IQMO)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C     A twin brother of MS4IN1. Don't get confused by the name of the
C     indices. We DO transform the second index pair.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               Second index-pair transformation.
C   
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (XX|LL), 2 : (XX|SS)
C     - ICL          Last class of integrals
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK        Number of spinor pairs for each compound symmetry
C     - IFPCK        Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPAIRS      Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C
C     Output :
C
C     - HMAT         Symmetry packed half-transformed integrals
C     - NZMO         Number of non-zero quaternion components of HMAT
C     - IQMO         Quaternion packing of HMAT
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(*),HMAT(*)
      DIMENSION NFPCK(2),IFPCK(2,2)
      DIMENSION IJPAIRS(2)
      DIMENSION IQAO(4),IQMO(4,0:7),NZMO(0:7)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN2')
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C
      NZAO = 1
      IQAO(1) = 1
      IQAO(2) = 0
      IQAO(3) = 0
      IQAO(4) = 0
C
C     Transform righthand indices to molecular spinor basis
C
      IOFF = 1
      JOFF = 1
      DO IREPIJ = 1, NFSYM
        CALL MEMGET ('REAL',KHTMP,NFPCK(IREPIJ)*NZ,WORK,KFREE,LFREE)
        CALL DZERO (WORK(KHTMP),NFPCK(IREPIJ)*NZ)
        DO IJ = 1, IJPAIRS(IREPIJ)
          DO IREPPQ = 0, NBSYM-1
          IF (IREPIJ.EQ.JBTOF(IREPPQ,1)) THEN
            DO IZ = 1, NZ
               DO IC = ICS, ICF
                  IF (ICS.EQ.2) IOFF = IOFF + NSPCK(IREPPQ,1)
                  DO IREPQ = 0, NBSYM-1
                     IREPP = IBTXOR(IREPPQ,IREPQ)
                     IREPI =  JBTOF(IREPP,IC)
                     IREPJ =  JBTOF(IREPQ,IC)
                     KOFF = ICMOQR(IREPI,1) + IBBAS(IREPP,IC) 
     &                    - IBAS(IREPI)
                     LOFF = ICMOQR(IREPJ,2) + IBBAS(IREPQ,IC) 
     &                    - IBAS(IREPJ)
                     NP = NBBAS(IREPP,IC)
                     NQ = NBBAS(IREPQ,IC)
                     NI = NSTR1(IREPI)
                     NJ = NSTR2(IREPJ)
                     NRQ1 = NDMOQR(1,IREPI,1)
                     NCQ1 = NDMOQR(2,IREPI,1)
                     NRQ2 = NDMOQR(1,IREPJ,2)
                     NCQ2 = NDMOQR(2,IREPJ,2)
                     IF (IC.EQ.1) THEN
                        IREPPI = IREPP
                        IREPQI = IREPQ
                     ELSE
                        IREPPI = IBTXOR(IXYZ,IREPP)
                        IREPQI = IBTXOR(IXYZ,IREPQ)
                     ENDIF
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                      GMAT(IOFF),NP,NQ,NZAO,IQAO,
     &                      WORK(KHTMP),NI,NJ,NZMO(IREPPQ),
     &                                            IQMO(1,IREPPQ),
     &                      Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                      Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                      WORK(KFREE),LFREE,IPRINT)
C
C                      Order the integrals with NKL as first and IZ2
C                      as second index
C
                       DO IZ2 = 1, NZ
                          JOFF1 = JOFF + (IZ2-1)*NFPCK(IREPIJ) +
     &                                   IFPCK(IREPI,IREPJ)
                          JOFF2 = KHTMP + (IZ2-1)*NI*NJ
                          CALL DAXPY(NI*NJ,D1,WORK(JOFF2),1,
     &                               HMAT(JOFF1),1)
                       ENDDO
                     ENDIF
                     IOFF = IOFF + NP * NQ
                  ENDDO
               ENDDO
               JOFF = JOFF + NFPCK(IREPIJ)*NZ
            ENDDO
          ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      CALL QEXIT('MS4IN2')
C
      RETURN
      END
