C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : tradr4.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR4T */
      SUBROUTINE TRDR4T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                  NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2)
C
C     Written by Luuk Visscher, Februari 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include <fpvm3h.h>
#endif
#include <infpar.h>
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <dgroup.h>
      LOGICAL ANTIS,LMP2,TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION IQMO(4,0:7)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDR4T')
      KFRSAV = KFREE
C
C     Get hold of the slaves if we run in parallel.
C
#if defined (VAR_PVM) || defined (VAR_MPI)
      IF (PARCAL) THEN
         ITASK = 2
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
         CALL TRAPARI(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,IDUM,IDUM,
     &                IDUM)
      ENDIF
#endif
C
      TRIAN(1) = ISAME(2).EQ.ISAME(1)
      TRIAN(2) = .FALSE.
C
C     Define the packing of the 2-index transformed integrals
C
      NFPCK12T = 0
      NFPCK34T = 0
      DO IREPIJ =1, NFSYM
         NFPCK12(IREPIJ) = 0
         NFPCK34(IREPIJ) = 0
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
C
C           Packing can be triangular (i=>j) or square.
C
            IF (.NOT.TRIAN(1).OR.(TRIAN(1).AND.IREPI.GT.IREPJ)) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ) = NFPCK12(IREPIJ) +
     &                           NSTR(IREPI,0,1)*NSTR(IREPJ,0,2)
            ELSEIF (TRIAN(1).AND.IREPI.EQ.IREPJ) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ) = NFPCK12(IREPIJ) +
     &                           NSTR(IREPI,0,1)*(NSTR(IREPI,0,1)+1)/2
            ELSE
               IFPCK12(IREPI,IREPJ) = - 1
            ENDIF
            IF (.NOT.TRIAN(2).OR.(TRIAN(2).AND.IREPI.GT.IREPJ)) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ) = NFPCK34(IREPIJ) +
     &                           NSTR(IREPI,0,3)*NSTR(IREPJ,0,4)
            ELSEIF (TRIAN(2).AND.IREPI.EQ.IREPJ) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ) = NFPCK34(IREPIJ) +
     &                           NSTR(IREPI,0,3)*(NSTR(IREPI,0,3)+1)/2
            ELSE
               IFPCK34(IREPI,IREPJ) = - 1
            ENDIF
         ENDDO
         NFPCK12T = NFPCK12T + NFPCK12(IREPIJ)
         NFPCK34T = NFPCK34T + NFPCK34(IREPIJ)
      ENDDO
C
      IF (PARCAL) THEN
C
C     We send out all the tasks (shells in this case) to the slaves
C
         DO IC = 1, 2
            CALL SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
            DO ISHLA = IASTRT, IASMAX
#if defined (VAR_MPI)
               CALL MPI_RECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                       MPI_COMM_WORLD,ISTAT,IERR)
               CALL MPI_SEND(ISHLA,1,MPI_INTEGER,NODE,30,
     &                       MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(IC,1,MPI_INTEGER,NODE,30,
     &                       MPI_COMM_WORLD,IERR)
#endif
            ENDDO
         ENDDO
C
C     Get the slaves out of the transformation loop
C
         NTEST = -1
         DO I = 1,NODES
#if defined (VAR_MPI)
            CALL MPI_RECV(NODE,1,MPI_INTEGER,MPI_ANY_SOURCE,20,
     &                    MPI_COMM_WORLD,ISTAT,IERR)
            CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                    MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(NTEST,1,MPI_INTEGER,NODE,30,
     &                    MPI_COMM_WORLD,IERR)
#endif
         ENDDO
C
C     Serial code 
C
      ELSE
C
C     Initialize scratch output file for final integrals
C
         NODE = 0
         CALL MDSCRI(WORK,KFREE,LFREE,NFPCK12,NFPCK34,NODE)
C
         DO IC = 1, 2
            CALL SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
            DO ISHLA = IASTRT, IASMAX
               CALL TRDRV4(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                     NSTR,ANTIS,LMP2,INDX,TRIAN,
     &                     Q1,Q2,Q3,Q4,E1,E2,IC,IQMO)
            ENDDO
         ENDDO
      ENDIF
C
C     Initialize MOLFDIR integral file
C
      CALL MEMGET('INTE',KIKRA,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIKRB,NORBT,WORK,KFREE,LFREE)
      CALL MDINTI(NSTR,WORK(KIKRA),WORK(KIKRB),IPRINT)
      CALL MEMREL('TRDR4T',WORK,1,KIKRA,KFREE,LFREE)
C
C     Complete the integrals to the desired classes
C
      CALL MEMGET('INTE',KDXKR12,NFPCK12T*2,WORK,KFREE,LFREE)
      CALL MKINDXKR (NSTR,TRIAN(1),WORK(KDXKR12))
C
      NIJOFF = 0
      DO IREPIJ = 1, NFSYM
         NIJ = NFPCK12(IREPIJ)
         CALL DRV4SYM(WORK,KFREE,LFREE,IPRINT,NSTR(1,0,3),NSTR(1,0,4),
     &                IREPIJ,NIJ,IQMO,WORK(KDXKR12),NIJOFF)
         NIJOFF = NIJOFF + NIJ
      ENDDO
c     do ibsym = 0, nbsym - 1 
c        do iz = 1, nz
c           print*,ibsym,iz,ipqtoq(iz,ibsym),iqmo(iz,ibsym)
c        enddo
c     enddo
C
C     Close MOLFDIR integral file and delete direct access file
C
      IF (.NOT.PARCAL) CALL MDSCRF
      CALL MDINTF
C
C     Release the slaves if we run in parallel.
C
      IF (PARCAL) THEN
         ITASK = -1
#if defined (VAR_MPI)
         CALL MPI_BCAST(ITASK,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMFPACK(INTEGER4,ITASK,1,1,INFO)
#endif
      ENDIF
C
      CALL MEMREL('TRDR4T',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDR4T')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDRV4 */
      SUBROUTINE TRDRV4(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                  NSTR,ANTIS,LMP2,INDX,TRIAN,Q1,Q2,Q3,Q4,E1,E2,
     &                  IC,IQMO)
C
C     Written by Luuk Visscher, december 1996
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxorb.h>
#include <dcbgen.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAL ANTIS, LMP2, NOPV, NODV, TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION IQMO(4,0:7),NZMO(0:7)
      DIMENSION IJP12(0:7)
      DIMENSION IJP34(2)
      DIMENSION IHTM(0:7),IHM(0:7),IFM(0:7),MHTM(2),NHTM(2)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*),E1(*),E2(*)
      DIMENSION IRECC(2),NBFSZ(2),NBUF(2),NIJBUF(2),NHTIM(2),IBOFF(2)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDRV4')
      KFRSAV = KFREE
C
C     We transform all shells in the first halftransformation
C
      CALL ICOPY (24,NSPCK,1,NSPCK12,1)
      CALL ICOPY (128,ISPCK,1,ISPCK12,1)
      CALL ICOPY (24,NBBAS,1,NBBAS1,1)
      CALL ICOPY (24,NBBAS,1,NBBAS2,1)
      CALL ICOPY (24,IBBAS,1,IBBAS1,1)
      CALL ICOPY (24,IBBAS,1,IBBAS2,1)
      CALL ICOPY (2,IBAS,1,IBAS1,1)
      CALL ICOPY (2,IBAS,1,IBAS2,1)
C
C     Copy the first half of the NDMOQ and ICMOQ arrays
C
      CALL ICOPY (8,NDMOQR,1,NDMOQS,1)
      CALL ICOPY (8,NDMOQR,1,NDMOQT,1)
      CALL ICOPY (4,ICMOQR,1,ICMOQS,1)
      CALL ICOPY (4,ICMOQR,1,ICMOQT,1)
C
C     Set up symmetry packing for reduced shell space in the second
C     half transformation. The fourth index is taken completely, the
C     third has only one shell active.
C
      CALL ICOPY (24,NBBAS,1,NBBAS4,1)
      CALL ICOPY (24,IBBAS,1,IBBAS4,1)
      CALL ICOPY (2,IBAS,1,IBAS4,1)
      CALL GISPCK(NSPCK34,ISPCK34,ISHLA,0,NBBAS3,IDUM,IBBAS3,IDUM,
     &            IBAS3,IDUM)
C
C     Generate pointers for reduced coefficient matrices
C
      NQ3S = 0
      NQ4T = 0
      DO IFSYM = 1, NFSYM
         NROW = 0
         DO IREPR = 0, NBSYM-1 
            IF (IFSYM.EQ.JBTOF(IREPR,IC)) THEN
               NROW = NROW + NBBAS3(IREPR,IC)
            ENDIF
         ENDDO
         NDMOQS(1,IFSYM,3) = NROW
         NDMOQS(2,IFSYM,3) = NDMOQR(2,IFSYM,3)
         NDMOQS(1,IFSYM,4) = NDMOQR(1,IFSYM,4)
         NDMOQS(2,IFSYM,4) = NDMOQR(2,IFSYM,4)
         NDMOQT(1,IFSYM,3) = NDMOQR(1,IFSYM,3)
         NDMOQT(2,IFSYM,3) = NDMOQR(2,IFSYM,3)
         NDMOQT(1,IFSYM,4) = NROW
         NDMOQT(2,IFSYM,4) = NDMOQR(2,IFSYM,4)
         ICMOQS(IFSYM,3) = NQ3S + 1
         ICMOQS(IFSYM,4) = ICMOQR(IFSYM,4)
         ICMOQT(IFSYM,3) = ICMOQR(IFSYM,3)
         ICMOQT(IFSYM,4) = NQ4T + 1
         NQ3S = NQ3S + NROW * NDMOQR(2,IFSYM,3) * NZ
         NQ4T = NQ4T + NROW * NDMOQR(2,IFSYM,4) * NZ
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
      NSIZET = 0
      NBUFM = 0
      DO IREPIJ = 1, NFSYM
C
         IF (NHTM(IREPIJ).EQ.0) THEN
            NIJBUF(IREPIJ) = 0 
            NHTIM (IREPIJ) = 0 
            NBUF  (IREPIJ) = 0
            NBFSZ (IREPIJ) = 0
            GOTO 90 
         ENDIF
C
C        Maximum size of the number of half-transformed integrals
C        that we can keep in core. Fix it to 1 million.
C
         NHTIM(IREPIJ) = 1 000 000
         NHTIM(IREPIJ) = MIN0(NHTM(IREPIJ),NHTIM(IREPIJ))
C
C        Calculate the maximum number of ij-pairs that we can handle
C
         NIJBUF(IREPIJ) = NHTIM(IREPIJ) / MHTM(IREPIJ)
C
C        We can now deduce the actual size of NHTIM(IREPIJ)
C
         NHTIM(IREPIJ) = NIJBUF(IREPIJ) * MHTM(IREPIJ)
C
C        Calculate the number of buffers
C
         NBUF(IREPIJ) = NFPCK12(IREPIJ) / NIJBUF(IREPIJ)
         IF (MOD(NFPCK12(IREPIJ),NIJBUF(IREPIJ)).NE.0) THEN
            NBUF(IREPIJ) = NBUF(IREPIJ) + 1
         ENDIF
         NBUFM = MAX0(NBUFM,NBUF(IREPIJ)+1)
C
C        Choose the size of the buffer array : it should be a multiple
C        (named NB) of the number of ij-pairs in one batch multiplied 
C        by NZ
C
         NSIZE = 100 000
         NSIZE = MIN0(NHTM(IREPIJ),NSIZE)
         NB = NSIZE / (NIJBUF(IREPIJ)*NBUF(IREPIJ)*NZ)
         IF (NB.LT.1) NB = 1
         NSIZE = NB * NIJBUF(IREPIJ) * NBUF(IREPIJ) * NZ
C
C        Pointer to the begin of buffer blocks for this irrep
C        All buffers are stored in one 1-dimensional array, but you may
C        think of it as NFSYM sets of 2-dimensional matrices, with the
C        first dimension the size of the buffer and the second the
C        number of buffers for this fermion irrep
C
         IBOFF(IREPIJ) = NSIZET+1
         NSIZET = NSIZET + NSIZE
         NBFSZ(IREPIJ) = NB * NIJBUF(IREPIJ) * NZ
  90     CONTINUE
         IF (IPRINT .GE. 2) WRITE (LUPRI,'(//,4(1X,A,I8,/))')
     &              'Buffer information for parity :',irepij,
     &              'Number of buffers :',nbuf(irepij),
     &              'Pairs in  buffers :',nijbuf(irepij),
     &              'Size   of buffers :',nbfsz(irepij) 
      ENDDO
C
C     Get the memory for the buffers
C
      CALL MEMGET('REAL',KJBOFF,2*NBUFM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRBUF,NSIZET,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLBUF1,NSIZET,WORK,KFREE,LFREE)
C
      CALL INITBUF (IRECC,NBFSZ,NBUF,WORK(KJBOFF))
      CALL MEMCHK('TRDRV4.1',WORK,1)
C
      CALL SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
      IF (I2TYP.LT.0) THEN
         RETURN
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
C     Loop over the shells
C
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
C           First a dummy run to calculate the memory requirements
C           for GMAT. We assume that HMAT is not more than twice GMAT,
C           and reserve the same memory also for use in HERMIT.
C
         NSIZG = LFREE / 3
         NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &                    .FALSE.,INDX,WORK(KDXAB),IJP12,NSIZG,IDUM)
C
         LDXAB = 5*NINSHA*NINSHB
         CALL MEMGET('INTE',KDXAB  ,LDXAB  ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
C        Set INDXAB and INDPASS
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
         ENDIF
C
C        Get distributions
C
         NODV = NASHT.EQ.0
         NOPV = NASHT.LT.2
         CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &               NODV,NOPV,WORK(KFREE),LFREE,IPRINT)
C
C        Do the transformation of the first index pair and sort the
C        half-transformed integrals to file.
C
C        De-allocate buffers and allocate GMAT if necessary.
C
         IF (NPASS.GT.1) THEN
            CALL MEMREL('TRDRV4.1',WORK,1,KALGBUF,KFREE,LFREE)
            CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
            CALL DZERO(WORK(KGMAT),NSIZG)
         ENDIF
C
         CALL MS4IN1D(WORK,KFREE,LFREE,IPRINT,.TRUE.,.FALSE.,
     &                WORK(KALGREC),ICS12,ICF12,IC,
     &                TRIAN(1),NSTR(1,0,1),NSTR(1,0,2),
     &                WORK(KIJPASS),
     &                Q1,Q2,WORK(KGMAT),NSIZG,NZMO,IQMO,
     &                MHTM,IHTM,INDX,WORK(KDXAB),IRECC,WORK(KJBOFF),
     &                IBOFF,NBUF,NIJBUF,NBFSZ,WORK(KRBUF),
     &                WORK(KLBUF1))
C
         CALL MEMREL('TRDRV4.2',WORK,1,KFRSAV3,KFREE,LFREE)
C
C        Delete scalar integral buffer file
C
         IF (NPASS.GT.1) THEN
            CALL DELGBUF(LGFIL)
         ENDIF
C
      ENDDO
      CALL MEMCHK('TRDRV4.2',WORK,1)
C
C     Write out last set of buffers
C     De-allocate buffer space
C
      CALL FINABUF (IRECC,NBFSZ,NBUF,IBOFF,
     &              WORK(KJBOFF),WORK(KRBUF),WORK(KLBUF1))
      CALL MEMREL('TRDRV4.3',WORK,1,KJBOFF,KFREE,LFREE)
C
C     Second half-transformation
C
C
C     Get the reduced coefficient matrices (only active shells)
C
      CALL MEMGET('REAL',KQ3S,NQ3S,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ4T,NQ4T,WORK,KFREE,LFREE)
C
      DO IFSYM = 1, NFSYM
         CALL QREDUC (ISHLA,IFSYM,
     &                Q3(ICMOQR(IFSYM,3)),
     &                NDMOQR(1,IFSYM,3),NDMOQR(2,IFSYM,3),
     &                WORK(KQ3S+ICMOQS(IFSYM,3)-1),
     &                NDMOQS(1,IFSYM,3),NDMOQS(2,IFSYM,3))
c        print*,' q3 '
c        CALL PRQMAT(Q3(ICMOQR(IFSYM,3)),
c    &                NDMOQR(1,IFSYM,3),NDMOQR(2,IFSYM,3),
c    &                NDMOQR(1,IFSYM,3),NDMOQR(2,IFSYM,3),
c    &                NZ,IPQTOQ(1,0),LUPRI)
c        print*,' q3s '
c        CALL PRQMAT(work(kQ3S+ICMOQS(IFSYM,3)-1),
c    &                NDMOQS(1,IFSYM,3),NDMOQS(2,IFSYM,3),
c    &                NDMOQS(1,IFSYM,3),NDMOQS(2,IFSYM,3),
c    &                NZ,IPQTOQ(1,0),LUPRI)
         CALL QREDUC (ISHLA,IFSYM,
     &                Q4(ICMOQR(IFSYM,4)),
     &                NDMOQR(1,IFSYM,4),NDMOQR(2,IFSYM,4),
     &                WORK(KQ4T+ICMOQT(IFSYM,4)-1),
     &                NDMOQT(1,IFSYM,4),NDMOQT(2,IFSYM,4))
      ENDDO
C
      CALL MEMCHK('TRDRV4.3',WORK,1)
C
      ICS34 = IC
      ICF34 = IC
C
      CALL MEMGET('INTE',KDXKR12,NFPCK12T*2,WORK,KFREE,LFREE)
      CALL MKINDXKR (NSTR,TRIAN(1),WORK(KDXKR12))
C
      CALL MS4IN2D(WORK,KFREE,LFREE,IPRINT,ICS34,ICF34,
     &             NSTR(1,0,3),NSTR(1,0,4),
     &             WORK(KQ3S),Q4,Q3,WORK(KQ4T),NZMO,
     &             IQMO,IRECC,NBFSZ,NBUF,NIJBUF,MHTM,WORK(KDXKR12))
C
C     Delete half-transformed integral files and intermediary file
C
      CALL DELEBUF
C
C     Release all memory and exit
C
      CALL MEMREL('TRDRV4.5',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDRV4')
C
      RETURN
C
 1000 FORMAT (/' Storing all half-transformed integrals to disk',
     & /' Disk requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DRV4SYM*/
      SUBROUTINE DRV4SYM(WORK,KFREE,LFREE,IPRINT,NSTR3,NSTR4,
     & IREPIJ,NIJ,IQMO,INDXKR,NIJOFF)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               Symmetrization of integrals.
C
C     Input :
C
C     - IPRINT       Print flag
C     - NSTR3        Number of active spinors for index 3
C     - NSTR4        Number of active spinors for index 4
C     - IREPIJ       Fermion symmetry (parity)
C     - NIJ          Number of blocks for this fermion symmetry
C     - Q3           Coefficients for index 3
C     - Q4           Coefficients for index 4
C
C     Output is written directly in the MOLFDIR-type file MDCINT
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
#include <dcbgen.h>
#include <infpar.h>
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION INDXKR(2,*)
      DIMENSION IQMO(4,0:7)
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbtr3.h>
C
      KFRSAV = KFREE
C
C     Number of boson symmetry for each fermion irrep (parity)
C
      NBSYMP = NBSYM/NFSYM
C
C     Size of integral batches.
C
      NKL = NFPCK34(IREPIJ)
      CALL MEMGET('REAL',KHMAT,NKL*NZ*NZ*NBSYMP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHMATS,NKL*NZ*NZ*NBSYMP,WORK,KFREE,LFREE)
      CALL DZERO (WORK(KHMAT),NKL*NZ*NZ*NBSYMP)
C
      IF (PARCAL) THEN
        DO IJ = 1, NIJ
          IKR = INDXKR(1,IJ+NIJOFF)
          JKR = INDXKR(2,IJ+NIJOFF)
          IF (MYTID.NE.MPARID)
     &    CALL MDSCRR(IPRINT,IREPIJ,IJ,NKL,WORK(KHMAT),WORK,KFREE,LFREE)
#if defined (VAR_MPI)
          CALL MPI_REDUCE(WORK(KHMAT),WORK(KHMATS),NKL*NZ*NZ*NBSYMP,
     &                    MPI_DOUBLE_PRECISION,MPI_SUM,MPARID,
     &                    MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
          write pvm code
#endif
          IF (MYTID.EQ.PARID) THEN
             CALL SYMFINT(IPRINT,IREPIJ,IJ,IKR,JKR,NSTR3,NSTR4,NKL,
     &                     WORK(KHMATS),IPQTOQ,WORK,KFREE,LFREE)
ctest&                     WORK(KHMATS),IQMO,WORK,KFREE,LFREE)
          ENDIF
        ENDDO
      ELSE
        DO IJ = 1, NIJ
          IKR = INDXKR(1,IJ+NIJOFF)
          JKR = INDXKR(2,IJ+NIJOFF)
          CALL MDSCRR(IPRINT,IREPIJ,IJ,NKL,WORK(KHMAT),WORK,KFREE,LFREE)
          CALL SYMFINT(IPRINT,IREPIJ,IJ,IKR,JKR,NSTR3,NSTR4,NKL,
     &                  WORK(KHMAT),IPQTOQ,WORK,KFREE,LFREE)
ctest&                  WORK(KHMAT),IQMO,WORK,KFREE,LFREE)
        ENDDO
      ENDIF
      CALL MEMREL('DRV4SYM',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C Start of parallel section : activitate compilation only when necessary
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (VAR_PVM) || defined (VAR_MPI)
C/* Deck TRAPARI*/
      SUBROUTINE TRAPARI(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                   NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,KINDX,
     &                   KQ,KE)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher (based on HERNOD) August 1997
C
C    *******************************************************************
C    *   Initializes the nodes for the parallel 4-index transformation *
C    *******************************************************************
C
#include <implicit.h>
C
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
C
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <dcbgen.h>
C
#include <dorps.h>
#include <dcborb.h>
#include <infpar.h>
C
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION KQ(4),KE(2)
      DIMENSION WORK(*),DINTSKP(2,4)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH
C
      CALL QENTER('TRAPARI')
C
      KFRSAV = KFREE
      PARHER = .TRUE.
      PARCAL = .TRUE.
C
      FINISH = .FALSE.
C
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
C
#endif
      IF (MYTID.EQ.0) THEN
C
C        **********************************************************
C        ***** Send HERMIT initialization data to the other nodes *
C        **********************************************************
C
C----------------------------------------------------------------
C        Setup information for the two-electron integralroutines.
C----------------------------------------------------------------
C
         NODV = .FALSE.
         NOPV = .FALSE.
         IF (NASHT .EQ. 0) NODV = .TRUE.
         IF (NASHT .LT. 2) NOPV = .TRUE.
         IPRALL =  0
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
         NOCONT = .FALSE.
C
         MAXDIF = 0
         ITYPE  = 4
C
         CALL SDINIT(DUMMY,0,IDUMMY,IDUMMY,ITYPE,MAXDIF,IDUMMY,NODV,
     &               NOPV,NOCONT,TKTIME,RTNTWO,I2TYP,IDUMMY,IDUMMY,
     &               DUMMY,DUMMY,DUMMY)
C
      ELSE 
C
C        **********************************************************
C        ***** Receive HERMIT initialization data from master *****
C        **********************************************************
C
         CALL RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &               IATOM,MAXDIF,MYTID,NODV,NOPV,NOCONT,RETUR,
     &               TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,IFTHRS,
     &               KFMAT,KDMAT,KIFC,KIRD,KDMRAO,KGABAO,KDMRSO)
C
C     We call MEMREL because we do not need any of the variable size
C     arrays in the 4-index transformation (these are used in the
C     parallel HERMIT driver)
C
         CALL MEMREL('TRAPARI',WORK,1,KFRSAV,KFREE,LFREE)
C
      ENDIF
C
C        *****************************************************
C        ***** Send or receive DIRAC initialization data *****
C        *****************************************************
C
      CALL TRAPRD(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &            NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,KINDX,KQ,KE)
C
C     Do not call MEMREL since the variable size arrays are represented
C     as pointers when the routine is called from a slave
C
      CALL QEXIT('TRAPARI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck traprd */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TRAPRD (WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                   NSTR,ANTIS,LMP2,INDX,Q1,Q2,Q3,Q4,E1,E2,
     &                   KINDX,KQ,KE)
C
C     Communicate DIRAC common blocks for parallel 4-index transform.
C     Communicate arguments received from main program
C
C     Common blocks initialized : dcborb, dcbbas, dgroup, dcbtri,
C                                 dcbtr3, comdis and symmet
C
C     Written by Luuk Visscher, august 1997.
C
#include <implicit.h>
C
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <dcborb.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <symmet.h>
#include <dgroup.h>
#include <infpar.h>
#include <parint.h>
C
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION Q1(*),Q2(*),Q3(*),Q4(*)
      DIMENSION E1(*),E2(*)
      DIMENSION KQ(4),KE(2)
      DIMENSION WORK(*)
C
      DIMENSION NQSZ(4),NESZ(2)
C
#if defined (VAR_PVM)
      CALL PVMFMYTID(MYTID)
      IF (MYTID .LT. 0) CALL PFAIL(MYTID)
C
      CALL PVMFPARENT(MPARID)
      IF (MPARID .LT. 0) CALL PFAIL(MPARID)
C
#endif
#if defined (VAR_MPI)
      NUMELM = NDCBORB
      CALL MPI_BCAST(NORB,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NDCBBAS
      CALL MPI_BCAST(NBBAS,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NDGROUP
      CALL MPI_BCAST(NZ,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NDCBTRI
      CALL MPI_BCAST(IPRTRA,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NCOMDIS
      CALL MPI_BCAST(NSPCK,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTI
      CALL MPI_BCAST(MAXREP,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = NSYMMTR
      CALL MPI_BCAST(FMULT,NUMELM,MPI_DOUBLE_PRECISION,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 1
      CALL MPI_BCAST(INTFLG,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 1
      CALL MPI_BCAST(ANTIS,NUMELM,MPI_LOGICAL,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 1
      CALL MPI_BCAST(LMP2,NUMELM,MPI_LOGICAL,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 24
      CALL MPI_BCAST(NSTR,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 16
      CALL MPI_BCAST(NDMOQR,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
      NUMELM = 8
      CALL MPI_BCAST(ICMOQR,NUMELM,MPI_INTEGER,MPARID,
     &               MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
      NUMELM = NDCBORB
      CALL PVMFPACK(INTEGER4,NORB,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NDCBBAS
      CALL PVMFPACK(INTEGER4,NBBAS,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NDGROUP
      CALL PVMFPACK(INTEGER4,NZ,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NDCBTRI
      CALL PVMFPACK(INTEGER4,IPRTRA,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NCOMDIS
      CALL PVMFPACK(INTEGER4,NSPCK,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NSYMMTI
      CALL PVMFPACK(INTEGER4,MAXREP,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = NSYMMTR
      CALL PVMFPACK(REAL8,FMULT,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 1
      CALL PVMFPACK(INTEGER4,INTFLG,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 1
      CALL PVMFPACK(LOGICAL,ANTIS,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 1
      CALL PVMFPACK(LOGICAL,LMP2,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 24
      CALL PVMFPACK(INTEGER4,NSTR,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 16
      CALL PVMFPACK(INTEGER4,NDMOQR,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
      NUMELM = 8
      CALL PVMFPACK(INTEGER4,ICMOQR,NUMELM,1,INFO)
      IF (INFO .LT. 0) CALL PFAIL(INFO)
#endif
C
C     Calculate size of variable size arrays, when called from
C     the master we send out the data, when called from a slave
C     we first allocate the array and then receive. The pointers
C     will be sent up.
C
      IF (MYTID.NE.MPARID) THEN
         CALL MEMGET('INTE',KINDX,3*NTBAS(0),WORK,KFREE,LFREE)
      ENDIF
C
      DO I = 1, 4
         NQSZ(I) = 0
         DO J = 1, NFSYM
            NQSZ(I) = NQSZ(I) + NDMOQR(1,J,I)*NDMOQR(2,J,I)*NZ
         ENDDO
         IF (MYTID.NE.MPARID) THEN
            CALL MEMGET('REAL',KQ(I),NQSZ(I),WORK,KFREE,LFREE)
         ENDIF
      ENDDO
C
      DO I = 1, 2
         NESZ(I) = 0
         DO J = 1, NFSYM
            NESZ(I) = NESZ(I) + NSTR(J,0,I)
         ENDDO
         IF (MYTID.NE.MPARID) THEN
            CALL MEMGET('REAL',KE(I),NESZ(I),WORK,KFREE,LFREE)
         ENDIF
      ENDDO
C
#if defined (VAR_MPI)
      IF (MYTID.EQ.MPARID) THEN
         CALL MPI_BCAST(INDX,NTBAS(0)*3,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(Q1,NQSZ(1),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(Q2,NQSZ(2),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(Q3,NQSZ(3),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(Q4,NQSZ(4),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(E1,NESZ(1),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(E2,NESZ(2),MPI_DOUBLE_PRECISION,MPARID,
     &                     MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_BCAST(WORK(KINDX),NTBAS(0)*3,MPI_INTEGER,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KQ(1)),NQSZ(1),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KQ(2)),NQSZ(2),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KQ(3)),NQSZ(3),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KQ(4)),NQSZ(4),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KE(1)),NESZ(1),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(WORK(KE(2)),NESZ(2),MPI_DOUBLE_PRECISION,MPARID,
     &                  MPI_COMM_WORLD,IERR)
      ENDIF
#endif
#if defined (VAR_PVM)
      write pvm code for this part
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRANOD*/
      SUBROUTINE TRANOD(WORK,LWORK)
C
C     Written by Luuk Visscher, August 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#else
#include <fpvm3h.h>
#endif
#include <infpar.h>
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbbas.h>
#include <dgroup.h>
      LOGICAL ANTIS,LMP2,TRIAN(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION IQMO(4,0:7)
      DIMENSION KQ(4),KE(2)
      DIMENSION WORK(*)
C
      CALL QENTER('TRANOD')
#include <memint.h>
      KFRSAV = KFREE
C
C     Initialize the common blocks.
C
      CALL TRAPARI(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &             NSTR,ANTIS,LMP2,IDUM,DUM,DUM,DUM,DUM,DUM,DUM,KINDX,
     &             KQ,KE)
C
      TRIAN(1) = ISAME(2).EQ.ISAME(1)
      TRIAN(2) = .FALSE.
C
C     Define the packing of the 2-index transformed integrals
C     This should be modified so that we ony pack them for a specific
C     value of ISHLA (define NSPCK inside the loop and split to NSPCK12
C     and NSPCK34). To use this packing we
C     also need to modify the second halftransformation and we need to
C     loop over squares of shells rather than triangles. This cding is
C     postponed for the moment. The code now becomes n6 in the second
C     part (with the extra n the number of shells).
C
      NFPCK12T = 0
      NFPCK34T = 0
      DO IREPIJ =1, NFSYM
         NFPCK12(IREPIJ) = 0
         NFPCK34(IREPIJ) = 0
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
C
C           Packing can be triangular (i=>j) or square.
C
            IF (.NOT.TRIAN(1).OR.(TRIAN(1).AND.IREPI.GT.IREPJ)) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ) = NFPCK12(IREPIJ) +
     &                           NSTR(IREPI,0,1)*NSTR(IREPJ,0,2)
            ELSEIF (TRIAN(1).AND.IREPI.EQ.IREPJ) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ) = NFPCK12(IREPIJ) +
     &                           NSTR(IREPI,0,1)*(NSTR(IREPI,0,1)+1)/2
            ELSE
               IFPCK12(IREPI,IREPJ) = - 1
            ENDIF
            IF (.NOT.TRIAN(2).OR.(TRIAN(2).AND.IREPI.GT.IREPJ)) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ) = NFPCK34(IREPIJ) +
     &                           NSTR(IREPI,0,3)*NSTR(IREPJ,0,4)
            ELSEIF (TRIAN(2).AND.IREPI.EQ.IREPJ) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ) = NFPCK34(IREPIJ) +
     &                           NSTR(IREPI,0,3)*(NSTR(IREPI,0,3)+1)/2
            ELSE
               IFPCK34(IREPI,IREPJ) = - 1
            ENDIF
         ENDDO
         NFPCK12T = NFPCK12T + NFPCK12(IREPIJ)
         NFPCK34T = NFPCK34T + NFPCK34(IREPIJ)
      ENDDO
C
C     Initialize output file for final integrals
C
      CALL MDSCRI(WORK,KFREE,LFREE,NFPCK12,NFPCK34,MYTID)
C
C     Wait for tasks until we receive the end message
C
  100 CONTINUE
#if defined (VAR_MPI)
         CALL MPI_SEND(MYTID,1,MPI_INTEGER,MPARID,20,
     &                 MPI_COMM_WORLD,IERR)
         CALL MPI_RECV(ISHLA,1,MPI_INTEGER,MPARID,30,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
         CALL MPI_RECV(IC,1,MPI_INTEGER,MPARID,30,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
#else
#endif
         IF (ISHLA.EQ.-1) GOTO 200
         IPRINT = -1
         CALL TRDRV4(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &               NSTR(1,0,1),ANTIS,LMP2,WORK(KINDX),TRIAN,
     &               WORK(KQ(1)),WORK(KQ(2)),WORK(KQ(3)),WORK(KQ(4)),
     &               WORK(KE(1)),WORK(KE(2)),IC,IQMO)
         GOTO 100
  200 CONTINUE
C
C     Complete the integrals to the desired classes
C
      CALL MEMGET('INTE',KDXKR12,NFPCK12T*2,WORK,KFREE,LFREE)
      CALL MKINDXKR (NSTR,TRIAN(1),WORK(KDXKR12))
C
      NIJOFF = 0
      DO IREPIJ = 1, NFSYM
         NIJ = NFPCK12(IREPIJ)
         CALL DRV4SYM(WORK,KFREE,LFREE,IPRINT,NSTR(1,0,3),NSTR(1,0,4),
     &                IREPIJ,NIJ,IQMO,WORK(KDXKR12),NIJOFF)
         NIJOFF = NIJOFF + NIJ
      ENDDO
C
C     Delete direct access file
C
      CALL MDSCRF
C
C     We check whether all the work is done and return to the 
C     main slave menu.
C
#if defined (VAR_MPI)
      CALL MPI_BCAST(NTEST,1,MPI_INTEGER,MPARID,MPI_COMM_WORLD,IERR)
#else
#endif
      IF (NTEST.NE.-1) CALL QUIT('Received wrong message')
C
      CALL MEMREL('TRANOD',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRANOD')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C& End of parallel section : activitate compilation only when necessary
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qreduc */
      SUBROUTINE QREDUC (ISHELA1,IFSYM,QM,NDIM1,NDIM2,
     &                                 QMR,NDIMR1,NDIMR2)
C
C     Takes out active rows of coefficient matrix for a given irrep.
C     NDIM and NDIMR is now NDIM1,NDIM2,NDIMR1 and NDIMR2
C
C     Luuk Visscher, 25-8-97
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <aovec.h>
#include <blocks.h>
#include <dcbtra.h>
#include <dgroup.h>
#include <dcbbas.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
C
      DIMENSION QM(NDIM1,NDIM2,NZ),QMR(NDIMR1,NDIMR2,NZ)
C
#include <ibtfun.h>
C
      ICA    = LCLASH(ISHELA1)
      NHKTA1 = NHKTSH(ISHELA1)
      KHKTA1 = KHKTSH(ISHELA1)
      MULA1  = ISTBSH(ISHELA1)
      NORBA1 = NORBSH(ISHELA1)
      NSTRA1 = IORBSB(IORBSH(ISHELA1,1))
C
      INDA2 = 0
      DO IREPA1 = 0, MAXREP
C
C        Loop over components
C
         DO NA = 1,KHKTA1
            NSTRNA = NSTRA1 + NA
            ITYNA  = ISYMAO(NHKTA1,NA)
            IF ( JBTOF(IREPA1,ICA).EQ.IFSYM .AND.
     &           IBTAND(MULA1,IBTXOR(IREPA1,ITYNA)) .EQ. 0) THEN
C
C              Loop over contracted orbitals
C
               DO IA = 1,NORBA1
                  INDA = IPTSYM(NSTRNA + KHKTA1*(IA-1),IREPA1)
                  INDA1 = INDA - ICOS(IREPA1+1,ICA) 
     &                  + IBBAS(IREPA1,ICA) - IBAS(IFSYM)
                  INDA2 = INDA2 + 1
C
C                 This function is active : copy the row to the reduced
C                 matrix
C
                  DO IZ = 1, NZ
                     CALL DCOPY (NDIMR2,QM(INDA1,1,IZ),NDIM1,
     &                           QMR(INDA2,1,IZ),NDIMR1)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
