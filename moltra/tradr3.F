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
C FILE    : tradr3.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR3T */
      SUBROUTINE TRDR3T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                  NSTR,ANTIS,LMP2,INDX,KQ,KE,KIBE,DINTSKP)
C
C     Written by Luuk Visscher, december 1996
C     Screening introduced by T.Saue Sep 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
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
      DIMENSION KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION IJP12(0:7)
      DIMENSION IJP34(2)
      DIMENSION IHTM(0:7),IHM(0:7),IFM(0:7),MHTM(2),NHTM(2)
      DIMENSION IRECC(2),NBFSZ(2),NBUF(2),NIJBUF(2),NHTIM(2),IBOFF(2)
      DIMENSION WORK(*)
C
      CALL QENTER('TRDRV3')
      KFRSAV = KFREE
C
C     Define the packing of the 2-index transformed integrals
C
      TRIAN(1) = ISAME(2).EQ.ISAME(1)
      TRIAN(2) = .FALSE.
      CALL PCK2IN(NSTR,TRIAN,IPRINT)
C
      CALL ICOPY (24,NSPCK,1,NSPCK12,1)
      CALL ICOPY (24,NSPCK,1,NSPCK34,1)
      CALL ICOPY (128,ISPCK,1,ISPCK12,1)
      CALL ICOPY (128,ISPCK,1,ISPCK34,1)
      CALL ICOPY (24,NBBAS,1,NBBAS1,1)
      CALL ICOPY (24,NBBAS,1,NBBAS2,1)
      CALL ICOPY (24,NBBAS,1,NBBAS3,1)
      CALL ICOPY (24,NBBAS,1,NBBAS4,1)
      CALL ICOPY (24,IBBAS,1,IBBAS1,1)
      CALL ICOPY (24,IBBAS,1,IBBAS2,1)
      CALL ICOPY (24,IBBAS,1,IBBAS3,1)
      CALL ICOPY (24,IBBAS,1,IBBAS4,1)
      CALL ICOPY (2,IBAS,1,IBAS1,1)
      CALL ICOPY (2,IBAS,1,IBAS2,1)
      CALL ICOPY (2,IBAS,1,IBAS3,1)
      CALL ICOPY (2,IBAS,1,IBAS4,1)
      CALL ICOPY (16,NDMOQR,1,NDMOQS,1)
      CALL ICOPY (16,NDMOQR,1,NDMOQT,1)
      CALL ICOPY (8,ICMOQR,1,ICMOQS,1)
      CALL ICOPY (8,ICMOQR,1,ICMOQT,1)
C
      CALL MEMGET('INTE',KDXKR12,NFPCK12T*2,WORK,KFREE,LFREE)
      CALL MKINDXKR (NSTR,TRIAN(1),WORK(KDXKR12))
      CALL MEMGET('INTE',KDXB12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MKINDXB (NSTR,WORK(KIBE(1,1)),WORK(KIBE(2,1)),
     &              WORK(KIBE(1,2)),WORK(KIBE(2,2)),
     &              TRIAN(1),WORK(KDXB12))
      CALL MEMGET('INTE',KDXB34,2*NFPCK34T,WORK,KFREE,LFREE)
      CALL MKINDXB (NSTR,WORK(KIBE(1,3)),WORK(KIBE(2,3)),
     &              WORK(KIBE(1,4)),WORK(KIBE(2,4)),
     &              TRIAN(2),WORK(KDXB34))
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
               MHTM(IREPIJ) = MHTM(IREPIJ) + NSPCK(IREPAB,0)*NZ
            ENDIF
            NHTM(IREPIJ) = MHTM(IREPIJ)*NFPCK12(IREPIJ)
         ENDDO
         NHTMT = NHTMT + NHTM(IREPIJ)
      ENDDO
C
C     Calculate disk space for half transformed integrals
C
      IF (IPRINT .GE. 0) WRITE (LUPRI,1000) REAL(NHTMT)*16./(1024.**2)
C
      NSIZET = 0
      NBUFM = 0
      DO IREPIJ = 1, NFSYM
C
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
C
C     Prepare for screening
C
      IF(SCRTRA.GT.D0) THEN
        CALL PR4SCR(KGAB,KDRIJ,WORK,KFREE,LFREE,DINTSKP,
     &              WORK(KQ(1,1)),WORK(KQ(1,2)),
     &              WORK(KQ(1,3)),WORK(KQ(1,4)),
     &              NDMOQR,ICMOQR,NSTR,ISAME,IPRINT)
      ELSE
        KDRIJ = KFREE
        KGAB  = KFREE
      ENDIF
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
C           First a dummy run to calculate the memory requirements
C           for GMAT. We assume that HMAT is not more than twice GMAT,
C           and reserve the same memory also for use in HERMIT.
C
            NSIZG = LFREE / 3
            LDXAB = 5*NINSHA*NINSHB
            CALL MEMGET('INTE',KDXAB  ,LDXAB  ,WORK,KFREE,LFREE)
            NUMDIS = NDISTRN(I2TYP,IPRINT,.FALSE.,.FALSE.,
     &                       .FALSE.,INDX,WORK(KDXAB),
     &                       IJP12,NSIZG,IDUM)
C
C
            CALL MEMGET('INTE',KIJPASS,8*NPASS,WORK,KFREE,LFREE)
C
C           Set INDXAB and INDPASS
C
            NUMDIS = NDISTRN(I2TYP,IPRINT,.TRUE.,.TRUE.,.FALSE.,
     &               INDX,WORK(KDXAB),IJP12,NSIZG,WORK(KIJPASS))
C
            IF (IPRINT .GE. 3) THEN
               WRITE (LUPRI,'(//,4(1X,A,2I4,/))')
     &         ' Calculation of integral distribution set:',ISHLA,ISHLB
               WRITE (LUPRI,'(//,4(1X,A,I8,/))')
     &         ' Number of distributions in this set:     ',NUMDIS,
     &         ' Maximum size of scalar integral batch:   ',NSIZG,
     &         ' Number of passes :                       ',NPASS
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
            CALL CALDIS(I2TYP,WORK(KGMAT),INDX,WORK(KDXAB),
     &               NODV,NOPV,WORK(KGAB),WORK(KDRIJ),DINTSKP,SCRTRA,
     &               WORK(KFREE),LFREE,IPRINT)
C
C           Do the transformation of the first index pair and sort the
C           half-transformed integrals to file.
C
C           De-allocate buffers and allocate GMAT if necessary.
C
            IF (NPASS.GT.1) THEN
               CALL MEMREL('TRDRV3.1',WORK,1,KALGBUF,KFREE,LFREE)
               CALL MEMGET('REAL',KGMAT,NSIZG,WORK,KFREE,LFREE)
               CALL DZERO(WORK(KGMAT),NSIZG)
            ENDIF
C
            CALL MS4IN1D(WORK,KFREE,LFREE,IPRINT,.FALSE.,.FALSE.,
     &                  WORK(KALGREC),ICS,ICF,
     &                  IC,TRIAN(1),NSTR(1,0,1),NSTR(1,0,2),
     &                  WORK(KIJPASS),WORK(KQ(1,1)),WORK(KQ(1,2)),
     &                  WORK(KGMAT),NSIZG,MHTM,IHTM,
     &                  INDX,WORK(KDXAB),IRECC,WORK(KJBOFF),
     &                  IBOFF,NBUF,NIJBUF,NBFSZ,WORK(KRBUF),
     &                  WORK(KLBUF1))
C
            CALL MEMREL('TRDRV3.2',WORK,1,KFRSAV3,KFREE,LFREE)
C
C           Delete scalar integral buffer file
C
            IF (NPASS.GT.1) THEN
               CALL DELGBUF(LGFIL)
            ENDIF
C
         ENDDO
         ENDDO
C
C----------------
 100     CONTINUE
C----------------
C
C     Screening statistics
C
      IF(SCRTRA.GT.0) THEN
        CALL ST4SCR(DINTSKP)
      ENDIF
C
C     Write out last set of buffers
C     De-allocate buffer space
C
      CALL FINABUF (IRECC,NBFSZ,NBUF,IBOFF,
     &              WORK(KJBOFF),WORK(KRBUF),WORK(KLBUF1))
      CALL MEMREL('TRDRV3.3',WORK,1,KJBOFF,KFREE,LFREE)
C
C     Initialize MOLFDIR integral file
C
      CALL MEMGET('INTE',KIKRA,NORBT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIKRB,NORBT,WORK,KFREE,LFREE)
      CALL MDINTI(NSTR,WORK(KIKRA),WORK(KIKRB),IPRINT)
      CALL MEMREL('TRDRV3.4',WORK,1,KIKRA,KFREE,LFREE)
C
C     Second half-transformation
C
      ICS = 1
      ICF = 2
      CALL MS4IN2D(WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &             NSTR(1,0,3),NSTR(1,0,4),
     &             WORK(KQ(1,3)),WORK(KQ(1,4)),
     &             WORK(KQ(1,3)),WORK(KQ(1,4)),
     &             IRECC,NBFSZ,NBUF,NIJBUF,MHTM,
     &             WORK(KDXKR12),WORK(KDXB12),WORK(KDXB34))
C
C     Delete half-transformed integral files
C
      CALL DELEBUF
C
C     Close integral file and delete intermediary file
C
      CALL MDINTF
C
C     Release all memory and exit
C
      CALL MEMREL('TRDRV3.5',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDRV3')
C
      RETURN
C
 1000 FORMAT (/' Storing all half-transformed integrals to disk',
     & /' Disk requirements ',F10.3,' Megabytes')
 1001 FORMAT (/' Storing all fully-transformed integrals to memory',
     & /' Memory requirements ',F10.3,' Megabytes')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MS4INW*/
      SUBROUTINE MS4INW (WORK,KFREE,LFREE,IPRINT,
     &                   SUB3,SUB4,IC,IHM,MHTM,IHTM,INDX,INDXAB,IPASS,
     &                   HMAT,IRECC,JBOFF,IBOFF,NBUF,NIJBUF,NBFSZ,
     &                   RBUF,LBUF1)

C
C     Sort half-transformed integrals
C     This version writes them to disk
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
#include <dcbtr3.h>
#include <dcbbas.h>
C
      PARAMETER (HALF=0.5D0)
C
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5)
      DIMENSION IHTM(0:7),IHM(0:7),MHTM(2)
      DIMENSION HMAT(*)
      DIMENSION WORK(*)
      DIMENSION IND(2)
      LOGICAL SUB3,SUB4
C
C     For the buffered output
C
      DIMENSION IRECC(2),IBOFF(2),NBFSZ(2),NBUF(2),NIJBUF(2)
      DIMENSION JBOFF(2,0:*),RBUF(*),LBUF1(*)
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
            IF (IPASS.NE.INDXAB(INDAR,INDBR,5)) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
            INDA1 = INDA - ICOS(IREPA+1,IC) ! Absolute position
            INDB1 = INDB - ICOS(IREPB+1,IC) ! Absolute position
            IF (SUB3) INDA1 = INDX(3,INDA)  ! Relative position
            IF (SUB4) INDB1 = INDX(3,INDB)  ! Relative position
C
C           The offset within the untransformed indexes. These were
C           the right-hand side (slowest varying) and are sorted to
C           become the left hand side. Symmetrization is done
C           in the second half-transformation (the function indices are
C           left triangular), but we scale the diagonal elements to 
C           avoid double counting.
C
            JOFFAB = ISPCK34(IREPA,IREPB,IC)
     &             + (INDB1-1)*NBBAS3(IREPA,IC)+INDA1
            JOFF1 = IHTM(IREPAB) + JOFFAB
            IF (INDA.EQ.INDB) 
     &         CALL DSCAL (NZ*NFPCK12(IREPIJ),HALF,HMAT(IOFF),1)
C
C           Number of pairs in one buffer
C
            NIJ  = NIJBUF(IREPIJ)
            NIJL = NFPCK12(IREPIJ) - NIJ*(NBUF(IREPIJ)-1)
            IF(NIJL.EQ.0) GOTO 10
C
            DO IZ = 1, NZ
               IJOFF = 0
C
C              Distribute the integrals over the output buffers
C              The buffers are kept in 1-dimensional arrays with
C              the offsets stored in JBOFF. This allows a variable
C              buffer length for the different irreps. The amount of
C              integrals in each buffer will be the same because we
C              do not check for zeroes.
C              The number of integrals in the buffers is JBOFF(IREPIJ,0)
C
               DO IBUF = 1, NBUF(IREPIJ)-1
                  JBUF = JBOFF(IREPIJ,IBUF)
                  CALL DCOPY(NIJ,HMAT(IOFF+IJOFF),1,RBUF(JBUF),1)
C
C                 Write the absolute addresses of the integrals
C
                  DO IJ = 1, NIJ
                     LBUF1(JBUF) = JOFF1+(IJOFF+IJ-1)*MHTM(IREPIJ)
                     JBUF = JBUF + 1
                  ENDDO
                  JBOFF(IREPIJ,IBUF) = JBUF
                  IJOFF = IJOFF + NIJ
               ENDDO
               JBOFF(IREPIJ,0) = JBOFF(IREPIJ,0) + NIJ
C
C              Fill the last buffer : the number of integrals may be
C              different
C
               IBUF = NBUF(IREPIJ)
               JBUF = JBOFF(IREPIJ,IBUF)
               CALL DCOPY(NIJL,HMAT(IOFF+IJOFF),1,RBUF(JBUF),1)
               DO IJ = 1, NIJL
                  LBUF1(JBUF) = JOFF1+(IJOFF+IJ-1)*MHTM(IREPIJ)
                  JBUF = JBUF + 1
               ENDDO
               JBOFF(IREPIJ,IBUF) = JBUF
C
               IOFF  = IOFF  + NFPCK12(IREPIJ)
               JOFF1 = JOFF1 + NSPCK34(IREPAB,0)
            ENDDO
C
C           When the buffers are full we write them to disk
C           In this case JBUF has become equal to the length
C           of the buffer array for this irrep.
C
            IF (JBOFF(IREPIJ,0).GT.NBFSZ(IREPIJ)) THEN
               CALL QUIT ('Buffer overflow for half-transformed ints')
            ELSEIF (JBOFF(IREPIJ,0).EQ.NBFSZ(IREPIJ)) THEN
C
C              Now write the buffers : write always the full length,
C              but write also the actual length in the first spot
C
               IREC = IRECC(IREPIJ)
               JBUFF = IBOFF(IREPIJ)
               DO IBUF = 1, NBUF(IREPIJ)
                  NBUFT = JBOFF(IREPIJ,IBUF) - JBUFF
                  JBUFL = JBUFF + NBFSZ(IREPIJ) - 1
                  WRITE (LUTRA1+IREPIJ,REC=IREC) NBUFT,
     &                  (RBUF(JBUF),JBUF=JBUFF,JBUFL),
     &                  (LBUF1(JBUF),JBUF=JBUFF,JBUFL) 
C
C                 Reset offsets 
C
                  JBOFF(IREPIJ,IBUF) = JBUFF
                  JBUFF = JBUFL + 1
                  IREC = IREC + 1
               ENDDO
               JBOFF(IREPIJ,0) = 0
               IRECC(IREPIJ) = IREC
            ENDIF
  10     CONTINUE
  20  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN2D*/
      SUBROUTINE MS4IN2D(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,
     & Q3S,Q4S,Q3T,Q4T,IRECC,NBFSZ,NBUF,NIJBUF,MHTM,INDXKR,
     & INDXB12,INDXB34)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C
C     Driver for the second index transformation
C     Get part of the half-transformed integrals from memory, call
C     MS4IN2
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - NSTR3        Number of active spinors for index 3
C     - NSTR4        Number of active spinors for index 4
C     - Q3S          Coefficients for index 3
C     - Q4S          Coefficients for index 4
C     - Q3T          Coefficients for index 3 but for transposed AO ints
C     - Q4T          Coefficients for index 4 but for transposed AO ints
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      DIMENSION IRECC(2),MHTM(2),NBUF(2),NIJBUF(2),NBFSZ(2)
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN2D')
C
C     The integrals are stored with their absolute address, IOFF is 
C     the start of the current batch. The relative address is obtained 
C     by subtracting IOFF.
C
C     NIJOFF is the offset for INDXKR (one list for all irreps)
C     IJS is the offset for the start of the batch within an irrep
C
      IOFF = 0
      NIJOFF = 1
      NKLOFF = 1
      DO IREPIJ = 1, NFSYM
         IJS = 0
         NIJ  = NIJBUF(IREPIJ)
         NIJL = NFPCK12(IREPIJ) - NIJ*(NBUF(IREPIJ)-1)
         NKL  = NFPCK34(IREPIJ)
         CALL MEMGET('REAL',KHMAT,NIJ*MHTM(IREPIJ),WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KRBUF,NBFSZ(IREPIJ)*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KLBUF1,NBFSZ(IREPIJ)*NZ,WORK,KFREE,LFREE)
         DO IBUF = 1, NBUF(IREPIJ)
            IF (IBUF.EQ.NBUF(IREPIJ)) NIJ = NIJL
            CALL DZERO(WORK(KHMAT),NIJ*MHTM(IREPIJ))
            CALL READBUF(IREPIJ,IBUF,IOFF,IRECC(IREPIJ),NBFSZ(IREPIJ),
     &                   NBUF(IREPIJ),WORK(KRBUF),WORK(KLBUF1),
     &                   WORK(KHMAT))
            CALL MS4IN2P(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,
     &                   IREPIJ,IJS,NIJ,Q3S,Q4S,Q3T,Q4T,
     &                   WORK(KHMAT),INDXKR(1,NIJOFF),
     &                   INDXB12(1,NIJOFF),INDXB34(1,NKLOFF))
            IOFF = IOFF + NIJ * MHTM(IREPIJ)
            IJS = IJS + NIJ
            NIJOFF = NIJOFF + NIJ
         ENDDO
         NKLOFF = NKLOFF + NIJ
         CALL MEMREL('MS4IN2D',WORK,1,KHMAT,KFREE,LFREE)
      ENDDO
C
      CALL QEXIT('MS4IN2D')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN2P*/
      SUBROUTINE MS4IN2P(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,
     & IREPIJ,IJS,NIJ,Q3S,Q4S,Q3T,Q4T,HMAT,INDXKR,INDXB12,
     & INDXB34)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               Second index-pair transformation.
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - NSTR3        Number of active spinors for index 3
C     - NSTR4        Number of active spinors for index 4
C     - IREPIJ       Fermion symmetry (parity)
C     - IJS          First blocks in this pass
C     - NIJ          Number of blocks in this pass
C     - Q3S          Coefficients for index 3
C     - Q4S          Coefficients for index 4
C     - Q3T          Coefficients for index 3 but for transposed AO ints
C     - Q4T          Coefficients for index 4 but for transposed AO ints
C     - HMAT         Symmetry packed half-transformed integrals
C
C     Output is written directly in the MOLFDIR-type file MDCINT
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION HMAT(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN2P')
      KFRSAV = KFREE
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C
C     Number of boson symmetry for each fermion irrep (parity)
C
      NBSYMP = NBSYM/NFSYM
C
C     Transform indices 3 and 4 to molecular spinor basis
C
      NKL = NFPCK34(IREPIJ)
      CALL MEMGET('REAL',KHMAT,NKL*NZ*NZ*NBSYMP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHTMP,NKL*NZ,WORK,KFREE,LFREE)
      CALL DZERO (WORK(KHTMP),NKL*NZ)
C
      IOFF = 1
C
C     Loop over the number of right hand side (12) batches.
C
      DO IJ = 1, NIJ
        JOFF = 0
        CALL DZERO(WORK(KHMAT),NKL*NZ*NZ*NBSYMP)
C
C       Loop over all boson irreps. On the next line we pick out the
C       ones with the right parity.
C
        DO IREPRS = 0, NBSYM-1
          IF (IREPIJ.EQ.JBTOF(IREPRS,1)) THEN
C
C           Loop over the quaternion units of the first transformed
C           density (this is IZ1 in my notes).
C
            DO IZ = 1, NZ
C
C              Loop over the compound component (LL/SS) of the second
C              untransformed density.
C
               DO IC = ICS, ICF
                  IF (ICS.EQ.2) IOFF = IOFF + NSPCK34(IREPRS,1)
C
C                 Loop over the boson symmetries of the fourth
C                 (untransformed) index. Get the fermion symmetry
C                 of that index and get the boson and fermion symmetry
C                 of the third index.
C
                  DO IREPS = 0, NBSYM-1
                     IREPR = IBTXOR(IREPRS,IREPS)
                     IREPK =  JBTOF(IREPR,IC)
                     IREPL =  JBTOF(IREPS,IC)
                     KOFF = ICMOQS(IREPK,3) + IBBAS3(IREPR,IC)
     &                    - IBAS3(IREPK)
                     LOFF = ICMOQS(IREPL,4) + IBBAS4(IREPS,IC)
     &                    - IBAS4(IREPL)
                     NR = NBBAS3(IREPR,IC)
                     NS = NBBAS4(IREPS,IC)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     NRQ3 = NDMOQS(1,IREPK,3)
                     NCQ3 = NDMOQS(2,IREPK,3)
                     NRQ4 = NDMOQS(1,IREPL,4)
                     NCQ4 = NDMOQS(2,IREPL,4)
                     IF (IC.EQ.1) THEN
                        IREPRI = IREPR
                        IREPSI = IREPS
                     ELSE
                        IREPRI = IBTXOR(IXYZ,IREPR)
                        IREPSI = IBTXOR(IXYZ,IREPS)
                     ENDIF
C
C                    Do the actual transformation. 

                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        CALL QTRANS('AOMO','S',D0,NR,NS,NK,NL,
     &                       HMAT(IOFF),NR,NS,1,IPQTOQ(1,0),
     &                       WORK(KHTMP),NK,NL,
     &                       NZ,IPQTOQ(1,IREPRS),
     &                       Q3S(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPRI),
     &                       Q4S(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPSI),
     &                       WORK(KFREE),LFREE,IPRINT)
C
C                       Order the integrals with NKL as first and IZ2
C                       as second index
C
                        DO IZ2 = 1, NZ
                           JOFF1 = JOFF+(IZ2-1)*NKL+IFPCK34(IREPK,IREPL)
                           JOFF2 = KHTMP + (IZ2-1)*NK*NL
                           CALL DAXPY(NK*NL,D1,WORK(JOFF2),1,
     &                                WORK(KHMAT+JOFF1),1)
                        ENDDO
                     ENDIF
C
C                    Repeat this transformation with the transposed
C                    AO-matrix (the half-transformed integrals were
C                    triangular in the AO-indices). Diagonal elements
C                    are scaled before storage to avoid double-counting.
C                    
                     IREPK =  JBTOF(IREPS,IC)
                     IREPL =  JBTOF(IREPR,IC)
                     KOFF = ICMOQT(IREPK,3) + IBBAS4(IREPS,IC)
     &                    - IBAS4(IREPK)
                     LOFF = ICMOQT(IREPL,4) + IBBAS3(IREPR,IC)
     &                    - IBAS3(IREPL)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     NRQ3 = NDMOQT(1,IREPK,3)
                     NCQ3 = NDMOQT(2,IREPK,3)
                     NRQ4 = NDMOQT(1,IREPL,4)
                     NCQ4 = NDMOQT(2,IREPL,4)
C
                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        CALL MEMGET('REAL',KTRMAT,NR*NS,
     &                              WORK,KFREE,LFREE)
                        CALL MTRSP (NR,NS,HMAT(IOFF),NR,WORK(KTRMAT),NS)
                        CALL QTRANS('AOMO','S',D0,NS,NR,NK,NL,
     &                       WORK(KTRMAT),NS,NR,1,IPQTOQ(1,0),
     &                       WORK(KHTMP),NK,NL,NZ,IPQTOQ(1,IREPRS),
     &                       Q3T(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPSI),
     &                       Q4T(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPRI),
     &                       WORK(KFREE),LFREE,IPRINT)
                        CALL MEMREL('MS4IN2P',WORK,1,KTRMAT,KFREE,LFREE)
C
C                       Order the integrals with NKL as first and IZ2
C                       as second index
C
                        DO IZ2 = 1, NZ
                           JOFF1 = JOFF+(IZ2-1)*NKL+IFPCK34(IREPK,IREPL)
                           JOFF2 = KHTMP + (IZ2-1)*NK*NL
                           CALL DAXPY(NK*NL,D1,WORK(JOFF2),1,
     &                                WORK(KHMAT+JOFF1),1)
                        ENDDO
                     ENDIF
C
                     IOFF = IOFF + NR * NS
                  ENDDO
               ENDDO
               JOFF = JOFF + NKL*NZ
            ENDDO
C
          ENDIF
        ENDDO
C
C       We can symmetrize the integral here or write it to a
C       scratch file.
C
        IF (ISTRAT.EQ.3) THEN
           IKR = INDXKR(1,IJ)
           JKR = INDXKR(2,IJ)
           CALL SYMFINT(IPRINT,IREPIJ,INDXB12(1,IJ),INDXB34,IJ,IKR,JKR,
     &                  NSTR3,NSTR4,NKL,WORK(KHMAT),
     &                  WORK,KFREE,LFREE)
        ELSEIF (ISTRAT.EQ.4) THEN
           CALL MDSCRU(IPRINT,IREPIJ,IJS+IJ,NKL,WORK(KHMAT)
     &                 ,WORK,KFREE,LFREE)
        ENDIF
      ENDDO
C
      CALL MEMREL('MS4IN2P',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('MS4IN2P')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INITBUF */
      SUBROUTINE INITBUF(IRECC,NBFSZ,NBUF,JBOFF)
C
C     Initialize the buffer file and index arrays for the half-
C     transformed integrals
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION NBFSZ(2),NBUF(2),JBOFF(2,0:*),IRECC(2)
      CHARACTER*10 LFNAM,FNODE,FNAME*5
C
      DO IRREP = 1, NFSYM
         NWORDS = NBFSZ(IRREP)
C
C     The record size needs to be nwords real + nwords + 1 integer
C
         IRECL = IRECLEN(NWORDS,NWORDS+1,0)
         LUTRA11 = LUTRA1 + IRREP
         WRITE (FNAME,'(A4,I1)') 'HTIN',IRREP
         FNODE = LFNAM(FNAME)
         OPEN (UNIT=LUTRA11,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
         IRECC(IRREP) = 1
      ENDDO
C
      JBOFF1 = 1
      DO IRREP = 1, NFSYM
         DO IBUF = 1, NBUF(IRREP)
            JBOFF(IRREP,IBUF) = JBOFF1
            JBOFF1 = JBOFF1 + NBFSZ(IRREP)
         ENDDO
         JBOFF(IRREP,0) = 0
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DELEBUF */
      SUBROUTINE DELEBUF
C
C     Delete the buffer files for the half-transformed integrals
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
C
      DO IRREP = 1, NFSYM
         LUTRA11 = LUTRA1 + IRREP
         CLOSE (UNIT=LUTRA11,STATUS='DELETE')
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck FINABUF */
      SUBROUTINE FINABUF(IRECC,NBFSZ,NBUF,IBOFF,JBOFF,RBUF,LBUF1)
C
C     Write the last set of (half-filled) buffers
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION IRECC(2),NBFSZ(2),NBUF(2),JBOFF(2,0:*),IBOFF(2)
      DIMENSION RBUF(*),LBUF1(*)
C
      DO IREPIJ = 1, NFSYM
         IF (JBOFF(IREPIJ,0).NE.0) THEN
C
            IREC = IRECC(IREPIJ)
            JBUFF = IBOFF(IREPIJ)
            DO IBUF = 1, NBUF(IREPIJ)
               JBUFT = JBOFF(IREPIJ,IBUF) - JBUFF
               JBUFL = JBUFF + NBFSZ(IREPIJ) - 1
C
               WRITE (LUTRA1+IREPIJ,REC=IREC) JBUFT,
     &               (RBUF(JBUF),JBUF=JBUFF,JBUFL),
     &               (LBUF1(JBUF),JBUF=JBUFF,JBUFL)
               JBOFF(IREPIJ,IBUF) = JBUFF
               JBUFF = JBUFL + 1
               IREC = IREC + 1
            ENDDO
            IRECC(IREPIJ) = IREC
            JBOFF(IREPIJ,0) = 0
         ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck READBUF */
      SUBROUTINE READBUF(IREPIJ,IBUF,IOFF,IRECC,NBFSZ,NBUF,
     &                   RBUF,LBUF1,HMAT)
C
C     Read a subset of the halftransformed integrals into memory
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
#include <dcbtra.h>
#include <dgroup.h>
      DIMENSION RBUF(NBFSZ),LBUF1(NBFSZ),HMAT(*)
C
      DO IREC = IBUF, IRECC-1, NBUF
         READ (LUTRA1+IREPIJ,REC=IREC) NABFSZ,RBUF,LBUF1
C
C        Distribute the integrals. 
C
         DO JBUF = 1, NABFSZ
            HMAT(LBUF1(JBUF)-IOFF) = RBUF(JBUF)
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN1D*/
      SUBROUTINE MS4IN1D(WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,LGREC,
     &                   ICS,ICF,IC,
     &                   TRIAN,NSTR1,NSTR2,
     &                   IJPASS,Q1,Q2,
     &                   GMAT,NSIZG,MHTM,IHTM,INDX,
     &                   INDXAB,IRECC,JBOFF,IBOFF,NBUF,
     &                   NIJBUF,NBFSZ,RBUF,LBUF1)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher May 1997.
C
C     PURPOSE : Driver for the first index transformation.
C               Loop over batches for a given shell combination and
C               buffer the half-transformed integrals to file
C
C     Input :
C
C     - IPRINT       Print flag
C     - NPASS        Number of passes through the scalar integrals
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
C     - ICL          Last class of integrals
C     - IC           Component of the right hand : (XX|LL) or (XX|SS)
C     - TRIAN        Store the transformed integrals as a lower triangle
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
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
C     - INDXAB       Information for a pair of boson function belonging
C                    to the current shell pair, 
C     - IRECC        Current record on the half-transformed integr files
C     - JBOFF        Number of integrals in the buffers
C     - IBOFF        Starting adress for buffers of a fermion ircop
C     - NBUF         Number of buffers for half-transformed integrals
C     - NIJBUF       Size of buffers for half-transformed integrals
C     - NBFSZ        Size buffer arrays for half-transformed integrals
C     - RBUF         Buffer arrays with value of half-transformed integ
C     - LBUF1        Buffer arrays with label of half-transformed integ
C
C     Output :
C
C     -              Updated values of buffer blocks
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL   TRIAN
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(NSIZG)
      DIMENSION IJPASS(0:7,NPASS),INDX(3,*)
      LOGICAL   SUB3,SUB4
C
C     For the buffered input
C
      DIMENSION LGREC(NPASS)
C
C     For the buffered output
C
      DIMENSION IRECC(2),IBOFF(2),NBFSZ(2),NBUF(2),NIJBUF(2)
      DIMENSION JBOFF(2,0:*),RBUF(*),LBUF1(*)
C
      DIMENSION IHM(0:7)
C
#include <dgroup.h>
#include <dcbtra.h>
#include <dcbtr3.h>
#include <dcbibt.h>
C
      CALL QENTER('MS4IN1D')
      KFRSAV = KFREE
C
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
             CALL MEMREL('MS4IN1D',WORK,KFRSAV,KRBUF,KFREE,LFREE)
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
         CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KHMAT),NSIZH*NZ)
C
C        Do first step of 4-index transformation :
C        Transform first pair of indices
C 
         CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,
     &                TRIAN,NSTR1,NSTR2,NDMOQR,
     &                ICMOQR,IJPASS(0,IPASS),Q1,Q2,
     &                GMAT,WORK(KHMAT))
C
C        Scatter half-transformed integrals to the right position
C        Write them out to disk 
C
         CALL MS4INW (WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,
     &                IC,IHM,MHTM,IHTM,INDX,INDXAB,IPASS,
     &                WORK(KHMAT),IRECC,JBOFF,IBOFF,NBUF,NIJBUF,
     &                NBFSZ,RBUF,LBUF1)
         CALL MEMREL('MS4IN1D',WORK,KFRSAV,KHMAT,KFREE,LFREE)
C
      ENDDO
C
      CALL QEXIT('MS4IN1D')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck READGBF */
      SUBROUTINE READGBF(LGFIL,IREC,RBUF,LBUF,NGBFSZ,GMAT)
C
C     Read a batch of scalar integrals into memory
C
C     Written by Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      DIMENSION RBUF(NGBFSZ),LBUF(NGBFSZ),GMAT(*)
C
      IF (IREC.EQ.0) RETURN
C
C     Read a batch and the link to the previous record.
C
    1 READ (LGFIL,REC=IREC) IRECN,NABFSZ,RBUF,LBUF
      IREC = IRECN
C
C     Distribute the integrals. 
C
      DO JBUF = 1, NABFSZ
         GMAT(LBUF(JBUF)) = RBUF(JBUF)
      ENDDO
C
C     We read in reverse order, if we have the first block we're done
C
      IF (IREC.NE.0) GOTO 1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck FINAGBF */
      SUBROUTINE FINAGBF (LGFIL,NPASS,IREC,LGREC,LGBUF,NGBFSZ,
     &                    RGBUF,IGBUF)
C
C     Write the last set of (half-filled) buffers of scalar integrals
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      DIMENSION LGREC(NPASS),LGBUF(NPASS)
      DIMENSION RGBUF(NGBFSZ,NPASS),IGBUF(NGBFSZ,NPASS)
C
      DO IPASS = 1, NPASS
         IF (LGBUF(IPASS).NE.0) THEN
C
            IREC = IREC + 1
            WRITE (LGFIL,REC=IREC) LGREC(IPASS),LGBUF(IPASS),
     &            (RGBUF(JBUF,IPASS),JBUF=1,NGBFSZ),
     &            (IGBUF(JBUF,IPASS),JBUF=1,NGBFSZ)
            LGREC(IPASS) = IREC
            LGBUF(IPASS) = 0
         ENDIF
      ENDDO
c     print*,'Number of buffers written', IREC
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INITGBF */
      SUBROUTINE INITGBF(LGFIL,NGBFSZ,NPASS,IREC,LGREC,LGBUF)
C
C     Initialize the buffer file and index arrays for the scalar
C     integrals
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
      DIMENSION LGREC(NPASS),LGBUF(NPASS)
      CHARACTER*10 FNODE,LFNAM
C
      NWORDS = NGBFSZ
c     print*,'NGBFSZ',NGBFSZ
C
C     The record size needs to be nwords real + (nwords + 2) integer
C
      IRECL = IRECLEN(NWORDS,NWORDS+2,0)
C
      FNODE = LFNAM('SCLIN')
      OPEN (UNIT=LGFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
C
      IREC = 0
      CALL IZERO(LGREC,NPASS)
      CALL IZERO(LGBUF,NPASS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DELGBUF */
      SUBROUTINE DELGBUF(LGFIL)
C
C     Delete the buffer files for the scalar integrals
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <implicit.h>
#include <priunit.h>
C
      CLOSE (UNIT=LGFIL,STATUS='DELETE')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LFNAM */
      CHARACTER*10 FUNCTION LFNAM (FILNAM)
C
C     Add the node identification to the file name for shared disk systems
C     Written by Luuk Visscher, February 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <dcbgen.h>
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#else
#include <fpvm3h.h>
#endif
      CHARACTER*5  FILNAM
      CHARACTER*10 LFILNAM
#if defined (VAR_MPI)
      IF (PARCAL) THEN
         NODE = MYTID
      ELSE
         NODE = 0
      ENDIF
#else
      NODE = 0
#endif
C
 1000 FORMAT (A5,I5)
C
      WRITE (LFILNAM,1000) FILNAM,NODE
      DO I = 6, 10
         IF (LFILNAM(I:I).EQ.' ') LFILNAM(I:I) = 'X'
      ENDDO
C
      LFNAM = LFILNAM
C
      RETURN
      END
