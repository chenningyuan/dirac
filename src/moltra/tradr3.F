!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRDR3T */
      SUBROUTINE TRDR3T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                  NSTR,ANTIS,LMP2,INDX,KQ,KE,KIBE,DINTSKP)
C
C     Written by Luuk Visscher, december 1996
C     Removed by Luuk Visscher, july 2002
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      LOGICAL ANTIS, LMP2, NOPV, NODV, TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION WORK(*)
C
      CALL QUIT ("ASKED FOR OBSOLETE STRATEGY 3")
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MS4INW*/
      SUBROUTINE MS4INW (WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,
     &                   IC,IHM,MHTM,IHTM,INDX,INDXAB,IPASS,
     &                   HMAT,IRECC,JBOFF,IBOFF,NBUF,NIJBUF,NBFSZ,
     &                   RBUF,LBUF1)

C
C     Sort half-transformed integrals
C     This version writes them to disk
C
C     Luuk Visscher
C
C
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
C
      PARAMETER (HALF=0.5D0)
C
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5),IOFF4(0:7)
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
#include "ibtfun.h"
C
C     Find out the component type of the functions
C
      CALL ICTYPES(IC,IC1,IC2)
C
C     Loop over the functions inside the shell-block.
C     We do it this way rather than just over the boson irreps of the
C     density and then over the pairs within that irrep, because we
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
            IF (IC.EQ.3) IREPIJ = JBTOF(IREPAB,2)
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
            INDA1 = INDA - ICOS(IREPA+1,IC1) ! Absolute position
            INDB1 = INDB - ICOS(IREPB+1,IC2) ! Absolute position
            IF (SUB3) INDA1 = INDX(3,INDA)  ! Relative position
            IF (SUB4) INDB1 = IOFF4(IREPB) + INDX(3,INDB)  ! Relative position
C
C           The offset within the untransformed indexes. These were
C           the right-hand side (slowest varying) and are sorted to
C           become the left hand side. Symmetrization is done
C           in the second half-transformation (the function indices are
C           left triangular), but we scale the diagonal elements to 
C           avoid double counting.
C
            JOFFAB = ISPCK34(IREPA,IREPB,IC)
     &             + (INDB1-1)*NBBAS3(IREPA,IC1)+INDA1
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
               JOFF1 = JOFF1 + NSPCK34(IREPAB,IC)
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      DIMENSION IRECC(2),MHTM(2),NBUF(2),NIJBUF(2),NBFSZ(2)
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION HMAT(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      LOGICAL CORRECT_PARITY
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
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
      LGMAT = NKL*NZ*NZ*NBSYMP
      CALL MEMGET('REAL',KGMAT,LGMAT,WORK,KFREE,LFREE)
      IF (ISTRAT.EQ.4) THEN
C        for MDSCRU and for NFSYM.gt.1
         CALL MEMGET('REAL',KGTMP,LGMAT,WORK,KFREE,LFREE)
      ELSE IF (NFSYM.GT.1) THEN
         CALL MEMGET('REAL',KGTMP,NKL*NZ,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KGTMP,     0,WORK,KFREE,LFREE)
      END IF
      NR = 0
      NS = 0
      DO IC = ICS,ICF
         CALL ICTYPES(IC,IC3,IC4)
         DO IREPRS = 0,NBSYM-1
            NR = MAX(NR,NBBAS3(IREPRS,IC3))
            NS = MAX(NS,NBBAS4(IREPRS,IC4))
         END DO
      END DO
      CALL MEMGET('REAL',KTRMAT,NR*NS,WORK,KFREE,LFREE)
C
      IOFF = 1
C
C     Loop over the number of right hand side (12) batches.
C
      DO IJ = 1, NIJ
        JGMAT = KGMAT
        IF (NFSYM.EQ.1) CALL DZERO(WORK(KGMAT),NKL*NZ*NZ*NBSYMP)
C
C       Loop over all boson irreps. On the next line we pick out the
C       ones with the right parity. Note that for Gaunt we need to take
C       into account the difference in parity between large and small
C       components by using JBTOF(..,2) instead of JBTOF(..,1).
C
        DO IREPRS = 0, NBSYM-1
          CORRECT_PARITY = IREPIJ.EQ.JBTOF(IREPRS,1)
          IF (ICS.EQ.3) CORRECT_PARITY = IREPIJ.EQ.JBTOF(IREPRS,2)
          IF (CORRECT_PARITY) THEN
C
C           Loop over the quaternion units of the first transformed
C           density (this is IZ1 in my notes).
C
            DO IZ = 1, NZ
               IF (NFSYM.GT.1) CALL DZERO (WORK(KGTMP),NKL*NZ)
C
C              Loop over the compound component (LL/SS) of the second
C              untransformed density.
C
               IF (ICS.EQ.2) IOFF = IOFF + NSPCK34(IREPRS,1)
               DO IC = ICS, ICF
                  CALL ICTYPES(IC,IC3,IC4)
C
C                 Loop over the boson symmetries of the fourth
C                 (untransformed) index. Get the fermion symmetry
C                 of that index and get the boson and fermion symmetry
C                 of the third index.
C
                  DO IREPS = 0, NBSYM-1
                     IREPR = IBTXOR(IREPRS,IREPS)
                     NR = NBBAS3(IREPR,IC3)
                     NS = NBBAS4(IREPS,IC4)
                     IF (IC3.EQ.1) THEN
                        IREPRI = IREPR
                     ELSE
                        IREPRI = IBTXOR(IXYZ,IREPR)
                     ENDIF
                     IF (IC4.EQ.1) THEN
                        IREPSI = IREPS
                     ELSE
                        IREPSI = IBTXOR(IXYZ,IREPS)
                     ENDIF
C
                     IREPK =  JBTOF(IREPR,IC3)
                     IREPL =  JBTOF(IREPS,IC4)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        KOFF = ICMOQS(IREPK,3) + IBBAS3(IREPR,IC3)
     &                       - IBAS3(IREPK)
                        LOFF = ICMOQS(IREPL,4) + IBBAS4(IREPS,IC4)
     &                       - IBAS4(IREPL)
                        NRQ3 = NDMOQS(1,IREPK,3)
                        NCQ3 = NDMOQS(2,IREPK,3)
                        NRQ4 = NDMOQS(1,IREPL,4)
                        NCQ4 = NDMOQS(2,IREPL,4)
C
C                       Do the actual transformation. 
C
C                       Because QTRANS will generate them in a different order
C                       for cases with inversion symmetry, we write to the
C                       temporary array KGTMP and sort later on
C                       ELSE: place them directly in array KGMAT
C
                        IF (NFSYM.GT.1) THEN
                           KLOFF = KGTMP + IFPCK34(IREPK,IREPL)*NZ
                        ELSE
                           KLOFF = JGMAT
                        END IF
C
                        CALL QTRANS('AOMO','S',D1,NR,NS,NK,NL,
     &                       HMAT(IOFF),NR,NS,1,IPQTOQ(1,0),
     &                       WORK(KLOFF),NK,NL,
     &                       NZ,IPQTOQ(1,IREPRS),
     &                       Q3S(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPRI),
     &                       Q4S(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPSI),
     &                       WORK(KFREE),LFREE,IPRINT)
#if 0
                        if ( ij .eq. 1 ) then
                           call header('ij = 1',-1 )
                           write(6,*) 'hmat'
                           call prqmat(hmat(ioff),nr,ns,nr,ns,
     $                          nz,ipqtoq(1,0),6)
                           write(6,*) 'khtmp'
                           call prqmat(work(khtmp),nk,nl,nk,nl,
     $                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q4s'
                           call prqmat(q3s(koff),nr,nl,nrq3,ncq3,
     $                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q3s'
                           call prqmat(q4s(loff),ns,nk,nrq4,ncq4,
     $                          nz,ipqtoq(1,0),6)
                        end if
#endif
                     ENDIF
C
C                    Repeat this transformation with the transposed
C                    AO-matrix (the half-transformed integrals were
C                    triangular in the AO-indices). Diagonal elements
C                    are scaled before storage to avoid double-counting.
C                    
                     IREPK =  JBTOF(IREPS,IC4)
                     IREPL =  JBTOF(IREPR,IC3)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        KOFF = ICMOQT(IREPK,3) + IBBAS4(IREPS,IC4)
     &                       - IBAS4(IREPK)
                        LOFF = ICMOQT(IREPL,4) + IBBAS3(IREPR,IC3)
     &                       - IBAS3(IREPL)
                        NRQ3 = NDMOQT(1,IREPK,3)
                        NCQ3 = NDMOQT(2,IREPK,3)
                        NRQ4 = NDMOQT(1,IREPL,4)
                        NCQ4 = NDMOQT(2,IREPL,4)
C
                        IF (NFSYM.GT.1) THEN
                           KLOFF = KGTMP + IFPCK34(IREPK,IREPL)*NZ
                        ELSE
                           KLOFF = JGMAT
                        END IF
C
                        CALL MTRSP (NR,NS,HMAT(IOFF),NR,WORK(KTRMAT),NS)
                        CALL QTRANS('AOMO','S',D1,NS,NR,NK,NL,
     &                       WORK(KTRMAT),NS,NR,1,IPQTOQ(1,0),
     &                       WORK(KLOFF),NK,NL,NZ,IPQTOQ(1,IREPRS),
     &                       Q3T(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPSI),
     &                       Q4T(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPRI),
     &                       WORK(KFREE),LFREE,IPRINT)
#ifdef UNDEF
                        if ( ij .eq. 1 ) then
                           call header('ij = 1',-1 )
                           write(6,*) 'ktrmat'
                           call prqmat(work(ktrmat),ns,nr,ns,nr,
     $                          1,ipqtoq(1,0),6)
                           write(6,*) 'khtmp'
                           call prqmat(work(khtmp),nk,nl,nk,nl,
     $                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q3t'
                           call prqmat(q3t(koff),ns,nk,nrq3,ncq3,
     $                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q4t'
                           call prqmat(q4t(loff),nr,nl,nrq4,ncq4,
     $                          nz,ipqtoq(1,0),6)
                        end if
#endif
                     ENDIF
C
                     IOFF = IOFF + NR * NS
                  ENDDO
C                 ... ENDDO IREPS loop
               ENDDO
C              ... ENDDO IC loop
C
               IF (NFSYM.GT.1) THEN
C
C                 Order the integrals with NKL as first and IZ2
C                 as second index
C
                  DO IZ2 = 1, NZ
                     DO IREPK = 1, NFSYM
                        IREPL = MOD(IREPK+IREPIJ,2) + 1
                        NK = NSTR3(IREPK)
                        NL = NSTR4(IREPL)
                        JOFF1 = JGMAT+(IZ2-1)*NKL+IFPCK34(IREPK,IREPL)
                        JOFF2 = KGTMP + IFPCK34(IREPK,IREPL)*NZ
     &                                + (IZ2-1)*NK*NL
                        CALL DCOPY(NK*NL,WORK(JOFF2),1,WORK(JOFF1),1)
                     ENDDO
                  ENDDO
               ENDIF
C
               JGMAT = JGMAT + NKL*NZ
            ENDDO
C           ... ENDDO IZ loop
C
          ENDIF
        ENDDO
C       ... ENDDO IREPRS loop
C
C       We can symmetrize the integral here or write it to a
C       scratch file.
C
        IF (ISTRAT.EQ.3) THEN
           IKR = INDXKR(1,IJ)
           JKR = INDXKR(2,IJ)
           CALL SYMFINT(IPRINT,IREPIJ,INDXB12(1,IJ),INDXB34,IJ,IKR,JKR,
     &                  NSTR3,NSTR4,NKL,WORK(KGMAT))
        ELSEIF (ISTRAT.EQ.4) THEN
           CALL MDSCRU(IREPIJ,IJS+IJ,LGMAT,WORK(KGMAT),WORK(KGTMP))
        ENDIF
      ENDDO
C     ... ENDDO IJ loop
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
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
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
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
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
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
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
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
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
      SUBROUTINE MS4IN1D(WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,LGREC,
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL   TRIAN
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(NSIZG)
      DIMENSION IJPASS(0:7,NPASS),INDX(3,*),IOFF4(0:7)
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
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
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
             CALL READGBF(LGFIL,IREC,NGBFSZ,GMAT)
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
         CALL MS4INW (WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,
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
      SUBROUTINE READGBF(LGFIL,IREC,NGBFSZ,GMAT)
C
C     Read a batch of scalar integrals into memory
C
C     Written by Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
C
      Real(8)              :: GMAT(*)
      Real(8), allocatable :: RBUF(:)
      Integer, allocatable :: LBUF(:)
	  Integer, intent(in)  :: NGBFSZ
	  Integer :: istat
	  
      IF (IREC.EQ.0) RETURN
! miro: provide with checking parameter
      Allocate (RBUF(NGBFSZ),stat=istat)
	  If (istat.ne.0) then
	    Write(LUPRI,*) 
     & "Real*8 RBUF(NGBFSZ) allocation failed, NGBFSZ=",NGBFSZ
		Call Quit('Real*8 RBUF(NGBFSZ) allocation failed !')
	  EndIf
      Allocate (LBUF(NGBFSZ),stat=istat)
	  If (istat.ne.0) then
	    Write(LUPRI,*)
     &  "Integer LBUF(NGBFSZ) allocation failed, NGBFSZ=",NGBFSZ
		Call Quit('Integer LBUF(NGBFSZ) allocation failed !')
	  EndIf	  
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
      DeAllocate (RBUF,stat=istat)
      DeAllocate (LBUF,stat=istat)
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
#include "implicit.h"
#include "priunit.h"
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
#include "implicit.h"
#include "priunit.h"
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
#include "implicit.h"
#include "priunit.h"
C
      CLOSE (UNIT=LGFIL,STATUS='DELETE')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ICTYPES */
      SUBROUTINE ICTYPES(IC,IC1,IC2)
C
C     Get the individual components from the compound IC label
C     Written by Luuk Visscher, July 2002
C
      IF (IC.LE.2) THEN
         IC1 = IC
         IC2 = IC
      ELSEIF (IC.EQ.3) THEN
         IC1 = 2
         IC2 = 1
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gauntr */
      SUBROUTINE GAUNTR(FADD,NRAO,NCAO,NRMO,NCMO,
     &                  FAO,LRAO,LCAO,NZAO,IQAO,
     &                  FMO,LRMO,LCMO,NZMO,IQMO,
     &                  TM1,LR1,LC1,NZTM1,IQTM1,
     &                  TM2,LR2,LC2,NZTM2,IQTM2,
     &                  BUF,LBUF,IPRINT)
C*****************************************************************************
C
C     Dedicated routine based on QTRANS to do transformation of a vector
C     of spin matrices as required in the transformation of the Gaunt
C     integrals. Only AO-to-MO-transformation.
C
C     ADDF gives the possibility of adding the results
C
C     The unitary C-matrix has a time-symmetric structure.
C
C     The AO-to-MO-transformation can be set up quaternionically as
C
C       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
C
C       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
C
C     Written by L.Visscher : July 2002
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),
     &          FMO(LRMO,LCMO,NZMO,3),IQMO(*),
     &          TM1(LR1,LC1,*),IQTM1(*),
     &          TM2(LR2,LC2,*),IQTM2(*),
     &          BUF(LBUF)
      DIMENSION IQW(4),IQWA(4)
C
      CALL QENTER('GAUNTR')
C
C     *******************************************************
C     *****  AO - to - MO -transformation : F' = (C+)FC *****
C     *******************************************************
C
C
C     Look for the most efficient transformation
C       ------------------------------------------
C
      IF(NCAO.GE .NRAO) THEN
C
C       Perform the two steps:
C         1. W  = FC
C         2. F' = (C+)W
C
        CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
        NBUF = NRAO*NCMO*NZW
        IF(NBUF.GT.LBUF) THEN
           WRITE (LUPRI,'(/A,I10/A,I10)')
     & ' >>> GAUNTR error, need work space    ',NBUF,
     & '                   current work space ',LBUF
        ENDIF
C
C       First part of transformation: W = FC
C       ------------------------------------
C
        CALL QGEMM(NRAO,NCMO,NCAO,D1,
     &           'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &           'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &      D0,IQW,BUF,NRAO,NCMO,NZW)
C
C       Second part of transformation: H = (C+)W
C       ----------------------------------------
C
        DO IALPHA = 2, 4
C
C         Multiply the one-index transformed matrix
C         by the quaternion unit corresponding to the
C         alpha matrix.
C
          CALL IQPACK(IALPHA,1,IQW,NZW,IQWA,NZWA)
          CALL QGEMM(NRMO,NCMO,NRAO,D1,
     &       'H','N',IQTM1,TM1 ,LR1,LC1,NZTM1,
     &       'N','N',IQWA,BUF,NRAO,NCMO,NZWA,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
        ENDDO
      ELSE
C
C       Perform the two steps:
C         1. W  = (C+)F
C         2. F' = WC
C
        CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
        NBUF = NRMO*NCAO*NZW
        IF(NBUF.GT.LBUF) THEN
           WRITE (LUPRI,'(/A,I10/A,I10)')
     & ' >>> GAUNTR error, need work space    ',NBUF,
     & '                   current work space ',LBUF
        ENDIF
C
C       First part of transformation: W = (C+)F
C       ----------------------------------------
C
        CALL QGEMM(NRMO,NCAO,NRAO,D1,
     &     'H','N',IQTM1,TM1,LR1,LC1,NZTM1,
     &     'N','N',IQAO,FAO,LRAO,LCAO,NZAO,
     &      D0,IQW,BUF,NRMO,NCAO,NZW)
C
C       Second part of transformation: H = WC
C       -------------------------------------
C
        DO IALPHA = 2, 4
          CALL IQPACK(IQW,NZW,IALPHA,1,IQWA,NZWA)
          CALL QGEMM(NRMO,NCMO,NCAO,D1,
     &       'N','N',IQWA,BUF,NRMO,NCAO,NZWA,
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
C
        ENDDO
      ENDIF
C
      CALL QEXIT('GAUNTR')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
