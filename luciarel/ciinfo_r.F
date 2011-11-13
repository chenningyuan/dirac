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

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
* Lucia.f : GAS implementing (no pair) relativistic Theory (YES !!!!)
*
      SUBROUTINE DBG_INT_LISTS(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
     &                         ISPINFREE)
*
* Number of double group integral lists to be constructed
*
*     Jeppe Olsen, Feb. 98
*
* Reordered lists for convenience: We only need the list
* (uu|uu) in an initial or plain SOCI calculation. So this
* is going to be list number 1. Allocation modified as well.
* Modifications in GETINT_DBG incorporated (and elsewhere).
* We are furthermore using real wave functions. This implies
* That complex conjugation symmetry can also be used in the
* quaternion matrix groups C1 and Ci, because their irrep
* basis functions are not generally complex. But if we
* generalize the code to e.g. atomic complex 4- (or 2-)spinors
* or the like, these things have to be distinguished.
* So we will have 3 cases:
*
*   1) Initial or plain SOCI: 1 list
*   2) SOCI + spinor optimization in some way (INS or MCSCF): 3 lists
*   3) Complex wave function: 6 lists
*
* Here, time-reversal symmetry has already been accounted for in
* all cases, so these are non-redundant integral lists.
* At the moment, the first two cases are assumed to be implemented,
* so we will distinguish them by the input keyword DENSI, i.e.
* whether density matrices are to be calculated or not.
*
*     Timo Fleig, June 1999
*
* revised for DIRAC environment
*     Timo Fleig, December 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.output
      INTEGER IUB1(2,10),IUB2(4,10),ISM2(3,10)
*
*. ======================================
*. 1e Lists : all four lists constructed
*. ======================================
*
      NL1 = 4
*. 1 : unbarred, unbarred
      IUB1(1,1) = 1
      IUB1(2,1) = 1
*. 2 : unbarred, barred
      IUB1(1,2) = 1
      IUB1(2,2) = 2
*. 3 : barred, unbarred
      IUB1(1,3) = 2
      IUB1(2,3) = 1
*. 4 : barred, barred
      IUB1(1,4) = 2
      IUB1(2,4) = 2
*
*. =========
*. 2e Lists
*. =========
*
*  For the moment, no complex wave function:
      IWFCOMP = 0
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
*=====================
* List 1 ( u u ! u u )
*=====================
*
      IUB2(1,1) = 1
      IUB2(2,1) = 1
      IUB2(3,1) = 1
      IUB2(4,1) = 1
*. No Symmetry between 12 and 34
      ISM2(1,1) = 0
      ISM2(2,1) = 0
      ISM2(3,1) = 0
*
      if (IDENSI.ne.0.or.ISPINFREE.eq.0) then
*
*=====================
* List 2 ( u b ! u u )
*=====================
*
        IUB2(1,2) = 1
        IUB2(2,2) = 2
        IUB2(3,2) = 1
        IUB2(4,2) = 1
*. No Symmetry between 12 and 34
        ISM2(1,2) = 0
        ISM2(2,2) = 0
        ISM2(3,2) = 0
*
*=====================
* List 3 ( b u ! u b )
*=====================
*
        IUB2(1,3) = 2
        IUB2(2,3) = 1
        IUB2(3,3) = 1
        IUB2(4,3) = 2
*. No Symmetry between 12 and 34
        ISM2(1,3) = 0
        ISM2(2,3) = 0
        ISM2(3,3) = 0
*
*=====================
* List 4 ( b u ! b u )
*=====================
*
        IUB2(1,4) = 2
        IUB2(2,4) = 1
        IUB2(3,4) = 2
        IUB2(4,4) = 1
*. No Symmetry between 12 and 34
        ISM2(1,4) = 0
        ISM2(2,4) = 0
        ISM2(3,4) = 0
*
        if (IWFCOMP.eq.1) then
*
*=====================
* List 5 ( u b ! u b )
*=====================
*
          IUB2(1,5) = 1
          IUB2(2,5) = 2
          IUB2(3,5) = 1
          IUB2(4,5) = 2
*. Symmetry between 12 and 34
          ISM2(1,5) = 0
          ISM2(2,5) = 0
          ISM2(3,5) = 0
CTEMP     ISM2(3,5) = 1
*
*=====================
* List 6 ( b u ! u u )
*=====================
*
          IUB2(1,6) = 2
          IUB2(2,6) = 1
          IUB2(3,6) = 1
          IUB2(4,6) = 1
*. No Symmetry between 12 and 34
          ISM2(1,6) = 0
          ISM2(2,6) = 0
          ISM2(3,6) = 0
*
        end if
      end if
*
      if (IDENSI.eq.0.and.ISPINFREE.eq.1) then
        NL2 = 1
      else
        NL2 = 4
        if (IWFCOMP.eq.1) NL2 = 6
      end if
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DETTYP(IPRNT)
*
* Allowed set of determinants : Distribution of
* electrons in alpha(up) -and beta(down) spinors
*
* Defined in terms of
* NMS2VAL : Number of MS2 values (MK2 values )
* MS2VAL  : The allowed MS2 values
*
* Jeppe Olsen, July 97
* Timo Fleig, August '97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
#include "mxpdim.inc"
#include "crun.inc"
#include "lucinp.inc"
*. Output
#include "cstate.inc"
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
        write(6,*)
        write(6,*) ' -------------------------------- '
        write(6,*) '   DETTYP speaking  '
        write(6,*) ' -------------------------------- '
        write(6,*)
        write(6,'(A17,I4,I4)') 'MK2REF_CI,MK2DEL_CI are',
     &                          MK2REF_CI,MK2DEL_CI
      end if
#endif
*
      IF(IDBGCI.EQ.0.AND.IRELACI.EQ.0) THEN
*. Single MS2 space
        NMS2VAL = 1
        MS2VAL(1) = MK2REF_CI
      ELSE
*. Potentially several MS2/MK2 values
        MK2MAX = MK2REF_CI + MK2DEL_CI
*. MK2MAX gives highest MK2 value for possibly coupled strings (?)
        if (MK2MAX.gt.NACTEL) then
          MK2MAX = NACTEL
        else
          call evenodd(IEVODMX,MK2MAX)
          call evenodd(IEVODN,NACTEL)
          if (IEVODN.eq.1.and.IEVODMX.eq.2) MK2MAX = MK2MAX - 1
        end if
        MK2MIN = MK2REF_CI - MK2DEL_CI
        IF(MK2MIN.LT.-NACTEL) MK2MIN = -NACTEL
* Pairing ? Ensure symmetry between +MK and - MK spaces
        IPAIR = 1
        call evenodd(IEVOD,NACTEL)
        IF(IPAIR.EQ.1) THEN
          IF(MK2MIN.GT.0) THEN
            MK2MINP = -MK2MIN
          ELSE
            if (IEVOD.eq.1) then
               MK2MIN = 1
               MK2MINP = -1
            else if (IEVOD.eq.2) then
               MK2MIN = 0
               MK2MINP = -2
            else
               write(6,*) 'Problem in evenodd'
               Call Abend2( 'Quitting in dettyp.' )
            end if
          END IF
          MK2MAXP = -MK2MAX
        ELSE
          if (IEVOD.eq.1) then
             MK2MINP = -1
             MK2MAXP = 1
          else if (IEVOD.eq.2) then
             MK2MINP = -2
             MK2MAXP = 0
          else
             write(6,*) 'Problem in evenodd'
             Call Abend2( 'Quitting in dettyp.' )
          end if
        END IF
*. And construct the MS2 spaces
        NMS2VAL = 0
        DO MS2 = MK2MIN,MK2MAX,2
          NMS2VAL = NMS2VAL + 1
          MS2VAL(NMS2VAL) = MS2
        END DO
        DO MS2 = MK2MAXP,MK2MINP,2                           !sic
          NMS2VAL = NMS2VAL+1
          MS2VAL(NMS2VAL) = MS2
        END DO
      END IF
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of MS2/MK2 values', NMS2VAL
        WRITE(6,*)
        WRITE(6,'(A)') ' Included MS2/MK2 values '
        WRITE(6,'(A)') ' ======================= '
        DO IMS2 = 1, NMS2VAL
          WRITE(6,'(1H ,I3)') MS2VAL(IMS2)
        END DO
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE EXCCLS(NCLS,IACTIN,IACTOUT,IEXC)
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* Master routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*. From the common blocks
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cgas.inc"
*
      CALL EXCCLS_S(NGAS,WORK(KLOCCLS),NCLS,IACTIN,IACTOUT,IEXC)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE EXCCLS_S(NGAS,ICLS,NCLS,IACTIN,IACTOUT,IEXC)
*
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* Slave routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ICLS(NGAS,NCLS)
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*
C     write(6,*) ' ICLS in EXCCLS '
C     call iwrtma(icls,ngas,ncls,ngas,ncls)

      IZERO = 0
      CALL ISETVC(IACTOUT,IZERO,NCLS)
      DO ICLSIN = 1, NCLS
        IF(IACTIN(ICLSIN).NE.0) THEN
        DO ICLSOUT = 1, NCLS
*. Number of anihilations and creations required to connect classes
          NANNI = 0
          NCREA = 0
          DO IGAS = 1, NGAS
*
            IDIF = ICLS(IGAS,ICLSOUT)-ICLS(IGAS,ICLSIN)
            IF(IDIF .GT. 0 ) THEN
              NCREA = NCREA + IDIF
            ELSE IF (IDIF .LT. 0 ) THEN
              NANNI = NANNI - IDIF
            END IF
          END DO
*
          IF(NCREA.LE.IEXC) THEN
            IACTOUT(ICLSOUT) = 1
          END IF
*
        END DO
        END IF
      END DO
*
      NACTOUT = 0
      DO ICLSOUT = 1, NCLS
        IF(IACTOUT(ICLSOUT).GT.0) NACTOUT = NACTOUT + 1
      END DO
*
      NTEST = 1
#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Output from EXCCLS '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' Number of active output classes ',NACTOUT
         WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Active output classes '
         CALL IWRTMA(IACTOUT,NCLS,1,NCLS,1)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Dummy routine as long as development goes on!
      SUBROUTINE FINMO_REL
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
      SUBROUTINE GASANA_REL(C,NBLOCK,IBLOCK,IBLTP,LUC,LUSCR,ICISTR,IRC)
*
*
* Analyze CI vector
*
* Jeppe Olsen, August 1995
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
*
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "cands.inc"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      dimension C(*)
*
      CALL QENTER('ANACI')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRWFC)
*
      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call itods(-1,1,-1,LUSCR)
*
* Loop over real and imaginary parts
      do IRIR = 1,IRC,1
        if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
        if (IRIR.EQ.2) then
*
*. Skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if
        do IMK2_R = 1,NMS2VAL,1
          MK2_R = MS2VAL(IMK2_R)
          IACTP = IST_FOR_DT(1,IMK2_R)
          IBCTP = IST_FOR_DT(2,IMK2_R)
*
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2(IMK2_R-1)
          end if
*
          call REWINE(LUSCR,-1)
          NBL_C = NBLK_MS2(IMK2_R)
          CALL COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          CALL ITODS(-1,1,-1,LUSCR)
#if defined LUCI_DEBUG
          if (NTEST.ge.50) then
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' * Start of contribution from given MK_2     * '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
            WRITE(6,'(4X,A,2I3)')
     &      ' MK2/IRI vector of C :  MK2_R, IRIR =', MK2_R,IRIR
            WRITE(6,*)
          end if
#endif
          write(6,'(/A/A,I4,I3/A)')
     &         ' ----------------------------------------------',
     &         ' Contrib. Kramers proj., real/im  ',MK2_R,IRIR,
     &         ' ----------------------------------------------'
          NAEL = NELEC(IACTP)
          NBEL = NELEC(IBCTP)
*
          NOCTPA = NOCTYP(IACTP)
          NOCTPB = NOCTYP(IBCTP)
*
          IOCTPA = IBSPGPFTP(IACTP)
          IOCTPB = IBSPGPFTP(IBCTP)
*
          IASPGPTP = IACTP
          IBSPGPTP = IBCTP
*
#if defined LUCI_DEBUG
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' ================'
            WRITE(6,*) ' GASANA speaking '
            WRITE(6,*) ' ================'
            WRITE(6,'(A26)') ' IACTP   IBCTP   NAEL   NBEL'
            WRITE(6,'(3X,4I6)')  IACTP,IBCTP,NAEL,NBEL
            WRITE(6,'(A14,2I6)') ' NOCTPA NOCTPB ', NOCTPA,NOCTPB
          END IF
#endif
*
**. Info on block structure of space
*
*
CTF       IF (ICISTR.NE.1) CALL REWINE(LUC,-1)
*. Number of terms to be printed
          IF (NTEST.GE.100) THEN
            THRES = 0.0D0
            MAXTRM = 100000
          ELSE
            THRES = 0.02
            MAXTRM = 500
          END IF
*
          IUSLAB = 0
*. Number of occupation classes
          IWAY = 1
          NEL = NAEL + NBEL
          CALL OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &                   IGSOCC(1,1),IGSOCC(1,2))
*
*. MARK for local allocations
          KDUM = 0
          call memmar(KDUM,DUMMY,'MARK  ',DUMMY,'GASANA')
*. and then the occupation classes
          call memmar(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
          IWAY = 2
          CALL OCCLS_REL(IWAY,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &                   IGSOCC(1,1),IGSOCC(1,2))
*
          call memmar(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
          call memmar(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
          LENGTH = NOCCLS*10
          call memmar(KNCPMT,LENGTH,'ADDL  ',1,'KNCPMT')
          call memmar(KWCPMT,LENGTH,'ADDL  ',2,'KWCPMT')
*
*. Occupation of strings of given sym and supergroup
          CALL GASANAS_REL(C,LUSCR,NBL_C,WORK(KNSTSO(IACTP)),
     &                     WORK(KNSTSO2(IBCTP)),
     &                     NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &                     NBLOCK,IBLOCK(1,IBOFF),
     &                     THRES,MAXTRM,NAEL,NBEL,
     &                     WORK(KLASTR),WORK(KLBSTR),
     &                     IBLTP,NSMST,IUSLAB,
     &                     IDUMMY,WORK(KNCPMT),WORK(KWCPMT),
     &                     NELFSPGP,
     &                     NOCCLS,NGAS,WORK(KLOCCLS),
     &                     IASPGPTP,IBSPGPTP,
     &                     ICISTR,NTOOB,IPRWFC)
*
          KDUM = 0
          call memmar(KDUM,KDUM,'FLUSM',KDUM,'GASANA')
        end do
      end do
*
      CALL QEXIT('ANACI')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GASANAS_REL(C,LUC,NBL_C,NSSOA,NSSOB,
     &                       NOCTPA,NOCTPB,
     &                       MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                       THRES,MAXTRM,NAEL,NBEL,
     &                       IASTR,IBSTR,IBLTP,NSMST,IUSLAB,
     &                       IOBLAB,NCPMT,WCPMT,NELFSPGP,NOCCLS,
     &                       NGAS,IOCCLS,
     &                       IASPGPTP,IBSPGPTP,
     &                       ICISTR,NORB,IPRNT)
*
* Analyze CI vector :
*
*      1) Print atmost MAXTRM  combinations with coefficients
*         larger than THRES
*         Currently the corresponding dets are not GIVEN !!
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given
*         occupation class
*
* Jeppe Olsen , Jan. 1989 ,
*               Aug 1995 : GAS version
*               May 1997 : BLOCK driven
*
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
*
*
*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION C(*)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      CHARACTER*6 IOBLAB(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBL_C)
*. Output
      DIMENSION NCPMT(10,NOCCLS)
      DIMENSION WCPMT(10,NOCCLS)
*
      NTESTL = 0000
      NTEST = max(IPRNT,NTESTL)
*
      IF(IUSLAB.NE.0) THEN
         WRITE(6,*)
         WRITE(6,*)
     & ' Labels for orbitals are of the type n l ml starting with n = 1'
         WRITE(6,*)
     & ' so the user should not be  alarmed by labels like 1 f+3 '
         WRITE(6,*)
      END IF

      MINPRT = 0
      ITRM = 0
      IDET = 0
      IIDET = 0
      ILOOP = 0
      NCIVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
 2001 CONTINUE
      IF( ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
         XMAX = 1.0D0
         XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
         XMAX = XMIN
         XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNT.EQ.2) THEN
*. Print in one shot
         XMAX = 3006.1956
         XMIN = 0.0D0
      END IF
      IDET = 0
*
      WRITE(6,'(/A,E10.4,A,E10.4/A/)')
     & '  Printout of coefficients in interval  ',XMIN,' to ',XMAX,
     & '  ========================================================='
*
#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
         WRITE(6,'(A,I6)') ' GASANAS : NBL_C = ',NBL_C
         WRITE(6,'(A,I6)') ' GASANAS : IASPGPTP = ',IASPGPTP
         WRITE(6,'(A,I6)') ' GASANAS : IBSPGPTP = ',IBSPGPTP
      end if
#endif
*
      DO JBLOCK = 1, NBL_C
*
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
#if defined LUCI_DEBUG
         if (NTEST.ge.5) write(6,'(A,1X,I6,I6,I6,I6)')
     &        ' IATP,IBTP,IASM,IBSM = ',IATP,IBTP,IASM,IBSM
#endif
*. Obtain alpha strings of sym IASM and type IATP
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,
     &        NASTR1,IASTR,
     &        NORB,0,IDUM,IDUM)
*
*. Obtain beta strings of sym IBSM and type IBTP
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,
     &        NBSTR1,IBSTR,
     &        NORB,0,IDUM,IDUM)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
#if defined LUCI_DEBUG
         if (NTEST.ge.5) write(6,'(A,1X,I9,I9)')
     &        ' NIA,NIB = ',NIA,NIB
#endif
*
         IMZERO = 0
         if (ICISTR.ge.2) then
*. Read in a Type-Type-symmetry block
*? Perhaps use GSTTBL here?
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
#if defined LUCI_DEBUG
            if (NTEST.ge.10) then
              write(6,'(A,I8)') ' Number of elements read in ',IDET
            end if
            if (NTEST.ge.1000) then
              write(6,*)
              write(6,*) 'Partition of C vector is: '
              call wrtmat(C,IDET,1,IDET,1)
            end if
#endif
            IDET = 0
         END IF
C?       write(6,*) ' IMZERO in gasanas:  ',IMZERO
         IF(IMZERO.NE.1) THEN
*
            IBBAS = 1
            IABAS = 1
*
            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO  IA = MINIA,IABAS+NIA-1
*
                  IF (ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
                  IDET = IDET + 1
                  IF (XMAX .GE. ABS(C(IDET)) .AND.
     &                 ABS(C(IDET)).GT. XMIN ) THEN
                     ITRM = ITRM + 1
                     IIDET = IIDET + 1
                     IF( ITRM .LE. MAXTRM ) THEN
                        CNORM = CNORM + C(IDET) ** 2
                        WRITE(6,'(A,I8,A,E14.8/A)')
     &                   '  Coefficient of combination ',IDET,' is ',
     &                       C(IDET),
     &                   '  Corresponding alpha- and beta-strings'
                        IF(IUSLAB.EQ.0) THEN
                           WRITE(6,'("  alpha",20I4)')
     &                          (IASTR(IEL,IA),IEL = 1, NAEL )
                           WRITE(6,'("   beta",20I4)')
     &                          (IBSTR(IEL,IB),IEL = 1, NBEL )
                        ELSE
                           WRITE(6,'("  alpha",10(1X,A6))')
     &                      (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
                           WRITE(6,'("   beta",10(1X,A6))')
     &                      (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
                        END IF
                     END IF
                  END IF
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of if statement for nonvanishing blocks
      END DO
*     ^ End of loop over blocks
      IF(IIDET .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
*
      WRITE(6,'(A,E15.8)')
     &     '  Norm of printed CI vector .. ', CNORM
*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      CALL ISETVC(NCPMT,0    ,10*NOCCLS)
      CALL SETVEC(WCPMT,0.0D0,10*NOCCLS)
      DO JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*. Occupation class corresponding to given occupation
         JOCCLS = 0
         DO JJOCCLS = 1, NOCCLS
            IM_THE_ONE = 1
            DO IGAS = 1, NGAS
               IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &          NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &          IM_THE_ONE = 0
            END DO
            IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
         END DO
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
*
         IMZERO = 0
         IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         END IF
         IF(IMZERO.EQ.0) THEN

            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO IA = MINIA,IABAS+NIA-1
                  IDET = IDET + 1
                  DO IPOT = 1, 10
                     IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
     &                 ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
                       NCPMT(IPOT,JOCCLS)=
     &                       NCPMT(IPOT,JOCCLS)+ 1
                       WCPMT(IPOT,JOCCLS)=
     &                       WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
                     END IF
                  END DO
*     ^ End of loop over powers of ten
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of CI coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
         W = 0.0D0
         N = 0
         DO 290 JOCCLS = 1, NOCCLS
            N = N + NCPMT(IPOT,JOCCLS)
            W = W + WCPMT(IPOT,JOCCLS)
 290     CONTINUE
         WACC = WACC + W
         NACC = NACC + N
         WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 300  CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &     ' IS  ',NCIVAR - NACC
*
      IF(NOCCLS.NE.1) THEN
         WRITE(6,'(A)')
         WRITE(6,'(A)')
     &    '   Magnitude of CI coefficients for each excitation level '
         WRITE(6,'(A)')
     &    '  ========================================================='
         WRITE(6,'(A)')
         DO 400 JOCCLS = 1, NOCCLS
            N = 0
            DO 380 IPOT = 1, 10
               N = N + NCPMT(IPOT,JOCCLS)
 380        CONTINUE
            IF(N .NE. 0 ) THEN
               WRITE(6,*)
               WRITE(6,'(A,15I3)')'       Occupation of active sets :',
     &              (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
               WRITE(6,'(A,I9)')
     &          '         Number of coefficients larger than 10-11 ', N
               WRITE(6,*)
               WACC = 0.0D0
               DO 370 IPOT = 1, 10
                  N =  NCPMT(IPOT,JOCCLS)
                  W =  WCPMT(IPOT,JOCCLS)
                  WACC = WACC + W
                  WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &                 '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 370           CONTINUE
            END IF
 400     CONTINUE
*
*. Total weight and number of dets per excitation level
*
         WRITE(6,'(A)')
         WRITE(6,'(A)')
     &        '   Total weight and number of SD''s (> 10 ** -11 )  : '
         WRITE(6,'(A)')
     &        '  ================================================='
         WRITE(6,'(A)')
       WRITE(6,*) '        N      Weight      Acc. Weight   Occupation '
       WRITE(6,*) ' ==================================================='
         WACC = 0.0D0
         DO 500 JOCCLS = 1, NOCCLS
            N = 0
            W = 0.0D0
            DO 480 IPOT = 1, 10
               N = N + NCPMT(IPOT,JOCCLS)
               W = W + WCPMT(IPOT,JOCCLS)
 480        CONTINUE
            WACC = WACC + W
            IF(N .NE. 0 ) THEN
               WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,16(1X,I2))')
     &              N,W,WACC,(IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
            END IF
 500     CONTINUE
      END IF
*
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C**********************************************************************
C
C     allocate two blocks: VEC1 and VEC2
C     dimension for resolution block is computed:
C                                    - two partitioned blocks of
C                                      length LSCR2 can be stored
C
C      Jeppe Olsen, Jan 1998 - Relativistic version, Feb. 98
C
C      S. Knecht - Jul 2007: memory saving version enabled
C      S. Knecht - Nov 2008: modified for CI property runs
C
C**********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      character*6 CIRUN
C
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "krciprop.h"
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
C
#include "cgas.inc"
#include "lucinp.inc"
#include "cands.inc"
C
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',KDUM,'GET_3B')
C
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDIA)
C
      ISPC = ISSPC
      LSCR2 = 0
C
      IF( CIRUN .ne. 'PROP1 ' ) IXSYMACT = 1
C      
      DO I = 1, IXSYMACT
C
        IJSYM = ISSM
        IF( CIRUN .eq. 'PROP1 ' ) IJSYM = IXSYMLIST(I)
C
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
C
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NAEL = NELEC(IATP)
          NBEL = NELEC(IBTP)
C
          NOCTPA = NOCTYP(IATP)
          NOCTPB = NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
          call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
          call memmar(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
C
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
          CALL ZBLTP_REL(ISMOST(1,IJSYM),NSMST,IDC,WORK(KLCBLTP))
C         largest block of strings in zero order space
          MXSTBL0 = MXNSTR
C         this must be generalized ! alpha and beta strings with an 
C         electron removed
          IATPM1 = ITPM1(1,IMS2)
          IBTPM1 = ITPM1(2,IMS2)
C         alpha and beta strings with two electrons removed
          IATPM2 = ITPM2(1,IMS2)
          IBTPM2 = ITPM2(2,IMS2)
C         largest number of strings of given symmetry and type
          MAXA = MXNSTR
          IF(NAEL.GE.1) THEN
            MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),
     &                        2)
            MAXA = MAX(MAXA,MAXA1)
          END IF
          IF(NAEL.GE.2) THEN
            MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),
     &                        2)
            MAXA = MAX(MAXA,MAXA1)
          END IF
          MAXB = 0
          IF(NBEL.GE.1) THEN
            MAXB1 =IMNMX_REL(WORK(KNSTSO2(IBTPM1)),NSMST*NOCTYP(IBTPM1),
     &                       2)
            MAXB = MAX(MAXB,MAXB1)
          END IF
          IF(NBEL.GE.2) THEN
            MAXB1 =IMNMX_REL(WORK(KNSTSO2(IBTPM2)),NSMST*NOCTYP(IBTPM2),
     &                       2)
            MAXB = MAX(MAXB,MAXB1)
          END IF
          MXSTBL = MAX(MAXA,MAXB)
C
C         largest number of resolution strings and spectator strings
C         that can be treated simultaneously
          MAXK = MIN(MXINKA,MXSTBL)
C
C         scratch space for CJKAIB resolution matrices
C         size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
          CALL MXRESC_REL(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                    NSMST,NSTFSMSPGP,NSTFSMSPGP2,MXPNSMST,NSMOB,
     &                    MXPNGAS,NGAS,NOBPTS,NOBPTS2,IPRCIX,MAXK,
     &                    NELFSPGP,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,
     &                    MXADKBLK)
          LSCR2 = MAX(LSCR2,MXCJ,MXCIJA,MXCIJB)
        END DO
C       ^ end loop over MK2 values
      END DO
C     ^ end loop over active symmetry irreps
C       
C     L2BLOCK and LBLOCK are stored in crun.inc; dimensioning in Z_BLK_*
      L0BLOCK = MXSOOB
CSK...6-jul-07      LBLOCK = MAX(L2BLOCK,L0BLOCK)
      L3BLOCK = MAX(L2BLOCK,L0BLOCK)
CSK...6-jul-07 new
      L3BLOCK = MIN(LBLOCK,L3BLOCK)
CSK...6-jul-07      LSCR12 = MAX(LBLOCK,2*LSCR2)
      LSCR12 = MAX(LBLOCK,2*LSCR2)
      LCJRES_SAVE = LSCR12
C     
#if defined LUCI_DEBUG
      write(luwrt,*) 'L3BLOCK is',L3BLOCK
      write(luwrt,*) 'LBLOCK is',LBLOCK
      write(luwrt,*) 'LSCR2 is',LSCR2
      write(luwrt,*) 'LSCR12 is',LSCR12
      write(luwrt,*) 'LCJRES_SAVE is',LCJRES_SAVE
#endif
C
      IDUM = 0
      call memmar(KDUM,IDUM,'FLUSM ',KDUM,'GET_3B')
C
      WRITE(LUWRT,'(/3X,A)') 
     & '=============================================================='
      WRITE(LUWRT,'(3X,A)') 
     & '==> allocation of two CI vectors and one resolution vector <=='
      write(LUWRT,*) ' '
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'current available free memory in double words:',LMEMFREE_PTR
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'allocate two CI vectors each of length:       ',LBLOCK
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'allocate resolution vector of length:         ',LCJRES_SAVE
      WRITE(LUWRT,'(3X,A/)') 
     & '=============================================================='
C
C     static allocation
      call memmar(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      call memmar(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_SPGP_INF_REL(ISPGP,ITP,IGRP)
*
* Obtain groups defining supergroup ISPGP,ITP
*
* Jeppe Olsen, November 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Output
      DIMENSION IGRP(*)
*
      NTEST = 000
*. Absolute group number
      ISPGPABS = ISPGP + IBSPGPFTP(ITP) -1
      CALL ICOPY(NGAS,ISPGPFTP(1,ISPGPABS),1,IGRP,1)
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_SPGP_INF_REL : ISPGP ITP ISPGPABS',
     &              ISPGP, ITP, ISPGPABS
        WRITE(6,*) ' Groups of supergroups'
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  This subroutine generates all group info needed for
*   - integral transformation
*   - CI
*   - Spinor rotation
*
      subroutine grinfo(PNTGRP)
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      integer, intent(inout) :: PNTGRP
*
#include "mxpdim.inc"
#include "symm.inc"
! stefan: test for spinfree run in D2h symmetry
#include "dcbham.h"
*
      if (PNTGRP.eq.9) DOUGRP = 4
      if (PNTGRP.eq.9.and.spinfr) DOUGRP = 9
      if (PNTGRP.eq.8) DOUGRP = 5
      if (PNTGRP.eq.7) DOUGRP = 6
      if (PNTGRP.eq.6) DOUGRP = 7
      if (PNTGRP.eq.5) DOUGRP = 8
*
      if (DOUGRP.eq.1) then
         NIRR_DG = 10
         NIRR_PN = 8
         IHASDG = 4
         PNTGRP = 1
      else if (DOUGRP.eq.2) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.3) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.4) then
         NIRR_DG = 8
         NIRR_PN = 4
         IHASDG = DOUGRP
         PNTGRP = 9
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         if (DOUGRP.eq.5) PNTGRP = 8
         if (DOUGRP.eq.6) PNTGRP = 7
      else if (DOUGRP.eq.7) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         PNTGRP = 6
      else if (DOUGRP.eq.8) then
         NIRR_DG = 2
         NIRR_PN = 1
         IHASDG = DOUGRP
         PNTGRP = 5
      else if (DOUGRP.eq.9) then
         NIRR_DG = 12
         NIRR_PN = 8
         IHASDG = DOUGRP
         PNTGRP = 9
         print *, 
     & '*** error in grinfo: no spinfree formalism yet in LUCIAREL. ***'
         call quit(
     & '*** error in grinfo: no spinfree formalism in LUCIAREL.!***')
      else
         write(*,*) 'Double group out of range!'
         write(*,*) 'Double group is  ',DOUGRP
         write(*,*) 'Picasso will put down his brush.'
         call quit('*** error in grinfo: Double group out of range!***')
      end if
*
*  Set up info for highest abelian double subgroup
*
      if (DOUGRP.ne.IHASDG) then
         if (DOUGRP.eq.1) then
            NSDGIRR = 8
            NSPGIRR = 4
         else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
            NSDGIRR = 4
            NSPGIRR = 2
         else if (DOUGRP.eq.9) then
            NSDGIRR = 12
            NSPGIRR = 8
         end if
      end if
 
!     stefan: generalization of dimensions
      MXPDBGDIM  = max(8,NIRR_DG)
      MXPPGPDIM  = max(4,NIRR_PN)
      MXPDBGDIM2 = MXPDBGDIM ** 2
   
!     initialize
      CALL IZERO(IDBGMULT,MXPDBGDIM2)
      CALL IZERO(IMOSF_SP, MXPDBGDIM)
      CALL IZERO(ISPSF_MO, MXPDBGDIM)

!     stefan: adapt here for straightforward implementation of linear symmetry 
!             and spinfree LUCIAREL. SK - Sep 2010
      MXPDBGDIM  = max(8,NIRR_DG)
      IF(PNTGRP.EQ.5) THEN
*. C1 point group
        DBGSYM(1:6) ='C1    '
*
        IPGMULT(1,1) = 1
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
*. All spinors have symmetry 2, all orbitals symmetry 1
        IMOSF_SP(1,1) = 2
        IMOSF_SP(1,2) = 2
*
        ISPSF_MO(2,1) = 1
        ISPSF_MO(2,2) = 1
*. Spinors and their cc have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
      ELSE IF(PNTGRP.EQ.6 ) THEN
*. Ci pointgroup
        DBGSYM(1:6) ='Ci    '
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 1
        IDBGMULT(3,4) = 2
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 2
        IDBGMULT(4,4) = 1
*. MO have symmetries 1,2 spinors symmetries 3,4 alpha,beta have symm 3
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 3
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 4
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 1
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 2
*. Functions and adjoints have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
      ELSE IF(PNTGRP.EQ.7.OR.PNTGRP.EQ.8) THEN
*. CS and C2 point groups
        IF(PNTGRP.EQ.7) THEN
          DBGSYM(1:6) ='Cs    '
        ELSE
          DBGSYM(1:6) ='C2    '
        END IF
*
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 2
        IDBGMULT(3,4) = 1
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 1
        IDBGMULT(4,4) = 2
*. Orbitals have sym 1,2, spinors 3,4, alpha 3, beta 4
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 4
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 3
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 2
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 1
*. Adjoints
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 4
        IADJSYM(4) = 3
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 4
        INVELM(4) = 3
      ELSE IF(PNTGRP.EQ.9) THEN
*. C2H
        DBGSYM(1:6) ='C2h   '
*
*. Adjoints
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
        IADJSYM(5) = 6
        IADJSYM(6) = 5
        IADJSYM(7) = 8
        IADJSYM(8) = 7
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
        INVELM(5) = 6
        INVELM(6) = 5
        INVELM(7) = 8
        INVELM(8) = 7
*
        ISPSF_MO(5,1) = 1
        ISPSF_MO(6,2) = 1
        ISPSF_MO(6,1) = 2
        ISPSF_MO(5,2) = 2
        ISPSF_MO(7,1) = 3
        ISPSF_MO(8,2) = 3
        ISPSF_MO(8,1) = 4
        ISPSF_MO(7,2) = 4
*
        IMOSF_SP(1,1) = 5
        IMOSF_SP(1,2) = 6
        IMOSF_SP(2,1) = 6
        IMOSF_SP(2,2) = 5
        IMOSF_SP(3,1) = 7
        IMOSF_SP(3,2) = 8
        IMOSF_SP(4,1) = 8
        IMOSF_SP(4,2) = 7
*
        do JJ=1,MXPDBGDIM,1
           do II=1,MXPDBGDIM,1
              if (II.le.4.and.JJ.le.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 1
                 else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.5) then
                    IDBGMULT(II,JJ) = 4
                 else
                    IDBGMULT(II,JJ) = 3
                 end if
              else if (II.gt.4.and.JJ.gt.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.13) then
                    IDBGMULT(II,JJ) = 3
                 else if ((II+JJ).eq.11.or.(II+JJ).eq.15) then
                    IDBGMULT(II,JJ) = 1
                 else
                    IDBGMULT(II,JJ) = 4
                 end if
              else
                 if ((II+JJ).eq.9) then
                    IDBGMULT(II,JJ) = 8
                 else if ((II+JJ).eq.7.or.(II+JJ).eq.11) then
                    IDBGMULT(II,JJ) = 6
                 else if ((II+JJ).eq.6.or.(II+JJ).eq.12) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.6.and.JJ.eq.2) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.2.and.JJ.eq.6) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.3.and.JJ.eq.7) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.7.and.JJ.eq.3) then
                    IDBGMULT(II,JJ) = 5
                 else
                    IDBGMULT(II,JJ) = 7
                 end if
              end if
           end do
        end do
        do JJ=1,MXPPGPDIM,1
           do II=1,MXPPGPDIM,1
              if (II.eq.JJ) then
                 IPGMULT(II,JJ) = 1
              else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                 IPGMULT(II,JJ) = 2
              else if ((II+JJ).eq.5) then
                 IPGMULT(II,JJ) = 4
              else
                 IPGMULT(II,JJ) = 3
              end if
           end do
        end do
      END IF
*
CSK   NTEST = 100
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Symmetry Information from GRINFO'
        WRITE(6,*) ' ================================'
        WRITE(6,*)
        WRITE(6,'(A,A)') ' Double group : ', DBGSYM
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of irreps ',NIRR_DG
        WRITE(6,*)
        WRITE(6,*) ' Multiplication table '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        DO IREP = 1, NIRR_DG
          WRITE(6,'(10X,8I3)') (IDBGMULT(IREP,JREP),JREP=1,NIRR_DG)
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' orbital X spin => spinor '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO IOBSM = 1, NIRR_DG/2
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    IOBSM, ISF, IMOSF_SP(IOBSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Spinor  X spin => orbital'
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO ISPSM =  NIRR_DG/2+1, NIRR_DG
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    ISPSM, ISF, ISPSF_MO(ISPSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Adjoint symmetry array '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,'(10X,8I3)') (IADJSYM(IREP),IREP=1,NIRR_DG)
*
        WRITE(6,*) ' Inverse elements '
        WRITE(6,*) ' ================='
        WRITE(6,'(10X,8I3)') (INVELM(IREP),IREP=1,NIRR_DG)
      END IF
#endif
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE IAIBCM_REL(ICISPC,IATP_TF,IBTP_TF,IAIB)
*
* obtain allowed combination of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*
*  Corrected for rel. GASCI
*       Timo Fleig, 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Output
      INTEGER IAIB(*)
*
      NAEL = NELEC(IATP_TF)
      NBEL = NELEC(IBTP_TF)
*
      NOCTPA = NOCTYP(IATP_TF)
      NOCTPB = NOCTYP(IBTP_TF)
*
      IOCTPA = IBSPGPFTP(IATP_TF)
      IOCTPB = IBSPGPFTP(IBTP_TF)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &                IGSOCCX,NOCTPA,
     &                NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NELFGP,MXPNGAS,NGAS,IAIB,IPRDIA)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTDIM_REL(IPRNT)
*
* Number of integrals and storage mode
*
*. Last modifications : Jan 98, relativistic integrals added
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "crun.inc"
#include "cgas.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),
     &              SXDXSX(2*MXPOBS,4*MXPOBS)
*. Local arrays
      INTEGER NNORB(MXNDGIRR,4)
*
* =======
*. Output
* =======
*. Atmost 10 double group integral lists
C     COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2,
C    &              NL1D,IUB1D(2,10),L1D(10),
C    &              NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10)
#include "cintfo.inc"
C     initial definitions
#if defined (VAR_MPI2)
      MASTER = MPARID
      MYPROC = MYTID
*     Add the master node
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
#if defined LUCI_DEBUG
      if (NTEST.ge.2)
     &    WRITE(6,'(A25,I3)') ' INTDIM speaking, IDBG = ', IDBG
#endif
*
* =================================
* integrals in point group symmetry
* =================================
*
      IPNTGRP_SAVE = PNTGRP
*. Force the use of point group symmetry
      IF(IDBG.EQ.1) PNTGRP = 1
*.1 : Number of one-electron integrals
C     NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
      NINT1 =  NSXFSM2(NSMSH,NGSSHT,NGSSHT,ITSSX,1,IPRNT)
C              NSXFSM2(NSMOB,NO1PS,NO2PS,ISXSM,ISYM,IPRNT)
*.2 : Number of two-electron integrals
*. Full eightfold symmetry can be used
      I12S = 1
      I34S = 1
      I1234S = 1
      NINT2 = NDXFSM2(NSMSH,NGSSHT,NGSSHT,NGSSHT,NGSSHT,
     &                I12S,I34S,I1234S,IPRNT )
C     NDXFSM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
C    &         IS12,IS34,IS1234,IPRNT)
*
* =======================
* Double group integrals
* =======================
*
      IF(IDBG.EQ.0) THEN
        NL1D = 0
        NL2D = 0
      ELSE
        PNTGRP = IPNTGRP_SAVE
*. Number of double group integral lists
        CALL DBG_INT_LISTS(NL1D,IUB1D,NL2D,IUB2D,ISM2D,IDENSI,
     &                     ISPINFREE)
C            DBG_INT_LISTS(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
C     &                    ISPINFREE)
*
        DO IL1 = 1, NL1D
          DO INDEX = 1, 2
            IF(IUB1D(INDEX,IL1).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L1D(IL1) = NSXFSM2(NSMOB,NNORB(1,1),NNORB(1,2),ITSSX,1,IPRNT)
        END DO
*
        DO IL2 = 1, NL2D
          DO INDEX = 1, 4
            IF(IUB2D(INDEX,IL2).EQ.1) THEN
              CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
            ELSE
              CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
            END IF
          END DO
          L2D(IL2) = NDXFSM2(NSMOB,
     &               NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &               ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &               IPRNT )
        END DO
      END IF
*. Largest block of two-electron integrals
      LMAXE2LST = NINT2
      DO IL2 = 1, NL2D
        LMAXE2LST = MAX(LMAXE2LST,L2D(IL2))
      END DO
*
      IF( MYPROC .eq. MASTER .and. NTEST .gt. 0) THEN
        write(6,'(/A,I10)')
     &   '   The largest block of two-electron integrals is',LMAXE2LST
      END IF
*
*
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTPNT_DBG
*
* Pointers to symmetry blocks of integrals for DOUBLE GROUP
*
* Jeppe Olsen, Jan 98

      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
C     COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2,
C    &              NL1D,IUB1D(2,10),L1D(10),
C    &              NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10),
C    &              LMAXE2LST
*
* =======
*. Output
* =======
*
*. Local scratch
      DIMENSION NNORB(MXNDGIRR,4)
*
*.1 : Offset for one-electron integrals
      DO IL1 = 1, NL1D
        DO INDEX = 1,2,1
         IF(IUB1D(INDEX,IL1).EQ.1) THEN
           CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
         ELSE
           CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
         END IF
        END DO
*
        CALL PNT2DM2(1,NSMOB,NNORB(1,1),NNORB(1,2),ITSSX,
     &       WORK(KPINT1D(IL1)) )
C            PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)
      END DO
*
*.2 : two-electron integrals
      DO IL2 = 1, NL2D
        DO INDEX = 1,4,1
         IF(IUB2D(INDEX,IL2).EQ.1) THEN
           CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
         ELSE
           CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
         END IF
        END DO
*
        CALL PNT4DM2(NSMOB,
     &       NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &       ITSDX,ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &       WORK(KPINT2D(IL2)))
      END DO
*
      CALL RMEMCHK
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTPNT_PNG(IPNT1,IPNT2)
*
* Pointers to symmetry blocks of integrals for POINT GROUP
* IPNT1 : Pointer to given one-electron block
* IPNT2 : Pointer to given two-electron block
*
* Note : Generation of symmetry of last index (ISL1, ISL2) eliminated

      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "cgas.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
* =======
*. Output
* =======
*
      INTEGER IPNT1(NSMOB)
      INTEGER IPNT2(NSMOB,NSMOB,NSMOB)
*
* if the actual group is a double group, hide this for the moment
* and use the corresponding pointgroup
      IPNTGRP_SAVE = PNTGRP
      IF(IDBG.EQ.1) THEN
        PNTGRP = 1
      END IF
*
*.1 : Number of one-electron integrals
      CALL PNT2DM2(1,NSMSH,NGSSHT,NGSSHT,ITSSX,IPNT1)
C          PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)

*.2 : two-electron integrals
      I12S = 1
      I34S = 1
      I1234S = 1
      CALL PNT4DM2(NSMSH,NGSSHT,NGSSHT,NGSSHT,NGSSHT,
     &            ITSDX,I12S,I34S,I1234S,IPNT2)
C          PNT4DM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
C    &           IDXSM,IS12,IS34,IS1234,IPNTR)

      PNTGRP = IPNTGRP_SAVE
*
      CALL RMEMCHK
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION IOFF_DBG_LIST(ILIST,ISYM,JSYM,KSYM,LSYM)
*
* Offset to given integral block of double group integrals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "lucinp.inc"
*
      IOFF_DBG_LIST =
     &IOFF_DBG_LIST_S(WORK(KPINT2D(ILIST)),ISYM,JSYM,KSYM,LSYM,
     &                NSMOB)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION IOFF_DBG_LIST_S(IPNT,ISYM,JSYM,KSYM,LSYM,NSMOB)
*
* Read offset to given symmetry block from IPNT
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IPNT(NSMOB,NSMOB,NSMOB)
*
      IOFF_DBG_LIST_S = IPNT(ISYM,JSYM,KSYM)
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE KRA_SYM(ISYM,IBU,IKRSM,IBUK)
* A spinor symmetry ISYM is given for a spinor, Barred/Unbarred
* given by IBU
*
* Find corresponding symmetry of spinor with opposite Barred/Unbarred symmetry
* as well as opposite Kramers index IBUK
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
*
      IPNTGP_SYM = ISPSF_MO(ISYM,IBU)
      IF(IBU.EQ.1) THEN
        IKRSM = IMOSF_SP(IPNTGP_SYM,2)
        IBUK = 2
      ELSE
        IKRSM = IMOSF_SP(IPNTGP_SYM,1)
        IBUK = 1
      END IF
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Kramers connected symmetries',ISYM,IKRSM
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE LCISPC_REL(LCSBLK,IPRNT)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               July 1997, multiple Ms spaces allowed
*
* GAS VERSION
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* ===================
*.Input common blocks
* ===================
*
*./LUCINP/
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "strbas.inc"
*./CSM/ : NSMST is used
#include "csm.inc"
#include "stinf.inc"
*./CSTATE/ : IDC is used
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRNT)
*
      CALL QENTER('LCISP')
*. Number of spaces
      NICISP = NCMBSPC
*
*.    set mark for local memory
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'LCISP ')
*
*.  Largest NOCTPA*NOCTPB
      MXNOCAB = 0
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
        NOCAB = NOCTYP(IATP)*NOCTYP(IBTP)
        MXNOCAB = MAX(MXNOCAB,NOCAB)
      END DO
C     WRITE(6,*) ' MXNOCAB = ',MXNOCAB
      call memmar(KLIOIO,MXNOCAB,   'ADDL  ',2,'KLIOIO')
*
      call memmar(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP')
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &  call memmar(KLCVST,NSMST,'ADDL  ',2,'KLCVST')
*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
      CALL SMOST_REL(NSMST,NSMCI,MXPCSM,ISMOST)
*
      MXSB = 0
      MXSOOB = 0
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*
C?      WRITE(6,*) ' IMS2 IATP IBTP NOCTPA NOCTPB '
C?      WRITE(6,*)   IMS2,IATP,IBTP,NOCTPA,NOCTPB
*
        DO ICI = 1, NICISP
*. allowed combination of types
          CALL IAIBCM_GAS(LCMBSPC(ICI),ICMBSPC(1,ICI),
     &                    IGSOCCX,NOCTPA,
     &                    NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                    NELFGP,MXPNGAS,NGAS,WORK(KLIOIO),IPRNT)
*
          KLCVST = 1
          DO  ISYM = 1, NSMCI
            CALL ZBLTP_REL(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP))
            CALL NGASDT_REL(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),
     &                      NGAS,ISYM,NSMST,NOCTPA,NOCTPB,
     &                      WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &                      ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                      ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &                      LCOL,WORK(KLIOIO))
*
            IF(IMS2.EQ.1) THEN
              XISPSM(ISYM,ICI) = XNCOMB
              NBLKIC(ISYM,ICI) = NTTSBL
              LCOLIC(ISYM,ICI) = LCOL
            ELSE
              XISPSM(ISYM,ICI) = XISPSM(ISYM,ICI) + XNCOMB
              NBLKIC(ISYM,ICI) = NBLKIC(ISYM,ICI) + NTTSBL
              LCOLIC(ISYM,ICI) = LCOLIC(ISYM,ICI) + LCOL
            END IF
            MXSOOB = MAX(MXSOOB,MXSOO)
            MXSB = MAX(MXSB,MXS)
          END  DO
*         ^ End of loop over symmetries
        END DO
*       ^ End of loop over CI spaces
      END DO
*     ^ End of loop over MS2 spaces
*
*     Setting size of resolution matrices
      LCSBLK = IMAXLBLKSZ
*      
#if defined LUCI_DEBUG
      IF(NTEST.ge.1) write(LUWRT,'(/A,I18/)') ' Set LCSBLK to ',LCSBLK
#endif
*
*     Check for minimum size of resolution matrices:
      if (MXSB.gt.LCSBLK) then
        write(LUWRT,'(/A,I18/)') ' WARNING: Resetting LCSBLK to: ',MXSB
        LCSBLK = MXSB
      end if
*
      NTEST = 0
*
#if defined LUCI_DEBUG
      if (NTEST.ge.1) then
        WRITE(LUWRT,*)
     &  ' Number of internal combinations per symmetry '
        WRITE(LUWRT,*)
     &  ' =========================================== '
*
        DO 200 ICI = 1, NCMBSPC
          WRITE(LUWRT,*) ' CI space ', ICI
          CALL WRTMATMN(XISPSM(1,ICI),1,NSMCI,1,NSMCI,LUWRT)
  200   CONTINUE
        WRITE(LUWRT,*) ' Largest symmetry block           ',MXSB
        WRITE(LUWRT,*) ' Largest Symmetry-type-type block ',MXSOOB
*
        WRITE(LUWRT,*) ' Number of TTS subblocks per CI expansion '
        WRITE(LUWRT,*) ' ======================================== '
*
        DO ICI = 1,  NCMBSPC
          WRITE(LUWRT,*) ' Internal CI space ', ICI
          CALL IWRTMAMN(NBLKIC(1,ICI),1,NSMCI,1,NSMCI,LUWRT)
        END DO
      end if
#endif
*
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(LUWRT,*) ' Largest number of blocks in CI expansion',
     &                   MXNTTS
        WRITE(LUWRT,*) ' Number of columns per CI expansion '
        WRITE(LUWRT,*) ' =================================== '
*
        DO ICI = 1,NCMBSPC
          WRITE(LUWRT,*) ' Internal CI space ', ICI
          CALL IWRTMAMN(LCOLIC(1,ICI),1,NSMCI,1,NSMCI,LUWRT)
        END DO
      END IF
#endif
*
C
C     eliminate local memory
      IDUM = 0
      call memmar(KDUM ,IDUM,'FLUSM ',2,'LCISP ')

      CALL QEXIT('LCISP')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MINMAX_FOR_SYM_DIST_DBG(NIGRP,IGRP,MNVAL,MXVAL,NDIST,
     &                                   IAB)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Include blocks     
#include "mxpdim.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "wrkspc.inc"
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
#endif
*
      DO JGRP = 1, NIGRP
C Requires setting up of  MINMAX_SM_GP and  MINMAX_SM_GP2 
        IF(IAB.EQ.1) THEN
          MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
        ELSE
          MNVAL(JGRP) = MINMAX_SM_GP2(1,IGRP(JGRP))
          MXVAL(JGRP) = MINMAX_SM_GP2(2,IGRP(JGRP))
        END IF
      END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
*
* Transfer between ML,IPARI notation and compound notation ISM
*
* IWAY = 1 : IML,IPARI => Compound
* IWAY = 2 : IML,IPARI <= Compound
*
* TYPE : 'SX','OB','ST','DX','CI'
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
      CHARACTER*2 TYPE
*./NONAB/
      LOGICAL INVCNT
C     PARAMETER (MXPOBS = 20 )
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*./CSM/
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
*.(Tired of warnings from 3090 Compiler so )
* (
      NML = 0
      MXML= 0
      MNML= 0
*             )
      IF(TYPE.EQ.'OB') THEN
        NML = NMLOB
        MXML = MXMLOB
        MNML = MNMLOB
      ELSE IF(TYPE.EQ.'SX') THEN
        NML = NMLSX
        MXML = MXMLSX
        MNML = MNMLSX
      ELSE IF(TYPE.EQ.'DX') THEN
        NML = NMLDX
        MXML = MXMLDX
        MNML = MNMLDX
      ELSE IF(TYPE.EQ.'ST') THEN
        NML = NMLST
        MXML = MXMLST
        MNML = MNMLST
      ELSE IF(TYPE.EQ.'CI') THEN
        NML = NMLCI
        MXML = MXMLCI
        MNML = MNMLCI
      END IF
*
      IF(IWAY.EQ.1) THEN
C        ISM = (IPARI-1)*NML + MNML - 1
         ISM = (IPARI-1)*NML + IML - MNML + 1
      ELSE IF(IWAY.EQ.2) THEN
        IF(ISM.GT.NML) THEN
          IPARI = 2
          IML = ISM - NML + MNML - 1
        ELSE
          IPARI = 1
          IML = ISM       + MNML - 1
        END IF
      ELSE
        WRITE(6,*) ' Error in MLSM , IWAY = ' ,IWAY
        WRITE(6,*) ' MLSM stop !!! '
        Call Abend1( 20 )
      END IF
*
#if defined LUCI_DEBUG
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,A)') ' MLSM speaking ,type= ',TYPE
        WRITE(6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NATORB_LCR(RHO1,NSMOB,NTOPSM,NACPSM,NINPSM,
     &                      ISTOB,XNAT,RHO1SM,OCCNUM,
     &                      NACOB,SCR,IPRDEN)
*
* Obtain natural orbitals in symmetry blocks
*
* Jeppe Olsen, June 1994
*              Modification, Oct 94
*              Last modification, Feb. 1998 (reorder deg eigenvalues)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION RHO1(NACOB,NACOB)
      DIMENSION ISTOB(*)
      DIMENSION NTOPSM(NSMOB), NACPSM(NSMOB),NINPSM(NSMOB)
*. Output
      DIMENSION RHO1SM(*),OCCNUM(*),XNAT(*)
*. Scratch ( Largest symmetry block )
      DIMENSION SCR(*)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRDEN)
*. IOBOFF : Offset for active orbitals in symmetry order
      DO ISMOB = 1, NSMOB
        IF(ISMOB.EQ.1) THEN
          IOBOFF = NINPSM(1)+1
          IMTOFF = 1
        ELSE
          IOBOFF =
     &    IOBOFF + NTOPSM(ISMOB-1)-NINPSM(ISMOB-1)+NINPSM(ISMOB)
          IMTOFF = IMTOFF + NACPSM(ISMOB-1)**2
        END IF
        LOB = NACPSM(ISMOB)
*
*. Extract symmetry block of density matrix
*
        DO IOB = IOBOFF,IOBOFF + LOB-1
           DO JOB = IOBOFF,IOBOFF + LOB-1
*. Corresponding type indeces
             IOBP = ISTOB(IOB)
             JOBP = ISTOB(JOB)
             RHO1SM(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1)
     &     = RHO1(IOBP,JOBP)
           END DO
        END DO
*
        IF(NTEST.GE.2 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' Density matrix for symmetry  = ', ISMOB
          WRITE(6,*) ' ======================================='
          WRITE(6,*)
          CALL WRTMAT(RHO1SM(IMTOFF),LOB,LOB,LOB,LOB)
        END IF
*. Pack and diagonalize
        SIGNTP = 1.0
        CALL TRIPAK(RHO1SM(IMTOFF),SCR,1,LOB,LOB,SIGNTP)
        ONEM = -1.0D0
*. scale with -1 to get highest occupation numbers as first eigenvectors
        CALL SCALVE(SCR,ONEM,LOB*(LOB+1)/2)
        CALL EIGEN(SCR,XNAT(IMTOFF),LOB,0,1)
*
        DO  I = 1, LOB
          OCCNUM(IOBOFF-1+I) = - SCR(I*(I+1)/2)
        END DO
*. Order the degenerate eigenvalues so diagonal terms are maximized
        TESTY = 1.0D-11
        DO IOB = 2, LOB
          IF(ABS(OCCNUM(IOBOFF-1+IOB)-OCCNUM(IOBOFF-2+IOB))
     &       .LE.TESTY) THEN
            XII   = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB  ))
            XI1I1 = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB-1))
            XII1  = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB  ))
            XI1I  = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB-1))
*
            IF( XI1I.GT.XII.AND.XII1.GT.XI1I1 ) THEN
*. interchange orbital IOB and IOB -1
              CALL SWAPVE(XNAT(IMTOFF+(IOB-1)*LOB),
     &                    XNAT(IMTOFF+(IOB-1-1)*LOB),LOB)
              SS = OCCNUM(IOBOFF-1+IOB-1)
              OCCNUM(IOBOFF-1+IOB-1) = OCCNUM(IOBOFF-1+IOB)
              OCCNUM(IOBOFF-1+IOB)   = SS
              write(6,*) ' Orbitals interchanged ',
     &        IOBOFF-1+IOB,IOBOFF-2+IOB
            END IF
          END IF
        END DO
*
        IF(NTEST.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*)
     &    ' Natural occupation numbers for symmetry = ', ISMOB
          WRITE(6,*)
     &    ' ==================================================='
          WRITE(6,*)
          CALL WRTMAT(OCCNUM(IOBOFF),1,LOB,1,LOB)
          IF(NTEST.GE.2 ) THEN
            WRITE(6,*)
            WRITE(6,*) ' Corresponding Eigenvectors '
            WRITE(6,*)
            CALL WRTMAT(XNAT(IMTOFF),LOB,LOB,LOB,LOB)
          END IF
        END IF
      END DO
*. ( End of loop over orbital symmetries )
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION NDXFSM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
     &         IS12,IS34,IS1234,IPRNT)
*
* Number of double excitations with total symmetry IDXSM
*
* IS12 (0,1,-1)   : Permutational symmetry between index 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between index 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between index 12 and 34
*
* Compared to NDXFSM : Reference to symmetry arrays eliminated by
*                      calls to SYMCOM ( easier when several groups are
*                      used )
*. General input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*
      MDX = 0
      IDXSM = 1
      DO 200 I12SM = 1, NSMOB
        DO 190 I1SM = 1, NSMOB
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL SYMCOM_REL(2,1,I1SM,I2SMA,I12SM)
          I2SM = I_GET_ADJSYM(I2SMA)
C         I2SM = ADSXA(I1SM,I12SM)
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 190
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMSX+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
C              I34SM = SXDXSX(I12SM,IDXSM)
C              SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL SYMCOM_REL(2,1,I12SM,I34SM,IDXSM)
          DO 90 I3SM = 1, NSMOB
C           I4SM = ADSXA(I3SM,I34SM)
            CALL SYMCOM_REL(2,1,I3SM,I4SMA,I34SM)
            I4SM = I_GET_ADJSYM(I4SMA)
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 90
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMSX+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 90
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              MDX = MDX + N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12+1)/2
              ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12-1)/2
            END IF
C?          WRITE(6,*) ' I1SM I2SM I3SM I4SM MDX '
C?          WRITE(6,*)   I1SM,I2SM,I3SM,I4SM,MDX
   90       CONTINUE
C 100     CONTINUE
  190   CONTINUE
  200 CONTINUE
*
      NDXFSM2 = MDX
*
      NTEST = 0
#if defined LUCI_DEBUG
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.NE.0) THEN
         WRITE(6,'(A39,I10)')
     &   ' Number of double excitations obtained ', MDX
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NEWTYP_REL(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*
* an input supergroup is given.
* apply a string of elementary operators to this supergroup and
* obtain supergroup number of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95
*
* ------
* Input
* ------
*
* INSPGP  : input super group ( given occupation in each AS )
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER ITPOP(*),IACOP(*)
*. Number of active spaces  (NGAS )
#include "mxpdim.inc"
#include "cgas.inc"
#include "strbas.inc"
#include "wrkspc.inc"
*. Local scratch
      DIMENSION IEL(MXPNGAS)
*. output
      INTEGER OUTSPGP
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.1) THEN
        CALL NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &       NGAS,WORK(KSPGPAN),WORK(KSPGPCR),OUTSPGP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
      ELSE IF(INEW_OR_OLD.EQ.2) THEN
*. Number of electrons in AS1,AS2, .... for input supergroup
        CALL GTSPGP_REL(IEL,INSPGP,2)
*
        IDELTA = 0
        DO IOP = 1, NOP
*. Change in number of orbitals
          IF(IACOP(IOP).EQ.1) THEN
            IDELTA = IDELTA - 1
            IDELTAI = -1
          ELSE
            IDELTA = IDELTA + 1
            IDELTAI = +1
          END IF
          IEL(ITPOP(IOP)) = IEL(ITPOP(IOP)) + IDELTAI
        END DO
*. output supergroup
        CALL GTSPGP_REL(IEL,OUTSPGP,1)
      END IF
*     /\ End  of NEW_OR_OLD switch
*
#if defined LUCI_DEBUG
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP ,  OUTSPGP ', OUTSPGP
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GTSPGP_REL(IEL,ISPGP,IWAY)
*
*
* Relation between number of electrons in AS1, AS2 ... and
* supergoup number
*
* IWAY = 1 :
* Get ISPGP : Supergroup of strings
*             with IEL(*)  electrons in each AS
* IWAY = 2 :
* GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP
*
*
* Jeppe Olsen, Another lonely night in Lund
*               GAS version July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Generel input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. input(IWAY = 2 ), output (IWAY = 1 )
      INTEGER IEL(*)
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Number
        ISPGP = -1
        DO JSPGP = 1, NTSPGP
          IF(ISPGP.EQ.-1) THEN
            IEQUAL = 1
            DO IGAS = 1, NGAS
              IF(NELFSPGP(IGAS,JSPGP).NE.IEL(IGAS))  IEQUAL= 0
            END DO
            IF(IEQUAL.EQ.1) ISPGP = JSPGP
          END IF
        END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Number => Occupation
        DO IGAS = 1, NGAS
         IEL(IGAS) = NELFSPGP(IGAS,ISPGP)
        END DO
      END IF
*
#if defined LUCI_DEBUG
      NTEST  = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from GTSPGP_REL '
        WRITE(6,*)
     &   ' IWAY ISPGP IEL ', IWAY,ISPGP,(IEL(IGAS),IGAS = 1, NGAS)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NEXT_SYM_DISTR_REL(NGAS,MINVAL,MAXVAL,
     &                              ISYM,ISYM_TOT,IFIRST,NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*. Symmetry of first NGAS -1 spaces
C     write(6,*) ' next_sym...:nonew,ifirst',nonew,ifirst
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS-1
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        IFIRST = 0
        NONEW = 0
      ELSE
C            NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
      END IF
      KSYM = -1
*. Symmetry of last space
      IF(NONEW.EQ.0) THEN
        JSYM = 1
CSK     write(6,*)'before SYMCOM_REL #1: NGAS-1 is',NGAS-1
        DO IGAS = 1, NGAS-1
          CALL SYMCOM_REL(3,0,JSYM,ISYM(IGAS),KSYM)
          JSYM = KSYM
        END DO
        CALL SYMCOM_REL(2,0,JSYM,ISYM(NGAS),ISYM_TOT)
C      write(6,*) ' after SYMCOM : JSYM, ISYM(NGAS), ISYM_TOT'
C      write(6,*) JSYM,ISYM(NGAS),ISYM_TOT
*
C       IF(MINVAL(NGAS).GT.ISYM(NGAS).OR.
C    &     MAXVAL(NGAS).LT.ISYM(NGAS)    )NONEW = 1
      END IF
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE NGASDT_REL(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                      NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                      IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXSB,MXSOOB,
     &                      IBLTP,NTTSBL,LCOL,IOCOC)
*
* Number of combinations (NCOMB) with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM_REL(2,1,IASM,IBSM,ITOTSM)
C       write(6,*) ' Iasm ibsm ', iasm,ibsm
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           LTSSA = NSSOA(IASM,IATP)
           XLTSSA = LTSSA
           DO 100 IBTP = 1, MXBTP
*
#if defined LUCI_DEBUG
             IF(NTEST.GE.10) THEN
               WRITE(6,*) ' Alpha super group and beta super group'
               CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
               CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
             END IF
#endif
*
             IF(IOCOC(IATP,IBTP).EQ.1) THEN
*. Size of unpacked block
               LTTSUP =  LTSSA*NSSOB(IBSM,IBTP)
*. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = LTSSA*NSSOB(IBSM,IBTP)
                 XNCOMB = XNCOMB + XLTSSA*FLOAT(NSSOB(IBSM,IBTP))
               ELSE
                 LTTSBL = LTSSA *(LTSSA +1)/2
                 XNCOMB = XNCOMB + XLTSSA*(XLTSSA+1.0D0)*0.5
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               NTTSBL = NTTSBL + 1
               LCOL = LCOL + NSSOB(IBSM,IBTP)
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION NSXFSM2(NSMOB,NO1PS,NO2PS,ISXSM,ISYM,IPRNT)
*
* Number of single excitations of symmetry ISXSM
*
* ISYM = 0 : All symmetry allowed excitations
* ISYM = 1 : Only excitations a+iaj with I.ge.J
* ISYM =-1 : Only excitations a+iaj with I.gt.J
*
* Compared to NSXFSM : Reference to symmetries (ADSXA removed,
*                      replaced with call to SYMCOM
*
COLD  INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER NO1PS(*),NO2PS(*)
*
      MSXFSM = 0
C     WRITE(6,*) ' NO1PS in NSXFSM2 '
C     CALL IWRTMA(NO1PS,1,NSMOB,1,NSMOB)
C     WRITE(6,*) ' NSMOB ',NSMOB
      DO 100 IO1SM = 1,NSMOB
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM_REL(2,1,IO1SM,IO2SMA,ISXSM)
COLD    IO2SMA = (IO1SM,ISXSM)
*. IO2SMA is symmetry of annihilation operator, find symmetry of creation op
        IO2SM = I_GET_ADJSYM(IO2SMA)
C               I_GET_ADJSYM(ISYM)
C?      WRITE(6,*) ' IO1SM,IO2SM',IO1SM,IO2SM
        IF(ISYM.EQ.0.OR.IO1SM.GT.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*NO2PS(IO2SM)
        ELSE IF( ISYM.EQ. 1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
        ELSE IF( ISYM.EQ.-1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)-1)/2
        END IF
  100 CONTINUE
*
      NSXFSM2 = MSXFSM
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A32,I8,A1,I8)')
     &  ' Number of single excitations of symmetry ',
     &    ISXSM,',',NSXFSM2
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = allowed distributions of electrons
*
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION  IOCCLS(NGAS,*)
*. Local scratch
#include "mxpdim.inc"
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
#endif
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
#if defined LUCI_DEBUG
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : '
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
#endif
            CALL ICOPY(NGAS,IOC,1,IOCCLS(1,NOCCLS),1)
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE OCC_INCLUDED(IOCC,INCLUDED)
*
* A supergroup is given as occupation for each gas space
*
* Check if all these groups are included, and return 
* answer in INCLUDED
*
* Jeppe and Timo, Oct 2001 
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
*. Specific input
      INTEGER IOCC(NGAS)
#include "ipoist8.inc"
*
      NTEST = 000
*
      INCLUDED = 1
      DO IOBTP = 1, NGAS
        JJGRP = 0
        DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
          IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
        END DO
        IF(JJGRP.EQ.0 ) INCLUDED = 0
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation  '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Included = ',INCLUDED
      END IF
#endif
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ORBINF_REL(LUOUT,IPRNT)
*
* Obtain information about orbitals from shell information
*
* =====
* input
* =====
* Shell and symmetry information in /LUCINP/
*
* ======
* Output
* ======
* Orbital information in /ORBINP/
*
* Jeppe Olsen, Winter of 1991
*              Updated, July 97 : Double group symmetry added

*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      INTEGER CITYP
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cgas.inc"
*
#include "orbinp.inc"
#include "symm.inc"
#include "csm.inc"
C
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
*. Double group
      IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        IDBG = 1
*. number of symmetries in point group sym
        NSMSH = NSMOB/2
      ELSE
        IDBG = 0
        NSMSH = NSMOB
      END IF
*. Total number of shells per symmetry
      DO ISHSM = 1, NSMSH
        LORB = 0
        DO JGAS = 1, NGAS
          LORB = LORB + NGSSH(ISHSM,JGAS)
        END DO
        NGSSHT(ISHSM) = LORB
      END DO
*
CTF   DO ISMSH = 1, NSMSH
CTF     IF(ISMSH.EQ.1) THEN
CTF      IBSH(ISMSH) = 1
CTF     ELSE
CTF      IBSH(ISMSH) = IBSH(ISMSH-1)+NGSSHT(ISMSH)
CTF     END IF
CTF   END DO
* There is something fishy here; try this:
      IBSH(1) = 1
      do ISMSH = 2,NSMSH,1
         IVAL = 1
         do I = 1,ISMSH-1,1
            IVAL = IVAL + NGSSHT(I)
         end do
         IBSH(ISMSH) = IVAL
      end do
*
C?    WRITE(6,*) ' IBSH in ORBINP '
C?    CALL IWRTMA(IBSH,1,NSMSH,1,NSMSH)

*
* Point group => double group information if required
* already done. Info on symm.inc
      call grinfo(PNTGRP)
*                                              *
* 1 : From shell format to orbital format *
*
      CALL OSPIR_REL(NOSPIR,IOSPIR,PNTGRP,NIRREP,IPRNT,
     &               NORB_SET)
*
* 2 : Shell information to orbital information for each group of orbital
*
*
* ===============
*     GAS case
* ===============
*
      DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
        CALL SHTOOB_REL(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
     &              IMOSF_SP(1,1),NGSOB(1,IGAS),NGSOBT(IGAS) )
      END DO
*
      IF(IDBG.EQ.1) THEN
*. A second orbital set must be generated corresponding to beta
* spinors. The corresponding orbital information is stored in
* NGSOBT2 NGSOBT2
*. Orbital => spinors for beta spintors
C?      WRITE(6,*) ' Generate beta spinors '
C?      WRITE(6,*) ' NIRR_DG ', NIRR_DG
C       DO IRREP = 1, NIRR_DG
C         IOSPIR(1,IRREP) = IOSPIR(2,IRREP)
C       END DO
        DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
          CALL SHTOOB_REL(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
     &         IMOSF_SP(1,2),NGSOB2(1,IGAS),NGSOBT2(IGAS))
C    &                IOSPIR,NGSOB2(1,IGAS),NGSOBT2(IGAS))
        END DO
      END IF
*
*  ========================================================
*. Number of inactive, active, occupied , deleted orbitals
*  ========================================================
*
*
* current inactive and deleted orbitals are not identified so
      IGSINA = 0
      IGSDEL = 0
*
      CALL IZERO(NTOOBS,NSMOB)
      CALL IZERO(NOCOBS,NSMOB)
      CALL IZERO(NACOBS,NSMOB)
*
      NTOOB = 0
      NACOB = 0
      NOCOB = 0
      DO IGAS = 1, NGAS
*. Inactive orbitals
        IF(IGAS.EQ.IGSINA) THEN
          CALL ICOPVE(NGSOB(1,IGAS),NINOBS,NSMOB)
          NINOB = NGSOBT(IGAS)
        END IF
*. Deleted orbitals
        IF(IGAS.EQ.IGSDEL) THEN
          CALL ICOPVE(NGSOB(1,IGAS),NDEOBS,NSMOB)
          NDEOB = NGSOBT(IGAS)
        END IF
*. Add to total number of orbitals
        CALL IVCSUM(NTOOBS,NTOOBS,NGSOB(1,IGAS),1,1,NSMOB)
        NTOOB = NTOOB + NGSOBT(IGAS)
*. Add to occupied orbitals
        IF(IGAS.NE.IGSDEL) THEN
          CALL IVCSUM(NOCOBS,NOCOBS,NGSOB(1,IGAS),1,1,NSMOB)
          NOCOB = NOCOB + NGSOBT(IGAS)
        END IF
*. Add to active orbitals
        IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
          CALL IVCSUM(NACOBS,NACOBS,NGSOB(1,IGAS),1,1,NSMOB)
          NACOB = NACOB + NGSOBT(IGAS)
        END IF
      END DO
* ===============================================
*. Well, report back
* ===============================================
CSK   NTEST = 100
#if defined LUCI_DEBUG
      IF(NTEST.GT.0) THEN
        IF(IDBG.EQ.0) THEN
          WRITE(LUOUT,'(/A)') ' Number of orbitals per symmetry :'
          WRITE(LUOUT,'( A)') ' ================================='
        ELSE
          WRITE(LUOUT,'(/A)') ' Number of up spinors per symmetry :'
          WRITE(LUOUT,'( A)') ' =================================='
        END IF
        WRITE(LUOUT,'(/A,10I4,A)')
     &  '             Symmetry  ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(A,2X,10A)')
     &  '             ========== ',('====',I = 1,NSMOB)
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(         A,I3,7X,A,10I4,8X,I3)')
     &    '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB)
        END DO
*
        WRITE(LUOUT,*) ' Total number of orbitals ', NTOOB
        WRITE(LUOUT,*) ' Total number of occupied orbitals ', NOCOB
      END IF
#endif
*. Offsets for orbitals of given symmetry
      ITOOBS(1) = 1
      DO  ISMOB = 2, NSMOB
        ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Offsets for orbital of given symmetry '
        CALL IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
      END IF
#endif
*
      IF(IDBG.EQ.1) THEN
*.     Repeat the above for construction of beta spinors
*      current inactive and deleted orbitals are not identified so
       IGSINA = 0
       IGSDEL = 0
*
       CALL IZERO(NTOOBS2,NSMOB)
       CALL IZERO(NOCOBS2,NSMOB)
       CALL IZERO(NACOBS2,NSMOB)
*
       NTOOB2 = 0
       NACOB2 = 0
       NOCOB2 = 0
       DO IGAS = 1, NGAS
*. Inactive orbitals
         IF(IGAS.EQ.IGSINA) THEN
           CALL ICOPVE(NGSOB2(1,IGAS),NINOBS2,NSMOB)
           NINOB2 = NGSOBT2(IGAS)
         END IF
*. Deleted orbitals
         IF(IGAS.EQ.IGSDEL) THEN
           CALL ICOPVE(NGSOB2(1,IGAS),NDEOBS2,NSMOB)
           NDEOB2 = NGSOBT(IGAS)
         END IF
*. Add to total number of orbitals
         CALL IVCSUM(NTOOBS2,NTOOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
         NTOOB2 = NTOOB2 + NGSOBT2(IGAS)
*. Add to occupied orbitals
         IF(IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NOCOBS2,NOCOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
           NOCOB2 = NOCOB2 + NGSOBT2(IGAS)
         END IF
*. Add to active orbitals
         IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NACOBS2,NACOBS2,NGSOB2(1,IGAS),1,1,NSMOB)
           NACOB2 = NACOB2 + NGSOBT2(IGAS)
         END IF
       END DO
* ===============================================
*. Well, report back
* ===============================================
#if defined LUCI_DEBUG
       IF(NTEST.GT.0) THEN
         WRITE(LUOUT,*)
         WRITE(LUOUT,*) ' Number of down spinors per symmetry :'
         WRITE(LUOUT,*) ' ===================================='
         WRITE(LUOUT,*)
         WRITE(LUOUT,'(1H ,A,10I4,A)')
     &   '            Symmetry  ',(I,I = 1,NSMOB)
         WRITE(LUOUT,'(1H ,A,2X,10A)')
     &   '           ========== ',('====',I = 1,NSMOB)
         DO IGAS = 1, NGAS
           WRITE(LUOUT,'(1H      ,A,I3,7X,A,10I4,8X,I3)')
     &     '   GAS',IGAS,'      ',(NGSOB2(I,IGAS),I=1,NSMOB)
         END DO
         WRITE(LUOUT,*) ' Total number of orbitals ', NTOOB2
         WRITE(LUOUT,*) ' Total number of occupied orbitals ', NOCOB2
       END IF
#endif
*. Offsets for orbitals of given symmetry
       ITOOBS2(1) = 1
       DO  ISMOB = 2, NSMOB
         ITOOBS2(ISMOB) = ITOOBS2(ISMOB-1)+NTOOBS2(ISMOB-1)
       END DO
*
#if defined LUCI_DEBUG
       IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Offsets for orbital of given symmetry '
         CALL IWRTMA(ITOOBS2,1,NSMOB,1,NSMOB)
       END IF
#endif
      END IF
*.Information about spinor symmetries
#if defined LUCI_DEBUG
      if (NTEST.ge.1) then
        WRITE(LUOUT,'(/A)')
     &  '  Spin-up Spinors in Double group symmetry  '
        WRITE(LUOUT,'(A)')
     &  ' =========================================== '
        WRITE(LUOUT,'(/A,10I4,A)')
     &    '                 Irrep ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(A,2X,10A,A)')
     &    '                 ===== ',('====',I = 1,NSMOB )
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &    '        GAS',IGAS,'          ',
     &    (NGSOB(IRREP,IGAS),IRREP = 1, NSMOB )
        END DO
*
        WRITE(LUOUT,'(/A)')
     &  '  Spin-down Spinors in Double group symmetry  '
        WRITE(LUOUT,'(A)')
     &  ' =========================================== '
        WRITE(LUOUT,'(/A,10I4,A)')
     &    '                 Irrep ',(I,I = 1,NSMOB)
        WRITE(LUOUT,'(A,2X,10A,A)')
     &    '                 ===== ',('====',I = 1,NSMOB )
        DO IGAS = 1, NGAS
          WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &    '        GAS',IGAS,'          ',
     &    (NGSOB2(IRREP,IGAS),IRREP = 1, NSMOB )
        END DO
      end if
#endif
*
********************************************
*                                          *
* Part 2 : Reordering arrays for orbitals  *
*                                          *
********************************************
      IPRNT = NTEST
        CALL ORBORD_GAS_REL(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &       NOCOBS,NTOOBS,NTOOB,
     &       IREOST,IREOTS,ISMFTO,ITPFSO,
     &       IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &       ISMFSO,ITPFTO,NOBPT,IPRNT,IOSPIR,1)
      IF(IDBG.EQ.1) THEN
*. Reorder array for down spinors
        CALL ORBORD_GAS_REL(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB2,NGSOBT2,
     &       NOCOBS2,NTOOBS2,NTOOB2,
     &       IREOST2,IREOTS2,ISMFTO2,ITPFSO2,
     &       IBSO2,NTSOB2,IBTSOB2,ITSOB2,NOBPTS2,IOBPTS2,
     &       ISMFSO2,ITPFTO2,NOBPT,IPRNT,IOSPIR,2)
        END IF
        NTEST = 0
        IPRNT = 0
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ORBORD_GAS_REL(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &                  NOCOBS,NTOOBS,NTOOB,
     &                  IREOST,IREOTS,ISFTO,ITFSO,IBSO,
     &                  NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &                  ISFSO,ITFTO,NOBPT,IPRNT,IORD,ICOMP)
*
* Obtain Reordering arrays for orbitals
* ( See note below for assumed ordering )
*
*
* GAS version
*
* =====
* Input
* =====
*  NSMOB  : Number of orbital symmetries
*  MXPOBS : Max number of orbital symmetries allowed by program
*  MXPNGAS: Max number of GAS spaces allowed by program
*  NGAS   : Number of GAS spaces
*  NGSOB  : Number of GAS orbitals per symmetry and space
*  NGSOBT : Number of GAS orbitals per space
*  NOCOBS : Number of occupied orbitals per symmetry
*  NTOOBS : Number of orbitals per symmetry,all types
*  IORD   : Order in which orbital symmetries are filled
*
* ======
* Output
* ======
*  IREOST : Reordering array symmetry => type
*  IREOTS : Reordering array type     => symmetry
*  ISFTO  : Symmetry array for type ordered orbitals
*  ITFSO  : Type array for symmetry ordered orbitals( not activated )
*  IBSO   : First orbital of given symmetry ( symmetry ordered )
*  NOBPTS : Number of orbitals per subtype and symmetry
*  IOBPTS : Offsets for orbitals of given subtype and symmetry
*           ordered according to input integrals
*
* ISFSO  : Symmetry of orbitals, symmetry ordereing
* ITFTO  : Type of orbital, type ordering
*
* Jeppe Olsen, Winter 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION NGSOB(MXPOBS,MXPNGAS),NOCOBS(*),NTOOBS(*)
      DIMENSION NGSOBT(MXPNGAS)
      DIMENSION IORD(MXPOBS,*)
*. Output
      DIMENSION IREOST(*),IREOTS(*),ISFTO(*),ITFSO(*),IBSO(*)
      DIMENSION ISFSO(*),ITFTO(*)
      DIMENSION NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
      DIMENSION NOBPT(MXPNGAS)

* ==========================
* Note on order of orbitals
* ==========================
*
* The orbitals are supposed to be imported ordered symmetry-type
* ordered as
*
* Loop over symmetries of orbitals
*  Loop over GAS spaces
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of Loop over Gas spaces
* End of loop over symmetries
*
* Internally the orbitals are reordered to type symmetry order
* where the outer loop is over types and the inner loop is
* over symmetries, i.e.
*
* Loop over GAS spaces
*  Loop over symmetries of orbitals
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of loop over symmetries
* End of Loop over Gas spaces
*
*. 1:  Construct ISFTO, ITFTO, IREOST,IREOTS,NOBPTS,IOBPTS
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
*
#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
        WRITE(6,*) ' ICOMP = ', ICOMP
        WRITE(6,*) ' Four elements of relevant IORD row '
        WRITE(6,*)
     &   IORD(ICOMP,1),IORD(ICOMP,2),
     &   IORD(ICOMP,3),IORD(ICOMP,4)
      end if
#endif
*
      ITSOFF = 1
      DO IGAS = 1, NGAS
        DO IISYM = 1, NSMOB
          ISYM = IORD(ICOMP,IISYM)
          IF(IISYM.EQ.1) THEN
            IBSSM = 1
          ELSE
            IBSSM = IBSSM + NTOOBS(IORD(ICOMP,IISYM-1))
          END IF
          NPREV = 0
          DO JGAS = 1, IGAS-1
            NPREV = NPREV + NGSOB(ISYM,JGAS)
          END DO
          IADD = 0
          NOBPTS(IGAS,ISYM) = NGSOB(ISYM,IGAS)
          IOBPTS(IGAS,ISYM) = ITSOFF
          DO IORB = ITSOFF,ITSOFF+NGSOB(ISYM,IGAS)-1
            IADD = IADD + 1
            IREOTS(IORB) = IBSSM-1+NPREV+IADD
            IREOST(IBSSM-1+NPREV+IADD) = IORB
            ITFTO(IORB) = IGAS
            ISFTO(IORB) = ISYM
          END DO
          ITSOFF = ITSOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 2 : ISFSO,ITFSO
*
      ISTOFF = 1
      DO IISYM = 1, NSMOB
        ISYM = IORD(ICOMP,IISYM)
        DO IGAS = 1, NGAS
          DO IORB = ISTOFF,ISTOFF+NGSOB(ISYM,IGAS)-1
            ISFSO(IORB) = ISYM
            ITFSO(IORB) = IGAS
          END DO
          ISTOFF = ISTOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 3 IBSO, NOBPT
*
      IOFF = 1
      DO IISM = 1, NSMOB
       ISM = IORD(ICOMP,IISM)
       IBSO(ISM) = IOFF
       IOFF = IOFF + NTOOBS(ISM)
      END DO
      DO IGAS = 1, NGAS
        NOBPT(IGAS) = NGSOBT(IGAS)
      END DO
*
#if defined LUCI_DEBUG
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Output from ORBORD '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry of orbitals , type ordered '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Symmetry => type reordering array '
        CALL IWRTMA(IREOST,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Type => symmetry reordering array '
        CALL IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' IBSO array '
        CALL IWRTMA(IBSO,1,NSMOB,1,NSMOB)
*
        WRITE(6,*) ' NOBPTS '
        CALL IWRTMA(NOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
        WRITE(6,*) ' NOBPT '
        CALL IWRTMA(NOBPT,NGAS,1,MXPNGAS,1)
        WRITE(6,*) ' IOBPTS '
        CALL IWRTMA(IOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
*
        WRITE(6,*) ' ISFTO array : '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFSO array : '
        CALL IWRTMA(ITFSO,1,NTOOB,1,NTOOB)
*
        WRITE(6,*) ' ISFSO array : '
        CALL IWRTMA(ISFSO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFTO array : '
        CALL IWRTMA(ITFTO,1,NTOOB,1,NTOOB)
      END IF
#endif
*

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE OSPIR_REL(NOSPIR,IOSPIR,PNTGRP,NIRREP,IPRNT,
     &                     NORB_SET)
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   PNTGRP  : type of pointgroup
*         = 1 => D2h or a subgroup of D2H
*         = 2 => C inf v
*         = 3 => D inf h
*         = 4 => O 3
*         = 5 to 9 => different double groups :
*             9 :   C2H
*             8 :   C2
*             7 :   CS
*             6 :   CI
*             5 :   C1
*   NIRREP : Number of irreducible representations per point group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPOBS : Largest allowed number of orbital symmetries
*
* ======
* Output
* ======
*
*   NOSPIR : Number of orbital symmetries per irrep
*   IOSPIR : Orbital symmetries corresponding to a given irrep
*
* Jeppe Olsen , Winter of 1991
*              Updated, July 97 : Double group symmetry added
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
      INTEGER PNTGRP
*. Double group info
#include "symm.inc"
*. Output
      DIMENSION NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR)
*
*. In general -except for double groups we have just a single
*  orbitals set
      NORB_SET = 1
      IF(PNTGRP.GE.5) NORB_SET = 2
      if (PNTGRP.eq.1) then
*=====
*.D2h
*=====
        NSMOB = 0
        DO 10 IRREP = 1, 8
          NOSPIR(IRREP) = 1
          IOSPIR(1,IRREP) = IRREP
   10   CONTINUE
      ELSE IF(PNTGRP.EQ.2) THEN
* =========
*. C inf V
* =========
* orbital symmetry is numbered as IML - MNMLOB + 1
        MNMLOB = -(NIRREP-1)
        DO 20 IRREP = 1, NIRREP
*.Irrep I contains orbitals with ML = -(IRREP-1),+(IRREP-1)
          IF(IRREP.EQ.1) THEN
            NOSPIR(IRREP) = 1
            IOSPIR(1,IRREP) = IRREP - 1 - MNMLOB + 1
          ELSE
            NOSPIR(IRREP) = 2
            IOSPIR(1,IRREP) = -(IRREP - 1) - MNMLOB + 1
            IOSPIR(2,IRREP) =  (IRREP - 1) - MNMLOB + 1
          END IF
   20   CONTINUE
      ELSE IF(PNTGRP.EQ.3) THEN
* ========
*. D inf H
* ========
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        IRREP = 0
        DO 35 IPARI = 1, 2
          IADD = (IPARI-1)*NMLOB
          DO 30 ML = 0,MXMLOB
            IRREP = IRREP + 1
            IF(ML.EQ.0) THEN
              NOSPIR(IRREP) = 1
              IOSPIR(1,IRREP) = IADD + ML - MNMLOB + 1
            ELSE
              NOSPIR(IRREP) = 2
              IOSPIR(1,IRREP) = IADD - ML - MNMLOB + 1
              IOSPIR(2,IRREP) = IADD + ML - MNMLOB + 1
            END IF
   30     CONTINUE
   35   CONTINUE

      ELSE IF(PNTGRP.EQ.4) THEN
* =====
*. O 3
* =====
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        DO 45 L = 0, NIRREP - 1
          IF(MOD(L,2).EQ.0) THEN
            IPARI = 1
          ELSE
            IPARI = 2
          END IF
          IADD = (IPARI-1)*NMLOB
          IRREP = L + 1
          NOSPIR(IRREP) = 2 * L + 1
          ICOMP = 0
          DO 40 ML = MNMLOB,MXMLOB
            ICOMP = ICOMP + 1
            IOSPIR(ICOMP,IRREP) = IADD + ML - MNMLOB + 1
   40     CONTINUE
   45   CONTINUE
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ========================
*. Double group symmetry
* ========================
*
* We are handling the transformation from orbitals to spinors
* by having two sets of orbitals. This turns out to be
* significant more convenient ( we hope ) than having
* a single set with double dimension
        NORB_SET = 2
        DO  IRREP = 1, NIRR_DG
          NOSPIR(IRREP) = 1
        END DO
*. IOSPIR become the symmetry of the spinors
        DO  IRREP = 1, NIRR_DG
*. The first set corresponds to alpha spin functions
          IOSPIR(1,IRREP) = IMOSF_SP(IRREP,1)
*. The second set corresponds to beta spin functions
          IOSPIR(2,IRREP) = IMOSF_SP(IRREP,2)
        END DO
*. The first NSMOB/2 symmetries are now the fermion irreps, whereas
*  the boson irreps are missing. Add these for completeness
        DO ISMOB = 1, NIRR_PN
         IOSPIR(1,NIRR_PN+ISMOB) = ISMOB
         IOSPIR(2,NIRR_PN+ISMOB) = ISMOB
        END DO
      ELSE
        WRITE(6,*) ' Sorry  PNTGRP out of range , PNTGRP = ', PNTGRP
        WRITE(6,*) ' OSPIR fatally wounded '
        Call Abend1( 5 )
      END IF
*
CSK   NTEST = 10
      NTEST = MAX(IPRNT,NTEST)
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' OSPIR speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Number of orbital sets ',NORB_SET
        WRITE(6,*) ' Number of orbitals per irrep '
*
        IF(.NOT.(PNTGRP.GE.5.AND.PNTGRP.LE.9)) THEN
*. Info for orbitals
          CALL IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
          WRITE(6,*) ' Orbital symmetries per irrep '
          DO 100 IRREP = 1, NIRREP
          CALL IWRTMA(IOSPIR(1,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
  100     CONTINUE
        ELSE
*. Info for spinors
          WRITE(6,*) ' IOSPIR for alpha spinors '
          WRITE(6,'(10I4)') (IOSPIR(1,I),I=1,NIRR_DG)
          WRITE(6,*) ' IOSPIR for beta spinors '
          WRITE(6,'(10I4)') (IOSPIR(2,I),I=1,NIRR_DG)
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV3     '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(6,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      IASM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, IASM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,IASM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(IASM).EQ.0) GOTO 1000
      IF(IBLTP(IASM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(IASM)
*. can this block be included
      IBSM = ISMOST(IASM)
      NSTA = NSSOA(IASM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(IASM).EQ.1.OR.(IBLTP(IASM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(IASM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?    write(6,*) ' IA IB IASM LBLOCK ', IA,IB,IASM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = IASM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV3'
        WRITE(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        stop ' In PART_CIV3. '
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 'Output from PART_CIV3'
        WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV4(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
*
* IMPORTANT change compared to PART_CIV3:
*
* MXLNG is no longer assumed to be equals to the total number 
* of dets. Therefore we have to partition the CI vector being aware  
* of putting symmetry blocks of same TT into different batches.
*
*  
* REASON: SIGMA-vector / DENSITY matrix calculation is done wrt loops
* over NSMST! 
* 
* S. Knecht - July 02 2007
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      I_SAVE_IBLOCK = 0
      N_TRY = 0
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' =================='
        WRITE(LUWRT,*) '     PART_CIV4     '
        WRITE(LUWRT,*) ' =================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMAMN(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA,LUWRT)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMAMN(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB,LUWRT)
      END IF
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV4'
        WRITE(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        stop ' In PART_CIV4. '
      ELSE
*
*     this batch is finished, goto next batch
*
*       is this the TT-1 block?
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include'//
     &                 ' two TT-blocks of different symmetry in'//
     &                   'a batch'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV4'
            WRITE(LUWRT,*) ' block LENGTH, actual batch length ',
     &                       LENGTH, LEBATCH(NBATCH)
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            stop ' In PART_CIV4 because of N_TRY .ge. 1. '
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,6)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV4'
        WRITE(LUWRT,*) '====================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '  Number of blocks included ', LBATCH(JBATCH)
          WRITE(LUWRT,*) '  TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV3_PAR(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) '     PART_CIV3_PAR    '
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
C       only s-blocks corresponding to distribution
C*causes error?!?        IF( IBLOCKD(IBLOCK) .ne. MYPROC ) IBATCH(1,IBLOCK) = 0
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV3'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH, LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        CALL QUIT( ' In PART_CIV3_PAR. ')
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV3_PAR'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lenghts of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV3_SPC(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ========================='
        WRITE(LUWRT,*) '     PART_CIV3_PAR_SPC    '
        WRITE(LUWRT,*) ' ========================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
      I_SAVE_IBLOCK = 0
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
      N_TRY = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV3_SPC'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        stop ' In PART_CIV3_SPC. '
      ELSE
*. This batch is finished, goto next batch
*
*       is this the TT-1 block?
*
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include a '//
     &                     ' single Block'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV3_SPC'
            WRITE(LUWRT,*) '  LENGTH, LBLOCK ',LENGTH, LBLOCK
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            CALL RMEMCHK
            CALL QUIT( ' In PART_CIV3_SPC because of N_TRY .ge. 1. ')
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,LUWRT)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV3_SPC'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lenghts of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE PNT2DM2(I12SM,NSMOB,IPSM,JPSM,IJSM,IPNTR)
*
* Pointer to two dimensional array
*
* Jeppe Olsen, Jan. 98
*
* Modified from PNT2DM2 : Reference to symmetry products referred to
*                         SYMCOM
*                         Generation of symmetry of last index eliminated
*
* =====
* Input
* =====
* I12SM  : ne.0 => restrict to lower half
*          eq.0 => completer matrix
* NSMOB : Number of orbital symmetries
* IPSM : Number of orbitals per symmetry for index 1
* JPSM : Number of orbitals per symmetry for index 2
* IJSM  : Symmetry of two index array
*
* =======
* Output
* =======
* IPNTR : Pointer to block with first index of given symmetry
*         = 0 indicates forbidden block
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER IPSM(*),JPSM(*)
*.Output
      DIMENSION IPNTR(*)
*
      CALL IZERO(IPNTR,NSMOB)
COLD  CALL ISETVC(ISM2 ,IZERO,NSMOB)
*
      IOFF = 1
      DO 100 ISM = 1,NSMOB
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM_REL(2,1,ISM,JSMA,IJSM)
        IF(JSMA.EQ.0) GOTO 100
        JSM = I_GET_ADJSYM(JSMA)
        IF(I12SM.EQ.0.OR.ISM.GE.JSM) THEN
*. Allowed block
          IPNTR(ISM) = IOFF
          IF(I12SM.GT.0.AND.ISM.EQ.JSM) THEN
            IOFF = IOFF + IPSM(ISM)*(IPSM(ISM)+1)/2
          ELSE
            IOFF = IOFF + IPSM(ISM)*JPSM(JSM)
          END IF
        END IF
  100 CONTINUE
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' dimension of two-dimensional array ',IOFF-1
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Pointer '
        CALL IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE PNT4DM2(NSMOB,NO1PS,NO2PS,NO3PS,NO4PS,
     &                   IDXSM,IS12,IS34,IS1234,IPNTR)
*
* Pointer for 4 dimensional array with total symmetry IDXSM
* Pointer is given as 3 dimensional array corresponding
* to the first 3 indeces
*
* IS12 (0,1,-1)   : Permutational symmetry between indeces 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between indeces 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between indeces 12 and 34
*
* Jeppe Olsen, Feb. 98
*
* Derived from PNT4DM : Symmetry handling transferred to SYMCOM
*
* note ISM4A is not active
*
*.input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
#include "ipoist8.inc"
*.Output
      INTEGER IPNTR(NSMOB,NSMOB,NSMOB)
*
      NTEST = 000
      CALL IZERO(IPNTR,NSMOB ** 3 )
*
      IOFF= 1
      DO 10 I1SM = 1, NSMOB
        DO 20 I2SM= 1, NSMOB
*. I2SM is symmetry of creation operator, symmetry of annihilation op
          I2SMA = I_GET_ADJSYM(I2SM)
          CALL SYMCOM_REL(3,1,I1SM,I2SMA,I12SM)
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
C         I12SM = ADASX(I1SM,I2SM)
          CALL SYMCOM_REL(2,1,I12SM,I34SM,IDXSM)
C         I34SM = SXDXSX(I12SM,IDXSM)
          IF(I34SM.EQ.0) GOTO 20
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 20
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMOB+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
          DO 30 I3SM = 1, NSMOB
            CALL SYMCOM_REL(2,1,I3SM,I4SMA,I34SM)
            I4SM = I_GET_ADJSYM(I4SMA)
C           I4SM = ADSXA(I3SM,I34SM)
            IF(I4SM.EQ.0) GOTO 30
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 30
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMOB+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 30
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF+ N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF + N12*(N12+1)/2
            ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
C             ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF=  IOFF+ N12*(N12-1)/2
            END IF
#if defined LUCI_DEBUG
            if (NTEST.ge.100) then
              WRITE(6,*) ' I1SM I2SM I3SM I4SM    IOFF'
              WRITE(6,'(1H ,4I4,I9)')   I1SM,I2SM,I3SM,I4SM,IOFF
            end if
#endif
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*
*
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Length of 4 dimensional array ', IOFF - 1
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE RSMXMN_REL(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                      MIN1_T,MAX1,MIN3,MAX3_T,NTEST)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION MINEL(*),MAXEL(*)
*
      NORB = NORB1 + NORB2 + NORB3
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
COLD  MIN1A = MIN1_T
      MIN1A = MAX(MIN1_T,NEL-MAX3_T-NORB2)
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3_T
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF (NORB1.GT.0) MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF (NORB1+NORB2.GT.0)
     &        MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
*
      IF( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from RSMXMN '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SHTOOB_REL(NSHPIR,NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,
     &                      NOBPS,NOB)
*
* Number of shells per irrep => Number of orbitals per symmetry
*
* =====
* Input
* =====
*
*  NSHPIR : Number of shells per irrep
*  NIRREP : Number of irreps
*  MXPOBS : Largest allowed number of orbitals symmetries
*  NSMOB  : Number of orbital symmetries
*  NOSPIR : Number of orbital symmetries per irrep
*  IOSPIR : Orbital symmetries per irrep
*
* ======
* Output
* ======
*  NOBPS  : Number of orbitals per symmetry
*  NOB    : Number of orbitals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
C
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
C
*. Input
      DIMENSION NSHPIR(*),NOSPIR(*),IOSPIR(*)
*. Output
      DIMENSION NOBPS(*)
*
CSK      IF( MYPROC .le. MASTER + 1 ) THEN
CSK        WRITE(6,*) 'THESE ARE MY IOSPIR and NSHPIR matrix',myproc
CSK        CALL IWRTMA(IOSPIR,1,NIRREP,1,NIRREP)
CSK        CALL IWRTMA(NSHPIR,1,NIRREP,1,NIRREP)
CSK      END IF
*
      CALL IZERO(NOBPS,NSMOB)
      NOB = 0
      JSM = 0
      DO 100 IRREP = 1, NIRREP
C         IISM = IOSPIR(ISM,IRREP)
          IISM  = IOSPIR(IRREP)
          if (IISM.ne.0) then
             JSM = JSM + 1
             NOBPS(IISM) = NOBPS(IISM) + NSHPIR(IRREP)
          end if
          NOB = NOB + NSHPIR(IRREP)
  100 CONTINUE
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' SHTOOB Speaking '
         WRITE(6,*) ' =============== '
         WRITE(6,*) ' Number of orbitals obtained ', NOB
         WRITE(6,*) ' Number of orbitals per symmetry '
         CALL IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SIGVST_REL(ISGVST,NSMST)
*
* Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      INTEGER ISGVST(*)
*
      DO 100 ISM = 1, NSMST
C            MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISM,'ST',2)
        MIML = - IML
        CALL MLSM_REL(MIML,IPARI,MISM,'ST',1)
        ISGVST(ISM) = MISM
  100 CONTINUE
*
      NTEST = 1
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' ISGVST array '
        WRITE(6,*) ' ============ '
        CALL IWRTMA(ISGVST,1,NSMST,1,NSMST)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SMOST_REL(NSMST,NSMCI,MXPCSM,ISMOST)
*
* ISMOST(ISYM,ITOTSM) : Symmetry of an internal state is ITOTSM
*                       if symmetry of 1 string is ISYM, the
*                       symmetry of the other string is
*                       ISMOST(ISYM,ITOTSM)
*
* Jeppe Olsen , Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION ISMOST(MXPCSM,MXPCSM)
*
      DO 1000 ITOTSM = 1, NSMCI
       DO 900 ISTSM  = 1, NSMST
C            SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM_REL(2,1,ISTSM,JSTSM,ITOTSM)
        ISMOST(ISTSM,ITOTSM) = JSTSM
  900  CONTINUE
 1000 CONTINUE
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF( NTEST.NE. 0 ) THEN
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Info from SMOST '
        WRITE(6,*) ' ==============='
        DO 1010 ITOTSM = 1, NSMCI
          WRITE(6,*) ' ISMOST array for ITOTSM = ', ITOTSM
          CALL IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
 1010   CONTINUE
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE STSTSM_REL(STSTSX,STSTDX,NSMST)
*
* construct  STSTSX and STSTDX giving
* symmetry of sx (dx) connecting two given string symmetries
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      INTEGER STSTSX(NSMST,NSMST),STSTDX(NSMST,NSMST)
*
      DO 100 ILSTSM = 1, NSMST
        DO 50 IRSTSM = 1, NSMST
          CALL SYMCOM_REL(1,5,ISXSM,IRSTSM,ILSTSM)
          CALL SYMCOM_REL(1,6,IDXSM,IRSTSM,ILSTSM)
          STSTSX(ILSTSM,IRSTSM) = ISXSM
          STSTDX(ILSTSM,IRSTSM) = IDXSM
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' STSTSM : STSTSX, STSTDX '
        CALL IWRTMA(STSTSX,NSMST,NSMST,NSMST,NSMST)
        CALL IWRTMA(STSTDX,NSMST,NSMST,NSMST,NSMST)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SXTYP_GAS_REL(NSXTYP,ITP,JTP,NGAS,MXPNGAS,ILTP,IRTP)
*
* Two supergroups are given. Find single excitations that connect
* these supergroups
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS)
*. Output
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one single excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SYMCMDG(ITASK,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1*I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* Double group version, written for compatibility with general symmetry
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "symm.inc"
*
      IF(ITASK.EQ.1) THEN
        I1  = IDBGMULT(INVELM(I2),I12)
      ELSE IF(ITASK.EQ.2) THEN
        I2  = IDBGMULT(INVELM(I1),I12)
      ELSE IF (ITASK.EQ.3) THEN
        I12 = IDBGMULT(I1,I2)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = I12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I2
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* IOBJ = 1 : I1,I2 are strings I12 determinant
* ( Other things can follow )
* IOBJ = 2 : I1,I2,I3 are externals
* IOBJ = 3 : I1 is an external, I2,I3 are dets
* IOBJ = 4 : I1 is orbital, I2 is string,l, I12 is string
* IOBJ = 5 : I1 is single excitation, I2 is string,l, I12 is string
* IOBJ = 6 : I1 is orbital, I2 is Orbital I12 is single excitation
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
* ================
*. Driver routine
* ================
*.LUCINP (PNTGRP is used )
C     INTEGER PNTGRP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "implicit.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
*
*
#if defined LUCI_DEBUG
      NTEST = 0000
      if (NTEST.ge.2000) then
        write(6,*)
        write(6,*) 'symcom test:'
        write(6,*) 'ITASK : ',ITASK
        write(6,*) 'IOBJ  : ',IOBJ
        write(6,*) 'I1    : ',I1
        write(6,*) 'I2    : ',I2
        write(6,*) 'I12   : ',I12
      end if
#endif
*
      IF(PNTGRP.EQ.1) THEN
        CALL SYMCM1(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
        CALL SYMCM2(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
        CALL SYMCMDG(ITASK,I1,I2,I12)
      ELSE
        WRITE(6,*) ' PNTGRP parameter out of bounds ', PNTGRP
        WRITE(6,*) ' Enforced stop in SYMCOM '
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE SYMINF_REL(IPRNT)
*
* Information about number of symmetries
*
* Input : /LUCINP/,/ORBINP
* Output : /CSM/,/CSMPRO/
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     INTEGER PNTGRP,EXTSPC
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*

*. Output
* NSMSX : number of symmetries of single excitations
* NSMDX : Number of symmetries of double excitations
* NSMST : Number of symmetries of strings
* NSMCI : NUmber of symmetries of CI spaces
* ITSSX : Total symmetrix single excitation
* ITSDX : Total symmetrix double excitation
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
* ADASX : symmetry of orbs i and i => symmetry of a+iaj
* ASXAD : symmetry of orb j and excit a+iaj => symmetry of i
* ADSXA : symmetry of orb i and excit a+iaj => symmetry of j
*
* SXSXDX : Symmetry of two single excitations
*          => symmetry of double  excitation
* SXDXSX : Symmetry of single excitation and double excitation
*          => symmetry of single  excitation

*.
      IF(PNTGRP.EQ.1) THEN
* =====
* D 2 h
* =====
        CALL ZSYM1_REL(NIRREP,IPRNT)
      ELSE IF(PNTGRP.EQ.2) THEN
* ========
* C inf V
* ========
C            ZNONAB_REL(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB_REL(0,MAXML,NSMOB,IPRNT)
        CALL ZSYM2_REL(IPRNT)
      ELSE IF(PNTGRP.EQ.3.OR.PNTGRP.EQ.4) THEN
* ===========
* D inf H O3
* ===========
C            ZNONAB_REL(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB_REL(1,MAXML,NSMOB,IPRNT)
        CALL ZSYM2_REL(IPRNT)
      ELSE IF(PNTGRP.GE.5.AND.PNTGRP.LE.9) THEN
* ===========
* Double group
* ============
*
        CALL ZSYMDG(IPRNT)
      ELSE
        WRITE(6,*) ' You are too early , sorry '
        WRITE(6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE Z_BLKFO_REL(IDC,NMS2VAL,ISPC,ISM,KPCLBT,KPCLEBT,
     &                       KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                       NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,
     &                       NBLK_MS2_C,CALC_MS2_C,I_SET_L2BLOCK,
     &                       I_USE_PC,NPARBLOCK)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : Block type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#include "krciprop.h"
*
      LOGICAL CALC_MS2_C
*.output
      INTEGER NBLK_MS2(*), IBLK_MS2(*), NBLK_MS2_C(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
*. Pointers to output arrays
      NTTS = MXNTTS
      call memmar(KPCLBT ,MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,MXNTTS,'ADDL  ',1,'CI1BT ')
      call memmar(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      call memmar(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*.    ^ These should be preserved after exit so put mark for 
*       flushing after
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLKF')
*
      IDUM = 0
*
*     allowed length of each batch
CSK   IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .ge. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
csk        WRITE(LUWRT,'(1X,A,1X,I20)')
csk     & '  Current available free memory in double words:',L2BLOCK
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*
csk        WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk     &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1000 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
      LTEST_BLOCK = L_COMBI_MAX
      IF( I_USE_PC .ne. 0)THEN
CSK     WRITE(LUWRT,*) ' LTEST_BLOCK, L2BLOCK, L_COMBI',
CSK  &                   LTEST_BLOCK, L2BLOCK, L_COMBI
*
        IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*         reset to L_COMBI because that is already enough
          WRITE(LUWRT,*) ' *** restricting memory usage ***'
          WRITE(LUWRT,*) ' reset L2BLOCK =',L2BLOCK,' to LTEST_BLOCK =',
     &                     LTEST_BLOCK
          WRITE(LUWRT,*) ' *** restricting memory usage ***'
          L2BLOCK = LTEST_BLOCK
*
        END IF
        LTEST_BLOCK = MIN(L2BLOCK,L_COMBI_MAX)
      END IF
*
*     ... set LBLOCK value
      LBLOCK  = MAX(L0BLOCK,L2BLOCK,LTEST_BLOCK)
      IF( I_SET_L2BLOCK .gt. 1 .or. I_USE_PC .eq. 0 ) THEN
       LBLOCK = MAX(LBLOCK,LCSBLK)
      END IF
      IF( I_SET_L2BLOCK .ge. 1 )
     &  WRITE(LUWRT,*) '  LBLOCK set to:',LBLOCK
#if defined LUCI_DEBUG
      WRITE(LUWRT,*) ' inside Z_BLKFO_REL: reset LBLOCK to',LBLOCK
      WRITE(LUWRT,*) ' LBLOCK, L0BLOCK, L2BLOCK, LCSBLK, LTEST_BLOCK',
     &                 LBLOCK, L0BLOCK, L2BLOCK, LCSBLK, LTEST_BLOCK
#endif
C
C     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
C
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
C
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
C       info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KPCBLTP))
C
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
C       number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
C
C       batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
          IDEBUGPRNT = 000
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
          IDEBUGPRNT = 000
*
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
       WRITE(LUWRT,*) ' Total number of batches', NBATCH
       WRITE(LUWRT,*) ' Total number of blocks ', NBLOCK
*
       WRITE(LUWRT,*) ' Number of blocks and batches per MS2 value '
       WRITE(LUWRT,*) ' ==========================================='
       WRITE(LUWRT,*)
       WRITE(LUWRT,*) ' IMS2    MS2      NBLOCK      NBATCH '
       WRITE(LUWRT,*) ' ===================================='
       DO IMS2 = 1, NMS2VAL
         IATP = IST_FOR_DT(1,IMS2)
         IBTP = IST_FOR_DT(2,IMS2)
         IIMS2 = NELEC(IATP)-NELEC(IBTP)
         WRITE(LUWRT,'(4(3X,I4))') 
     &              IMS2,IIMS2,NBLK_MS2(IMS2),NBAT_MS2(IMS2)
       END DO
*
      END IF
#endif
*
* Transfer information to common block for KRMC - LUCIAREL
* communication:
      if (INTIMP.eq.6.or.INTIMP.eq.7) then
        IF( MYPROC .eq. MASTER ) THEN
          call civcinf_mc(WORK(KPCIBT),WORK(KPCLBT),WORK(KPCI1BT),
     &                    WORK(KPCLEBT),NBAT_MS2,NMS2VAL,NBLOCK)
        END IF
      end if
C
C     this part is only relevant for sigma calculations as part of 
C     the property module
C
      IF( .NOT. CALC_MS2_C ) GOTO 9999
C
C     if we run a property calculation and the call to this 
C     routine occurs from inside SIGDEN_REL we need to define 
C     atleast NBLK_MS2_C which may differ from NBLK_MS2 above 
C     since either ISSM == ICSM (fine) or ISSM != ICSM (this loop
C     is needed)
C
C     hidden information: IABS(I_SET_L2BLOCK) = ICSM
C
      ICSM_SCR = IABS(I_SET_L2BLOCK)
C     WRITE(LUWRT,*) ' ISM, ICSM_SCR is',ISM, ICSM_SCR
C
      CALL MEMMAR(KPCLBTX ,  MXNTTS,'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEBTX,  MXNTTS,'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1BTX,  MXNTTS,'ADDL  ',1,'CI1BTX')
      CALL MEMMAR(KPCIBTX ,8*MXNTTS,'ADDL  ',1,'CIBTX ')
      CALL MEMMAR(KPCBLTPX,  NSMST,'ADDL  ',2,'CBLTPX')
C
C     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
C
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
C
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
C       info needed for generation of block info
        call memmar(KLCIOIOX,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIOX))
C
        CALL ZBLTP_REL(ISMOST(1,ICSM_SCR),NSMST,IDC,WORK(KPCBLTPX))
C
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
C       number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
C
C       batches of C vector
        ITTSS_ORD = 2
        CALL PART_CIV4(IDC,WORK(KPCBLTPX),WORK(KNSTSO(IATP)),
     &                 WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                 IOFFBLK,IOFFBTC,
     &                 NSMST,LBLOCK,WORK(KLCIOIOX),
     &                 ISMOST(1,ICSM_SCR),
     &                 NBATCH,WORK(KPCLBTX),WORK(KPCLEBTX),
     &                 WORK(KPCI1BTX),WORK(KPCIBTX),0,ITTSS_ORD)
C
C       number of BLOCKS so far
        NBLOCK = IFRMR_REL(WORK(KPCI1BTX),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBTX),1,NBATCH) - 1
*
C       WRITE(LUWRT,*) ' IMS2, NBLK_MS2_C(IMS2) is',
C    &                   IMS2,NBLOCK - NBLOCKP
        NBLK_MS2_C(IMS2) = NBLOCK - NBLOCKP
*
      END DO
*     ^ End of loop over MS2 values
C
9999  call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLKF')
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE Z_BLKFO_REL_1(IDC,NMS2VAL,ISPC,ISM,KPCIBT,KPCBLTP,
     &                         KNODELIST,KBLOCKLIST,KRCCTOS,
     &                         NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                         NBAT_MS2,IBAT_MS2,I_SET_L2BLOCK,
     &                         I_USE_PC,NPARBLOCK,IRC,IT_TTPL,
     &                         IGROUPLIST)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT  : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
*
*.output
      DIMENSION IT_TTPL(*)
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*) ' Output from Z_BLKFO_REL_1  '
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
*     allowed length of each batch
CSK   IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .ge. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
CSK     WRITE(LUWRT,'(1X,A,1X,I20)')
CSK  & '  Current available free memory in double words:',L2BLOCK
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK_SAVE_V = 0
        L2BLOCK_SAVE_V = L2BLOCK
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*       check if L2BLOCK is a positive integer
        IF( L2BLOCK .le. 0 )THEN
*         set back to default...
          WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
          WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
          ISMEMFAC = 9
          L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4
        END IF
*
csk     WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk  &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
*
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI)
*
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*       reset to L_COMBI because that is already enough
        L2BLOCK = LTEST_BLOCK
*
      END IF
*
*     ... set LBLOCK value
*
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
*
      IF( I_SET_L2BLOCK .eq. 2 ) RETURN
*
*     pointers to output arrays
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTP ')
*.    ^ These should be preserved after exit so put mark for 
*       flushing after
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLKS')
*
*     moved here since they are not in use for GASCI_REL
      call memmar(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BT ')
*
      IDUM = 0
*
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*. Info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
*. Number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLKS')
*
*     start second part...
*
      NUM_BLOCKS  = NBLOCK
      NUM_BLOCKS2 = NBLOCK * IRC
*
      CALL MEMMAR(KNODELIST,  NUM_BLOCKS,'ADDS  ',1,'NODELS')
      CALL MEMMAR(KBLOCKLIST, NUM_BLOCKS,'ADDS  ',1,'BLCKLS')
      CALL MEMMAR(KRCCTOS,   NUM_BLOCKS2,'ADDS  ',1,'RCCTOS')
*.    ^ These should be preserved after exit so put mark for
*       flushing after
*
      MIN2  = - 2
*
      CALL ISETVC(WORK(KNODELIST),MIN2,NUM_BLOCKS)
      CALL IZERO(WORK(KBLOCKLIST),NUM_BLOCKS)
      CALL IZERO(WORK(KRCCTOS),NUM_BLOCKS2)
*     initialize block array
      CALL INI_BLOCKL(WORK(KBLOCKLIST),NBLOCK,WORK(KPCIBT))
*
      
#if defined (VAR_MPI2)
      COMPDISTL     = .TRUE.
      COMPRCCTOS    = .TRUE.
      CALL BLOCK_DISTR_DRV_REL(WORK(KNODELIST),WORK(KBLOCKLIST),
     &                         WORK(KBLOCKLIST),WORK(KRCCTOS),
     &                         IT_TTPL,IGROUPLIST)
#endif
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE Z_BLKFO_REL_2(IDC,NMS2VAL,ISPC,ISM,KPCIBT,KPCBLTP,
     &                         NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                         NBAT_MS2,IBAT_MS2,I_SET_L2BLOCK,
     &                         I_USE_PC,NPARBLOCK,IRC)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT  : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
*
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*) ' Output from Z_BLKFO_REL_2  '
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLK2')
*     moved here since they are not in use for GASCI_REL
      call memmar(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BT ')
*
      IDUM = 0
*
*     allowed length of each batch
CSK   ELSE IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .eq. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK_SAVE_V = 0
        L2BLOCK_SAVE_V = L2BLOCK
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*       check if L2BLOCK is a positive integer
        IF( L2BLOCK .le. 0 )THEN
*         set back to default...
          WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
          WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
          ISMEMFAC = 9
          L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4 
        END IF
*
csk        WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk     &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
*
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI)
*
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*       reset to L_COMBI because that is already enough
        L2BLOCK = LTEST_BLOCK
*
      END IF
*
*     ... set LBLOCK value
*
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
*
CACTU      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
CTEST      LBLOCK  = MIN(L0BLOCK,L2BLOCK)
*          old versions works without L0BLOCK and L2BLOCK
CVOLD     LBLOCK = MAX(LBLOCK,LCSBLK)
*
      NBATCH = 0
      NBLOCK = 0
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*. Info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
*. Number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
       WRITE(LUWRT,*) ' Total number of batches', NBATCH
       WRITE(LUWRT,*) ' Total number of blocks ', NBLOCK
*
       WRITE(LUWRT,*) ' Number of blocks and batches per MS2 value '
       WRITE(LUWRT,*) ' ==========================================='
       WRITE(LUWRT,*)
       WRITE(LUWRT,*) ' IMS2    MS2      NBLOCK      NBATCH '
       WRITE(LUWRT,*) ' ===================================='
       DO IMS2 = 1, NMS2VAL
         IATP = IST_FOR_DT(1,IMS2)
         IBTP = IST_FOR_DT(2,IMS2)
         IIMS2 = NELEC(IATP)-NELEC(IBTP)
         WRITE(LUWRT,'(4(3X,I4))') 
     &              IMS2,IIMS2,NBLK_MS2(IMS2),NBAT_MS2(IMS2)
       END DO
*
      END IF
#endif
*
*     transfer information to common block for KRMC - LUCIAREL
*     communication:
      if (INTIMP.eq.6.or.INTIMP.eq.7) then
        IF( MYPROC .eq. MASTER ) THEN
          call civcinf_mc(WORK(KPCIBT),WORK(KPCLBT),WORK(KPCI1BT),
     &                    WORK(KPCLEBT),NBAT_MS2,NMS2VAL,NBLOCK)
        END IF
      end if
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLK2')
*
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE Z_BLKFO_XPROP(IDC,NMS2VAL,ISPC,KNODELIST,KBLOCKLIST,
     &                         KRCCTOS,KXSYMDISTR,IXSYMBLK,IXSYMPAIRS,
     &                         IRC,IGROUPLIST)
C***********************************************************************
C
C     construct information about batch and block structure of all 
C     CI spaces defined by ISPC(=1) and ISM used in property run.
C
C     Output is given in the form of pointers to vectors in WORK
C
C     based on Z_BLKFO_REL. 
C
C     Written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "parluci.h"
#include "infpar.h"
#include "krciprop.h"
#if defined (VAR_MPI2)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
      DIMENSION IXSYMBLK(*)
C
C     min length of each batch
      L0BLOCK = MXSOOB
C
C     L2BLOCK = max memory for c-vec and sigma-vec from
C     ( curr. total mem - 3 mio for string info etc.) / 4
C
      L2BLOCK = 0
C
      CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
      L2BLOCK        = LMEMFREE_PTR
C
C     we want to keep three blocks in memory at the same time
C     CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
C     division by a factor of 4 = safety!
C
      L2BLOCK_SAVE_V = 0
      L2BLOCK_SAVE_V = L2BLOCK
      L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
C     check if L2BLOCK is a positive integer
      IF( L2BLOCK .le. 0 )THEN
C       set back to default...
        WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
        WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
        ISMEMFAC = 9
        L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4
      END IF
C
C
      L2BLOCK = MIN(LCSBLK,L2BLOCK)
C
C     test if L2BLOCK has become too large...
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI_MAX)
C
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
C       reset to L_COMBI because that is already enough
C       WRITE(LUWRT,*) '  L2BLOCK,LTEST_BLOCK,L_COMBI_MAX ',
C    &                    L2BLOCK,LTEST_BLOCK,L_COMBI_MAX
        L2BLOCK = LTEST_BLOCK
      END IF
C
C     ... set LBLOCK value
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
C     WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK,LBLOCK ',
C    &                  L0BLOCK,L2BLOCK,LCSBLK,LBLOCK
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'ZBXPRP')
C
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBTX ')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTPX')
      CALL MEMMAR(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BTX')
C
      NBLOCK_MAX = 0
      IXSYMACT   = 0
C
C     WRITE(LUWRT,*) ' ... before IJSYM loop: LBLOCK = ', LBLOCK
      DO 100 IJSYM = 1, MXPROPKRCI_SYM
C
        IF( IXSYMLIST(IJSYM) .eq. 0 ) GOTO 100
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NOCTPA =  NOCTYP(IATP)
          NOCTPB =  NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
C         info needed for generation of block info
          CALL MEMMAR(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CXIOIO')
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
          CALL ZBLTP_REL(ISMOST(1,IJSYM),NSMST,IDC,WORK(KPCBLTP))
C
          IF(IMS2.EQ.1) THEN
            IOFFBLK = 1
            IOFFBTC = 1
          ELSE
            IOFFBLK = NBLOCK + 1
            IOFFBTC = NBATCH + 1
          END IF
C         number of batches and blocks obtained until now
          IF(IMS2.EQ.1) THEN
            NBATCHP = 0
            NBLOCKP = 0
          ELSE
            NBATCHP = NBATCH
            NBLOCKP = NBLOCK
          END IF
C         batches of C vector
          ITTSS_ORD = 2
C
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,IJSYM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
C
C         number of BLOCKS so far
          NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &           + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
C
        END DO
C       ^ End of loop over MS2 values
        IXSYMBLK(IJSYM)     = NBLOCK
        NBLOCK_MAX          = MAX(NBLOCK_MAX,NBLOCK)
        IXSYMACT            = IXSYMACT + 1
        IXSYMLIST(IXSYMACT) = IJSYM
 100  CONTINUE
C     ^ End of loop over different symmetries
C     WRITE(LUWRT,*) ' ... after IJSYM loop: LBLOCK = ', LBLOCK

      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ZBXPRP')
C
C     start second part...
C
      IF( NBLOCK_MAX .le. 0 )THEN
          WRITE(LUWRT,*) ' *** ERROR in Z_BLKFO_XPROP: max number of'// 
     &                   ' TTSS blocks is zero.'
          CALL QUIT('*** Z_BLKFO_XPROP: max number of TTSS blocks is
     &               zero ***')
      END IF
C
      NUM_BLOCKS  = NBLOCK_MAX
      NUM_BLOCKS2 = NBLOCK_MAX * IRC
C
      CALL MEMMAR(KNODELIST, NUM_BLOCKS *IXSYMACT  ,'ADDS  ',1,'NODELS')
      CALL MEMMAR(KBLOCKLIST,NUM_BLOCKS *IXSYMACT  ,'ADDS  ',1,'BLCKLS')
      CALL MEMMAR(KXSYMDISTR,            IXSYMACT  ,'ADDS  ',1,'IXSYMD')
      CALL MEMMAR(KRCCTOS,   NUM_BLOCKS2*NXSYMPAIRS,'ADDS  ',1,'RCCTOS')
C     ^ keep them after exit
C
      CALL ISETVC(WORK(KNODELIST),-2,NUM_BLOCKS *IXSYMACT)
      CALL IZERO(WORK(KBLOCKLIST),   NUM_BLOCKS *IXSYMACT)
      CALL IZERO(WORK(KXSYMDISTR),               IXSYMACT)
      CALL IZERO(WORK(KRCCTOS),      NUM_BLOCKS2*NXSYMPAIRS)
C
C     mark for local memory
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'ZBXPR2')
C
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBTX2')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTX2')
      CALL MEMMAR(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBTX2')
      CALL MEMMAR(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBX2')
      CALL MEMMAR(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BX2')
C
C     WRITE(LUWRT,*) ' ... before IXSYMACT loop: LBLOCK = ', LBLOCK
      DO 200 I = 1, IXSYMACT
C
        IJSYM = IXSYMLIST(I)
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NOCTPA =  NOCTYP(IATP)
          NOCTPB =  NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
C         info needed for generation of block info
          CALL MEMMAR(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
          CALL ZBLTP_REL(ISMOST(1,IJSYM),NSMST,IDC,WORK(KPCBLTP))
C
          IF(IMS2.EQ.1) THEN
            IOFFBLK = 1
            IOFFBTC = 1
          ELSE
            IOFFBLK = NBLOCK + 1
            IOFFBTC = NBATCH + 1
          END IF
C         number of batches and blocks obtained until now
          IF(IMS2.EQ.1) THEN
            NBATCHP = 0
            NBLOCKP = 0
          ELSE
            NBATCHP = NBATCH
            NBLOCKP = NBLOCK
          END IF
C         batches of C vector
          ITTSS_ORD = 2
C
C         WRITE(LUWRT,*) 'IJSYM is',IJSYM
C
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,IJSYM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
C
C         number of BLOCKS so far
          NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &           + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
C
        END DO
C       ^ End of loop over MS2 values
C
C       initialize symmetry irrep block array
        CALL INI_BLOCKL_PRP(WORK(KBLOCKLIST),NUM_BLOCKS,I,
     &                      NBLOCK,WORK(KPCIBT))
 200  CONTINUE
C     ^ End of loop over different symmetries
C     WRITE(LUWRT,*) ' ... after IXSYMACT loop: LBLOCK = ', LBLOCK

      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ZBXPR2')
C
#if defined (VAR_MPI2)
C     WRITE(LUWRT,*) ' ... before BLOCK_DISTR_PROP: LBLOCK = ', LBLOCK
      CALL BLOCK_DISTR_PROP(WORK(KNODELIST),WORK(KBLOCKLIST),
     &                      WORK(KRCCTOS),IXSYMPAIRS,
     &                      WORK(KXSYMDISTR),IGROUPLIST)
C     WRITE(LUWRT,*) ' ... after BLOCK_DISTR_PROP: LBLOCK = ', LBLOCK
#endif
      END
#if defined (VAR_MPI2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE BLOCK_DISTR_PROP(NPARBLOCK_X,IBLOCKL_X,RCCTOS_X,
     &                            IXSYMPAIRS,ISYM_DISTR,IGROUPLIST)
C***********************************************************************      
C
C     driver routine for symmetry irrep dependent block distribution 
C
C     written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "krciprop.h"
#include "mpif.h"
#include "parluci.h"
#include "infpar.h"
#include "cands.inc"
#include "mxpdim.inc"
#include "symm.inc"
      DIMENSION ISTAT(MPI_STATUS_SIZE), ISYM_DISTR(IXSYMACT)
      DIMENSION NPARBLOCK_X(NUM_BLOCKS,*),IBLOCKL_X(NUM_BLOCKS,*)
      DIMENSION IXSYMPAIRS(MXPROPKRCI_SYM,*)
      INTEGER   RCCTOS_X(NUM_BLOCKS2,*)
C
      IXPAIR = 0
C
C     loop over all "active" symmetry irrep pairs
C     for simplicity we assume the same even distribution for a given 
C     sigma irrep no matter from which T_X x C pair it will be generated
C
C     COMPDISTL and COMPRCCTOS on common block in include/parluci.h.
C
C     WRITE(LUWRT,*) ' block distribution for sum_(act. irreps) =',
C    &                 IXSYMACT
C
      DO I = 1, IXSYMACT
        IISYM = IXSYMLIST(I)
        DO J = 1, IXSYMACT
           IJSYM = IXSYMLIST(J)
           IF( IXSYMPAIRS(IISYM,IJSYM) .eq. 1 )THEN
             IXPAIR = IXPAIR + 1
C
C            reset IXSYMPAIRS to current pair number
C            ---------------------------------------
             IXSYMPAIRS(IISYM,IJSYM) = IXPAIR
             IF( ISYM_DISTR(J) .eq. 0 ) THEN
               COMPDISTL     = .TRUE.
               COMPRCCTOS    = .TRUE.
               ISYM_DISTR(J) = 1
             ELSE
               COMPDISTL  = .FALSE.
               COMPRCCTOS = .TRUE.
             END IF
C            set correct symmetry for C and sigma
             ICSM = IISYM
             ISSM = IJSYM
C
C            set T operator arrays
             ISYM_T    = IDBGMULT(INVELM(ICSM),ISSM)
             IPRHAM_CI = 0
             NO_TTWO   = 1
C	 WRITE(LUWRT,*) ' inside loop I,J (1): LBLOCK',I,J,LBLOCK
             CALL SET_HOP_DBG(IHTYPE_X,NO_TTWO,ISYM_T,IPRHAM_CI)
C	 WRITE(LUWRT,*) ' inside loop I,J (2): LBLOCK',I,J,LBLOCK
C
             CALL BLOCK_DISTR_DRV_REL(NPARBLOCK_X(1,J),IBLOCKL_X(1,J),
     &                                IBLOCKL_X(1,I),RCCTOS_X(1,IXPAIR),
     &                                IDUMMY,IGROUPLIST)
C	 WRITE(LUWRT,*) ' inside loop I,J (3): LBLOCK',I,J,LBLOCK
           END IF
        END DO
      END DO
      END
#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE GET_SC_PAIRS(IXSYMPAIRS,IXSYM,NX_OP,ISYMOPX,
     &                        IXROOTS,IXROOTSYM)
C***********************************************************************      
C
C     obtain all relevant sigma/C symmetry pairs
C
C     written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************      
      IMPLICIT REAL*8(A-H,O-Z)
#include "parluci.h"
#include "mxpdim.inc"
#include "symm.inc"
#include "krciprop.h"
      DIMENSION IXSYMPAIRS(IXSYM,IXSYM), ISYMOPX(*),IXROOTSYM(*)
C
      IONE = 1
C
      DO 100 I = 1, NX_OP
         IF( ISYMOPX(I) .eq. 0 ) GOTO 100
         DO J = 1, IXROOTS
           ICSM = IXROOTSYM(J)
           ISSM = IDBGMULT(ICSM,ISYMOPX(I))
           IXSYMPAIRS(ICSM,ISSM) = IONE
         END DO
 100  CONTINUE
C
C     number of symmetry pairs
      NXSYMPAIRS = 0
      DO I = 1, IXSYM
        DO J = 1, IXSYM
          IF( IXSYMPAIRS(I,J).eq. IONE ) NXSYMPAIRS = NXSYMPAIRS + 1
        END DO
      END DO
c     WRITE(LUWRT,*) ' GET_SC_PAIRS: number of symmetry pairs = ',
c    &                 NXSYMPAIRS
C
       END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
* Codes for general symmetry handling
*
*                - ZSTINF : generate /STINF/ info on strings and mapping
*                - MEMSTR : allocates memory for string information
*                - WEIGHT : Weights for strings
*                - NSTRSO : Number of strings per sym and class
*                - ZBASE  : offset arrays for strings
*                - ZSMCL  : symmetry and class for each string
*                - GENSTR : Generate strings ordered by sym and class
*                - MEMEXT : Memory for external blocks
*
      SUBROUTINE ZBLTP_REL(ISMOST,MAXSYM,IDC,ICBLTP)
*
* Generate vector ICBLTP giving type of each block
*
* Simplified version : IDC = 3, IDC = 4 path eliminated
*
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
#include "ipoist8.inc"
      DIMENSION ISMOST(*)
*. Output
      DIMENSION ICBLTP(*)
*
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF(IBSYM .EQ. 0 ) GOTO 100
        IF(IDC.EQ.2.AND.IBSYM.GT.IASYM) THEN
*.Symmetry block excluded
          ICBLTP(IASYM) = 0
        ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM).OR.IDC.EQ.1) THEN
*.Complete symmetry block included
          ICBLTP(IASYM) = 1
        ELSE
*.Lower half  symmetry block included
          ICBLTP(IASYM) = 2
        END IF
  100 CONTINUE
*
      NTEST = 0
#if defined LUCI_DEBUG
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Block type of symmetry blocks '
         CALL IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZIRAT_REL
*
* Ratio between real and integer
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
#include "SysDef.inc"
#include "machine.inc"
*. Output
#include "irat.inc"
*. For CRAY or explicit 64 bit architecture : IRAT = 1
c     IF(MACHINE(1:4).EQ.'CRAY'.OR.MACHINE(1:2).EQ.'64') THEN
c       IRAT = 1
c     ELSE
*. For all other architectures : 2 integers per real
c       IRAT = 2
c     END IF
*
CTF No variability in thsi version.
*   We always use RTOI = 1, but in some MOLCAS interface routines,
*   RTOI is set to 2.
c     IRAT = RTOI
      IRAT = 1
*
C     WRITE(6,*)
C     WRITE(6,'(1H ,6X,A,I2)')
C    &'Ratio between Integer and Real word length ', IRAT
C     WRITE(6,*)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZNELFSPGP_REL(NTESTG)
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
#include "gasstr.inc"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZNONAB_REL(INVCEN,MAXMLO,NSMOB,IPRNT)
*
*
* ============================
* Set up common block /NONAB/
* ============================
*
*========
* Input :
*========
*      INVCNT :inversion center is present(1), absent(0)
*      MAXMLO : Largest ML value of orbitals
*      NSMOB  : Number of symmetries of orbitals
*      Contents of common block /STRINP/,/ORBINP/
*=========
* output :
*=========
*======================
* Orbital Information
*======================
*      NORASM : Number of orbitals per abelian symmetry
*      MNMLOB : Smallest ML of orbitals
*      MXMLOB : largest ML of orbitals
*      NMLOB  : number of ML values for orbitals
*
*======================
* String Information
*======================
*      MNMLST : smallest ML of any strings
*      MXMLST : largest ML of any strings
*      NMLST  : Number of ML values of strings
*      NSMST  : Number of symmetries of strings
*
*==============================
* Single excitation Information
*==============================
*      MNMLSX : SMALLEST ML OF SINGLE EXCITATION
*      MXMLSX : LARGEST ML OF SINGLE EXCITATIONS
*      NMLSX  : NUMBER OF ML VALUES FOR SINGLE EXCITATIONS
*      NSMSX  : NUMBER OF SYMMETRIES FOR SINGLE EXCITATIONS
*=============================================
* External configurations (upto 4 electrons )
*=============================================
*      MNMLXT : SMALLEST ML OF External configurations
*      MXMLSX : LARGEST ML OF external configurations
*      NMLXT  : NUMBER OF ML VALUES FOR ext. configurations
*      NSMXT  : NUMBER OF SYMMETRIES FOR ext. configurations
*
* =============
* General input
* =============
*
*./ORBINP/
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     PARAMETER(MXPORB=500)
C     INTEGER ORBSM
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
*./CSM/
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*./STRINP/
C     PARAMETER(MXPSTT=40)
C     COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),
C    &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),
C    &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),
C    &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
#include "strinp.inc"
* =======
*. Output
* =======
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
*. Inversion symmetry
      IF( INVCEN .EQ. 0 ) THEN
        INVCNT = .FALSE.
        NIG = 1
      ELSE
        INVCNT = .TRUE.
        NIG = 2
      END IF
*
** 1 : Information about orbitals
*
      MXMLOB = MAXMLO
      MNMLOB =-MAXMLO
      NMLOB = MXMLOB - MNMLOB + 1
*. Number of orbitals per symmetry
      DO 10 ISYM = 1, NSMOB
        NORASM(ISYM) = IFREQ(ISMFTO,ISYM,NACOB)
   10 CONTINUE
#if defined LUCI_DEBUG
      IF( NTEST.GE. 2 ) THEN
        WRITE(6,*) ' NORASM '
        CALL IWRTMA(NORASM,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' MNMLOB,MXMLOB ',MNMLOB,MXMLOB
        WRITE(6,*) ' NMLOB, NSMOB ',NMLOB,NSMOB
      END IF
#endif
*
**  2. Information about strings
*
      MXMLST = 0
      MNMLST = 0
      DO 50 ITYPE = 1, NSTTYP
        IEL = NELEC(ITYPE)
*
        MXMLTP = 0
        DO 40 IML = MXMLOB,MNMLOB,-1
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MXMLTP = MXMLTP + IEL2*IML
          IEL = IEL - IEL2
   40   CONTINUE
        MXMLST = MAX(MXMLST,MXMLTP)
*
        MNMLTP = 0
        IEL = NELEC(ITYPE)
        DO 45 IML = MNMLOB,MXMLOB
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MNMLTP = MNMLTP + IEL2*IML
          IEL = IEL - IEL2
   45   CONTINUE
        MNMLST = MIN(MNMLST,MNMLTP)
   50 CONTINUE
*
      NMLST  = MXMLST - MNMLST + 1
      NSMST  = NIG * NMLST
*
#if defined LUCI_DEBUG
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' MXMLST,MNMLST,NSMST'
        WRITE(6,*)   MXMLST,MNMLST,NSMST
      END IF
#endif
*
** 3. Information about single excitations
*
      MNMLSX = MNMLOB - MXMLOB
      MXMLSX = MXMLOB - MNMLOB
      NMLSX  = MXMLSX - MNMLSX +1
      NSMSX  = NIG * NMLSX

#if defined LUCI_DEBUG
      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLSX,NSMSX,MNMLSX ',NMLSX,NSMSX,MNMLSX
      END IF
#endif
*
** 4 : External configurations(double excitations)
*
      MXMLDX = 4*MAXMLO
      MNMLDX = -4*MAXMLO
      NMLDX  = MXMLDX - MNMLDX + 1
      NSMDX  = NIG * NMLDX
#if defined LUCI_DEBUG
      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLDX,NSMDX,MNMLDX ',NMLDX,NSMDX,MNMLDX
      END IF
#endif
*
** 5 : Determinants
*
      MXMLCI =  2*MXMLST + MXMLDX
      MNMLCI = - MXMLCI
      NMLCI = 2 * MXMLCI + 1
      NSMCI = NIG * NMLCI
*
*.6 Total symmetrix single excitation and external
*
      ITSSX = 0 - MNMLSX + 1
      ITSDX = 0 - MNMLDX + 1

#if defined LUCI_DEBUG
      IF ( NTEST .GE. 1 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of orbitals     .. ', NSMOB
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of strings      .. ', NSMST
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of single excit. . ', NSMSX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of double excit. . ', NSMDX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of determinants .. ', NSMCI
        WRITE(6,*)
*
        WRITE(6,*) ' Total symmetric single excitation .. ',ITSSX
        WRITE(6,*) ' Total symmetric double excitation .. ',ITSDX
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSYM1_REL(NIRREP,IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
* ( trivial, written for compatibility with higher point groups)
*
#include "ipoist8.inc"
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
C     PARAMETER ( MXPOBS = 20 )
#include "mxpdim.inc"
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      NSMXT = NIRREP
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1

*
C     COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
      CALL ICPMT2(SYMPRO,ADASX,8,8,MXPOBS,MXPOBS,1)
      CALL ICPMT2(SYMPRO,ADSXA,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,ASXAD,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSYM2_REL(IPRNT)
*
* Symmetry connecting arrays
*
* ======
*. Input
* ======
*
*./LUCINP/
C     INTEGER PNTGRP,CITYP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20 )
C     PARAMETE
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONABR/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*./CSM
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
* ======
*.Output
* ======
*
*./CSMPRD/
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
**. ADASX,ASXAD,ADSXA
      CALL ISETVC(ADASX,0,MXPOBS**2)
      CALL ISETVC(ASXAD,0,2*MXPOBS**2)
      CALL ISETVC(ADSXA,0,2*MXPOBS**2)
*
      DO 100 ISM = 1, NSMOB
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISM,'OB',2)
        DO 90 JSM = 1, NSMOB
          CALL MLSM_REL(JML,JPARI,JSM,'OB',2)
*.a+ i a j symmetry
          IJML = IML - JML
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJSM = (IJPARI-1)*NMLSX + IJML - MNMLSX + 1
          ADASX(ISM,JSM) = IJSM
          ASXAD(JSM,IJSM) = ISM
          ADSXA(ISM,IJSM) = JSM
   90   CONTINUE
  100 CONTINUE
*.SXSXDX,SXDXSX
      DO 200 ISX = 1, NSMSX
C       MLSM_REL(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM_REL(IML,IPARI,ISX,'SX',2)
        DO 190 JSX = 1, NSMSX
          CALL MLSM_REL(JML,JPARI,JSX,'SX',2)
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJML = IML + JML
          IJSM = (IJPARI-1)*NMLDX+IJML - MNMLDX + 1
          SXSXDX(ISX,JSX) = IJSM
          SXDXSX(ISX,IJSM) = JSX
  190   CONTINUE
  200 CONTINUE
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ADASX '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADASX,NSMOB,NSMOB,MXPOBS,MXPOBS)
         WRITE(6,*) ' ASXAD '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ASXAD,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' ADSXA '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADSXA,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXSXDX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXSXDX,NSMSX,NSMSX,2*MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXDXSX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXDXSX,NSMSX,NSMDX,2*MXPOBS,4*MXPOBS)
      END IF
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSYMDG(IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
*
*
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "csm.inc"
#include "symm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
      NSMSX = NIRR_DG
      NSMDX = NIRR_DG
      NSMST = NIRR_DG
      NSMCI = NIRR_DG
      NSMXT = NIRR_DG
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1

*
      IZERO = 0
      CALL ISETVC(ADASX,IZERO,MXPOBS**2)
      CALL ISETVC(ASXAD,IZERO,2*MXPOBS**2)
      CALL ISETVC(ADSXA,IZERO,2*MXPOBS**2)
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
*A+ A => SX
          ADASX(IREP,JREP) = IDBGMULT(IREP,IADJSYM(JREP))
        END DO
      END DO
*
      DO IREP = 1, NIRR_DG
        DO JREP = 1, NIRR_DG
          DO IJSM = 1, NIRR_DG
            IF(ADASX(IREP,JREP).EQ.IJSM) THEN
              ADSXA(IREP,IJSM) = JREP
              ASXAD(JREP,IJSM) = IREP
            END IF
          END DO
        END DO
      END DO

      CALL ICPMT2(IDBGMULT,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(IDBGMULT,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
