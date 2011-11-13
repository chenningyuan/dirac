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
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAS_REL(IMK2NUM,NMS2VAL,MS2VAL,NAEL,
     &           IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,NIRREP,
     &           H1A,H1B,XA,XB,SCR,RJ,RK,
     &           RK2,NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,IASPGPTP,IBSPGPTP,
     &           ISPINFREE,NPARBLOCK,ILASTBLK,IDIA_OFFSET)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* Fixed bug in first exchange term.    Timo Fleig, 07.Aug.98
* Fixed bug in third integral term.    Timo and Jeppe, 25.Aug.98
* Correct sym-type ordering for all integrals implemented.
*                                      Timo Fleig, fall 2000
*
* IHAM12 = 1 => only one-body part
*        = 2 =>      one+two-body part
*
* separate alpha and beta 1-electron integrals  May 1998
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H1A(NORB), H1B(NORB), NPARBLOCK(*)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
      dimension MS2VAL(NMS2VAL)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION RK2(NTOOB,NTOOB)
*
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
#if defined (VAR_MPI2)
      INTEGER(KIND=MPI_OFFSET_KIND) IDIA_OFFSET
#endif
*. Output
      DIMENSION DIAG(*)
*
      call qenter('GASDIAS')
*
      NTEST =  0000
      NTEST = MAX(NTEST,IPRNT)
*
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*
      IF (NTEST.GE.20) THEN
        if (IMK2NUM.eq.1) then
          WRITE(6,*)
     &    ' Diagonal one electron integrals, unbarred'
          CALL WRTMAT(H1A,1,NORB,1,NORB)
          WRITE(6,*)
     &    ' Diagonal one electron integrals,   barred'
          CALL WRTMAT(H1B,1,NORB,1,NORB)
          WRITE(6,*) ' Core energy ', ECORE
          IF(IHAM12.EQ.2) THEN
            WRITE(6,*) ' Coulomb and exchange integrals '
            CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
            WRITE(6,*)
            CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
            if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
              WRITE(6,*)
              CALL WRTMAT(RK2,NTOOB,NTOOB,NTOOB,NTOOB)
            end if
          END IF
        end if
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO

        WRITE(6,*) ' Number of alpha strings , type 1 '
        CALL IWRTMA(NSSOA,1,NSMST,1,NSMST)
        WRITE(6,*) ' Number of beta strings , type 1 '
        CALL IWRTMA(NSSOB,1,NSMST,1,NSMST)
*
      END IF
*
*  Diagonal elements according to Timo's formula
*
*   DIAG(IDET)
*          = sum_I (NIU*HII + NIB*HI_I_)
*            + sum_KM 1/2 * NKU*NMU [(KK|MM)     - (KM|MK)    ]
*                   + 1/2 * NKB*NMB [(K_K_|M_M_) - (K_M_|M_K_)]
*                   +       NKU*NMB [(KK|M_M_)   - (KM_|M_K)  ]
*
      IDET = 0
      ITDET = 0
*
      DO IBLK = 1, NBLOCK
#if defined (VAR_MPI2)
        JBLOCK = IBLK + ILASTBLK
        IF( NPARBLOCK( JBLOCK ) .ne. MYPROC ) GOTO 2425 
#endif
        if (NTEST.ge.50) then
          write(6,*) '-----------------------------------'
          write(6,'(A14,I8,A13)') 
     &               ' Block number ',IBLK,' in diagonal '
          write(6,*) '-----------------------------------'
        end if
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
        if (NTEST.ge.10)
     &    WRITE(6,'(A21,4I4)') 
     &    ' IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
*
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
*
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
C?        WRITE(6,*)  ' IB = ', IB
          IBREL = IB - IBSTRT + 1
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
C?        WRITE(6,*)' IASTRT,IASTOP', IASTRT,IASTOP
*
          DO IA = IASTRT,IASTOP
C?          WRITE(6,*)  ' IA = ', IA
            IDET = IDET + 1
            ITDET = ITDET + 1
*
            X = ECORE
*
            DO IEL = 1, NAEL
* Occupation of alpha spinors
              X = X + H1A(IASTR(IEL,IA))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1A(',IASTR(IEL,IA),')',
     &                       H1A(IASTR(IEL,IA))
              end if
            END DO
*
            DO IEL = 1, NBEL
* Occupation of beta spinors
              X = X + H1B(IBSTR(IEL,IB))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1B(',IBSTR(IEL,IB),')',
     &                       H1B(IBSTR(IEL,IB))
              end if
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NAEL
                X = X + 0.5*
     &            RJ(IASTR(IEL,IA),IASTR(JEL,IA))
     &                - 0.5*
     &            RK(IASTR(IEL,IA),IASTR(JEL,IA))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IASTR(IEL,IA),
     &                          IASTR(JEL,IA),
     &                          ') = ',
     &            0.5*RJ(IASTR(IEL,IA),IASTR(JEL,IA)),
     &                          '- 0.5 K (',
     &                          IASTR(IEL,IA),
     &                          IASTR(JEL,IA),') =',
     &            0.5*RK(IASTR(IEL,IA),IASTR(JEL,IA))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NBEL
              DO JEL = 1, NBEL
               X = X + 0.5*
     &            RJ(IBSTR(IEL,IB),IBSTR(JEL,IB))
     &               - 0.5*
     &            RK(IBSTR(IEL,IB),IBSTR(JEL,IB))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IBSTR(IEL,IB),
     &                          IBSTR(JEL,IB),
     &                          ') = ',
     &            0.5*RJ(IBSTR(IEL,IB),IBSTR(JEL,IB)),
     &                          '- 0.5 K (',
     &                          IBSTR(IEL,IB),
     &                          IBSTR(JEL,IB),') =',
     &            0.5*RK(IBSTR(IEL,IB),IBSTR(JEL,IB))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NBEL
                X = X + 
     &            RJ(IASTR(IEL,IA),IBSTR(JEL,IB))
                if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
                  X = X - 
     &            RK2(IASTR(IEL,IA),IBSTR(JEL,IB))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding J(',
     &                          IASTR(IEL,IA),
     &                          -IBSTR(JEL,IB),
     &                          ') = ',
     &            RJ(IASTR(IEL,IA),IBSTR(JEL,IB)),
     &                          '- K2 (',
     &                          IASTR(IEL,IA),
     &                          -IBSTR(JEL,IB),') =',
     &            RK2(IASTR(IEL,IA),IBSTR(JEL,IB))
                  end if
                end if
              END DO
            END DO
*
            if (NTEST.ge.100) then
               write(6,*) 'Element ITDET and value X:'
               write(6,*) ITDET,X
               write(6,*) 'offset IDET:'
               write(6,*) IDET
               WRITE(6,*) ' alpha and beta string : '
               CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
               CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
               WRITE(6,*)
            end if
            DIAG(IDET) = X
            IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
2345    continue
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMAT(DIAG,1,IDET,1,IDET)
          END IF
#if defined (VAR_MPI2)
CSK          IF( IDET .ne. 0 ) THEN
CSK          write(LUWRT,*) ' number of diagonal elements to disc ',IDET
CSK          CALL WRTMATMN(DIAG,1,IDET,1,IDET,LUWRT)
          CALL MPI_FILE_WRITE_AT(IDIA,IDIA_OFFSET,DIAG,IDET,
     &                         MPI_REAL8,ISTAT,IERR)
*         new offset
          IDIA_OFFSET = IDIA_OFFSET + IDET
CSK          END IF
#endif
          IF( .NOT. NO_CDIAF )THEN
            CALL ITODS(IDET,1,-1,LUDIA)
            CALL TODSC(DIAG,IDET,-1,LUDIA)
          END IF
          IDET = 0
        END IF
2425    CONTINUE
      END DO
*     ^ End of loop over blocks
      if (NTEST.ge.500) then
         WRITE(6,*)
     &   ' Number of diagonal elements generated ',ITDET
      end if
*
      IF (NTEST.GE.1000) THEN
        WRITE(6,*) ' CIDIAGONAL FIXME. (write statement)'
c       CALL WRTMAT(DIAG(1),ITDET,1,ITDET,1)
      END IF
*
      call qexit('GASDIAS')
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAT_REL(DIAG,LUDIA,ECORE,ICISTR,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,KPICI1,INTSPACE,T_BUFF,T_BUFF_D,
     &           IDENSI,NIRREP,ISPINFREE,ENVIRO,SIGDEN_ROUTE,IRECO,
     &           NPARBLOCK,CPIC_PASS_D)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
* For DBG-RELCI, March21 98
* Modifed fetch of 1 electron integrals, May 98
* Modified for parallel runs, S. Knecht, May 2007
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "symm.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
      LOGICAL CPIC_PASS_D
      DIMENSION IBLTP(*), T_BUFF_D(*)
      DIMENSION IBLKFO(8,NBLOCK), NPARBLOCK(*)
#if defined (VAR_MPI2)
      INTEGER(KIND=MPI_OFFSET_KIND) IDIA_OFFSET
#endif
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('CIDIA')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*. Largest Number of electrons around
      MAXNEL = 0
      DO IMK2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
        MAXNEL = MAX(MAXNEL,NELEC(IATP),NELEC(IBTP))
      END DO
*. Local memory
      IDUM = 0
      call memmar(KDUM,  IDUM,    'MARK  ',IDUM,'GASDIA')
*
CTF  Total number of orbitals and active orbitals (spinors, actually)
CTF  from BOTH sets, barred and unbarred !
      NTOOBTOT = NTOOB + NTOOB2
      NACOBTOT = NACOB + NACOB2
*
*
      call memmar(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      call memmar(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLK2  ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLSCR2,2*NTOOBTOT**2,'ADDL  ',2,'KLSC2 ')
      call memmar(KLXA  ,NACOBTOT,   'ADDL  ',2,'KLXA  ')
      call memmar(KLXB  ,NACOBTOT,   'ADDL  ',2,'KLXB  ')
      call memmar(KLSCR ,2*NACOBTOT, 'ADDL  ',2,'KLSCR ')
C     CALL MEMMAR(KLH1D ,NACOBTOT, 'ADDL  ',2,'KLH1D ')
      call memmar(KLH1DA,NTOOB, 'ADDL  ',2,'KLH1DA')
      call memmar(KLH1DB,NTOOB, 'ADDL  ',2,'KLH1DB')
*. Space for blocks of strings
      call memmar(KLASTR,MXNSTR*MAXNEL,'ADDL  ',1,'KLASTR')
      call memmar(KLBSTR,MXNSTR*MAXNEL,'ADDL  ',1,'KLBSTR')
      KLSVST = 1
      call memmar(KLRJKA,MXNSTR, 'ADDL  ',2,'KLRJKA')
*
      if (IRECO.eq.1) then
        if( CPIC_PASS_D )then
C
          IOFF_T_D = 1
C
          CALL DCOPY(NTOOB,T_BUFF_D,1,WORK(KLH1DA),1)
          IOFF_T_D = IOFF_T_D + NTOOB
          CALL DCOPY(NTOOB,T_BUFF_D(IOFF_T_D),1,WORK(KLH1DB),1)
          IOFF_T_D = IOFF_T_D + NTOOB
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLJ),1)
          IOFF_T_D = IOFF_T_D + NTOOB**2
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLK),1)
          IOFF_T_D = IOFF_T_D + NTOOB**2
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLK2),1)
C
CCC       WRITE(LUWRT,*) ' T_BUFF_D:'
CCC       CALL WRTMATMN8(T_BUFF_D,1,MY_T_LEN_D,1,MY_T_LEN_D,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLH1DA):'
CCC       CALL WRTMATMN(WORK(KLH1DA),1,NTOOB,1,NTOOB,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLH1DB):'
CCC       CALL WRTMATMN(WORK(KLH1DB),1,NTOOB,1, NTOOB,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLJ):'
CCC       CALL WRTMATMN(WORK(KLJ),1, NTOOB**2,1, NTOOB**2,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLK):'
CCC       CALL WRTMATMN(WORK(KLK),1, NTOOB**2,1, NTOOB**2,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLK2):'
CCC       CALL WRTMATMN(WORK(KLK2),1, NTOOB**2,1, NTOOB**2,LUWRT)
C
        else
          if (IHAM12.eq.1.or.IHAM12.eq.2) then
            call geth1dia4(WORK(KLH1DA),WORK(KLH1DB),WORK(KPICI1),
     &                     INTSPACE,T_BUFF,ENVIRO,SIGDEN_ROUTE)
          end if
*         obtain 2e integrals needed for diagonal
          if (IHAM12.EQ.2.or.IHAM12.eq.3) then
            call get_dia_2e(WORK(KLJ),WORK(KLK),WORK(KLK2),T_BUFF,
     &                      IDENSI,ISPINFREE,NGAS,NGSOB,NGSOB2,
     &                      NGSOBT,NGSOBT2,ISPSF_MO,ENVIRO,
     &                      SIGDEN_ROUTE,IPRDIA)
          end if
        end if
C       ^  IT_SHL?
      end if
C
#if defined (VAR_MPI2)
      IDIA_OFFSET = 0
*
*     calculate offset IDIA_OFFSET for LUDIA file
*
      IDIA_OFFSET = IDIA_OFFSET + MY_DIA_OFF + 
     &              MY_VEC1_IOFF  * ( IRECO - 1 )
CSK      WRITE(LUWRT,*)' Start writing on LUDIA at',IDIA_OFFSET
*
#endif
*
      DO IMK2 = 1, NMS2VAL
        if (NTEST.ge.10) then
          write(LUWRT,*) '=============================='
          write(LUWRT,*) 'Treating MK space in diagonal '
          write(LUWRT,'(A15,I4)') '         MK2 = ',MS2VAL(IMK2)
          write(LUWRT,*) '=============================='
        end if
        IF(IMK2.EQ.1) THEN
          IBOFF = 1
          ILASTBLK = 0
        ELSE
          IBOFF = IBOFF + NBLK_MS2(IMK2-1)
          ILASTBLK = ILASTBLK + NBLK_MS2(IMK2-1)
        END IF
        NNBLK = NBLK_MS2(IMK2)
CSK        WRITE(LUWRT,*)' IBOFF, NNBLK, ILASTBLK',IBOFF, NNBLK,ILASTBLK
*
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
*
        NAEL = NELEC(IATP)
        NBEL = NELEC(IBTP)
*
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
*
        IASPGPTP = IATP
        IBSPGPTP = IBTP
*
*. Offsets for alpha and beta supergroups
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
        IF(NTEST.GE.10) THEN
          WRITE(LUWRT,*) ' ================'
          WRITE(LUWRT,*) ' GASDIA speaking '
          WRITE(LUWRT,*) ' ================'
          WRITE(LUWRT,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
          write(LUWRT,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
          write(LUWRT,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        END IF
*
** Info on block structure of internal state
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
*
        CALL GASDIAS_REL(IMK2,NMS2VAL,MS2VAL,NAEL,WORK(KLASTR),NBEL,
     &         WORK(KLBSTR),NTOOB,DIAG,NSMST,NIRREP,
     &         WORK(KLH1DA),WORK(KLH1DB),
     &         WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &         WORK(KLK),WORK(KLK2),
     &         WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &         LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &         WORK(KLRJKA),IHAM12,IBLTP,NNBLK,IBLKFO(1,IBOFF),
     &         IASPGPTP,IBSPGPTP,ISPINFREE,NPARBLOCK,ILASTBLK,
     &         IDIA_OFFSET)
*
      END DO
*     ^End of loop over MK2 spaces
*
*.Flush local memory
      call memmar(KDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('CIDIA')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_DIA_2E(UUUUC,UUUUE,UBBUE,T_BUFF,IDENSI,ISPINFREE,
     &                      NGAS,NGSOB,NGSOB2,NGSOBT,NGSOBT2,ISPSF_MO,
     &                      ENVIRO,SIGDEN_ROUTE,IPRNT)
*
* Obtain 2electron integrals needed for
* constructing the determinant diagonal
*
*===============================================
* Jeppe Olsen, Timo Fleig
*  SIGDEN route implemented, 
*     Timo Fleig, August 2001
*===============================================
*
* UUUUC(I,J) = ( I I  ! J  J)
* UUUUE(I,J) = ( I J  ! J  I)
* UBBUE(I,J) = ( I J_ ! J_ I)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION UUUUC(NTOOB,NTOOB),UUUUE(NTOOB,NTOOB)
      DIMENSION UBBUE(NTOOB,NTOOB)
*
      character*6 ENVIRO,SIGDEN_ROUTE
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "ctcc.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "integrals_off.inc"
*
      dimension INDGATH(MXPORB),ISPSF_MO(MXNDGIRR,2),IKPERM(4),
     &          NGSOB(MXPOBS,MXPNGAS),NGSOB2(MXPOBS,MXPNGAS),
     &          NGSOBT(MXPNGAS),NGSOBT2(MXPNGAS)
      DIMENSION T_BUFF(*)
*     Input list of integrals (or in SIGDEN on WORK(KT_CC))
*     Input list of integrals (or in SIGDEN in "shared memory"
*     mode on T_BUFF (for "node master" == N_MASTER)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.2) then
        write(6,*)
      ZERO = 0.D0
        write(6,*) ' //////////////////////////////////// '
        write(6,*) '  Entering get_dia_2e                 '
        write(6,*) ' //////////////////////////////////// '
      end if
C     WRITE(LUWRT,*) ' T_BUFF print test'
C     CALL WRTMATMN8(T_BUFF,1,12,1,12,6)
*
* Read in indices from INDICES file
      IDUM = 1
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'GTDIA2')
      NINT12 = N1ELINT + N2ELINT
*
* UUUUC and UUUUE obtained from integral list 3
*
* ( I I ! J J )
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal Coulomb ints. (II|JJ):'
      end if
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        DO ISM = NSMOB/2+1, NSMOB
          DO JSM = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM ISM ! JSM JSM)
            IOFF_DBG = IOFF_DBG_LIST(1,ISM,ISM,JSM,JSM)
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            if (NTEST.ge.50) then
              write(6,*) 'I symmetry is ',ISM
              write(6,*) 'J symmetry is ',JSM
              write(6,*) 'LIST offset is ',IOFF_DBG
              write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
              write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
              write(6,*) 'I offset for this symmetry ',IOFF
              write(6,*) 'J offset for this symmetry ',JOFF
            end if
            DO I = 1, NI
              DO J = 1, NJ
*. address of (II!JJ)
                IJ = IOFF_DBG-1 + ((J-1)*NJ+(J-1))*NI*NI
     &                          +  (I-1)*NI+I
                UUUUC(I-1+IOFF,J-1+JOFF) =
     &          WORK(KINT2D(1,1)-1+IJ)
                IADRIN = KINT2D(1,1)-1+IJ
                if (NTEST.ge.100) then
                   write(6,*) 'selected address:  ',IADRIN
                   write(6,*) '  selected integral: ',WORK(IADRIN)
                end if
              END DO
            END DO
*           ^ End of loop over orbitals
          END DO
        END DO
*       ^ End of loop over symmetries
*
      else if (SIGDEN_ROUTE.eq.'SIGDEN') then
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
* Storage mode for integrals contains "-" sign for beta/barred indices.
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if  
        IGSOFF = 0
        JGSOFF = 0
        ISMOFF = 0
        JSMOFF = 0
        do IGS = 1,NGAS,1
          if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
          do JGS = 1,NGAS,1
            if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
            do ISM = NSMOB/2+1, NSMOB
              if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
              do JSM = NSMOB/2+1, NSMOB
                if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
                NI = NGSOB(ISM,IGS)
                NJ = NGSOB(JSM,JGS)
                IOFF = IGSOFF + ISMOFF
                JOFF = JGSOFF + JSMOFF
                if (NTEST.ge.10) then
                  write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                  write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                  write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                end if
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
                    call geth2dia_sigden(UUUUC(I,J),
     &                                   WORK(KT_CC),T_BUFF,
     &                                   I,J,
     &                                   'COUL',NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS,
     &                                   WORK(KMDIASV),              
     &                                   N1DIAINT,NDIAINT,NTEST)
                  end do
                end do
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
          end do
          JGSOFF = 0
        end do
      end if
*
*
* ( I J ! J I )
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ|JI):'
      end if
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        DO ISM =  NSMOB/2+1, NSMOB
          DO JSM  = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM JSM ! JSM ISM)
            IOFF_DBG = IOFF_DBG_LIST(1,ISM,JSM,JSM,ISM)
            NI = NTOOBS(ISM)
            NJ = NTOOBS(JSM)
            IOFF = IBSO(ISM)
            JOFF = IBSO(JSM)
            if (NTEST.ge.50) then
              write(6,*) 'I symmetry is ',ISM
              write(6,*) 'J symmetry is ',JSM
              write(6,*) 'LIST offset is ',IOFF_DBG
              write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
              write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
              write(6,*) 'I offset for this symmetry ',IOFF
              write(6,*) 'J offset for this symmetry ',JOFF
            end if
            DO I = 1, NI
             DO J = 1, NJ
*. address of (IJ!JI)
               IJ = IOFF_DBG-1 + ((I-1)*NJ+(J-1))*NJ*NI
     &                         +  (J-1)*NI+I
               UUUUE(I-1+IOFF,J-1+JOFF) =
     &         WORK(KINT2D(1,1)-1+IJ)
               IADRIN = KINT2D(1,1)-1+IJ
               if (NTEST.ge.100) then
                  write(6,*) 'selected address:  ',IADRIN
                  write(6,*) '  selected integral: ',WORK(IADRIN)
               end if
             END DO
            END DO
*           ^ End of loop over orbitals
          END DO
        END DO
*
      else if (SIGDEN_ROUTE.eq.'SIGDEN') then
        IKFAC = 1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4
        if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
          IKFAC = -1
* Use Kramers symmetry on second index pair.
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 4
          IKPERM(4) = 3
        end if  
        IGSOFF = 0
        JGSOFF = 0
        ISMOFF = 0
        JSMOFF = 0
        do IGS = 1,NGAS,1
          if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
          do JGS = 1,NGAS,1
            if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
            do ISM = NSMOB/2+1, NSMOB
              if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
              do JSM = NSMOB/2+1, NSMOB
                if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
                NI = NGSOB(ISM,IGS)
                NJ = NGSOB(JSM,JGS)
                IOFF = IGSOFF + ISMOFF
                JOFF = JGSOFF + JSMOFF
                if (NTEST.ge.10) then
                  write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                  write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                  write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                end if
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
                    call geth2dia_sigden(UUUUE(I,J),
     &                                   WORK(KT_CC),T_BUFF,
     &                                   I,J,
     &                                   'EXC1',NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS,
     &                                   WORK(KMDIASV),              
     &                                   N1DIAINT,NDIAINT,NTEST)
                  end do
                end do
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
          end do
          JGSOFF = 0
        end do
      end if
c     stop 'After 1st exchange type diagonal 2-particle integrals.'
*
* ( I J_ ! J_ I)
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ_|J_I):'
      end if
*
*( obtained from list 3 )
      if ((ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0).or.
     &    (ENVIRO(1:6).eq.'MOLCAS'.and.IDENSI.ge.1)) then
        if (SIGDEN_ROUTE.eq.'TRADIT') then
          DO ISM =  NSMOB/2+1, NSMOB
            DO JSM  = NSMOB/2+1, NSMOB
*. Pointer to first integral in list ( ISM JSM ! JSM ISM)
              IOFF_DBG = IOFF_DBG_LIST(3,ISM,JSM,JSM,ISM)
              NI = NTOOBS(ISM)
              NJ = NTOOBS(JSM)
              IOFF = IBSO(ISM)
              JOFF = IBSO(JSM)
              if (NTEST.ge.50) then
                write(6,*) 'I symmetry is ',ISM
                write(6,*) 'J symmetry is ',JSM
                write(6,*) 'LIST offset is ',IOFF_DBG
                write(6,*) 'I Orbitals in this symmetry (all GAS) ',NI
                write(6,*) 'J Orbitals in this symmetry (all GAS) ',NJ
                write(6,*) 'I offset for this symmetry ',IOFF
                write(6,*) 'J offset for this symmetry ',JOFF
              end if
              DO I = 1, NI
               DO J = 1, NJ
*. address of (IJ_!J_I)
                 IJ = IOFF_DBG-1 + ((I-1)*NJ+(J-1))*NJ*NI
     &                           +  (J-1)*NI+I
                 UBBUE(I-1+IOFF,J-1+JOFF) =
     &           WORK(KINT2D(1,3)-1+IJ)
                 IADRIN = KINT2D(1,3)-1+IJ
                 if (NTEST.ge.100) then
                   write(6,*) 'selected address:  ',IADRIN
                   write(6,*) '  selected integral: ',WORK(IADRIN)
                 end if
               END DO
              END DO
*             ^ End of loop over orbitals
            END DO
          END DO
*         ^ end loop over symmetries
*
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
          IKFAC = 1
          IKPERM(1) = 1
          IKPERM(2) = 2
          IKPERM(3) = 3
          IKPERM(4) = 4
          if (ENVIRO(1:6).eq.'DIRAC '.and.ISPINFREE.eq.0) then
            IKFAC = -1
            IKPERM(1) = 1
            IKPERM(2) = 2
            IKPERM(3) = 3
            IKPERM(4) = 4
          end if  
          IGSOFF = 0
          JGSOFF = 0
          ISMOFF = 0
          JSMOFF = 0
          do IGS = 1,NGAS,1
            if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
            do JGS = 1,NGAS,1
              if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT2(JGS-1)
              do ISM = NSMOB/2+1, NSMOB
                if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
                do JSM = NSMOB/2+1, NSMOB
                  if (JSM.gt.NSMOB/2+1) 
     &                       JSMOFF = JSMOFF + NGSOB2(JSM-1,JGS)
                  NI = NGSOB(ISM,IGS)
                  NJ = NGSOB2(JSM,JGS)
                  IOFF = IGSOFF + ISMOFF
                  JOFF = JGSOFF + JSMOFF
                  if (NTEST.ge.10) then
                    write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                    write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                    write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                  end if
                  do I = 1+IOFF,NI+IOFF,1
                    do J = 1+JOFF,NJ+JOFF,1
                      call geth2dia_sigden(UBBUE(I,J),
     &                                   WORK(KT_CC),T_BUFF,
     &                                   I,J,
     &                                   'EXC2',NINT12,
     &                                   I2DIATP,IKFAC,IKPERM,
     &                                   IREOTS,IREOTS2,
     &                                   WORK(KMDIASV),              
     &                                   N1DIAINT,NDIAINT,NTEST)
                    end do
                  end do
                end do
                JSMOFF = 0
              end do
              ISMOFF = 0
            end do
            JGSOFF = 0
          end do
        end if
      end if
*     ^ end if non-vanishing integral lists
*
csk      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Two-electron integrals for CI diagonal'
        WRITE(6,*)
        WRITE(6,*) ' (I I ! J J) matrix '
        WRITE(6,*)
        CALL WRTMAT(UUUUC,NTOOB,NTOOB,NTOOB,NTOOB)
        WRITE(6,*)
        WRITE(6,*) ' (I J ! J I) matrix '
        WRITE(6,*)
        CALL WRTMAT(UUUUE,NTOOB,NTOOB,NTOOB,NTOOB)
        if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
          WRITE(6,*)
          WRITE(6,*) ' (I J_! J_I) matrix '
          WRITE(6,*)
          CALL WRTMAT(UBBUE,NTOOB,NTOOB,NTOOB,NTOOB)
        end if
      END IF
csk      NTEST = 000
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'GTDIA2')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GETH1DIA4(H1DIA_U,H1DIA_B,SPININTS,INTSPACE,T_BUFF,
     &                     ENVIRO,SIGDEN_ROUTE)
*
* Fetch unbarred and barred diagonal one-electron integrals
*
*-----------------------------------
*   Timo Fleig
*    SIGDEN route incorporated, August 2001
*-----------------------------------
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "integrals_off.inc"
*
      character*6 SIGDEN_ROUTE
*
*     Input list of integrals (or in SIGDEN on WORK(KT_CC))
*     Input list of integrals (or in SIGDEN in "shared memory"
*     mode on T_BUFF (for "node master")  
      DIMENSION SPININTS(INTSPACE)
*. Output
      DIMENSION H1DIA_U(*),H1DIA_B(*),T_BUFF(*)
*. Just locally
      dimension IOBPTS_PIC(MXPNGAS+MXPR4T,MXPOBS)
      dimension IOBPTS2_PIC(MXPNGAS+MXPR4T,MXPOBS)
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDIA)
*
      ZERO = 0.D0
      if (NTEST.ge.50) then
         write(LUWRT,*) '----------------------------------------'
         write(LUWRT,*) ' GETH1DIA4 speaking '
         write(LUWRT,*) '----------------------------------------'
      end if
*
      IDUM = 1
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'GT1DIA')
      call memmar(KLH1,NTOOB**2,'ADDL  ',2,'H1BL  ')
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
*
*  Check full integral list
        if (NTEST.ge.1000) then
          write(6,*)
          write(6,*) 'Checking full integral list in geth1dia4 '
          write(6,*) '---------------------------------------- '
          ITRALINE = INTSPACE / 32
          do LIST=1,32,1
            write(6,'(A18,I2)') 'Integral list no. ',LIST
            IOFF = 0
            do IL=1,LIST-1,1
              IOFF = IOFF + ITRALINE
            end do
            do INT=1,ITRALINE,1
              write(6,'(I7,1F18.14)') INT,SPININTS(IOFF+INT)
            end do
          end do
        end if
*
*  Do some local reordering. Nothing else is affected by this.
        ICOMP = 1
        call reo1fdia(IOBPTS_PIC,IOSPIR,NGSOB,NSMOB,NGAS,ICOMP)
        ICOMP = 2
        call reo1fdia(IOBPTS2_PIC,IOSPIR,NGSOB2,NSMOB,NGAS,ICOMP)
*
*. unbarred-unbarred diagonal integrals
*
        DO IGAS = 1, NGAS
          DO ISMOB = NSMOB/2+1, NSMOB
*. Fetch one-electron integrals,
            call geth1rel3(SPININTS,WORK(KLH1),
     &           ISMOB,IGAS,ISMOB,IGAS,1,1,0,1,IRETURN,SCALIN,
     &           INTSPACE,ENVIRO)
C           subroutine geth1rel3(SPININTS,
C    &                           H,ISM,ITYP,JSM,JTYP,IRISIG,IRICO,IPM1,
C    &                           IJUB,IRETURN,SCALIN)
            IOFF = IOBPTS_PIC(IGAS,ISMOB)
            LORB = NOBPTS(IGAS,ISMOB)
            if (NTEST.ge.100) then
               write(6,*) 'IGAS is  ',IGAS
               write(6,*) 'Symm is  ',ISMOB
               write(6,*) 'IOFF is  ',IOFF
            end if
            CALL COPDIA(WORK(KLH1),H1DIA_U(IOFF),LORB,0)
C                COPDIA(A,VEC,NDIM,IPACK)
          END DO
        END DO
*
*. barred-barred diagonal integrals
*
        DO IGAS = 1, NGAS
          DO ISMOB = NSMOB/2+1, NSMOB
*. Fetch one-electron integrals,
            call geth1rel3(SPININTS,WORK(KLH1),
     &           ISMOB,IGAS,ISMOB,IGAS,1,1,0,2,IRETURN,SCALIN,
     &           INTSPACE,ENVIRO)
            IOFF = IOBPTS2_PIC(IGAS,ISMOB)
            LORB = NOBPTS2(IGAS,ISMOB)
            if (NTEST.ge.100) then
               write(6,*) 'IGAS is  ',IGAS
               write(6,*) 'Symm is  ',ISMOB
               write(6,*) 'IOFF is  ',IOFF
            end if
            CALL COPDIA(WORK(KLH1),H1DIA_B(IOFF),LORB,0)
          END DO
        END DO
*
      else if (SIGDEN_ROUTE.eq.'SIGDEN') then
        IOFFSD = 1
C
C       GAS-SYMM ordering (GAS is outer loop!)
C
        do IGAS=1,NGAS,1
          do ISMOB=1,NSMOB/2,1
* Determine index range for this symmetry/GAS block:
            NORBSM = 0
            do IGASS=1,IGAS-1,1
              do ISM=1,NSMOB/2,1
                NORBSM = NORBSM + NGSSH(ISM,IGASS)
              end do
            end do
            do JSM=1,ISMOB-1,1
              NORBSM = NORBSM + NGSSH(JSM,IGAS)
            end do
            ISTART = NORBSM + 1
            IEND = NORBSM + NGSSH(ISMOB,IGAS)
            if (NTEST.ge.50) then
              write(LUWRT,'(2X,A,5I4)') 'IGAS,ISMOB,NORBSM,ISTART,IEND',
     &                               IGAS,ISMOB,NORBSM,ISTART,IEND
            end if
* Check for occupation of current space
            if (MXNELGS(IGAS).ne.0) then
              IUB = 1
              call geth1dia_sigden(IUB,WORK(KLH1),WORK(KT_CC),T_BUFF,
     &                             WORK(KMDIASV),              
     &                             N1ELINT,N1DIAINT,NGSSH(ISMOB,IGAS),
     &                             ISTART,IEND,IREOTS,NTEST)
              call dcopy(NGSSH(ISMOB,IGAS),WORK(KLH1),1,
     &                   H1DIA_U(IOFFSD),1)
              IUB = 2
              call geth1dia_sigden(IUB,WORK(KLH1),WORK(KT_CC),T_BUFF,
     &                             WORK(KMDIASV),              
     &                             N1ELINT,N1DIAINT,NGSSH(ISMOB,IGAS),
     &                             ISTART,IEND,IREOTS2,NTEST)
              call dcopy(NGSSH(ISMOB,IGAS),WORK(KLH1),1,
     &                   H1DIA_B(IOFFSD),1)
            else
              call dzero(H1DIA_U(IOFFSD),NGSSH(ISMOB,IGAS))
              call dzero(H1DIA_B(IOFFSD),NGSSH(ISMOB,IGAS))
            end if
            IOFFSD = IOFFSD + NGSSH(ISMOB,IGAS)
          end do
        end do
      end if
*
*
csk      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Diagonal integrals, unbarred '
        WRITE(LUWRT,*)
        CALL WRTMATMN(H1DIA_U,1,NTOOB,1,NTOOB,LUWRT)
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Diagonal integrals, barred '
        WRITE(LUWRT,*)
        CALL WRTMATMN(H1DIA_B,1,NTOOB,1,NTOOB,LUWRT)
      END IF
*
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'GT1DIA')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Get diagonal one-electron integrals in SIGDEN formalism.
*
      subroutine geth1dia_sigden(IUB,XOUT,XINT,T_BUFF,
     &                           MDIASV,N1E,N1DIAINT,NUMINT,
     &                           ISTART,IEND,IREOTS,IPRNT)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "krmc_shmem.h"
*
      DIMENSION IREOTS(*),MDIASV(4,*)
      DIMENSION XINT(N1E),XOUT(*),T_BUFF(*)
      REAL*8 SXINT
      COMMON/H1DIA_INT/SXINT
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
csk      NTEST = 100
*
      IZERO = 0
      ZERO = 0.D0
      ICT = 1
      do ISEA=ISTART,IEND,1
*
* Rename to SYMM-TYPE ordering on external list:
        ISEA_REO = IREOTS(ISEA)
C       do INT=1,N1E,1
        do INT=1,N1DIAINT,1
C         if (INDEX(IUB,INT).eq.ISEA_REO.and.INDEX(IUB+2,INT).eq.
C    &        ISEA_REO.and.INDEX(IUB+1,INT).eq.IZERO) then
          if (MDIASV(IUB,INT).eq.ISEA_REO.and.MDIASV(IUB+2,INT).eq.
     &        ISEA_REO.and.MDIASV(IUB+1,INT).eq.IZERO) then
C
            if (NTEST.ge.100) then
              write(LUWRT,'(A,I5,A,I5,A,1F18.12)') 
     &                'Writing integral no. ',IPOSDIA(INT),'to ',ICT,
     &                '  Integral ',XINT(IPOSDIA(INT))
            end if
C           XOUT(ICT) = XINT(INT)
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
C           write(LUWRT,*) ' fetch integral from pos',IPOSDIA(INT)
            XOUT(ICT) = T_BUFF(IPOSDIA(INT))
#else
            if( lowsrt_ijkl )then
              XOUT(ICT) = T_BUFF(IPOSDIA(INT))
            else
              XOUT(ICT) = XINT(IPOSDIA(INT))
            end if
#endif
            ICT = ICT + 1
            goto 10
          end if
        end do
        write(LUWRT,*) 'No integral found in geth1dia_sigden.'
        write(LUWRT,*) 'Index pair ',ISEA_REO,ISEA_REO
        CALL QUIT( '*** Error in geth1dia_sigden.No integral found ***')
  10    continue
      end do
      if ((ICT-1).ne.NUMINT) then
        write(LUWRT,*) 'Number of processed integrals: ',ICT-1
        write(LUWRT,*) 'Number of diagonal elements  : ',NUMINT
        write(LUWRT,*) 'Problem in geth1dia_sigden.'
        CALL Abend2( ' Quitting in geth1dia_sigden.')
      end if
*
csk      NTEST = 000
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *

*                                                                     *
***********************************************************************
      subroutine geth2dia_sigden(XOUT,XINT,T_BUFF,ISEA,JSEA,TYPE,
     &                           N12,I2DIATP,IKFAC,IKPERM,
     &                           IREOTS,IREOTS2,MDIASV,N1DIAINT,
     &                           NDIAINT,IPRNT)
*
*     Get diagonal two-electron integrals in SIGDEN formalism.
*     Timo Fleig, August 2001
*
*     Modified to account for integrals over true spinors, i.e.
*     explicit reference to unbarred (+) and barred (-) indices (IKFAC).
*
*     Revised routine - Timo Fleig, February 2002
*
*     "shared memory" MPI implementation by Stefan Knecht - Nov. 2007
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "krmc_shmem.h"
*
      character*4 TYPE
*
      dimension IAC(2),JAC(2),         IKPERM(4),
     &          IREOTS(*),IREOTS2(*),MDIASV(4,*)           
      dimension XINT(NINT12), T_BUFF(*)
      REAL*8 SXINT2
      COMMON/H2DIA_INT/SXINT2
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
csk      NTEST  = 100
*
      IERROR = 1
      N2B = IBTSOSO_TP(I2DIATP)
C     do INT = N2B,N12,1
      do INT = N1DIAINT+1,NDIAINT,1
        if (TYPE.eq.'COUL') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(2),INT)
C         K = INDEX(IKPERM(3),INT) * IKFAC
C         L = INDEX(IKPERM(4),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(2),INT)
          K = MDIASV(IKPERM(3),INT) * IKFAC
          L = MDIASV(IKPERM(4),INT) * IKFAC
        else if (TYPE.eq.'EXC1') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(4),INT) * IKFAC
C         K = INDEX(IKPERM(2),INT)
C         L = INDEX(IKPERM(3),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT) * IKFAC
          K = MDIASV(IKPERM(2),INT)
          L = MDIASV(IKPERM(3),INT) * IKFAC
        else if (TYPE.eq.'EXC2') then
C         I = INDEX(IKPERM(1),INT)
C         J = INDEX(IKPERM(4),INT)
C         K = INDEX(IKPERM(2),INT) * IKFAC
C         L = INDEX(IKPERM(3),INT) * IKFAC
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT)
          K = MDIASV(IKPERM(2),INT) * IKFAC
          L = MDIASV(IKPERM(3),INT) * IKFAC
        end if
*
* Rename: External ST ordering, internal TS ordering:
        if (I.eq.IREOTS(ISEA).and.
     &      J.eq.IREOTS(ISEA).and.
     &      K.eq.IREOTS2(JSEA).and.
     &      L.eq.IREOTS2(JSEA)) then
C
          if (NTEST.ge.100) then
C           write(6,'(A,X,A,I5,3X,4I3,1F18.10)') 
C    &            'Diagonal int type ',TYPE,
C    &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            write(LUWRT,'(A,1X,A,1I4,3X,4I3,1F18.10)') 
     &            'Diagonal int type ',TYPE,
     &             IPOSDIA(INT),(MDIASV(II,INT),II=1,4),
     &             XINT(IPOSDIA(INT))
          end if
C         XOUT = XINT(INT)
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
C         write(6,*) ' fetch integral (2e) from pos',IPOSDIA(INT)
          XOUT = T_BUFF(IPOSDIA(INT))
#else
          if( lowsrt_ijkl )then 
            XOUT = T_BUFF(IPOSDIA(INT))
          else
            XOUT = XINT(IPOSDIA(INT))
          end if
#endif
          IERROR = 0
* and exit.
          goto 10
        end if
      end do
      if (IERROR.eq.1) then
        write(6,*) 'No integral found in geth2dia_sigden.'
        write(6,*) 'Type of integral = ',TYPE
        write(6,*) 'Index pair ',IREOTS(ISEA),IREOTS2(JSEA)
        write(6,*) 'Indices ',I,J,K,L
        stop 'Quitting.'
      end if
*
  10  continue
csk      NTEST = 0
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* Find range of active orbital indices for a given symmetry
* INSYM from GAS/SYM ordering.
*   Timo Fleig, August 2001
* Indices of unoccupied spaces are set to zero, the corresponding
* integrals as well.
*   Timo Fleig, October 2001
*
      subroutine orbind_gssm(INSYM,NSM,INDGATH,NELM,IPRNT)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      dimension INDGATH(NELM)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      IGSOFF = 0
      do IGAS = 1,NGAS,1
        NORBSM = 0
        do JGAS = 1,IGAS-1,1
          do ISM = 1,NSM,1
            NORBSM = NORBSM + NGSSH(ISM,JGAS)
          end do
        end do
        do ISM = 1,INSYM-1,1
          NORBSM = NORBSM + NGSSH(ISM,IGAS)
        end do
        ISTA = NORBSM + 1
        IEND = NORBSM + NGSSH(INSYM,IGAS)
* write to array
        do I = 1,NGSSH(INSYM,IGAS),1
          if (MXNELGS(IGAS).ne.0) then
            INDGATH(I+IGSOFF) = ISTA + I - 1
          else
            INDGATH(I+IGSOFF) = IZERO
          end if
        end do
        IGSOFF = IGSOFF + NGSSH(INSYM,IGAS)
      end do
*
      if (NTEST.ge.5) then
        write(6,*) 'Output from orbind_gssm:'
        write(6,*) 'Symmetry ',INSYM
        write(6,*) 'INDGATH array :'
        call iwrtma(INDGATH,1,NELM,1,NELM)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*  Create new arrays IOBPTS and IOBPTS2 for reordering diagonal
*  one-electron integrals (symmetry/type ordering) to internal
*  LUCIA ordering type/symmetry.
*
      subroutine reo1fdia(IOBPTS,IORD,NGSOB,NSMOB,NGAS,ICOMP)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*
      dimension IOBPTS(MXPNGAS,*),IORD(MXPOBS,*)
      dimension NGSOB(MXPOBS,MXPNGAS)
*
      ITSOFF = 1
      DO IISYM = 1, NSMOB
        DO IGAS = 1, NGAS
          ISYM = IORD(ICOMP,IISYM)
          IOBPTS(IGAS,ISYM) = ITSOFF
          ITSOFF = ITSOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
      return
      end
