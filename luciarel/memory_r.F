***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ALLOC_REL
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
* /GLBBAS/,/CDIM/
*
*
* =======
* Version
* =======
*
* Modified Jan 1997
* SIGDEN_ROUTE SIGDEN added by Timo Fleig, May-June 2001

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "csm.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cintfo.inc"
#include "ctcc.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

*.Output
#include "glbbas.inc"
CTF adding global info
      common/glbref/ISPC/ISM/ITERGLB
*.1 : One electron integrals( Complete matrix allocated )
      call memmar(KINT1,NTOOB ** 2,'ADDS  ',2,'INT1   ')
*.1.1 : Inactive fock matrix
      call memmar(KFI  ,NTOOB ** 2,'ADDS  ',2,'FI     ')
*.2 : Two electron integrals
      IF(NOINT.EQ.0.AND.INCORE.EQ.1) THEN
      call memmar(KINT2,NINT2,'ADDS  ',2,'INT2   ')
      END IF
*. Pointers to symmetry block of integrals
      call memmar(KPINT1,NBINT1,'ADDS  ',2,'PINT1   ')
      call memmar(KPINT2,NBINT2,'ADDS  ',2,'PINT2   ')
*. Symmetry of last index as a function of initial index
      call memmar(KLSM1,NBINT1,'ADDS  ',2,'LSM1   ')
      call memmar(KLSM2,NBINT2,'ADDS  ',2,'LSM2   ')
*
* =============================
*. For double group integrals
* =============================
* Caution!
* Double group calculation can also be non-relativistic!
* (Well, at least technically!)
      IF(IDBG.EQ.1) THEN
        if (SIGDEN_ROUTE.eq.'TRADIT') then
*. 1e integrals
          DO IL1 = 1, NL1D
*. Real part
            call memmar(KINT1D(1,IL1),NTOOB ** 2,'ADDS  ',2,'INT1D ')
*. Imaginary part
            call memmar(KINT1D(2,IL1),NTOOB ** 2,'ADDS  ',2,'INT1D ')
*. and a pointer
            call memmar(KPINT1D(IL1),NBINT1,'ADDS  ',2,'PINT1D')
          END DO
          call memmar(KPINT1,NBINT1,'ADDS  ',2,'PINT1   ')
          call memmar(KPINT2,NBINT2,'ADDS  ',2,'PINT2   ')
*. 2e integrals
          DO IL2 = 1, NL2D
*. Real part
            LENGTH = L2D(IL2)
C?          write(6,*) 'Current LENGTH is  ',LENGTH,'IL2 is  ',IL2
            call memmar(KINT2D(1,IL2),LENGTH,'ADDS  ',2,'INT2D ')
*. Imaginary part
            call memmar(KINT2D(2,IL2),LENGTH,'ADDS  ',2,'INT2D ')
*. and a pointer
            call memmar(KPINT2D(IL2),NBINT2,'ADDS  ',2,'PINT2D')
          END DO
*
*  Allocation of density matrices
*  See readin.f for explanations
          if (IDENSI.eq.1.or.IDENSI.eq.2) then
            call memmar(KRHO1,NTOOB ** 2,'ADDS  ',2,'RHO1  ')
            if (IDENSI.eq.2) then
              LRHO2 = NTOOB**2*(NTOOB**2+1)/2
              call memmar(KRHO2,LRHO2     ,'ADDS  ',2,'RHO2  ')
              write(6,*) 'Message from alloc:'
              write(6,*) 'These density routines are not operating yet.'
              write(6,*) 'Change DENSI in input or ask the experts.'
              Call Abend2( 'LUCIAREL is quitting.' )
            end if
          else if (IDENSI.eq.3.or.IDENSI.eq.4) then
* We will calculate UU and UB, real and imaginary parts
            LRHO1 = NTOOB ** 2
            do IROOTCT = 1,NROOT,1
              do IRI = 1,2,1
                do JBU = 1,2,1
                  do IBU = 1,2,1
* Redundant density blocks!
c                   IBU = 2
                    call memmar(KRHO1_DBG(IBU,JBU,IRI,IROOTCT),
     &                          LRHO1,'ADDS  ',2,'RHO_DB')
                  end do
                end do
              end do
            end do
            if (IDENSI.eq.4) then
* Still to be allocated !!
              write(6,*) 'Spin-dependent two-particle densities '
              write(6,*) 'not yet debugged. '
              Call Abend2( 'I would better quit.' )
            end if
          end if
*
* We set the pointer for the T array here as well. There is no
* explicit access to the array in this case, so no memory faults
* should occur.
          call memmar(KT_CC,1,'ADDS  ',2,'KT_CC ')
        else if (SIGDEN_ROUTE.eq.'SIGDEN') then
*
* Allocate the amplitude array T only with general sigden
* structure for integrals or densities.
          NINT_CC = 0
          do IL1=1,NL1D,1
            NINT_CC = NINT_CC + 2 * NTOOB**2
          end do
*
          do IL2=1,NL2D,1
            NINT_CC = NINT_CC + 2 * L2D(IL2)
          end do
*
          call memmar(KT_CC,NINT_CC,'ADDS  ',2,'KT_CC ')
C FIXME: Add symmetry savings
          write(6,*) 'NOTE:'
          write(6,*) 'Allocation should be modified to account '
          write(6,*) 'for possible symmetry savings.'
        end if
      end if
*
*.4. Integrals (ij!kk) -(ik!kj), pointer and space
      call memmar(KPNIJ,NTOOB ** 2,'ADDS  ',1,'KPNIJ ')
      call memmar(KIJKK,NTOOB**2 *(NTOOB+1) / 2, 'ADDS  '
     &            ,2,'KIJKK ')
CTF  Allocate memory for initial guess for CI vector
      IF(MXP1+MXP2+MXQ. NE. 0 .OR. IPROCC.NE.0 ) THEN
*. Allocate memory for explicit hamiltonian and roots
        NSBDET = MXP1 + MXP2 + MXQ
        NSBDETP = MAX(NSBDET,IPROCC)
        MXP = MXP1 + MXP2
        call memmar(KSBEVC,NSBDETP*NROOT,'ADDS  ',2,'KSBEVC')
        call memmar(KSBEVL,       NROOT,'ADDS  ',2,'KSBEVL')
*. KSBIDT must be able to hold list of subspace dets / print dets
        LSCR = NSBDETP+1
        call memmar(KSBIDT,LSCR,'ADDS  ',1,'KSBIDT')
*. ( NSBDET + 1 elements allocated since NSBDET+1 values are
* obtained in order to check for degenerencies )
CADDB
        call memmar(KSBCNF,LSCR        ,'ADDS  ',1,'KSBCNF')
CADDE
        call memmar(KSBIA ,LSCR        ,'ADDS  ',1,'KSBIA ')
        call memmar(KSBIB ,LSCR        ,'ADDS  ',1,'KSBIB ')
*
        LH0 = MXP*(MXP+1)/2 + MXP1*MXQ + MXQ
        call memmar(KH0,LH0,'ADDS  ',2,'KH0   ')
*. Scratch space for manipulating with H0
        LH0SCR =
     &  MAX(6*NSBDETP,4*NSBDETP+4*NOCOB,
     &      MXP1*(MXP1+1)/2+MXP1**2)
        call memmar(KH0SCR,LH0SCR      ,'ADDS  ',2,'KH0SCR')
      ELSE
        KSBEVC = 1
        KSBEVL = 1
        KSBIDT = 1
        KSBCNF = 1
        KH0    = 1
        KH0SCR = 1
      END IF
*. indices for pair of orbitals symmetry ordered
*
      call memmar(KINH1,NTOOB*NTOOB,'ADDS  ',2,'KINH1  ')
*
      IF (NOMOFL.EQ.0) THEN
        write(6,*)
        WRITE(6,*)
     &  ' Alloc : General allocation for SPAO trafo '
        WRITE(6,*)
     &  '         still to be done !! (now locally)'
*. MO-AO transformation
CTF     CALL MEMMAR(KMOAOIN,NTOOB ** 2,'ADDS  ',2,'MOAOIN')
CTF     CALL MEMMAR(KMOAOUT,NTOOB ** 2,'ADDS  ',2,'MOAOUT')
*. Spinor-spin atomic orbital transformation
CTF     CALL MEMMAR(KSPSAOIN,4*NTOOB**2,'ADDS  ',2,'SPOAIN')
CTF     CALL MEMMAR(KSPSAOUT,4*NTOOB**2,'ADDS  ',2,'SPOAUT')
*. and their offsets
CTF     CALL MEMMAR(KISPSAOIN,4,'ADDS  ',1,'ISPOAI')
CTF     CALL MEMMAR(KISPSAOUT,4,'ADDS  ',1,'ISPOAU')
*
      ELSE
        KMOAOIN = 0
        KMOAOUT = 0
        KSPSAOIN = 0
        KSPSAOUT = 0
        KISPSAOIN = 0
        KISPSAOUT = 0
      END IF
*. (more to come )
      call memmar(KDUMMY,IDUMMY,'CHECK',IDUMMY,'Dummy ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE RMEMCHK
#include "ipoist8.inc"
*
* Check memory allocated  with the memory manager
*
      call memmar(KDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MEMMAR(KBASE,KADD,TASK,IR,IDENT)
*
* Memory manager routine for LUCIAREL "R"
*
* KBASE : New base address
*         If TASK = INI, KBASE is offset for memory to be controlled
*         by MEMMAR
* KADD  : Dimension of array to be added
*         If TASK = INI, KADD is total length of array
* TASK  : = INI  : Initialize                 Character*6
*         = ADDS : Add static memory
*         = ADDL : Add Local memory
*         = FLUSH : Flush local memory
*         = CHECK : Check memory paddings
*         = FREE  : Return first Free word in KBASE
*         = MARK  : Set a mark at current free adress
*         = FLUSM : Flush local memory to previous mark
*         = PRINT : Print memory allocations
* IR    : 1 => integer , 2 => real,
*         ratio between integer and real is IRAT
* IDENT : identifier of memory slice,Character*6
*
* Local Memory not flushed before allocation of additional static memory
* is tranferred to static memory
*
      IMPLICIT REAL*8(A-H,O-Z)
#if defined (IA64)
      integer*8 KBASE,KFREES,KFREEL,IBASE,NWORD,JBASE,JBASEN,IMARK
#endif
      CHARACTER*6 TASK,IDENT,IIDENT,MARKC
      PARAMETER(NPAD = 1 )
      PARAMETER(MAXLVL = 4096 )
      PARAMETER(MAXMRK = 1024)
      COMMON/CMEMOR/NWORD,KFREES,KFREEL,IBASE(MAXLVL),IMARK(MAXMRK),
     &              LENGTH(MAXLVL),IIDENT(MAXLVL),NS,NL,NM,
     &              MARKL(MAXMRK),MARKS(MAXMRK),MARKC(MAXMRK)
*. Two real*8 words, one added NPAD times before each array, another
*. added NPAD times after each array
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
*
#include "mxpdim.inc"
#include "irat.inc"
#include "wrkspc.inc"
*
      ILPRNT = 0
      IF(ILPRNT .NE. 0 ) THEN
        WRITE(6,*) ' KBASE KADD, IR ', KBASE,KADD,IR
        WRITE(6,'(3A)') ' TASK , IDENT ', TASK,IDENT
        WRITE(6,*) '  NS NL NM ',NS,NL,NM
      END IF
*
      ISTOP = 0
      ITSOK  = 1
      IPRNTMP = 0
*
      IF(TASK(1:5).EQ.'PRINT') THEN
        IPRNTMP = 1
        NSNLI= NS + NL
*
      ELSE IF(TASK(1:3).EQ.'INI') THEN
*
**.Initialize
*
        NS = 0
        NL = 0
        NSNLI = 0
        NM = 0
        KFREES = KBASE
        KFREEL = KBASE
*       NWORD = MXPWRD
*
*       Compute pointer to end element
*
        NWORD = KBASE + MXPWRD-1
        IPRNTMP = 0
        ISTOP = 0
*
** First free word
*
      ELSE IF (TASK(1:4).EQ.'FREE') THEN
       KBASE = KFREEL
       NSNLI = NS+NL
*
**. Static memory
*
      ELSE IF(TASK(1:4).EQ.'ADDS') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1)/IRAT + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(6,*)
          WRITE(6,*) ' You can''t always get what you want'
          WRITE(6,*) ' No, you can''t always get what you want'
          WRITE(6,*) ' But if you try sometime, you may find '
          WRITE(6,*) ' you get what you need '
          WRITE(6,*) '                       Jagger/Richard '
*
          WRITE(6,*) ' MEMMAR : work array too short '
          WRITE(6,*) ' current and required length ',NWORD,KFREES-1
*
          WRITE(6,*) ' Trying to allocate : identifer,length'
          WRITE(6,'(20X,A,I9)')IDENT,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(6,*) ' Too many levels in MEMMAR '
          write(6,*) ' NS has value  ',NS
          WRITE(6,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1( 11 )
        END IF
        IIDENT(NS) = IDENT
        KFREEL = KFREES
        IBASE(NS) = KBASE
        DO 10 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   10   CONTINUE
*
**. Local memory
*
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1)/IRAT + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(6,*)
          WRITE(6,*) ' You can''t always get what you want'
          WRITE(6,*) ' No, you can''t always get what you want'
          WRITE(6,*) ' But if you try sometime, you may find '
          WRITE(6,*) ' you get what you need '
          WRITE(6,*) '                       Jagger/Richard '

          WRITE(6,*) ' MEMMAR : work array too short '
          WRITE(6,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(6,*) ' Trying to allocate : identifer,offset,  length'
          WRITE(6,'(24X,A,2I9)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(6,*) ' Too many levels in MEMMAR '
          WRITE(6,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1( 11 )
        END IF
        IIDENT(NS+NL) = IDENT
        IBASE(NS+NL) = KBASE
        DO 20 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   20   CONTINUE
*
** Flush local memory
*
      ELSE IF(TASK(1:5).EQ.'FLUSH') THEN
        NSNLI = NS+NL
        KFREEL = KFREES
        NL = 0
*. Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF(TASK(1:4).EQ.'MARK') THEN
*. Set a mark at current free address
        NM = NM + 1
        IF(NM.GT.MAXMRK) THEN
          WRITE(6,*) ' Too many marks  in MEMMAR '
          WRITE(6,*) ' Increase MAXMRK from ', MAXMRK
          Call Abend1( 11 )
        END IF
        MARKC(NM) = IDENT
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
        NSNLI = NS + NL
      ELSE IF (TASK(1:5).EQ.'FLUSM') THEN
        NSNLI = NS+NL
        KFREEL = IMARK(NM)
        IF(KFREES.GT.IMARK(NM)) KFREES = IMARK(NM)
        NL = MARKL(NM)
        NS = MARKS(NM)
        NM = NM - 1
*. Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF( TASK(1:5).EQ.'CHECK') THEN
        NSNLI = NS+ NL
      ELSE
          WRITE(6,*) ' MEMMAR : Unknown task parameter ',TASK
          WRITE(6,*) ' Too confused to continue  '
          Call Abend1( 11 )
      END IF
*
**. Check paddings
*
      ICHECK = 1
      ITSOK = 1
      IF(TASK(1:5).EQ.'CHECK'.OR.ICHECK.EQ.1) THEN
        DO 100 IL = 1, NSNLI
          JBASE = IBASE(IL)
          IF(IL.NE.NSNLI) THEN
           JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF
          L1OK = 1
          L2OK = 1
          DO IPAD = 1, NPAD
            IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1 .OR.
     &         WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) THEN
               ITSOK = 0
               WRITE(6,*) ' Memory problem for : '
               WRITE(6,*) '   Level (IL) ',IL
               IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
               IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
               IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
                WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
               ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
                WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY       -  '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &            IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
                WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'     -       -    '
               END IF
            END IF
          END DO
  100   CONTINUE
      END IF
 1001 CONTINUE
*
        IF(ITSOK.EQ.0.OR.IPRNTMP.NE.0) THEN
          WRITE(6,'(A,A)') ' Current task : ', TASK
          WRITE(6,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          WRITE(6,'(A,A)') ' Current Identifier : ', IDENT
          IF (ITSOK.EQ.0) THEN
            WRITE(6,*) '  Sorry to say it , but memory is CORRUPTED '
	  END IF
          WRITE(6,*) '  Memory map : '
          WRITE(6,*) '  Identifier   Offset    Pad1 okay Pad2 okay '
          WRITE(6,*) '  ========== ==========  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            DO 40 IPAD = 1, NPAD
              IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
   40       CONTINUE
            DO 50 IPAD = 1, NPAD
              IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
   50       CONTINUE
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(6,'(4X,A6,2X,I20,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       -    '
            END IF
  200     CONTINUE
*
* Marks
*
          WRITE(6,*)
          WRITE(6,*) '======='
          WRITE(6,*) ' Marks '
          WRITE(6,*) '======='
          WRITE(6,*)
*
          WRITE(6,*) ' Identifier  Start of free memory '
          WRITE(6,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(6,'(3X,A6,10X,I20)') MARKC(JMARK),IMARK(JMARK)
          END DO
*
        IF(ITSOK.EQ.0) 
     &     Call Abend2( ' Error observed by  memory manager ' )
        END IF
*
      IF(ISTOP.NE.0) call Abend2
     &               ( ' Error observed by  memory manager ' )
      RETURN
      END
