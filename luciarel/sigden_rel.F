      SUBROUTINE SIGDEN_REL(C,HC,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                      LLUC,LLUHC,T,ISIGDEN)
*
* Outer routine for SIGMA and Density calculation for 
* general relativistic code
*
* For Sigma T is the input set of coefficients
* For Densi T is the output set of coefficients
*
* Jeppe Olsen, Summer of 99, March 2001
* 
* Modifications for vector storage and sigma loops
* Timo Fleig, Halloween 2001
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
#include "cands.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#include "clunit.inc"
*
* Batches of sigma
      integer LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      dimension HC(*),C(*)
*
      CALL QENTER('SIGDE')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDEN')
*
      NTESTL = 0000
      NTEST = max(NTESTL,max(IPRSIG,IPRDEN))
*
      IF(NTEST.GE.100) THEN 
        write(6,*)
        write(6,*)
        write(6,*) '++++++++++++++++++++++++++++++++++++++'
        write(6,*) '++++++++++++++++++++++++++++++++++++++'
        WRITE(6,*) '       SIGDEN_REL in action           '
        write(6,*) '++++++++++++++++++++++++++++++++++++++'
        write(6,*) '++++++++++++++++++++++++++++++++++++++'
        WRITE(6,*) '          ICSM, ISSM = ', ICSM, ISSM
      END IF
*
*=====================================================
* Implementing time-reversal symmetry 
* in the even-particle case, according to
*   Cm_ = Cm*
*   Cm  = Cm_*
* No signs here, phase factors always 1. (In odd-particle
* case -1 for Nbeta odd, Nalpha even!)
*
*  Timo Fleig, 31.07.2001
*
* define symmetry, factors, and loop restrictions
      call evenodd(IEONELEC,NACTEL)
      if (IEONELEC.eq.1) then
        IKRA_SYM = 0
      else if (IEONELEC.eq.2) then
        IKRA_SYM = 1
      end if
      IKRA_SYM = 0         ! FOR THE MOMENT !
* 
*=====================================================
*
      LBLK = -1
*
      if (NTEST.ge.1000) then
        write(6,*) '--------------------------------------'
        write(6,*) ' Printing vectors in SIGDEN_REL.      '
        write(6,*) '--------------------------------------'
        write(6,*) '  Right vector: '
        call wrtvcdc(C,LLUC,1,-1)
        if (ISIGDEN.eq.2) then
          write(6,*) '  Left vector: '
          call wrtvcdc(HC,LLUHC,1,-1)
        end if
      end if
*
      CALL REWINO(LLUC)
      CALL REWINO(LLUHC)
* FIXME: Only real density working currently.
*        Problem with sigma vector reading of imaginary part.
*        ???? Still a problem ???? Check !!!!
*        Probably not.
      IRILP = 2
* LUCIAREL is a relativistic CI program. Spinfree formalism not
*          supported, only for testing purposes. Therefore:
*          ALWAYS real AND imaginary!
*
*. Loop over vectors of sigma and C : MK2 space + Real/Imag
      DO IRIS = 1,IRILP,1
        if (NTEST.ge.50) then
          write(6,*)
          write(6,*) '-----------------------'
          write(6,*) ' SIGMA LOOP , IRIS = ',IRIS
          write(6,*) '-----------------------'
        end if
        JBATABS = 0
        DO IMK2_S = 1, NMS2VAL
          MK2_S = MS2VAL(IMK2_S)
          if (NTEST.ge.50) then
            write(6,*) '- - - - - - - - - - - -'
            write(6,*) ' MK2_S LOOP , MK2_S = ',MK2_S
            write(6,*) '- - - - - - - - - - - -'
          end if
*
          IASTP = IST_FOR_DT(1,IMK2_S)
          IBSTP = IST_FOR_DT(2,IMK2_S)
          NOCTPA_S = NOCTYP(IASTP)
          NOCTPB_S = NOCTYP(IBSTP)
*. Arrays giving allowed type combinations (of alpha and beta)
          CALL MEMMAN(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
          CALL IAIBCM(ISSPC,IASTP,IBSTP,WORK(KSIOIO))
*. Arrays giving block type
          KSVST = 1
          CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
          CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP))
*
* Loop over batches of sigma for this MK2_S projection value
          do JBATS = 1,NBAT_MS2(IMK2_S)
            JBATABS = JBATABS + 1
*
            if (ISIGDEN.eq.1) then
              if (NTEST.ge.500) then
                 WRITE(6,*)
                 WRITE(6,*)
     &           ' =============================================='
                 WRITE(6,*)
     &           ' Start of construction of Sigma batch ', JBATABS
                 write(6,'(2X,A,I3)') ' for MK2_S value ',MK2_S
                 WRITE(6,*) 
     &           ' =============================================='
                 WRITE(6,*)
              end if
              LS = LEBATS(JBATABS)
              ZERO = 0.D0
              call setvec(HC,ZERO,LS)
            end if
*
            call rewino(LLUC)
            DO IRIC = 1,IRILP,1
            if (NTEST.ge.50) then
              write(6,*)
              write(6,*) '-----------------------'
              write(6,*) ' C     LOOP , IRIC = ',IRIC
              write(6,*) '-----------------------'
            end if
              IF(IRIC.EQ.2) THEN
*. Skip EOV mark between real and imaginary part
                CALL IFRMDS(IONEM,1,-1,LLUC)
              END IF
*
              DO IMK2_C = 1, NMS2VAL
                MK2_C = MS2VAL(IMK2_C)
                if (NTEST.ge.50) then
                  write(6,*) '- - - - - - - - - - - -'
                  write(6,*) ' MK2_C LOOP , MK2_C = ',MK2_C
                  write(6,*) '- - - - - - - - - - - -'
                end if
*
* Kramers symmetry restriction
                if (IKRA_SYM.eq.1.and.MK2_C.lt.0) then
                  if (IRIC.eq.1) IKSFAC = 1
                  if (IRIC.eq.2) IKSFAC = -1
* copy relevant block of sigma vector to position and 
* multiply with IKSFAC.
* FIXME !!
                  stop 'KS not implemented yet.'
                else
                  IACTP = IST_FOR_DT(1,IMK2_C)
                  IBCTP = IST_FOR_DT(2,IMK2_C)
*. Copy blocks with given MK2_C to LUSC37
                  CALL REWINO(LUSC37)
                  NBL_C = NBLK_MS2(IMK2_C)
                  CALL COPNBLKD(LLUC,LUSC37,C,NBL_C,0,LBLK)
*. End of vector mark
                  CALL ITODS(-1,1,-1,LUSC37)
*
* Analyze T coefficients before doing the job
                  call flush(6)
                  if (NTEST.ge.100) then
c                   call ana_gencc(T,ISYM_T)
c                   stop
                  end if
*
                  call sigden_rel2(C,HC,JBATABS,NBATCH_MK2,
     &                             LBATS,LEBATS,I1BATS,IBATS,
     &                             LUSC37,LLUHC,T,
     &                             ISIGDEN,WORK(KVEC3),IACTP,IBCTP,
     &                             IRIS,IRIC,IMK2_C,LUSC41,IPRDEN)
                end if
*               ^ End if Kramers restriction
              END DO
*             ^ End of loop over IMK2_C
            END DO
*           ^ End of loop over IRIC
*
            if (ISIGDEN.eq.1) then
              if (NTEST.ge.100) then
                 WRITE(6,*)
                 WRITE(6,*) ' **********************************'
                 WRITE(6,*) ' *                                *'
                 WRITE(6,*) ' * Array containing final sbatch  *'
                 WRITE(6,*) ' *                                *'
                 WRITE(6,*) ' **********************************'
                 WRITE(6,*)
                 CALL WRTMAT(HC,1,LS,1,LS)
              end if
              KBATSEND = I1BATS(JBATABS) + LBATS(JBATABS)
              do ISBLK = I1BATS(JBATABS),KBATSEND-1,1
                IOFF = IBATS(6,ISBLK)
                LEN_T = IBATS(8,ISBLK)
                call itods(LEN_T,1,-1,LLUHC)
                call todsc(HC(IOFF),LEN_T,-1,LLUHC)
              end do
            end if
*
          end do
*         ^ End loop over S batches of given MK2_S value
        END DO
*       ^ End of loop over IMK2_S
*
        if (ISIGDEN.eq.1) then
* Write end of vector mark between real and imaginary part.
          call itods(-1,1,-1,LLUHC)
        else if (ISIGDEN.eq.2) then
* Skip end of vector mark between real and imaginary part.
          call ifrmds(IONEM,1,-1,LLUHC)
        end if
*
      END DO
*     ^ End of loop over IRIS
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'SIGDEN')
*
      CALL QEXIT('SIGDE')
*
* Close scratch unit for Right vector
      close(unit=LUSC37,status='DELETE')
*
      RETURN
      END
