      SUBROUTINE GETH1DIA4(H1DIA_U,H1DIA_B,SPININTS,INTSPACE,
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
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
      character*6 SIGDEN_ROUTE
*
*. Input list of integrals (or in SIGDEN on WORK(KT_CC))
      DIMENSION SPININTS(INTSPACE)
*. Output
      DIMENSION H1DIA_U(*),H1DIA_B(*)
*. Just locally
      dimension IOBPTS_PIC(6+MXPR4T,MXPOBS)
      dimension IOBPTS2_PIC(6+MXPR4T,MXPOBS)
*
      NTESTL = 0010
      NTEST = MAX(NTESTL,IPRDIA)
*
      ZERO = 0.D0
      if (NTEST.ge.50) then
         write(6,*) '----------------------------------------'
         write(6,*) ' GETH1DIA4 speaking '
         write(6,*) '----------------------------------------'
      end if
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT1DIA')
      CALL MEMMAN(KLH1,NTOOB**2,'ADDL  ',2,'H1BL  ')
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
        do ISMOB=1,NSMOB/2,1
          do IGAS=1,NGAS,1
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
              write(6,'(2X,A,5I4)') 'IGAS,ISMOB,NORBSM,ISTART,IEND',
     &                               IGAS,ISMOB,NORBSM,ISTART,IEND
            end if
* Check for occupation of current space
            if (MXNELGS(IGAS).ne.0) then
              IUB = 1
              call geth1dia_sigden(IUB,WORK(KLH1),WORK(KT_CC),
     &                             WORK(KINDEX),
     &                             N1ELINT,NGSSH(ISMOB,IGAS),
     &                             ISTART,IEND,NTEST)
              call copvec(WORK(KLH1),H1DIA_U(IOFFSD),NGSSH(ISMOB,IGAS))
              IUB = 2
              call geth1dia_sigden(IUB,WORK(KLH1),WORK(KT_CC),
     &                             WORK(KINDEX),
     &                             N1ELINT,NGSSH(ISMOB,IGAS),
     &                             ISTART,IEND,NTEST)
              call copvec(WORK(KLH1),H1DIA_B(IOFFSD),NGSSH(ISMOB,IGAS))
            else
              call setvec(H1DIA_U(IOFFSD),ZERO,NGSSH(ISMOB,IGAS))
              call setvec(H1DIA_B(IOFFSD),ZERO,NGSSH(ISMOB,IGAS))
            end if
            IOFFSD = IOFFSD + NGSSH(ISMOB,IGAS)
          end do
        end do
      end if
      call memman(IDUM,IDUM,'PRINT ',IDUM,'      ')
*
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Diagonal integrals, unbarred '
        WRITE(6,*)
        CALL WRTMAT(H1DIA_U,1,NTOOB,1,NTOOB)
        WRITE(6,*)
        WRITE(6,*) ' Diagonal integrals, barred '
        WRITE(6,*)
        CALL WRTMAT(H1DIA_B,1,NTOOB,1,NTOOB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT1DIA')
      RETURN
      END
