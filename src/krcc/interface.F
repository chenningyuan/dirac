      SUBROUTINE TRANSFER_KRCC_TO_KRMC(NKRCCGAS,NKRCC_MK2REF,
     &                           NKRCC_MK2DEL,NKRCCAELEC,NFSYM,
     &                           NKRCC_MINMK2,NKRCC_MAXMK2,NKRCCISH,
     &                           NKRCCASH,NKRCCGSH,NKRCCSSH,NKRCCFRO,
     &                           NKRCCGSP,MAXCIT_KRCC)
*
* Transfer common block info from KRCC to KRMC
* Lasse 2010   
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "crun.inc"
#include "lucinp.inc"
*
      INTEGER NKRCCISH(2),NKRCCASH(2),NKRCCSSH(2)
      INTEGER NKRCCFRO(2),NKRCCGSP(2,MXPNGAS),NKRCCGSH(2,MXPNGAS)
*
      NTEST = 100
*
* Initialize
      DO I =1,NFSYM
        NKRCCISH(I) = 0
        NKRCCASH(I) = 0
        NKRCCSSH(I) = 0
        NKRCCFRO(I) = 0
        DO J = 1,NGAS
          NKRCCGSH(I,J) = 0
          NKRCCGSP(J,I) = 0
        END DO
      END DO
*
      NKRCCGAS     = NGAS
      NKRCC_MK2REF = MK2REF
      NKRCC_MK2DEL = MK2DEL
      NKRCC_MINMK2 = MK2REF - MK2DEL
      NKRCC_MAXMK2 = MK2REF + MK2DEL
      NKRCCAELEC   = NACTEL
      MAXCIT_KRCC  = MAXIT
*
      DO I =1,NFSYM
        NKRCCISH(I) = KRCC_NISH(I)
        DO J = 1,NGAS
          NKRCCASH(I) = NKRCCASH(I) + NGSSH(I,J)
          NKRCCGSH(I,J) = NGSSH(I,J)
        END DO
        NKRCCSSH(I) = 0
        NKRCCFRO(I) = 0
      END DO
*
      DO I =1,NGAS
        DO J =1,2
          NKRCCGSP(J,I) = IGSOCCX(J,I,1)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ALL INFORMATION TRANSFERED '
        WRITE(6,*) ' NKRCCGAS ',NKRCCGAS
        WRITE(6,*) ' NKRCC_MK2REF ',NKRCC_MK2REF
        WRITE(6,*) ' NKRCC_MK2DEL ',NKRCC_MK2DEL
        WRITE(6,*) ' NKRCC_MINMK2 ',NKRCC_MINMK2
        WRITE(6,*) ' NKRCC_MAXMK2 ',NKRCC_MAXMK2
        WRITE(6,*) ' NKRCCAELEC   ',NKRCCAELEC
        WRITE(6,*) ' MAXCIT_KRCC  ',MAXCIT_KRCC
        DO I =1,NFSYM
          WRITE(6,*) ' NKRCCISH ',NKRCCISH(I),' I ',I
        END DO
        DO I =1,NFSYM
          WRITE(6,*) ' NKRCCASH ',NKRCCASH(I),' I ',I
        END DO
        DO I =1,NFSYM
          DO J = 1,NGAS
            WRITE(6,*) ' NKRCCGSH ',NKRCCGSH(I,J),' I ',I,' J ',J
          END DO
        END DO
        DO I =1,NFSYM
          DO J = 1,NGAS
            WRITE(6,*) ' NKRCCGSP ',NKRCCGSP(I,J),' I ',I,' J ',J
          END DO
        END DO
        DO I =1,NFSYM
          WRITE(6,*) ' NKRCCSSH ',NKRCCSSH(I),' I ',I
        END DO
        DO I =1,NFSYM
          WRITE(6,*) ' NKRCCFRO ',NKRCCFRO(I),' I ',I
        END DO
      END IF
*
      RETURN
      END
