      SUBROUTINE SPOBEX_TO_ABOBEX(ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &                            IFLAG,NAOBEX_TP,NBOBEX_TP,IAOBEX_TP,
     &                            IBOBEX_TP)
*
* Split spin-orbital excitations into alpha and beta-orbital excitations
*
* IFLAG = 1 : Find only number of alpha- and beta- orbital excitations
*
* Jeppe Olsen, July 2000
*
#include "implicit.inc"
*
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
*. Output
      INTEGER IAOBEX_TP(2*NGAS,*),IBOBEX_TP(2*NGAS,*)
*
C?    WRITE(6,*) ' SPOBEX_TO_ABOBEX : NSPOBEX_TP = ', NSPOBEX_TP
      DO IAB = 1, 2
        LEX_TP = 0
        DO JSPOBEX_TP = 1, NSPOBEX_TP
*. Has this a or b excitation been observed before
          I_AM_OLD_HAT = 0
          DO KSPOBEX_TP = 1, LEX_TP
            I_DENTICAL = 1
            DO ICA = 1, 2
              IEXP_OFF = (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
              DO IGAS = 1, NGAS
                IF(ISPOBEX_TP(IGAS+IEXP_OFF-1,JSPOBEX_TP).NE.
     &             ISPOBEX_TP(IGAS+IEXP_OFF-1,KSPOBEX_TP)    )
     &             I_DENTICAL = 0
              END DO
            END DO
            IF(I_DENTICAL.EQ.1) I_AM_OLD_HAT = 1
          END DO
          IF(I_AM_OLD_HAT.EQ.0) THEN
            LEX_TP = LEX_TP + 1
            IF(IFLAG.EQ.0) THEN
              DO ICA = 1, 2
               IEXP_OFF =  (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
               IAB_OFF  =  (ICA-1)*NGAS + 1
               IF(IAB.EQ.1) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IAOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               ELSE IF (IAB.EQ.2) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IBOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               END IF
              END DO
            END IF
          END IF
*         ^ End of I_AM_OLD_HAT = 0
        END DO
*       ^ End of loop over JSPOBEX_TP
      IF(IAB.EQ.1) THEN
        NAOBEX_TP = LEX_TP
      ELSE
        NBOBEX_TP = LEX_TP
      END IF
      END DO
*     ^ End of loop over IAB
*
      NTEST = 000
      IF(NTEST.GE.3) THEN
       WRITE(6,*) ' Number of alpha-excitation operators ', NAOBEX_TP
       WRITE(6,*) ' Number of beta-excitation operators ',  NBOBEX_TP
      END IF
      IF(NTEST.GE.5.AND.IFLAG.EQ.0) THEN
        WRITE(6,*) ' Alpha-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NAOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' alphaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IAOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IAOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' beta-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NBOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' betaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IBOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IBOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
      END IF
*     ^ End of print is active
*
      RETURN
      END
