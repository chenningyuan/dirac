      SUBROUTINE GET_INT2(XOUT,XINT,ISEA,JSEA,KSEA,LSEA,
     &                                   NALPHA)
*
* Routine to fetch integral one by one from CC code.
* Done by calculating an offset on XINT, accieved by arranging indices
* Integral storage mode is assumed the default for luciarel hence integral indexing
* 
C  Plan is as follows
C 1. Identify which GAS space each index belong to
C 2. Make sure that this integral type will be there (only for NALPHA = 0 or 4)
C 3. Identify TYPE of integral
C 4. Calculate integral offset 
C 5. Fetch integral
C
C Could perhaps just pass down IHTYPE instead of finding it?
C ISUMORB and ISUBLENGHT could be replaced by Gauss's formula
C
C ISEA = C1
C JSEA = A1
C KSEA = C2
C LSEA = A2
C
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Local scratch
      INTEGER I_INDX_TO_GAS(4), ISUBLENGHT(NGAS+1)
      INTEGER IORB(NGAS+1), ITYPE2ALPHA(3), IPOS2ALPHA(3)
      DIMENSION XINT(*)
*
* Set some local variables
*
C     print*,'ISEA,JSEA,KSEA,LSEA',ISEA,JSEA,KSEA,LSEA
      NTEST = 000
      SIGN2 = 1.0D0
* 1. Identify the GAS and mark those with one orbital
      IORBCOUNT = 0
      ILENGHT = 0
      IEND = 0
      JEND = 0
      KEND = 0
      LEND = 0
      ISWAP = 0
      ISUBLENGHT(1) = 0
      IORB(1) = 0
      DO IGAS = 1,NGAS
        IORBCOUNT = IORBCOUNT + NGSOBT(IGAS)
        IF(IEND.EQ.0) THEN
          IF(IORBCOUNT.GE.ISEA) THEN
            I_INDX_TO_GAS(1) = IGAS
            IEND = 1
          END IF
        END IF
        IF(JEND.EQ.0) THEN
          IF(IORBCOUNT.GE.JSEA) THEN
            I_INDX_TO_GAS(2) = IGAS
            JEND = 1
          END IF
        END IF
        IF(KEND.EQ.0) THEN
          IF(IORBCOUNT.GE.KSEA) THEN
            I_INDX_TO_GAS(3) = IGAS
            KEND = 1
          END IF
        END IF
        IF(LEND.EQ.0) THEN
          IF(IORBCOUNT.GE.LSEA) THEN
            I_INDX_TO_GAS(4) = IGAS
            LEND = 1
          END IF
        END IF
        ILENGHT = ILENGHT + IGAS
        ISUBLENGHT(IGAS+1) =  ILENGHT
        IORB(IGAS+1) = IORBCOUNT
      END DO
*
* 2. Finding the index of this type of integral in luciarel format
*    We here want the following condition on indices l>j and k>i for NALPHA 4 or 0
*    NALPHA 2 is not sorted like this
*
C      WRITE(6,*) 'ISEA,JSEA,KSEA,LSEA,start = ',ISEA,JSEA,KSEA,LSEA
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
        IF(ABS(ISEA).GT.ABS(KSEA)) THEN
          IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap all 4 indices and I_INDX_TO_GAS
C           print*,'enter 1'
            ISWAP = I_INDX_TO_GAS(1)
            I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
            I_INDX_TO_GAS(3) = ISWAP
            ISWAP = I_INDX_TO_GAS(2)
            I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
            I_INDX_TO_GAS(4) = ISWAP
            ISWAP = ISEA
            ISEA = KSEA
            KSEA = ISWAP
            ISWAP = JSEA
            JSEA = LSEA
            LSEA = ISWAP
            SIGN2 = -1.0D0 
          ELSE
C swap only 2 indices and I_INDX_TO_GAS
C           print*,'enter 2'
            ISWAP = I_INDX_TO_GAS(1)
            I_INDX_TO_GAS(1) = I_INDX_TO_GAS(3)
            I_INDX_TO_GAS(3) = ISWAP
            ISWAP = ISEA
            ISEA = KSEA
            KSEA = ISWAP
            SIGN2 = 1.0D0 !?
          END IF
        ELSE IF(ABS(JSEA).GT.ABS(LSEA)) THEN
C swap only 2 indices and I_INDX_TO_GAS
C           print*,'enter 3'
            ISWAP = I_INDX_TO_GAS(2)
            I_INDX_TO_GAS(2) = I_INDX_TO_GAS(4)
            I_INDX_TO_GAS(4) = ISWAP
            ISWAP = JSEA
            JSEA = LSEA
            LSEA = ISWAP
            SIGN2 = 1.0D0
        ELSE
C           print*,'enter 4'
            SIGN2 = -1.0D0
        END IF
      END IF
*
* 3. Identify relative type of integral and find offset of type and integral from lists
* Remember to check if there is only one kramers pair in gas shell. Done in step 1
* Should also add relative offset of type. Done
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
        ITYPE = ISUBLENGHT(I_INDX_TO_GAS(4)) + I_INDX_TO_GAS(2)-1 +
     & ILENGHT*(ISUBLENGHT(I_INDX_TO_GAS(3)) + I_INDX_TO_GAS(1)-1)
*
        IF(NALPHA.EQ.4) THEN
          ITYPE = ITYPE + IRELTYPE(3)
C        print*,'ITYPE ABS',ITYPE
        ELSE IF(NALPHA.EQ.0) THEN
C          print*,'ITYPE,IRELTYPE(4)',ITYPE,IRELTYPE(4)
          ITYPE = ITYPE + IRELTYPE(4)
        ELSE
          ITYPE = ITYPE + IRELTYPE(8)
        END IF
      ELSE
C
C Here we need coulomb exchange2 and spin-orbit. Spin-orbit is like above
C Spin-orbit is here a special case since it is sorted like for NALPHA 4 or 0
C Will just make a small subroutine that goes like above to find it (if it is there!)
C Placed at the very end since no duplication of indices is then needed
        ITYPE =           (I_INDX_TO_GAS(2)-1) + 
     &          NGAS    * (I_INDX_TO_GAS(4)-1) + 
     &          NGAS**2 * (I_INDX_TO_GAS(3)-1) +
     &          NGAS**3 * (I_INDX_TO_GAS(1)-1)
        ITYPE2ALPHA(1) = ITYPE + IRELTYPE(5)
C       ITYPE =           (I_INDX_TO_GAS(4)-1) +
C    &          NGAS    * (I_INDX_TO_GAS(2)-1) +
C    &          NGAS**2 * (I_INDX_TO_GAS(3)-1) +
C    &          NGAS**3 * (I_INDX_TO_GAS(1)-1)
C       ITYPE2ALPHA(2) = ITYPE + IRELTYPE(5) +
C    &                  (IRELTYPE(6)-IRELTYPE(5))/2
C       print*,'ITYPE2ALPHA',ITYPE2ALPHA(1)!,ITYPE2ALPHA(2)
C       print*,'IRELTYPE(5),IRELTYPE(6)',IRELTYPE(5)!,IRELTYPE(6)
      END IF
* 4. Calculate position of integral
* First calculate position of integral in type then add offset for type
C Number of integrals can be found from array
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
C Four possibilties
        IF(I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3).AND.
     &     I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
          IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &           ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &          (JSEA - IORB(I_INDX_TO_GAS(2)) - 1) * 
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &           ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C          print*,'JSEA,LSEA,ISEA,KSEA',JSEA,LSEA,ISEA,KSEA
C          print*,'IORB(I_INDX_TO_GAS(2)-1)',IORB(I_INDX_TO_GAS(2)-1)
C          print*,'IPOS for 1',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 1',IPOS,IOFFINTTYPE(ITYPE)
        ELSE IF (I_INDX_TO_GAS(1).NE.I_INDX_TO_GAS(3).AND.
     &     I_INDX_TO_GAS(2).NE.I_INDX_TO_GAS(4)) THEN
          IPOS = (JSEA - IORB(I_INDX_TO_GAS(2))-1) * 
     &                   NGSOBT(I_INDX_TO_GAS(1)) * 
     &                   NGSOBT(I_INDX_TO_GAS(3)) * 
     &                   NGSOBT(I_INDX_TO_GAS(4)) +
     &           (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &           (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &            KSEA - IORB(I_INDX_TO_GAS(3))-1
C          print*,'JSEA,LSEA,ISEA,KSEA',JSEA,LSEA,ISEA,KSEA
C          print*,'IORB(I_INDX_TO_GAS(2)-1)',IORB(I_INDX_TO_GAS(2)-1)
C          print*,'IPOS',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 2',IPOS,IOFFINTTYPE(ITYPE)
C          print*,'
        ELSE IF (I_INDX_TO_GAS(1).EQ.I_INDX_TO_GAS(3)) THEN
          IPOS = ISEA - IORB(I_INDX_TO_GAS(1)) - 1 +
     &           ISUMORB(KSEA - IORB(I_INDX_TO_GAS(3)) - 1) +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(4))  *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3))) +
     &           (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &           ISUMORB(NGSOBT(I_INDX_TO_GAS(3)))
C          print*,'IPOS for 3',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) for 3',IPOS,IOFFINTTYPE(ITYPE)
        ELSE IF (I_INDX_TO_GAS(2).EQ.I_INDX_TO_GAS(4)) THEN
          IPOS = (ISEA - IORB(I_INDX_TO_GAS(1))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &            KSEA - IORB(I_INDX_TO_GAS(3))-1 +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) +
     &           ISUMORB(LSEA - IORB(I_INDX_TO_GAS(4)) - 1) *
     &                   NGSOBT(I_INDX_TO_GAS(1)) *
     &                   NGSOBT(I_INDX_TO_GAS(3)) 
C          print*,'IPOS 4',IPOS
          IPOS = IPOS + IOFFINTTYPE(ITYPE)
C          print*,'IPOS,IOFFINTTYPE(ITYPE) 4',IPOS,IOFFINTTYPE(ITYPE)
        END IF
      ELSE
          IPOS = (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) *
     &                        NGSOBT(I_INDX_TO_GAS(2)) +
     &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) *
     &                        NGSOBT(I_INDX_TO_GAS(3)) +
     &           (KSEA - IORB(I_INDX_TO_GAS(3))-1) *
     &                        NGSOBT(I_INDX_TO_GAS(1)) +
     &            ISEA - IORB(I_INDX_TO_GAS(1))-1
C        print*,'ISEA,JSEA,KSEA,LSEA',ISEA,JSEA,KSEA,LSEA
        IPOS2ALPHA(1) = IPOS + IOFFINTTYPE(ITYPE2ALPHA(1))
C         IPOS = (LSEA - IORB(I_INDX_TO_GAS(4))-1) *
C    &                        NGSOBT(I_INDX_TO_GAS(1)) *
C    &                        NGSOBT(I_INDX_TO_GAS(3)) *
C    &                        NGSOBT(I_INDX_TO_GAS(2)) +
C    &           (JSEA - IORB(I_INDX_TO_GAS(2))-1) *
C    &                        NGSOBT(I_INDX_TO_GAS(1)) *
C    &                        NGSOBT(I_INDX_TO_GAS(3)) +
C    &           (KSEA - IORB(I_INDX_TO_GAS(3))-1) *
C    &                        NGSOBT(I_INDX_TO_GAS(1)) +
C    &            ISEA - IORB(I_INDX_TO_GAS(1))-1
C       IPOS2ALPHA(2) = IPOS + IOFFINTTYPE(ITYPE2ALPHA(2))
C        print*,'IPOS',IPOS
C see if there is a spin-orbit contribution and find it!
C       IF(ISPINORBIT.EQ.2) THEN
C       IF(ISEA.NE.KSEA.AND.JSEA.NE.LSEA) THEN
C         IPOS2ALPHA(3) = 0
C       END IF
C       ELSE
C         IPOS2ALPHA(3) = 0
C       END IF
C        print*,'IPOS2ALPHA',IPOS2ALPHA(1),IPOS2ALPHA(2),IPOS2ALPHA(3)
      END IF
* 5. Fetch integral from list
C
C      print*,'ipos',ipos
C      print*,'i am here'
      IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
        XOUT = SIGN2*XINT(IPOS)
C       print*,'IPOS = ',IPOS
      ELSE
C       IF(IPOS2ALPHA(3).EQ.0) THEN
          XOUT = XINT(IPOS2ALPHA(1)) !-
C         print*,'IPOS2ALPHA(1) = ',IPOS2ALPHA(1)
C    &           XINT(IPOS2ALPHA(2))
C       ELSE
C         XOUT = XINT(IPOS2ALPHA(1)) - 
C    &           XINT(IPOS2ALPHA(2)) !+
C    &    SIGNS* XINT(IPOS2ALPHA(3))
C       END IF
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'IORBCOUNT = ', IORBCOUNT 
        WRITE(6,*) 'ILENGHT = ',ILENGHT
        DO IGAS = 1,NGAS
          WRITE(6,*) 'ISUBLENGHT(IGAS) = ', ISUBLENGHT(IGAS)
        END DO
        DO IGAS = 1,NGAS
          WRITE(6,*) 'IORB(IGAS) = ', IORB(IGAS)
        END DO
        DO I = 1,4
          WRITE(6,*) 'I_INDX_TO_GAS(I) = ',I_INDX_TO_GAS(I)
        END DO
        WRITE(6,*) 'ISEA,JSEA,KSEA,LSEA = ',ISEA,JSEA,KSEA,LSEA
        IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
          WRITE(6,*) 'ITYPE = ',ITYPE
        ELSE
          DO I=1,1!3
            WRITE(6,*) 'ITYPE2ALPHA(I) = ',ITYPE2ALPHA(I)
          END DO
        END IF
        IF(NALPHA.EQ.4.OR.NALPHA.EQ.0) THEN
          WRITE(6,*) 'IPOS = ',IPOS
        ELSE
          DO I=1,1!3
            WRITE(6,*) 'IPOS2ALPHA(I) = ',IPOS2ALPHA(I)
          END DO
        END IF
        WRITE(6,*) 'XOUT = ',XOUT
      END IF
      WRITE(6,*) 'ISEA,JSEA,KSEA,LSEA = ',ISEA,JSEA,KSEA,LSEA
*
      RETURN
      END
*
      SUBROUTINE GET_INTERMEDIATE(XOUT,XINT)
*
* Routine to fetch intermediates on an array
*
#include "implicit.inc"
*
      DIMENSION XINT(*)
*
      RETURN
      END
