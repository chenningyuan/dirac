      SUBROUTINE C_TO_CKK(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                    IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
*
* Version where matrix notation is eliminated
*
* Notice : K and T is in general batched. This is not visible here
*
*. Input
#include "implicit.inc"
      INTEGER    IKJA(NTA*NKA), IKJB(NTB,NKB)             
      DIMENSION  XKJA(NTA*NKA), XKJB(NTB,NKB)             
*. Input/Output 
C     DIMENSION C(NJA,NJB),CKK(NKA,NKB,*)
      DIMENSION C(NJA*NJB),CKK(*)
*
      DO ITB = 1, NTB
      DO ITA = 1, NTA
       ITAB = (ITB-1)*NTA+ITA-1
*
       KAKBITAB = ITAB*NKA*NKB
       IF(IWAY.EQ.1) THEN
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
          IF(JB.NE.0) THEN
           SIGNB = XKJB(ITB,KB)
           JAJB0 = (JB-1)*NJA 
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
            JA = IKJA(ITAKA)
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
                JAJB = JAJB0 + JA
                SIGN = XKJA(ITAKA)*SIGNB
                CKK(KAKBITAB) = SIGN*C(JAJB)
            ELSE 
               CKK(KAKBITAB) = 0.0D0
            END IF
           END DO
          ELSE IF (JB.EQ.0) THEN
           DO KAKBITAB_EFF =  KAKBITAB+1,  KAKBITAB + NKA  
               CKK(KAKBITAB_EFF) = 0.0D0
           END DO
           KAKBITAB  =  KAKBITAB + NKA
          END IF
*         ^ End if JB .NE. 0
         END DO
       ELSE 
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
          IF(JB.NE.0) THEN
           JAJB0 = (JB-1)*NJA
           SIGNB = XKJB(ITB,KB)
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
            JA = IKJA(ITAKA)
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
             JAJB = JAJB0 + JA
             SIGN = XKJA(ITAKA)*SIGNB
             C(JAJB) = C(JAJB) + SIGN*CKK(KAKBITAB)
            END IF
           END DO
          ELSE
            KAKBITAB =  KAKBITAB + NKA
          END IF
*         ^ End if JB = 0
         END DO
       END IF
*      ^ End of IWAY switch
      END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Output from C_TO_CKK, IWAY =  ', IWAY
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) ' C(Ja,Jb) '
        CALL WRTMAT(C,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' C(Ka,Kb,ITa, Itb) '
        ITAB = 0 
        DO ITB = 1, NTB
        DO ITA = 1, NTA
          ITAB = ITAB + 1
          WRITE(6,*) ' C(Ka,Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
          CALL WRTMAT(CKK(1+ITAB*NKA*NKB),NKA,NKB,NKA,NKB)
        END DO
        END DO
      END IF
*
      RETURN
      END 
