      FUNCTION ITDIANUM(ICA,ICB,IAA,IAB,NC,NA)
*
* Find adress of element T(ICA,ICB,IAA,IAB) in diagonal 
* block
* It is assumed that (ICA,IAA) .ge. (ICB,IAB) has been checked outside
*
* Jeppe Olsen, July 2000, HNIE
*
#include "implicit.inc"
*. T(1,1,Iab,Iab)
      I11IABIAB = ((IAB-1)*NA - IAB*(IAB-1)/2)*NC*NC+
     &             (IAB-1)*NC*(NC+1)/2
*. T(1,1,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       I11IAAIAB =  I11IABIAB + NC*(NC+1)/2 + (IAA-IAB-1)*NC*NC
      ELSE
       I11IAAIAB =  I11IABIAB 
      END IF
*. T(Ica,Icb,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA
      ELSE
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA - ICB*(ICB-1)/2
      END IF
*
      ITDIANUM = IADR
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Address of diagonal T-element : '
        WRITE(6,*) ' ICA, ICB, IAA, IAB, NC,NA, IADR = ',
     &               ICA, ICB, IAA, IAB, NC,NA, IADR
      END IF
*
      RETURN
      END
