       SUBROUTINE ASSIGN_POINTERS77(IW,RW,LW,IVW,N1,N2,NNZ)
C
C
C in Fortran77 style !!!
C
       use pointers

#include "implicit.h"
#include "priunit.h"

       INTEGER, TARGET:: IW(N1,N2,NNZ),IVW(N1)
       REAL*8,TARGET :: RW(N1,N2,NNZ)
       LOGICAL, TARGET :: LW(N1,N2,NNZ)


C .. fill arrays !!!

       DO I=1,N1
            IVW(I) = I-10
       DO J=1,N2
       DO IZ=1,NNZ
         IW(I,J,IZ) = I*2+J+2*IZ-100
         RW(I,J,IZ) = DFLOAT(I+J+2*IZ)+0.250
         IF (I.NE.J) THEN
          LW(I,J,IZ) = .TRUE.
         ELSE
          LW(I,J,IZ) = .FALSE.
         ENDIF
       ENDDO
       ENDDO
       ENDDO

C       CALL ARRANGE_POINTERS90(IW,RW,LW)
C      ...finally  assign pointers here
         IA => IW
         IVA => IVW
         RA => RW
         LA => LW

       RETURN
       END


