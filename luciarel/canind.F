      INTEGER FUNCTION CANIND(I,J)
C
      IF(I.GT.J) THEN
       CANIND=I*(I-1)/2 + J
      ELSE
       CANIND=J*(J-1)/2 + I
      END IF
      RETURN
      END
