      SUBROUTINE INPCTL(LUIN,LUOUT)
*
* Read and check input for LUCIA
*
*.1 : Position input file after  line containing &LUCIA &END
C      CALL SLASK
*.2 : read input in
      CALL LREADIN(LUIN,LUOUT)
      RETURN
      END
