      SUBROUTINE GXFTP(NOP,IAC,ILOCC,IROCC,NGAS,NX,IX,MXNX)
*
* Occupation of left and right string given by ILOCC,IROCC
* These are connected with a string containing NOP operators.
*
* <ILOCC! IOP ! IROCC>
*
* Find possible gasspaces of the each operator
*
* Jeppe Olsen, Winter of 1998
*
*      IMPLICIT REAL*8(A-H,O-Z)
*
* ======
*. Input
* ======
*
* NOP : Number of operators
* IACI(I) : =1 => operator I is annihilation
*           =2 => operator I is creation
* ILOCC   : Occupation of left  string
* IROCC   : Occupation of right string
*
* ======
* Output
* ======
*
* NX : Number of connecting excitations
* IX : Connecting excitations
* MXNX : Max number of connecting excitations
*
*. Input
CTI  This array has to be passed from above or put on a common block
CTI      DIMENSION IOP(NOP)
      INTEGER ILOCC(NGAS),IROCC(NGAS)
*. output
      DIMENSION IX(MXNX,NOP)
*
      RETURN
      END
