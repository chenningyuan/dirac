*
*   Coefficient common block especially for Jeppe's routines
*
      integer NBSPFSM,NUSPFSM,IBSPEXP
      real*8 CSPEXP
      dimension CSPEXP(MXNSPIPDMGRP**2)
      dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
      dimension IBSPEXP(2,2)
      common /KRCC_COEFF/IBSPEXP,CSPEXP,NBSPFSM,NUSPFSM
*   
