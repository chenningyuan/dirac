*  Import integrals from MOLCAS/SEWARD environment
*  For run with externally generated Fock matrix.
*  (Driver routine)
*
      subroutine rd_traone(ONEINT,LU2,FNONEMO,IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
*
      character*6 FNONEMO
      dimension ONEINT(ITRIASUM)
*
      call geth0(ONEINT(1),LU2,NIRR_PN,SYDI(1),
     &           NFRO_L(1),FNONEMO,IPRINT)
*
      return
      end
