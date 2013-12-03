       module qcorr_interface

       implicit none

public init_qcorr


contains 

      subroutine init_qcorr()


       
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"




! integer block
      NGAS_q      = NGAS
      IBSPGPFTP_q = IBSPGFTP
      MXNSTR_q    = MXNSTR
      NMS2VAL_q   = NMS2VAL

      allocate(IGSOCC_q(MXPNGAS))
      IGSOCC_q(1:MXPNGAS) = IGSOCC(1:MXPNGAS)

      allocate(NELFSPGP_q(MXPNGAS,MXPSTT))
      NELFSPGP(1:MXPNGAS,1:MXPSTT) = NELFSPGP(1:MXPNGAS,1:MXPSTT)

      allocate(KNSTSO_q(MXPSTT))
      KNSTSO_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
      
      allocate(KNSTSO2_q(MXPSTT))
      KNSTSO2_q(1:MXPSTT) = KNSTSO2(1:MXPSTT) 
     
      allocate(MS2VAL_q(100))
      MS2VAL_q(1:100) = MS2VAL(1:100)

      allocate(NBLK_MS2_q(100))
      NBLK_MS2_q(1:100) = NBLK_MS2(1:100)

      allocate(NELEC_q(MXPSTT))
      NELEC_q(1:MXPSTT) = NELEC(1:MXPSTT)

      allocate(NOCTYP_q(MXPSTT))
      NOCTYP_q(1:MXPSTT) = NOCTYP(1:MXPSTT)

      allocate(IST_FOR_DT_q(MXPSTT))
      IST_FOR_DT_q(2,MXPSTT) = IST_FOR_DT(2,MXPSTT)

      end subroutine init_qcorr

      end module

