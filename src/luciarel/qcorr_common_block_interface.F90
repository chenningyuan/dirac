       module qcorr_interface

       use qcorr_cfg

       implicit none

       public init_qcorr
       public exit_qcorr


contains

      subroutine init_qcorr(relQC)

      logical, intent(in) :: relQC


      select case(relQC)

        case(.true.)
#ifdef PRG_DIRAC
          call init_qcorr_rel()
#else
          call quit('no relativistic CI in this program environment')
#endif
        case(.false.)
          call init_qcorr_nr()

      end select
      
      end subroutine init_qcorr

      subroutine exit_qcorr(relQC)

      logical, intent(in) :: relQC

      deallocate(KNSTSO_q)
      deallocate(KNSTSO2_q)
      deallocate(NELEC_q)
      if(relQC)then
        deallocate(MS2VAL_q)
        deallocate(NBLK_MS2_q)
        deallocate(IST_FOR_DT_q)
      end if
      
      end subroutine exit_qcorr

#ifdef PRG_DIRAC
      subroutine init_qcorr_rel()

 use symmetry_setup_krci
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"

! integer block
      MXNSTR_q    = MXNSTR
      NMS2VAL_q   = NMS2VAL
      NACTEL_q    = real(NACTEL)
      NORB_q      = 0
      NSMST_q     = NSMST
      

      allocate(KNSTSO_q(MXPSTT))
      KNSTSO_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
      
      allocate(KNSTSO2_q(MXPSTT))
      KNSTSO2_q(1:MXPSTT) = KNSTSO2(1:MXPSTT) 
     
      allocate(NELEC_q(MXPSTT))
      NELEC_q(1:MXPSTT) = NELEC(1:MXPSTT)

      allocate(MS2VAL_q(100))
      MS2VAL_q(1:100) = MS2VAL(1:100)

      allocate(NBLK_MS2_q(100))
      NBLK_MS2_q(1:100) = NBLK_MS2(1:100)

      allocate(IST_FOR_DT_q(2,MXPSTT))
      IST_FOR_DT_q(1:2,1:mxpstt) = IST_FOR_DT(1:2,1:MXPSTT)

      end subroutine init_qcorr_rel
#endif

      subroutine init_qcorr_nr()

#include "../lucita/mxpdim.inc"
#include "../lucita/lucinp.inc" ! contains NACTEL
#include "../lucita/strbas.inc"
#include "../lucita/cstate.inc"
#include "../lucita/strinp.inc"
#include "../lucita/stinf.inc"
#include "../lucita/cgas.inc"
#include "../lucita/gasstr.inc"
#include "../lucita/orbinp.inc"
#include "../lucita/csm.inc"


! integer block
      MXNSTR_q    = MXNSTR
      NMS2VAL_q   = 1 
      NACTEL_q    = real(NACTEL)
      NORB_q      = ntoob
      NSMST_q     = NSMST
      

      allocate(KNSTSO_q(MXPSTT))
      KNSTSO_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
      
      allocate(KNSTSO2_q(MXPSTT))
      KNSTSO2_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
     
      allocate(NELEC_q(MXPSTT))
      NELEC_q(1:MXPSTT) = NELEC(1:MXPSTT)

      end subroutine init_qcorr_nr

      end module

