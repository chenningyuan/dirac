*
*  Arrange density matrices for diagonalization
*  Input matrices are GAS-symmetry ordered.
*
*   1) Loop over output spinor symmetries (done outside!)
*    2) Loop over real/imaginary part (done outside!)
*     3) Loop over unbarred/barred type of left index
*      4) Loop over unbarred/barred type of right index
*         Pick relevant symmetry blocks from GAS blocks
*         and arrange them according to the spinor symmetry.
*
*   Timo Fleig
*    fixed bug in UB handling (05.11.2001)
*
      subroutine arr_dens(RHO1UU,RHO1UB,RHO1BU,RHO1BB,
     &                    NGSOB,NGSOB2,
     &                    NGSOB_HE,NGSOB2_HE,BLKOUT,
     &                    READEN,IDBG,ISMDIM,
     &                    ITDIM,ISPSM,NIRR,IGRP,
     &                    ISPSF_MO,ITOTDM,NGAS,IPRDEN)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
*
      dimension RHOSYM(ISMDIM,ISMDIM,NIRR)
      dimension BLKOUT(*),READEN(ISMDIM**2)
      dimension NGSOB_HE(2,NGAS),NGSOB2_HE(2,NGAS),
     &          NGSOB(MXPOBS,MXPNGAS),NGSOB2(MXPOBS,MXPNGAS)
      dimension ISPSF_MO(MXNDGIRR,2),ITOTDM(MXNSPIR)
*
      NTESTL = 20
      NTEST = max(IPRDEN,NTESTL)
*
*  Copy matrix with number of unbarred spinors and number of
*  unbarred spinors (2) per (dgsym,GAS) to help array.
      do ISYM=1,2,1
        do IGAS=1,NGAS,1
          NGSOB_HE(ISYM,IGAS) = NGSOB(ISYM+NIRR+2*(IGRP-1),IGAS)
          NGSOB2_HE(ISYM,IGAS) = NGSOB2(ISYM+NIRR+2*(IGRP-1),IGAS)
        end do
      end do
*
      if (NTEST.ge.20) then
        write(6,*) 'Checking help matrices NGSOB(2)_HE: '
        call iwrtma(NGSOB_HE,NIRR,NGAS,NIRR,NGAS)
        call iwrtma(NGSOB2_HE,NIRR,NGAS,NIRR,NGAS)
      end if
*
      if (IDBG.ge.4.and.IDBG.le.6) then
        do IUB=1,2,1
          do JUB=1,2,1
* UB = 1 means unbarred index
* UB = 2 means barred index
            if (NTEST.ge.50) write(6,*) ' Processing UB block ',IUB,JUB
            IDIMR = ITOTDM(ISPSF_MO(ISPSM,IUB))
            IDIMC = ITOTDM(ISPSF_MO(ISPSM,JUB))
            if (IUB.eq.1.and.JUB.eq.1) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB_HE,
     &                   NGSOB_HE,IUB,JUB,TRDUM,TRDUM,RHO1UU,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = 1
                ICOLOFF = 1
              else if (IEOSS.eq.2) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.1.and.JUB.eq.2) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB_HE,
     &                   NGSOB2_HE,IUB,JUB,TRDUM,TRDUM,RHO1UB,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              else if (IEOSS.eq.2) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.2.and.JUB.eq.1) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB2_HE,
     &                   NGSOB_HE,IUB,JUB,TRDUM,TRDUM,RHO1BU,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = 1
              else if (IEOSS.eq.2) then
                IROWOFF = 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.2.and.JUB.eq.2) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB2_HE,
     &                   NGSOB2_HE,IUB,JUB,TRDUM,TRDUM,RHO1BB,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              else if (IEOSS.eq.2) then
                IROWOFF = 1
                ICOLOFF = 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            end if
          end do
        end do
      else if (IDBG.ge.7.and.IDBG.le.8) then
        write(6,*) ' Double group ',IDBG
        write(6,*) ' not inplemented yet in density calculation. '
        stop 'Quitting.'
      end if
*
      if (NTEST.ge.10) then
        write(6,*)
        write(6,*) 'Checking output block READEN in arr_dens.'
        call wrtmat(READEN,ISMDIM,ISMDIM,ISMDIM,ISMDIM)
      end if
*
      return
      end
