*  This little routine transfers input information
*  from LUCIA to PICASSO, so that PICASSO no longer
*  needs an input of his own ...
*
      subroutine luc2pic(ISOOP,IPRTRA)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cgas.inc"
#include "symm.inc"
#include "input.inc"
#include "crun.inc"
*
      NGASPI = NGAS
      if (NGAS.gt.MXPNGAS) then
        write(6,*) 'Too many GAS spaces.'
        write(6,*) 'Number of GAS spaces is ',NGAS
        write(6,*) 'Increase MXPNGAS in mxpdim.inc or reconsider.'
        Call Abend2( 'Quitting.' )
      end if
*
      if (PNTGRP.eq.9) DOUGRP = 4
      if (PNTGRP.eq.8) DOUGRP = 5
      if (PNTGRP.eq.7) DOUGRP = 6
      if (PNTGRP.eq.6) DOUGRP = 7
      if (PNTGRP.eq.5) DOUGRP = 8
*
      NOBPIRR = 0
      IOFPGSM(1) = 0
      do IRR=1,NIRR_PN,1
         do IGAS=1,NGASPI,1
            NOBPIRR = NOBPIRR + NGSSH(IRR,IGAS)
         end do
         if (IRR.lt.NIRR_PN) IOFPGSM(IRR+1) = IOFPGSM(IRR) + NOBPIRR
         ITOTDM(IRR) = NOBPIRR
         SYDI(IRR) = ITOTDM(IRR)
         NOBPIRR = 0
      end do
*
* Determine largest symmetry block
* Determine triangle lengths
* ... and the sum of these
      J = 1
      IMAXDIM = J
      ITRIASUM = 0
      ITRIASUM_F = 0
      do I=1,NIRR_PN,1
         if (SYDI(I).gt.SYDI(J)) then
            IMAXDIM = I
            J = I
         end if
*
         ITRIA(I) = (ITOTDM(I)**2 + ITOTDM(I)) / 2
         ITRIA_F(I) = (NORB_F(I)**2 + NORB_F(I)) / 2
         ITRIASUM = ITRIASUM + ITRIA(I)
         ITRIASUM_F = ITRIASUM_F + ITRIA_F(I)
*
         do J=1,NGASPI,1
            NGSOBPI(I,J) = NGSSH(I,J)
         end do
      end do
*
* ISO1 = ISOOP in readin.
* '0' means inactive, '1' means active SO operator of this type.
* ISO2 not implemented yet.
      if (ISOOP.eq.0) IOPER = 1   ! Only sf one-electron operator
      if (ISOOP.eq.1) IOPER = 3   ! sf+SO one-electron operator
C     Only SO operator doesn't make any sense!
*
      IPRINT = IPRTRA
*
      return
      end
