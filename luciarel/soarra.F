* Arrange vector HELPMAT2 with spin-orbit integrals in analogy
* to spin-free integrals: two symmetry blocks on a triangle
* matrix
*
      subroutine soarra(SOINT,HELPMAT2,NIRR_PN,IDIMGRP,IDGIDIM,
     &                  NTRIAEL
*
      implicit real*8 (A-H,O-Z)
*
      dimension SOINT(IDGIDIM(IDIMGRP)**2)
      dimension HELPMAT2(NTRIAEL)
*
      if (IDBG.eq.4) then
*
*  determine offset for this run
         if (IDIMGRP.eq.1) then
            IOFF = 1
         else if (IDIMGRP.eq.2) then
            IOFF = 1 + ITRIADB(1)
            J = 0
            do I=IOFF,ITRIADB(2),1
               SOINT(J) = HELPMAT2(I)
               J = J+1
            end do

         end if




      else
         continue
      end if
      return
      end
