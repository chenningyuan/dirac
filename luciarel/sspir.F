      subroutine sspir(ISOI)
*
      implicit real*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "coeff.inc"
#include "crun.inc"
*
*   Coefficient common block especially for Jeppe's routines
*
*     integer NBSPFSM,NUSPFSM,IBSPEXP
*     real*8 CSPEXP
*     dimension CSPEXP(MXNSPIPDMGRP**2)
*     dimension NBSPFSM(MXNSPIR),NUSPFSM(MXNSPIR)
*     dimension IBSPEXP(2,2)
*     common /COEFF/CSPEXP,NBSPFSM,NUSPFSM,IBSPEXP
*
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   DOUGRP  : type of double group
*   NIRREP : Number of irreducible representations per double group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPSBS : Largest allowed number of spinor symmetries
*   MXNSPIR : Maximum number of spatial irreps
*   MXNSPIPDMGRP: Maximum number of spinors per dimension group
*
* ======
* Output
* ======
*
*
* Jeppe Olsen , Winter of 1991
*
* modified: Timo Fleig, Fall of 1996
*
      integer ISSPIR,NSSPIR
*. Output
*
*  Scratch
      do IRREP=1,NIRR_DG,1
         NBSPFSM(IRREP) = 0
         NUSPFSM(IRREP) = 0
      end do
*
      if (DOUGRP.eq.1) then
*=======
* D2h(*)
*=======
        NSMOB = 0
      else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
         write(*,*) 'These cases have not yet been implemented.'
         Call Abend1( 5 )
*=======
* C2v(*)
*=======
      else if (DOUGRP.eq.4) then
*=======
* C2h(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM(2) = ITOTDM(3) + ITOTDM(4)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         SYMDIM_F(2) = NORB_F(3) + NORB_F(4)
         DDIM = (SYMDIM(1))*2 + (SYMDIM(2))*2
         NBSPFSM(5) = ITOTDM(2)
         NBSPFSM(6) = ITOTDM(1)
         NBSPFSM(7) = ITOTDM(4)
         NBSPFSM(8) = ITOTDM(3)
         NUSPFSM(5) = ITOTDM(1)
         NUSPFSM(6) = ITOTDM(2)
         NUSPFSM(7) = ITOTDM(3)
         NUSPFSM(8) = ITOTDM(4)
*
         do I=1,2,1
            ITRIADB(I) = (SYMDIM(I)**2 + SYMDIM(I)) / 2
            ITRIADB_F(I) = (SYMDIM_F(I)**2 + SYMDIM_F(I)) / 2
         end do
         do I=1,2,1
            ITRIADBSUM_F = ITRIADBSUM_F + ITRIADB_F(I)
         end do
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
*=======
* C2(*) and Cs(*)
*=======
         SYMDIM(1) = ITOTDM(1) + ITOTDM(2)
         SYMDIM_F(1) = NORB_F(1) + NORB_F(2)
         DDIM = (SYMDIM(1))*2
         NBSPFSM(3) = ITOTDM(2)
         NBSPFSM(4) = ITOTDM(1)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.7) then
*=======
* Ci(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM(2) = ITOTDM(2)*2
         SYMDIM_F(1) = NORB_F(1)*2
         SYMDIM_F(2) = NORB_F(2)*2
         DDIM = SYMDIM(1) + SYMDIM(2)
         NBSPFSM(3) = ITOTDM(1)
         NBSPFSM(4) = ITOTDM(2)
         NUSPFSM(3) = ITOTDM(1)
         NUSPFSM(4) = ITOTDM(2)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else if (DOUGRP.eq.8) then
*=======
* C1(*)
*=======
         SYMDIM(1) = ITOTDM(1)*2
         SYMDIM_F(1) = NORB_F(1)*2
         NBSPFSM(2) = ITOTDM(1)
         NUSPFSM(2) = ITOTDM(1)
         ITRIADB(1) = (SYMDIM(1)**2 + SYMDIM(1)) / 2
         ITRIADBSUM = ITRIADB(1)
         ITRIADB_F(1) = (SYMDIM_F(1)**2 + SYMDIM_F(1)) / 2
         ITRIADBSUM_F = ITRIADB_F(1)
      else
         write(*,*) ' Sorry  DOUGRP out of range , DOUGRP = ', DOUGRP
         write(*,*) ' SSPIR fatally wounded '
         Call Abend1( 5 )
      end if
*
*  Determine spinor info and symmetry within the double group
*
Cjan      call spiinf(LU1)
*
Cjan      write(*,*) 'This is sspir speaking:'
Cjan      do I=1,DDIM,1
Cjan         write(*,*) 'Spinor ',I,' type is   ',SPINOR(I)
Cjan         write(*,*) 'Spinor symmetry is        ',SPISYM(I)
Cjan      end do
*
      return
      end
