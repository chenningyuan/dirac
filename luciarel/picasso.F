****************************************************************
***                     PICASSO                              ***
***                                                          ***
***   P rogram						     ***
***     for 						     ***
***   I ntelligent                                           ***
***   C omplete                                              ***
***   A ctive                                                ***
***   S pace                                                 ***
***   S pin                                                  ***
***   O rbit                                                 ***
***     calculations                                         ***
***         (at least what it was once planned to become...) ***
***                                                          ***
****************************************************************
***
*** Written by Timo Fleig, Bonn
***        and Jeppe Olsen, Lund
***
*** October 1996 - February 1998
***
*** (poor picasso is merely one of LUCIAREL's subroutines ... :-(
***
      subroutine picasso(INTSPACE,ISOOP,IPRTRA,MCSCF,
     &                   NL2D,IUB2D,ISM2D,F1,F2,IBOSYM,CIRUN)
*
*******************************************************************
*  !!!  INATIT states iteration point for natural spinor method !!!
*******************************************************************
*
      implicit real*8 (A-H,O-Z)
*
      integer POI(30)
      logical MCSCF
      common /WORKPIC/IPOI
*
#include "mxpdim.inc"
#include "input.inc"
#include "symm.inc"
#include "onel.inc"
#include "crun.inc"
*
* for integral summation
      ISUMNEW = 1
* for global purposes
      INTITER = 1
      N_SPINR_EXP = 2
*
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *  This is PICASSO,                       *'
C     write(6,*) '     *  (Transformation section for            *'
C     write(6,*) '     *   one-electron sf and SO integrals      *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *        0.1   stoneage version           *'
C     write(6,*) '     *        0.2   GAS version                *'
C     write(6,*) '     *        0.2.1  AO ints from MOLECULE     *'
C     write(6,*) '     *        0.2.2  SOMF ints from AMFI       *'
C     write(6,*) '     *        0.2.2.1 GAS-DG dependent         *'
C     write(6,*) '     *        0.2.3  general trafo CI/MCSCF    *'
C     write(6,*) '     *        0.2.3.1 SGI portable version     *'
C     write(6,*) '     *        0.3   MEMMAN allocation included *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
*
C     write(6,*) '     *******************************************'
C     write(6,*) '     *                                         *'
C     write(6,*) '     * written by                              *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *      TIMO FLEIG, University of Bonn     *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *           and                           *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *      JEPPE OLSEN, University of Lund    *'
C     write(6,*) '     *                                         *'
C     write(6,*) '     *******************************************'
*
*  Determine double group information
      call sspir(ISOI)
*
*  Move to integral transformation section:
*  Define type of transformation
*  Here for now: complete two-component trans for CI
*
      if (MCSCF) then
        ITRA_TYPE = 4
      else
        if (INATIT.eq.1) then
          ITRA_TYPE = 2
        else if (INATIT.ge.2) then
          ITRA_TYPE = 4
          ORBFILE = .true.
        end if
      end if
*
      call itrctl(INTSPACE,ITRA_TYPE,INTITER,ENVIRO,
     &            ISPINFREE,IFOCKEX,NORB_F,NL2D,IUB2D,ISM2D,ISYM_T,
     &            SIGDEN_ROUTE,F1,F2,IBOSYM,MCSCF,CIRUN)
*
      return
      end
