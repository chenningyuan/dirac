      subroutine tweltra(I12O,I34O,I1234O,ICCSMO,
     &                   IRIO,ISPITRA_RC,ISPITRA,NAOFSM,
     &                   LSPINOR,ITRA_TYPE,ITER)
*
* Test general spinor transformation program
*
* Jeppe + Timo, April-July 1997
*
* ITRA_TYPE Defines the type of transformation
* ITRA_TYPE = 1 : Normal nonrelativistic transformation
* ITRA_TYPE = 2 : 2-component relativistic transformation for CI,
*                 all nonredundant integral lists
* ITRA_TYPE = 3 : 2-component relativistic transformation for CI,
*                 Excluding integrals with an odd number of bars
* ITRA_TYPE = 4 : 2-component relativistic transformation for MCSCF,
*                 all nonredundant integral lists
* ITRA_TYPE = 5 : 2-component relativistic transformation for MCSCF,
*                 Excluding integrals with an odd number of bars
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
*
*. output lists of integrals
      dimension I12O(*),I34O(*),I1234O(*)
      dimension ICCSMO(*),IRIO(MXBLKO)
      dimension ISPITRA(4,*),ISPITRA_RC(*)
*
*. size of spinor expansions
      dimension LSPINOR(N_SPINR_EXP)
*
#include "clunit.inc"
*
777   format(1X,4I3)
888   format(1X,2I3)
999   format(1X,8I3)
700   format(6I3)
*
*
      WRITE(6,*) ' ***************************** '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' *     Welcome to TWELTRA    * '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' ***************************** '
*
*  Test
*
      if (DOUGRP.eq.4) then
         write(*,*) 'Point group mult. matrix'
         do I=1,NIRR_PN,1
            write(6,'(10X,4I3)') (IPGMULT(I,J),J=1,NIRR_PN,1)
         end do
         write(*,*) 'Double group mult. matrix'
         do I=1,NIRR_DG,1
            write(6,'(10X,8I3)') (IDBGMULT(I,J),J=1,NIRR_DG,1)
         end do
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
         do I=1,NIRR_PN,1
            write(6,'(10X,2I3)') (IPGMULT(I,J),J=1,NIRR_PN,1)
         end do
         write(*,*) 'Double group mult. matrix'
         do I=1,NIRR_DG,1
            write(6,'(10X,4I3)') (IDBGMULT(I,J),J=1,NIRR_DG,1)
         end do
      end if
*
*. Normal list of input integrals
      I12 = 1
      I34 = 1
      I1234 = 1
      ICCSM = 1
      IOPSMPN = 1
      IOPSMDG = 1
*
*  Read input from file 40:
      read(LU40,'(A)')
      read(LU40,*) NLISTO
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,700) I12O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) I34O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) I1234O(I)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) ICCSMO(I)
      end do
*
      read(LU40,'(A)')
      do I=1,4,1
         read(LU40,700) (ISPITRA(I,J),J=1,NLISTO,1)
      end do
*
      read(LU40,'(A)')
      do I=1,NLISTO,1
         read(LU40,*) ISPITRA_RC(I)
      end do
*
      write(6,*) ' tweltra: ITRA_TYPE = ' ,ITRA_TYPE
      if (ITRA_TYPE.eq.1) then
*  Good old non-relativistic integral transformation
*. A normal integral transformation, spinors have just
* pn symmetry, and symmetry effects of spin functions
* are neglected
         NSPINOR = 1
C         NLISTO = 1
C         I12O(1) = 1
C         I34O(1) = 1
C         I1234O(1) = 1
C         ICCSMO(1) = 1
C         ISPITRA(1,1) = 1
C         ISPITRA(2,1) = 1
C         ISPITRA(3,1) = 1
C         ISPITRA(4,1) = 1
C         ISPITRA_RC(1) = 1
      else if (ITRA_TYPE.eq.2) then
         NSPINOR = 2
*  Relativistic transformation to double group spinors,
*  CI with all classes of integrals included.
*  Unbarred trans is 1, bar trans is 2.
*
* First list : (K L_ | M N_) where _ = bar
* Second list : (K L_ | M N)
* Third list : (K L | M N)
* Fourth list : (K_ L | M N_)
* Fifth list : (K_ L | M N)
* Sixth list : (K_ L | M_ N)
*
      else if(ITRA_TYPE.EQ.3) then
*. relativistic transformation to double group spinors,
*. CI with integrals including an even number of bars
         NSPINOR = 2
*
* First list : (K L_ | M N_) Where _ = bar
* second list : (K L | M N)
* Third list : (K_ L | M N_)
* Fourth list : (K_ L | M_ N)
      else if (ITRA_TYPE.eq.4) then
         write(6,*) ' Add types for MCSCF transformation'
      ELSE IF (ITRA_TYPE.EQ.5) THEN
         WRITE(6,*) ' Add types for MCSCF transformation'
      END IF
*
*. Relative size of spinor expansions.
*  As an example, in the case of C2(*) there are two different spinor
*  expansions, one for alpha and one for beta spinors. The size of the
*  expansions of two particular spatial symmetries might be different, so that
*  there is an overall difference of the expansion sizes of two indices in
*  question. LSPINOR assumes the value 1 for a given expansion and another
*  1 if the expansion compared to it has the same size, a value of 2 if the
*  second expansion is smaller than the first. The expansion compared to
*  each other are characterized by the value of ISPITRA(index,list), which
*  in this case is either 1 or two, i.e. there are 2 expansions, and in this
*  case they are of equal size.
*
      DO JSPINOR = 1,NSPINOR,1
         LSPINOR(JSPINOR) = 1
      END DO

      CALL SPITRA2(I12O,I34O,I1234O,
     &             ICCSMO,ISPITRA,ISPITRA_RC,
     &             LSPINOR,IRIO,ICCSM,
     &             NLISTO,I1234,I12,I34,IOPSMPN,IOPSMDG)
*
      WRITE(6,*) ' ***************************** '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' *     TWELTRA says goodbye  * '
      WRITE(6,*) ' *                           * '
      WRITE(6,*) ' ***************************** '
*
      return
      END
