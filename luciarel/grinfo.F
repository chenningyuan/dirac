*  This subroutine generates all group info needed for
*   - integral transformation
*   - CI
*   - Spinor rotation
*
      subroutine grinfo(PNTGRP)
      implicit real*8 (A-H,O-Z)
*
      integer PNTGRP
*
#include "mxpdim.inc"
#include "symm.inc"
*
      if (PNTGRP.eq.9) DOUGRP = 4
      if (PNTGRP.eq.8) DOUGRP = 5
      if (PNTGRP.eq.7) DOUGRP = 6
      if (PNTGRP.eq.6) DOUGRP = 7
      if (PNTGRP.eq.5) DOUGRP = 8
*
      if (DOUGRP.eq.1) then
         NIRR_DG = 10
         NIRR_PN = 8
         IHASDG = 4
         PNTGRP = 1
      else if (DOUGRP.eq.2) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.3) then
         NIRR_DG = 5
         NIRR_PN = 4
         IHASDG = 5
         PNTGRP = 0
         write(6,*) 'Point group not defined in LUCIA '
         write(6,*) 'watch out!'
      else if (DOUGRP.eq.4) then
         NIRR_DG = 8
         NIRR_PN = 4
         IHASDG = DOUGRP
         PNTGRP = 9
      else if (DOUGRP.eq.5.or.DOUGRP.eq.6) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         if (DOUGRP.eq.5) PNTGRP = 8
         if (DOUGRP.eq.6) PNTGRP = 7
      else if (DOUGRP.eq.7) then
         NIRR_DG = 4
         NIRR_PN = 2
         IHASDG = DOUGRP
         PNTGRP = 6
      else if (DOUGRP.eq.8) then
         NIRR_DG = 2
         NIRR_PN = 1
         IHASDG = DOUGRP
         PNTGRP = 5
      else
         write(*,*) 'Double group out of range!'
         write(*,*) 'Double group is  ',DOUGRP
         write(*,*) 'Picasso will put down his brush.'
         Call Abend()
      end if
*
*  Set up info for highest abelian double subgroup
*
      if (DOUGRP.ne.IHASDG) then
         if (DOUGRP.eq.1) then
            NSDGIRR = 8
            NSPGIRR = 4
         else if (DOUGRP.eq.2.or.DOUGRP.eq.3) then
            NSDGIRR = 4
            NSPGIRR = 2
         end if
      end if
*
      MXPDBGDIM = 8
      MXPDBGDIM2 = MXPDBGDIM ** 2
*
      IZERO = 0
      CALL ISETVC(IDBGMULT,IZERO,MXPDBGDIM2)
      CALL ISETVC(IMOSF_SP,IZERO,MXPDBGDIM)
      CALL ISETVC(ISPSF_MO,IZERO,MXPDBGDIM)
      IF(PNTGRP.EQ.5) THEN
*. C1 point group
        DBGSYM(1:6) ='C1    '
*
        IPGMULT(1,1) = 1
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
*. All spinors have symmetry 2, all orbitals symmetry 1
        IMOSF_SP(1,1) = 2
        IMOSF_SP(1,2) = 2
*
        ISPSF_MO(2,1) = 1
        ISPSF_MO(2,2) = 1
*. Spinors and their cc have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
      ELSE IF(PNTGRP.EQ.6 ) THEN
*. Ci pointgroup
        DBGSYM(1:6) ='Ci    '
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 1
        IDBGMULT(3,4) = 2
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 2
        IDBGMULT(4,4) = 1
*. MO have symmetries 1,2 spinors symmetries 3,4 alpha,beta have symm 3
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 3
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 4
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 1
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 2
*. Functions and adjoints have identical symmetries
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
      ELSE IF(PNTGRP.EQ.7.OR.PNTGRP.EQ.8) THEN
*. CS and C2 point groups
        IF(PNTGRP.EQ.7) THEN
          DBGSYM(1:6) ='Cs    '
        ELSE
          DBGSYM(1:6) ='C2    '
        END IF
*
        IPGMULT(1,1) = 1
        IPGMULT(1,2) = 2
        IPGMULT(2,1) = 2
        IPGMULT(2,2) = 1
*
        IDBGMULT(1,1) = 1
        IDBGMULT(1,2) = 2
        IDBGMULT(1,3) = 3
        IDBGMULT(1,4) = 4
        IDBGMULT(2,1) = 2
        IDBGMULT(2,2) = 1
        IDBGMULT(2,3) = 4
        IDBGMULT(2,4) = 3
        IDBGMULT(3,1) = 3
        IDBGMULT(3,2) = 4
        IDBGMULT(3,3) = 2
        IDBGMULT(3,4) = 1
        IDBGMULT(4,1) = 4
        IDBGMULT(4,2) = 3
        IDBGMULT(4,3) = 1
        IDBGMULT(4,4) = 2
*. Orbitals have sym 1,2, spinors 3,4, alpha 3, beta 4
        IMOSF_SP(1,1) = 3
        IMOSF_SP(1,2) = 4
        IMOSF_SP(2,1) = 4
        IMOSF_SP(2,2) = 3
*
        ISPSF_MO(3,1) = 1
        ISPSF_MO(3,2) = 2
        ISPSF_MO(4,1) = 2
        ISPSF_MO(4,2) = 1
*. Adjoints
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 4
        IADJSYM(4) = 3
* inverse elements
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 4
        INVELM(4) = 3
      ELSE IF(PNTGRP.EQ.9) THEN
*. C2H
        DBGSYM(1:6) ='C2h   '
*
        IADJSYM(1) = 1
        IADJSYM(2) = 2
        IADJSYM(3) = 3
        IADJSYM(4) = 4
        IADJSYM(5) = 6
        IADJSYM(6) = 5
        IADJSYM(7) = 8
        IADJSYM(8) = 7
*
        INVELM(1) = 1
        INVELM(2) = 2
        INVELM(3) = 3
        INVELM(4) = 4
        INVELM(5) = 6
        INVELM(6) = 5
        INVELM(7) = 8
        INVELM(8) = 7
*
        ISPSF_MO(5,1) = 1
        ISPSF_MO(6,2) = 1
        ISPSF_MO(6,1) = 2
        ISPSF_MO(5,2) = 2
        ISPSF_MO(7,1) = 3
        ISPSF_MO(8,2) = 3
        ISPSF_MO(8,1) = 4
        ISPSF_MO(7,2) = 4
*
        IMOSF_SP(1,1) = 5
        IMOSF_SP(1,2) = 6
        IMOSF_SP(2,1) = 6
        IMOSF_SP(2,2) = 5
        IMOSF_SP(3,1) = 7
        IMOSF_SP(3,2) = 8
        IMOSF_SP(4,1) = 8
        IMOSF_SP(4,2) = 7
*
        do JJ=1,8,1
           do II=1,8,1
              if (II.le.4.and.JJ.le.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 1
                 else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.5) then
                    IDBGMULT(II,JJ) = 4
                 else
                    IDBGMULT(II,JJ) = 3
                 end if
              else if (II.gt.4.and.JJ.gt.4) then
                 if (II.eq.JJ) then
                    IDBGMULT(II,JJ) = 2
                 else if ((II+JJ).eq.13) then
                    IDBGMULT(II,JJ) = 3
                 else if ((II+JJ).eq.11.or.(II+JJ).eq.15) then
                    IDBGMULT(II,JJ) = 1
                 else
                    IDBGMULT(II,JJ) = 4
                 end if
              else
                 if ((II+JJ).eq.9) then
                    IDBGMULT(II,JJ) = 8
                 else if ((II+JJ).eq.7.or.(II+JJ).eq.11) then
                    IDBGMULT(II,JJ) = 6
                 else if ((II+JJ).eq.6.or.(II+JJ).eq.12) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.6.and.JJ.eq.2) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.2.and.JJ.eq.6) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.3.and.JJ.eq.7) then
                    IDBGMULT(II,JJ) = 5
                 else if (II.eq.7.and.JJ.eq.3) then
                    IDBGMULT(II,JJ) = 5
                 else
                    IDBGMULT(II,JJ) = 7
                 end if
              end if
           end do
        end do
        do JJ=1,4,1
           do II=1,4,1
              if (II.eq.JJ) then
                 IPGMULT(II,JJ) = 1
              else if ((II+JJ).eq.3.or.(II+JJ).eq.7) then
                 IPGMULT(II,JJ) = 2
              else if ((II+JJ).eq.5) then
                 IPGMULT(II,JJ) = 4
              else
                 IPGMULT(II,JJ) = 3
              end if
           end do
        end do
        continue
      END IF
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Symmetry Information from GET_DBGSYM'
        WRITE(6,*) ' ===================================='
        WRITE(6,*)
        WRITE(6,'(A,A)') ' Double group : ', DBGSYM
        WRITE(6,*)
        WRITE(6,'(A,I2)') ' Number of irreps ',NIRR_DG
        WRITE(6,*)
        WRITE(6,*) ' Multiplication table '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        DO IREP = 1, NIRR_DG
          WRITE(6,'(10X,8I3)') (IDBGMULT(IREP,JREP),JREP=1,NIRR_DG)
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' orbital X spin => spinor '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO IOBSM = 1, NIRR_DG/2
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    IOBSM, ISF, IMOSF_SP(IOBSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Spinor  X spin => orbital'
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        DO ISPSM =  NIRR_DG/2+1, NIRR_DG
         DO ISF = 1, 2
          WRITE(6,'(3X,I3,3X,I3,5X,I3)')
     &    ISPSM, ISF, ISPSF_MO(ISPSM,ISF)
         END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Adjoint symmetry array '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,'(10X,8I3)') (IADJSYM(IREP),IREP=1,NIRR_DG)
*
        WRITE(6,*) ' Inverse elements '
        WRITE(6,*) ' ================='
        WRITE(6,'(10X,8I3)') (INVELM(IREP),IREP=1,NIRR_DG)
      END IF
*
      return
      end
