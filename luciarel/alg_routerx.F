      SUBROUTINE ALG_ROUTERX(ILOC,IROC,NOP,IOP_TYP,IOP_AC,IOP_REO,
     &                      SIGN)
*
* Decide route for calculating <ILOC! operator string !ROC>
* with smallest amount of operations/storage
*
* The operator string contains NOP elementary operators
* defined by
*
* IOP_TYP : Orbital type of each operator
* IOP_AC  : Creation/annihilation operator
*
*. Output :
* ==========
*
* IOP_REO : New to old order
* sign    : sign of full rank operator
*
* Method : use IPHGAS , and move creation of holes
*                       and annihilation of particles to the left
*
* Jeppe Olsen, October 1997
* version  of : dec 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
*. Specific input, I*OC : Contains number of electrons in each space
*  (not used !!!)
       INTEGER ILOC(NGAS),IROC(NGAS)
       INTEGER IOP_TYP(*), IOP_AC(*)
*. Output, new to old order :
       INTEGER IOP_REO(*)
*
*. The loops are assumed to be of the resolution type, so
*. we will insert a resolution of the identity after NOP/2
*  operators(we assume that NOP is even ) :
*
* <ILOC! IOPL !K><K! IOPR ! IROC >
*=<ILOC! IOPL !K><IROC! (IOPR)+ ! K >
*
*
* and the question is now how to order the operators.
* The operation count can in general be mininized by
*. According to above : Should operators stand left or
*                       right of the resolution K
*
      NTEST = 00
*. First time around : not very general
      NMOVE = 0
      IF(NOP.EQ.2) THEN
*
* One-Electron operator
*
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.1) THEN
* a+1 a 2 : Move around only if both spaces are hole spaces
          NMOVE = 0
          IF(IPHGAS(IOP_TYP(1)).EQ.2) NMOVE = 1
          IF(IPHGAS(IOP_TYP(2)).EQ.2) NMOVE = NMOVE + 1
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) '  IOP_TYP(1), IOP_TYP(2), NMOVE ',
     &                    IOP_TYP(1), IOP_TYP(2), NMOVE
          END IF
        ELSE
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*
        IF(NMOVE.EQ.2) THEN
*. Move the operators around in expansion
          IOP_REO(1) = 2
          IOP_REO(2) = 1
          SIGN = -1.0D0
        ELSE
*. No reorganization adviced
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          SIGN = 1.0D0
        END IF
      ELSE IF(NOP.EQ.4) THEN
*
* Two-electron operator
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.2.AND.
     &     IOP_AC(3).EQ.1.AND.IOP_AC(4).EQ.1) THEN
* a+1 a+2 a3 a4
          IF(IPHGAS(IOP_TYP(1)).EQ.2) THEN
            IMOVE1 = 1
          ELSE
            IMOVE1 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(2)).EQ.2) THEN
            IMOVE2 = 1
          ELSE
            IMOVE2 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(3)).EQ.2) THEN
            IMOVE3 = 1
          ELSE
            IMOVE3 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(4)).EQ.2) THEN
            IMOVE4 = 1
          ELSE
            IMOVE4 = 0
          END IF
        ELSE
*. Not yet implemented
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*. Number of left and right operators that would like to be moved
        IMOVEL = IMOVE1 + IMOVE2
        IMOVER = IMOVE3 + IMOVE4
        IMOVE = MIN(IMOVEL,IMOVER)
        IF(IMOVE.EQ.2) THEN
*. Well, we are all moving
          IOP_REO(1) = 3
          IOP_REO(2) = 4
          IOP_REO(3) = 1
          IOP_REO(4) = 2
          SIGN = 1.0D0
        ELSE IF(IMOVE.EQ.0) THEN
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
        ELSE IF(IMOVEL.EQ.1.AND.IMOVER.EQ.1) THEN
* a+ a a+ a
* 1 : Position the two operators to be moved as operators 2, 3
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
          IF(IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -SIGN
          END IF
          IF(IMOVE4.EQ.1) THEN
           IOP_REO(3) = 4
           IOP_REO(4) = 3
           SIGN = - SIGN
          END IF
* 2 : and interchange operators 2 and 3
          ISAVE = IOP_REO(2)
          IOP_REO(2) = IOP_REO(3)
          IOP_REO(3) = ISAVE
          SIGN = - SIGN
        ELSE IF (IMOVEL.EQ.2.AND.IMOVER.EQ.1) THEN
* Final operator is a a+ a+ a
          IF(IMOVE4.EQ.1) THEN
            IOP_REO(3) = 4
            IOP_REO(4) = 3
            SIGN = -1.0D0
          ELSE
            IOP_REO(3) = 3
            IOP_REO(4) = 4
            SIGN = 1.0D0
          END IF
          IOP_REO(1) = IOP_REO(3)
          IOP_REO(2) = 1
          IOP_REO(3) = 2
          SIGN = SIGN
        ELSE IF( IMOVEL.EQ.1 .AND. IMOVER .EQ. 2 ) THEN
* Final operator is a+ a a a+
          IF( IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -1.0D0
          ELSE
            IOP_REO(1) = 1
            IOP_REO(2) = 2
            SIGN = 1.0D0
          END IF
          IOP_REO(4) = IOP_REO(2)
          IOP_REO(2) = 3
          IOP_REO(3) = 4
          SIGN = SIGN
        END IF
*       ^ End of number of pairs to be moved around
      ELSE
*. Not one- or two- electron operator
        WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
        WRITE(6,*) ' NOP  = ', NOP
        WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
        Call Abend2( ' ALG_ROUTER : Path not implemented ' )
      END IF
*

      IF(NOP.EQ.4.AND.NTEST.GE.100) THEN
        WRITE(6,*) ' Information from ALG_ROUTER '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
     &  ' Input : anni- or crea- operators + orbital types'
        CALL IWRTMA(IOP_AC,1,NOP,1,NOP)
        CALL IWRTMA(IOP_TYP,1,NOP,1,NOP)
        WRITE(6,*)
        WRITE(6,*) ' Suggested output order '
        CALL IWRTMA(IOP_REO,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
