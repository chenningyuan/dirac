      SUBROUTINE WRTHD_KRCC(LUOUT)
*
* Introduce progam and author
*
      WRITE(6,'(1H ,8X,A)')
     &'****************************************************************'
      WRITE(6,'(1H ,8X,A)')
     &'*                                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'*           Welcome to KRCC, a CC program                      *'
      WRITE(6,'(1H ,8X,A)')
     &'*                                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'* Written by Lasse K. Soerensen, University of Aarhus, Denmark *'
      WRITE(6,'(1H ,8X,A)')
     &'* Version of 2010                                              *'
      WRITE(6,'(1H ,8X,A)')
     &'****************************************************************'
*
      WRITE(6,*)
      WRITE(6,'(1H ,8X,A)')
     &' In case of trouble please contact : '
      WRITE(6,'(1H ,12X,2A)')
     &' Lasse K. Soerensen, Dept. of Chemistry, University of Aarhus,',
     &' Aarhus, DK-8000 Denmark '
      WRITE(6,'(1H ,12X,A)')
     &' Telephone :+45 '
      WRITE(6,'(1H ,12X,A)')
     &' E-mail : lks@chem.au.dk   '
      WRITE(6,*)
      RETURN
      END
*
      SUBROUTINE WRT_TCC_BLK_CC_REL(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
*
*. Write TCC block containing all symmetries, total sym of T is ITCC_SM
*
* Jeppe Olsen, summer of 99
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "multd2h.inc"
#include "symm.inc"
*. Input
       DIMENSION TCC(*)
       INTEGER NCA(*), NCB(*), NAA(*), NAB(*)
*
       WRITE(6,*) ' Block of Coupled cluster vector '
       WRITE(6,*) ' ================================'
       IOFF = 1
       DO ISM_C = 1, NSMST
         ISM_A = IDBGMULT(ITCC_SM,INVELM(ISM_C))
         DO ISM_CA = 1, NSMST
           ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
           DO ISM_AA = 1, NSMST
             ISM_AB =  IDBGMULT(ISM_A,INVELM(ISM_AA))
C             WRITE(6,*) ' ISM_AB = ', ISM_AB
             LCA = NCA(ISM_CA)
             LCB = NCB(ISM_CB)
             LAA = NAA(ISM_AA)
             LAB = NAB(ISM_AB)
             LENGTH = LCA*LCB*LAA*LAB
C             print*,'LENGTH',LENGTH
*
             IF(LENGTH.NE.0) THEN
               WRITE(6,'(A,4I4)') ' Sym of CA, CB, AA, AB ',
     &         ISM_CA, ISM_CB, ISM_AA, ISM_AB
               CALL WRTMAT(TCC(IOFF),1,LENGTH,1,LENGTH)
             END IF
             IOFF = IOFF + LENGTH
           END DO
         END DO
       END DO
*
       RETURN
       END
