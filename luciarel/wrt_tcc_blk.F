      SUBROUTINE WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
*
*. Write TCC block containing all symmetries, total sym of T is ITCC_SM
*
* Jeppe Olsen, summer of 99
*
#include "implicit.inc"
#include "multd2h.inc"
*. Input
       DIMENSION TCC(*)
       INTEGER NCA(*), NCB(*), NAA(*), NAB(*)
*
       WRITE(6,*) ' Block of Coupled cluster vector '
       WRITE(6,*) ' ================================'
       IOFF = 1
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ITCC_SM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
             ISM_AB =  MULTD2H(ISM_A,ISM_AA)
             WRITE(6,*) ' ISM_AB = ', ISM_AB
             LCA = NCA(ISM_CA)
             LCB = NCB(ISM_CB)
             LAA = NAA(ISM_AA)
             LAB = NAB(ISM_AB)
             LENGTH = LCA*LCB*LAA*LAB
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
