      SUBROUTINE WRT_SPOX_TP(IEX_TP,NEX_TP)
*
* Print types of spin-orbital excitations
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER IEX_TP(4*NGAS,NEX_TP)
*
      WRITE(6,*)
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*) ' Information about spinorbital excitations '
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*)
*
      DO JEX_TP = 1, NEX_TP
        WRITE(6,*)
        WRITE(6,*) ' Included spinorbitalexcitation ', JEX_TP
        WRITE(6,'(A,16I4)')
     &  ' Creation of alpha     :',
     &  (IEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Creation of beta      :',
     &  (IEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of alpha :',
     &  (IEX_TP(I+2*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)')
     &  ' Annihilation of beta  :',
     &  (IEX_TP(I+3*NGAS,JEX_TP),I=1,NGAS)
      END DO
*
      RETURN
      END
