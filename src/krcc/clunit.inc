* CLUNIT
*
*   Assigning logical units
*
      character*6 DIRACTWO
      character*7 DIRACONE
      COMMON/KRCC_CLUNIT/LUIN,LUOUT,
     &              LUDIA,LUC,LUHC,
     &              LUSC1,LUSC2,LUSC3,LUSC4,LUSC5,LUSC6,
     &              LUSC41,
     &              LUSC81,LUSC82,LU_ARDUCCA,
*   Added for MRCC
     &              LU_CCVEC,LU_CCVECT,LU_CCAMP,
     &              LU_CCVECF,LU_CCVECL,LU_CCVECFL,LU_CCEXC_OP,
*   Intermediates to lrcc
     &              LU_M10,LU_M02,LU_M11,LU_M01,LU_M12,LU_M11_TO_M01
      COMMON/KRCC_CLNAME/DIRACONE,DIRACTWO
*
*   Logical units are: See diskun_krcc for information
