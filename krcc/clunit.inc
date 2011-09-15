* CLUNIT
*
*   Assigning logical units
*
      character*6 DIRACTWO
      character*7 DIRACONE
      COMMON/KRCC_CLUNIT/LUIN,LUOUT,
     &              LUDIA,LUC,LUHC,
     &              LUSC1,LUSC2,LUSC3,
     &              LUSC41,
     &              LUSC81,LUSC82,
*   Added for MRCC
     &              LU_CCVEC,LU_CCVECT,LU_CCAMP,
     &              LU_CCVECF,LU_CCVECL,LU_CCVECFL,LU_CCEXC_OP
      COMMON/KRCC_CLNAME/DIRACONE,DIRACTWO
*
*   Logical units are: See diskun.f for information
