* CLUNIT
*
*   Assigning logical units
*
      integer   LU2,LU5,LU9,LU19,LU40,LU69,LU70,LU71,
     &          LUSGI1,LUSGI2,LUSGI3,LUSGI4
      character*6 FNONEAO,FNIPH,FNONEMO,DIRACTWO
      character*7 DIRACONE
      COMMON/KRCC_CLUNIT/LUIN,LUOUT,LUCIVI,LU1INT,LU1SO,LU2INT,
     &              LUPRP,
     &              LUDIA,LUC,LUHC,
     &              LUSC1,LUSC2,LUSC3,
     &              LUSC34,LUSC35,LUSC36,LUSC37,LUSC38,LUSC39,
     &              LUSC40,LUSC41,LUSC61,LUSC62,
     &              LUCIVO,LUORB,LUMOIN,LUMOUTR,LUMOUTI,
     &              LU2,LU5,LU9,LU19,LU40,LU69,LU70,LU71,
     &              LUSGI1,LUSGI2,LUSGI3,LUSGI4,
     &              LUSC81,LUSC82,
     &              LUSC92,LUSC93,LUSC94,LUSC95,LUSC96,
     &              LUSC97,LUSC98,LUSC99,
*   Added for MRCC
     &              LU_CCVEC,LU_CCVECT,LU_CCAMP,
     &              LU_CCVECF,LU_CCVECL,LU_CCVECFL,LU_CC_FROM_CI,
     &              LU_CCEXC_OP
      COMMON/KRCC_CLNAME/FNONEAO,
     &              FNIPH,FNONEMO,DIRACONE,DIRACTWO
*
*   Logical units are: See diskun.f for information
