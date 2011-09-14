!
!    DIRAC - KRCC memory pointers
!
      INTEGER*8
! ctcc.inc start
! Hamiltonian
     &             KLSOBEX,KLLSOBEX,KLIBSOBEX,
     &             KLIHIND,KSIGNNHX,KLABEXTP,
     &             KNH22,KNH21,KNH12,KNH20,
     &             KNH11,KNH02,KNH10,KNH01,KNH00,
! Cluster operator
     &             KLSOBEX_CC,KLIBSOBEX_CC,KLLSOBEX_CC,
     &             KT_IDXS,KT_IDXF, 
! Intermediates
     &             KINTM10 ,KINTM01 ,KINTM11,
     &             KINTM02 ,KINTM20 ,KINTM12,
     &             KINT01_IDXS,KINT01_IDXF,KINT10_IDXS,KINT10_IDXF,
     &             KINT11_IDXS,KINT11_IDXF,KINT02_IDXS,KINT02_IDXF,
     &             KINT20_IDXS,KINT20_IDXF,KINT12_IDXS,KINT12_IDXF,
! Help arrays for intermediates
     &             KLAOBEX_CC,KLBOBEX_CC,
     &             KLAOBEX_KRCC,KLBOBEX_KRCC,
! Contraction mappings
     &             KLLSOBEX_M01,KLIBSOBEX_M01,KM01,
     &             KLLSOBEX_M10,KLIBSOBEX_M10,KM10,
     &             KLLSOBEX_M11,KLIBSOBEX_M11,KM11,
     &             KLLSOBEX_M02,KLIBSOBEX_M02,KM02,
     &             KLLSOBEX_M12,KLIBSOBEX_M12,KM12,
! Interface integrals
     &             KFCACM,KH2AC,KINC_HX,
! Others 
     &             KLSOX_TO_OX,
! Only used once in strings
     &             KLOP_REO,KSIOPREO,
! from here only in program.F
     &             KOBEX_TP,KLCOBEX_TP,KLAOBEX_TP,
     &             KIBSOX_FOR_OX, KNSOX_FOR_OX,KISOX_FOR_OX,
     &             KLSPOBEX_AB, KIBSOX_FOR_OCCLS,KNSOX_FOR_OCCLS,
     &             KISOX_FOR_OCCLS,
! ctcc.inc end
! strbas.inc start
     & KSTINF, KSTINE, KOCSTR, KNSTSO,
     & KISTSO, KSTSTM, KZ,
     & KSTREO, KNSTSGP,KISTSGP ,
     & KIOCLS, KSPGPAN, KSPGPCR,
     & KNSTSGP2, KISTSGP2,
     & KOCSTR2, KSTREO2,
     & KSTSTM2, KNSTSO2, KISTSO2,
! strbas.inc end
! glbbas.inc start
     & KRHO1, KVEC3, KFI, KLOOK,
     & KRHO2, KLOCCLS, KKR,
     & KINVREO, KT_CC, KBOSYM, KISREO,
     & KFDIA_AL, KFDIA_BE, KFDIA_AL_IM, KFDIA_BE_IM,
     & KCC1, KCC2, KCC3, KCC4, KCC5, KDIA,
     & KSBSPJA, KSRHO1, KBMAT, KSCDIIS,
! glbbas.inc end
! cc_scr2.inc start
     &KDUM, KFREEL, 
     & KLZ, KLZSCR, KLSTOCC1,
     & KLSTOCC2, KLSTOCC3, KLSTOCC4, KLSTREO,
     & KIX1_CA, KSX1_CA, KIX1_CB, KSX1_CB,
     & KIX1_AA, KSX1_AA, KIX1_AB, KSX1_AB,
     & KIX2_CA, KSX2_CA, KIX2_CB, KSX2_CB,
     & KIX2_AA, KSX2_AA, KIX2_AB, KSX2_AB,
     & KLTSCR1, KLTSCR2, KLOPSCR,
     & KLSMD1, KLSMEX, KLSMK1, KLSML1,
     & KLNMD1, KLNMEX, KLNMK1, KLNML1,
     & KLOCK1, KLOCL1,
! cc_scr2.inc end
! Not on other common block
! Only used for integrals
     & KISTR_AA, KISTR_AB,
     &KISTR_CA, KISTR_CB,

     &KC2,  
     &KLBLTP,  
     &KLCBLTP, KLCIOIO,  
     &KLIOIO,  
     &KLSCR, KLSCR2, 
     &KLSOBEX_TEMP, KLSTR1_OCC, 
     &KLSTR2_OCC, KLSTR3_OCC, KLSTR4_OCC, 
     &KNCPMT_OX, KNCPMT_SOX,  
     &KPCBLTP, KPCI1BT, KPCIBT, KPCLBT, KPCLEBT, 
     & KSB,
     &KVEC1, KVEC2, KWCPMT_OX, 
     &KWCPMT_SOX, 
     &KOCSTR_ADD,
     &KMUB
