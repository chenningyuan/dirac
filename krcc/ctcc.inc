      COMMON/KRCC_CTCC/
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
! Not pointers
! max something
     &             MX_ST_TSOSO_MX,MX_ST_TSOSO_BLK_MX,
     &             LEN_T_VEC,LEN_T_VEC_MX,LEN_T_VEC_CC,MX_SBSTR,NALLINT,
! Not really used but kept
! look at IDIM_TCC_KRCC
     &             N1ELINT,N2ELINT,ILTTWO,NHOL_SPC,
! Dimensions for H 
     &             NSPOBEX_TP,
     &             NH22,NH21,NH12,NH20,
     &             NH11,NH02,NH10,NH01,NH00,
! T
     &             NOBEX_TP,NSPOBEX_TPE,NSPOBEX_TP_CC,
! M
     &             NINTER10,NINTER01,NINTER11,
     &             NINTER02,NINTER20,NINTER12,
! For help to M
     &             LENA_KRCC,LENB_KRCC,
! Contraction mappings
     &             LEN_M01_VEC_CC,LEN_M10_VEC_CC,
     &             LEN_M11_VEC_CC,LEN_M02_VEC_CC,LEN_M12_VEC_CC,
! Other
     &             NAOBEX_TP,NBOBEX_TP
