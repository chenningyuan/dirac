*/CRUN/ : Defines RUN 
* NOTE : Differs from OCT9 version of CRUN
      COMMON/KRCC_CRUN/
     &            THRES_E,XLAMBDA,E_THRE,C_THRE,E_CONV,C_CONV,E0READ,
     &            RESP_W(20),RESP_OP(2,20),AVE_OP(20),CCCONV    ! added by lasse for gen_cc

* !!! NOTE : only integers in next section !!! ( real*8 above, logicals at bottom)
      COMMON/KRCC_CRUNI/
     &            MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,
     &            ICISTR,NOCSF,IDIAG,NOINT,IDMPIN,MXINKA,MXINKA_CC,
     &            ICJKAIB,MXCCV,NTOTORB,
     &            INIREF,IRESTRF,MOCAA,MOCAB,IPERT,NPERT,
     &            IAPRREF,IAPRZER,IIDUM,
     &            NSEQCI(MXCMBSPCK),ISEQCI(MXPICI,MXPICI),
* added by lasse for sequen
     &            IE0AVEX,IC1DSC,
     &            IH0SPC,NPTSPC,IPTSPC(MXPPTSPC),
     &            IOCPTSPC(2,MXPNGAS,MXPPTSPC),IH0INSPC(MXPPTSPC),
     &            IRFROOT,NH0EXSPC,IH0EXSPC(MXPPTSPC),INIDEG,
     &            LCSBLK,NOMOFL,IFINMO,NPSSPC,
     &            NPSSH(MXPIRR,MXPNGAS),
     &            ICLSSEL,IDENSI,IPTEKT,
     &            NPTEKT,LPTEKT(20,2),IH0ROOT,IRST2,ISKIPEI,IPAD,
     &            MK2REF,MK2DEL,IRELACI,IDBGCI,
     &            ISPINFREE,ISO1,ISO2,
     &            IXYZSYM(3),NPROP,
     &            ITRANS,ISPEXP_INI,
     &            NATITER,INATIT,ISYM_T,
     &            IDCOMH,IFOCKEX,NORB_F(MXPIRR),
     &            ISIMTRH,
* more
     &            MXSPOX,IAAEXC_TYP,
* added by lasse for gen_cc
     &            ITRAPRP,NEXCSTATE,IEXCSYM,IRESPONS,NRESP,
     &            MAXORD_OP(2,20),
     &            MXITLE,  N_AVE_OP, IROOTHOMING,
     &            IUSE_PH, IADVICE,ITRACI,ITRACI_CR,ITRACI_CN,
     &            IUSE_PA,IPTFOCK,NSXE,NDXEE,N_CC_AMP,
     &            I_RESTRT_CC,LAST_CC_SPC,LAST_CC_RUN,ITRA_FI,
     &            ITRA_IN,MULSPC,LPAT,IPAT(100),IFMULSPC, IRELAX,
     &            ICC_CUR, NCNV_RT, IDOQD,
     &            I_RE_MS2_SPACE,I_RE_MS2_VALUE, IPRECOND,
     &            ISIMSYM, IUSE_HW, 
     &            IUSEH0P, I_DO_MCSCF, MX_EXC_LEVEL,
     &            IMAXKRFLIP,IMUBMAX,
     &            I_DO_LZ2, ICCSOLVE,I_DO_CCN, I_DO_SBSPJA,
     &            IDUM1, N_HOLE_ORBSPACEX, I_DO_CC3,
     &            I_DO_CI_TO_CC,LCCBD12,LCCB,
     &            I_DO_CC_EXC_E, NEXC_PER_SYM(8),IRES_EXC,
     &            IFRZ_CC, NFRZ_CC,IFRZ_CC_AR(100),
     &            I_DO_CC_EXP, I_DO_GENTRD,
     &            ISM_AT, IMS2_AT, NSTATE_AT, I_DO_MASK_CC,
     &            MASK_SD(100,2),MSK_AEL,MSK_BEL, NOAAEX, MXCIV_ORIG,
     &            N_TDL_MAX,IUSE_TR, ISPIN_RESTRICTED, IGENTRD,
     &            IGST_OCC(MXPNGAS,2),IGST_SM,IGST_MS2, 
     &            I_DO_REO_ORB,NSWITCH,IREO_ORB(3,MXPORB), I_DO_CC,
     &            ICEXC_RANK_MIN,ICEXC_RANK_MAX,ICEXC_INT,
     &            I_DO_CMPCCI, I_DO_CC_TO_CI,
     &            IRECOM,    ! for real/complex algebra
     &            I_EXPT_REF ! to control ISIMSYM_LOC
* !!! and now all logicals:
      logical KRASYM
      COMMON/KRCC_CRUNL/
     &            KRASYM

      character*72 ENVIRO
      CHARACTER*6 PROPER,SIGDEN_ROUTE,CSEQCI,CCFORM
      COMMON/KRCC_CRUNC/ENVIRO,SIGDEN_ROUTE,PROPER(20),CCFORM,
     &            CSEQCI(MXPICI,MXCMBSPCK)
