!
!     KR-CI common block /CRUN/
!      
!     this common block contains all information that is 
!     needed to determine KR-CI runtype, convergence information 
!     as well as some symmetry variables. 
!
!     ISYM_T:    Hamiltonian operator symmetry
!     MXCIV_CI:  max. dimension of CI subspace in Davidson diagonalizer
!     E_THRE:    convergence control wrt to energy (real algebra)
!     C_THRE:    convergence control wrt to energy (complex algebra)
!     LBLOCK:    max. length of vector block to be allocated for
!                Davidson diagonalizer - dynamically determined
!     MAXIT:     max. number of Davidson iterations
!
      logical KRASYM
      COMMON/CRUNR/THRES_E,XLAMBDA,E_THRE,C_THRE,E_CONV,C_CONV,E0READ,
     &            MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV_CI,
     &            ICISTR,NOCSF,IDIAG,NOINT,IDMPIN,MXINKA,ICJKAIB,
     &            INIREF,IRESTRF,MOCAA,MOCAB,IPERT,NPERT,
     &            IAPRREF,IAPRZER,IIDUM,
     &            NSEQCI(10),ISEQCI(10,10),IEXTKOP,
     &            IE0AVEX,IC1DSC,
     &            IH0SPC,NPTSPC,IPTSPC(MXPPTSPC),
     &            IOCPTSPC(2,MXPNGAS,MXPPTSPC),IH0INSPC(MXPPTSPC),
     &            IRFROOT,NH0EXSPC,IH0EXSPC(MXPPTSPC),INIDEG,
     &            LCSBLK,NOMOFL,IFINMO,NPSSPC,
     &            NPSSH(MXPIRR,MXPNGAS),
     &            ICLSSEL,IDENSI,IPTEKT,
     &            NPTEKT,LPTEKT(20,2),IH0ROOT,IRST2,ISKIPEI,IPAD,
     &            MK2REF_CI,MK2DEL_CI,IRELACI,IDBGCI,
     &            ISPINFREE,ISO1,ISO2,
     &            IXYZSYM(3),NPROP,
     &            ITRANS,ISPEXP_INI,
     &            NATITER,INATIT,ISYM_T,
     &            IDCOMH,IFOCKEX,NORB_F(MXPIRR),
     &            KRASYM,L2BLOCK,LBLOCK
      character*72 ENVIRO
      CHARACTER*6 PROPER,SIGDEN_ROUTE
      COMMON/CRUNCR/ENVIRO,SIGDEN_ROUTE,PROPER(20)
