!
!     KR-CI common block /CRUN/
!      
!     this common block contains all information that is 
!     needed to determine KR-CI runtype, convergence information 
!     as well as some symmetry variables. 
!
!     THRES_G:   CI gradient convergence threshold
!     THRES_E:   CI energy change convergence threshold
!
!     ISYM_T:    Hamiltonian operator symmetry
!     MXCIV_CI:  max. dimension of CI subspace in Davidson diagonalizer
!     LBLOCK:    max. length of vector block to be allocated for
!                Davidson diagonalizer - dynamically determined
!     MAXIT:     max. number of Davidson iterations
!
!  Not used in this version:
!     E_THRE:    individual second order selection of included csf in trial vector
!     C_THRE:    individual first order selection of included csf in trial vector
!
      logical KRASYM
      COMMON/CRUNR/THRES_G,THRES_E,XLAMBDA,E_THRE,C_THRE,E_CONV,C_CONV,
     &            MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV_CI,
     &            ICISTR,IDIAG,NOINT,MXINKA,ICJKAIB,
     &            IPERT,NPERT,IIDUM,NSEQCI(10),ISEQCI(10,10),
     &            INIDEG,LCSBLK,NPSSPC,
     &            ICLSSEL,IDENSI,
     &            IH0ROOT,ISKIPEI,MK2REF_CI,MK2DEL_CI,
     &            ISPINFREE,NPROP,NATITER,ISYM_T,
     &            IDCOMH,L2BLOCK,LBLOCK,KRASYM
      character*72 ENVIRO
      CHARACTER*6 PROPER,SIGDEN_ROUTE
      COMMON/CRUNCR/ENVIRO,SIGDEN_ROUTE,PROPER(20)
