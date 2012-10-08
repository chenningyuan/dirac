*/CRUN/ : Defines RUN 
      COMMON/KRCC_CRUN/CCCONV

* !!! NOTE : only integers in next section !!! ( real*8 above, logicals at bottom)
      COMMON/KRCC_CRUNI/
     &            MAXIT,
! See if next line is needed
     &            IDIAG,
! use MXCCV
     &            MXCCV,NTOTORB,
     &            LCSBLK,
! Add idensi to input
     &            IDENSI,
     &            MK2REF,MK2DEL,
! Set ISPINFREE or remove it
     &            ISPINFREE,
! For properties
     &            NPROP,
! Add symmetry of T to input and use it
     &            ISYM_T,N_CC_AMP,
     &            NORB_F(MXPIRR),
! Note sure to keep IAAEXC_TYP
     &            MXSPOX,IAAEXC_TYP,
     &            MX_EXC_LEVEL,
     &            IMAXKRFLIP,IMUBMAX,
! add ICCSOLVE and I_DO_SBSPJA to input for different optimization
     &            ICCSOLVE, I_DO_SBSPJA,
! add an option to set LCCB from input
     &            LCCB,
! For block sizes in read and write
     &            I_READ_BLOCK,
! For CC excitation energies
     &            I_DO_CC_EXC_E, NEXC_PER_SYM(8),
! set and add input for MXCIV_ORIG
     &            MXCIV_ORIG,MXCIV,MXITLR,CCCONV_EX,
! add input for I_DO_CC_EXP
     &            I_DO_CC_EXP,
     &            NO4IDX, ! No 4-index transformation
! logicals
     &            IRECOM,    ! for real/complex algebra
     &            IFORM ! restart (set automatically)
