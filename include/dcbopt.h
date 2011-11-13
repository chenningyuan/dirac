!
! VERSION : $Revision: 7428 $
! DATE    : $Date: 2008-08-25 15:22:41 +0200 (Mon, 25 Aug 2008) $
! FILE    : dcbopt.h
!
!*** 2nd order opt. ***
!     IPROPT    - print level
!     IOPT_ITRINT   - LS/SS integrals from iteration ITRINT(1)/ITRINT(2)
!     IOPT_INTDEF    - default INTFLG
!     IOPT_INTFLG    - LL/SL/SS integrals (packed in bits 0-2)
!     IOPT_INTBUF    - value of INTFLG from last macro iteration
!     OPT_INTTYP     - string with integral types
!     IOPT_IHMROP    - Hermitian symmetry of Fock matrices
!     IOPT_ISYMOP    - Symmetry of Fock matrices
!     IOPT_IFCKOP    - Type of Fock matrices
!     
!     NZVAR     - total number of variables
!     NZCONF    - number of configurational parameters
!     NZXOPE    - number of orbital electron-electron parameters
!     NZXOPP    - number of orbital electron-positron parameters
!     NZXOPT    - total number of orbital parameters
!     JOPSY     - symmetry of wave-function
!
!     OPT_INDSTR - inactive/active/secondary orbitals
!     
!     KZCONF    - pointer to configurational parameters
!     KZXOPE    - pointer to orbital e-e parameters
!     KZXOPP    - pointer to orbital e-p parameters
!
!     MOPT_MXMICRO   - maximum micro iterations
!     MOPT_MXMACRO   - maximum macro iterations
!
!     IOPT_ISTATE    - converge to state no. ISTATE (0=ground state)
!
!     FLAG:       (1)  Always Newton-Raphson iterations
!                 (2)  Always NEO iterations
!                 (19) Tight step contol for ground state also
!                 (20) Skip extra termination tests
!
!     IDENSLR_STATE  - natural orbital occupation numbers for state xxx
!     ISHMEM_TYPE    - type of MPI-2 based shared memory level
!
!     JKRRUNTYPE     - determines KR run type:
!                                             0: second-order SCF
!                                             1: second-order MCSCF
!                                             2: KR-CI
!
!     Parameters for MCTYPE:
!
      INTEGER JDHF, JODHF, JCAS, JMCMIN, JGAS
      PARAMETER (JDHF = 0, JODHF = 1, JCAS = 2, JMCMIN = JCAS, JGAS =   &
     &           3)

      INTEGER IPROPT, IOPT_ITRINT(2), IOPT_INTDEF, IOPT_INTFLG,         &
     &        IOPT_INTBUF, IOPT_IHRMOP(2), IOPT_IFCKOP(2),              &
     &        IOPT_ISYMOP(2), NZVAR, NZVARQ, NZCONF, NZCONFQ, NZXOPT,   &
     &        NZXOPTQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ, NZHOPT,        &
     &        NZHOPTQ, NZHOPE, NZHOPEQ, MOPT_MXMICRO, MOPT_MXMACRO,     &
     &        MOPT_MAXBCK, NCIROOT, MXCIV, IRESTRK, IOPT_ISTATE, IANACI,&
     &        IOPT_SYMMETRY, KTRLVL, JTRLVL, LENPV, ITMAC, ITMICT,      &
     &        ITBCK, ITABS, ITMACN, MAXCIT, IOPTST, MK2DEL, MK2REF,     &
     &        MCTYPE, JHESSIAN, IKRMCCNO, IDOPARIO, IMEMFAC,            &
     &        IDENSLR_STATE, ISHMEM_TYPE, NDSH(2), NDSHT, N2DSHX,       &
     &        N2DSHXQ, NVSH(2), NVSHT, N2VSHX, N2VSHXQ, JKRRUNTYPE
      COMMON /DCIOPT/ IPROPT, IOPT_ITRINT, IOPT_INTDEF, IOPT_INTFLG,    &
     &                IOPT_INTBUF, IOPT_IHRMOP, IOPT_IFCKOP,            &
     &                IOPT_ISYMOP, NZVAR, NZVARQ, NZCONF, NZCONFQ,      &
     &                NZXOPT, NZXOPTQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ,&
     &                NZHOPT, NZHOPTQ, NZHOPE, NZHOPEQ, MOPT_MXMICRO,   &
     &                MOPT_MXMACRO, MOPT_MAXBCK, NCIROOT, MXCIV, IANACI,&
     &                IRESTRK, IOPT_ISTATE, IOPT_SYMMETRY, KTRLVL,      &
     &                JTRLVL, LENPV, ITMAC, ITMICT, ITBCK, ITABS,       &
     &                ITMACN, MAXCIT, IOPTST, MK2DEL, MK2REF, MCTYPE,   &
     &                JHESSIAN, IKRMCCNO, IDOPARIO, IMEMFAC,            &
     &                IDENSLR_STATE, ISHMEM_TYPE, NDSH, NDSHT, N2DSHX,  &
     &                N2DSHXQ, NVSH, NVSHT, N2VSHX, N2VSHXQ, JKRRUNTYPE

      REAL*8 OPT_CNVINT(2), GNORM(5), EMCSCF, EMCOLD, DEPRED, BETA,     &
     &       BETMIN, BETMAX, RTRUST, RTTOL, RATMIN, RATGOD, EMY, EACTIV,&
     &       RATREJ, STPLEN, STPMAX, STPINC, STPRED, THQMIN, THQLIN,    &
     &       THQKVA, GAMMA, RBETA, STPC, STPO, ECORR, OPT_THRGRD,       &
     &       STPLCI, STPLEE, STPLEP, DLVLSH, THRPCI, COMPFAC, ECORE_LR
      COMMON /DCROPT/ OPT_CNVINT, GNORM, EMCSCF, EMCOLD, DEPRED, BETA,  &
     &                BETMIN, BETMAX, RTRUST, RTTOL, RATMIN, RATGOD,    &
     &                EMY, EACTIV, RATREJ, STPLEN, STPMAX, STPINC,      &
     &                STPRED, THQMIN, THQLIN, THQKVA, GAMMA, RBETA,     &
     &                STPC, STPO, ECORR, OPT_THRGRD, STPLCI, STPLEE,    &
     &                STPLEP, DLVLSH, THRPCI, COMPFAC, ECORE_LR

      INTEGER NOPTFLAGS
      PARAMETER (NOPTFLAGS = 20)

      LOGICAL FLAG(NOPTFLAGS)
      LOGICAL OPT_SKIPEE, OPT_SKIPEP, OPT_NOPFQ, OPT_NOFQX
      LOGICAL OPT_NOCI, OPT_UCIBOS, OPT_MKCNST
      LOGICAL DONR, FROZCI, LDIRNR, NATONL, FOCKON, COMPROT, NATOLCR
      LOGICAL CRDFO_MAT, CWRTFO_MAT, CSHMEMO, CINT_REORD, CINT_SPLIT
      LOGICAL GENFOCK
      COMMON /DCLOPT/ OPT_SKIPEE, OPT_SKIPEP, OPT_NOPFQ, OPT_NOFQX,     &
     &                OPT_NOCI, OPT_UCIBOS, OPT_MKCNST, DONR, FROZCI,   &
     &                LDIRNR, FLAG, NATONL, FOCKON, COMPROT, NATOLCR,   &
     &                CRDFO_MAT, CWRTFO_MAT, CSHMEMO, CINT_REORD,       &
     &                CINT_SPLIT, GENFOCK

      CHARACTER*72 OPT_INDSTR(3, 2)
      CHARACTER*9 OPT_INTTYP
      CHARACTER*8 OPT_CIPROGRAM
      CHARACTER*72 OPT_FRZSTR(2), OPT_DELSTR(2)
      COMMON /DCCOPT/ OPT_INTTYP, OPT_INDSTR, OPT_CIPROGRAM, OPT_FRZSTR,&
     &                OPT_DELSTR

      INTEGER KZCONF, KZXOPE, KZXOPP, LZCONF, LZXOPE, LZXOPP
      COMMON /DCMOPT/ KZCONF, KZXOPE, KZXOPP, LZCONF, LZXOPE, LZXOPP

!     lzconf is integer length of WORK(KZCONF) etc.
!
!
      INTEGER NOPTTIM
      PARAMETER (NOPTTIM = 20)

!     Detailed timings are saved in CPUOPT and WLLOPT.
!
!     01: time in ROPTST     (wf. start guess)
!     02: time in RCIST      (CI start guess)
!     03: time in RGETH2     (FQ and H2AC)
!     04: time in RGETH2TX   (FQX, FQT, and H2ACX)
!     05: time in GMOLITX    (FCX, FVX, and FVT)
!     06: time in FMOLI      (one index transforms)
!     07: time in RTRACTL    (MO 4-index transformation)
!     08: time in RMAKDM     (density matrices)
!     09: time in RGRAD      (gradient)
!     10: time in RCIGRAD    (CI gradient)
!     11: time for CI sigma vectors
!     12: time for RFCKMAT   (FC and FV)
!     13: time for XRSSEP    (orbital parts of orbital and conf. sigma vec)
!     14: time for XRSTDM    (transition density matrices)
!     15: time in RTR1H1     (one index transforms)
!     16: time for RSIGOC, RSIGOO (orbital sigma vectors)
!     17: time for RSIGCO    (conf part of orbital sigma vectors)
!     18: time for XRSSVC    (sigma vectors)
!
      REAL*8 CPUOPT(NOPTTIM), WLLOPT(NOPTTIM)
      COMMON /DCTOPT/ CPUOPT, WLLOPT
