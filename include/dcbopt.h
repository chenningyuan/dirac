C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbopt.h
C SOURCE  : $Source$
C
C*** 2nd order opt. ***
C     IPROPT    - print level
C     IOPT_ITRINT    - LS/SS integrals from iteration ITRINT(1)/ITRINT(2)
C     IOPT_INTDEF    - default INTFLG
C     IOPT_INTFLG    - LL/SL/SS integrals (packed in bits 0-2)
C     IOPT_INTBUF    - value of INTFLG from last macro iteration
C     OPT_INTTYP     - string with integral types
C     IHMROP    - Hermitian symmetry of Fock matrices
C     ISYMOP    - Symmetry of Fock matrices
C     IFCKOP    - 
C     
C     NZVAR     - total number of variables
C     NZCONF    - number of configurational parameters
C     NZXOPE    - number of orbital electron-electron parameters
C     NZXOPP    - number of orbital electron-positron parameters
C     NZXOPT    - total number of orbital parameters
C     JOPSY     - symmetry of wave-function
C
C     OPT_INDSTR - inactive/active/secondary orbitals
C     
C     KZCONF    - pointer to configurational parameters
C     KZXOPE    - pointer to orbital e-e parameters
C     KZXOPP    - pointer to orbital e-p parameters
C
C     MOPT_MXMICRO   - maximum micro iterations
C     MOPT_MXMACRO   - maximum macro iterations
C
C     IOPT_ISTATE    - converge to state no. ISTATE (0=ground state)
C
C     FLAG:       (1)  Always Newton-Raphson iterations
C                 (2)  Always NEO iterations
C                 (19) Tight step contol for ground state also
C                 (20) Skip extra termination tests
C
      PARAMETER (NOPTFLAGS = 20)
      LOGICAL FLAG(NOPTFLAGS)
      CHARACTER OPT_INTTYP*9, OPT_CIPROGRAM*8
      COMMON /DCLOPT/ FLAG
      COMMON /DCCOPT/ OPT_INTTYP
      COMMON /DCIOPT/ IPROPT,IOPT_ITRINT(2),IOPT_INTDEF,
     &                IOPT_INTFLG,IOPT_INTBUF,
     &                IOPT_IHRMOP(2),IOPT_IFCKOP(2),IOPT_ISYMOP(2),
     &                NZVAR, NZVARQ ,NZCONF,NZCONFQ,
     &                NZXOPT,NZXOPTQ,NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,
     &                NZHOPT,NZHOPTQ,NZHOPE,NZHOPEQ,
     &                MOPT_MXMICRO,MOPT_MXMACRO,IOPT_ISTATE,
     $     IOPT_SYMMETRY,
     $     KTRLVL,JTRLVL,
     $     LENPV
C
      COMMON /DCROPT/ OPT_CNVINT(2),GNORM(5),EMCSCF,EMCOLD,DEPRED,BETA
     &                BETMIN,BETMAX,RTRUST,RTTOL,RATMIN,RATGOD,EMY,
     &                EACTIV,RATREJ,STPLEN,STPMAX,STPINC,STPRED,
     &                THQMIN,THQLIN,THQKVA,GAMMA,RBETA,
     &                ITMAC,ITMICT,STPC,STPO,ECORR,OPT_THRGRD
C
      LOGICAL OPT_SKIPEE,OPT_SKIPEP,OPT_NOPFQ
      COMMON /DCLOPT/ OPT_SKIPEE,OPT_SKIPEP,OPT_NOPFQ
C
      CHARACTER*72 OPT_INDSTR(3,2)
      COMMON /DCCOPT/
     $     OPT_INDSTR,
     $     OPT_CIPROGRAM
C
      COMMON /DCMOPT/ KZCONF,KZXOPE,KZXOPP
