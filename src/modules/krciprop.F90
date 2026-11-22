!-----------------------------------------
! Module with variables previously defined
! at krciprop.h include file (COMMON /LUCIIPRP/, /LUCILPRP/, /LUCICPRP/, /LUCIPROPNUC/)
!-----------------------------------------
!     property calculations inside the KR-CI module in
!     DIRAC.
!     DOPROP_KRCI: activate property calculation with LUCI*
!                  modules
!     DOOSCILLST:  calculate oscillator strength for electronic
!                  transitions
!     dodipmom_krci: calculate permanent electronic dipole moments of
!                    electronic ground and excited states.
!     DOJZEXP:     calculate <j_z> expectation value
!                  (omega quantum number) for electronic states
!     NPROP_KRCI:  total number of one-electron operators
!     LPROP_KRCI:  pointer list to general common block list
!                  of one-electron operators
!     NPROP_ROOTS_KRCI:
!                  total number of states (roots) in property
!                  calculation
!     ISYMOPPRP_KRCI:
!                  symmetry of one-electron operators for the
!                  corresponding double group that is used in
!                  KR-CI:
!                  DIRAC: D2h --> KR-CI: C2h
!                  DIRAC: C2v --> KR-CI: C2
!     IN1ELPRP_KRCI:
!                  total number of integrals in KR-CI order
!                  for this one-electron operator
!     ISYMEIG_KRCI:
!                  symmetry irrep for each eigenstate
!     IXSYMACT  :  total number of "active" symmetry irreps
!     NXSYMPAIRS:  number of C/sigma symmetry pairs
!     IHTYPE_X  :  X operator type
!     LXPRPKRCI :  total dimension of KR-CI property matrix
!     ISCAL_T   :  total symmetric operator?
!     IPROPPRINT:  print level
!     FILE    : luciprop.h (now krciprop_mod)
!     ================================================================
   module krciprop_mod
   implicit none
   integer, parameter :: MXPROP_KRCI = 60
!  max. total number of roots in property run - value is 4 x MXROOT
!  defined in luciarel/mxpdim.inc which is the max. number of roots
!  per irrep.
   integer, parameter :: MXPRPE_KRCI = 750
   integer, parameter :: MAXNUCLEI = 15
   integer, parameter :: IPROPPRINT = 30
!  max symmetry irreps - this value has to be equal to
!  PARAMETER (MAX_NKRCI_MAX_SYM = 128) in dcbkrci.h - SK Oct 2008
   integer, parameter :: MXPROPKRCI_SYM = 128
   real*8 :: VKRCI_NUCSPIN(MAXNUCLEI), VKRCI_NUCMAGMOM(MAXNUCLEI)
   integer :: NPROP_KRCI, LPROP_KRCI(MXPROP_KRCI), NPROP_ROOTS_KRCI, &
              ISYMOPPRP_KRCI(MXPROP_KRCI), LXPRPKRCI, IN1ELPRP_KRCI(MXPROP_KRCI), &
              ISYMEIG_KRCI(MXPRPE_KRCI), IXSYMLIST(MXPROPKRCI_SYM), IXSYMACT, NXSYMPAIRS, &
              IHTYPE_X
   logical :: DOPROP_KRCI, TWOEPROP, DOOSCILLST, DOJZEXP, DOPROPREOD, DOSIGPROP, RUNXPROP, &
              DODIPMOM_KRCI, DOEEDM_KRCI, DOMHYP_KRCI, DOGENP_KRCI, DOENSPS_KRCI, &
              DONMQM_KRCI
   character(len=4) :: XREPEIG(MXPROPKRCI_SYM)
   end module krciprop_mod
