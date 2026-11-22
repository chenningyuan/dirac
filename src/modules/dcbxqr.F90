!-----------------------------------------
! Module with variables previously defined
! at dcbxqr.h include file (COMMON /XCBLQR/, /XCBRQR/, /XCBIQR/, /XCBCQR/)
!-----------------------------------------
! dcbxqr.h - control common block for pamxqr.h and pamtpa.h
! **************************************************************
! ***********    Quadratic response  <<A;B,C>>    **************
! **************************************************************
! MAXQR  - max number of allowed QR functions
! MAXLQR - max number of operators specified
! MAXFQR - max number of frequencies specified
! MZYEE  - max length of electronic response vectors
! MZYEP  - max length of positronic response vectors
! Used to read in USER specified requests
!----------------------------------------
! LAQROP - list of A-operator pointers
! LBQROP - list of B-operator pointers
! LCQROP - list of C-operator pointers
! IPRXQR - general print level
! ITRXQR - max number of iterations
! NBQRFR - number of B frequencies
! BQRFR  - B frequencies
! NCQRFR - number of C frequencies
! CQRFR  - C frequencies
! XQR_SKIPEE - do not include e-e rotations
! XQR_SKIPEP - do not include e-p rotations
! XQR_ALLCMB - calculate all response functions disregarding
!              permutational symmetries (debug option)
! XQR_LEVICI - Include only functions with 3 different cartesian
!              indices (for .VERDET)
! XQR_INDSTR - list of orbitals to include in rotation vectors
! Used to store the necessary linear response equations
!------------------------------------------------------
! NQROP  - number of unique operator pointers which are to
!          to be solved linear response for
! LQROP  - list of operator pointers, and number of
!          frequencies for that operator
! QRFREQ - list of frequencies for operators
! Used to define nontrivial QR response functions
!------------------------------------------------
! NQRHYP  - number of nontrivial QR functions to be computed
! LQRHYP  - list of A/B/C operator pointers
! QRFRHYP - list of A/B/C frequencies
! Used in TPA calculations (EXCIT is also used by EXCPRP)
!--------------------------------------------------------
! IEXCOFF - offset when creating direct access labels for excitations
! TPACNV  - number of states per irrep
! TPASIM  -
! TPASTV  -
! EXCIT   - list of excitation frequencies in each irrep
! Used in EXCPRP calculations
! ---------------------------
! LSTATE - list of left excited states
! RSTATE - list of right excited states
   module dcbxqr_mod
   implicit none
   integer, parameter :: MAXLQR = 10
   integer, parameter :: MAXFQR = 15
   integer, parameter :: MAXQR = 70
   integer, parameter :: IEXCOFF = 90
   integer, parameter :: LUQRINFO = 40
   integer, parameter :: LUQREE = 41
   integer, parameter :: LUQREP = 42
   real*8 :: THCQR, BQRFR(MAXFQR), CQRFR(MAXFQR), QRFREQ(3*MAXFQR,3*MAXLQR), &
             QRFRHYP(MAXQR,3), CNVXQR(2), RESXQR, EXCIT(MAXFQR,8), DAMPFREQ_QR
   integer :: IPRXQR, ITRXQR, NBQRFR, NCQRFR, NQROP, MAXQRM, LAQROP(MAXLQR), NAQROP, &
              LBQROP(MAXLQR), NBQROP, LCQROP(MAXLQR), NCQROP, LQROP(3*MAXLQR,2), &
              LQRHYP(MAXQR,3), NQRHYP, INTXQR, ITRIQR(2), MZYEE, MZYEP, TPACNV(8), &
              TPASIM(8), TPASTV(8), LSTATE(8), RSTATE(8), NDAMP_QR
   logical :: XQRNRM, XQRDIH, XQR_SKIPEE, XQR_SKIPEP, XQR_ALLCMB, XQR_LEVICI
   character(len=72) :: XQR_INDSTR(3,2)
   end module dcbxqr_mod
