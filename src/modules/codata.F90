!-----------------------------------------
! Module with variables previously defined
! at codata.h include file (COMMON /CODVAL/)
!-----------------------------------------
! ** COMMON block for CODATA
!     CODSET - Codata set of values.
!              Available as: CODATAXX, where XX possible values are
!                            86, 98, 02, 06, 10, 14, 18, 22
!     Fundamental and derived constants declared as real*8 values
   module codata_mod
   implicit none
   real*8 :: XTANG, ECHARGE, HBAR, XFMOL, UMASS, PMASS, EMASS, CCM, FERMICC, S2THETAW, &
             AUTK, PLANCK, XTANGM10, CVEL, ALPHAC, ALPHA2, XTJ, XTKAYS, XTHZ, XTEV, XKJMOL, &
             XKCMOL, XTNM, XAJOUL, AUTIME, TESLA, DEBYE, XTKMML, XFAMU, XFMP, NMAGN, &
             NMAGNAU, EFAUMKSA, CMINV, NM, GFERMI, GFAC, G1FAC
   character(len=8) :: CODSET
   end module codata_mod
