!-----------------------------------------
! Module with variables previously defined
! at huckel.h include file (COMMON /CBHUCKEL/)
!-----------------------------------------
!  HUCCNT    : Huckel factor, normally 1.75
!  HUCEXC(i) : atomic orbital energies for Huckel SOs organized after
!              symmetry , i=1..NHUCBA
!  NHUCAO(8) : number of Huckel SOs per symmetry
!  NHUCBA    : total number of Huckel SOs (e.g. 3 p-functions)
!  NHUCOCC(8): Huckel orbital occupation for NSYM symmetries 
!  IHUCPT(j) : pointer to Huckel AO "i" from total basis fu list
!              (used in her1*.F when Huckel matrix is constructed)
!              j=1..KMAX (KMAX is total number of symmetry independent
!                         AO shells, )
!  ISETHUCKEL : which basis set number is Huckel basis set
!  EWMO      : do EWMO instead of EHT
   module huckel_mod
   use maxorb_mod, only: MXSHEL
   implicit none
   real*8 :: HUCCNT, HUCEXC(MXSHEL)
   integer :: NHUCAO(8), NHUCBA, IHUCPT(MXSHEL), IPRHUC, ISETHUCKEL
   logical :: EWMO, HUCPROJCMO
   end module huckel_mod
