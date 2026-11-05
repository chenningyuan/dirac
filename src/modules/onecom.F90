!-----------------------------------------
! Module with variables previously defined
! at onecom.h include file (COMMON /ONECOM/)
!-----------------------------------------
   module onecom_mod
   implicit none
   real*8 :: CORAX, CORAY, CORAZ, CORBX, CORBY, CORBZ, SIGNBX, SIGNBY, SIGNBZ, HKAB
   integer :: NHKTA, NHKTB, KHKTA, KHKTB, KHKTAB, KCKTA, KCKTB, KCKTAB, IDENA, IDENB, &
              NCENTA, NCENTB, ICENTA, ICENTB, NUCA, NUCB, JSTA, JSTB, MULA, MULB, MAB, &
              JMAX, ISTEPU, ISTEPV, NAHGTF, NUMCFA, NUMCFB
   logical :: ONECEN, LDIAG, SPHRAB, SPHRA, SPHRB
   end module onecom_mod
