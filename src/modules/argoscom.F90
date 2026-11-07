!-----------------------------------------
! Module with variables previously defined
! at argoscom.h include file (COMMON /ARGOSCOM/, /ECPINT/)
!-----------------------------------------
!    This is for the ECP routine
   module argoscom_mod
   implicit none
   integer :: IGRPAG, ECPNUC, ECPNONT(1000), NDPTAG, MAXANG(200), IROWMAX(200,200), &
              IREPARG(200,200,100)
   integer :: ECPINT_PRINT, ECPINT_DEBUG
   ! formerly a DATA-initialized local array of every unit including argoscom.h (read only)
   integer, parameter :: NDPTAGS(0:7) = [0, 0, 0, 0, 1, 1, 1, 7]
   end module argoscom_mod
