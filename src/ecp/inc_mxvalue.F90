!-----------------------------------------
! Module with variables previously defined
! at inc_mxvalue.h include file
!-----------------------------------------
!  msu    = max number of symmetry inequivalent types of atoms.
!  msup increased to 1000 to allow for large number of PCs (PSB 05/03/20)
!  also change recp_lnk.F line 145 for local definition of msup
!  number of centers after symmetry operations
   module inc_mxvalue_mod
   implicit none
   integer, parameter :: MSUP = 203
   integer, parameter :: MSTUP = 8
   integer, parameter :: MRCRUP = 50
   integer, parameter :: MCONUP = 50
   integer, parameter :: MCUP = 50
   integer, parameter :: MCONSP = 200
   integer, parameter :: MSFUP = 200
   integer, parameter :: MNSFUP = (msfup*(msfup+1))/2
   integer, parameter :: MGUP = 49
   integer, parameter :: MCCUP = 182
   integer, parameter :: KAORDP = 240
   integer, parameter :: MGCSUP = 240
   integer, parameter :: MRUP = 100
   integer, parameter :: MCRUP = 140
   integer, parameter :: MSFRUP = 1000
   integer, parameter :: MNRUP = 5000
   end module inc_mxvalue_mod
