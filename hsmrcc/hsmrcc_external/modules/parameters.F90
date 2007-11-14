module Parameters

   implicit none

   ! number of integers in floating point number (machine dependent)

   integer, parameter :: INTSFP = 1

   ! number of logical files and subroutines

   integer, parameter :: NVFIL  = 14
   integer, parameter :: NSUBR  = 17

   integer, parameter :: V2H2P  = 01
   integer, parameter :: VPPPP  = 02
   integer, parameter :: VHHHH  = 03
   integer, parameter :: VHHPP  = 04
   integer, parameter :: V2HPP  = 05
   integer, parameter :: V2PHH  = 06
   integer, parameter :: VPHHP  = 07
   integer, parameter :: VPH2P  = 08
   integer, parameter :: VPHPP  = 09
   integer, parameter :: V2HHP  = 10
   integer, parameter :: VHHHP  = 11
   integer, parameter :: FHH    = 12
   integer, parameter :: FHP    = 13
   integer, parameter :: FPP    = 14

   integer, parameter :: HEFFEC = 15
   integer, parameter :: RENORM = 16
   integer, parameter :: INTGEN = 17

   ! define alpha spin and beeta spin

   integer, parameter :: ALPHA_SPIN = +1
   integer, parameter :: BEETA_SPIN = -1

   ! define core, active and virtual

   integer, parameter :: CORE_ORB = +1
   integer, parameter :: ACTV_ORB =  0
   integer, parameter :: VIRT_ORB = -1

   ! define Q, P-M and P-I spaces

   integer, parameter :: Q_SPACE  = -1
   integer, parameter :: PM_SPACE = +1
   integer, parameter :: PI_SPACE =  0

   ! some double precision parameters

   double precision, parameter :: zero = 0.0d00
   double precision, parameter :: one  = 1.0d00
   double precision, parameter :: half = 0.5d00
   double precision, parameter :: NaN    = Z'7FF8000000000000'
   double precision, parameter :: NegInf = Z'FFF0000000000000'
   double precision, parameter :: PosInf = Z'7FF0000000000000'

   ! HSMRCC parameters

   integer :: NCMO
   integer :: NHOLE
   integer :: NPTCL
   integer :: NEMO
   integer :: NELE

   integer :: NVMO
   integer :: NHMO
   integer :: NPMO
   integer :: NTMO

   integer :: NHMOSQ
   integer :: NPMOSQ
   integer :: NVMOSQ

   integer :: NSING
   integer :: NDBLE
   integer :: NTSIZE

   integer :: MAXEXC

   integer :: NMAINORB
   integer :: NINTMORB
   integer :: NINTMORB_CORE
   integer :: NINTMORB_VIRT

   integer :: NSTATE
   integer :: NSTATE_MAIN
   integer :: NSTATE_INTM

   ! general print level parameter
   integer :: IPRINT = 0

end module Parameters
