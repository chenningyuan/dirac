!-----------------------------------------
! Module with variables previously defined
! at luctinp.h include file (COMMON /LUCTINFC/, /LUCTINFI/, /LUCTINFR/)
!-----------------------------------------
!     Common block for LUCITA input information
!     read in dirrdn.F/PSIINP/LUCTINP
!     Information is transferred to LUCITA routine dirluct.F
!     for further processing.
   module luctinp_mod
   implicit none
   integer, parameter :: NTABLE = 30
   integer, parameter :: MXNGAS = 16
   real*8 :: CTRUNC_FAC, MY_CONVERGENCE
   integer :: IMOKW(NTABLE), NROOTD, ISSYMD, NACTED, IMULTD, IPRNGD, IPRNLD, IDENSD, &
              MXHL1D, MXEL3D, INGASD, NSEQCD, IRSTLT, MXCIVE, IPARMODEL, ICIMAXITER, &
              IMAXBLKSIZE, I_USE_DIST_ROUTE, IN_MEMFAC
   character(len=72) :: CRDINA, CRDCCLR, CRDGAS(MXNGAS), CRDGOC(MXNGAS), CRDFRO, CRDRS1, &
                        CRDRS2, CRDRS3, WAFFCD, CALCTP, SZCALD, TITLUC, CALCTYP(MXNGAS*4), &
                        IUSEQ, UCMBCC, UNEWCCV, URES_CC
   end module luctinp_mod
