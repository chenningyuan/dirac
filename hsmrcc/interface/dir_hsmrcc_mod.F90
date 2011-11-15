!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module Dirac_HSMRCC
!----------------------------------------------------------------------------------
!
!     handles all necessary data from DIRAC for the standalone HSMRCC code 
!
!----------------------------------------------------------------------------------
  !  ... lu units of dirac prepared files for HSMRCC
!integer :: NTMO=0, NCMO=0, NVMO=0, NEMO=0, NHOLE=0, NPTCL=0, NELE=0
!integer :: NHMO=0, NPMO=0
!integer :: NFRZV, NACT

!real(kind=8) :: TotEn, NucEn, TrfTol

integer, public :: NCMO_DIRAC=0, NACT_DIRAC=0, NFRZV_DIRAC=-1, NEMO_DIRAC=-1

integer, public :: iprint_hscc = 0

logical, public :: nohsccfiles=.false., nohsccdirac=.false.

! ... stuff from HSMRCC code for preparing VINF,VIND,VINT, FHH, FHP, FPP files ...
!integer, parameter :: NVFIL  = 14
!integer, dimension(NVFIL) :: NumRec, RecLen
!integer, dimension(NVFIL) :: StartRec, EndngRec
!integer, dimension(NVFIL) :: TotNZInt, NumBlock

! ... all is in HSMRCC !

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_VINF_file
!
!            prepare *VINF - file for DIRAC
! 
!    Needs data from already prepared VINT,VIND integral files !
!
!    content of VINF file  ....
!
   ! ...............................................................................................
   ! ...  Rec     Module    : Variables present in the record                                    ...
   ! ...............................................................................................
   ! ...  01 [ Parameters ] : NCMO, NHOLE, NPTCL, NEMO                                           ...
   ! ...  02 [ Options    ] : OPTION                                                             ...
   ! ...  03 [ ConvMod    ] : MAXITER, CONVTOL                                                   ...
   ! ...  04 [ IntGenMod  ] : CUTOFF                                                             ...
   ! ...  05 [ SymInf     ] : point_group, naxis, norder, numirrep, irrep_product,               ...
   ! ...                      totally_symmetric_irrep, irrep_label                               ...
   ! ...  06 [ Shift      ] : CoreShiftVal, VirtShiftVal                                         ...
   ! ...  07 [ GMSMod     ] : NAO, NSO, NFRZV, TotEn, NucEn,                                     ...
   !                          TrfTol, IUNIQ, LABELSIZE, NINTMX                                   ...
   ! ...  08 [ GMSMod     ] : SpinOfOrb, OccGMSRef, IrrepOrb, IrpLabOrb,                         ...
   !                          OrbType, OrbSpace, OrbEnergy, OrbCoeff, HCoreMatrix                ...
   ! ...  09 [ IOModule   ] : NumRec, RecLen                                                     ...
   ! ...  10 [ IOModule   ] : StartRec(V2H2P), EndngRec(V2H2P), TotNZInt(V2H2P), NZIntRec_V2H2P  ...
   ! ...  11 [ IOModule   ] : StartRec(VPPPP), EndngRec(VPPPP), TotNZInt(VPPPP), NZIntRec_VPPPP  ...
   ! ...  12 [ IOModule   ] : StartRec(VHHHH), EndngRec(VHHHH), TotNZInt(VHHHH), NZIntRec_VHHHH  ...
   ! ...  13 [ IOModule   ] : StartRec(VHHPP), EndngRec(VHHPP), TotNZInt(VHHPP), NZIntRec_VHHPP  ...
   ! ...  14 [ IOModule   ] : StartRec(V2HPP), EndngRec(V2HPP), TotNZInt(V2HPP), NZIntRec_V2HPP  ...
   ! ...  15 [ IOModule   ] : StartRec(V2PHH), EndngRec(V2PHH), TotNZInt(V2PHH), NZIntRec_V2PHH  ...
   ! ...  16 [ IOModule   ] : StartRec(VPHHP), EndngRec(VPHHP), TotNZInt(VPHHP), NZIntRec_VPHHP  ...
   ! ...  17 [ IOModule   ] : StartRec(VPH2P), EndngRec(VPH2P), TotNZInt(VPH2P), NZIntRec_VPH2P  ...
   ! ...  18 [ IOModule   ] : StartRec(VPHPP), EndngRec(VPHPP), TotNZInt(VPHPP), NZIntRec_VPHPP  ...
   ! ...  19 [ IOModule   ] : StartRec(V2HHP), EndngRec(V2HHP), TotNZInt(V2HHP), NZIntRec_V2HHP  ...
   ! ...  20 [ IOModule   ] : StartRec(VHHHP), EndngRec(VHHHP), TotNZInt(VHHHP), NZIntRec_VHHHP  ...
   ! ...  21 [ IOModule   ] : StartRec(FHH)  , EndngRec(FHH)  , TotNZInt(FHH)                    ...
   ! ...  22 [ IOModule   ] : StartRec(FHP)  , EndngRec(FHP)  , TotNZInt(FHP)                    ...
   ! ...  23 [ IOModule   ] : StartRec(FPP)  , EndngRec(FPP)  , TotNZInt(FPP)                    ...
   ! ...............................................................................................
!   
!
!modules/iomodule.F90:   integer, dimension(NVFIL) :: NumRec, RecLen
!modules/iomodule.F90:   NumRec(V2H2P) = NHMO*NPMO
!modules/iomodule.F90:   NumRec(VPPPP) = NPMOSQ
!modules/iomodule.F90:   NumRec(VHHHH) = NHMOSQ
!modules/iomodule.F90:   NumRec(VHHPP) = NHMOSQ
!modules/iomodule.F90:   NumRec(V2HPP) = NHMO*NHMO
!modules/iomodule.F90:   NumRec(V2PHH) = NPMO*NPMO
!modules/iomodule.F90:   NumRec(VPHHP) = NHMO*NPMO
!modules/iomodule.F90:   NumRec(VPH2P) = NPMO*NPMO
!modules/iomodule.F90:   NumRec(VPHPP) = NPMO*NHMO
!modules/iomodule.F90:   NumRec(V2HHP) = NHMO*NHMO
!modules/iomodule.F90:   NumRec(VHHHP) = NHMOSQ
!modules/iomodule.F90:   NumRec(FHH)   = NHMO
!modules/iomodule.F90:   NumRec(FHP)   = NHMO
!modules/iomodule.F90:   NumRec(FPP)   = NPMO
!
!modules/iomodule.F90:   integer, dimension(NVFIL) :: NumRec, RecLen
!modules/iomodule.F90:   integer, dimension(NVFIL) :: StartRec, EndngRec
!modules/iomodule.F90:   integer, dimension(NVFIL) :: TotNZInt, NumBlock
!modules/parameters.F90:   integer, parameter :: NVFIL  = 14
!
!----------------------------------------------------------------------------------
use Memory
use Vectors
!use io_moltra
#include "implicit.h"

call Header('Output from "prepare_dirac_HSMRCC_VINF_file"',-1)


end subroutine prepare_dirac_HSMRCC_VINF_file


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_int_files
!-------------------------------------------------------------------------
!
!
   ! ... Generate Integrals ...

   ! ... The order of integral
   ! ... generation is the order
   ! ... in which NZIntRec info
   ! ... is written into FileVINF
   ! ... is the same as order in which
   ! ... Integrals Info is written to
   ! ... FileVIND & FileVINT
   ! ... FileFHH, FileFHP & FileFPP
   !
   ! ... The order of generation
   !
   !   Two-body integral files
   !
   !       01. V2H2P  <IJ||AB>  all I,J; all A,B
   !       02. VPPPP  <AB||CD>  A>B, C>D
   !       03. VHHHH  <IJ||KL>  I>J,K>L
   !       04. VHHPP  <IJ||AB>  I>J,A>B
   !       05. V2HPP  <IJ||AB>  all I,J ; A>B
   !       06. V2PHH  <AB||IJ>  all A,B; I>J
   !       07. VPHHP
   !       08. VPH2P
   !       09. VPHPP
   !       10. V2HHP
   !       11. VHHHP
   !
   !   One-body itegral files
   !
   !       01. FHH
   !       02. FHP
   !       03. FPP
   !     
!-------------------------------------------------------------------------
end subroutine prepare_dirac_HSMRCC_int_files

end module Dirac_HSMRCC
