!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
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
!#include "implicit.h"

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


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_HSMRCC_Ffiles
!----------------------------------------------------------------
! prepare FHH, FPH, FPP HSMRCC files
!----------------------------------------------------------------
Call Header('Output from prepare_HSMRCC_Ffiles',-1)
end subroutine prepare_HSMRCC_Ffiles

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_hscc_VINF_data
!----------------------------------------------------------
!
!    set important HSCC parameters
!  for VINF file
!
!----------------------------------------------------------
use io_moltra
use Dirac_HSMRCC
use Parameters
use GMSMod
implicit none
#include "priunit.h"
integer :: numirrep_dirac,i,ii
integer, external :: GetNBSYM

call header('Output from set_hscc_data',-1)

!  NCMO_DIRAC,  NACT_DIRAC
 NAO = NORBT

!
!NEMO = NAO ! WRONG !!!
!NEMO=NEMO_DIRAC

! ... active space - double excitations
NPTCL = 2
NHOLE = 2

!NACT = NACT_DIR
!NCMO = NCORE_DIR

 NCMO = NCMO_DIRAC

NFRZV = NFRZV_DIRAC

NSO = NSP

! ...  prepare variable  ...
NVMO = NHOLE + NPTCL

call count_nelec_tot
write(LUPRI,*) 'nelec_tot=',nelec_tot
numirrep_dirac = GetNBSYM()
write(LUPRI,*) 'numirrep_dirac=',numirrep_dirac

if (get_iprint() >=8) then
 write(LUPRI,*)  "going to call alloc_VINF_arrays"
endif

! ... set options ...
OPTION(HSCC_DIRAC_MODE)=.false.
OPTION(DIRAC_MODE)     =.true. ! for further "alloc_VINF_arrays"

call alloc_VINF_arrays

! ... fill these arrays ...
call fill_VINF_arrays

call GenHSMRpars

end subroutine set_hscc_VINF_data

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_VINF_arrays
!
! fill SpinOf, OccGMSRef, IrrepOrb, IrpLabOrb, OrbType, OrbSpace, OrbEnergy
!
! IrpLabOrb is integer for char ! all arrays are in OrbMod module
!
use io_moltra
use Dirac_HSMRCC
use Parameters
! use OrbMod ! TODO: move here print_gms_data ?
use GMSMod

implicit none
#include "priunit.h"
character(len=4) :: repna_temp_a, repna_temp_b
integer :: multbs(1:8,1:8)
integer :: numirrep_dirac,i,ii

if (get_iprint() >=8) then
  write(LUPRI,'(/,4x,a,/)') "Output from fill_VINF_arrays"
endif

write(LUPRI,*) 'numirrep_dirac=',numirrep_dirac

write(LUPRI,'(/,2x,a,i4,a)') 'Irreps (1-',NSP,' )'
do i=1, NSP
  write(LUPRI,'(2x,i3,a,2i2,a,a4,a,f12.7)')   &  
  i,': irrep no ',IREPSP(i),IREPSP(i+NSP),'- irrep label:',REPNA(IREPSP(i+NSP)),' orb.energy:',EPS1(i) 
enddo
write(LUPRI,*)

! once again...
do i=1, NSP/2
 write(LUPRI,'(2x,i3,1x,a4,1x,a4)') i,REPNA(IREPSP(i+NSP)),REPNA(IREPSP(i+NSP/2+NSP))
enddo

NEMO_DIRAC=0

! fill orb.en.
ii = 1

do i=1, NSP/2
! write(LUPRI,*) 'EPS1(',i,')=',EPS1(i),EPS1(i+NSP/2)
 OrbEnergy(ii)   = EPS1(i)
 OrbEnergy(ii+1) = EPS1(i+NSP/2)
 repna_temp_a = REPNA(IREPSP(i+NSP))
 repna_temp_b = REPNA(IREPSP(i+NSP/2+NSP))
 write(LUPRI,'(2x,i3,1x,a4,i2,1x,a4,i2)') i, repna_temp_a,IREPSP(i), repna_temp_b, IREPSP(i+NSP/2)
 IrrepOrb(ii)    = IREPSP(i+NSP)
 read(repna_temp_a,'(a4)') IrpLabOrb(ii)
 IrrepOrb(ii+1)  = IREPSP(i+NSP) ! previously saved irrep...
 read(repna_temp_b,'(a4)') IrpLabOrb(ii+1)
 SpinOfOrb(ii)   = +1
 SpinOfOrb(ii+1) = -1
! write(LUPRI,'(a,2i4,2f12.7)') 'orb_en:',ii,ii+1,OrbEnergy%p_rA1D(ii),OrbEnergy%p_rA1D(ii+1)
 MOLTRA_indx(ii)   = i
 MOLTRA_indx(ii+1) = i+(NSP/2)
 ii = ii + 2
enddo

write(LUPRI,'(/,2x,a,i3)') 'IrrepOrb, SpinOfOrb ... arrays filled, NSP=',NSP

NEMO_DIRAC=0

do i=1,NSP

 if (i<=nelec_tot) then
  OccGMSRef(i)=1
 else
  OccGMSRef(i)=0
 endif
!
 if (i<=NCMO_DIRAC) then
  OrbType(i) = -1
  OrbSpace(i) = -1
 else if ( i>NCMO_DIRAC .and. i<=NCMO_DIRAC+NACT_DIRAC ) then
  OrbType(i) = 0
  OrbSpace(i) = 1
 else if ( i>NCMO_DIRAC+NACT_DIRAC ) then
  OrbType(i) = 1
  NEMO_DIRAC=NEMO_DIRAC+1
  OrbSpace(i) = -1
 else
  call quit('dead branch !')
 endif

 if (iprint_hscc >= 5) then
   call print_gms_data_line(i)
 endif

enddo

! .. update
SpinOfActvOrb = SpinOfOrb(NCMO+1:NCMO+NVMO)
OccRef = OccGMSRef
totally_symmetric_irrep = 1
numirrep = numirrep_dirac ! .... NOT initialized !
NEMO= NEMO_DIRAC

if (iprint_hscc >= 5) then
  write(LUPRI,'(/,2x,a,/)') '***  HSMRCC active space info ***'
  write(LUPRI,'(4x,a,i3)') 'NCMO_DIRAC=',NCMO_DIRAC
  write(LUPRI,'(4x,a,i3)') 'NACT_DIRAC=',NACT_DIRAC
  write(LUPRI,'(4x,a,i3)') 'NEMO_DIRAC=',NEMO_DIRAC
endif

if (NEMO_DIRAC == -1) call quit('NEMO_DIRAC bad number !')

call print_gms_data

end subroutine fill_VINF_arrays

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine make_VIND_VINT_files
!----------------------------------------------------------------
!
!    prepare two-electron integral&indexes files : VINT & VIND
!
!----------------------------------------------------------------
use Memory
use Vectors

use io_moltra
use Dirac_HSMRCC
use Parameters
use GMSMod
use IntGenMod
implicit none
#include "priunit.h"

integer, external :: get_moltra_index
integer :: StVRec,StRec, EnRec, NxtRec

!integer         , dimension(NZVINTBUFSIZE) :: NZVPkIBuf
!double precision, dimension(NZVINTBUFSIZE) :: NZVIntBuf

real*8 :: get_moltra_integral
real*8 :: zerothrs=1.0d-14 ! to be replaced byu CUTOFF

integer :: IA, IPQ, inz, IJ, IB, ii_mo, ia_mo,NZVIntBufInd, icount, II, ij_mo, ib_mo
integer :: IC, ID, ic_mo, id_mo, IK, IL, ik_mo, il_mo
real*8 :: SInt_ijab,SInt_ijba,ASInt_ijaba,SInt_abcd,SInt_abdc,ASInt_abcd
real*8 :: ASInt_ijab, SInt_ijlk, SInt_ijkl, ASInt_ijkl, ASInt_ijb
real*8 :: SInt_abij,SInt_abji,ASInt_abij,ASInt_ajib,ASInt_aibc,SInt_ajib,SInt_ajbi,SInt_aibc,SInt_aicb
real*8 :: SInt_ijka,SInt_ijak,ASInt_ijka

call print_basic_params
call print_gms_parameters

NZVIntBufInd=0
NZVIntBuf    = DBLE(0)
NZIntRec_V2H2P = 0 ! array
NZVPkIBuf    = 0
TotNZInt(V2H2P) = 0

StVRec=1
StartRec(V2H2P) = StVRec
EndngRec = DBLE(0)
StartRec = StVRec
NxtRec   = StVRec

call header('two-electron MOLTRA transformed integrals sorted for HSMRCC',-1)

write(LUPRI,'(/,4x,a,/)') '01. V2H2P, 1324, <IJ||AB>'
icount=0 ! nonzero integrals counter
do II=1,NHMO
do IA=1,NPMO
 IPQ=RecNumHP(II,IA)
 NZIntRec_V2H2P(IPQ)=0 
 inz=0
do IJ=1,NHMO
do IB=1,NPMO
 ii_mo=MOLTRA_indx(II);      ij_mo=MOLTRA_indx(IJ); 
 ia_mo=MOLTRA_indx(IA+NCMO); ib_mo=MOLTRA_indx(IB+NCMO)
 SInt_ijab=get_moltra_integral(ii_mo,ia_mo,ij_mo,ib_mo)
 SInt_ijba=get_moltra_integral(ii_mo,ib_mo,ij_mo,ia_mo)
 ASInt_ijab=SInt_ijab-SInt_ijba
 if (DABs(ASInt_ijab).gt.zerothrs) then
  icount = icount + 1; inz = inz + 1
  write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <IJ||AB>, <',II,' ',IJ,'||',  & 
       IA,' ',IB,'>=',ASInt_ijab,' (ia|jb):',ii_mo,ia_mo,ij_mo,ib_mo
  ! save
  NZVIntBufInd=NZVIntBufInd+1
  NZIntRec_V2H2P(IPQ)=NZIntRec_V2H2P(IPQ)+1
  if ( NZVIntBufInd > NZVINTBUFSIZE ) then
     if (.not.nohsccfiles) then
       write(UVINT,rec=NxtRec) NZVIntBuf
       write(UVIND,rec=NxtRec) NZVPkIBuf
     endif
     NxtRec = NxtRec + 1
     NZVIntBufInd = 1
     NZVPkIBuf = 0
     NZVIntBuf = DBLE(0) 
  endif 
  ! ... add nonzero integral into integral-buffer
  NZVIntBuf(NZVIntBufInd) = ASInt_ijab
  ! ... add packed (R,S)-indexes into index-buffer
  NZVPkIBuf(NZVIntBufInd) = PackRSHP(IJ,IB)
 endif
enddo
enddo
   if (inz.gt.0) then
     write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,A=',II,IA,' is :',inz
   endif
enddo
enddo

TotNZInt(V2H2P) = TotNZInt(V2H2P) + NZIntRec_V2H2P(IPQ)
EndngRec(V2H2P) = NxtRec
if ( TotNZInt(V2H2P) /= 0 ) then
  if (.not.nohsccfiles) then
    write(UVINT,rec=NxtRec) NZVIntBuf
    write(UVIND,rec=NxtRec) NZVPkIBuf
  endif
  StVRec = EndngRec(V2H2P) + 1
endif


write(LUPRI,'(/,4x,a,/)') '02. VPPPP  <AB||CD>'
icount=0
do IA=2,NPMO
do IB=1,IA-1
 inz=0
do IC=2,NPMO
do ID=1,IC-1

ia_mo=MOLTRA_indx(IA+NCMO); ib_mo=MOLTRA_indx(IB+NCMO)
ic_mo=MOLTRA_indx(IC+NCMO); id_mo=MOLTRA_indx(ID+NCMO)
SInt_abcd=get_moltra_integral(ia_mo,ic_mo,ib_mo,id_mo) ! (ac|bd) = <ab|cd>
SInt_abdc=get_moltra_integral(ia_mo,id_mo,ib_mo,ic_mo) ! (ad|bc) = <ab|dc>
ASInt_abcd=SInt_abcd-SInt_abdc ! <AB||CD>=<ab|cd>-<ab|dc>
if (DABs(ASInt_abcd).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <AB||CD>, <',IA,' ',IB,'||',  & 
       IC,' ',ID,'>=',ASInt_abcd,' (ac|bd):',ia_mo,ic_mo,ib_mo,id_mo
endif
enddo
enddo
   write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair A,B=',IA,IB,' is :',inz
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '03. VHHHH  <IJ||KL>'
icount=0
do II=2,NHMO
do IJ=1,II-1
 inz=0
do IK=2,NHMO
do IL=1,IK-1
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ik_mo=MOLTRA_indx(IK);il_mo=MOLTRA_indx(IL)
 SInt_ijkl=get_moltra_integral(ii_mo,ik_mo,ij_mo,il_mo) 
 SInt_ijlk=get_moltra_integral(ii_mo,il_mo,ij_mo,ik_mo)
 ASInt_ijkl=SInt_ijkl-SInt_ijlk
 if (DABs(ASInt_ijkl).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <IJ||KL>, <',II,' ',IJ,'||',  & 
       IK,' ',IL,'>=',ASInt_ijkl,' (ij|kl):',ii_mo,ij_mo,ik_mo,il_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,J=',II,IJ,' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '04. VHHPP  <IJ||AB>'
icount=0
do II=2, NHMO
do IJ=1, II-1
 inz=0
do IA=2, NPMO
do IB=1, IA-1
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ia_mo=MOLTRA_indx(NCMO+IA);ib_mo=MOLTRA_indx(NCMO+IB)
 SInt_ijab=get_moltra_integral(ii_mo,ia_mo,ij_mo,ib_mo) ! (ia|jb)=<ij|ab>
 SInt_ijba=get_moltra_integral(ii_mo,ib_mo,ij_mo,ia_mo) ! (ib|ja)=<ij|ab>
 ASInt_ijab=SInt_ijab-SInt_ijba ! <ij||ab>=<ij|ab>-<ij|ba>
 if (DABs(ASInt_ijab).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <IJ||AB>, <',II,' ',IJ,'||',  & 
       IA,' ',IB,'>=',ASInt_ijab,' (ij|ab):',ii_mo,ij_mo,ia_mo,ib_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,J=',II,IJ,' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '05. V2HPP  <IJ||AB>'
icount=0
do II=1,NHMO
do IJ=1,NHMO
 inz=0
do IA=2,NPMO
do IB=1,IA-1
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ia_mo=MOLTRA_indx(NCMO+IA);ib_mo=MOLTRA_indx(NCMO+IB)
 SInt_ijab=get_moltra_integral(ii_mo,ia_mo,ij_mo,ib_mo)  ! (ia|jb)=<ij|ab>
 SInt_ijba=get_moltra_integral(ii_mo,ib_mo,ij_mo,ia_mo)  ! (ib|ja)=<ij|ba>
 ASInt_ijab=SInt_ijab-SInt_ijba ! <ij||ab>=<ij|ab>-<ij|ba>
 if (DABs(ASInt_ijab).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <IJ||AB>, <',II,' ',IJ,'||',  & 
       IA,' ',IB,'>=',ASInt_ijb,' (ij|ab):',ii_mo,ij_mo,ia_mo,ib_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,J=',II,IJ,' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '06. V2PHH, 1234, <AB||IJ>'
icount=0
do IA=1,NPMO
do IB=1,NPMO
 inz=0
do II=2,NHMO
do IJ=2,II-1
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ia_mo=MOLTRA_indx(NCMO+IA);ib_mo=MOLTRA_indx(NCMO+IB)
 SInt_abij=get_moltra_integral(ia_mo,ii_mo,ib_mo,ij_mo) ! (ai|bj)=<ab|ij>
 SInt_abji=get_moltra_integral(ia_mo,ij_mo,ib_mo,ii_mo) ! (aj|bi)=<ab|ji>
 ASInt_abij=SInt_abij-SInt_abji ! <AB||IJ>
 if (DABs(ASInt_abij).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <AB||IJ>, <',IA,' ',IB,'||',  & 
       II,' ',IJ,'>=',ASInt_abij,' (ab|ij):',ia_mo,ib_mo,ii_mo,ij_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair A,B=',IA,IB,' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '07. VPHHP, 3124,<AJ||IB>'
icount=0
do II=1, NHMO
do IA=1, NPMO
 inz=0
do IJ=1, NHMO
do IB=1, NPMO
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ia_mo=MOLTRA_indx(NCMO+IA);ib_mo=MOLTRA_indx(NCMO+IB)
 SInt_ajib=get_moltra_integral(ia_mo,ii_mo,ij_mo,ib_mo) 
 SInt_ajbi=get_moltra_integral(ia_mo,ib_mo,ij_mo,ii_mo) 
 ASInt_ajib=SInt_ajib-SInt_ajbi
 if (DABs(ASInt_ajib).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, & 
       '. 3124, <AJ||IB>, <',IA,' ',IJ,'||',   & 
       II,' ',IB,'>=',ASInt_ajib,' (ij|ab):',ii_mo,ij_mo,ia_mo,ib_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,A=',II,IA, & 
    ' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '08. VPH2P, 1324, <AI||BC>'
icount=0
do IA=1,NPMO
do IB=1,NPMO
 inz=0
do II=1,NHMO
do IC=1,NPMO
 ii_mo=MOLTRA_indx(II);ia_mo=MOLTRA_indx(NCMO+IA);ib_mo=MOLTRA_indx(NCMO+IB);ic_mo=MOLTRA_indx(NCMO+IC)
 SInt_aibc=get_moltra_integral(ia_mo,ib_mo,ii_mo,ic_mo) 
 SInt_aicb=get_moltra_integral(ia_mo,ic_mo,ii_mo,ib_mo) 
 ASInt_aibc=SInt_aibc-SInt_aicb
 if (DABs(ASInt_aibc).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1324, <AI||BC>, <',IA,' ',II,'||',  & 
       IB,' ',IC,'>=',ASInt_aibc,' (ab|ic):',ia_mo,ib_mo,ii_mo,ic_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair A,B=',IA,IB, & 
    ' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '09. VPHPP, 1234, <AI||BC>'
icount=0
do IA=1,NPMO
do II=1,NHMO
 inz=0
do IB=2,NPMO
do IC=1,IB-1
 ii_mo=MOLTRA_indx(II);ia_mo=MOLTRA_indx(NCMO+IA);
 ib_mo=MOLTRA_indx(NCMO+IB);ic_mo=MOLTRA_indx(NCMO+IC)
 SInt_aibc=get_moltra_integral(ia_mo,ib_mo,ii_mo,ic_mo) 
 SInt_aicb=get_moltra_integral(ia_mo,ic_mo,ii_mo,ib_mo) 
 ASInt_aibc=SInt_aibc-SInt_aicb
 if (DABs(ASInt_aibc).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <AI||BC>, <',IA,' ',II,'||',  & 
       IB,' ',IC,'>=',ASInt_aibc,' (ab|ic):',ia_mo,ib_mo,ii_mo,ij_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair A,I=',IA,II, & 
    ' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '10. V2HHP, 1324, <IJ||KA>'
icount=0
do II=1,NHMO
do IK=1,NHMO
 inz=0
do IJ=1,NHMO
do IA=1,NPMO
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ia_mo=MOLTRA_indx(NCMO+IA);ik_mo=MOLTRA_indx(IK)
 SInt_ijka=get_moltra_integral(ii_mo,ik_mo,ij_mo,ia_mo) 
 SInt_ijak=get_moltra_integral(ii_mo,ia_mo,ij_mo,ik_mo)
 ASInt_ijka=SInt_ijka-SInt_ijak
 if (DABs(ASInt_ijka).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1324, <IJ||KA>, <',II,' ',IJ,'||',  & 
       IK,' ',IA,'>=',ASInt_ijka,' (ik|ja):',ii_mo,ik_mo,ij_mo,ia_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,K=',II,IK, & 
    ' is :',inz
  endif
enddo
enddo

write(LUPRI,'(/,4x,a,/)') '11. VHHHP, 1234, <IJ||KA>'
icount=0
do II=2,NHMO
do IJ=1,II-1
 inz=0
do IK=1,NHMO
do IA=1,NPMO
 ii_mo=MOLTRA_indx(II);ij_mo=MOLTRA_indx(IJ);ik_mo=MOLTRA_indx(IK);ia_mo=MOLTRA_indx(NCMO+IA)
 SInt_ijka=get_moltra_integral(ii_mo,ik_mo,ij_mo,ia_mo) 
 SInt_ijak=get_moltra_integral(ii_mo,ia_mo,ij_mo,ik_mo)
 ASInt_ijka=SInt_ijka-SInt_ijak
 if (DABs(ASInt_ijka).gt.zerothrs) then
   icount = icount + 1; inz = inz + 1
   write(LUPRI,'(2x,i5,a,i3,a,i3,a,i3,a,i3,a,f18.12,a,4i3)') icount, &  
       '. 1234, <IJ||KA>, <',II,' ',IJ,'||',  & 
       IK,' ',IA,'>=',ASInt_ijka,' (ik|ja):',ii_mo,ik_mo,ij_mo,ia_mo
 endif
enddo
enddo
  if (inz.gt.0) then
    write(LUPRI,'(1x,a,2i4,a,i4)') 'Number of nonzero integrals for pair I,J=',II,IJ,' is :',inz
  endif
enddo
enddo

end subroutine make_VIND_VINT_files

subroutine PAMHSMRCC
! called from main/dirac.F
use memory_allocator
use allocator_cfg

implicit none

real(8), allocatable :: WORK(:)
integer :: lwork,  WORK_memsize

LWORK = WORK_memsize
call alloc(WORK,LWORK,id='WORK in PAMHSMRCC')
call HSCC_STUFF(WORK,LWORK)
call dirprg_hscc
call dealloc(WORK, id='WORK in PAMHSMRCC')
end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCC_STUFF(WORK,LWORK)
!----------------------------------------------------------------------------------
!
!  Called from PAMTRA, if defined (MOD_HSCC)
!
!  Written by Miro Ilias, Tel Aviv/Prievidza 2007,2008,2009
!
!----------------------------------------------------------------------------------
use Memory
use Dirac_HSMRCC
use io_moltra

implicit none
#include "priunit.h"
#include "dcbpsi.h"

integer,intent(in) :: LWORK
real*8, intent(in) :: WORK(LWORK)
logical :: check_point = .false.
integer :: ipos1,ipos2,ipos3,ipos4
logical :: check_MOints = .true.

! ... first run test_oof features
if (iprint_hscc >= 10) then
  call TEST_OOF(WORK,LWORK)
  call TEST_OO_MEMORY(WORK,LWORK)
endif

if (.not.DOHSCC) then
 write(LUPRI,'(2x,a)') 'DOHSCC=.false., no HSCC data preparation'
 return
endif

write(LUPRI,*)
write(LUPRI,*) 'in HSCC_STUFF !! LWORK=',LWORK

call print_HSMRCC_Dirac_info

call StickMemory(WORK,LWORK)

call read_moltra_files 

call UnstickMemory
end subroutine HSCC_STUFF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine read_moltra_files
!
! read MRCONEE & MDCINT; prepare sorted integrals for HSMRCC
!
  use Memory
  use Vectors
  use io_moltra
  use GMSMod
  use VINFMod
  use Dirac_HSMRCC
implicit none
#include "priunit.h"
call read_mrconee;  call print_mrconee_data
call read_mdcint
call set_dirac_hscc_file_names
call set_hscc_VINF_data
call InitIOMod
if (.not.nohsccfiles) then
 call OpenIntFiles
endif
! 2-el transformed integrals
call make_VIND_VINT_files
! ... deallocate big arrays for 2-el MO integrals & their indeces (allocation was done in read_mdcint)
call DeAllocObject(pos_kl)
call DeAllocObject(pos_ij)
call DeAllocObject(cbuf)
! 1-el. transformed integrals
if (.not.nohsccfiles) then
call CloseIntFiles ! close HSMRCC integral files 
! ... finally update VINF file
call OpenVINF
call UpdateVINF
call CloseVINF
endif
! ... deallocate VINF arrays ...
call EndGMSMod
call reset_Init
if (get_iprint() >=6) then
  write(LUPRI,*) "...end routine >>>read_moltra_files<<<"
endif
end subroutine read_moltra_files

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine HSCCINP(WORD,reset)
!-----------------------------------------------------------------------------
!
!  Read parameters for the external HSMRCC program, which are passed 
! through the formatted input file, together with other DIRAC input
! 
! Called from TRAINP if DOHSCC=.true., i.e. at the beginning
!
! Written by Miro ILIAS, Tel Aviv, Prievidza 2007, 2008
!
!-----------------------------------------------------------------------------

 use Dirac_HSMRCC  ! ... parameters for the external HSMRCC code
 use Parameters ! .. from hsmrcc_external/modules

implicit none
#include "priunit.h"

Logical, save :: set = .false.
Logical, intent(in) :: reset

CHARACTER:: PROMPT*1=' '
CHARACTER(len=7) :: WORD, KEYWORD=' '

! ... nice F90 version 
character(len=7) :: TABLE_TEMP, TABLE(10) =  (/'.PRINT ','.NACT  ','.NCORE ','.NOHSD ','.DRYRUN', & 
                                              '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX' /)

 if (set) then
  if (reset) set = .false.
  write(lupri,*) 'HSCCINP was called .... leaving at the beginning ... '
  Return
 endif
 if (reset) then
   set = .false.
 else
   set = .true.
 endif

 write(LUPRI,*) 'in HSCCINP, size of TABLE=',size(TABLE)

! REWIND (LUCMD,IOSTAT=IOS); print *,'input file rewind...' 

! Do While (KEYWORD /= '*HSCC  ')
!   READ (LUCMD,'(A7)',ERR=910, END=910) KEYWORD
!   print *,'read KEYWORD=',KEYWORD
! EndDo

 call FLSHFO(LUPRI)

!
!  Initializations...
!

Do 

  Read(LUCMD,'(A7)',ERR=910, END=930) WORD; PROMPT = WORD(1:1);  print *,'WORD=',WORD
  KEYWORD = WORD 

! Now decide !

  If (PROMPT == '*' ) Then
    print *,'WORD with "*" catched !', WORD
    Exit
  EndIf

  If ((PROMPT == '#') .or. (PROMPT == '!')) Then
     print *,'comment catched !', WORD
     Cycle ! ... continue the cycle..
  EndIf
  
  If (PROMPT == '.') Then
    ! ... keywords are beginning ...
    If (Any(TABLE(1:size(TABLE)) == WORD )) Then
     Select Case (WORD)
      Case ('.OPTION')
        Call PRTAB(size(TABLE), TABLE, KEYWORD//' input keywords',LUPRI)
      Case ('.NCORE')
        READ(LUCMD,*) NCMO_DIRAC
        write(LUPRI,*) ' read NCMO_DIRAC=',NCMO_DIRAC
      Case ('.NACT')
        READ(LUCMD,*) NACT_DIRAC
        write(LUPRI,*) ' read NACT_DIRAC=',NACT_DIRAC
      Case ('.NOHSDIR')
        nohsccdirac=.true. 
      Case ('.PRINT')
        READ(LUCMD,*) iprint_hscc
        write(LUPRI,*) ' read IPRINT_hscc=',iprint_hscc
      Case ('.DRYRUN')
        nohsccfiles=.true.
      Case Default
        Write(LUPRI,'(/,3A,/)') ' Keyword "',KEYWORD, &
        '" is not yet implemented in HSCCINP.'
      End Select
    EndIf

  Else
     WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,'" not recognized in HSCCINP.'
     CALL PRTAB(size(TABLE),TABLE,WORD//' input keywords',LUPRI)
     CALL QUIT('Illegal keyword in HSCCINP.')

  EndIf

EndDo

!  ... print out read parameters ... 

 WRITE (LUPRI,'(//)')
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(1X,A)') 'HSCCINP: Set-up for the external standalone HSMRCC code'
 CALL PRSYMB(LUPRI,'=',75,0)
 WRITE(LUPRI,'(/1X,A,I3)') '* General print level   : ',iprint_hscc
! WRITE(LUPRI,'( 1X,A,I3)') '*                NELE   : ',NELE

 WRITE(LUPRI,'(1X,A,I3)')  '*             NCMO_DIRAC=',NCMO_DIRAC
 WRITE(LUPRI,'(1X,A,I3)')  '*             NACT_DIRAC=',NACT_DIRAC
! WRITE(LUPRI,'(1X,A,I3)')  '*            NFRZV_DIRAC=',NFRZV_DIRAC

 if (nohsccfiles) then
  WRITE(LUPRI,'(1X,A)') '* Do not prepare HSCC files (dry run). Take them from GAMESS.'
 else
  WRITE(LUPRI,'(1X,A)') '* Prepare HSCC files by DIRAC.'
 endif

 if (nohsccdirac) then
  WRITE(LUPRI,'(1X,A)') '* Only prepare HSCC files by DIRAC, do not continue with "afterDIRAC" HSCC.'
 endif

 call set_dirac_UOUT

 call set_dirac_hscc_interface_mode

 call check_dirac_hscc_interface_mode

 if (.not.get_dirac_hscc_interface_mode()) call quit('NOT in dirac_hscc_interface_mode !')

 call set_hscc_iprint(iprint_hscc) ! call hsmrcc_external routine 

 WRITE(LUPRI,*) 'HSCC get_iprint()=',get_iprint()
 WRITE(LUPRI,*) 'HSCC variable UOUT=',UOUT
 
 print *,' returning from HSCCINP ...'

 Return

910  CONTINUE
        CALL QUIT('Error reading LUCMD, no "*HSCC" input found !')
920  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input some error occured...')
930  CONTINUE
        CALL QUIT('Error reading LUCMD, after "*HSCC" input end of file reached...')

end subroutine HSCCINP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_HSMRCC_Dirac_info
!-----------------------------------------------------------------------------
!
! Print-out basic info on the Dirac_HSMRCC interface run
!
! Called from:
!
! Written by Miro ILIAS, Tel Aviv 2007, Prievidza 2008
!
!-----------------------------------------------------------------------------
#include "priunit.h"

write(LUPRI,9001)

9001 format('                                                                 ',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '|        I N T E R F A C E                                      |',/, &
            '|                                                               |',/, &
            '| FOR THE STANDALONE PROGRAM FOR THE EFFECTIVE AND INTERMEDIATE |',/, &
            '| HAMILTONIAN HILBERT SPACE MULTI-REFERENCE COUPLED CLUSTER     |',/, &
            '| METHODS                                                       |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
            '|                                                               |',/, &
            '| Written by Miro ILIAS, Shamasunar K.R.                        |',/, &
            '|                                                               |',/, &
            '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' )

end subroutine print_HSMRCC_Dirac_info

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7
function get_moltra_index(is_hole,i) result(moi)
!
! returns MOLTRA spinor index of HSMRCC NHMO/NPMO  index i
!
use io_moltra
use Parameters
use GMSMod

implicit none
#include "priunit.h"

integer, intent(in) :: i
logical, intent(in) :: is_hole
integer :: moi

!  write(LUPRI,*) 'NCMO:',NCMO
!  write(LUPRI,*) 'NVMO:',NVMO
!  write(LUPRI,*) 'NEMO:',NEMO

if (is_hole) then ! ... is hole (1-NHMO)
 ! write(LUPRI,*) 'hole index, i=',i
 !call print_gms_data_line(i)
 moi=MOLTRA_indx(i)
else ! .. is particle (index 1-NPMO)
 ! write(LUPRI,*) 'particle index, i=',i
 ! call print_gms_data_line(i+NCMO)
  moi=MOLTRA_indx(i+NCMO)
endif

!write(LUPRI,*) 'moi=',moi

end function get_moltra_index

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_moltra_integral(i,j,k,l) result(moint)
!-----------------------------------------------------------------------------
!
!  get nonzero MOLTRA integral (IJ|KL) from the MOLTRA integral list
!
!  if none found, simply returns zero
!
! Called from hsmrcc/dir_hscc_mods.F90:make_VIND_VINT_files
!
!------------------------------------------------------------------------------
!C     Real groups       : classes 1-4 are used
!C     Complex groups    : classes 1-8 are used
!C     Quaternion groups : classes 1-16 are used
!C
!C     Class 1  : Real part of integrals      (i    j    |k    l   )
!C     Class 2  : Real part of integrals      (i    j    |kbar lbar)
!C     Class 3  : Real part of integrals      (i    jbar |k    lbar)
!C     Class 4  : Real part of integrals      (i    jbar |kbar l   )
!C     Class 5  : Imaginary part of integrals (i    j    |k    l   )
!C     Class 6  : Imaginary part of integrals (i    j    |kbar lbar)
!C     Class 7  : Imaginary part of integrals (i    jbar |k    lbar)
!C     Class 8  : Imaginary part of integrals (i    jbar |kbar l   )
!C     Class 9  : Real part of integrals      (i    j    |kbar l   )
!C     Class 10 : Real part of integrals      (i    j    |k    lbar)
!C     Class 11 : Real part of integrals      (ibar j    |k    l   )
!C     Class 12 : Real part of integrals      (i    jbar |k    l   )
!C     Class 13 : Imaginary part of integrals (i    j    |kbar l   )
!C     Class 14 : Imaginary part of integrals (i    j    |k    lbar)
!C     Class 15 : Imaginary part of integrals (ibar j    |k    l   )
!C     Class 16 : Imaginary part of integrals (i    jbar |k    l   )
!-----------------------------------------------------------------------------
use io_moltra
use Parameters
use GMSMod
implicit none
#include "priunit.h"
integer, intent(in) :: i,j,k,l ! KU indexes
integer :: i_kr,j_kr,k_kr,l_kr
integer :: ix,jx,kx,lx,ipos_ijkl,ipos_jilk,ipos_klij,ipos_lkji,int_jilk
real*8 :: moint

moint=0.0d0

!write(LUPRI,'(4x,a,4i3)') 'get_moltra_integral: (i,j,k,l):',i,j,k,l

i_kr=+1;j_kr=+1;k_kr=+1;l_kr=+1

if (i>NSP/2) i_kr=-1
if (j>NSP/2) j_kr=-1
if (k>NSP/2) k_kr=-1
if (l>NSP/2) l_kr=-1


! TODO: print-out 
!write(LUPRI,'(4x,a,a,i3,a,i2,a,i3,a,i2,a,i3,a,i2,a,i3,a,i2,a)')  & 
!          'get_moltra_integral: (i,j,k,l):',                     & 
!          '(',i,'[',i_kr,']',j,'[',j_kr,']|',k,'[',k_kr,']',l,'[',l_kr,'])'

! ... take into account Kramers symmetries of integrals ...
if (i_kr.eq.+1.and.j_kr.eq.+1.and.k_kr.eq.+1.and.l_kr.eq.+1) then      ! (aa|aa)
 ix=i; jx=j; kx=k;     lx=l
else if (i_kr.eq.+1.and.j_kr.eq.+1.and.k_kr.eq.-1.and.l_kr.eq.-1) then ! (aa|bb)
 ix=i; jx=j; kx=k-NSP/2; lx=l-NSP/2
else if (i_kr.eq.-1.and.j_kr.eq.-1.and.k_kr.eq.+1.and.l_kr.eq.+1) then ! (bb|aa)
 ix=i-NSP/2; jx=j-NSP/2; kx=k; lx=l
else if (i_kr.eq.-1.and.j_kr.eq.-1.and.k_kr.eq.-1.and.l_kr.eq.-1) then ! (bb|bb)
 ix=i-NSP/2; jx=j-NSP/2; kx=k-NSP/2; lx=l-NSP/2
else
 return ! no suitable combination of 'bars' for nonzero integral
endif

! write(LUPRI,'(4x,a,4i3)') '...new indexes, i,j,k,l:',ix,jx,kx,lx 

!
! ... take into account general symmetry properties of (IJ|KL) 
!
!  (IJ|KL) = (JI|LK)=(KL|IJ)=(LK|JI) for Coulomb integrals
!

! (IJ|KL)

  ipos_ijkl=get_nonzero_integral_position(ix,jx,kx,lx)
  ipos_jilk=get_nonzero_integral_position(jx,ix,lx,kx)
  ipos_klij=get_nonzero_integral_position(kx,lx,ix,jx)
  ipos_lkji=get_nonzero_integral_position(lx,kx,jx,ix)
 
 if (ipos_ijkl.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_ijkl)
 else if (ipos_jilk.ne.-1) then
    moint=cbuf%p_rA2D(1,int_jilk)
 else if (ipos_klij.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_klij)
 else if (ipos_lkji.ne.-1) then
    moint=cbuf%p_rA2D(1,ipos_lkji)
 endif

 return

end function get_moltra_integral


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine make_dirac_syminf
!-----------------------------------------------------------------------------
! prepare file for the gamess-like input file
!
! $GMS/"syminf.f90"
!-----------------------------------------------------------------------------
!  end subroutine make_dirac_syminf
