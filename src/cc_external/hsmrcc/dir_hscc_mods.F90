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
