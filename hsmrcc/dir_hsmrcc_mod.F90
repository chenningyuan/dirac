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
use io_moltra
#include "implicit.h"

call Header('Output from "prepare_dirac_HSMRCC_VINF_file"',-1)


end subroutine prepare_dirac_HSMRCC_VINF_file

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine prepare_dirac_HSMRCC_file
!----------------------------------------------------------------------------------
!
!  prepare files with DIRAC data (orbitals, symmetry, number of electrons..etc) for the 
! standalone HSMRCC run
!
!----------------------------------------------------------------------------------
use Memory
use Vectors
use io_moltra
#include "implicit.h"
#include "priunit.h"
!character(len=22) :: dirac_file='Dirac_for_HSMRCC.txt'
!character(len=*) :: dirac_file='Dirac_HSMRCC'
!character :: dirac_file='Dirac_HSMRCC'
character(len=12) :: dirac_file='Dirac_HSMRCC'
character(len=4) :: repna_temp_a, repna_temp_b
integer :: multbs(1:8,1:8)
logical :: prepare_file=.false.

type(Vect) :: SpinOfOrb,OccGMSRef,IrrepOrb,IrpLabOrb,OrbType,OrbSpace,OrbEnergy
 
if (iprint_hscc >= 5) then
  call Header('Output from "prepare_dirac_HSMRCC_file"',-1)
endif

! ... check if memory block is present
if (.not.IsMemorySticked() ) then
 write(LUPRI,*) 'memory array not present !' 
 call Quit(': memory array not present !  ')
else
 write(LUPRI,*) ' memory sticked !' 
endif

SpinOfOrb  = MakeVector('INTE',NSP,'GAMESS-like SpinOfOrb')
OccGMSRef  = MakeVector('INTE',NSP,'GAMESS-like OccGMSRef')
IrrepOrb   = MakeVector('INTE',NSP,'GAMESS-like IrrepOrb')
IrpLabOrb  = MakeVector('INTE',NSP,'GAMESS-like IrpLabOrb')
OrbType    = MakeVector('INTE',NSP,'GAMESS-like OrbType')
OrbSpace   = MakeVector('INTE',NSP,'GAMESS-like OrbSpace')

OrbEnergy  = MakeVector('REAL',NSP,'GAMESS-like OrbEnergy')

write(LUPRI,*)

nelec_tot = 0
do irp=1,NSYMRP
  write(LUPRI,'(2x,i2,a,i2,a,i2,a)') irp,'/',NSYMRP,' -> ',NELEC(irp),' electrons'
  nelec_tot = nelec_tot + NELEC(irp)
enddo

write(LUPRI,*) 'nelec_tot=',nelec_tot
write(LUPRI,*)

numirrep_dirac = GetNBSYM()
write(LUPRI,*) 'numirrep_dirac=',numirrep_dirac

! ... fill irreps now ...
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

write(LUPRI,*)

! fill orb.en.
ii = 1
do i=1, NSP/2

! write(LUPRI,*) 'EPS1(',i,')=',EPS1(i),EPS1(i+NSP/2)
 OrbEnergy%p_rA1D(ii)   = EPS1(i) 
 OrbEnergy%p_rA1D(ii+1) = EPS1(i+NSP/2)
 
 repna_temp_a = REPNA(IREPSP(i+NSP))
 repna_temp_b = REPNA(IREPSP(i+NSP/2+NSP))

 !write(LUPRI,'(2x,i3,1x,a4,1x,a4)') i,REPNA(IREPSP(i+NSP)),REPNA(IREPSP(i+NSP/2+NSP))
 write(LUPRI,'(2x,i3,1x,a4,i2,1x,a4,i2)') i, repna_temp_a,IREPSP(i), repna_temp_b, IREPSP(i+NSP/2)

 IrrepOrb%p_iA1D(ii)    = IREPSP(i+NSP)
 read(repna_temp_a,'(a4)') IrpLabOrb%p_iA1D(ii)

! IrrepOrb%p_iA1D(ii+1)  = IREPSP(i+NSP/2+NSP)
 IrrepOrb%p_iA1D(ii+1)  = IREPSP(i+NSP) ! previously saved irrep...

 read(repna_temp_b,'(a4)') IrpLabOrb%p_iA1D(ii+1)

 SpinOfOrb%p_iA1D(ii)   = +1
 SpinOfOrb%p_iA1D(ii+1) = -1

! write(LUPRI,'(a,2i4,2f12.7)') 'orb_en:',ii,ii+1,OrbEnergy%p_rA1D(ii),OrbEnergy%p_rA1D(ii+1)

   ii = ii + 2

enddo

do i=1,NSP
 if (i<=nelec_tot) then
  OccGMSRef%p_iA1D(i) = 1
 else
  OccGMSRef%p_iA1D(i) = 0
 endif
!
 if (i<=NCMO) then
  OrbType%p_iA1D(i) = -1
  OrbSpace%p_iA1D(i) = -1
 else if (i>NCMO .and. i<=NCMO+NACT) then
  OrbType%p_iA1D(i) = 0
  OrbSpace%p_iA1D(i) = 1
 else if (i>NCMO+NACT) then
  OrbType%p_iA1D(i) = 1
  OrbSpace%p_iA1D(i) = -1
 else
  call quit('dead branch !')
 endif

 if (iprint_hscc >= 5) then
   write(LUPRI,876) i,SpinOfOrb%p_iA1D(i),OccGMSRef%p_iA1D(i),  & 
                    OrbType%p_iA1D(i),OrbSpace%p_iA1D(i),       & 
                    IrrepOrb%p_iA1D(i),IrpLabOrb%p_iA1D(i),OrbEnergy%p_rA1D(i)

   876 format(i3,1x,'spin=',i2,1x,'occ=',i2,1x,'orbtype=',i2,1x,'orbspace=',i2,1x, &
                   'irp=',i2,1x,' irplab=',a4,1x,'orb.en=',f12.7)
 endif

enddo

If (prepare_file) then

!open(ludh,file='Dirac_for_HSMRCC.txt',form='formatted',status='unknown', access='sequential') 
!open(ludh,file=dirac_file,form='formatted',status='unknown', access='sequential') 
open(ludh,file=dirac_file,form='unformatted',status='unknown', access='sequential') 
write(ludh) iprint_hscc
write(ludh) nelec_tot
write(ludh) SPFR
write(ludh) ECORE
write(ludh) NSP
write(ludh) NORBT
write(ludh) NCMO
write(ludh) NACT
write(ludh) numirrep_dirac
write(ludh) irrep_product_dirac
!write(ludh) ( SpinOfOrb%p_iA1D(i), i=1, NSP )
write(ludh) ( SpinOfOrb%p_iA1D(1:NSP) )
write(ludh) ( OccGMSRef%p_iA1D(1:NSP) )
write(ludh) ( OrbType%p_iA1D(1:NSP) )
write(ludh) ( OrbSpace%p_iA1D(1:NSP) )
write(ludh) ( IrrepOrb%p_iA1D(1:NSP) )
write(ludh) ( IrpLabOrb%p_iA1D(1:NSP) )
!write(ludh) ( OrbEnergy%p_rA1D(i), i=1, NSP )
write(ludh) ( OrbEnergy%p_rA1D(1:NSP) )
close(ludh,status="keep")

if (iprint_hscc >= 5) then
    write(LUPRI,'(/,2x,a,a)')   & 
    '... (unformatted) Dirac input file for the HSMRCC run ready under name:',dirac_file
endif

write(LUPRI,*) 
write(LUPRI,*) 'saved NELEC_TOT=',nelec_tot
write(LUPRI,*) 'saved NCORE=',NCORE
write(LUPRI,*) 'saved NACT=',NACT
write(LUPRI,*) 'saved NORBT=',NACT
write(LUPRI,*) 'saved NSP=',NSP
write(LUPRI,*) 'saved numirrep_dirac=',numirrep_dirac

endif

call DeAllocObject(OrbEnergy)
call DeAllocObject(OrbSpace)
call DeAllocObject(OrbType)
call DeAllocObject(IrpLabOrb)
call DeAllocObject(IrrepOrb)
call DeAllocObject(OccGMSRef)
call DeAllocObject(SpinOfOrb)

!call Quit('prepare_dirac_HSMRCC_file')
end subroutine prepare_dirac_HSMRCC_file


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
use Memory
use Vectors
!use io_moltra
#include "implicit.h"
#include "priunit.h"
!character(len=22) :: dirac_file='Dirac_for_HSMRCC.txt'
!character(len=*) :: dirac_file='Dirac_HSMRCC'
!character :: dirac_file='Dirac_HSMRCC'
!character(len=17) :: dirac_file_VIND='Dirac_HSMRCC.VIND'
!character(len=17) :: dirac_file_VINT='Dirac_HSMRCC.VINT'
!character(len=17) :: dirac_file_FHH='Dirac_HSMRCC.FHH'
!character(len=17) :: dirac_file_FHP='Dirac_HSMRCC.FHP'
!character(len=17) :: dirac_file_FPP='Dirac_HSMRCC.FPP'

logical :: doprepare=.false.
!logical:: doprepare=.true.


if (doprepare) then

  write(LUPRI,*) 'preparing DIRAC integral files for HSMRCC !!!'

!open(lu_vind,file=dirac_file_VIND,form='unformatted',status='unknown', access='sequential') 
!open(lu_vint,file=dirac_file_VINT,form='unformatted',status='unknown', access='sequential') 
!open(lu_fhh,file=dirac_file_FHH,form='unformatted',status='unknown', access='sequential') 
!open(lu_fhp,file=dirac_file_FHP,form='unformatted',status='unknown', access='sequential') 
!open(lu_fpp,file=dirac_file_FPP,form='unformatted',status='unknown', access='sequential') 
!close(lu_vind)
!close(lu_vint)
!close(lu_fhh)
!close(lu_fhp)
!close(lu_fpp)

endif

end subroutine prepare_dirac_HSMRCC_int_files

end module Dirac_HSMRCC
