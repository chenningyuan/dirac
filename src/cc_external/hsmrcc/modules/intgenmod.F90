!
! Module dealing with the generation of transformed two-electron antisymmetric integrals
!
module IntGenMod

   use Parameters
   use IOModule
   use GMSMod
   use VacMod
   implicit none
   private
   integer, parameter :: ASINTBUFSIZE  = 5460
   double precision, public :: CUTOFF
   ! NZVINTBUFSIZE in IOModule
   integer         , dimension(NZVINTBUFSIZE) :: NZVPkIBuf
   double precision, dimension(NZVINTBUFSIZE) :: NZVIntBuf
   double precision, allocatable :: GenIntBuf(:)
   interface IntFileGen ! common name for 2 procedures
      module procedure IntFileGen_OneBody
      module procedure IntFileGen_TwoBody
   end interface IntFileGen

  ! ... this specification is needed for showing out internal routines 
   public :: DrvIntGen, DrvIntGen_DIRAC,NZVPkIBuf,NZVIntBuf 
   public :: check_gms_integrals, control_print_VIND_VINT, get_gms_spacial_integrals
   public :: PackRSHP,RecNumHP
 
contains

subroutine DrvIntGen
!--------------------------------------------------------------------------------
!
!         called from main F90 'HSMRCC' routine (source/hsmrcc.F90)
!
!--------------------------------------------------------------------------------
   use TimingMod
   implicit none
   integer :: StVRec
   integer :: GenIntBufSize
   integer :: i
   logical :: isMR2EINT
   character*4 :: integ_label

   ! ... Progress Report on Integral Generation ...

   write(UOUT,9011)
   write(UOUT,9012) CUTOFF
   write(UOUT,9013,advance='no')

9011 format('******************************',/,&
            'HSMRCC INTEGRAL FILE GENERATOR',/,&
            '******************************'/)  

9012 format('Integral CUTOFF is ',E13.8)
9013 format('Generating 14 Integral Files .. ')

   ! ... Allocate Integral Generation Buffer, first get it's size ...

   GenIntBufSize = MAXVAL(RecLen)
   GenIntBufSize = MAX(GenIntBufSize,RecLen(FHH)*NumRec(FHH))
   GenIntBufSize = MAX(GenIntBufSize,RecLen(FHP)*NumRec(FHP))
   GenIntBufSize = MAX(GenIntBufSize,RecLen(FPP)*NumRec(FPP))
   GenIntBufSize = 100*GenIntBufSize
   allocate(GenIntBuf(GenIntBufSize))

   if ( get_iprint() >= 7) then
     write(UOUT,'(/2x,a,i12)')  &
     'DrvIntGen: buffer for integrals allocated, size=',GenIntBufSize
   endif

   ! ... check the present of Gamess integral file 
   inquire(file=FileMR2EINT,exist=isMR2EINT)

   if (.not.isMR2EINT) then
     write(UOUT,*) 'DrvIntGen: Gamess integral file <FileMR2EINT> not present ! Error exit !'
     call flushout
     call quit('DrvIntGen: missing Gamess integral file FileMR2EINT !')
   endif

   ! ... Open Gamess integral FileMR2EINT file ...
   open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
   rewind(UMR2EINT)

   if ( get_iprint() >= 7) then
     write(UOUT,'(/2x,a,a/)')  &
     'DrvIntGen: gamess integral file open and rewinded, FileMR2EINT=',FileMR2EINT
   endif

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
   !       01. V2H2P
   !       02. VPPPP
   !       03. VHHHH
   !       04. VHHPP
   !       05. V2HPP
   !       06. V2PHH
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

   ! ... Two-Body integrals ...

   call StartTiming
   StVRec = 1 ! ... starting position of integral record...
 
   integ_label='2H2P'
   call IntFileGen(TimeTaken(:,V2H2P),RecLen(V2H2P),NumRec(V2H2P),StVRec,          &
                   StartRec(V2H2P),EndngRec(V2H2P),TotNZInt(V2H2P),NZIntRec_V2H2P, &
                   SIB_V2H2P,ResRecHP,RecNumHP,RecNumHP,                           &
                   SetQRangeHP,SetRRangeHP,SetSRangeHP,PackRSHP,integ_label)
   
   integ_label='PPPP'
   call IntFileGen(TimeTaken(:,VPPPP),RecLen(VPPPP),NumRec(VPPPP),StVRec,          &
                   StartRec(VPPPP),EndngRec(VPPPP),TotNZInt(VPPPP),NZIntRec_VPPPP, &
                   SIB_VPPPP,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangePP,SetSRangeMM,PackRSMM,integ_label)

   integ_label='HHHH'
   call IntFileGen(TimeTaken(:,VHHHH),RecLen(VHHHH),NumRec(VHHHH),StVRec,          &
                   StartRec(VHHHH),EndngRec(VHHHH),TotNZInt(VHHHH),NZIntRec_VHHHH, &
                   SIB_VHHHH,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangeHH,SetSRangeMM,PackRSMM,integ_label)

   integ_label='HHPP'
   call IntFileGen(TimeTaken(:,VHHPP),RecLen(VHHPP),NumRec(VHHPP),StVRec,          &
                   StartRec(VHHPP),EndngRec(VHHPP),TotNZInt(VHHPP),NZIntRec_VHHPP, &
                   SIB_VHHPP,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangePP,SetSRangeMM,PackRSMM,integ_label)

   integ_label='2HPP'
   call IntFileGen(TimeTaken(:,V2HPP),RecLen(V2HPP),NumRec(V2HPP),StVRec,          &
                   StartRec(V2HPP),EndngRec(V2HPP),TotNZInt(V2HPP),NZIntRec_V2HPP, &
                   SIB_V2HPP,ResRec2H,RecNum2H,RecNumMM,                           &
                   SetQRange2H,SetRRangePP,SetSRangeMM,PackRSMM,integ_label)

   integ_label='2PHH'
   call IntFileGen(TimeTaken(:,V2PHH),RecLen(V2PHH),NumRec(V2PHH),StVRec,          &
                   StartRec(V2PHH),EndngRec(V2PHH),TotNZInt(V2PHH),NZIntRec_V2PHH, &
                   SIB_V2PHH,ResRec2P,RecNum2P,RecNumMM,                           &
                   SetQRange2P,SetRRangeHH,SetSRangeMM,PackRSMM,integ_label)

   integ_label='PHHP'
   call IntFileGen(TimeTaken(:,VPHHP),RecLen(VPHHP),NumRec(VPHHP),StVRec,          &
                   StartRec(VPHHP),EndngRec(VPHHP),TotNZInt(VPHHP),NZIntRec_VPHHP, &
                   SIB_VPHHP,ResRecHP,RecNumHP,RecNumHP,                           &
                   SetQRangeHP,SetRRangeHP,SetSRangeHP,PackRSHP,integ_label)

   integ_label='PH2P'
   call IntFileGen(TimeTaken(:,VPH2P),RecLen(VPH2P),NumRec(VPH2P),StVRec,          &
                   StartRec(VPH2P),EndngRec(VPH2P),TotNZInt(VPH2P),NZIntRec_VPH2P, &
                   SIB_VPH2P,ResRec2P,RecNum2P,RecNumHP,                           &
                   SetQRange2P,SetRRangeHP,SetSRangeHP,PackRSHP,integ_label)

   integ_label='PHPP'
   call IntFileGen(TimeTaken(:,VPHPP),RecLen(VPHPP),NumRec(VPHPP),StVRec,          &
                   StartRec(VPHPP),EndngRec(VPHPP),TotNZInt(VPHPP),NZIntRec_VPHPP, &
                   SIB_VPHPP,ResRecPH,RecNumPH,RecNumMM,                           &
                   SetQRangePH,SetRRangePP,SetSRangeMM,PackRSMM,integ_label)

   integ_label='2HHP'
   call IntFileGen(TimeTaken(:,V2HHP),RecLen(V2HHP),NumRec(V2HHP),StVRec,          &
                   StartRec(V2HHP),EndngRec(V2HHP),TotNZInt(V2HHP),NZIntRec_V2HHP, &
                   SIB_V2HHP,ResRec2H,RecNum2H,RecNumHP,                           &
                   SetQRange2H,SetRRangeHP,SetSRangeHP,PackRSHP,integ_label)

   integ_label='HHHP'
   call IntFileGen(TimeTaken(:,VHHHP),RecLen(VHHHP),NumRec(VHHHP),StVRec,          &
                   StartRec(VHHHP),EndngRec(VHHHP),TotNZInt(VHHHP),NZIntRec_VHHHP, &
                   SIB_VHHHP,ResRecMM,RecNumMM,RecNumHP,                           &
                   SetQRangeMM,SetRRangeHP,SetSRangeHP,PackRSHP,integ_label)

   if ( get_iprint() >= 7) then
     write(UOUT,'(2x,a)') 'DrvIntGen: all two-el.integrals done'
   endif

   ! ... One-Body integrals ...

   call IntFileGen(TimeTaken(:,FHH),UFHH,RecLen(FHH),NumRec(FHH),  &
                   StartRec(FHH),EndngRec(FHH),TotNZInt(FHH),      &
                   SIB_FHH,RecNum2H)
   call IntFileGen(TimeTaken(:,FHP),UFHP,RecLen(FHP),NumRec(FHP),  &
                   StartRec(FHP),EndngRec(FHP),TotNZInt(FHP),      &
                   SIB_FHP,RecNumHP)
   call IntFileGen(TimeTaken(:,FPP),UFPP,RecLen(FPP),NumRec(FPP),  &
                   StartRec(FPP),EndngRec(FPP),TotNZInt(FPP),      &
                   SIB_FPP,RecNum2P)

   if ( get_iprint() >= 7) then
     write(UOUT,'(2x,a)') 'DrvIntGen: one-el.integrals done'
   endif

   ! ... Integral Generation Done ...

   write(UOUT,9021) ToTTime()
   write(UOUT,9022)
   call flushout

9021 format('Done in ',F6.1,' seconds.'/)
9022 format('**************************************',/,&
            'END OF HSMRCC INTEGRAL FILE GENERATION',/,&
            '**************************************'/)

   deallocate(GenIntBuf) ! ... deallocate big buffer memory ...

   close(UMR2EINT) ! ... close Gamess 2-el. integrals file

   return

end subroutine DrvIntGen

#include "maxmin.h"

subroutine SIB_V2H2P(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   integer :: I, J, A, B
   integer :: RecNum, BufInd

   ! For contributing to V2H2P, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <I,J||A,B>

   if ( P  > NHMO ) return
   if ( Q  > NHMO ) return
   if ( R <= NCMO ) return
   if ( S <= NCMO ) return

   I = P
   J = Q
   A = R - NCMO
   B = S - NCMO

   RecNum = (I-1)*NPMO + A   ! IA
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   ! ... calculate buffer index
   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + B  ! JB
   ! ... store the <IJ||AB> integral
   GenIntBuf(BufInd) = AsyInt

   return
end subroutine SIB_V2H2P

subroutine SIB_V2HPP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: I, J, A, B
   integer :: RecNum, BufInd

   ! For contributing to V2HPP, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <I,J||A,B>

   if ( P  > NHMO ) return

   if ( Q  > NHMO ) return

   if ( R <= NCMO ) return

   if ( S <= NCMO ) return

   I = P
   J = Q
   A = R - NCMO
   B = S - NCMO

   RecNum = (I-1)*NHMO + J  ! IJ
   if ( RecNum < StRec ) return

   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( A < B ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(A,B)
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_V2HPP

subroutine SIB_V2PHH(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: A, B, I, J
   integer :: RecNum, BufInd

   ! For contributing to V2PHH, first index
   ! pair should be particle-particle type and second
   ! index pair should be hole-hole type.
   !
   ! <P,Q||R,S> = <A,B||I,J>

   if ( P <= NCMO ) return

   if ( Q <= NCMO ) return

   if ( R  > NHMO ) return

   if ( S  > NHMO ) return

   A = P - NCMO
   B = Q - NCMO
   I = R
   J = S

   RecNum = (A-1)*NPMO + B ! AB
   if ( RecNum < StRec ) return

   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(I,J) ! IJ
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_V2PHH

subroutine SIB_VHHHH(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: I, J, K, L
   integer :: RecNum, BufInd


   ! For contributing to VHHHH, all
   ! indices should be hole type.
   !
   ! <P,Q||R,S> = <I,J||K,L>

   if ( P > NHMO ) return

   if ( Q > NHMO ) return

   if ( R > NHMO ) return

   if ( S > NHMO ) return

   I = P
   J = Q
   K = R
   L = S

   RecNum = MAXMIN(I,J) ! IJ

   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal
   IF ( K < L ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(K,L) ! KL
   GenIntBuf(BufInd) = IntVal
   
   return

end subroutine SIB_VHHHH

subroutine SIB_VHHPP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: I, J, A, B
   integer :: RecNum, BufInd


   ! For contributing to VHHPP, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <I,J||A,B>

   if ( P  > NHMO ) return

   if ( Q  > NHMO ) return

   if ( R <= NCMO ) return

   if ( S <= NCMO ) return

   I = P
   J = Q
   A = R - NCMO
   B = S - NCMO

   RecNum = MaxMin(I,J)  ! IJ
   if ( RecNum < StRec ) return

   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal
   IF ( A < B ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(A,B) ! AB
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_VHHPP

subroutine SIB_VPPPP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: A, B, C, D
   integer :: RecNum, BufInd


   ! For contributing to VPPPP, all
   ! indices should be particle type.
   !
   ! <P,Q||R,S> = <A,B||C,D>

   if ( P <= NCMO ) return

   if ( Q <= NCMO ) return

   if ( R <= NCMO ) return

   if ( S <= NCMO ) return

   A = P - NCMO
   B = Q - NCMO
   C = R - NCMO
   D = S - NCMO

   RecNum = MaxMin(A,B)
   if ( RecNum < StRec ) return

   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( A < B ) IntVal = -IntVal
   IF ( C < D ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(C,D)
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_VPPPP

subroutine SIB_VPHHP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   integer :: A, J, I, B
   integer :: RecNum, BufInd

   ! For contributing to VPHHP, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-particle type.
   !
   if ( P <= NCMO ) return
   if ( Q  > NHMO ) return
   if ( R  > NHMO ) return
   if ( S <= NCMO ) return

   A = P - NCMO
   J = Q
   I = R
   B = S - NCMO ! <PQ||RS>=<AJ||IB>

   RecNum = (I-1)*NPMO + A   ! IA
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + B ! JB
   GenIntBuf(BufInd) = AsyInt


   return

end subroutine SIB_VPHHP

subroutine SIB_V2HHP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   integer :: I, J, K, A
   integer :: RecNum, BufInd

   ! For contributing to V2HHP, first three indices
   ! should be of hole type and fourth index should
   ! be of particle type.
   !
   ! <P,Q||R,S> = <I,J||K,A>

   if ( P  > NHMO ) return
   if ( Q  > NHMO ) return
   if ( R  > NHMO ) return
   if ( S <= NCMO ) return

   I = P
   J = Q
   K = R
   A = S - NCMO

   RecNum = (I-1)*NHMO + K   ! IK
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + A ! JA
   GenIntBuf(BufInd) = AsyInt

   return
end subroutine SIB_V2HHP

subroutine SIB_VHHHP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: I, A, J, K
   integer :: RecNum, BufInd

   ! For contributing to VHHHP, first three indices
   ! should be of hole type and fourth index should
   ! be of particle type.
   !
   ! <P,Q||R,S> = <I,J||K,A>

   if ( P  > NHMO ) return
   if ( Q  > NHMO ) return
   if ( R  > NHMO ) return
   if ( S <= NCMO ) return

   I = P
   J = Q
   K = R
   A = S - NCMO ! substract occup.shells

   RecNum = MAXMIN(I,J) ! IJ
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + (K-1)*NPMO + A ! KA
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_VHHHP

subroutine SIB_VPH2P(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   integer :: A, I, B, C
   integer :: RecNum, BufInd

   ! For contributing to VPH2P, first, third & fourth
   ! indices should be of particle type and fourth
   ! index should be of particle type.
   !
   if ( P <= NCMO ) return
   if ( Q  > NHMO ) return
   if ( R <= NCMO ) return
   if ( S <= NCMO ) return

   A = P - NCMO
   I = Q
   B = R - NCMO
   C = S - NCMO !  <P,Q||R,S> = <A,I||B,C>

   RecNum = (A-1)*NPMO + B ! AB
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (I-1)*NPMO + C ! IC
   GenIntBuf(BufInd) = AsyInt

   return

end subroutine SIB_VPH2P

subroutine SIB_VPHPP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   double precision :: IntVal
   integer :: A, I, B, C
   integer :: RecNum, BufInd

   ! For contributing to VPHPP, second index should
   ! be hole type and first, third & fourth indices
   ! should be of particle type.
   !
   ! <P,Q||R,S> = <A,I||B,C>

   if ( P <= NCMO ) return
   if ( Q  > NHMO ) return
   if ( R <= NCMO ) return
   if ( S <= NCMO ) return

   A = P - NCMO
   I = Q
   B = R - NCMO
   C = S - NCMO

   RecNum = (A-1)*NHMO + I ! AI
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( B < C ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(B,C) ! BC
   GenIntBuf(BufInd) = IntVal

   return

end subroutine SIB_VPHPP

subroutine SIB_FHH(NumRec,RecLen,P,Q,R,S,VInt)

   implicit none
   integer, intent(in) :: NumRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: VInt
   integer :: I, J, Lamda
   integer :: RecNum, BufInd

   if ( P > NHMO ) return
   if ( R > NHMO ) return
   if ( Q /= S )   return

   ! <PQ||RQ> , Q = S CASE

   if ( Q <= NCMO ) RETURN
   if ( Q >  NHMO ) RETURN

   ! <P,Lamda||R,Lamda> = <I,Lamda||J,Lamda> CASE

   I = P
   J = R
   Lamda = Q - NCMO

   RecNum = I

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum-1)*RecLen + (J-1)*NVMO + Lamda
   GenIntBuf(BufInd) = VInt

   return

end subroutine SIB_FHH

subroutine SIB_FHP(NumRec,RecLen,P,Q,R,S,VInt)

   implicit none
   integer, intent(in) :: NumRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: VInt
   integer :: I, A, Lamda
   integer :: RecNum, BufInd

   if ( P  > NHMO ) return
   if ( R <= NCMO ) return
   if ( Q /= S )   return

   ! <PQ||RQ> , Q = S CASE

   if ( Q <= NCMO ) RETURN
   if ( Q >  NHMO ) RETURN

   ! <P,Lamda||R,Lamda> = <I,Lamda||J,Lamda> CASE

   I = P
   A = R - NCMO
   Lamda = Q - NCMO

   RecNum = I

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum-1)*RecLen + (A-1)*NVMO + Lamda
   GenIntBuf(BufInd) = VInt

   return

end subroutine SIB_FHP

subroutine SIB_FPP(NumRec,RecLen,P,Q,R,S,VInt)

   implicit none
   integer, intent(in) :: NumRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: VInt
   integer :: A, B, Lamda
   integer :: RecNum, BufInd

   if ( P <= NCMO ) return
   if ( R <= NCMO ) return
   if ( Q /= S )   return

   ! <PQ||RQ> , Q = S CASE

   if ( Q <= NCMO ) RETURN
   if ( Q >  NHMO ) RETURN

   ! <P,Lamda||R,Lamda> = <I,Lamda||J,Lamda> CASE

   A = P - NCMO
   B = R - NCMO
   Lamda = Q - NCMO

   RecNum = A

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum-1)*RecLen + (B-1)*NVMO + Lamda
   GenIntBuf(BufInd) = VInt

   return

end subroutine SIB_FPP

function FINTGET(PQ,vac)

   implicit none
   double precision :: FINTGET
   integer, intent(in) :: PQ, vac
   integer :: Lamda

   FIntGet = dble(0) 

   do Lamda = 1, NVMO
      if (IsActHole(Lamda,vac)) FINTGET = FINTGET + GenIntBuf(PQ+Lamda)
      if (IsActHole(Lamda,1  )) FINTGET = FINTGET - GenIntBuf(PQ+Lamda)
   enddo

   return

end function FINTGET

subroutine ResRecMM(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! Perform Reverse of MaxMin
   ! Get P and Q from MAXMIN(P,Q)
   ! P,Q are either hole-hole (HH) or particle-particle indexes (PP)

   P = 1
   do
     P = P + 1
     Q = PQ - (((P-2)*(P-1))/2)
     if (Q<P) exit
   enddo

   return

end subroutine ResRecMM

subroutine ResRecHP(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NPMO + Q
   ! finding P & Q from PQ, where P is hole(H), Q particle (P) index

   P = (PQ-1)/NPMO + 1
   Q = PQ - (P-1)*NPMO

   return

end subroutine ResRecHP

subroutine ResRecPH(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NHMO + Q
   ! finding P & Q from PQ, where P is particle and Q hole index

   P = (PQ-1)/NHMO + 1
   Q = PQ - (P-1)*NHMO

   return

end subroutine ResRecPH

subroutine ResRec2H(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NHMO + Q
   ! finding P & Q from PQ, where both P,Q are hole indexes

   P = (PQ-1)/NHMO + 1
   Q = PQ - (P-1)*NHMO

   return

end subroutine ResRec2H

subroutine ResRec2P(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NPMO + Q
   ! finding P & Q from PQ. Both P,Q are particle type indexes (P)

   P = (PQ-1)/NPMO + 1
   Q = PQ - (P-1)*NPMO

   return

end subroutine ResRec2P

function RecNumMM(P,Q)
   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNumMM

   !if (P==Q) stop 'Error in RecNumMM'
   if (P==Q) then
     write(UOUT,*) 'RecNumMM: P,Q:',P,Q
     call quit('P ==Q => Error in RecNumMM!')
   endif
   RecNumMM = MAXMIN(P,Q)
   return
end function RecNumMM

function RecNumHP(P,Q)
! returns position of PQ common index; P is hole index, Q particle index
   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNumHP

   RecNumHP = (P-1)*NPMO + Q

   return
end function RecNumHP

function RecNumPH(P,Q)
   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNumPH

   RecNumPH = (P-1)*NHMO + Q

   return
end function RecNumPH

function RecNum2H(P,Q)
   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNum2H

   RecNum2H = (P-1)*NHMO + Q

   return
end function RecNum2H

function RecNum2P(P,Q)
   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNum2P

   RecNum2P = (P-1)*NPMO + Q

   return
end function RecNum2P

subroutine SetQRangeMM(P,SP,SQ,EP,EQ,QL,QH)
   implicit none
   integer, intent(in)  :: P, SP, SQ, EP, EQ
   integer, intent(out) :: QL, QH

   QL = 1
   QH = (P-1)

   !if (QH<QL) stop 'Error in SetQRangeMM'
   if (QH<QL) then 
     write(UOUT,*) 'Error:QH<QL!; QH,QL=',QH,QL
     Call Quit('QH<QL...error in SetQRangeMM!')
   endif

   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   return
end subroutine SetQRangeMM

subroutine SetQRangeHP(P,SP,SQ,EP,EQ,QL,QH)

   implicit none
   integer, intent(in)  :: P, SP, SQ, EP, EQ
   integer, intent(out) :: QL, QH

   QL = 1
   QH = NPMO
   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   return
end subroutine SetQRangeHP

subroutine SetQRangePH(P,SP,SQ,EP,EQ,QL,QH)

   implicit none
   integer, intent(in)  :: P, SP, SQ, EP, EQ
   integer, intent(out) :: QL, QH

   QL = 1
   QH = NHMO
   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   return

end subroutine SetQRangePH

subroutine SetQRange2H(P,SP,SQ,EP,EQ,QL,QH)

   implicit none
   integer, intent(in)  :: P, SP, SQ, EP, EQ
   integer, intent(out) :: QL, QH

   QL = 1
   QH = NHMO
   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   return

end subroutine SetQRange2H

subroutine SetQRange2P(P,SP,SQ,EP,EQ,QL,QH)

   implicit none
   integer, intent(in)  :: P, SP, SQ, EP, EQ
   integer, intent(out) :: QL, QH

   QL = 1
   QH = NPMO
   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   return

end subroutine SetQRange2P

subroutine SetRRangeHH(RL,RH)

   implicit none
   integer, intent(out) :: RL, RH

   RL = 2
   RH = NHMO

   return

end subroutine SetRRangeHH

subroutine SetRRangePP(RL,RH)

   implicit none
   integer, intent(out) :: RL, RH

   RL = 2
   RH = NPMO

   return

end subroutine SetRRangePP

subroutine SetRRangeHP(RL,RH)

   implicit none
   integer, intent(out) :: RL, RH

   RL = 1
   RH = NHMO

   return

end subroutine SetRRangeHP

subroutine SetSRangeMM(R,SL,SH)

   implicit none
   integer, intent(in)  :: R
   integer, intent(out) :: SL, SH

   SL = 1
   SH = (R-1)

   if (SH<SL) call quit('SH<SL => Error in SetSRangeMM!')

   return

end subroutine SetSRangeMM

subroutine SetSRangeHP(R,SL,SH)

   implicit none
   integer, intent(in)  :: R
   integer, intent(out) :: SL, SH

   SL = 1
   SH = NPMO

   return

end subroutine SetSRangeHP

function PackRSMM(R,S)

   implicit none
   integer, intent(in) :: R, S
   integer :: PackRSMM

   if (R==S) call quit('R==S => Error in PackRSMM')

   PackRSMM = MAXMIN(R,S)

   return

end function PackRSMM

function PackRSHP(R,S)

   implicit none
   integer, intent(in) :: R, S
   integer :: PackRSHP

   PackRSHP = (R-1)*NPMO+S

   return

end function PackRSHP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,from_where)
!----------------------------------------------------------------------------
!
!     check P,Q,R,S,LABELSIZE,LABEL,NPACK parameters of Gamess integrals
!
!  text string "from_where" specifies exact place of call in the code
!
!  Written by Miro ILIAS, Prievidza, 2008
!
!----------------------------------------------------------------------------
integer :: P, Q, R, S
integer :: LABELSIZE, LABEL, NPACK
character(len=*) :: from_where

if (P==0 .or. Q==0 .or. R==0 .or. S==0) then
 write(UOUT,'(/,2x,a)') '*** check_PQRS error catch ***'
 write(UOUT,'(2x,a,a)')  &
 'entering text label:',from_where
 write(UOUT,'(2x,3(a,i5))')  &
 'entering LABELSIZE=',LABELSIZE,'  NPACK=',NPACK,'  LABEL=',LABEL
 if (LABEL == 0) then
   write(UOUT,'(2x,a)') 'LABEL .eq. 0 ! wrong IXPREAD(:) before !'
 endif
 write(UOUT,'(2x,4(a,i2))') '...thus P=',P,' Q=',Q,' R=',R,' S=',S
 call Quit('check_PQRS: wrong indexes P,Q,R,S !')
endif

end subroutine check_PQRS


subroutine  written_int_ind(NxtRec)
!
! prints out what was written into HSMRCC integral file
!
integer, intent(in) :: NxtRec ! # of record
integer :: i

if (get_iprint() < 20) return ! returns due to print level

write(UOUT,'(/,2x,i4,a)') NxtRec,'.  block written:'
!do i = 1, NZVINTBUFSIZE
! write(UOUT,'(2x,i,a,f16.8,a,i8)') i,' integral:', NZVIntBuf(i),' packed index:',NZVPkIBuf(i) 
!enddo

end subroutine  written_int_ind

subroutine IntFileGen_TwoBody(TimeTaken,RecLen,NumRec,StVRec,              &
                              StartRec,EndngRec,TotNZInt,NZIntRec,         &
                              SetIntBuf,ResolvRec,RecordNumPQ,RecordNumRS, &
                              SetQRange,SetRRange,SetSRange,PackRS,integ_label)
!--------------------------------------------------------------------------------------------------------
!
!  routine is repeatingly reading Gamess integral file  FileMR2EINT (TRIM(IFile)//'.F296'; lu=UMR2EINT)
!  and is preparing 2-electron integral/indexes file for the HSMRCC run :
!
!  FileVIND(TRIM(IFile)//'.VIND'; lu=UVIND) & FileVINT(TRIM(IFile)//'.VINT'; lu=UVINT)
!  These are to be produced by DIRAC !
!
!  Entering: RecLen
!            NumRec
!
!            RecordNumPQ
!            RecordNumRS
!            PackRS
!            SetIntBuf - as SIB_Vxxxx (xxxx=2H2P,PPPP...all integral types)
!
!  
!  SetIntBuf - SIB_Vxxxx routine which places integral on the buffer array  
!  
!   
!  Called from: DrvIntGen ("intgenmod.F90")  
!
!---------------------------------------------------------------------------------------------------------
   implicit none
   ! ... Dummy Argument Variables ...
   real,    intent(inout) :: TimeTaken(2)
   integer, intent(in)    :: RecLen, NumRec
   integer, intent(out)   :: StartRec, EndngRec, TotNZInt
   integer, intent(inout) :: StVRec
   integer, intent(out)   :: NZIntRec(NumRec)
   character*4, intent(in) :: integ_label
   ! ... Dummy Procedure and Function Interfaces
   interface
      ! ... called as SIB_Vxxxx 
      subroutine SetIntBuf(St,En,RL,P,Q,R,S,VInt)
         implicit none
         integer, intent(in) :: St, En, RL
         integer, intent(in) :: P, Q, R, S
         double precision, intent(in) :: VInt
      end subroutine SetIntBuf
      subroutine ResolvRec(Record,P,Q)
         implicit none
         integer, intent(in)  :: Record
         integer, intent(out) :: P, Q
      end subroutine ResolvRec
      function RecordNumPQ(P,Q)
         implicit none
         integer, intent(in) :: P, Q
         integer :: RecordNumPQ
      end function RecordNumPQ
      function RecordNumRS(R,S)
         implicit none
         integer, intent(in) :: R, S
         integer :: RecordNumRS
      end function RecordNumRS
      subroutine SetQRange(P,SP,SQ,EP,EQ,QL,QH)
         implicit none
         integer, intent(in)  :: P, SP, SQ, EP, EQ
         integer, intent(out) :: QL, QH
      end subroutine SetQRange
      subroutine SetRRange(RL,RH)
         implicit none
         integer, intent(out) :: RL, RH
      end subroutine SetRRange
      subroutine SetSRange(R,SL,SH)
         implicit none
         integer, intent(in)  :: R
         integer, intent(out) :: SL, SH
      end subroutine SetSRange
      function PackRS(R,S)
         implicit none
         integer, intent(in) :: R, S
         integer :: PackRS
      end function PackRS
   end interface
   ! ... Local Variables ...
   integer :: P, Q, R, S, R1, S1
   integer :: StP, EnP, StQ, EnQ
   integer :: QLow, QHig
   integer :: RLow, RHig
   integer :: SLow, SHig
   integer :: PQ, RS, PQRS
   integer :: ASBufInd, NZVIntBufInd
   integer :: StRec, EnRec, NxtRec
   integer :: BufSize, NumRecBuf, NumBlocks, IndBlock
   integer :: NumIntRec, EnPkInd
   double precision :: ASInt
   integer :: IM56=-56,IM48=-48,IM40=-40,IM32=-32,IM16=-16

   integer NXX
   integer NPACKIC, NPACK, LABEL
   real :: dtim, dtime, tarray(2)
!   external :: dtime

   if (get_iprint() >= 5) then
     call header('Output from IntFileGen_TwoBody',-1)
     write(UOUT,'(2x,4(a,i7))') '...StVRec=',StVRec,'...StartRec=',StartRec,  & 
                                '...EndngRec=',EndngRec,'...RecLen=',RecLen
     write(UOUT,'(5x,a,a)') 'processing integral type:',integ_label
   endif

   dtim = dtime(tarray)

   ! Calculate number of records that fit
   ! into the buffer.
   BufSize = SIZE(GenIntBuf)
   if ( BufSize < RecLen ) then
      write(UOUT,*) 'BufSize >= RecLen!  BufSize,RecLen:',BufSize, RecLen 
      call quit('IntFileGen_TwoBody: Buffer must be big enough to contain a record !')
   endif
   NumRecBuf = BufSize/RecLen
   NumBlocks = NumRec/NumRecBuf
   if ( NumBlocks*NumRecBuf /= NumRec ) NumBlocks = NumBlocks + 1

   StartRec = StVRec
   NxtRec   = StVRec
   TotNZInt = 0
   NZVIntBuf    = DBLE(0)
   NZVPkIBuf    = 0
   NZVIntBufInd = 0

   if (.not.allocated(XXPREAD)) allocate(XXPREAD(NINTMX))
   if (.not.allocated(IXPREAD)) allocate(IXPREAD(NINTMX))

   dtim = dtime(tarray)

   do IndBlock = 1, NumBlocks

      if (get_iprint() >= 20) then
        write(UOUT,'(2x,a,i4,a,i4)') 'IntFileGen_TwoBody: IndBlock/NumBlocks=',  &
        IndBlock,'/',NumBlocks
      endif

      StRec = (IndBlock-1)*NumRecBuf + 1
      EnRec = MIN(IndBlock*NumRecBuf,NumRec)

      ! Construct the integrals for the current block.
      ! Initialize the buffer before constructing each block.

      GenIntBuf = DBLE(0)

      ! Rewind MR2EINT file and
      ! Read first three records

      rewind(UMR2EINT)
      read(UMR2EINT)
      read(UMR2EINT)
      read(UMR2EINT)

      NPACKIC = 0
      do
         call PREAD(LABELSIZE,INTSFP,UMR2EINT,XXPREAD,IXPREAD,NXX,NINTMX)
         NumIntRec = NXX
         IF (LABELSIZE.EQ.2) THEN
            do ASBufInd = 1, IABS(NumIntRec)
               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(2*NPACK-1)
               P = ISHFT( LABEL, -16 )
               Q = IAND(  LABEL, 65535 )
               LABEL = IXPREAD(2*NPACK)
               R = ISHFT( LABEL, -16 )
               S = IAND(  LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, aa ')
#else
               LABEL = IXPREAD(NPACK)
              ! P =       ISHFT( LABEL, -48 )
               P =       ISHFT( LABEL, IM48 )
              ! Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               Q = IAND( ISHFT( LABEL, IM48 ), 65535 )
              ! R = IAND( ISHFT( LABEL, -16 ), 65535 )
               R = IAND( ISHFT( LABEL, IM16 ), 65535 )
               S = IAND( LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, ab')
#endif
               ! Now all permutations are
               ! generated & added to buffer.
               call SetIntBuf(StRec,EnRec,RecLen,P,Q,R,S,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,P,Q,S,R,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,Q,P,R,S,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,Q,P,S,R,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,R,S,P,Q,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,R,S,Q,P,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,S,R,P,Q,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,S,R,Q,P,+ASInt)
            enddo
         ELSE IF (LABELSIZE.EQ.1) THEN
            do ASBufInd = 1, IABS(NumIntRec)
               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(NPACK)
               P =       ISHFT( LABEL, -24 )
               Q = IAND( ISHFT( LABEL, -16 ), 255 )
               R = IAND( ISHFT( LABEL,  -8 ), 255 )
               S = IAND(        LABEL,        255 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, ba')
#else
               IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                  LABEL = IXPREAD(NPACK/2)
                  P = IAND( ISHFT( LABEL, -24 ), 255 )
                  Q = IAND( ISHFT( LABEL, -16 ), 255 )
                  R = IAND( ISHFT( LABEL,  -8 ), 255 )
                  S = IAND( LABEL, 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bb')
               ELSE
                  LABEL = IXPREAD((NPACK/2)+1)
                  !P = ISHFT( LABEL, -56 )
                  P = ISHFT( LABEL, IM56 )
                  !Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  Q = IAND( ISHFT( LABEL, IM48 ), 255 )
                  !R = IAND( ISHFT( LABEL, -40 ), 255 )
                  R = IAND( ISHFT( LABEL, IM40 ), 255 )
                 ! S = IAND( ISHFT( LABEL, -32 ), 255 )
                  S = IAND( ISHFT( LABEL, IM32 ), 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bc')  
               ENDIF
#endif
               ! Now all permutations are
               ! generated & added to buffer.
               call SetIntBuf(StRec,EnRec,RecLen,P,Q,R,S,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,P,Q,S,R,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,Q,P,R,S,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,Q,P,S,R,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,R,S,P,Q,+ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,R,S,Q,P,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,S,R,P,Q,-ASInt)
               call SetIntBuf(StRec,EnRec,RecLen,S,R,Q,P,+ASInt)
            enddo
         ENDIF
         if (get_iprint() >= 22) then
            write(UOUT,'(2x,a,i8)')   &
            'IntFileGen_TwoBody: read NumIntRec=',NumIntRec
            if (NumIntRec <= 0) write(UOUT,'(2x,a)') &
            '...negative, going to exit cycle of reading Gamess file !'
            call flushout
         endif
         if (NumIntRec <= 0) exit
      enddo

      !---------------------------------------------------------
      ! Now the buffer (GenIntBuf) is completely
      ! constructed. Collect NON-ZERO
      ! integrals and their packing
      ! indices into buffer and write
      ! to the disk.
      !---------------------------------------------------------

      call ResolvRec(StRec,StP,StQ)
      call ResolvRec(EnRec,EnP,EnQ)

      do P = StP, EnP
         call SetQRange(P,StP,StQ,EnP,EnQ,QLow,QHig)
         do Q = QLow, QHig
            PQ = RecordNumPQ(P,Q)
            NZIntRec(PQ) = 0
            call SetRRange(RLow,RHig)
            do R = RLow, RHig
               call SetSRange(R,SLow,SHig)
               do S = SLow, SHig
                  RS = RecordNumRS(R,S)
                  PQRS = (PQ - StRec)*RecLen + RS
                  !  ... write if integral is larger than threshold ...
                  if ( DABS(GenIntBuf(PQRS)) > CUTOFF ) then
                     ! ... increment index - array position
                     NZVIntBufInd = NZVIntBufInd + 1
                     ! ... increment counter of nonzero, PQ-based integrals (first two-indexes) 
                     NZIntRec(PQ) = NZIntRec(PQ) + 1

                     ! ... write block of integrals/indices of NZVINTBUFSIZE parameters
                     if ( NZVIntBufInd > NZVINTBUFSIZE ) then
                        ! ...  write buffer with integrals 
                        write(UVINT,rec=NxtRec) NZVIntBuf
                        ! ...  write buffer with indexes 
                        write(UVIND,rec=NxtRec) NZVPkIBuf

                        NxtRec = NxtRec + 1 ! shift to the next rectord
                        NZVIntBufInd = 1    ! reinitialize # of integrals per buffer 
                        NZVPkIBuf = 0       ! zero buffer with packed indexes
                        NZVIntBuf = DBLE(0) ! zero buffer with integrals
                        if (get_iprint() >= 40) then
                         write(UOUT,'(4X,A)') '...buffer flushed out 1'
                        endif
                     endif

                   ! ... add nonzero integral into integral-buffer
                     NZVIntBuf(NZVIntBufInd) = GenIntBuf(PQRS)
                   ! ... add packed (R,S)-indexes into index-buffer
                     NZVPkIBuf(NZVIntBufInd) = PackRS(R,S)

                     if (get_iprint() >= 40) then
                       call print_integral(P,Q,R,S,PQ,RS,PackRS,NZVIntBufInd,  & 
                          NZVIntBuf(NZVIntBufInd), NZVPkIBuf(NZVIntBufInd),integ_label)
                     endif

                  endif

               enddo ! end of S
            enddo ! end of R
            ! ... increase counter of nonzero integrals of PQ index
            TotNZInt = TotNZInt + NZIntRec(PQ) 

            if (get_iprint() >= 30) then
             if (NZIntRec(PQ).ne.0) then
              write(UOUT,'(1x,3(a,i5),a,i8)')     & 
              'number of nonzero <PQ||RS> integrals for P=',P,',Q=',Q,  &
              ' NZIntRec(PQ=',PQ,')=',NZIntRec(PQ)
             endif
            endif

         enddo ! end of Q
      enddo ! end of P
   enddo

   ! End of all blocks

   EndngRec = NxtRec
   if ( TotNZInt /= 0 ) then
      write(UVINT,rec=NxtRec) NZVIntBuf
      write(UVIND,rec=NxtRec) NZVPkIBuf
      StVRec = EndngRec + 1
      if (get_iprint() >= 40) then
           write(UOUT,'(4X,A)') '...buffer flushed out 2'
      endif
   endif
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   if (get_iprint() >=8) then
     write(UOUT,'(2x,a)') 'IntFileGen_TwoBody over...'
     write(UOUT,'(2x,2(a,i8))') '...# of nonzero integrals :',TotNZInt, '...set StVRec= ',StVRec
   endif

   return

end subroutine IntFileGen_TwoBody

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_integral(P,Q,R,S,PQ,RS,PackRS,ArrayPos, Integral, IntIndex, integ_label)
!-----------------------------------------------------------------------------------------------------
!
!     print the integral with proper indexes
!
!    ArrayPos : 1-NZVINTBUFSIZE
!
!-----------------------------------------------------------------------------------------------------
   implicit none
   integer, intent(in)    :: P,Q,R,S, PQ, RS
   integer, intent(in)    :: IntIndex, ArrayPos
   integer    :: R1, S1, I,J,K,L,A,B,C,D
   double precision, intent(in) :: Integral
   character*4, intent(in) :: integ_label

   interface
      function PackRS(R,S)
         implicit none
         integer, intent(in) :: R, S
         integer :: PackRS
      end function PackRS
   end interface

   if (.not.PackRS(R,S)==IntIndex) then
     write(UOUT,*) 'print_integral: RS=',RS
     write(UOUT,*) 'print_integral: R,S=',R,S
     call quit('print_integral: RS different from IntIndex !')
   endif

! -----------------------------------------------------------
!             ... print-out own integrals
! -----------------------------------------------------------
   selectcase (integ_label)

    case('2H2P') !  <PQ||RS>=<IJ||AB>
 
       call ResRecHP(PQ,I,A) ! P,Q = I,A
       call ResRecHP(RS,J,B) ! R,S = J,B

       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')        & 
         ArrayPos,'. V2H2P-1,1324,<IJ||AB>=<',I,' ',J,'||',A,' ',B,'>=',Integral  

       if (P.ne.I.or.Q.ne.A.or.R.ne.J.or.S.ne.B) then
         write(UOUT,'(8x,a,4i4)') '---  I,A,J,B:',I,A,J,B
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('2H2P: no correspondence of indexes !')
       endif

    case('PPPP')  ! <PQ||RS>=<AB||CD>

       call ResRecMM(PQ,A,B) ! PQ is A,B 
       call ResRecMM(RS,C,D) ! RS is C,D

       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VPPPP-2,1234,<AB||CD>=<',A,' ',B,'||',C,' ',D,'>=',Integral  

       if (P.ne.A.or.Q.ne.B.or.R.ne.C.or.S.ne.D) then
         write(UOUT,'(8x,a,4i4)') '---  A,B,C,D:',A,B,C,D
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VPPPP: no correspondence of indexes !')
       endif

    case('HHHH') ! <PQ||RS>=<IJ||KL>

       call ResRecMM(PQ,I,J) ! PQ is I,J
       call ResRecMM(RS,K,L) ! RS is K,L
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VHHHH-3,1234, <IJ||KL>=<',I,' ',J,'||',K,' ',L,'>=',Integral  

       if (P.ne.I.or.Q.ne.J.or.R.ne.K.or.S.ne.L) then
         write(UOUT,'(8x,a,4i4)') '---  I,J,K,L:',I,J,K,L
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VHHHH: no correspondence of indexes !')
       endif

    case('HHPP') ! <PQ||RS>=<IJ||AB>

       call ResRecMM(PQ,I,J) ! PQ is I,J
       call ResRecMM(RS,A,B) ! RS is A,B
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VHHPP-4,1234,<IJ||AB>=<',I,' ',J,'||',A,' ',B,'>=',Integral 

       if (P.ne.I.or.Q.ne.J.or.R.ne.A.or.S.ne.B) then
         write(UOUT,'(8x,a,4i4)') '---  I,J,A,B:',I,J,A,B
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VHHHH: no correspondence of indexes !')
       endif

    case('2HPP') ! <PQ||RS>=<IJ||AB>

       call ResRec2H(PQ,I,J) ! PQ is I,J
       call ResRecMM(RS,A,B) ! RS is A,B
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. V2HPP-5,1234, <IJ||AB>=<',I,' ',J,'||',A,' ',B,'>=',Integral  

       if (P.ne.I.or.Q.ne.J.or.R.ne.A.or.S.ne.B) then
         write(UOUT,'(8x,a,4i4)') '---  I,J,A,B:',I,J,A,B
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('V2HPP: no correspondence of indexes !')
       endif

    case('2PHH')  ! <PQ||RS>=<A,B||I,J>

       call ResRec2P(PQ,A,B) ! PQ is A,B
       call ResRecMM(RS,I,J) ! RS is I,J
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. V2PHH-6,1234,<AB||IJ>=<',A,' ',B,'||',I,' ',J,'>=',Integral 

       if (P.ne.A.or.Q.ne.B.or.R.ne.I.or.S.ne.J) then
         write(UOUT,'(8x,a,4i4)') '---  A,B,I,J:',A,B,I,J
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('V2PHH: no correspondence of indexes !')
       endif

    case('PHHP') ! <PQ||RS>=<AJ||IB>

       call ResRecHP(PQ,I,A) ! PQ is I,A
       call ResRecHP(RS,J,B) ! RS is J,B
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VPHHP-7,1324, <AJ||IB>=<',A,' ',J,'||',I,' ',B,'>=',Integral  

       if (P.ne.I.or.Q.ne.A.or.R.ne.J.or.S.ne.B) then
         write(UOUT,'(8x,a,4i4)') '---  A,I,J,B:',I,A,J,B
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VPHHP: no correspondence of indexes !')
       endif

    case('PH2P') !<PQ||RS>=<AI||BC>

       call ResRec2P(PQ,A,B) ! PQ is A,B
       call ResRecHP(RS,I,C) ! RS is I,C
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VPH2P-8,1324, <AI||BC>=<',A,' ',I,'||',B,' ',C,'>=',Integral 

       if (P.ne.A.or.Q.ne.B.or.R.ne.I.or.S.ne.C) then
         write(UOUT,'(8x,a,4i4)') '---  A,B,I,C:',A,B,I,C
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VPH2P: no correspondence of indexes !')
       endif

    case('PHPP') ! <PQ||RS>=<AI||BC>

       call ResRecPH(PQ,A,I) ! PQ is A,I
       call ResRecMM(RS,B,C) ! RS is B,C
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VPHPP-9,1234,<AI||BC>=<',A,' ',I,'||',B,' ',C,'>=',Integral  

       if (P.ne.A.or.Q.ne.I.or.R.ne.B.or.S.ne.C) then
         write(UOUT,'(8x,a,4i4)') '---  A,I,B,C:',A,I,B,C
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VPHPP: no correspondence of indexes !')
       endif


    case('2HHP') ! <P,Q||R,S>=<I,J||K,A>

       call ResRec2H(PQ,I,K) ! PQ is I,K
       call ResRecHP(RS,J,A) ! RS is J,A
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. V2HHP-10,1324,<IJ||KA>=<',I,' ',J,'||',K,' ',A,'>=',Integral 

       if (P.ne.I.or.Q.ne.K.or.R.ne.J.or.S.ne.A) then
         write(UOUT,'(8x,a,4i4)') '---  I,K,J,A:',I,K,J,A
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('V2HHP: no correspondence of indexes !')
       endif

    case('HHHP') ! <P,Q||R,S> = <I,J||K,A>, 1234

       call ResRecMM(PQ,I,J) ! PQ is I,J
       call ResRecHP(RS,K,A) ! RS is K,A 
       write(UOUT,'(2x,i4,a,i3,a,i3,a,i3,a,i3,a,f16.10)')  & 
         ArrayPos,'. VHHHP-11,1234,<IJ||KA>=<',I,' ',J,'||',K,' ',A,'>=',Integral  

       if (P.ne.I.or.Q.ne.J.or.R.ne.K.or.S.ne.A) then
         write(UOUT,'(8x,a,4i4)') '---  I,J,K,A:',I,J,K,A
         write(UOUT,'(8x,a,4i4)') '---  P,Q,R,S:',P,Q,R,S
         call quit('VHHHP: no correspondence of indexes !')
       endif
    
    case default
        write(UOUT,*) 'nothing, integ_label=',integ_label
        call quit('print_integral - error brach !')
   end select

end subroutine print_integral

subroutine IntFileGen_OneBody(TimeTaken,UF,RecLen,NumRec, &
                              StartRec,EndngRec,TotNZInt, &
                              SetIntBuf,RecordNumPQ)

   implicit none
   ! ... Dummy Arguments ...
   real, intent(inout)  :: TimeTaken(2)
   integer, intent(in)  :: UF, RecLen, NumRec
   integer, intent(out) :: StartRec, EndngRec, TotNZInt
   ! ... Interface for Dummy Procedures ...
   interface
      subroutine SetIntBuf(NR,RL,P,Q,R,S,VInt)
         implicit none
         integer, intent(in) :: NR, RL
         integer, intent(in) :: P, Q, R, S
         double precision, intent(in) :: VInt
      end subroutine SetIntBuf
      function RecordNumPQ(P,Q)
         implicit none
         integer, intent(in) :: P, Q
         integer :: RecordNumPQ
      end function RecordNumPQ
   end interface
   ! ... Local variables ...
   integer :: vac
   integer :: NumIntRec, EnPkInd
   integer :: P, Q, R, S
   integer :: PQ, IPQ
   integer :: ASBufInd, BufSize
   double precision :: ASInt, FInt
   integer :: IM56=-56,IM48=-48,IM40=-40,IM32=-32,IM16=-16

   integer NXX
   integer NPACKIC, NPACK, LABEL
   real :: dtim, dtime, tarray(2)
!   external :: dtime
   ! ... Automatic Array ...
   double precision :: IntF(NumRec*(RecLen/NVMO))

   dtim = dtime(tarray)
   BufSize = SIZE(GenIntBuf)
   if (BufSize < NumRec*RecLen) then
      stop 'Buffer too small'
   endif
   GenIntBuf(1:BufSize) = DBLE(0)

   ! Rewind MR2EINT file and
   ! Read first three records

   rewind(UMR2EINT)
   read(UMR2EINT)
   read(UMR2EINT)
   read(UMR2EINT)

   NPACKIC = 0
   do
      call PREAD(LABELSIZE,INTSFP,UMR2EINT,XXPREAD,IXPREAD,NXX,NINTMX)
      NumIntRec = NXX
      IF (LABELSIZE.EQ.2) THEN
         do ASBufInd = 1, IABS(NumIntRec)
            NPACK = NPACKIC + ASBufInd
            ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(2*NPACK-1)
               P = ISHFT( LABEL, -16 )
               Q = IAND(  LABEL, 65535 )
               LABEL = IXPREAD(2*NPACK)
               R = ISHFT( LABEL, -16 )
               S = IAND(  LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, ca')
#else
               LABEL = IXPREAD(NPACK)
               ! P =       ISHFT( LABEL, -48 )
               P =       ISHFT( LABEL, IM48 )
              ! Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               Q = IAND( ISHFT( LABEL, IM32 ), 65535 )
              ! R = IAND( ISHFT( LABEL, -16 ), 65535 )
               R = IAND( ISHFT( LABEL, IM16 ), 65535 )
               S = IAND( LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, cb')
#endif
            ! Now all permutations are
            ! generated & added to buffer.
            call SetIntBuf(NumRec,RecLen,P,Q,R,S,+ASInt)
            call SetIntBuf(NumRec,RecLen,P,Q,S,R,-ASInt)
            call SetIntBuf(NumRec,RecLen,Q,P,R,S,-ASInt)
            call SetIntBuf(NumRec,RecLen,Q,P,S,R,+ASInt)
            call SetIntBuf(NumRec,RecLen,R,S,P,Q,+ASInt)
            call SetIntBuf(NumRec,RecLen,R,S,Q,P,-ASInt)
            call SetIntBuf(NumRec,RecLen,S,R,P,Q,-ASInt)
            call SetIntBuf(NumRec,RecLen,S,R,Q,P,+ASInt)
         enddo
      ELSE IF (LABELSIZE.EQ.1) THEN
         do ASBufInd = 1, IABS(NumIntRec)
            NPACK = NPACKIC + ASBufInd
            ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(NPACK)
               P =       ISHFT( LABEL, -24 )
               Q = IAND( ISHFT( LABEL, -16 ), 255 )
               R = IAND( ISHFT( LABEL,  -8 ), 255 )
               S = IAND(        LABEL,        255 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, da')
#else
               IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                  LABEL = IXPREAD(NPACK/2)
                  P = IAND( ISHFT( LABEL, -24 ), 255 )
                  Q = IAND( ISHFT( LABEL, -16 ), 255 )
                  R = IAND( ISHFT( LABEL,  -8 ), 255 )
                  S = IAND( LABEL, 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, db')
               ELSE
                  LABEL = IXPREAD((NPACK/2)+1)
                  !P = ISHFT( LABEL, -56 )
                  P = ISHFT( LABEL, IM56 )
                 ! Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  Q = IAND( ISHFT( LABEL, IM48 ), 255 )
                 ! R = IAND( ISHFT( LABEL, -40 ), 255 )
                  R = IAND( ISHFT( LABEL, IM40 ), 255 )
                 ! S = IAND( ISHFT( LABEL, -32 ), 255 )
                  S = IAND( ISHFT( LABEL, IM32 ), 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, dc')
               ENDIF
#endif
            ! Now all permutations are
            ! generated & added to buffer.
            call SetIntBuf(NumRec,RecLen,P,Q,R,S,+ASInt)
            call SetIntBuf(NumRec,RecLen,P,Q,S,R,-ASInt)
            call SetIntBuf(NumRec,RecLen,Q,P,R,S,-ASInt)
            call SetIntBuf(NumRec,RecLen,Q,P,S,R,+ASInt)
            call SetIntBuf(NumRec,RecLen,R,S,P,Q,+ASInt)
            call SetIntBuf(NumRec,RecLen,R,S,Q,P,-ASInt)
            call SetIntBuf(NumRec,RecLen,S,R,P,Q,-ASInt)
            call SetIntBuf(NumRec,RecLen,S,R,Q,P,+ASInt)
         enddo
      ENDIF
      if (NumIntRec <= 0) exit
   enddo

   ! One-body integral files are
   ! generated for all vaccums.

   do vac = 1, NSTATE
      do P = 1, NumRec
      do Q = 1, (RecLen/NVMO)
         PQ = RecordNumPQ(P,Q)
         IPQ = (P-1)*RecLen + (Q-1)*NVMO
         FInt = FINTGET(IPQ,vac)
         ! Diagonal part will change from iteration
         ! to iteration and hence will be added
         ! at the time of diagram calculation.
         ! if ( P == Q ) FInt = FInt + (OrbEn(P)-GenEn(P))
         IntF(PQ) = FInt
      enddo
      enddo
      ! Write the integrals
      write(UF,rec=vac) IntF
   enddo

   StartRec = 1
   EndngRec = NSTATE
   TotNZInt = NSTATE*SIZE(IntF)
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   return

end subroutine IntFileGen_OneBody

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE PREAD(LABSIZ,NWDVAR,IS,XX,IX,NX,NINTMX)
!--------------------------------------------------------------
!
!    Driver for reading Gamess integral files
!
!    On input: IS - lu of gamess integral file
!
!
! Called from : IntFileGen_OneBody
!               IntFileGen_TwoBody
!
!   testing routine
!
!--------------------------------------------------------------

   implicit none
   INTEGER LABSIZ, NWDVAR, IS, NX, NINTMX
   DOUBLE PRECISION XX(NINTMX)
   INTEGER IX(*)
   INTEGER MAX2, MAX2_DIM

   IF (LABSIZ.EQ.1) THEN
      IF(NWDVAR.EQ.2) MAX2 = NINTMX
      IF(NWDVAR.EQ.1) MAX2 = (NINTMX+1)/2
   ELSE IF (LABSIZ.EQ.2) THEN
      IF(NWDVAR.EQ.2) MAX2 = NINTMX*2
      IF(NWDVAR.EQ.1) MAX2 = NINTMX
   ELSE
      WRITE(UOUT,*) 'PREAD: CONFUSION WITH -LABSIZ- !'
      WRITE(UOUT,*) 'PREAD: LABSIZ =',LABSIZ
      call quit('PREAD: CONFUSION WITH -LABSIZ- !')
   END IF

#if defined (ARCH32BIT)
   MAX2_DIM = MAX2
#else
   if (LABSIZ == 2 ) then
!   ... for 64 bit AND more than 256 spinorobitals
    ! MAX2_DIM = MAX2/2 ! this does not work !!!
     MAX2_DIM = MAX2 ! ... this works at least for short time ...
   else
     MAX2_DIM = MAX2
   endif
#endif

   if (get_iprint() >= 20) then
     write(UOUT,'(/,2x,a,i2,a,i2)')   &
     '** Output from PREAD: LABSIZ=',LABSIZ,' NWDVAR=',NWDVAR
     write(UOUT,'(2x,2(a,i7))')    &
     'MAX2=',MAX2,' and dimension in PRD, MAX2_DIM=',MAX2_DIM
     !write(UOUT,*) 'size(IX)=',size(IX)
 !    write(UOUT,*) 'NAO, NSO, NFRZV:',NAO, NSO, NFRZV
   endif

   CALL PRD(IS,XX,IX,NX,NINTMX,MAX2,MAX2_DIM)

   RETURN

END SUBROUTINE PREAD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE PRD(IS,XX,IX,NX,NINTMX,MAX2,MAX2_DIM)
!--------------------------------------------------------------
!
!        Own reading from the open Gamess file 
!
!  Called from subroutine PREAD.
!
!--------------------------------------------------------------

   implicit none
   INTEGER IS, NX, NINTMX,MAX2,MAX2_DIM,ios
   DOUBLE PRECISION XX(NINTMX)
   integer (kind=8), save :: ic = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  INTEGER IX(MAX2)   ! ... for less than 256 spinorobitals
 !  INTEGER IX(MAX2/2) ! ... for 64 bit AND more than 256 spinorobitals

 ! INTEGER IX(*)     ! ... perhaps general case - this does not work 
 !                   - compiler Error: Data
 ! transfer element at (1) cannot be a full reference to an assumed-size array

 !  INTEGER IX(MAX2_DIM) ! ... this is general I hope /MI,16.2.2008...

! ... for 64bit/LABELSIZE=2 error:
!  PGFIO-F-219/unformatted read/unit=11/attempt to read/write past end of record.
   READ(IS,iostat=ios) NX,IX,XX
  ! READ(IS,ERR=1200) NX,IX,XX

  if (ios.ne.0) then
    write(UOUT,*) 'PRD: error in reading from GAMESS file !'
    call quit('PRD: error in reading from GAMESS file !')
  endif

   if (get_iprint() >= 22) then
     ic = ic + 1
     write(UOUT,'(2x,2(a,i8))')   &
     'PREAD: read NX=',NX,' ic=',ic
     call flushout
   endif

   RETURN
 1200  Call Quit('PRD: Gamess unformatted/sequential file reading error exit !')
END SUBROUTINE PRD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine  check_gms_integrals
!------------------------------------------------------------------------------------
!
!             opens, and read GAMESS prepared integrals
!
!------------------------------------------------------------------------------------
   integer NXX
   integer NPACKIC, NPACK, LABEL
   integer :: ASBufInd, NZVIntBufInd
   integer :: P, Q, R, S
   integer :: NumIntRec, EnPkInd
   integer  ::   IUNIQ, LABELSIZE, NINTMX
   double precision :: ASInt
   integer :: IM56=-56,IM48=-48,IM40=-40,IM32=-32,IM16=-16

   double precision, allocatable :: XXPREAD(:)
   integer, allocatable :: IXPREAD(:)


 call Header('output from check_gms_integrals',-1)

 write(UOUT,*) 'check_gms_integrals...'

 open(UMR2EINT,file=FileMR2EINT,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
 if ( get_iprint() >= 7) then
   write(UOUT,'(/2x,a,a/)')  &
   '... gamess integral file open and rewinded, FileMR2EINT=',FileMR2EINT
 endif

! Rewind MR2EINT file and
! Read first three records - see "modules/gmsmod.F90"

  rewind(UMR2EINT)

  read(UMR2EINT)
  read(UMR2EINT)
  read(UMR2EINT)  IUNIQ, LABELSIZE, NINTMX

  ! print-out important integrals parameters
  write(UOUT,'(/,2x,a,i3,a,i3,a,i3)') 'IUNIQ=',IUNIQ,' LABELSIZE=',LABELSIZE
  write(UOUT,'(2x,a,i8,/)') ' NINTMX=',NINTMX


  if (.not.allocated(XXPREAD)) allocate(XXPREAD(NINTMX))
  if (.not.allocated(IXPREAD)) allocate(IXPREAD(NINTMX))


! ... own cycle of reading ...
      NPACKIC = 0
      do
        ! ... INTSFP set in Parameters
         call PREAD(LABELSIZE,INTSFP,UMR2EINT,XXPREAD,IXPREAD,NXX,NINTMX)
         NumIntRec = NXX
     
         write(UOUT,'(2x,a,i8)') 'PREAD done, NumIntRec=',NumIntRec
           
         IF (LABELSIZE.EQ.2) THEN
            do ASBufInd = 1, IABS(NumIntRec)

               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(2*NPACK-1)
               P = ISHFT( LABEL, -16 )
               Q = IAND(  LABEL, 65535 )
               LABEL = IXPREAD(2*NPACK)
               R = ISHFT( LABEL, -16 )
               S = IAND(  LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, aa ')
               call print_gms_int_ind(P,Q,R,S,ASInt)
#else
               LABEL = IXPREAD(NPACK)
              ! P =       ISHFT( LABEL, -48 )
               P =       ISHFT( LABEL, IM48 )
              ! Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               Q = IAND( ISHFT( LABEL, IM32 ), 65535 )
              ! R = IAND( ISHFT( LABEL, -16 ), 65535 )
               R = IAND( ISHFT( LABEL, IM16 ), 65535 )
               S = IAND( LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, ab')
               call print_gms_int_ind(P,Q,R,S,ASInt)
#endif
            enddo
         ELSE IF (LABELSIZE.EQ.1) THEN
            do ASBufInd = 1, IABS(NumIntRec)
               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(NPACK)
               P =       ISHFT( LABEL, -24 )
               Q = IAND( ISHFT( LABEL, -16 ), 255 )
               R = IAND( ISHFT( LABEL,  -8 ), 255 )
               S = IAND(        LABEL,        255 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, ba')
               call print_gms_int_ind(P,Q,R,S,ASInt)
#else
               IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                  LABEL = IXPREAD(NPACK/2)
                  P = IAND( ISHFT( LABEL, -24 ), 255 )
                  Q = IAND( ISHFT( LABEL, -16 ), 255 )
                  R = IAND( ISHFT( LABEL,  -8 ), 255 )
                  S = IAND( LABEL, 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bb')
                  call print_gms_int_ind(P,Q,R,S,ASInt)
               ELSE
                  LABEL = IXPREAD((NPACK/2)+1)
                 ! P = ISHFT( LABEL, -56 )
                  P = ISHFT( LABEL, IM56 )
                !  Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  Q = IAND( ISHFT( LABEL, IM48 ), 255 )
                 ! R = IAND( ISHFT( LABEL, -40 ), 255 )
                  R = IAND( ISHFT( LABEL, IM40 ), 255 )
                 ! S = IAND( ISHFT( LABEL, -32 ), 255 )
                  S = IAND( ISHFT( LABEL, IM32 ), 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bc')  
                  call print_gms_int_ind(P,Q,R,S,ASInt)
               ENDIF
#endif
            enddo
         ENDIF
         if (get_iprint() >= 18) then
            write(UOUT,'(2x,a,i8)')   &
            ': read NumIntRec=',NumIntRec
            if (NumIntRec <= 0) write(UOUT,'(2x,a)') &
            '...negative, going to exit cycle of reading Gamess file !'
            call flushout
         endif
         if (NumIntRec <= 0) exit
      enddo

     ! ... deallocate stuff ...
     if (.not.allocated(XXPREAD)) deallocate(XXPREAD)
     if (.not.allocated(IXPREAD)) deallocate(IXPREAD)

end subroutine check_gms_integrals

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine get_gms_spacial_integrals(lu_file,IUNIQ,LABELSIZE,NINTMX,v,ind,nnint)
!------------------------------------------------------------------------------------
!
!  Read GAMESS prepared integrals
!
!  Called from HighSect FSCC codes 
!
!------------------------------------------------------------------------------------
!#include "implicit.h"
#include "priunit.h"
   integer, intent(in) :: lu_file
   integer, intent(in)  ::   IUNIQ, LABELSIZE, NINTMX
   integer, intent(inout) :: nnint

   integer, intent(inout) :: ind(*) ! indexes
   double precision, intent(inout) :: v(*) ! integrals

   integer :: NXX, i
   integer :: NPACKIC, NPACK, LABEL
   integer :: ASBufInd, NZVIntBufInd
   integer :: P, Q, R, S
   integer :: NumIntRec, EnPkInd
   double precision :: ASInt
   integer :: IM56=-56,IM48=-48,IM40=-40,IM32=-32,IM16=-16

   double precision, allocatable :: XXPREAD(:)
   integer, allocatable :: IXPREAD(:)

 call Header('output from get_gms_spacial_integrals',-1)

 UOUT=LUPRI ! needed for other codes !

  ! print-out important integrals parameters
 write(UOUT,'(/,2x,a,i3,a,i3,a,i3)') 'IUNIQ=',IUNIQ,' LABELSIZE=',LABELSIZE
 write(UOUT,'(2x,a,i8,/)') ' NINTMX=',NINTMX

  if (.not.allocated(XXPREAD)) allocate(XXPREAD(NINTMX))
  if (.not.allocated(IXPREAD)) allocate(IXPREAD(NINTMX))

  nnint = 1

! ... own cycle of reading ...
      NPACKIC = 0
      do
        ! ... INTSFP set in Parameters
        ! call PREAD(LABELSIZE,INTSFP,UMR2EINT,XXPREAD,IXPREAD,NXX,NINTMX)
         call PREAD(LABELSIZE,INTSFP,lu_file,XXPREAD,IXPREAD,NXX,NINTMX)

         NumIntRec = NXX
     
!       .. copy onto big arrays 
         call dcopy(IABS(NXX),XXPREAD,1,v(nnint),1)

         call icopy(IABS(NXX),IXPREAD,1,ind(nnint),1)
        ! do i=0, IABS(NXX)-1,2
        !    ind(nnint+i)   =IXPREAD(i+1)
        !    ind(nnint+i+1) =IXPREAD(i+2)
        ! enddo    

         nnint = nnint + IABS(NXX)

         write(UOUT,'(2x,a,i8)') 'PREAD done, NumIntRec=',NumIntRec
         write(UOUT,'(2x,a,i8)') 'nnint=',nnint
           
         IF (LABELSIZE.EQ.2) THEN
            do ASBufInd = 1, IABS(NumIntRec)

               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(2*NPACK-1)
               P = ISHFT( LABEL, -16 )
               Q = IAND(  LABEL, 65535 )
               LABEL = IXPREAD(2*NPACK)
               R = ISHFT( LABEL, -16 )
               S = IAND(  LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, aa ')
               call print_gms_int_ind_mul(LABEL,P,Q,R,S,ASInt)
#else
               LABEL = IXPREAD(NPACK)
              ! P =       ISHFT( LABEL, -48 )
               P =       ISHFT( LABEL, IM48 )
              ! Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               Q = IAND( ISHFT( LABEL, IM32 ), 65535 )
              ! R = IAND( ISHFT( LABEL, -16 ), 65535 )
               R = IAND( ISHFT( LABEL, IM16 ), 65535 )
               S = IAND( LABEL, 65535 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, ab')
               call print_gms_int_ind_mul(LABEL,P,Q,R,S,ASInt)
#endif
            enddo

         ELSE IF (LABELSIZE.EQ.1) THEN
            do ASBufInd = 1, IABS(NumIntRec)
               NPACK = NPACKIC + ASBufInd
               ASInt = XXPREAD(NPACK)

! select 32/64-bit code by preprocessor statement
#if defined (ARCH32BIT) 
               LABEL = IXPREAD(NPACK)
               P =       ISHFT( LABEL, -24 )
               Q = IAND( ISHFT( LABEL, -16 ), 255 )
               R = IAND( ISHFT( LABEL,  -8 ), 255 )
               S = IAND(        LABEL,        255 )
               call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'32 bit, ba')
               call print_gms_int_ind_mul(LABEL,P,Q,R,S,ASInt)
#else
               IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                  LABEL = IXPREAD(NPACK/2)
                  P = IAND( ISHFT( LABEL, -24 ), 255 )
                  Q = IAND( ISHFT( LABEL, -16 ), 255 )
                  R = IAND( ISHFT( LABEL,  -8 ), 255 )
                  S = IAND( LABEL, 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bb')
                  call print_gms_int_ind_mul(LABEL,P,Q,R,S,ASInt)
               ELSE
                  LABEL = IXPREAD((NPACK/2)+1)
                 ! P = ISHFT( LABEL, -56 )
                  P = ISHFT( LABEL, IM56 )
                 ! Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  Q = IAND( ISHFT( LABEL, IM48 ), 255 )
                 ! R = IAND( ISHFT( LABEL, -40 ), 255 )
                  R = IAND( ISHFT( LABEL, IM40 ), 255 )
                 ! S = IAND( ISHFT( LABEL, -32 ), 255 )
                  S = IAND( ISHFT( LABEL, IM32 ), 255 )
                  call check_PQRS(P,Q,R,S,LABELSIZE,LABEL,NPACK,'64 bit, bc')  
                  call print_gms_int_ind_mul(LABEL,P,Q,R,S,ASInt)
               ENDIF
#endif
            enddo
         ENDIF

         if (get_iprint() >= 18) then
            write(UOUT,'(2x,a,i8)')   &
            ': read NumIntRec=',NumIntRec
            if (NumIntRec <= 0) write(UOUT,'(2x,a)') &
            '...negative, going to exit cycle of reading Gamess file !'
            call flushout
         endif
         if (NumIntRec <= 0) exit
      enddo

     ! ... deallocate stuff ...
     if (.not.allocated(XXPREAD)) deallocate(XXPREAD)
     if (.not.allocated(IXPREAD)) deallocate(IXPREAD)

     nnint = nnint - 1
 
  return

end subroutine get_gms_spacial_integrals

subroutine print_gms_int_ind_mul(PckIndx,P,Q,R,S,Integral)
! .. copied:
!
!C        IN THIS PASS, CONSTRUCT ALL MO INTEGRALS (IJ|KL) WITH ... 
!C                                                                       
!C             JMIN <= J <= JMAX                                      
!C             I >= J                                                   
!C             1 <= K <= I                                              
!C             IF (I.NE.K) 1 <= L <= K                                  
!C             IF (I.EQ.K) 1 <= L <= J                                  
!C                                      
integer, intent(in) :: P,Q,R,S,PckIndx
real*8, intent(in) :: Integral
if (get_iprint() >= 20) then
  write(UOUT,'(i22,a,2i4,a,2i4,a,f18.12)') PckIndx,'-pck index: (',P,Q,'|',R,S,')=',Integral
endif
end subroutine print_gms_int_ind_mul


subroutine print_gms_int_ind(P,Q,R,S,ASInt)
integer, intent(in) :: P,Q,R,S
real*8, intent(in) :: ASInt

if (get_iprint() >= 20) then
  write(UOUT,'(2x,a,2i4,a,2i4,a,f18.12)') '<',P,Q,'||',R,S,'>=',ASInt
endif

end subroutine print_gms_int_ind

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine DrvIntGen_DIRAC
logical :: isVIND,isVINT,isFPP,isFHP,isFHH

write(UOUT,*) 'In DrvIntGen_DIRAC!'

! ... check the present of prepared DIRAC for HSCC integral files 
   inquire(file=FileVINT,exist=isVINT)
   inquire(file=FileVIND,exist=isVIND)
   inquire(file=FileFHH,exist=isFHH)
   inquire(file=FileFHP,exist=isFHP)
   inquire(file=FileFPP,exist=isFPP)

 if (isVINT) write(UOUT,*) 'there is',FileVINT

call quit('DrvIntGen_DIRAC')

end subroutine DrvIntGen_DIRAC

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VIND_VINT
!----------------------------------------------------------------------------------
!
! make control print of all two-electron integrals when the file is open and read  
!
! called from PrintRunInfo (source/runinfo.F90)
!
!----------------------------------------------------------------------------------

call header('control printout of all two-electron transformed integrals and their indexes',-1)

call control_print_V2H2P ! 1.
call control_print_VPPPP ! 2. block
call control_print_VHHHH ! 3.
call control_print_VHHPP ! 4.
call control_print_V2HPP ! 5.
call control_print_V2PHH ! 6.
call control_print_VPHHP ! 7.
call control_print_VPH2P ! 8. block
call control_print_VPHPP ! 9. block
call control_print_V2HHP ! 10.
call control_print_VHHHP ! 11.

end subroutine control_print_VIND_VINT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_V2H2P
! 
! ... read and prints 2H2P integrals & indexes
!
implicit none
integer :: I,J,IJ,St,En,IA,A,B,JB,IJB,icount
integer,          dimension(:), pointer :: VIND
double precision, dimension(:), pointer :: VINT
double precision :: VIJ_AB

call check_VIND_VINT_readiness

call header('****  01. V2H2P, 1324, <IJ||AB>  ****',-1)

nullify(VINT); nullify(VIND)
call ReadFile(V2H2P)
En = 0 ; icount=0
DO I = 1, NHMO
DO A = 1, NPMO
   IA = RecNumHP(I,A) ! this calculates the record number of HP index pair
   St = En + 1
   En = En + NZIntRec(IA)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IJB = 1, NZIntRec(IA)
      JB = VIND(IJB)
      CALL ResRecHP(JB,J,B) ! obtain J,B
      VIJ_AB = VINT(IJB)
      icount = icount + 1
      if (.not.allocated(MOLTRA_indx)) then
        WRITE(UOUT,829) icount,I,J,A,B,VIJ_AB
      else
        WRITE(UOUT,2829) icount,I,J,A,B,VIJ_AB,  & 
        MOLTRA_indx(I),MOLTRA_indx(A+NCMO),MOLTRA_indx(J),MOLTRA_indx(B+NCMO)
      endif
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for IA pair, I=',I,' A=',A, &
   ' NZIntRec(',IA,')=',NZIntRec(IA)
ENDDO
ENDDO

 829 format(4x,i8,1x,'V2H2P-1, 1324, <IJ||AB>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)
2829 format(4x,i8,1x,'V2H2P-1, 1324, <IJ||AB>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12,   & 
            ' MOLTRA:(',i3,1x,i3,'||',i3,1x,i3,')')

end subroutine control_print_V2H2P

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VPPPP
!
! ...   
!
 implicit none
 integer :: St,En,IAB,A,B,AB,C,D,CD,Block,iblock,ii
 integer :: ICD,iblock_CD,StBl,EnBl
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision ::  VAB_CD

write(UOUT,'(/,4x,a,i3,/)')  & 
 '****  02. VPPPP  <AB||CD>  ****, # of integr.blocks:',NumBlock(VPPPP)
write(UOUT,*) 'NumBlock(VPPPP)=',NumBlock(VPPPP)

ii=0
do  Block = 1, NumBlock(VPPPP)

 nullify(VINT); nullify(VIND)
 call ReadFile(VPPPP,Block,StBl,EnBl)
 write(UOUT,'(4x,a,2i4)') 'StBl,EnBl:',StBl,EnBl

 En = 0
 do AB = StBl, EnBl
   CALL ResRecMM(AB,A,B) ! obtain a,B
   St = En + 1
   En = En + NZIntRec(AB)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO ICD = 1, NZIntRec(AB)
      CD = VIND(ICD)
      CALL ResRecMM(CD,C,D) ! obtain C,D
      VAB_CD = VINT(ICD); ii=ii+1
      if (.not.allocated(MOLTRA_indx)) then
        WRITE(UOUT,830) ii,A,B,C,D,VAB_CD,Block
      else
        WRITE(UOUT,2830) ii,A,B,C,D,VAB_CD,Block, & 
        MOLTRA_indx(A+NCMO),MOLTRA_indx(C+NCMO),MOLTRA_indx(B+NCMO),MOLTRA_indx(D+NCMO)
      endif
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for AB pair, A=',A,' B=',B, &
   ' NZIntRec(',AB,')=',NZIntRec(AB)

 enddo
enddo
 830 format(4x,i8,'.',1x,'VPPPP-2,1234,<AB||CD>=<',i3,' ',i3,'||',i3,' ',i3,'>=', & 
            f18.12,' block:',i3)
 2830 format(4x,i8,'.',1x,'VPPPP-2,1234,<AB||CD>=<',i3,' ',i3,'||',i3,' ',i3,'>=', & 
            f18.12,' block:',i3,' MOLTRA:(',i3,1x,i3,'||',i3,1x,i3,')' )

end subroutine control_print_VPPPP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VHHHH
implicit none
integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
integer,          dimension(:), pointer :: VIND
double precision, dimension(:), pointer :: VINT
double precision :: VIJ_KL

write(UOUT,'(/,4x,a,/)') '****  03. VHHHH, 1234, <IJ||KL>  ****'
nullify(VINT); nullify(VIND)
call ReadFile(VHHHH)
En = 0;ii=0
DO I = 2, NHMO
DO J = 1, I-1
   IJ = RecNumMM(I,J) ! IJ index pair
   St = En + 1
   En = En + NZIntRec(IJ)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IKL = 1, NZIntRec(IJ)
      KL = VIND(IKL)
      CALL ResRecMM(KL,K,L) ! obtain K,L
      VIJ_KL = VINT(IKL); ii=ii+1
      WRITE(UOUT,831) ii,I,J,K,L,VIJ_KL
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for first index pair, I=',I,' J=',J, &
   ' NZIntRec(',IJ,')=',NZIntRec(IJ)
ENDDO
ENDDO
 831 format(4x,i5,'. VHHHH-3, 1234, <IJ||KL>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_VHHHH

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VHHPP
! ... print out VHHPP integrals ...
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VIJ_AB

write(UOUT,'(/,4x,a,/)') '****  04. VHHPP, 1234, <IJ||AB> ****'
nullify(VINT)
nullify(VIND)
call ReadFile(VHHPP)
En = 0; ii=0
DO I = 2, NHMO
DO J = 1, I-1
   IJ = RecNumMM(I,J) ! 
   St = En + 1
   En = En + NZIntRec(IJ)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IAB = 1, NZIntRec(IJ)
      AB = VIND(IAB)
      CALL ResRecMM(AB,A,B) ! obtain A,B
      VIJ_AB = VINT(IAB); ii=ii+1
      WRITE(UOUT,8311) ii,I,J,A,B,VIJ_AB
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for first index pair, I=',I,' J=',J, &
   ' NZIntRec(',IJ,')=',NZIntRec(IJ)
ENDDO
ENDDO
 8311 format(4x,i5,'. VHHPP-4, 1234, <IJ||AB>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_VHHPP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_V2HPP
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VIJ_AB

write(UOUT,'(/,4x,a,/)') '**** 05. V2HPP, 1234, <IJ||AB> ****'
nullify(VINT)
nullify(VIND)
call ReadFile(V2HPP)
En = 0; ii=0
DO I = 1, NHMO
DO J = 1, NHMO
   IJ = RecNum2H(I,J) ! get common index, IJ
   St = En + 1
   En = En + NZIntRec(IJ)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IAB = 1, NZIntRec(IJ)
      AB = VIND(IAB)
      CALL ResRecMM(AB,A,B) ! obtain A,B
      VIJ_AB = VINT(IAB); ii=ii+1
      WRITE(UOUT,832) ii,I,J,A,B,VIJ_AB
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                            &
   'Number of nonzero integrals for first index pair, I=',I,' J=',J, &
   ' NZIntRec(',IJ,')=',NZIntRec(IJ)
ENDDO
ENDDO
 832 format(4x,i5,'. V2HPP-5, 1234, <IJ||AB>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_V2HPP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_V2PHH
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision ::  VAB_IJ

write(UOUT,'(/,4x,a,/)') '**** 06. V2PHH, 1234, <AB||IJ> ****'
nullify(VINT)
nullify(VIND)
call ReadFile(V2PHH)
En = 0; ii=0
DO A = 1, NPMO
DO B = 1, NPMO
   AB = RecNum2P(A,B) ! get A,B
   St = En + 1
   En = En + NZIntRec(AB)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IIJ = 1, NZIntRec(AB)
      IJ = VIND(IIJ)
      CALL ResRecMM(IJ,I,J) ! obtain I,J
      VAB_IJ = VINT(IIJ); ii=ii+1
      WRITE(UOUT,833) ii,A,B,I,J,VAB_IJ
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')  & 
   'Number of nonzero integrals for first index pair, A=',A,' B=',B, &
   ' NZIntRec(',AB,')=',NZIntRec(AB)
ENDDO
ENDDO
 833 format(4x,i5,'. V2PHH-6, 1234, <AB||IJ>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)
end subroutine control_print_V2PHH

subroutine control_print_VPHHP
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VAJ_IB

write(UOUT,'(/,4x,a,/)') '****  07. VPHHP, 3124 <AJ||IB>  ****'
nullify(VINT)
nullify(VIND)
call ReadFile(VPHHP)
En = 0; ii=0
DO I = 1, NHMO
DO A = 1, NPMO
   IA = RecNumHP(I,A) !
   St = En + 1
   En = En + NZIntRec(IA)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IJB = 1, NZIntRec(IA)
      JB = VIND(IJB)
      CALL ResRecHP(JB,J,B) ! obtain J,B
      VAJ_IB = VINT(IJB); ii=ii+1
      WRITE(UOUT,834) ii,A,J,I,B,VAJ_IB
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')  & 
   'Number of nonzero integrals for first index pair, I=',I,' A=',A, &
   ' NZIntRec(',IA,')=',NZIntRec(IA)
ENDDO
ENDDO
 834 format(4x,i5,'. VPHHP-7, 1324, <AJ||IB>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)
end subroutine control_print_VPHHP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VPH2P

 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision ::  VAI_BC

write(UOUT,'(/,4x,a,/)') '****  08. VPH2P, 1324, <AI||BC> ****'
nullify(VINT)
nullify(VIND)
En = 0; ii=0
do  Block = 1, NumBlock(VPH2P)
call ReadFile(VPH2P,block,StBl,EnBl)
write(UOUT,*) 'StBl,EnBl:',StBl,EnBl

En=0
DO A = StBl, EnBl
DO B = 1, NPMO
   AB =  (A-1)*NPMO + B
   St = En + 1
   En = En + NZIntRec(AB)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IIC = 1, NZIntRec(AB)
      IC = VIND(IIC)
      CALL ResRecHP(IC,I,C) ! obtain I,C
      VAI_BC = VINT(IIC); ii=ii+1
      WRITE(UOUT,835) ii,A,I,B,C,VAI_BC
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                            &
   'Number of nonzero integrals for first index pair, A=',A,' B=',B, &
   ' NZIntRec(',AB,')=',NZIntRec(AB)
ENDDO
ENDDO
enddo
 835 format(4x,i5,'.VPH2P-8, 1324, <AI||BC>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)
end subroutine control_print_VPH2P

subroutine control_print_VPHPP
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,BC,IBC,Block,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VAI_BC

write(UOUT,'(/,4x,a,/)') '****  09. VPHPP, 1234, <AI||BC> ****'

ii=0
do  Block = 1, NumBlock(VPHPP)
call ReadFile(VPHPP,block,StBl,EnBl)
write(UOUT,*) 'StBl,EnBl:',StBl,EnBl
nullify(VINT); nullify(VIND)
En = 0
DO A = StBl,EnBl
DO I = 1, NHMO
   AI = (A-1)*NHMO + I
   St = En + 1
   En = En + NZIntRec(AI)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IBC = 1, NZIntRec(AI)
      BC = VIND(IBC)
      CALL ResRecMM(BC,B,C) ! obtain B,C
      VAI_BC = VINT(IBC); ii=ii+1
      WRITE(UOUT,836) ii,A,I,B,C,VAI_BC
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                            &
   'Number of nonzero integrals for first index pair, A=',A,' I=',I, &
   ' NZIntRec(',AI,')=',NZIntRec(AI)
ENDDO
ENDDO
enddo
 836 format(4x,i8,'.',1x,'VPHPP-9, 1234, <AI||BC>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_VPHPP


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_V2HHP
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VIJ_KA

write(UOUT,'(/,4x,a,/)') '**** 10. V2HHP, 1324, <IJ||KA> ****'
nullify(VINT)
nullify(VIND)
call ReadFile(V2HHP)
En = 0; ii=0
DO I = 1, NHMO
DO K = 1, NHMO
   IK = RecNum2H(I,K) ! get I,K
   St = En + 1
   En = En + NZIntRec(IK)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IJA = 1, NZIntRec(IK)
      JA = VIND(IJA)
      CALL ResRecHP(JA,J,A) ! obtain J,A
      VIJ_KA = VINT(IJA); ii=ii+1
      WRITE(UOUT,837) ii,I,J,K,A,VIJ_KA
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for IK pair, I=',I,' K=',K, &
   ' NZIntRec(',IK,')=',NZIntRec(IK)
ENDDO
ENDDO
 837 format(4x,i5,'. V2HHP-10, 1324, <IJ||KA>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_V2HHP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine control_print_VHHHP
! ... print HHHP integrals from the file
 implicit none
 integer :: I,J,IJ,St,En,IAB,A,B,AB,C,D,CD,IB,IKA,KA,Block,iblock,ii
 integer :: IJB,IA,JB,ICD,K,L,IKL,KL,IIJ,IIC,AI,IC,StBl,EnBl,IK,IJA,JA
 integer,          dimension(:), pointer :: VIND
 double precision, dimension(:), pointer :: VINT
 double precision :: VIJ_KA

write(UOUT,'(/,4x,a,/)') '****  11. VHHHP, 1234, <IJ||KA>  ****'
nullify(VINT); nullify(VIND)
call ReadFile(VHHHP)
En = 0; ii=0
DO I = 2, NHMO
DO J = 1, I-1
   IJ = RecNumMM(I,J) ! get I,J
   St = En + 1
   En = En + NZIntRec(IJ)
   IF (En < St) cycle
   VINT => IntBuf(St:En)
   VIND => IndBuf(St:En)
   DO IKA = 1, NZIntRec(IJ)
      KA = VIND(IKA)
      CALL ResRecHP(KA,K,A) ! get K,A
      VIJ_KA = VINT(IKA); ii=ii+1
      WRITE(UOUT,838) ii,I,J,K,A,VIJ_KA
   ENDDO
   WRITE(UOUT,'(1X,A,I3,A,I3,A,I5,A,I6)')                   & 
   'Number of nonzero integrals for IA pair, I=',I,' J=',J, &
   ' NZIntRec(',IJ,')=',NZIntRec(IJ)
ENDDO
ENDDO
 838 format(4x,i5,'. VHHHP-11, 1234, <IJ||KA>=<',i3,' ',i3,'||',i3,' ',i3,'>=',f18.12)

end subroutine control_print_VHHHP

end module IntGenMod
