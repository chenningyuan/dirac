module IntGenMod

   use Parameters
   use IOModule
   use GMSMod
   use VacMod
   implicit none
   private
   integer, parameter :: ASINTBUFSIZE  = 5460
   double precision, public :: CUTOFF
   integer         , dimension(NZVINTBUFSIZE) :: NZVPkIBuf
   double precision, dimension(NZVINTBUFSIZE) :: NZVIntBuf
   double precision, allocatable :: GenIntBuf(:)
   interface IntFileGen ! common name for 2 procedures
      module procedure IntFileGen_OneBody
      module procedure IntFileGen_TwoBody
   end interface IntFileGen
   public :: DrvIntGen

contains

subroutine DrvIntGen

   use TimingMod
   implicit none
   integer :: StVRec
   integer :: GenIntBufSize
   integer :: i

   call QEnter('DrvIntGen')

   ! ... Progress Report on Integral Generation ...

   write(UOUT,9011)
   write(UOUT,9012) CUTOFF
   write(UOUT,9013,advance='no')

9011 format('******************************',/,&
            'HSMRCC INTEGRAL FILE GENERATOR',/,&
            '******************************'/)    ! MI: gfortran does not accept '/' at the end...

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

   ! ... Open Gamess FileMR2EINT file ...

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
   StVRec = 1
   call IntFileGen(TimeTaken(:,V2H2P),RecLen(V2H2P),NumRec(V2H2P),StVRec,          &
                   StartRec(V2H2P),EndngRec(V2H2P),TotNZInt(V2H2P),NZIntRec_V2H2P, &
                   SIB_V2H2P,ResRecHP,RecNumHP,RecNumHP,                           &
                   SetQRangeHP,SetRRangeHP,SetSRangeHP,PackRSHP)
   call IntFileGen(TimeTaken(:,VPPPP),RecLen(VPPPP),NumRec(VPPPP),StVRec,          &
                   StartRec(VPPPP),EndngRec(VPPPP),TotNZInt(VPPPP),NZIntRec_VPPPP, &
                   SIB_VPPPP,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangePP,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,VHHHH),RecLen(VHHHH),NumRec(VHHHH),StVRec,          &
                   StartRec(VHHHH),EndngRec(VHHHH),TotNZInt(VHHHH),NZIntRec_VHHHH, &
                   SIB_VHHHH,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangeHH,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,VHHPP),RecLen(VHHPP),NumRec(VHHPP),StVRec,          &
                   StartRec(VHHPP),EndngRec(VHHPP),TotNZInt(VHHPP),NZIntRec_VHHPP, &
                   SIB_VHHPP,ResRecMM,RecNumMM,RecNumMM,                           &
                   SetQRangeMM,SetRRangePP,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,V2HPP),RecLen(V2HPP),NumRec(V2HPP),StVRec,          &
                   StartRec(V2HPP),EndngRec(V2HPP),TotNZInt(V2HPP),NZIntRec_V2HPP, &
                   SIB_V2HPP,ResRec2H,RecNum2H,RecNumMM,                           &
                   SetQRange2H,SetRRangePP,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,V2PHH),RecLen(V2PHH),NumRec(V2PHH),StVRec,          &
                   StartRec(V2PHH),EndngRec(V2PHH),TotNZInt(V2PHH),NZIntRec_V2PHH, &
                   SIB_V2PHH,ResRec2P,RecNum2P,RecNumMM,                           &
                   SetQRange2P,SetRRangeHH,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,VPHHP),RecLen(VPHHP),NumRec(VPHHP),StVRec,          &
                   StartRec(VPHHP),EndngRec(VPHHP),TotNZInt(VPHHP),NZIntRec_VPHHP, &
                   SIB_VPHHP,ResRecHP,RecNumHP,RecNumHP,                           &
                   SetQRangeHP,SetRRangeHP,SetSRangeHP,PackRSHP)
   call IntFileGen(TimeTaken(:,VPH2P),RecLen(VPH2P),NumRec(VPH2P),StVRec,          &
                   StartRec(VPH2P),EndngRec(VPH2P),TotNZInt(VPH2P),NZIntRec_VPH2P, &
                   SIB_VPH2P,ResRec2P,RecNum2P,RecNumHP,                           &
                   SetQRange2P,SetRRangeHP,SetSRangeHP,PackRSHP)
   call IntFileGen(TimeTaken(:,VPHPP),RecLen(VPHPP),NumRec(VPHPP),StVRec,          &
                   StartRec(VPHPP),EndngRec(VPHPP),TotNZInt(VPHPP),NZIntRec_VPHPP, &
                   SIB_VPHPP,ResRecPH,RecNumPH,RecNumMM,                           &
                   SetQRangePH,SetRRangePP,SetSRangeMM,PackRSMM)
   call IntFileGen(TimeTaken(:,V2HHP),RecLen(V2HHP),NumRec(V2HHP),StVRec,          &
                   StartRec(V2HHP),EndngRec(V2HHP),TotNZInt(V2HHP),NZIntRec_V2HHP, &
                   SIB_V2HHP,ResRec2H,RecNum2H,RecNumHP,                           &
                   SetQRange2H,SetRRangeHP,SetSRangeHP,PackRSHP)
   call IntFileGen(TimeTaken(:,VHHHP),RecLen(VHHHP),NumRec(VHHHP),StVRec,          &
                   StartRec(VHHHP),EndngRec(VHHHP),TotNZInt(VHHHP),NZIntRec_VHHHP, &
                   SIB_VHHHP,ResRecMM,RecNumMM,RecNumHP,                           &
                   SetQRangeMM,SetRRangeHP,SetSRangeHP,PackRSHP)

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

   ! ... Integral Generation Done ...

   write(UOUT,9021) ToTTime()
   write(UOUT,9022)

9021 format('Done in ',F6.1,' seconds.'/)
9022 format('**************************************',/,&
            'END OF HSMRCC INTEGRAL FILE GENERATION',/,&
            '**************************************'/)

   deallocate(GenIntBuf) ! ... deallocate big buffer memory ...

   close(UMR2EINT)

   call QExit('DrvIntGen')

   return

end subroutine DrvIntGen

include "maxmin.h"

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

   RecNum = (I-1)*NPMO + A
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + B
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

   Call QEnter('SIB_V2HPP')

   ! For contributing to V2HPP, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <I,J||A,B>

   !if ( P  > NHMO ) return
   if ( P  > NHMO ) Goto 124

   !if ( Q  > NHMO ) return
   if ( Q  > NHMO ) Goto 124

   !if ( R <= NCMO ) return
   if ( R <= NCMO ) Goto 124

   !if ( S <= NCMO ) return
   if ( S <= NCMO ) Goto 124

   I = P
   J = Q
   A = R - NCMO
   B = S - NCMO

   RecNum = (I-1)*NHMO + J
   !if ( RecNum < StRec ) return
   if ( RecNum < StRec ) Goto 124

   !if ( RecNum > EnRec ) return
   if ( RecNum > EnRec ) Goto 124

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( A < B ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(A,B)
   GenIntBuf(BufInd) = IntVal

 124 Continue
   Call QExit('SIB_V2HPP')

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

   Call QEnter('SIB_V2PHH')

   ! For contributing to V2PHH, first index
   ! pair should be particle-particle type and second
   ! index pair should be hole-hole type.
   !
   ! <P,Q||R,S> = <A,B||I,J>

   !if ( P <= NCMO ) return
   if ( P <= NCMO ) Goto 125

   !if ( Q <= NCMO ) return
   if ( Q <= NCMO ) Goto 125

   !if ( R  > NHMO ) return
   if ( R  > NHMO ) Goto 125

   !if ( S  > NHMO ) return
   if ( S  > NHMO ) Goto 125

   A = P - NCMO
   B = Q - NCMO
   I = R
   J = S

   RecNum = (A-1)*NPMO + B
   !if ( RecNum < StRec ) return
   if ( RecNum < StRec ) Goto 125

   !if ( RecNum > EnRec ) return
   if ( RecNum > EnRec ) Goto 125

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(I,J)
   GenIntBuf(BufInd) = IntVal

 125 Continue
   Call QExit('SIB_V2PHH')
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

   Call QEnter('SIB_VHHHH')

   ! For contributing to VHHHH, all
   ! indices should be hole type.
   !
   ! <P,Q||R,S> = <I,J||K,L>

   !if ( P > NHMO ) return
   if ( P > NHMO ) Goto 126

   !if ( Q > NHMO ) return
   if ( Q > NHMO ) Goto 126

   !if ( R > NHMO ) return
   if ( R > NHMO ) Goto 126

   !if ( S > NHMO ) return
   if ( S > NHMO ) Goto 126

   I = P
   J = Q
   K = R
   L = S

   RecNum = MAXMIN(I,J)
   !if ( RecNum < StRec ) return
   if ( RecNum < StRec ) Goto 126

   !if ( RecNum > EnRec ) return
   if ( RecNum > EnRec ) Goto 126

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal
   IF ( K < L ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(K,L)
   GenIntBuf(BufInd) = IntVal

 126 Continue
   Call QExit('SIB_VHHHH')
   
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

   Call QEnter('SIB_VHHPP')

   ! For contributing to VHHPP, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <I,J||A,B>

   !if ( P  > NHMO ) return
   if ( P  > NHMO ) Goto 127

   !if ( Q  > NHMO ) return
   if ( Q  > NHMO ) Goto 127

   !if ( R <= NCMO ) return
   if ( R <= NCMO ) Goto 127

   !if ( S <= NCMO ) return
   if ( S <= NCMO ) Goto 127

   I = P
   J = Q
   A = R - NCMO
   B = S - NCMO

   RecNum = MaxMin(I,J)
   !if ( RecNum < StRec ) return
   if ( RecNum < StRec ) Goto 127

   !if ( RecNum > EnRec ) return
   if ( RecNum > EnRec ) Goto 127

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal
   IF ( A < B ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(A,B)
   GenIntBuf(BufInd) = IntVal

 127 Continue
   Call QExit('SIB_VHHPP')
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

   Call QEnter('SIB_VPPPP')

   ! For contributing to VPPPP, all
   ! indices should be particle type.
   !
   ! <P,Q||R,S> = <A,B||C,D>

   !if ( P <= NCMO ) return
   if ( P <= NCMO ) Goto 128

   !if ( Q <= NCMO ) return
   if ( Q <= NCMO ) Goto 128

   !if ( R <= NCMO ) return
   if ( R <= NCMO ) Goto 128

   !if ( S <= NCMO ) return
   if ( S <= NCMO ) Goto 128

   A = P - NCMO
   B = Q - NCMO
   C = R - NCMO
   D = S - NCMO

   RecNum = MaxMin(A,B)
   !if ( RecNum < StRec ) return
   if ( RecNum < StRec ) Goto 128

   !if ( RecNum > EnRec ) return
   if ( RecNum > EnRec ) Goto 128

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( A < B ) IntVal = -IntVal
   IF ( C < D ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(C,D)
   GenIntBuf(BufInd) = IntVal

  128 Continue
   Call QExit('SIB_VPPPP')

   return

end subroutine SIB_VPPPP

subroutine SIB_VPHHP(StRec,EnRec,RecLen,P,Q,R,S,AsyInt)

   implicit none
   integer, intent(in) :: StRec, EnRec, RecLen
   integer, intent(in) :: P, Q, R, S
   double precision, intent(in) :: AsyInt
   integer :: A, J, I, B
   integer :: RecNum, BufInd

   ! For contributing to V2H2P, first index
   ! pair should be hole-hole type and second
   ! index pair should be particle-paricle type.
   !
   ! <P,Q||R,S> = <A,J||I,B>

   if ( P <= NCMO ) return
   if ( Q  > NHMO ) return
   if ( R  > NHMO ) return
   if ( S <= NCMO ) return

   A = P - NCMO
   J = Q
   I = R
   B = S - NCMO

   RecNum = (I-1)*NPMO + A
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + B
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

   RecNum = (I-1)*NHMO + K
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (J-1)*NPMO + A
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
   A = S - NCMO

   RecNum = MAXMIN(I,J)
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( I < J ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + (K-1)*NPMO + A
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
   ! <P,Q||R,S> = <A,I||B,C>

   if ( P <= NCMO ) return
   if ( Q  > NHMO ) return
   if ( R <= NCMO ) return
   if ( S <= NCMO ) return

   A = P - NCMO
   I = Q
   B = R - NCMO
   C = S - NCMO

   RecNum = (A-1)*NPMO + B
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   BufInd = (RecNum - StRec)*RecLen + (I-1)*NPMO + C
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

   RecNum = (A-1)*NHMO + I
   if ( RecNum < StRec ) return
   if ( RecNum > EnRec ) return

   ! Now the given integral is
   ! in the current block.

   IntVal = AsyInt
   IF ( B < C ) IntVal = -IntVal

   BufInd = (RecNum - StRec)*RecLen + MAXMIN(B,C)
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
   ! finding P & Q from PQ.

   P = (PQ-1)/NPMO + 1
   Q = PQ - (P-1)*NPMO

   return

end subroutine ResRecHP

subroutine ResRecPH(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NHMO + Q
   ! finding P & Q from PQ.

   P = (PQ-1)/NHMO + 1
   Q = PQ - (P-1)*NHMO

   return

end subroutine ResRecPH

subroutine ResRec2H(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NHMO + Q
   ! finding P & Q from PQ.

   P = (PQ-1)/NHMO + 1
   Q = PQ - (P-1)*NHMO

   return

end subroutine ResRec2H

subroutine ResRec2P(PQ,P,Q)

   implicit none
   integer, intent(in)  :: PQ
   integer, intent(out) :: P, Q

   ! PQ = (P-1)*NPMO + Q
   ! finding P & Q from PQ.

   P = (PQ-1)/NPMO + 1
   Q = PQ - (P-1)*NPMO

   return

end subroutine ResRec2P

function RecNumMM(P,Q)

   implicit none
   integer, intent(in) :: P, Q
   integer :: RecNumMM

   if (P==Q) stop 'Error in RecNumMM'
   RecNumMM = MAXMIN(P,Q)

   return

end function RecNumMM

function RecNumHP(P,Q)

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

   call QEnter('SetQRangeMM')

   QL = 1
   QH = (P-1)

   !if (QH<QL) stop 'Error in SetQRangeMM'
   if (QH<QL) then 
     write(UOUT,*) 'Error:QH<QL!; QH,QL=',QH,QL
     Call Quit('QH<QL...error in SetQRangeMM!')
   endif

   if (P==SP) QL = SQ
   if (P==EP) QH = EQ

   call QExit('SetQRangeMM')

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
   if (SH<SL) stop 'Error in SetSRangeMM'

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

   if (R==S) stop 'Error in PackRSMM'
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
call QEnter('check_PQRS')

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

call QExit('check_PQRS')
end subroutine check_PQRS

subroutine IntFileGen_TwoBody(TimeTaken,RecLen,NumRec,StVRec,              &
                              StartRec,EndngRec,TotNZInt,NZIntRec,         &
                              SetIntBuf,ResolvRec,RecordNumPQ,RecordNumRS, &
                              SetQRange,SetRRange,SetSRange,PackRS)

   implicit none
   ! ... Dummy Argument Variables ...
   real,    intent(inout) :: TimeTaken(2)
   integer, intent(in)    :: RecLen, NumRec
   integer, intent(out)   :: StartRec, EndngRec, TotNZInt
   integer, intent(inout) :: StVRec
   integer, intent(out)   :: NZIntRec(NumRec)
   ! ... Dummy Procedure and Function Interfaces
   interface
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
   integer :: P, Q, R, S
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
   integer NXX
   integer NPACKIC, NPACK, LABEL
   real :: dtim, dtime, tarray(2)
   external :: dtime

   call QEnter('IntFileGen_TwoBody')

   dtim = dtime(tarray)
   ! Calculate number of records that fit
   ! into the buffer.
   BufSize = SIZE(GenIntBuf)
   if ( BufSize < RecLen ) then
      stop 'Buffer must be big enough to contain a record'
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
               P =       ISHFT( LABEL, -48 )
               Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               R = IAND( ISHFT( LABEL, -16 ), 65535 )
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
                  P = ISHFT( LABEL, -56 )
                  Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  R = IAND( ISHFT( LABEL, -40 ), 255 )
                  S = IAND( ISHFT( LABEL, -32 ), 255 )
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

      ! Now the buffer is completely
      ! constructed. Collect non-zero
      ! integrals and their packing
      ! indices into buffer and write
      ! to the disk.

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
                  if ( DABS(GenIntBuf(PQRS)) > CUTOFF ) then
                     NZVIntBufInd = NZVIntBufInd + 1
                     NZIntRec(PQ) = NZIntRec(PQ) + 1
                     if ( NZVIntBufInd > NZVINTBUFSIZE ) then
                        write(UVINT,rec=NxtRec) NZVIntBuf
                        write(UVIND,rec=NxtRec) NZVPkIBuf
                        NxtRec = NxtRec + 1
                        NZVIntBufInd = 1
                        NZVPkIBuf = 0
                        NZVIntBuf = DBLE(0)
                     endif
                     NZVIntBuf(NZVIntBufInd) = GenIntBuf(PQRS)
                     NZVPkIBuf(NZVIntBufInd) = PackRS(R,S)
                  endif
               enddo
            enddo
            TotNZInt = TotNZInt + NZIntRec(PQ)
         enddo
      enddo
   enddo

   ! End of all blocks

   EndngRec = NxtRec
   if ( TotNZInt /= 0 ) then
      write(UVINT,rec=NxtRec) NZVIntBuf
      write(UVIND,rec=NxtRec) NZVPkIBuf
      StVRec = EndngRec + 1
   endif
   dtim = dtime(tarray)
   TimeTaken = TimeTaken + tarray

   call QExit('IntFileGen_TwoBody')
   return

end subroutine IntFileGen_TwoBody

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
   integer NXX
   integer NPACKIC, NPACK, LABEL
   real :: dtim, dtime, tarray(2)
   external :: dtime
   ! ... Automatic Array ...
   double precision :: IntF(NumRec*(RecLen/NVMO))

   call QEnter('IntFileGen_OneBody')

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
               P =       ISHFT( LABEL, -48 )
               Q = IAND( ISHFT( LABEL, -32 ), 65535 )
               R = IAND( ISHFT( LABEL, -16 ), 65535 )
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
                  P = ISHFT( LABEL, -56 )
                  Q = IAND( ISHFT( LABEL, -48 ), 255 )
                  R = IAND( ISHFT( LABEL, -40 ), 255 )
                  S = IAND( ISHFT( LABEL, -32 ), 255 )
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

   call QExit('IntFileGen_OneBody')
   return

end subroutine IntFileGen_OneBody

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE PREAD(LABSIZ,NWDVAR,IS,XX,IX,NX,NINTMX)
!--------------------------------------------------------------
!
!    Driver for reading Gamess integral files
!
! Called from : IntFileGen_OneBody
!               IntFileGen_TwoBody
!
!--------------------------------------------------------------

   implicit none
   INTEGER LABSIZ, NWDVAR, IS, NX, NINTMX
   DOUBLE PRECISION XX(NINTMX)
   INTEGER IX(*)
   INTEGER MAX2, MAX2_DIM

   Call QEnter('PREAD')

   IF (LABSIZ.EQ.1) THEN
      IF(NWDVAR.EQ.2) MAX2 = NINTMX
      IF(NWDVAR.EQ.1) MAX2 = (NINTMX+1)/2
   ELSE IF (LABSIZ.EQ.2) THEN
      IF(NWDVAR.EQ.2) MAX2 = NINTMX*2
      IF(NWDVAR.EQ.1) MAX2 = NINTMX
   ELSE
      WRITE(6,*) 'PREAD: CONFUSION WITH -LABSIZ-'
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

   Call QExit('PREAD')
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
   INTEGER IS, NX, NINTMX,MAX2,MAX2_DIM
   DOUBLE PRECISION XX(NINTMX)
   integer (kind=8), save :: ic = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!  INTEGER IX(MAX2)   ! ... for less than 256 spinorobitals
 !  INTEGER IX(MAX2/2) ! ... for 64 bit AND more than 256 spinorobitals

 ! INTEGER IX(*)     ! ... perhaps general case - this does not work 
 !                   - compiler Error: Data
 ! transfer element at (1) cannot be a full reference to an assumed-size array

   INTEGER IX(MAX2_DIM) ! ... this is general I hope /MI,16.2.2008...

   Call QEnter('PRD')

! ... for 64bit/LABELSIZE=2 error:
!  PGFIO-F-219/unformatted read/unit=11/attempt to read/write past end of record.
   READ(IS) NX,IX,XX
  ! READ(IS,ERR=1200) NX,IX,XX

   if (get_iprint() >= 22) then
     ic = ic + 1
     write(UOUT,'(2x,2(a,i8))')   &
     'PREAD: read NX=',NX,' ic=',ic
     call flushout
   endif

   Call QExit('PRD')
   RETURN

 1200  Call Quit('PRD: Gamess unformatted/sequential file reading error exit !')

END SUBROUTINE PRD

end module IntGenMod
