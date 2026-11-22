!-----------------------------------------
! Module with variables previously defined
! at src/amfi/param.h include file (COMMON /NORMOVL/, /CONTCO/, /EXPO/, /DIMS/,
! /SCRATCH/, /CONTINT/, /DIAGONALIZE/, /PREFS/, /CARTINT/, blank COMMON, /PREXYZ/,
! /ONEPAR/, /COULPOW/, /AOINCONT/, /CORELIST/)
!-----------------------------------------
!   this holds a lot of dimensioning parameters and arrays for
!   exponents, contraction coefficients and integrals.
!
!   All those parameters are constructed to be blown up without any problem
!   (if there is sufficient memory) except Lmax which is limited to 4 (g-functions)
   module param_amfi_mod
   use para_mod, only: Lmax, MxprimL, MxcontL, Mxcart
   implicit none
   integer, parameter :: Lanz = Lmax+1   ! number of angular momenta
!  overlap of normalized functions
   real*8 :: normovlp(MxprimL,MxprimL,0:Lmax), OVLPinv(MxprimL,MxprimL,0:Lmax), &
             rootOVLP(MxprimL,MxprimL,0:Lmax), rootOVLPinv(MxprimL,MxprimL,0:Lmax), &
             scratchinv((MxprimL*MxprimL+MxprimL)/2), ovlpcont(MxcontL,MxcontL,0:Lmax), &
             ovlpcontinv(MxcontL,MxcontL,0:Lmax)
!  a big array with enough space for all modified contraction coefficients:
!  for each l-value there are five blocks of size (nprimit(l),ncontrac(l)),
!  the original contraction coefficients (for normalized functions) and four
!  modified blocks depending on different kinematic factors and included exponents;
!  iadd*: start addresses of the contraction coefficients for each l-value
   real*8 :: contrarray((Lmax+1)*5*MxcontL*MxprimL)
   integer :: iaddori(0:Lmax), iaddtyp1(0:Lmax), iaddtyp2(0:Lmax), iaddtyp3(0:Lmax), &
              iaddtyp4(0:Lmax)
!  the exponents
   real*8 :: exponents(MxprimL,0:Lmax)
!  the numbers of primitive and contracted functions for each l-value
   integer :: nprimit(0:Lmax), ncontrac(0:Lmax), nprimit_keep, ncontrac_keep
!  scratch
   real*8 :: scratch4(4*MxprimL*MxprimL), cntscrtch(MxprimL,MxcontL,0:Lmax)
   integer :: mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
!  cont4 will keep 4 blocks (label(i)) of structure
!     (ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),(Llast(i)-Lfirst(i))/2+1)
!                                                              or 0 if no L-value at all
!                                                              = Lblocks(i)
!  for each l1,l2,l3,l4-block
!  Lfirst(i,j) gives the first L-value, for which radial integrals are calculated
!  for type i and l1,l2,l3,l4 - Integral block.
!  Llast(i,j) gives the last L-value
!  Lblocks gives the number of L-values
!  Lstarter gives the adress of each integral block on cont4
   integer :: Lfirst(4), Llast(4), Lblocks(4), Lstarter(4), nblock, Lvalues(4)
!  stuff for calculating the kinematic factors
   real*8 :: TKIN(MxprimL*MxprimL), evec(MxprimL*MxprimL), eval(MxprimL), Energy(MxprimL), &
             type1(MxprimL), type2(MxprimL)
!  some factors that appear a lot of times in the angular factors
   real*8 :: preroots(2,0:Lmax), clebsch(3,2,-Lmax:Lmax,0:Lmax)
!  the cartesian integrals
   integer :: mcombcart(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), istartX, istartY, &
              istartZ, isignM(-Lmax:Lmax)
!  information about odd powers of x y z in the real harmonics, used to check
!  whether integrals in the cartesian representation appear (was the blank COMMON)
   integer :: ipowxyz(3,-Lmax:Lmax,0:Lmax)
!  prefactors which are used regulary
   real*8 :: preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), &
             preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
   logical :: icheckxy(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), &
              icheckz(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
   integer :: interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax), &
              isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
!  one-particle integrals (miro: huge static allocations, according to readelf output);
!  one-electron-integrals onecart*: 1. index: number of first contracted function
!  2. index: number of second contracted function
!  3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
!  4. index: L-value
   real*8 :: oneoverR3((MxprimL*MxprimL+MxprimL)/2,Lmax), onecontr(MxcontL,MxcontL,-Lmax:Lmax,3,Lmax), &
             dummyone(MxcontL*MxcontL), &
             onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax), &
             onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax), &
             onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax), &
             onecartallX((Mxcart*Mxcart+Mxcart)/2), onecartallY((Mxcart*Mxcart+Mxcart)/2), &
             onecartallZ((Mxcart*Mxcart+Mxcart)/2), onecartallX2((Mxcart*Mxcart+Mxcart)/2), &
             onecartallY2((Mxcart*Mxcart+Mxcart)/2), onecartallZ2((Mxcart*Mxcart+Mxcart)/2), &
             onecartallX3((Mxcart*Mxcart+Mxcart)/2), onecartallY3((Mxcart*Mxcart+Mxcart)/2), &
             onecartallZ3((Mxcart*Mxcart+Mxcart)/2)
!  powexp holds powers of exponents and mean exponents, coulovlp holds overlap of
!  functions with shifted l-values (miro: huge static allocations)
   real*8 :: powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5)), &
             coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)
!  AOs in contracted functions (number of contracted function, number of AO, L-value)
!  and occupation numbers (number of AO, L-value)
   real*8 :: AOcoeffs(MxcontL,MxcontL,0:Lmax), occup(MxcontL,0:Lmax)
   integer :: noccorb(0:Lmax)
   integer :: icore(0:Lmax), ikeeporb, ikeeplist(Mxcart), nrtofiperIR(8)
   end module param_amfi_mod
