module mconst
!
!*******************************************************************************
!  ----- This module contains the values of some widely used   -----
!  ----- values                                                -----
!===============================================================================
!  PI32 = sqrt(Pi**3)
!  PI12 = 1/sqrt(Pi)
!  PI212 = 2/sqrt(Pi)
!  A4SRPI = 4/sqrt(Pi)
!  PI54 = Pi**(5/4)
!  SQRPI = sqrt(Pi)
!  SQRFPI= = sqrt(4*Pi)
!  RLN10 = Log(10)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
   use kindvar
   implicit none
!
!     ----- parameters -----
!
   character (len=8), parameter :: BLANK8='        '
   character (len=6), parameter :: BLANK6='      '
   character (len=4), parameter :: IBLANK='    '
   double precision, parameter :: ZERO=0.0D+00,       &
                                  HALF=0.5D+00,       &
                                  ONE=1.0D+00,        &
                                  TWO=2.0D+00,        &
                                  THREE=3.0D+00,      &
                                  FOUR=4.0D+00,       &
                                  FIVE=5.0D+00,       &
                                  SIX=6.0D+00,        &
                                  SEVEN=7.0D00,       &
                                  EIGHT=8.0D+00,      &
                                  TEN=10.0D+00,       &
                                  TWELVE=12.0D+00,    &
                                  FORTIN=14.0D+00,    &
                                  FIFTIN=15.0D+00,    &
                                  SIXTIN=16.0D+00,    &
                                  EIGTIN=18.0D+00,    &
                                  FIFTY=50.0D+00,     &
                                  PT1=0.1D+00,        &
                                  PT187=1.875D+00,    &
                                  FOURTH=0.25D+00,    &
                                  PT25=0.25D+00,      &
                                  PT625=0.625D+00,    &
                                  PT75=0.75D+00,      &
                                  PT6562=6.5625D+00,  &
                                  TOLN=1.0D-10,       &
                                  H_TOL=1.0D-08,      &
                                  V15=1.5D+00,        &
                                  V25=2.5D+00,        &
                                  V375=3.75D+00,      &
                                  V125=1.25D+00,      &
                                  V225=2.25D+00,      &
                                  V1875=1.875D+00,    &
                                  HUNDRD=1.0D+02,     &
                                  HUND05=105.0D+00,   &
                                  TENM1=1.0D-01,      &
                                  TENM2=1.0D-02,      &
                                  TENM3=1.0D-03,      &
                                  TENM4=1.0D-04,      &
                                  TENM5=1.0D-05,      &
                                  TENM6=1.0D-06,      &
                                  TENM7=1.0D-07,      &
                                  TENM8=1.0D-08,      &
                                  TENM9=1.0D-09,      &
                                  TENM10=1.0D-10,     &
                                  TENM11=1.0D-11,     &
                                  TENM12=1.0D-12,     &
                                  TENM15=1.0D-15,     &
                                  TENM20=1.0D-20,     &
                                  TENM30=1.0D-30,     &
                                  THOUS=1.0D+03,      &
                                  TENP9=1.0D+09,      &
                                  DEGREE=360.0D+00,PIDEG=180.0D+00,  &
                                  PI=3.14159265358979D+00,       &
                                  A4SRPI=2.2567583341910251D+00, &
                                  PI32=5.56832799683170D+00,     &
                                  PI212=1.1283791670955D+00,     &
                                  PI12=0.564189583547756D+00,    &
                                  PI54=4.18251339837960D+00,     &
                                  SQRPI=1.772453850905D+00,      &
                                  SQRFPI=3.5449077018110D+00,    &
                                  SQRT2=1.41421356237310D+00,    &
                                  SQRT3=1.73205080756888D+00,    &
                                  SQRT5=2.23606797749979D+00,    &
                                  DSQRT6=0.244948974278318D+01,  &
                                  SQRT7=2.64575131106459D+00,    &
                                  DSQRT10=0.316227766016838D+01, &
                                  RLN10= 0.230258509299405D+01
   double precision, parameter :: SQRT53=SQRT5/SQRT3,TWOPI=TWO*PI,FOURPI=FOUR*PI,&
                                  M_PI=PI/PIDEG,FOURPISQ=FOURPI**2
#ifndef _SETCPL_
   double complex, parameter :: CZERO=dcmplx(ZERO,ZERO),CONE=dcmplx(ONE,ZERO),  &
                                C_I=dcmplx(ZERO,ONE)
#else
   complex (kind=K_C8), save :: C_I,CZERO,C_ONE
#endif
!
!     ----- machine dependent constants -----
!
   double precision, parameter :: MACHEP=epsilon(1.0D+00)
!  double precision, parameter :: MACHEP=2.220446049250313E-16
#ifdef _SETCPL_
   contains
   subroutine setcpl
!
!***********************************************************************
!    set double complex parameters (bloody pathscale)
!-----------------------------------------------------------------------
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- parameters -----
!
!     ----- arguments -----
!
!     ----- local variables -----
!
!     ----- variables in common -----
!
!     +++++
!     start
!     +++++
!
      C_I = dcmplx(ZERO,ONE)
      CZERO = dcmplx(ZERO,ZERO)
      CONE = dcmplx(ONE,ZERO)
      return
   end subroutine setcpl
#endif
end module mconst
