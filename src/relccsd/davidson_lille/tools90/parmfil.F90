module parmfil
!
!*******************************************************************************
!    ----- this module contains the global variables used -----
!    ----- in file handling                               -----
!===============================================================================
!     DO NOT modify this file
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!
!*******************************************************************************
!
!
   integer, parameter :: MAXFIL=30,MXFIL0=12,TOTFIL=14
!
   character (len=72), save, dimension (MAXFIL) :: dsn
   integer, save :: nfil0,nfile,nsup
   integer, save, dimension (MAXFIL) :: iunit
!  integer, parameter :: GEO=1,OMS=2,SYMO=3,IJKL=4,TCON=5,PCON=6,  &
!                        PSEU=7,ERML=8,SLSO=9,MFSO=10,INTE=11,DEN=12,  &
!                        CASM=13,ECAD=14,VEC=15,GUES=16,  &
!                        DIPX=25,DIPY=26,DIPZ=27
   integer, parameter :: GEO=1,OMS=2,SYMO=3,IJKL=4,SLSO=5,MFSO=6,INTE=7,  &
                         DEN=8,CASM=9,ECAD=10,VEC=11,GUES=12,TCON=13,PCON=14
  integer,parameter :: DEN_REF=2,L_MO=2,LX_MO=3,LY_MO=4,LZ_MO=5,R_MO=5,  &
                       RX_MO=6,RY_MO=7,RZ_MO=8,L2_MO=9,SIG_MO=10,RR_MO=10,  &
                       XX_MO=11,XY_MO=12,YY_MO=13,XZ_MO=14,YZ_MO=15,  &
                       ZZ_MO=16,NA_MO=16,NAX_MO=17,NAY_MO=18,NAZ_MO=19,  &
                       L_AO=20,LX_AO=21,LY_AO=22,LZ_AO=23,R_AO=23,RX_AO=24,  &
                       RY_AO=25,RZ_AO=26,RR_AO=26,XX_AO=27,XY_AO=28,YY_AO=29,  &
                       XZ_AO=30,YZ_AO=31,ZZ_AO=32,L2_AO=33,SOAA_MO=34,  &
                       SOAB_MO=35,NA_AO=35,NAX_AO=36,NAY_AO=37,NAZ_AO=38
   integer, parameter :: SR_CHIR=39,SC_CHIR=40,SC_CHIC=41,S_NAM=42,  &
                         DR_CHIR=43,DC_CHIR=44,DC_CHIC=45,D_NAM=46,  &
                         HR_CHIR=47,HC_CHIR=48,HC_CHIC=49,H_NAM=50,  &
                         S_OPE=51,S_SOPE=52,D_OPE=53,D_SOPE=54,  &
                         H_OPE=55,H_SOPE=56,S_SHL=57,S_MAT2=58,S_MAT3=59, &
                         S_TRM=60,S_DUM=61,D_TRM=62,D_DUM=63,O_SALC=64,  &
                         S_SYM=65,D_SYM=66,SD_OPE=67,D_MAT2=68,D_MAT3=69,  &
                         D_MAT4=70,SD_DEG=71
   integer, save :: prchir,pcchir,pcchic,pnam,popr,psopr
!
end module parmfil
