      SUBROUTINE DISKUN
      IMPLICIT REAL*8           (A-H,O-Z)
*
* Assign logical unit numbers for LUCIA:
*
*==================================================
*  Complete revision motivated by MCSCF development
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*
*  60 - 69 :  Important files retained (DIAG,C,HC etc.)
*  ---------
*       60 :    CI diagonal
*       61 :    CI reference vector
*       62 :    CI sigma vector
*
*  70 - 79 :  Scratch files
*  ---------
*
*  80 - 99 :  Files in use in different context 
*             (e.g. standalone, INS, MOLCAS etc.)
*  ---------
*
*  Timo Fleig    Mar 19, 2002
*==================================================
*
*  Old setup:
*
* All files with some kind of input information  :  10 - 19
* Internal files (retained through job)          :  20 - 29
* Scratch files                                  :  30 - 40
* Files for machine information			 :  50 - 59
* Additional scratch files		       	 :  41,61,62
* Files for relativistic business (SO ints etc.) :  70 - 79
* All files containing final results             :  90 - 99
*
#include "clunit.inc"
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
      LU5 = LUIN
*. Output file
      LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 80
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals
      LU2INT = 81
      LU40 = LU2INT
*. Input , one electron integrals
      LU1SO = 82
      LU1INT = 83
      LU2 = LU1INT        ! for MOLCAS environment
      LU9 = LU1INT        ! for MOLECULE environment
*. Input , property one-electron integral files
      LUPRP  = 84
*. Sirius interface file
      LUSIR1 = 85
*. Input MOAO coefficients
      LU19 = 86
*. Input Spin-orbit integrals (X)
      LU69 = 87
*. Input Spin-orbit integrals (X)
      LU70 = 88
*. Input Spin-orbit integrals (X)
      LU71 = 89
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
*. CI vector
      LUC =   61
*. Sigma vector file
      LUHC =  62
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC2 = 71
      LUSC3 = 72
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC35 = 74
      LUSC36 = 75
      LUSC37 = 76
      LUSC38 = 77
      LUSC39 = 78
      LUSC41 = 79
*. Reference vectors for (transition) densities
      LUSC61 = 65
C     LUSC62 = 66
C  avoid unit 66. Used in DALTON
      LUSC62 = 95    ! conflict might occur with LU18 in genh1.F
*. Some more general purpose scratch
      LUSC92 = 92
      LUSC93 = 93
      LUSC94 = 94
      LUSC95 = 95
      LUSC96 = 96
      LUSC97 = 97
      LUSC98 = 98
      LUSC99 = 99
* =================
* Output files
* =================
*. Natural spinors in terms of double group adapted atomic spin-orbitals
*.
      LUMOUTR = 90
      LUMOUTI = 91
*
      RETURN
      END
