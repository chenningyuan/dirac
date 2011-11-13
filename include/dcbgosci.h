C
C VERSION : $Revision: 5471 $
C DATE    : $Date: 2007-02-06 15:33:25 +0100 (Tue, 06 Feb 2007) $
C FILE    : dcbgosci.h
C
C     DEFINITION PARAMETERS GOSCIP :
C
C     MAXASH_GOSCI is the limit for the GOSCI CI program
C           (32 bit integers - 1 sign bit = 31 orbitals)
C     N5 is the max no of reps produced by GMULTSF (and GMULTA)
C
      PARAMETER (MAXASH_GOSCI = 31, N5 = 64)
C
C     /DCIGOS/
C        NREP  : number of symmetries (max 32)
C        IRP   : symmetry number (either boson or fermion)
C        IRRP  : symmetry number 
C                (1..NREP is boson, NREP+1..2*NREP is fermion)
C        MULTB : symmetry multiplication table.
C        IRPAMO: symmetry of spinors
C     
C     /DCCGOS/
C        REPNA : string identifying symmetry IRRP
C
      CHARACTER*4 REPNA
      COMMON /DCIGOS/ IRRP, IRP, NREP, MULTB(N5,0:N5),
     &     IRPAMO(2*MAXASH_GOSCI),
     &     IMK2(-MAXASH_GOSCI:MAXASH_GOSCI),
     &     NMK2(-MAXASH_GOSCI:MAXASH_GOSCI),
     &     IPTA2O(MAXASH_GOSCI),   IPTB2O(MAXASH_GOSCI), 
     &     IPTO2A(2*MAXASH_GOSCI), IPTO2B(2*MAXASH_GOSCI),
     &     MAXMK2,MINMK2
      COMMON /DCCGOS/ REPNA(N5)
