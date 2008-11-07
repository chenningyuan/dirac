C
C VERSION : $Revision$
C DATE    : $Date$
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
     &     IGOS_MAXMK2,IGOS_MINMK2,IPRGOS
      COMMON /DCCGOS/ REPNA(N5)
C
C     /GOSORB/
C     
      LOGICAL GOS_UCIBOS
      COMMON /GOSORB/ IGOS_NISH(2),IGOS_NASH(2),IGOS_NASHT,IGOSREP,
     &                IGOS_STATE,GOS_UCIBOS,
     &                IGOS_NGAS,IGOS_NAELEC,IGOS_MK2REF,IGOS_MK2DEL,
     &                IGOS_NGSH(2,MAXASH_GOSCI),
     &                IGOS_NGSP(2,MAXASH_GOSCI)

