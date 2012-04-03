C
C FILE    : dcbgascip.h
C
C     DEFINITION PARAMETERS GASCIP :
C
C     MAXASH_GASCIP is the limit for the GASCIP CI program
C           (32 bit integers - 1 sign bit = 31 orbitals)
C     N5 is the max no of reps produced by GMULTSF (and GMULTA)
C
      PARAMETER (MAXASH_GASCIP = 63, N5 = 64)
C
C     /DCI_GASCIP/
C        NREP  : number of symmetries (max 32)
C        IRP   : symmetry number (either boson or fermion)
C        IRRP  : symmetry number 
C                (1..NREP is boson, NREP+1..2*NREP is fermion)
C        MULTB : symmetry multiplication table.
C        IRPAMO: symmetry of spinors
C     
C     /DCC_GASCIP/
C        REPNA : string identifying symmetry IRRP
C
      CHARACTER*4 REPNA
      COMMON /DCI_GASCIP/ IRRP, IRP, NREP, MULTB(N5,0:N5),
     &     IRPAMO(2*MAXASH_GASCIP),
     &     IMK2(-MAXASH_GASCIP:MAXASH_GASCIP),
     &     NMK2(-MAXASH_GASCIP:MAXASH_GASCIP),
     &     IPTA2O(MAXASH_GASCIP),   IPTB2O(MAXASH_GASCIP), 
     &     IPTO2A(2*MAXASH_GASCIP), IPTO2B(2*MAXASH_GASCIP)
      COMMON /DCC_GASCIP/ REPNA(N5)
C
C ********************************************************************
C
C     The remaining common block variables are set in the
C     *GASCIP input module and used when DOGASCIP is true.
C
C     /DC_GASCIP_ORB/
C     
      LOGICAL GASCIP_UCIBOS
      COMMON /DC_GASCIP_ORB/ IGASCIP_NISH(2),IGASCIP_NASH(2),
     &     IGASCIP_NASHT,IGASCIP_REP,IGASCIP_STATE,GASCIP_UCIBOS,
     &     IGASCIP_NGAS,IGASCIP_NAELEC,
     &     IGASCIP_MK2REF,IGASCIP_MK2DEL,IGASCIP_MINMK2, IGASCIP_MAXMK2,
     &     IGASCIP_NGSH(2,MAXASH_GASCIP),IGASCIP_NGSP(2,MAXASH_GASCIP)
C --- end of dcbgascip.h ---
