C
C FILE    : dcbgascip.h
C
C     DEFINITION PARAMETERS GASCIP :
C
C     MAXASH_GASCIP is the limit for the GASCIP CI program
C           (64 bit integers = 64 orbitals)
C     N5 is the max no of reps produced by GMULTSF (and GMULTA)
C
      INTEGER MAXASH_GASCIP, N5
      PARAMETER (MAXASH_GASCIP = 64, N5 = 64)
C
C     /DCI_GASCIP/
C        NREP  : number of symmetries (max 32)
C        IRP   : symmetry number (either boson or fermion)
C        IRRP  : symmetry number
C                (1..NREP is boson, NREP+1..2*NREP is fermion)
C        MK2DEL, ... : MK2 limits
C        MULTB : symmetry multiplication table
C        IRPAMO: symmetry of spinors
C        IEL_STR(i)  : off-set in ISTR for strings with i electrons
C        NEL_STR(I)  : number of strings in ISTR with i electrons
C        NEL_STR_T = sum(nel_str(i), i=0..maxash_gascip)
C        IMK2(:,i), NMK2(:,i): i= 1:A-str, 2:B-str, 3:det
C
C     /DCC_GASCIP/
C        REPNA : string with irrep name for symmetry 1:2*NREP
C
      INTEGER I1_DCI_GASCIP,
     &     IRRP, IRP, NREP, MULTB,
     &     MK2DEL, MK2REF, MINMK2, MAXMK2,
     &     IRPAMO, NEL_STR_T, IEL_STR, NEL_STR,
     &     IMK2, NMK2,
     &     IPTA2O, IPTB2O, IPTO2A, IPTO2B,
     &     I2_DCI_GASCIP
      COMMON /DCI_GASCIP/ I1_DCI_GASCIP,
     &     IRRP, IRP, NREP, MULTB(N5,0:N5),
     &     MK2DEL, MK2REF, MINMK2, MAXMK2,
     &     IRPAMO(2*MAXASH_GASCIP), NEL_STR_T,
     &     IEL_STR(0:MAXASH_GASCIP), NEL_STR(0:MAXASH_GASCIP),
     &     IMK2(-MAXASH_GASCIP:MAXASH_GASCIP,3),
     &     NMK2(-MAXASH_GASCIP:MAXASH_GASCIP,3),
     &     IPTA2O(MAXASH_GASCIP),   IPTB2O(MAXASH_GASCIP),
     &     IPTO2A(2*MAXASH_GASCIP), IPTO2B(2*MAXASH_GASCIP),
     &     I2_DCI_GASCIP
!
!     IMPORTANT !!!
!     I1_DCI_GASCIP and I2_DCI_GASCIP are tags used to calculate the size of
!     this COMMON block for MPI transfers; they must remain at start/end
!

      CHARACTER*4 REPNA
      COMMON /DCC_GASCIP/ REPNA(N5)
C
C ********************************************************************
C
C     The remaining common block variables are set in the
C     *GASCIP input module and only used when DOGASCIP is true.
C
C     /DC_GASCIP_ORB/
C
      INTEGER IGASCIP_NISH,IGASCIP_NASH,
     &     IGASCIP_NASHT,IGASCIP_REP,IGASCIP_STATE,
     &     IGASCIP_NGAS,IGASCIP_NAELEC,
     &     IGASCIP_MK2REF,IGASCIP_MK2DEL,IGASCIP_MINMK2, IGASCIP_MAXMK2,
     &     IGASCIP_NGSH,IGASCIP_NGASSP
      LOGICAL GASCIP_UCIBOS
      COMMON /DC_GASCIP_ORB/ IGASCIP_NISH(2),IGASCIP_NASH(2),
     &     IGASCIP_NASHT,IGASCIP_REP,IGASCIP_STATE,
     &     IGASCIP_NGAS,IGASCIP_NAELEC,
     &     IGASCIP_MK2REF,IGASCIP_MK2DEL,IGASCIP_MINMK2, IGASCIP_MAXMK2,
     &     IGASCIP_NGSH(2,MAXASH_GASCIP),IGASCIP_NGASSP(2,MAXASH_GASCIP)
     &    ,GASCIP_UCIBOS
C --- end of dcbgascip.h ---
