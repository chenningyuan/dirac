C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbgrd.h
C
C Depencies: mxcent.h
C
C*** RELGRD is true if we want to calculate a relativistic gradient
C     IGRD_INTFLG: INTFLG for gradient
C     IGRD_INTBUF: buffer to INTFLG
C     DOTRCK    : Estimate LS and SS gradient and skip them if "small"
C                 Default for gradient calculation is: .FALSE.
C                 Default for geometry-optimization is: .TRUE.
C     DODERM    : used internally in ONEDRV
C     INTCLASS  : used internally in ONEDRV
C     GRD_NOSMLV: No small-small nuclear attraction 
C                 (eg. when using Levy-Leblond Hamiltonian)
C     GRD_NOSMLS: No small-small reorthonormalization term
C                 (eg. when using Levy-Leblond Hamiltonian)
Cjth: 1997/09/23> Added GRADLV
Cjth: 1997/06/20: Added GRADKN(0:2),GRADNU(0:2),GRADRO(0:2),GRADER(0:3)
Cjth              GRADKN is kinetic energy
Cjth              GRADNU is nuclear attraction
Cjth              GRADRO is reorthonormalization
Cjth              GRADER is electron-electron repulsion
Cjth              GRADLV is gradient of, and I quote: "Model inter-
Cjth                atomic SS-integral contribution by classical
Cjth                repulsion of small component atomic charges"
C
C     DMXWLT: Max wall time for slaves
C     DMXCPT: Max cpu time for slaves
C
      LOGICAL RELGRD,DOTRCK,DODERM,GRD_NOSMLV,
     &        GRD_NOSMLS,DONGRD
      COMMON /CBIGRD/ IPRGRD,INTCLASS,LUMGRD,
     &     IGRD_INTFLG, IGRD_INTBUF
      COMMON /CBLGRD/ RELGRD,DOTRCK,DODERM,GRD_NOSMLV,
     &                GRD_NOSMLS,DONGRD
      COMMON /CBRGRD/ GRADKN(MXCOOR,0:2),
     &                GRADNU(MXCOOR,0:2),
     &                GRADRO(MXCOOR,0:2),
     &                GRADER(MXCOOR,0:3),
     &                GRADLV(MXCOOR),
     &                GRADEE(MXCOOR),
     &                GRADNN(MXCOOR), SCRGRD
      COMMON /CBDGRT/ DMXWLT,DMXCPT
