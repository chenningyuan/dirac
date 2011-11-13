!
! VERSION : $Revision: 1.10 $
! DATE    : $Date: 2004/05/26 13:06:54 $
! FILE    : dcbgrd.h
!
! Depencies: mxcent.h
!
!*** RELGRD is true if we want to calculate a relativistic gradient
!     IGRD_INTFLG: INTFLG for gradient
!     IGRD_INTBUF: buffer to INTFLG
!     DOTRCK    : Estimate LS and SS gradient and skip them if "small"
!                 Default for gradient calculation is: .FALSE.
!                 Default for geometry-optimization is: .TRUE.
!     DODERM    : used internally in ONEDRV
!     INTCLASS  : used internally in ONEDRV
!     GRD_NOSMLV: No small-small nuclear attraction 
!                 (eg. when using Levy-Leblond Hamiltonian)
!     GRD_NOSMLS: No small-small reorthonormalization term
!                 (eg. when using Levy-Leblond Hamiltonian)
!jth: 1997/09/23> Added GRADLV
!jth: 1997/06/20: Added GRADKN(0:2),GRADNU(0:2),GRADRO(0:2),GRADER(0:3)
!jth              GRADKN is kinetic energy
!jth              GRADNU is nuclear attraction
!jth              GRADRO is reorthonormalization
!jth              GRADER is electron-electron repulsion
!jth              GRADLV is gradient of, and I quote: "Model inter-
!jth                atomic SS-integral contribution by classical
!jth                repulsion of small component atomic charges"
!
!     DMXWLT: Max wall time for slaves
!     DMXCPT: Max cpu time for slaves
!
      LOGICAL RELGRD, DOTRCK, DODERM, GRD_NOSMLV, GRD_NOSMLS, DONGRD
      INTEGER IPRGRD, INTCLASS, LUMGRD, IGRD_INTFLG, IGRD_INTBUF
      COMMON /CBIGRD/ IPRGRD, INTCLASS, LUMGRD, IGRD_INTFLG,            &
     &                IGRD_INTBUF

      COMMON /CBLGRD/ RELGRD, DOTRCK, DODERM, GRD_NOSMLV, GRD_NOSMLS,   &
     &                DONGRD

      REAL*8 GRADKN(MXCOOR, 0:2), GRADNU(MXCOOR, 0:2),                  &
     &       GRADRO(MXCOOR, 0:2), GRADER(MXCOOR, 0:3), GRADLV(MXCOOR),  &
     &       GRADEE(MXCOOR), GRADNN(MXCOOR), SCRGRD
      COMMON /CBRGRD/ GRADKN, GRADNU, GRADRO, GRADER, GRADLV, GRADEE,   &
     &                GRADNN, SCRGRD

      REAL*8 DMXWLT, DMXCPT
      COMMON /CBDGRT/ DMXWLT, DMXCPT
