!-----------------------------------------
! Module with variables previously defined
! at dcbgrd.h include file (COMMON /CBIGRD/, /CBLGRD/, /CBRGRD/, /CBDGRT/)
!-----------------------------------------
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
!     DMXWLT: Max wall time for slaves
!     DMXCPT: Max cpu time for slaves
   module dcbgrd_mod
   use mxcent_mod, only: MXCOOR
   implicit none
   real*8 :: GRADKN(mxcoor,0:2), GRADNU(mxcoor,0:2), GRADRO(mxcoor,0:2), &
             GRADER(mxcoor,0:3), GRADLV(mxcoor), GRADEE(mxcoor), GRADXC(mxcoor), &
             GRADNN(mxcoor), SCRGRD, DMXWLT, DMXCPT
   integer :: IPRGRD, INTCLASS, LUMGRD, IGRD_INTFLG, IGRD_INTBUF
   logical :: RELGRD, DOTRCK, DODERM, GRD_NOSMLV, GRD_NOSMLS, DONGRD
   end module dcbgrd_mod
