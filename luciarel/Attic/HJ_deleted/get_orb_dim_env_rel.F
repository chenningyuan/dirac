      SUBROUTINE GET_ORB_DIM_ENV_REL(ECORE_ENV)
*
* Obtain number of orbitals and basis functions from the
* programming environment.
* results stored in NAOS_ENV, NMOS_ENV
*
* Obtain environments CORE energy, ECORE_ENV
*
* Jeppe Olsen, December 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "crun.inc"
#include "lucinp.inc"
*
      IF(ENVIRO(1:6).EQ.'MOLCAS') THEN
CC        CALL GETOBS2(ECORE_ENV,NAOS_ENV,NMOS_ENV)
C            Getobs2(ECORE,NAOS_ENV,NMOS_ENV)
      ELSE IF(ENVIRO(1:6).EQ.'DALTON' ) THEN
CC        CALL GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV)
      ELSE IF(ENVIRO(1:4).EQ.'NONE') THEN
*. No environment,
        WRITE(6,*) 'GET_ORB_DIM_ENV  in problems '
        WRITE(6,*) 'No ENVIRO parameter defined '
        IZERO = 0
        CALL ISETVC(NAOS_ENV,IZERO,NSMOB)
        CALL ISETVC(NMOS_ENV,IZERO,NSMOB)
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' From GET_ORB_FROM_ENV : '
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NAOS_ENV'
        CALL IWRTMA(NAOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' NMOS_ENV'
        CALL IWRTMA(NMOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' ECORE_ENV=', ECORE_ENV
      END IF
*
      RETURN
      END
