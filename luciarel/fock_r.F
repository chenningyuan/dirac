!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE FIFAM_REL(FIFA)
*
*. Construct inactive + active fock matrix
*
*. On input FIFAM Should be the inactive Fock matrix, in symmetry packed form
*
* Jeppe Olsen
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
      DIMENSION FIFA(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "lucinp.inc"
*
*
      CALL FIFAMS_REL(FIFA,WORK(KRHO1),IBSO,NSMOB,
     &                NTOOBS,NACOB,NTOOB,IREOST)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE FIFAMS_REL(FIFA,RHO1,IOBSM,NSMOB,LOBSM,NACOB,
     &                      NORBT,ISTOB)
*
* Update inactive fock matrix with active contributions
*
*     FIFA(I,J) = FIFA(I,J) + sum(k,l) ((ij!ab)-0.5*(ib!ja))
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION FIFA(*),RHO1(NACOB,NACOB)
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
*
*.  Assume spatial symmetric fock matrix
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM_REL(2,6,ISM,JSM,IJSM)
C            SYMCOM_REL(ITASK,IOBJ,I1,I2,I12)
        IF(JSM.NE.0) THEN
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I                          
C?            write(6,*) ' I J ', I,J
              IP = ISTOB(I)
              JP = ISTOB(J)
C?            write(6,*) ' IP JP ', IP,JP
               IJ= IJ + 1
               DO IA = 1, NACOB
                 DO IB = 1, NACOB
                   FIFA(IJ) = FIFA(IJ) 
     &           + RHO1(IA,IB)
     &           *(GTIJKL_REL(IP,JP,IA,IB)-0.5*GTIJKL_REL(IP,IB,IA,JP))
                 END DO
               END DO
            END DO
          END DO
        END IF
      END DO
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' FI + FA in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 1
       CALL APRBLM2(FIFA,LOBSM,LOBSM,NSMOB,ISYM)
      END IF
* 
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Construct Fock matrices from one- and two-electron integral
*  files/arrays and one-particle densities
*
      subroutine fockmat(FMAT,OINTS,ITRAF)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      dimension OINTS(ITRAF)

      return
      end
