C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id$
C
#ifdef UNDEF
/* Comdeck log */
960818-vebjornb: New module taking care of redundant internal
                 coordinates. The procedures are called from abaopt.F
#endif
C  /* Deck inired */
      SUBROUTINE INIRED(WILBMT,BMTINV,PJINMT,TMPMAT,TMPMT2,TMPMT3,
     &                  TMPMT4)
C
C     Initializes everything that has to do with redundant internal
C     coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
      DIMENSION WILBMT(*), BMTINV(*), PJINMT(*)
      DIMENSION TMPMAT(2*MXCOOR,2*MXCOOR), TMPMT2(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMT3(2*MXCOOR,2*MXCOOR), TMPMT4(2*MXCOOR,2*MXCOOR)
C
C     First we have to determine the redundant internal coordinates.
C
      CALL FNDRED(TMPMAT,TMPMT2)
C
C     Then we construct Wilsons B matrix
C
      CALL GETWIL(TMPMAT,WILBMT)
C
C     ... and the derivative of the B matrix
C
      CALL GETDWL(TMPMAT,TMPMT2,TMPMT3,WILBMT)
C
C     We also need the inverse of the rectangular B matrix.
C
      CALL GTBINV(TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTINV,PJINMT)
C
      RETURN
      END

C  /* Deck fndred */
      SUBROUTINE FNDRED(ATMARR,IBOND)
C     
C     Finds natural redundant internal coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
      DIMENSION ATMARR(MXCENT,8), IBOND(MXCENT,MXCENT)
      CHARACTER CRDTXT*22, EMPTXT*14
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL DZERO(ATMARR,8*MXCENT)
      CALL IZERO(INTCRD,12*MXCOOR)
      CALL IZERO(IBOND,MXCENT*MXCENT)
      IINTCR = 0
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from FNDRED','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C      
      CALL ATMINI(ATMARR,IATOM,.FALSE.)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Expanded array of atoms (Aangstrom)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,5,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C
C     We find all bonds.
C
      CALL FNDBND(ATMARR,IBOND,IATOM)
C
C     Then all angles
C
      CALL FNDANG(ATMARR,IBOND,IATOM)
C
C     Finally all dihedral angles
C
      CALL FNDIHD(ATMARR,IBOND,IATOM)
C
C     Output
C
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Redundant internal coordinates',-1)
         EMPTXT = '              '
         DO 10 I = 1, IINTCR
            CRDTXT = '                      '
            IF (INTCRD(I,1) .LT. 10) THEN
               IF (INTCRD(I,1) .EQ. 1) CRDTXT='Regular bond'
               IF (INTCRD(I,1) .EQ. 3) CRDTXT='Interfragment bond'
               IF (INTCRD(I,1) .EQ. 5) CRDTXT='Hydrogen bond'
               WRITE(LUPRI,'(A,I3,A,A22,2I4)') EMPTXT,I,'   ',CRDTXT,
     &                                           INTCRD(I,2),INTCRD(I,3)
            ELSE IF (INTCRD(I,1) .GT. 20) THEN
               IF (INTCRD(I,1) .EQ. 21) CRDTXT='Dihedral angle'
               IF (INTCRD(I,1) .EQ. 23) CRDTXT='Extra dihedral angle'
               WRITE(LUPRI,'(A,I3,A,A22,4I4)') EMPTXT,I,'   ',CRDTXT,
     &                   INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
            ELSE
               IF (INTCRD(I,1) .EQ. 11) CRDTXT='Regular angle'
               IF (INTCRD(I,1) .EQ. 12) CRDTXT='Additional angle'
               WRITE(LUPRI,'(A,I3,A,A22,3I4)') EMPTXT,I,'   ',CRDTXT,
     &                               INTCRD(I,2),INTCRD(I,3),INTCRD(I,4)
            END IF
 10      CONTINUE
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,I5)')
     &      ' Total number of redundant internal coordinates: ',IINTCR
      END IF
      IF (IINTCR .GT. 2*MXCOOR) 
     &     CALL QUIT('Too many internal ccordinates!')
      RETURN
      END

C  /* Deck fndbnd */
      SUBROUTINE FNDBND(ATMARR,IBOND,IATOM)
C     
C     Finds all bonds.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), IBOND(MXCENT,MXCENT)
      LOGICAL POSSIB
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We have to find all regular bonds.
C
      IBEFOR = IINTCR
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Regular bonds',-1)
      END IF
      DO 10 I = 1, IATOM - 1
         DO 15 J = I + 1, IATOM
            RADI = ATMARR(I,5)
            RADJ = ATMARR(J,5)            
            VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
            VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
            VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
            DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
C
C     Regular bonds have dist. less than 1.3 times sum of cov. radii
C
            IF (DIST .LE. 1.3D0*(RADI+RADJ)) THEN
               IBOND(I,J) = 1
               IBOND(J,I) = 1
               IINTCR = IINTCR + 1
               INTCRD(IINTCR,1) = 1
               INTCRD(IINTCR,2) = I
               INTCRD(IINTCR,3) = J
               IF (IPRINT .GE. IPRMAX) THEN
                  WRITE(LUPRI,'(A,2I4)') ' Bond between atoms: ',i,j
               END IF
            END IF
 15         CONTINUE
 10   CONTINUE
      IF ((IPRINT .GE. IPRMAX) .AND. (IBEFOR .EQ. IINTCR)) THEN
         WRITE(LUPRI,'(A)') ' None were found.'
      END IF
C
C     Next we have to check if single atoms or fragments are too far
C     away to be bonded by normal bonds, in which case interfragment
C     bond(s) must be assigned.
C
C     The first step is to assign fragment numbers to all atoms, this
C     process is done iteratively.
C
      IBEFOR = IINTCR
      IFRAG = 1
      IAT = 1
 20   CONTINUE
      ATMARR(IAT,6) = IFRAG*1.0D0
 25   CONTINUE
      ICHANG = 0
      DO 30 II = 1,IATOM
         IF (NINT(ATMARR(II,6)) .EQ. IFRAG) THEN
            DO 35 JJ = 1, IATOM
               IF ((IBOND(II,JJ) .GE. 1) .AND.
     &              (NINT(ATMARR(JJ,6)) .EQ. 0)) THEN
                  ATMARR(JJ,6) = IFRAG*1.0D0
                  ICHANG = 1
               END IF
 35         CONTINUE
         END IF
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Assigning atoms to fragments',-1)
         CALL OUTPUT(ATMARR,1,IATOM,6,6,MXCENT,8,1,LUPRI)
      END IF
      IF (ICHANG .GT. 0) GOTO 25
C      
      II = 1
 32   CONTINUE
      IF ((II .LE. IATOM) .AND. (ATMARR(II,6) .GT. 0.5D0)) THEN
         II = II + 1
         GOTO 32
      END IF
      IF (II .LE. IATOM) THEN
         IFRAG = IFRAG + 1
         IAT = II
         ATMARR(IAT,6) = IFRAG*1.0D0
         GOTO 20
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,I4)') 
     &        ' Number of separate fragments found: ', IFRAG
      END IF
C     
C     Two fragments get connected through the shortest interfragment distance.
C     All distances shorter than the minimum of 2.0 Ang. and this distance
C     are also marked as interfragment bonds.
C
      IF (IFRAG .GT. 1) THEN
         DSTMIN = 1.0D10
         DO 40 II = 1, IATOM
            DO 45 JJ = 1, IATOM
               IF (NINT(ATMARR(II,6)) .NE. NINT(ATMARR(JJ,6))) THEN
                  VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                  VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                  VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                  DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                  IF (DIST .LT. DSTMIN) THEN
                     DSTMIN = DIST
                     IATM1 = II
                     IFRG1 = NINT(ATMARR(II,6))
                     IATM2 = JJ
                     IFRG2 = NINT(ATMARR(JJ,6))
                  END IF
               END IF
 45         CONTINUE
 40      CONTINUE
         IBOND(IATM1,IATM2) = 3
         IBOND(IATM2,IATM1) = 3
         IINTCR = IINTCR + 1
         INTCRD(IINTCR,1) = 3
         INTCRD(IINTCR,2) = IATM1
         INTCRD(IINTCR,3) = IATM2
         BNDLIM = MIN(1.3D0*DSTMIN, 2.0D0)
C
C     Additional interfragment bonds
C
         DO 50 II = 1, IATOM
            IF (NINT(ATMARR(II,6)) .EQ. IFRG1) THEN
               DO 55 JJ = 1, IATOM
                  IF ((IBOND(II,JJ) .EQ. 0) .AND.
     &                 (NINT(ATMARR(JJ,6)) .EQ. IFRG2)) THEN
                     VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                     VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                     VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     IF (DIST .LE. BNDLIM) THEN
                        IBOND(II,JJ) = 3
                        IBOND(JJ,II) = 3
                        IINTCR = IINTCR + 1
                        INTCRD(IINTCR,1) = 3
                        IF (II .LT. JJ) THEN
                           INTCRD(IINTCR,2) = II
                           INTCRD(IINTCR,3) = JJ
                        ELSE
                           INTCRD(IINTCR,2) = JJ
                           INTCRD(IINTCR,3) = II
                        END IF
                     END IF
                  END IF
 55            CONTINUE
            END IF
 50      CONTINUE
C
C     The process of connecting fragments runs iteratively, until all
C     atoms are connected, i.e. the number of fragments is one.
C
         IF (IFRAG .GT. 1) THEN
            DO 57 I = 1, IATOM
               ATMARR(I,6) = D0
 57         CONTINUE
            IAT = 1
            IFRAG = 1
            GOTO 20
         END IF
      END IF
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Interfragment bonds',-1)
         DO 60 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,2I4)')
     &           ' Interfragment bond between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3)
 60      CONTINUE
      END IF
C
C     Then we have to find all hydrogen bonds in the system.
C     X --- H ... Y (X,Y = N, O, F, P, S, Cl)
C
      IBEFOR = IINTCR
      DO 70 I = 1, IATOM
         IF (NINT(ATMARR(I,1)) .EQ. 1) THEN
            POSSIB = .FALSE.
            DO 72 J = 1, IATOM
               ID = NINT(ATMARR(J,1))
               IF ((IBOND(I,J) .EQ. 1) .AND. ((ID .EQ. 7) .OR.
     &              (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15)
     &              .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                  POSSIB = .TRUE.
                  VEC2(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC2(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC2(3) = ATMARR(J,4)-ATMARR(I,4)
               END IF
 72         CONTINUE
            IF (POSSIB) THEN
               VDWR1 = VDWRAD(1)
               DO 75 J = 1, IATOM
                  ID = NINT(ATMARR(J,1))
                  IF ((IBOND(I,J) .EQ. 0) .AND. ((ID .EQ. 7) .OR.
     &                 (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15)
     &                  .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                     VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                     VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                     VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     ICHG = NINT(ATMARR(J,1))
C
C     We require the distance to be less than 0.9 times the sum of the
C     van der Waals radii and larger than the covalent radii. Finally
C     we say that the angle (X-H-Y) has to be larger than 90 degrees.
C
                     VDWRI = VDWRAD(ICHG)
                     IF ((DIST .LE. 0.9D0*(VDWR1+VDWRI))
     &                .AND. (DIST .GE. (RADIUS(1)+RADIUS(ICHG)))) THEN
                        BNDANG = VECANG(VEC1,VEC2)
                        IF (ABS(BNDANG) .GE. 0.49D0*PI) THEN
                           IBOND(I,J) = 5
                           IBOND(J,I) = 5
                           IINTCR = IINTCR + 1
                           INTCRD(IINTCR,1) = 5
                           IF (I .LT. J) THEN
                              INTCRD(IINTCR,2) = I
                              INTCRD(IINTCR,3) = J
                           ELSE
                              INTCRD(IINTCR,2) = J
                              INTCRD(IINTCR,3) = I
                           END IF
                        END IF
                     END IF
                  END IF
 75            CONTINUE
            END IF
         END IF
 70   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Hydrogen bonds',-1)
         DO 77 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,2I4)')
     &           ' Hydrogen bond between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3)
 77      CONTINUE
      END IF
      RETURN
      END

C  /* Deck fndang */
      SUBROUTINE FNDANG(ATMARR,IBOND,IATOM)
C     
C     Finds all angles.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), IBOND(MXCENT,MXCENT)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We have to find all angles between bonded atoms
C
      IBEFOR = IINTCR
      DO 10 I = 1, IATOM
         NUMBND = 0
         DO 15 J = 1, IATOM
            IF (IBOND(I,J) .GT. 0) NUMBND = NUMBND + 1
 15      CONTINUE
         IF (NUMBND .GT. 1) THEN
            DO 20 J = 1, IATOM - 1
               IF (IBOND(I,J) .GT. 0) THEN
                  VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                  DO 25 K = J + 1, IATOM
                     IF (IBOND(I,K) .GT. 0) THEN
                        IINTCR = IINTCR + 1
                        INTCRD(IINTCR,1) = 11
                        INTCRD(IINTCR,2) = J
                        INTCRD(IINTCR,3) = I
                        INTCRD(IINTCR,4) = K
C
C     We have to check for linear or near linear bonds (these get two
C     internal orthogonal bond coordinates). The "critical" bond angle
C     is 175 degrees.
C
                        VEC2(1) = ATMARR(K,2)-ATMARR(I,2)
                        VEC2(2) = ATMARR(K,3)-ATMARR(I,3)
                        VEC2(3) = ATMARR(K,4)-ATMARR(I,4)
                        BNDANG = VECANG(VEC1,VEC2)
                        IF (ABS(BNDANG) .GE. 175D0*PI/180D0) THEN
                           IINTCR = IINTCR + 1
                           INTCRD(IINTCR,1) = 12
                           INTCRD(IINTCR,2) = J
                           INTCRD(IINTCR,3) = I
                           INTCRD(IINTCR,4) = K
                        END IF
                     END IF
 25               CONTINUE
               END IF
 20         CONTINUE
         END IF
 10   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Angles',-1)
         DO 30 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,3I4)')
     &           ' Angles between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4)
 30      CONTINUE
      END IF
      RETURN
      END

C  /* Deck fndihd */
      SUBROUTINE FNDIHD(ATMARR,IBOND,IATOM)
C     
C     Finds all dihedral angles.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMARR(MXCENT,8), IBOND(MXCENT,MXCENT)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175D0*PI/180D0)
C
C     We have to find all dihedral angles between bonded atoms
C
      IF (IATOM .LE. 3) RETURN
      IBEFOR = IINTCR
      DO 10 I = 1, IATOM
         DO 15 J = 1, IATOM
            IF (IBOND(I,J) .GT. 0) THEN
               DO 20 K = 1, IATOM
                  IF ((IBOND(J,K) .GT. 0) .AND. (K .NE. I)) THEN
                     DO 25 L = I, IATOM
                        IF ((IBOND(K,L) .GT. 0) .AND. (L .NE. J) .AND.
     &                       (L .NE. I)) THEN
C
C     Check for linearity.
C
                           VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
                           VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
                           VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
                           VEC2(1) = ATMARR(K,2)-ATMARR(J,2)
                           VEC2(2) = ATMARR(K,3)-ATMARR(J,3)
                           VEC2(3) = ATMARR(K,4)-ATMARR(J,4)
                           BNDAN1  = VECANG(VEC1,VEC2)
                           VEC1(1) = -VEC2(1)
                           VEC1(2) = -VEC2(2)
                           VEC1(3) = -VEC2(3)
                           VEC2(1) = ATMARR(L,2)-ATMARR(K,2)
                           VEC2(2) = ATMARR(L,3)-ATMARR(K,3)
                           VEC2(3) = ATMARR(L,4)-ATMARR(K,4)
                           BNDAN2  = VECANG(VEC1,VEC2)
                           IF ((ABS(BNDAN1) .LT. DEG175)
     &                          .AND. (ABS(BNDAN2) .LT. DEG175)) THEN
                              IINTCR = IINTCR + 1
                              INTCRD(IINTCR,1) = 21
                              INTCRD(IINTCR,2) = I
                              INTCRD(IINTCR,3) = J
                              INTCRD(IINTCR,4) = K
                              INTCRD(IINTCR,5) = L
                           END IF
                        END IF
 25                  CONTINUE
                  END IF
 20            CONTINUE
            END IF
 15      CONTINUE
 10   CONTINUE
C
      IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
         CALL HEADER('Dihedral angles',-1)
         DO 30 I = IBEFOR + 1, IINTCR
            WRITE(LUPRI,'(A,4I4)')
     &           ' Dihedral angles between atoms: ',
     &           INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
 30      CONTINUE
      END IF
C
C     Systems with more than three atoms where no regular dihedral
C     angles has been found, need addition of dihedral angle to
C     take care of out of plane bending.
C
      IF ((IATOM .GE. 4) .AND. (IBEFOR .EQ. IINTCR)) THEN
         ITRI = 0
         DO 50 I = 1, IATOM
            INUMBN = 0
            DO 55 J = 1, IATOM
               IF (IBOND(I,J) .GT. 0) INUMBN = INUMBN + 1
 55         CONTINUE
            IF (INUMBN .GE. 3) ITRI = I
 50      CONTINUE
      END IF
      IF (ITRI .GT. 0) THEN
         IBEFOR = IINTCR
         DO 60 I = IATOM, 1, -1
            IF (IBOND(I,ITRI) .GT. 0) II = I
 60      CONTINUE
         VEC1(1) = ATMARR(II,2)-ATMARR(ITRI,2)
         VEC1(2) = ATMARR(II,3)-ATMARR(ITRI,3)
         VEC1(3) = ATMARR(II,4)-ATMARR(ITRI,4)
         DO 65 J = II + 1, IATOM
            IF (IBOND(J,ITRI) .GT. 0) THEN
               VEC2(1) = ATMARR(J,2)-ATMARR(ITRI,2)
               VEC2(2) = ATMARR(J,3)-ATMARR(ITRI,3)
               VEC2(3) = ATMARR(J,4)-ATMARR(ITRI,4)
               BNDANG  = VECANG(VEC1,VEC2)
               IF (ABS(BNDANG) .LT. DEG175) JJ = J
            END IF
 65      CONTINUE
         DO 70 K = 1, IATOM
            IF ((IBOND(K,ITRI) .GT. 0) .AND. (K .NE. II)
     &                                   .AND. (K .NE. JJ)) THEN
               VEC3(1) = ATMARR(II,2)-ATMARR(K,2)
               VEC3(2) = ATMARR(II,3)-ATMARR(K,3)
               VEC3(3) = ATMARR(II,4)-ATMARR(K,4)
               VEC4(1) = ATMARR(JJ,2)-ATMARR(K,2)
               VEC4(2) = ATMARR(JJ,3)-ATMARR(K,3)
               VEC4(3) = ATMARR(JJ,4)-ATMARR(K,4)
               BNDAN1  = VECANG(VEC3,VEC1)
               BNDAN2  = VECANG(VEC4,VEC2)
               IF (ABS(BNDAN1) .GT. DEG175) THEN
                  KK = K
                  KBND = JJ
               ELSE IF (ABS(BNDAN2) .GT. DEG175) THEN
                  KK = K
                  KBND = II
               ELSE
                  XNRM3 = VEC3(1)*VEC3(1) + VEC3(2)*VEC3(2) +
     &                 VEC3(3)*VEC3(3)
                  XNRM4 = VEC4(1)*VEC4(1) + VEC4(2)*VEC4(2) +
     &                 VEC4(3)*VEC4(3)
                  IF (XNRM3 .LT. XNRM4) THEN
                     KK = K
                     KBND = II
                  ELSE
                     KK = K
                     KBND = JJ
                  END IF
               END IF
            END IF        
 70      CONTINUE
C
         IINTCR = IINTCR + 1
         INTCRD(IINTCR,1) = 23
         IF (KBND .EQ. II) THEN
            INTCRD(IINTCR,2) = JJ
            INTCRD(IINTCR,3) = ITRI
            INTCRD(IINTCR,4) = II
            INTCRD(IINTCR,5) = KK
         ELSE
            INTCRD(IINTCR,2) = II
            INTCRD(IINTCR,3) = ITRI
            INTCRD(IINTCR,4) = JJ
            INTCRD(IINTCR,5) = KK
         END IF
C
         IF ((IBEFOR .NE. IINTCR) .AND. (IPRINT .GE. IPRMAX)) THEN
            CALL HEADER('Extra dihedral angles',-1)
            DO 90 I = IBEFOR + 1, IINTCR
               WRITE(LUPRI,'(A,4I4)')
     &              ' Dihedral angles between atoms: ',
     &              INTCRD(I,2),INTCRD(I,3),INTCRD(I,4),INTCRD(I,5)
 90         CONTINUE
         END IF
      END IF
      RETURN
      END

C  /* Deck getwil */
      SUBROUTINE GETWIL(ATMARR,WILBMT)
C     
C     Constructs the Wilson B matrix used for transformations between
C     Cartesian and redundant natural internal coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMARR(MXCENT,8), WILBMT(2*MXCOOR,MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175D0*PI/180D0)
C
      CALL DZERO(WILBMT,2*MXCOOR*MXCOOR)
      CALL DZERO(ATMARR,8*MXCENT)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from GETWIL','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Expanded array of atoms (Bohr)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C      
      DO 20 IC = 1, IINTCR
C
C     The connection between Cartesian coordinates and bonds:
C     -------------------------------------------------------
C
         IF (INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            CALL NRMLVC(VEC1)
            DO 25 II = 1, 3 
               WILBMT(IC,(INTCRD(IC,2)-1)*3+II) =  VEC1(II)
               WILBMT(IC,(INTCRD(IC,3)-1)*3+II) = -VEC1(II)
 25         CONTINUE
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               DO 27 II = 1, IATOM
                  DO 28 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     BNDL2 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     CHANGE = (BNDL2-BNDL1)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,CHANGE
 28               CONTINUE
 27            CONTINUE
            END IF
C
C     The connection between Cartesian coordinates and angles:
C     --------------------------------------------------------
C
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
            VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
            VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            CALL NRMLVC(VEC1)
            CALL NRMLVC(VEC2)
C
C     All regular angles
C
            IF (INTCRD(IC,1) .EQ. 11) THEN
               CALL VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0D-8) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               CALL NRMLVC(VEC3)
C
C     Second coordinate of angles larger than 175 degrees.
C
            ELSE
               CALL VECPRD(VEC1,VEC3,VEC4)
               CALL NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
C
            CALL VECPRD(VEC1,VEC3,VEC4)
            CALL VECPRD(VEC3,VEC2,VEC5)
            CALL NRMLVC(VEC4)
            CALL NRMLVC(VEC5)
            DO 30 II = 1, 3
               WILBMT(IC,(INTCRD(IC,2)-1)*3+II) = VEC4(II)/BNDL1
               WILBMT(IC,(INTCRD(IC,3)-1)*3+II) =
     &                          -VEC4(II)/BNDL1 - VEC5(II)/BNDL2
               WILBMT(IC,(INTCRD(IC,4)-1)*3+II) = VEC5(II)/BNDL2
 30         CONTINUE
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
               VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
               VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
               CALL VECPRD(VEC3,VEC1,VEC4)
               CALL NRMLVC(VEC4)
               ANGLE = VECANG(VEC1,VEC2)
               IF ((INTCRD(IC,1) .EQ. 12) .OR. (INTCRD(IC+1,1) .EQ. 12))
     &              ANGLE = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               DO 47 II = 1, IATOM
                  DO 48 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     ANGLE2 = VECANG(VEC1,VEC2)
                     IF ((INTCRD(IC,1) .EQ. 12) .OR.
     &                  (INTCRD(IC+1,1) .EQ. 12))
     &                    ANGLE2 = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
                     CHANGE = (ANGLE2-ANGLE)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,CHANGE
 48               CONTINUE
 47            CONTINUE
            END IF
C
C     The connection between Cartesian coordinates and dihedral angles:
C     -----------------------------------------------------------------
C
         ELSE
            VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
            VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
            VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
            VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
            VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
            VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
            VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
            VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
            VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
            CALL NRMLVC(VEC1)
            CALL NRMLVC(VEC2)
            CALL NRMLVC(VEC3)
            UW =  DDOT(3,VEC1,1,VEC2,1)
            VW =  DDOT(3,VEC3,1,VEC2,1)
            CALL VECPRD(VEC1,VEC2,VEC4)
            CALL VECPRD(VEC3,VEC2,VEC5)
C
C     Check for linearity. The dihedral angle is undefined if three
C     or four of the atoms are linear.
C
            IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0D-16) .AND.
     &           (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0D-16)) THEN
               DO 50 II = 1, 3
C
C     We construct a few "building blocks".
C
                  CMPU1 = VEC4(II)/(BNDL1*(1.0D0-UW*UW))
                  CMPV1 = VEC5(II)/(BNDL3*(1.0D0-VW*VW))
                  CMPU2 = VEC4(II)*UW/(BNDL2*(1.0D0-UW*UW))
                  CMPV2 = VEC5(II)*VW/(BNDL2*(1.0D0-VW*VW))
                  WILBMT(IC,(INTCRD(IC,2)-1)*3+II) =  CMPU1
                  WILBMT(IC,(INTCRD(IC,3)-1)*3+II) = -CMPU1+CMPU2-CMPV2
                  WILBMT(IC,(INTCRD(IC,4)-1)*3+II) =  CMPV1+CMPV2-CMPU2
                  WILBMT(IC,(INTCRD(IC,5)-1)*3+II) = -CMPV1
 50            CONTINUE
            END IF
C
            IF (IPRINT .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
               VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
               VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
               VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
               VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
               VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
               VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
               VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
               VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
               CALL NRMLVC(VEC2)
               CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
               VEC1(1) = VEC1(1) - CMPNT1*VEC2(1) 
               VEC1(2) = VEC1(2) - CMPNT1*VEC2(2) 
               VEC1(3) = VEC1(3) - CMPNT1*VEC2(3) 
               CALL NRMLVC(VEC1)
               CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
               VEC3(1) = VEC3(1) - CMPNT2*VEC2(1) 
               VEC3(2) = VEC3(2) - CMPNT2*VEC2(2) 
               VEC3(3) = VEC3(3) - CMPNT2*VEC2(3) 
               CALL NRMLVC(VEC3)
               DIHED = VECANG(VEC1,VEC3)
               IF (ABS(DIHED) .GT. DEG175) THEN
                  CALL VECPRD(VEC1,VEC2,VEC4)
                  CALL NRMLVC(VEC4)
                  DIHED = VECANG(VEC1,VEC4)+VECANG(VEC4,VEC3)
               END IF
               VEC5(1) = VEC1(1)
               VEC5(2) = VEC1(2)
               VEC5(3) = VEC1(3)
               VEC6(1) = VEC3(1)
               VEC6(2) = VEC3(2)
               VEC6(3) = VEC3(3)
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               FAC = -D1
               DO 97 II = 1, IATOM
                  DO 98 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3) 
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4) 
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2) 
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3) 
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4) 
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(INTCRD(IC,5),2) 
     &                    - ATMARR(INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(INTCRD(IC,5),3) 
     &                    - ATMARR(INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(INTCRD(IC,5),4) 
     &                    - ATMARR(INTCRD(IC,4),4)
                     CALL NRMLVC(VEC2)
                     CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
                     VEC1(1) = VEC1(1) - CMPNT1*VEC2(1) 
                     VEC1(2) = VEC1(2) - CMPNT1*VEC2(2) 
                     VEC1(3) = VEC1(3) - CMPNT1*VEC2(3) 
                     CALL NRMLVC(VEC1)
                     CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
                     VEC3(1) = VEC3(1) - CMPNT2*VEC2(1) 
                     VEC3(2) = VEC3(2) - CMPNT2*VEC2(2) 
                     VEC3(3) = VEC3(3) - CMPNT2*VEC2(3) 
                     CALL NRMLVC(VEC3)
                     IF (FAC .LT. D0) THEN
                        FAC = D1
                        CALL VECPRD(VEC1,VEC3,VEC4)
                        IF (DDOT(3,VEC4,1,VEC2,1) .LT. D0) FAC = -D1
                     END IF
                     DIHED2 = VECANG(VEC1,VEC3)
                     CH1 = FAC*(VECANG(VEC1,VEC6)
     &                    -VECANG(VEC5,VEC6))/1.0D-6
                     CH2 = FAC*(VECANG(VEC3,VEC5)
     &                    -VECANG(VEC6,VEC5))/1.0D-6
                     IF (ABS(DIHED2) .GT. DEG175) THEN
                        CALL VECPRD(VEC1,VEC2,VEC4)
                        CALL NRMLVC(VEC4)
                        DIHED2 = VECANG(VEC1,VEC4)+VECANG(VEC4,VEC3)
                     END IF
                     CHANGE = FAC*(DIHED2-DIHED)/1.0D-6
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0D-6
                     WRITE(LUPRI,'(2I5,3F16.6)') II,JJ,CHANGE,CH1,CH2                     
 98               CONTINUE
 97            CONTINUE
            END IF
         END IF
C
         IF (IPRINT .GE. IPRMAX) THEN
            TXT = 'Internal coordinate number:     '
            WRITE(TXT(28:32),'(I5)') IC
            CALL HEADER(TXT,-1)
            DO 99 K = 1, IATOM
               WRITE(LUPRI,'(I25,A,F16.6)') 
     &              K,'x',WILBMT(IC,(K-1)*3+1)
               WRITE(LUPRI,'(I25,A,F16.6)') 
     &              K,'y',WILBMT(IC,(K-1)*3+2)
               WRITE(LUPRI,'(I25,A,F16.6)') 
     &              K,'z',WILBMT(IC,(K-1)*3+3)
 99         CONTINUE
         END IF
 20   CONTINUE
C
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Wilson B matrix [B(ij) = dq(i)/dx(j)]',-1)
         CALL OUTPUT(WILBMT,1,IINTCR,1,ICRTCR,2*MXCOOR,MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck getdwl */
      SUBROUTINE GETDWL(DWILBM,ATMARR,TMPMAT,WILBMT)
C     
C     Constructs the derivative of the Wilson B matrix used for transformations
C     between Cartesian and redundant natural internal coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
      DIMENSION DWILBM(MXCOOR,MXCOOR), WILBMT(2*MXCOOR,MXCOOR)
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL DZERO(ATMARR,8*MXCENT)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL TITLER('Output from GETDWL','*',103)
      END IF
C
C     First we initialize the ATMARR array.
C
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Expanded array of atoms (Bohr)',-1)
         CALL OUTPUT(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      ICRTCR = 3*IATOM
C
      DO 20 IC = 1, IINTCR
         CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
         CALL DZERO(DWILBM,MXCOOR*MXCOOR)
         IF (INTCRD(IC,1) .LT. 10) THEN
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 30 II = 2, 3
                  DO 35 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1)= ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2)= ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3)= ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     CALL NRMLVC(VEC1)
                     CALL DZERO(VEC2,3)
                     DO 40 J = 1, 3 
                        VEC2(J) = 1.0D0
                        COMPNT = DDOT(3,VEC1,1,VEC2,1)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,3)-1)*3+J) = -COMPNT
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,2)-1)*3+J) =  COMPNT
                        VEC2(J) = 0.0D0
 40                  CONTINUE
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 35               CONTINUE
 30            CONTINUE
               DO 50 II = 2, 3
                  DO 51 JJ = 1, 3
                     DO 53 KK = 2, 3
                        DO 54 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6                           
 54                     CONTINUE
 53                  CONTINUE
 51               CONTINUE
 50            CONTINUE
               DO 56 II = 1, ICRTCR
                  DO 58 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 58               CONTINUE
 56            CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C     
C     Analytically:
C
            CALL GTDWLM(IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 60 II = 2, 4
                  DO 65 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     CALL NRMLVC(VEC1)
                     CALL NRMLVC(VEC2)
C     
C     All regular angles
C     
                     IF (INTCRD(IC,1) .EQ. 11) THEN
                        CALL VECPRD(VEC1,VEC2,VEC3)
                        VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
                        IF (VNRM .LE. 1.0D-8) THEN
                           VEC3(1) =  VEC1(2)+VEC1(3)
                           VEC3(2) = -VEC1(1)+VEC1(3)
                           VEC3(3) = -VEC1(1)-VEC1(2)
                        END IF
                        CALL NRMLVC(VEC3)
C     
C     Second coordinate of angles larger than 175 degrees.
C     
                     ELSE
                        CALL VECPRD(VEC1,VEC3,VEC4)
                        CALL NRMLVC(VEC4)
                        VEC3(1) = VEC4(1)
                        VEC3(2) = VEC4(2)
                        VEC3(3) = VEC4(3)
                     END IF
C     
                     CALL VECPRD(VEC1,VEC3,VEC4)
                     CALL NRMLVC(VEC4)
                     CALL VECPRD(VEC3,VEC2,VEC5)
                     CALL NRMLVC(VEC5)
C     
                     DO 70 J = 1, 3
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,2)-1)*3+J) =  (VEC4(J)/BNDL1)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,3)-1)*3+J) = -(VEC4(J)/BNDL1)
     &                       - (VEC5(J)/BNDL2)
                        TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                       (INTCRD(IC,4)-1)*3+J) =  (VEC5(J)/BNDL2)
 70                  CONTINUE
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 65               CONTINUE
 60            CONTINUE
C     
               DO 90 II = 2, 4
                  DO 91 JJ = 1, 3
                     DO 93 KK = 2, 4
                        DO 94 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6
 94                     CONTINUE
 93                  CONTINUE
 91               CONTINUE
 90            CONTINUE
               DO 96 II = 1, ICRTCR
                  DO 98 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 98               CONTINUE
 96            CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C     
C     Analytically:
C
            CALL GTDWLM(IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
         ELSE
C
C     Numerically:
C
            IF (IPRINT .GE. IPRDBG) THEN
               DO 100 II = 2, 5
                  DO 105 JJ = 1, 3
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) + 1.0D-6
                     VEC1(1) = ATMARR(INTCRD(IC,2),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(INTCRD(IC,2),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(INTCRD(IC,2),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(INTCRD(IC,4),2)
     &                    - ATMARR(INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(INTCRD(IC,4),3)
     &                    - ATMARR(INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(INTCRD(IC,4),4)
     &                    - ATMARR(INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(INTCRD(IC,5),2) 
     &                    - ATMARR(INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(INTCRD(IC,5),3)
     &                    - ATMARR(INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(INTCRD(IC,5),4)
     &                    - ATMARR(INTCRD(IC,4),4)
                     VNRM1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     VNRM2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     VNRM3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
                     CALL NRMLVC(VEC1)
                     CALL NRMLVC(VEC2)
                     CALL NRMLVC(VEC3)
                     UW =  DDOT(3,VEC1,1,VEC2,1)
                     VW =  DDOT(3,VEC3,1,VEC2,1)
                     CALL VECPRD(VEC1,VEC2,VEC4)
                     CALL VECPRD(VEC2,VEC3,VEC5)
                     CALL VECPRD(VEC3,VEC2,VEC6)
C     
C     Check for linearity. The dihedral angle is undefined if three
C     or four of the atoms are linear.
C     
                     IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0D-16) .AND.
     &                    (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0D-16)) THEN
                        DO 120 J = 1, 3
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,2)-1)*3+J)
     &                          = VEC4(J)/(VNRM1*(1.0D0-UW*UW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,3)-1)*3+J) =
     &                          (VEC4(J)*(VNRM1*UW-VNRM2))/
     &                          (VNRM1*VNRM2*(1.0D0-UW*UW))
     &                          - (VW*VEC6(J))/(VNRM2*(1.0D0-VW*VW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,4)-1)*3+J) =
     &                          (VEC6(J)*(VNRM3*VW+VNRM2))/
     &                          (VNRM3*VNRM2*(1.0D0-VW*VW))
     &                          - (UW*VEC4(J))/(VNRM2*(1.0D0-UW*UW))
                           TMPMAT((INTCRD(IC,II)-1)*3+JJ,
     &                          (INTCRD(IC,5)-1)*3+J) =
     &                          VEC5(J)/(VNRM3*(1.0D0-VW*VW))
 120                    CONTINUE
                     END IF
                     ATMARR(INTCRD(IC,II),JJ+1) =
     &                    ATMARR(INTCRD(IC,II),JJ+1) - 1.0D-6
 105              CONTINUE
 100           CONTINUE
C     
               DO 200 II = 2, 5
                  DO 210 JJ = 1, 3
                     DO 220 KK = 2, 5
                        DO 230 LL = 1, 3
                           DWILBM((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) =
     &                          (TMPMAT((INTCRD(IC,KK)-1)*3+LL,
     &                          (INTCRD(IC,II)-1)*3+JJ) -
     &                          WILBMT(IC,(INTCRD(IC,II)-1)*3+JJ))
     &                          /1.0D-6
 230                    CONTINUE
 220                 CONTINUE
 210              CONTINUE
 200           CONTINUE
               DO 250 II = 1, ICRTCR
                  DO 255 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0D0
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0D-6)
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 255              CONTINUE
 250           CONTINUE
               CALL HEADER('Numerically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
C
C     Analytically:
C
            CALL GTDWLM(IC,DWILBM,ATMARR,WILBMT)
            IF (IPRINT .GE. IPRMAX) THEN
               CALL HEADER('Analytically diff. Wilson B matrix',-1)
               CALL OUTPUT(DWILBM,1,ICRTCR,1,ICRTCR,MXCOOR,
     &              MXCOOR,1,LUPRI)
            END IF
         END IF
 20   CONTINUE
      RETURN
      END

C  /* Deck gtdwlm */
      SUBROUTINE GTDWLM(IC,DWILBM,ATMARR,WILBMT)
C
C     Returns the derivative of the Wilson B matrix for _one_
C     redundant internal coordinate (IC).
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <priunit.h>
      DIMENSION DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR)
      DIMENSION VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(3)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      CALL DZERO(DWILBM,MXCOOR*MXCOOR)
      CALL ATMINI(ATMARR,IATOM,.TRUE.)
C
C     The derivative for bonds:
C     -------------------------
C
      IF (INTCRD(IC,1) .LT. 10) THEN
         VEC1(1)=ATMARR(INTCRD(IC,2),2)-ATMARR(INTCRD(IC,3),2)
         VEC1(2)=ATMARR(INTCRD(IC,2),3)-ATMARR(INTCRD(IC,3),3)
         VEC1(3)=ATMARR(INTCRD(IC,2),4)-ATMARR(INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         CALL NRMLVC(VEC1)
         DO 10 II = 1, 3
            DO 15 JJ = 1, 3
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms.
C
               FACIJ = D0
               IF (II .EQ. JJ) FACIJ = D1
               CMPUU = (FACIJ - VEC1(II)*VEC1(JJ))/BNDL1
               DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &              =  CMPUU
               DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &              = -CMPUU
               DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &          = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
               DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &              =  CMPUU
 15         CONTINUE
 10      CONTINUE
C
C     The derivative for angles:
C     --------------------------
C
      ELSE IF (INTCRD(IC,1) .LT. 20) THEN
         VEC1(1)=ATMARR(INTCRD(IC,2),2)-ATMARR(INTCRD(IC,3),2)
         VEC1(2)=ATMARR(INTCRD(IC,2),3)-ATMARR(INTCRD(IC,3),3)
         VEC1(3)=ATMARR(INTCRD(IC,2),4)-ATMARR(INTCRD(IC,3),4)
         VEC2(1)=ATMARR(INTCRD(IC,4),2)-ATMARR(INTCRD(IC,3),2)
         VEC2(2)=ATMARR(INTCRD(IC,4),3)-ATMARR(INTCRD(IC,3),3)
         VEC2(3)=ATMARR(INTCRD(IC,4),4)-ATMARR(INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         CALL NRMLVC(VEC1)
         CALL NRMLVC(VEC2)
         COSA = DDOT(3,VEC1,1,VEC2,1)
         SINA = SQRT(1-COSA*COSA)
C
C     We can only calculate the derivative, if sin(q) is non-zero, as
C     all expressions contains a division by sin(q).
C
         IF (ABS(SINA) .GT. 1.0D-6) THEN
            DO 25 II = 1, 3
               DADMI = WILBMT(IC,(INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(INTCRD(IC,4)-1)*3+II)
               DADOI = WILBMT(IC,(INTCRD(IC,3)-1)*3+II)
               DO 27 JJ = 1, 3
                  DADMJ = WILBMT(IC,(INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(INTCRD(IC,4)-1)*3+JJ)
                  DADOJ = WILBMT(IC,(INTCRD(IC,3)-1)*3+JJ)
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms.
C
                  FACIJ = D0
                  IF (II .EQ. JJ) FACIJ = D1
C     
C     We only need a few "building blocks" to make all expressions.
C
                  CMPUU  = (VEC1(II)*VEC2(JJ) + VEC1(JJ)*VEC2(II)
     &                 - 3.0D0*COSA*VEC1(II)*VEC1(JJ)
     &                 + FACIJ*COSA)/(SINA*BNDL1*BNDL1)
                  CMPVV  = (VEC2(II)*VEC1(JJ) + VEC2(JJ)*VEC1(II)
     &                 - 3.0D0*COSA*VEC2(II)*VEC2(JJ)
     &                 + FACIJ*COSA)/(SINA*BNDL2*BNDL2)
                  CMPWW1 = (VEC1(II)*VEC1(JJ) + VEC2(II)*VEC2(JJ)
     &                 - COSA*VEC1(II)*VEC2(JJ)
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
                  CMPWW2 = (VEC2(II)*VEC2(JJ) + VEC1(II)*VEC1(JJ)
     &                 - COSA*VEC2(II)*VEC1(JJ)
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
C     
C     All the different elements are constructed.
C
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &             = CMPUU - DADMI*DADMJ*COSA/SINA
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = -CMPUU - CMPWW1 - DADMI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &             = CMPWW1 - DADMI*DADNJ*COSA/SINA
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = CMPUU + CMPVV + CMPWW1 + CMPWW2
     &             - DADOI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &             = -CMPVV - CMPWW2 - DADNI*DADOJ*COSA/SINA
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,4)-1)*3+II)
     &             = DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &             = CMPVV - DADNI*DADNJ*COSA/SINA
 27            CONTINUE
 25         CONTINUE
C
C     If the second derivative is undefined (sin(q) ~ 0), we set the matrix equal to
C     the identity matrix.
C     
         ELSE
            DO 30 I = 1, ICRTCR
               DWILBM(I,I) = D1
 30         CONTINUE
         END IF
C
C     The derivative for dihedral angles:
C     -----------------------------------
C
      ELSE 
         VEC1(1) = ATMARR(INTCRD(IC,2),2) - ATMARR(INTCRD(IC,3),2)
         VEC1(2) = ATMARR(INTCRD(IC,2),3) - ATMARR(INTCRD(IC,3),3)
         VEC1(3) = ATMARR(INTCRD(IC,2),4) - ATMARR(INTCRD(IC,3),4)
         VEC2(1) = ATMARR(INTCRD(IC,4),2) - ATMARR(INTCRD(IC,3),2)
         VEC2(2) = ATMARR(INTCRD(IC,4),3) - ATMARR(INTCRD(IC,3),3)
         VEC2(3) = ATMARR(INTCRD(IC,4),4) - ATMARR(INTCRD(IC,3),4)
         VEC3(1) = ATMARR(INTCRD(IC,5),2) - ATMARR(INTCRD(IC,4),2)
         VEC3(2) = ATMARR(INTCRD(IC,5),3) - ATMARR(INTCRD(IC,4),3)
         VEC3(3) = ATMARR(INTCRD(IC,5),4) - ATMARR(INTCRD(IC,4),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
         CALL NRMLVC(VEC1)
         CALL NRMLVC(VEC2)
         CALL NRMLVC(VEC3)
         CALL VECPRD(VEC1,VEC2,VEC4)
         CALL VECPRD(VEC3,VEC2,VEC5)
         UW   = DDOT(3,VEC1,1,VEC2,1)
         VW   = DDOT(3,VEC3,1,VEC2,1)
         UWVW = DDOT(3,VEC4,1,VEC5,1)
         UVW  = DDOT(3,VEC1,1,VEC5,1)
C
C     We can only calculate the derivative, if sin(q) (or (u * (v x w)))
C     is non-zero, as all expressions contains a division by sin(q).
C
         IF (ABS(UVW) .GT. 1.0D-6) THEN
            DO 40 II = 1, 3
               DADMI = WILBMT(IC,(INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(INTCRD(IC,5)-1)*3+II)
               DO 45 JJ = 1, 3
                  DADMJ = WILBMT(IC,(INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(INTCRD(IC,5)-1)*3+JJ)
C
C     The variable FACIJ makes the distinction between "diagonal" and
C     "off-diagonal" elements, where "diagonal" refers to Cartesian
C     displacements (II = JJ) and not the atoms. The sign of this
C     factor follows this table:
C
C                                 i
C                              1  2  3
C                             ---------
C                           1| 0  +  - 
C                         j 2| -  0  + 
C                           3| +  -  0 
C
                  FACIJ = (JJ-II)*(-0.5D0)**(ABS(JJ-II))
C
C     For some of the elements we need the third Cartesian direction,
C     that is k =/ i,j.
C
                  KK = 6-II-JJ
C     
C     With the following "building blocks" all expressions can be constructed:
C
                  CMPUU  = (VEC4(II)*(UW*VEC2(JJ)-VEC1(JJ))
     &                   + VEC4(JJ)*(UW*VEC2(II)-VEC1(II)))/
     &                     (BNDL1*BNDL1*(1.0D0-UW*UW)**2)
                  CMPVV  = (VEC5(II)*(VW*VEC2(JJ)-VEC3(JJ))
     &                   + VEC5(JJ)*(VW*VEC2(II)-VEC3(II)))/
     &                     (BNDL3*BNDL3*(1.0D0-VW*VW)**2)
                  CMPUW  = 0.5D0*(VEC4(II)*(VEC2(JJ)-2.0D0*UW*VEC1(JJ)
     &                   + UW*UW*VEC2(JJ))+VEC4(JJ)*(VEC2(II)
     &                   - 2.0D0*UW*VEC1(II)+UW*UW*VEC2(II)))/
     &                     (BNDL1*BNDL2*(1.0D0-UW*UW)**2)
                  CMPVW  = 0.5D0*(VEC5(II)*(VEC2(JJ)-2.0D0*VW*VEC3(JJ)
     &                   + VW*VW*VEC2(JJ))+VEC5(JJ)*(VEC2(II)
     &                   - 2.0D0*VW*VEC3(II)+VW*VW*VEC2(II)))/
     &                     (BNDL3*BNDL2*(1.0D0-VW*VW)**2)
                  CMPWW1 = 0.5D0*(VEC4(II)*(VEC1(JJ)+UW*UW*VEC1(JJ)
     &                   - 3.0D0*UW*VEC2(JJ)+UW*UW*UW*VEC2(JJ))
     &                   + VEC4(JJ)*(VEC1(II)+UW*UW*VEC1(II)
     &                   - 3.0D0*UW*VEC2(II)+UW*UW*UW*VEC2(II)))/
     &                     (BNDL2*BNDL2*(1.0D0-UW*UW)**2)
                  CMPWW2 = 0.5D0*(VEC5(II)*(VEC3(JJ)+VW*VW*VEC3(JJ)
     &                   - 3.0D0*VW*VEC2(JJ)+VW*VW*VW*VEC2(JJ))
     &                   + VEC5(JJ)*(VEC3(II)+VW*VW*VEC3(II)
     &                   - 3.0D0*VW*VEC2(II)+VW*VW*VW*VEC2(II)))/
     &                     (BNDL2*BNDL2*(1.0D0-VW*VW)**2)
                  CMPKK1 = (UW*VEC2(KK)-VEC1(KK))/
     &                     (BNDL1*BNDL2*(1.0D0-UW*UW))
                  CMPKK2 = (VW*VEC2(KK)-VEC3(KK))/
     &                     (BNDL3*BNDL2*(1.0D0-VW*VW))
C
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,2)-1)*3+JJ)
     &                 = CMPUU
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = -CMPUU + CMPUW + FACIJ*CMPKK1
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = -CMPUW - FACIJ*CMPKK1
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,2)-1)*3+II)
     &             = DWILBM((INTCRD(IC,2)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,5)-1)*3+JJ)
     &                 = -CMPVV
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = -CMPVW - FACIJ*CMPKK2
                  DWILBM((INTCRD(IC,3)-1)*3+JJ,(INTCRD(IC,5)-1)*3+II)
     &             = DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = CMPVV + CMPVW + FACIJ*CMPKK2
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,5)-1)*3+II)
     &             = DWILBM((INTCRD(IC,5)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,3)-1)*3+JJ)
     &                 = CMPUU - 2.0D0*CMPUW - CMPWW1 + CMPWW2
                  DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = CMPUW + CMPVW + CMPWW1 - CMPWW2
     &                 + FACIJ*(CMPKK1-CMPKK2)
                  DWILBM((INTCRD(IC,4)-1)*3+JJ,(INTCRD(IC,3)-1)*3+II)
     &             = DWILBM((INTCRD(IC,3)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((INTCRD(IC,4)-1)*3+II,(INTCRD(IC,4)-1)*3+JJ)
     &                 = -CMPVV - 2.0D0*CMPVW + CMPWW2 - CMPWW1
 45            CONTINUE
 40         CONTINUE
C
C     If the second derivative is undefined (sin(q) ~ 0), we set the
C     matrix equal to the identity matrix.
C     
         ELSE
            DO 50 I = 1, ICRTCR
               DWILBM(I,I) = D1
 50         CONTINUE
         END IF
      END IF
      RETURN
      END

C  /* Deck gtbinv */
      SUBROUTINE GTBINV(UMAT,WMAT,VMAT,VAL,WILBMT,BMTINV,PJINMT)
C     
C     Decomposes the transpose of Wilson's B matrix, then one can
C     construct the inverse of the rectangulat matrix.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
      DIMENSION UMAT(MXCOOR,2*MXCOOR), WMAT(2*MXCOOR,2*MXCOOR)
      DIMENSION VMAT(2*MXCOOR,2*MXCOOR), VAL(2*MXCOOR)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR), BMTINV(2*MXCOOR,MXCOOR),
     &          PJINMT(2*MXCOOR,2*MXCOOR)
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      IF (IPRINT .GE. IPRMED) CALL TITLER('Output from GTBINV','*',103)
      IDIM1 = ICRTCR
      IDIM2 = IINTCR
      CALL DZERO(BMTINV,2*MXCOOR*MXCOOR)
      CALL DZERO(UMAT,2*MXCOOR*MXCOOR)
      CALL DZERO(WMAT,4*MXCOOR*MXCOOR)
      CALL DZERO(VMAT,4*MXCOOR*MXCOOR)
      CALL DZERO(VAL,2*MXCOOR)
C
C     The transpose is placed in UMAT
C
      DO 10 I = 1, IINTCR
         DO 12 J = 1, ICRTCR
            UMAT(J,I) = WILBMT(I,J)
 12      CONTINUE
 10   CONTINUE
C
C     The transpose of Wilsons B matrix is decomposed to:
C
C           B^t = U W V^t  ;   U,V orthogonal, W diagonal
C
      BMTNRM = D0
      SCL    = D0
      SGNFAC = D0
      DO 20 I = 1, IDIM2
         VAL(I) = SCL*SGNFAC
         SCL    = D0
         SGNFAC = D0
         SUM    = D0
         II = I + 1
         IF (I .LE. IDIM1) THEN 
            DO 25 J = I, IDIM1 
               SCL = SCL + ABS(UMAT(J,I))
 25         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 30 J = I, IDIM1
                  UMAT(J,I) = UMAT(J,I)/SCL
                  SUM = SUM + UMAT(J,I)*UMAT(J,I)
 30            CONTINUE
               DIAG = UMAT(I,I)
               SGNFAC = -SIGN(SQRT(SUM),DIAG)
               FAC = SGNFAC*DIAG - SUM
               UMAT(I,I) = DIAG - SGNFAC
               DO 35 J = II, IDIM2
                  SUM = D0
                  DO 37 K = I, IDIM1
                     SUM = SUM + UMAT(K,I)*UMAT(K,J)
 37               CONTINUE
                  FAC2 = SUM/FAC
                  DO 39 K = I, IDIM1
                     UMAT(K,J) = UMAT(K,J)+FAC2*UMAT(K,I)
 39               CONTINUE
 35            CONTINUE
               DO 40 J = I, IDIM1
                  UMAT(J,I) = SCL*UMAT(J,I)
 40            CONTINUE
            END IF
         END IF
         WMAT(I,I) = SCL*SGNFAC
         SUM    = D0
         SGNFAC = D0
         SCL    = D0
         IF ((I .LE. IDIM1) .AND. (I .NE. IDIM2)) THEN
            DO 45 J = II, IDIM2
               SCL = SCL + ABS(UMAT(I,J))
 45         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 50 J = II, IDIM2
                  UMAT(I,J) = UMAT(I,J)/SCL
                  SUM = SUM + UMAT(I,J)*UMAT(I,J)
 50            CONTINUE
               FAC = UMAT(I,II)
               SGNFAC = -SIGN(SQRT(SUM),FAC)
               FAC2 = SGNFAC*FAC - SUM
               UMAT(I,II) = FAC - SGNFAC
               DO 52 J = II, IDIM2
                  VAL(J) = UMAT(I,J)/FAC2
 52            CONTINUE
               DO 54 J = II, IDIM1
                  SUM = D0
                  DO 56 K = I+1, IDIM2
                     SUM = SUM + UMAT(J,K)*UMAT(I,K)
 56               CONTINUE
                  DO 58 K = II, IDIM2
                     UMAT(J,K) = UMAT(J,K) + SUM*VAL(K)
 58               CONTINUE
 54            CONTINUE
               DO 60 J = II, IDIM2
                  UMAT(I,J) = SCL*UMAT(I,J)
 60            CONTINUE
            END IF
         END IF
         BMTNRM = MAX(BMTNRM,(ABS(WMAT(I,I)) + ABS(VAL(I))))
 20   CONTINUE
      DO 70 I = IDIM2, 1, -1
         IF (I .LT. IDIM2) THEN
            IF (SGNFAC .NE. D0) THEN
               DO 72 J = II, IDIM2
                  VMAT(J,I) = UMAT(I,J)/(SGNFAC*UMAT(I,II))
 72            CONTINUE
               DO 74 J = II, IDIM2
                  SUM = D0
                  DO 76 K = II, IDIM2
                     SUM = SUM + UMAT(I,K)*VMAT(K,J)
 76               CONTINUE
                  DO 78 K = II, IDIM2
                     VMAT(K,J) = VMAT(K,J) + SUM*VMAT(K,I)
 78               CONTINUE
 74            CONTINUE
            END IF
            DO 79 J = II, IDIM2
               VMAT(I,J) = D0
               VMAT(J,I) = D0
 79         CONTINUE
         END IF
         VMAT(I,I) = 1.0D0
         SGNFAC = VAL(I)
         II = I
 70   CONTINUE
      DO 80 I = MIN(IDIM1,IDIM2), 1, -1
         II = I + 1
         SGNFAC = WMAT(I,I)
         DO 82 J = II, IDIM2
            UMAT(I,J) = D0
 82      CONTINUE
         IF (SGNFAC .NE. D0) THEN
            SGNFAC = 1.0D0/SGNFAC
            DO 84 J = II, IDIM2
               SUM = D0
               DO 86 K = II, IDIM1
                  SUM = SUM + UMAT(K,I)*UMAT(K,J)
 86            CONTINUE
               FAC = SGNFAC*(SUM/UMAT(I,I))
               DO 88 K = I, IDIM1
                  UMAT(K,J) = UMAT(K,J) + FAC*UMAT(K,I)
 88            CONTINUE
 84         CONTINUE
            DO 90 J = I, IDIM1
               UMAT(J,I) = SGNFAC*UMAT(J,I)
 90         CONTINUE
         ELSE
            DO 92 J = I, IDIM1
               UMAT(J,I) = D0
 92         CONTINUE
         END IF
         UMAT(I,I) = UMAT(I,I) + 1.0D0
 80   CONTINUE
      DO 100 J = IDIM2, 1, -1
         DO 105 ITRS = 1, 50
            DO 110 K = J, 1, -1
               KK = K - 1
               IF ((ABS(VAL(K))+BMTNRM) .EQ. BMTNRM) GOTO 902
               IF ((ABS(WMAT(KK,KK))+BMTNRM) .EQ. BMTNRM) GOTO 901
 110        CONTINUE
 901        CONTINUE
            FC   = D0
            SUM = 1.0D0
            DO 112  I = K, J
               FAC = SUM*VAL(I)
               VAL(I) = FC*VAL(I)
               IF ((ABS(FAC)+BMTNRM) .EQ. BMTNRM) GOTO 902
               SGNFAC = WMAT(I,I)
               FAC2 = SQRT(FAC*FAC + SGNFAC*SGNFAC)
               WMAT(I,I) = FAC2
               FAC2 = 1.0D0/FAC2
               FC  = SGNFAC*FAC2
               SUM =   -FAC*FAC2
               DO 114 L = 1, IDIM1
                  FY = UMAT(L,KK)
                  FZ = UMAT(L,I)
                  UMAT(L,KK) = FY*FC + FZ*SUM
                  UMAT(L,I)  = FZ*FC - FY*SUM  
 114           CONTINUE
 112        CONTINUE
 902        CONTINUE
            FZ = WMAT(J,J)
            IF (K .EQ. J) THEN
               IF (FZ .LT. D0) THEN
                  WMAT(J,J) = -FZ
                  DO 116 L = 1, IDIM2
                     VMAT(L,J) = -VMAT(L,J)
 116              CONTINUE
               END IF
               GOTO 903
            END IF
            IF (ITRS .EQ. 50) CALL QUIT('No convergence in GTBINV.')
            FX = WMAT(K,K)
            JJ = J - 1
            FY = WMAT(JJ,JJ)
            SGNFAC = VAL(JJ)
            FAC2 = VAL(J)
            FAC = ((FY-FZ)*(FY+FZ) + (SGNFAC-FAC2)*(SGNFAC+FAC2))/
     &           (2.0D0*FAC2*FY)
            SGNFAC = SQRT(FAC*FAC + 1.0D0)
            FAC = ((FX-FZ)*(FX+FZ)+FAC2*((FY/(FAC+SIGN(SGNFAC,FAC)))
     &           -FAC2))/FX
            FC  = 1.0D0
            SUM = 1.0D0
            DO 120 L = K, JJ
               I = L + 1
               SGNFAC = VAL(I)
               FY = WMAT(I,I)
               FAC2   = SGNFAC*SUM
               SGNFAC = SGNFAC*FC
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               VAL(L) = FZ
               FC  = FAC/FZ
               SUM = FAC2/FZ
               FAC = (FX*FC) + (SGNFAC*SUM)
               SGNFAC = -(FX*SUM) + (SGNFAC*FC)
               FAC2 = FY*SUM
               FY = FY*FC
               DO 122 M = 1, IDIM2
                  FX = VMAT(M,L)
                  FZ = VMAT(M,I)
                  VMAT(M,L) = (FX*FC) + (FZ*SUM)
                  VMAT(M,I) =-(FX*SUM) + (FZ*FC) 
 122           CONTINUE
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               WMAT(L,L) = FZ
               IF (FZ .NE. D0) THEN
                  FZ = 1.0D0/FZ
                  FC = FAC*FZ
                  SUM = FAC2*FZ
               END IF
               FAC = (FC*SGNFAC) + (SUM*FY)
               FX = -(SUM*SGNFAC) + (FC*FY)
               DO 124 M = 1, IDIM1
                  FY = UMAT(M,L)
                  FZ = UMAT(M,I)
                  UMAT(M,L) = (FY*FC) + (FZ*SUM)
                  UMAT(M,I) = -(FY*SUM) + (FZ*FC) 
 124           CONTINUE
 120        CONTINUE
            VAL(K) = D0
            VAL(J) = FAC
            WMAT(J,J) = FX
 105     CONTINUE
 903     CONTINUE
 100  CONTINUE
C
C     Some output
C
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('U matrix in GTBINV',-1) 
         CALL OUTPUT(UMAT,1,IDIM1,1,IDIM2,MXCOOR,2*MXCOOR,1,LUPRI)
         CALL HEADER('W matrix in GTBINV',-1) 
         CALL OUTPUT(WMAT,1,IDIM2,1,IDIM2,2*MXCOOR,2*MXCOOR,1,LUPRI)
         CALL HEADER('V matrix in GTBINV',-1) 
         CALL OUTPUT(VMAT,1,IDIM2,1,IDIM2,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
C
C     The inverse is constructed:
C
C           (B^t)^-1 = V (1/w_ii) U^t
C
C     All non-negative values of W are inverted, values close to zero
C     are set to zero.
C
      NPROJ = 0
      DO 130 I = 1, IDIM2
         IF (ABS(WMAT(I,I)) .LE. 1.0D-10) THEN
            WMAT(I,I) = D0
         ELSE
            WMAT(I,I) = 1.0D0/WMAT(I,I)
            NPROJ = NPROJ + 1
         END IF
 130  CONTINUE
C
C     A projection operator that will be used later,
C     is also constructed: P = V W W^-1 V^t
C
      CALL DZERO(PJINMT,4*MXCOOR*MXCOOR)
      DO 140 I = 1, IDIM2
         DO 142 J = 1, IDIM2
            DO 144 K = 1, IDIM2
               IF (ABS(WMAT(K,K)) .GE. 1.0D-10)
     &              PJINMT(I,J) = PJINMT(I,J) + VMAT(I,K)*VMAT(J,K)
 144        CONTINUE
 142     CONTINUE
 140  CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projection matrix',-1) 
         CALL OUTPUT(PJINMT,1,IDIM2,1,IDIM2,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
C
C     We construct (1/w_ii) U^t
C      
      DO 160 I = 1, IDIM2
         DO 162 J = 1, IDIM1
            DO 164 K = 1, IDIM2
               BMTINV(I,J) = BMTINV(I,J) + WMAT(I,K)*UMAT(J,K)
 164        CONTINUE
 162     CONTINUE
 160  CONTINUE
C
C     Then V (1/w_ii) U^t
C
      CALL DZERO(WMAT,4*MXCOOR*MXCOOR)
      DO 170 I = 1, IDIM2
         DO 172 J = 1, IDIM1
            DO 174 K = 1, IDIM2
               WMAT(I,J) = WMAT(I,J) + VMAT(I,K)*BMTINV(K,J)
 174        CONTINUE
 172     CONTINUE
 170  CONTINUE
      CALL DZERO(BMTINV,2*MXCOOR*MXCOOR)
      DO 180 J = 1, IDIM1
         DO 182 I = 1, IDIM2
            BMTINV(I,J) = WMAT(I,J)
 182     CONTINUE
 180  CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Inverse of B^t',-1) 
         CALL OUTPUT(BMTINV,1,IDIM2,1,IDIM1,2*MXCOOR,MXCOOR,1,LUPRI)
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Total number of internal coordinates: ',IINTCR
         WRITE(LUPRI,*) 'Number of non-redundant coordinates : ',NPROJ
         WRITE(LUPRI,*)
      END IF
C
C     Put number of coordinates to be projected away in NPROJ
C
      NPROJ = IINTCR - NPROJ
      RETURN
      END

C  /* Deck gq2gx */
      SUBROUTINE GQ2GX(GRADQ,GRADX,WILBMT)
C
C     Transforms the gradient in redundant internal coordinates to
C     Cartesian coordinates:
C                              g_x = B^t g_q
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
      DIMENSION GRADQ(2*MXCOOR), GRADX(MXCOOR), WILBMT(2*MXCOOR,MXCOOR)
      CALL DZERO(GRADX,MXCOOR)
      DO 10 I = 1, ICRTCR
         DO 12 J = 1, IINTCR
            GRADX(I) = GRADX(I) + WILBMT(J,I)*GRADQ(J)
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck gx2gq */
      SUBROUTINE GX2GQ(GRADX,GRADQ,BMTINV)
C
C     Transforms the gradient in Cartesian coordinates to redundant
C     internal coordinates:
C                              g_q = (B^t)^-1 g_x
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
      DIMENSION GRADX(MXCOOR), GRADQ(2*MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      CALL DZERO(GRADQ,2*MXCOOR)
      DO 10 I = 1, IINTCR
         DO 12 J = 1, ICRTCR
            GRADQ(I) = GRADQ(I) + BMTINV(I,J)*GRADX(J)
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck hq2hx */
      SUBROUTINE HQ2HX(ATMARR,TMPMAT,TMPMT2,HESSQ,GRADQ,HESSX,WILBMT)
C
C     Transforms the Hessian in redundant internal coordinates to
C     Cartesian coordinates:
C                              H_x = B^t H_q B + B'^t g_q
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
      DIMENSION ATMARR(MXCENT,8), TMPMAT(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMT2(2*MXCOOR,2*MXCOOR), HESSQ(2*MXCOOR,2*MXCOOR)
      DIMENSION GRADQ(2*MXCOOR), HESSX(MXCOOR,MXCOOR)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR)
      CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
      DO 10 IC = 1, IINTCR
         CALL GTDWLM(IC,HESSX,ATMARR,WILBMT)
         DO 12 I = 1, ICRTCR
            DO 14 J = 1, ICRTCR
               TMPMAT(I,J) = TMPMAT(I,J) + HESSX(I,J)*GRADQ(IC)
 14         CONTINUE
 12      CONTINUE
 10   CONTINUE
      CALL DZERO(TMPMT2,4*MXCOOR*MXCOOR)
      DO 20 I = 1, IINTCR
         DO 22 J = 1, ICRTCR
            DO 24 K = 1, IINTCR
               TMPMT2(I,J) = TMPMT2(I,J) + HESSQ(I,K)*WILBMT(K,J)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      CALL DZERO(HESSX,MXCOOR*MXCOOR)
      DO 30 I = 1, ICRTCR
         DO 32 J = 1, ICRTCR
            DO 34 K = 1, IINTCR
               HESSX(I,J) = HESSX(I,J) + WILBMT(K,I)*TMPMT2(K,J)
 34         CONTINUE
            HESSX(I,J) = HESSX(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
C
C     We make sure the Cartesian Hessian is symmetric.
C
      DO 40 I = 1, ICRTCR
         DO 42 J = 1, I
            HESSX(J,I) = HESSX(I,J)
 42      CONTINUE
 40   CONTINUE
      RETURN
      END

C  /* Deck hx2hq */
      SUBROUTINE HX2HQ(ATMARR,TMPMAT,TMPMT2,TMPMT3,HESSX,GRADQ,HESSQ,
     &                 WILBMT,BMTINV)
C
C     Transforms the Hessian in Cartesian coordinates to redundant
C     internal coordinates:
C                         H_q = (B^t)^-1 (H_x - B'^t g_q) B^-1
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      DIMENSION TMPMT2(2*MXCOOR,2*MXCOOR), TMPMT3(2*MXCOOR,2*MXCOOR)
      DIMENSION HESSX(MXCOOR,MXCOOR), GRADQ(2*MXCOOR)
      DIMENSION HESSQ(2*MXCOOR,2*MXCOOR)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      CALL DZERO(HESSQ,4*MXCOOR*MXCOOR)
      CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
      CALL DZERO(TMPMT2,4*MXCOOR*MXCOOR)
      CALL DZERO(TMPMT3,4*MXCOOR*MXCOOR)
      DO 10 IC = 1, IINTCR
         CALL GTDWLM(IC,TMPMAT,ATMARR,WILBMT)
         DO 12 I = 1, ICRTCR
            DO 14 J = 1, ICRTCR
               TMPMT2(I,J) = TMPMT2(I,J) + TMPMAT(I,J)*GRADQ(IC)
 14         CONTINUE
 12      CONTINUE
 10   CONTINUE
      DO 16 J = 1, ICRTCR
         DO 18 I = 1, ICRTCR
            TMPMT2(I,J) = HESSX(I,J) - TMPMT2(I,J)
 18      CONTINUE
 16   CONTINUE
      DO 20 I = 1, ICRTCR
         DO 22 J = 1, IINTCR
            DO 24 K = 1, ICRTCR
               TMPMT3(I,J) = TMPMT3(I,J) + TMPMT2(I,K)*BMTINV(J,K)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      DO 30 I = 1, IINTCR
         DO 32 J = 1, IINTCR
            DO 34 K = 1, ICRTCR
               HESSQ(I,J) = HESSQ(I,J) + BMTINV(I,K)*TMPMT3(K,J)
 34         CONTINUE
 32      CONTINUE
 30   CONTINUE
      RETURN
      END

C  /* Deck cghint */
      SUBROUTINE CGHINT(SCLVEC,HESTMP,TMPMAT,TMPMT2,TMPMT3,TMPMT4,
     &                  WILBMT,BMTINV,HESINT)
C
C     Transforms Cartesian gradient and Hessian to redundant internal
C     coordinates.
C
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <optinf.h>
#include <taymol.h>
#include <nuclei.h>
#include <symmet.h>
C     #include <tranuc.h>
      DIMENSION SCLVEC(MXCOOR), HESTMP(MXCOOR,MXCOOR)
      DIMENSION TMPMAT(2*MXCOOR,2*MXCOOR), TMPMT2(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMT3(2*MXCOOR,2*MXCOOR), TMPMT4(2*MXCOOR,2*MXCOOR)
      DIMENSION WILBMT(2*MXCOOR,MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      DIMENSION HESINT(2*MXCOOR,2*MXCOOR)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      CALL DZERO(SCLVEC,MXCOOR)
      DO 1 IREP = 0, MAXREP
         DO 2 ICENT = 1, NUCIND
#include <scalar.h>
            DO 3 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from CGHINT','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Scaling vector',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (MAXREP .GT. 0) THEN
         CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
         CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,NCRREP(0,1),3*NUCDEP)
         CALL GX2GQ(TMPMAT,GRDINT,BMTINV)
         GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(TMPMAT,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         CALL DZERO(HESTMP,MXCOOR*MXCOOR)
         CALL TRAHES(HESMOL,MXCOOR,HESTMP,TMPMT2,TMPMT3,MXCOOR,
     &               3*NUCDEP,1)
         CALL HX2HQ(TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESTMP,GRDINT,HESINT,
     &              WILBMT,BMTINV)
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian Hessian',-1)
            CALL OUTPUT(HESTMP,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,1,LUPRI)
         END IF
      ELSE
         CALL GX2GQ(GRDMOL,GRDINT,BMTINV)
         GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(GRDMOL,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         CALL HX2HQ(TMPMAT,TMPMT2,TMPMT3,TMPMT4,HESMOL,GRDINT,HESINT,
     &              WILBMT,BMTINV)
         IF (IPRINT .GE. IPRMED) THEN
            CALL HEADER('Cartesian Hessian',-1)
            CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Internal gradient',-1)
         CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Internal Hessian',-1)
         CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck prjint */
      SUBROUTINE PRJINT(TMPMAT,PJINMT,HESINT)
C
C     Because of possible redundancies, we have to project
C     both the gradient and Hessian.
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
      DIMENSION TMPMAT(2*MXCOOR,2*MXCOOR), PJINMT(2*MXCOOR,2*MXCOOR),
     &          HESINT(2*MXCOOR,2*MXCOOR)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from PRJINT','*',103)
         CALL HEADER('Projection matrix',-1) 
         CALL OUTPUT(PJINMT,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Unprojected gradient',-1) 
         CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Unprojected Hessian',-1) 
         CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
      CALL FLSHFO(LUPRI)
C
C     First we project the gradient
C
      CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
      DO 10 I = 1, IINTCR
         DO 12 J = 1, IINTCR
            TMPMAT(I,1) = TMPMAT(I,1) + PJINMT(I,J)*GRDINT(J)
 12      CONTINUE
 10   CONTINUE
      DO 15 I = 1, IINTCR
         GRDINT(I) = TMPMAT(I,1)
 15   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projected gradient',-1) 
         CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
C
C     Then the Hessian
C
      CALL DZERO(TMPMAT,4*MXCOOR*MXCOOR)
      DO 20 I = 1, IINTCR
         DO 22 J = 1, IINTCR
            DO 24 K = 1, IINTCR
               TMPMAT(I,J) = TMPMAT(I,J) + PJINMT(I,K)*HESINT(K,J)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      CALL DZERO(HESINT,4*MXCOOR*MXCOOR)
      DO 25 I = 1, IINTCR
         DO 27 J = 1, IINTCR
            DO 29 K = 1, IINTCR
               HESINT(I,J) = HESINT(I,J) + TMPMAT(I,K)*PJINMT(K,J)
 29         CONTINUE
 27      CONTINUE
 25   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Projected Hessian',-1) 
         CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck diaint */
      SUBROUTINE DIAINT(EVEC,EVCTMP,HESPCK,TMPHES,TMPMAT,THRIND,HESINT)
C
C     The Hessian in redundant internal coordinates is diagonalized.
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
      DIMENSION EVEC(2*MXCOOR,2*MXCOOR), EVCTMP(IINTCR,IINTCR)
      DIMENSION HESPCK(IINTCR*IINTCR), TMPHES(IINTCR,IINTCR)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), HESINT(2*MXCOOR,2*MXCOOR)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
      DO 10 J = 1, IINTCR
         DO 12 I = 1, IINTCR
            TMPHES(I,J) = HESINT(I,J)
 12      CONTINUE
 10   CONTINUE
      CALL DZERO(HESPCK,IINTCR*IINTCR)
      CALL DSITSP(IINTCR,TMPHES,HESPCK)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from DIAINT','*',103)
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Internal gradient',-1) 
            CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
            CALL HEADER('Internal Hessian',-1) 
            CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,
     &                                    2*MXCOOR,2*MXCOOR,1,LUPRI)
            CALL HEADER('Packed Hessian',-1) 
            CALL OUTPUT(HESPCK,1,1,1,(IINTCR*(IINTCR+1))/2,
     &                                      1,IINTCR*IINTCR,1,LUPRI)
         END IF
      END IF
      CALL DUNIT(EVCTMP,IINTCR)
      CALL JACO(HESPCK,EVCTMP,IINTCR,IINTCR,IINTCR,
     &                          TMPMAT(1),TMPMAT(2*MXCOOR*MXCOOR+1))
      INDTOT = 0
      DO 20 J = 1, IINTCR
         EVAL(J) = HESPCK(J*(J+1)/2)
         GRDDIA(J) = DDOT(IINTCR,GRDINT,1,EVCTMP(1,J),1)
         DO 22 I = 1, IINTCR
            EVEC(I,J) = EVCTMP(I,J)
 22      CONTINUE
         IF (EVAL(J) .LT. -THRIND) INDTOT = INDTOT + 1
 20   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Eigenvalues in DIAINT',-1)
         CALL OUTPUT(EVAL,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Eigenvectors in DIAINT',-1)
         CALL OUTPUT(EVEC,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
         CALL HEADER('Gradient (diagonal rep.) in DIAINT',-1)
         CALL OUTPUT(GRDDIA,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
C
C     The eigenvalues are sorted
C
      DO 30 I = 1, IINTCR
         JMIN = I
         EMIN = EVAL(I)
         DO 35 J = (I + 1), IINTCR
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
 35      CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
            CALL DSWAP(2*MXCOOR,EVEC(1,I),1,EVEC(1,JMIN),1)
            CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
 30   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Sorted eigenvalues in WLKEI1',-1)
         CALL OUTPUT(EVAL,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Sorted eigenvectors in WLKEI1',-1)
         CALL OUTPUT(EVEC,1,IINTCR,1,IINTCR,2*MXCOOR,2*MXCOOR,1,LUPRI)
         CALL HEADER('Gradient (sorted) in WLKEI1',-1)
         CALL OUTPUT(GRDDIA,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck fnstin */
      SUBROUTINE FNSTIN(SCLVEC,EVEC,TMPMAT,TMPMT2,TMPMT3,CSTEP,BMTINV)
C
C     We determine the step in redundant internal coordinates,
C     then we find the corresponding step vector in
C     Cartesian coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <optinf.h>
#include <pi.h>
      DIMENSION SCLVEC(MXCOOR), EVEC(2*MXCOOR,2*MXCOOR)
      DIMENSION TMPMAT(4*MXCOOR*MXCOOR), TMPMT2(4*MXCOOR*MXCOOR),
     &          TMPMT3(2*MXCOOR*MXCOOR)
      DIMENSION CSTEP(MXCOOR), BMTINV(*)
      EXTERNAL WLKSTL
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     First we find the internal step in diagonal representation.
C
      CALL DZERO(STPDIA,2*MXCOOR)
      DO 10 I = 1, IINTCR
         IF (EVAL(I) .LE. 1.0D-8) THEN
            STPDIA(I) = D0
         ELSE
            STPDIA(I) = -GRDDIA(I)/EVAL(I)
         END IF
 10   CONTINUE
      STPNRM = SQRT(DDOT(IINTCR,STPDIA,1,STPDIA,1))
C
C     If Newton step is larger that trust radius, we take a step
C     to the boundary.
C
      IF ((STPNRM .GT. TRSTRA) .AND. (.NOT. NOTRST)) THEN
         IF (IPRINT .GE. IPRMED) THEN
            WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', STPNRM
            WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', TRSTRA
         END IF
C
C     Simple rip-off of level-shift from ABAWLK.
C
         GRDMIN = ABS(GRDDIA(1))
         HESMIN = EVAL(1)
         XMIN = GRADNM/TRSTRA - MIN(HESMIN,D0)
         XMAX = MAX(D0, -HESMIN + DP5*GRDMIN/TRSTRA)
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,'(/A,2(1P,D12.5,2X))')
     &          ' XMIN AND XMAX BEFORE WLKBIS (WLKSTL):  ',XMIN,XMAX
         END IF
         CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA,EVAL,TRSTRA,IINTCR,
     &                                                  WLKSTL,IFAIL)
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,'(A,2(1P,D12.5,2X))')
     &          ' XMIN and XMAX after WLKBIS (WLKSTL):   ',XMIN,XMAX
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            WRITE(LUPRI,'(A,1P,D12.5)')' Level shift parameter: ',RNU
         END IF
         IF (IFAIL .EQ. 0) THEN
            WRITE(LUPRI,'(/A/,2(A,F10.6))')
     &        '  *** Wrong interval specified in WLKBIS (WLKSTL) ***',
     &        ' XMAX = ', XMAX, ' XMIN = ', XMIN
         ELSE IF (IFAIL .EQ. 1) THEN
            WRITE(LUPRI,'(3(/A))')
     &      '  *** WARNING WLKBIS (WLKSTL) ***',
     &      ' Desired accuracy not obtained in the specified maximum',
     &      ' number of iterations.'
         END IF
         DO 15 I = 1, IINTCR
            STPDIA(I) = -GRDDIA(I)/(EVAL(I) + RNU)
 15      CONTINUE
         STPNRM = SQRT(DDOT(IINTCR,STPDIA,1,STPDIA,1))
         IF (IPRINT .GE. IPRMED) THEN
            WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', STPNRM
         END IF
      END IF
      ERGPRD = DDOT(IINTCR,GRDDIA,1,STPDIA,1)
     &         + 0.5D0*DV3DOT(IINTCR,STPDIA,EVAL,STPDIA)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER('Output from FNSTIN','*',103)
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Internal diagonal gradient',-1) 
            CALL OUTPUT(GRDDIA,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
            CALL HEADER('Eigenvalues',-1) 
            CALL OUTPUT(EVAL,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal diagonal step',-1) 
            CALL OUTPUT(STPDIA,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         END IF
         WRITE (LUPRI,'(/A,F25.15)') ' Predicted energy change', ERGPRD
      END IF
C
C     The diagonal step is converted to ordinary internal
C     coordinates.
C
      CALL DZERO(STPINT,2*MXCOOR)
      DO 20 I = 1, IINTCR
         DO 22 J = 1, IINTCR
            STPINT(I) = STPINT(I) + STPDIA(J)*EVEC(I,J)
 22      CONTINUE
C
C     For angles and dihedral angles we have to avoid step components
C     giving multiples of 2*pi.
C
         IF (INTCRD(I,1) .GT. 10) STPINT(I) = MOD(STPINT(I),2*PI)
         IF (ABS(STPINT(I)) .LE. 1.0D-8) STPINT(I) = D0
 20   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Step in internal coordinates',-1) 
         CALL OUTPUT(STPINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
      STPNRM = SQRT(DDOT(IINTCR,STPINT,1,STPINT,1))
C
C     Find Cartesian step vector.
C
      CALL STPI2C(SCLVEC,TMPMAT(1),TMPMAT(8*MXCENT+1),TMPMT2(1),
     &     TMPMT2(2*MXCOOR+1),TMPMT2(4*MXCOOR+1),
     &     TMPMT2(6*MXCOOR+1),CSTEP,BMTINV,TMPMT3(1),
     &     TMPMT3(1+MXCOOR*MXCOOR))
      RETURN
      END
Cjth - SUB.. GETINT renamed to GETIN2 since GETINT exists in dirac/dirone.F
C  /* Deck getint */
      SUBROUTINE GETIN2(IATOM,ATMCRD,VALINT)
C     
C     Determines the value of all redundant internal coordinates.
C
#include <implicit.h>
#include <mxcent.h>
#include <optinf.h>
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMCRD(MXCENT,8), VALINT(2*MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3)
      PARAMETER (D0 = 0.0D0, DEG175 = 175D0*PI/180D0)
      CALL DZERO(VALINT,2*MXCOOR)
      DO 10 IC = 1, IINTCR
         IF (INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VALINT(IC) = SQRT(DDOT(3,VEC1,1,VEC1,1))
         ELSE IF (INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(INTCRD(IC,4),2) - ATMCRD(INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(INTCRD(IC,4),3) - ATMCRD(INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(INTCRD(IC,4),4) - ATMCRD(INTCRD(IC,3),4)
C     
C     Regular angles
C     
            IF (INTCRD(IC,1) .EQ. 11) THEN
               CALL VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0D-8) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               CALL NRMLVC(VEC3)
C     
C     Second coordinate of angles larger than 175 degrees.
C     
            ELSE
               CALL VECPRD(VEC1,VEC3,VEC4)
               CALL NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
            CALL VECPRD(VEC3,VEC1,VEC4)
            CALL NRMLVC(VEC4)
            VALINT(IC) = VECANG(VEC1,VEC2)
            IF ((INTCRD(IC,1) .EQ. 12) .OR. (INTCRD(IC+1,1) .EQ. 12))
     &           VALINT(IC) = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC2)
C     
C     Dihedral angles
C     
         ELSE
            VEC1(1) = ATMCRD(INTCRD(IC,2),2) - ATMCRD(INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(INTCRD(IC,2),3) - ATMCRD(INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(INTCRD(IC,2),4) - ATMCRD(INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(INTCRD(IC,4),2) - ATMCRD(INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(INTCRD(IC,4),3) - ATMCRD(INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(INTCRD(IC,4),4) - ATMCRD(INTCRD(IC,3),4)
            VEC3(1) = ATMCRD(INTCRD(IC,5),2) - ATMCRD(INTCRD(IC,4),2)
            VEC3(2) = ATMCRD(INTCRD(IC,5),3) - ATMCRD(INTCRD(IC,4),3)
            VEC3(3) = ATMCRD(INTCRD(IC,5),4) - ATMCRD(INTCRD(IC,4),4)
            CALL NRMLVC(VEC2)
            CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
            VEC1(1) = VEC1(1) - CMPNT1*VEC2(1) 
            VEC1(2) = VEC1(2) - CMPNT1*VEC2(2) 
            VEC1(3) = VEC1(3) - CMPNT1*VEC2(3) 
            CALL NRMLVC(VEC1)
            CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
            VEC3(1) = VEC3(1) - CMPNT2*VEC2(1) 
            VEC3(2) = VEC3(2) - CMPNT2*VEC2(2) 
            VEC3(3) = VEC3(3) - CMPNT2*VEC2(3) 
            CALL NRMLVC(VEC3)
            IF ((DDOT(3,VEC1,1,VEC1,1) .GT. 1.0D-16) .AND.
     &           (DDOT(3,VEC3,1,VEC3,1) .GT. 1.0D-16)) THEN
               VALINT(IC) = VECANG(VEC1,VEC3)
               IF (ABS(VALINT(IC)) .GT. DEG175) THEN
                  CALL VECPRD(VEC1,VEC2,VEC4)
                  CALL NRMLVC(VEC4)
                  VALINT(IC) = VECANG(VEC1,VEC4) + VECANG(VEC4,VEC3)
               END IF
               CALL VECPRD(VEC1,VEC3,VEC4)
               IF (DDOT(3,VEC2,1,VEC4,1) .LT. 0.0D0)
     &              VALINT(IC) = -VALINT(IC)
            ELSE
               VALINT(IC) = D0
            END IF
         END IF
 10   CONTINUE
      RETURN
      END

C  /* Deck stpi2c */
      SUBROUTINE STPI2C(SCLVEC,CRDOLD,CRDNEW,VALOLD,VALNEW,
     &                  DIFF,TMPVEC,CSTEP,BMTINV,CSTRA,SCTRA)
C     
C     Transforms internal step to Cartesian step iteratively.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <optinf.h>
C     #include <tranuc.h>
#include <pi.h>
      DIMENSION SCLVEC(MXCOOR), CRDOLD(MXCENT,8), CRDNEW(MXCENT,8)
      DIMENSION VALOLD(2*MXCOOR), VALNEW(2*MXCOOR), DIFF(2*MXCOOR)
      DIMENSION TMPVEC(MXCOOR), CSTEP(MXCOOR), BMTINV(2*MXCOOR,MXCOOR)
      DIMENSION CSTRA(MXCOOR*MXCOOR), SCTRA(MXCOOR*MXCOOR)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (ITRLIM = 25)
      CHARACTER STPTXT*16
C
C     We have to find the value of all internal coordinates for the
C     old geometry.
C
      CALL ATMINI(CRDOLD,IATOM,.TRUE.)
      CALL GETIN2(IATOM,CRDOLD,VALOLD)
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Previous coordinates',-1) 
         CALL OUTPUT(CRDOLD,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         CALL HEADER('Previous internal values',-1) 
         CALL OUTPUT(VALOLD,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
C
C     First estimate of Cartesian step.
C
      CALL DZERO(TMPVEC,MXCOOR)
      DO 10 I = 1, ICRTCR
         DO 12 J = 1, IINTCR
            TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*STPINT(J)
 12      CONTINUE
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Inverse of B^t',-1) 
         CALL OUTPUT(BMTINV,1,IINTCR,1,ICRTCR,2*MXCOOR,MXCOOR,1,LUPRI)
         CALL HEADER('First estimate of Cartesian step',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
      END IF
      ITRCRD = 1
      RMSOLD = D0
      RMS1ST = -D1
 123  CONTINUE
      DO 15 I = 1, IATOM
         CRDNEW(I,1) = CRDOLD(I,1)
         DO 17 J = 1, 3
            CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 17      CONTINUE
 15   CONTINUE
      CALL GETIN2(IATOM,CRDNEW,VALNEW)
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('New coordinates',-1) 
         CALL OUTPUT(CRDNEW,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         CALL HEADER('New internal values',-1) 
         CALL OUTPUT(VALNEW,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
      END IF
      CALL DZERO(DIFF,2*MXCOOR)
      RMSINT = D0
      DO 19 I = 1, IINTCR
         DIFF(I) = STPINT(I) - (VALNEW(I) - VALOLD(I))
         IF (INTCRD(I,1) .GT. 10) DIFF(I) = MOD(DIFF(I),2*PI)
         IF (ABS(STPINT(I)) .LT. 1.0D-14) DIFF(I) = D0
         RMSINT = RMSINT + DIFF(I)*DIFF(I)
 19   CONTINUE
      RMSINT = SQRT(RMSINT/(1.0D0*IINTCR))
      IF (RMS1ST .LT. D0) RMS1ST = RMSINT 
      DIFRMS = ABS(RMSOLD-RMSINT)
      RMSOLD = RMSINT
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Difference to wanted step',-1) 
         CALL OUTPUT(DIFF,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         CALL HEADER('Root-mean-square of difference',-1) 
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',RMSINT
         CALL HEADER('Change in root-mean-square of difference',-1) 
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',DIFRMS
      END IF
      IF ((RMSINT .GE. 1.0D-6) .AND. (RMSINT .LE. 1.0D2)
     &     .AND. (DIFRMS .GE. 1.0D-8) .AND. (ITRCRD .LE. ITRLIM)) THEN
         DO 20 I = 1, ICRTCR
            DO 22 J = 1, IINTCR
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*DIFF(J)
 22         CONTINUE
 20      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            STPTXT = 'Updated step #XX'
            WRITE(STPTXT(15:16),'(I2)') ITRCRD
            CALL HEADER(STPTXT,-1) 
            CALL OUTPUT(TMPVEC,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
         END IF
         ITRCRD = ITRCRD + 1
         GOTO 123
      ELSE IF (((ITRCRD .GE. ITRLIM) .AND. (RMS1ST .LE. RMSINT)) .OR.
     &        (RMSINT .GE. 1.0D2)) THEN
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
     &           'Step does not converge, reverting to first estimate.'
            WRITE(LUPRI,*)
         END IF
         CALL DZERO(TMPVEC,MXCOOR)
         DO 30 I = 1, ICRTCR
            DO 32 J = 1, IINTCR
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*STPINT(J)
 32         CONTINUE
 30      CONTINUE
         DO 40 I = 1, IATOM
            DO 42 J = 1, 3
               CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 42         CONTINUE
 40      CONTINUE
         CALL GETIN2(IATOM,CRDNEW,VALNEW)
         DO 45 I = 1, IINTCR
            STPINT(I) = VALNEW(I) - VALOLD(I)
            IF (INTCRD(I,1) .GT. 10) STPINT(I) = MOD(DIFF(I),2*PI)
 45      CONTINUE
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Adjusted step',-1) 
            CALL OUTPUT(STPINT,1,1,1,IINTCR,1,2*MXCOOR,1,LUPRI)
         END IF
      END IF
C
C     We symmetrize the step, and scale it appropriately.
C
      CALL DZERO(STPSYM,2*MXCOOR)
C
      CALL TRACOR(CSTRA,SCTRA,1,ICRTCR,0)
Cjth - bugfix      CALL MPAB(CSTRA ,ICRTCR,ICRTCR,ICRTCR,ICRTRC,
      CALL MPAB(CSTRA ,ICRTCR,ICRTCR,ICRTCR,ICRTCR,
     &          TMPVEC,ICRTCR,     1,MXCOOR,     1,
     &          STPSYM,              MXCOOR,     1)
      DO 90 I = 1, IINTCR
         STPSYM(I) = STPSYM(I)/SCLVEC(I)
 90   CONTINUE
      CALL WLKCOR(STPSYM,CSTEP,ICRTCR,MXCOOR,IPRINT)
      RETURN
      END

C  /* Deck nrmlvc */
      SUBROUTINE NRMLVC(VEC)
C     
C     Normalizes a three-dimensional vector.
C
#include <implicit.h>
      DIMENSION VEC(3)
      VECNRM = SQRT(DDOT(3,VEC,1,VEC,1))
      VEC(1) = VEC(1)/VECNRM
      VEC(2) = VEC(2)/VECNRM
      VEC(3) = VEC(3)/VECNRM
      RETURN
      END

C  /* Deck vecprd */
      SUBROUTINE VECPRD(VEC1,VEC2,VPRD)
C     
C     Finds the vector product of two three-dimensional vectors.
C
#include <implicit.h>
      DIMENSION VEC1(3),VEC2(3),VPRD(3)
      VPRD(1) =  VEC1(2)*VEC2(3) - VEC2(2)*VEC1(3)
      VPRD(2) = -VEC1(1)*VEC2(3) + VEC2(1)*VEC1(3)
      VPRD(3) =  VEC1(1)*VEC2(2) - VEC2(1)*VEC1(2)
      RETURN
      END
