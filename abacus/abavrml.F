C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : abavrml.F
C SOURCE  : $Source$
C
C
#ifdef UNDEF
/* Comdeck log */
960319-vebjornb: New module for visualization through VRML models.
#endif
C  /* Deck mkvrml */
      SUBROUTINE MKVRML(LAST,ATMARR,EVEC)
C     
C     Punch out geometry to VRML-file
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <cbiwlk.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL LAST
      DIMENSION ATMARR(MXCENT,8), EVEC(2*MXCOOR,2*MXCOOR)
      DIMENSION EVC1(MXCOOR)
      CHARACTER*12 FILENM
#include <ibtfun.h>
C
C     Initialize the ATMARR array. The first index runs over all
C     atoms, the second marks the following properties:
C
C               1 - Element number
C               2 - X coordinate of atom
C               3 - Y coordinate of atom
C               4 - Z coordinate of atom
C               5 - Covalent radius
C
      CALL ATMINI(ATMARR,IATOM,.FALSE.)
C
C     We proceed to open the output-file
C
      LUVRML = 66
      FILENM = 'first.wrl'
      IF (LAST) THEN
         LUVRML = 67
         FILENM = 'last.wrl'
      END IF
      CALL VRINI(LUVRML,FILENM)
C
C     We create all the atoms
C
      INDX = 1
      CALL DRWATM(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
C
C     Then we draw bonds between the atoms if this is requested.
C
      IF (VRBOND) CALL DRWBND(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
      CALL VREND(LUVRML)
C
C     Finally all eigenvectors are visualized
C
      IF (VREIGV) THEN
         LUVRML = 68
         IEIG = 1
         IEVEC = 1
         IF (IPRINT .GT. 0) THEN
            WRITE(LUPRI,*)
            CALL HEADER('VRML Visualization of Eigenvectors',-1)
            WRITE(LUPRI,'(A)')
     &        ' Eig.Vec.     Filename      Symmetry      Eig.Value '
            WRITE(LUPRI,'(A)')
     &        '----------------------------------------------------'
         END IF
C
C     We loop over all symmetries...
C
         DO 30 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               II = 0
               NCR = NCRREP(IREP,1)
               DO 35 I = 0, IREP - 1
                  II = II + NCRREP(I,1)
 35            CONTINUE
C
C     ... and all vectors in each symmetry
C
               DO 40 IVEC = 1, NCR
C
C     Only eigenvectors with a non-zero eigenvalue is visualized,
C     that is only the eigenvectors associated with.
C
                  IF (EVAL(IEVEC+IVEC-1) .GT. 1.0D-3) THEN
                     DO 50 I = 1, NCR
                        EVC1(I) = EVEC(II+I,II+IVEC)
 50                  CONTINUE
                     INDX = 1
                     FILENM = 'eigv_XXX.wrl'
                     WRITE(FILENM(6:8),'(I3)') IEIG
                     IF (IEIG .LT. 100) WRITE(FILENM(6:6),'(A1)') '0'
                     IF (IEIG .LT. 10) WRITE(FILENM(7:7),'(A1)') '0'
                     CALL VRINI(LUVRML,FILENM)
                     CALL DRWATM(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
                     IF (VRBOND)
     &                    CALL DRWBND(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
                     CALL DRWEIG(LUVRML,INDX,ATMARR,EVC1,IREP,NCR)
                     CALL VREND(LUVRML)
                     IF (IPRINT .GT. 0) THEN
                        WRITE(LUPRI,'(I5,A,A12,A,I1,A,F16.6)')
     &                       IEIG,'       ',FILENM,'       ',
     &                       IREP,'   ',EVAL(IEVEC+IVEC-1)
                     END IF
                     IEIG = IEIG + 1
                  END IF
 40            CONTINUE
            END IF
            IEVEC = IEVEC + NCRREP(IREP,1)
 30      CONTINUE
      END IF
      RETURN
      END

C  /* Deck vrini */
      SUBROUTINE VRINI(LUVRML,FILENM)
C     
C     Open unit LUVRML for output
C
#include <implicit.h>
      CHARACTER*12 FILENM
      CHARACTER*8 CHDATE, CHTIME
      CALL GETDAT(CHDATE,CHTIME)
#if defined (SYS_VAX)
      OPEN(LUVRML,STATUS='NEW',FORM='FORMATTED',
     &                                 CARRIAGECONTROL='LIST')
#endif
#if defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_CRAY) || defined (SYS_IRIX)
      OPEN(LUVRML,FILE=FILENM,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if defined (SYS_IBM)
      OPEN(LUVRML,FILE=FILENM,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
#if !defined (SYS_CRAY) && !defined (SYS_VAX) && !defined (SYS_ALLIANT) && !defined (SYS_CONVEX) && !defined (SYS_AIX) && !defined (SYS_PARAGON) && !defined (SYS_IBM) && !defined (SYS_IRIX)
      OPEN(LUVRML,STATUS='UNKNOWN',FORM='FORMATTED')
#endif
      WRITE(LUVRML,'(A)') '#VRML V1.0 ascii'
      WRITE(LUVRML,'(A)') 'Separator {'
      WRITE(LUVRML,'(A)') 'Info {'
      WRITE(LUVRML,'(A)') 'string "Created by DALTON'
      WRITE(LUVRML,'(A,A8,A8)') 'date: ',CHDATE,CHTIME
      WRITE(LUVRML,'(A)') '" }'
      WRITE(LUVRML,'(A)')
     &     'Complexity { fields [ SFFloat value ] value 0.5 }'
      WRITE(LUVRML,'(A)') 'Separator {'
      WRITE(LUVRML,'(A)') 'Separator {'
      WRITE(LUVRML,'(A)')
     &     'Complexity { fields [ SFFloat value ] value 0.5 }'
      RETURN
      END

C  /* Deck vrend */
      SUBROUTINE VREND(LUVRML)
C     
C     End of output, close unit LUVRML
C
#include <implicit.h>
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') '}'
      WRITE(LUVRML,'(A)') '}'
      CLOSE(LUVRML)
      RETURN
      END

C  /* Deck drwatm */
      SUBROUTINE DRWATM(LUVRML,INDX,IATOM,VIBRA,ATMARR)
C     
C     Draws all the atoms as spheres with appropriate colour
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL VIBRA
      DIMENSION ATMARR(MXCENT,8), RGBCOL(3)      
      DIMENSION COLOUR(60)
#include <ibtfun.h>
C
C     Array with colours for the first 20 elements
C
      DATA (COLOUR(I), I = 1, 60)/
C
C     H-He
C
     & 0.750D0, 0.750D0, 0.750D0, 0.950D0, 0.950D0, 0.600D0,
C
C     Li-Ne
C
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.400D0,
     & 0.400D0, 0.400D0, 0.300D0, 0.300D0, 0.300D0, 0.150D0, 0.150D0,
     & 1.000D0, 1.000D0, 0.100D0, 0.100D0, 0.100D0, 0.950D0, 0.100D0,
     & 0.700D0, 0.100D0, 0.100D0,
C
C     Na-Ar
C
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0,
     & 0.650D0, 0.650D0, 0.400D0, 0.400D0, 0.400D0, 0.760D0, 0.600D0,
     & 0.000D0, 0.960D0, 0.800D0, 0.200D0, 0.800D0, 0.960D0, 0.600D0,
     & 0.950D0, 0.600D0, 0.950D0,
C
C     K-Ca
C
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0/
C
C     We create all the atoms as spheres with appropriate colour
C
      LSTCHG = 0
      DO 10 INUC = 1, IATOM
         ICHARG = NINT(ATMARR(INUC,1))
         IF (ICHARG .NE. LSTCHG) THEN
            IF (INDX .LT. 10) THEN
               WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
            ELSE
               WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
            END IF
            IF (ICHARG .LE. 20) THEN
               RGBCOL(1) = COLOUR(3*(ICHARG-1)+1)
               RGBCOL(2) = COLOUR(3*(ICHARG-1)+2)
               RGBCOL(3) = COLOUR(3*(ICHARG-1)+3)
            ELSE
               RGBCOL(1) = 0.500D0
               RGBCOL(2) = 0.500D0
               RGBCOL(3) = 0.500D0
            END IF
            WRITE(LUVRML,'(A,3F7.3)') ' ambientColor',
     &           RGBCOL(1), RGBCOL(2), RGBCOL(3)
            WRITE(LUVRML,'(A,3F7.3)') ' diffuseColor',
     &           RGBCOL(1), RGBCOL(2), RGBCOL(3)
            WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
            WRITE(LUVRML,'(A)') ' shininess 0.750 }'
            INDX= INDX + 1
            LSTCHG = ICHARG
         END IF
         WRITE(LUVRML,'(A)') 'Separator {'
         WRITE(LUVRML,'(A,3F7.3,A)')
     &        ' Translation { translation ', ATMARR(INUC,2),
     &        ATMARR(INUC,3), ATMARR(INUC,4), ' }'
         RAD = 0.5D0*ATMARR(INUC,5)
         IF (VIBRA) RAD = 0.15D0
         WRITE(LUVRML,'(A,F5.3,A)') ' Sphere { radius ',RAD,' } }'
 10   CONTINUE
      RETURN
      END

C  /* Deck drwbnd */
      SUBROUTINE DRWBND(LUVRML,INDX,IATOM,VIBRA,ATMARR)
C
C     Draws bonds between the atoms
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL VIBRA
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), VEC3(3)
#include <ibtfun.h>
C
      RADIUS = 0.06D0
      IF (VIBRA) RADIUS = 0.015D0
C
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') 'Separator {'
      IF (INDX .LT. 10) THEN
         WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
      ELSE
         WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
      END IF
      WRITE(LUVRML,'(A)') ' ambientColor 0.600 0.600 0.600'
      WRITE(LUVRML,'(A)') ' diffuseColor 0.600 0.600 0.600'
      WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
      WRITE(LUVRML,'(A)') ' shininess 0.750 }'
      INDX = INDX + 1
      DO 10 I = 1, IATOM - 1
         DO 20 J = I + 1, IATOM
            RADI = ATMARR(I,5)
            RADJ = ATMARR(J,5)            
C     We find the bond vector ...
            VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
            VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
            VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
            DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
            IF (DIST .LE. 1.2D0*(RADI+RADJ)) THEN
C     ... the center of the bond ...
               VEC2(1) = 0.5D0*(ATMARR(I,2)+ATMARR(J,2))
               VEC2(2) = 0.5D0*(ATMARR(I,3)+ATMARR(J,3))
               VEC2(3) = 0.5D0*(ATMARR(I,4)+ATMARR(J,4))
C     ... the angle to rotate ...
               VEC3(1) = 0.0D0
               VEC3(2) = 1.0D0
               VEC3(3) = 0.0D0
               ANG = VECANG(VEC3,VEC1)
C     ... and finally a vector to rotate the bond around ... 
               VEC3(1) = VEC1(3)
               VEC3(2) = 0.0D0
               VEC3(3) = -VEC1(1)
C     For vibrational visualization, all radii are equal
               H = 0.25D0*(RADJ-RADI)
               IF (VIBRA) THEN
                  RADI = 0.300D0
                  RADJ = 0.300D0
                  H = 0.0D0
               END IF
C     We have to move the center of the bond in accordance with the radii
               VEC2(1) = VEC2(1) + H*(VEC1(1)/DIST)
               VEC2(2) = VEC2(2) + H*(VEC1(2)/DIST)
               VEC2(3) = VEC2(3) + H*(VEC1(3)/DIST)
C     We also have to shorten the bond, so that it only touches the spheres
               H = 1.05D0*(DIST - 0.5D0*(RADI + RADJ))
               IF (VIBRA) H = 0.97D0*H
C
               WRITE(LUVRML,'(A)') 'Separator {'
               WRITE(LUVRML,'(A,3F7.3,A)')
     &              ' Translation { translation ', 
     &              VEC2(1), VEC2(2), VEC2(3), ' }'
               WRITE(LUVRML,'(A,4F7.3,A)')
     &              ' Rotation { rotation ', 
     &              VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
               WRITE(LUVRML,'(A,F5.3,A,F5.3,A)')
     &              ' Cylinder { parts SIDES radius ',RADIUS,' height ',
     &              H,' } }'
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck drweig */
      SUBROUTINE DRWEIG(LUVRML,INDX,ATMARR,EVEC1,IREP,NCR)
C
C     Draws vectors to illustrate vibrational modes
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
#include <trkoor.h>
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION EVEC1(2*MXCOOR), EVEC2(2*MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3)
#include <ibtfun.h>
C
      CALL DZERO(EVEC2,2*MXCOOR)
      CALL DAXPY(NCR,1.0D0,EVEC1,1,EVEC2,1)
      CALL DZERO(EVEC1,2*MXCOOR)
      DO 10 I = 1, NCR
         EVEC1(I) = EVEC2(I)*1.0D0
 10   CONTINUE
C
C     The eigenvector in symmetry basis is transformed to
C     cartesian coordinates.
C
      CALL DZERO(EVEC2,2*MXCOOR)
      DO 20 IAT = 1, NUCIND
         DO 30 ICO = 1, 3
            ICCOOR = 3*(IAT - 1) + ICO
            ISCOOR = IPTCNT(ICCOOR,IREP,1)
            IF (ISCOOR .GT. 0) THEN
               EVEC2(ICCOOR)=EVEC1(ICCOOR)/SQRT(FMULT(ISTBNU(IAT)))
             END IF
 30       CONTINUE
 20   CONTINUE
C
C     We add necessary information to the ATMARR array:
C
C               6 - X component of eigenvector
C               7 - Y component of eigenvector
C               8 - Z component of eigenvector
C
      IATOM = 1
      DO 40 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 45 ISYMOP = 0, MAXOPR
C
C     This factor is a humble attempt to give the vectors a nice length
C
            FAC = 0.5D0*SQRT(1.0D0*NUCDEP)
C
C     Vectors on the symmetry dependent centres, should have a direction
C     according to the symmetry.
C
C     *****************************************************************
C     NOTE!!!!! This test is probably _NOT_ correct!
C     *****************************************************************
C
            IF (IBTAND(ISYMOP,IREP) .GT. 0) FAC = -FAC
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               ATMARR(IATOM,6) =
     &      FAC*PT(IBTAND(ISYMAX(1,1),ISYMOP))*EVEC2((ICENT-1)*3+1)
               ATMARR(IATOM,7) =
     &      FAC*PT(IBTAND(ISYMAX(2,1),ISYMOP))*EVEC2((ICENT-1)*3+2)
               ATMARR(IATOM,8) =
     &      FAC*PT(IBTAND(ISYMAX(3,1),ISYMOP))*EVEC2((ICENT-1)*3+3)
               IATOM = IATOM + 1
            END IF
 45      CONTINUE
 40   CONTINUE
      IATOM = IATOM - 1
C
C     We draw the white vectors...
C
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,.TRUE.)
C
C     ... then we turn all the vectors before we draw the black ones.
C
      DO 50 I = 1, IATOM
         ATMARR(I,6) = -ATMARR(I,6)
         ATMARR(I,7) = -ATMARR(I,7)
         ATMARR(I,8) = -ATMARR(I,8)
 50   CONTINUE
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,.FALSE.)
      RETURN
      END

C  /* Deck drwvec */
      SUBROUTINE DRWVEC(LUVRML,INDX,IATOM,ATMARR,WHITE)
C
C     Draws vectors
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <priunit.h>
      LOGICAL WHITE
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION VEC1(3), VEC2(3), VEC3(3)
C
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') 'Separator {'
      IF (INDX .LT. 10) THEN
         WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
      ELSE
         WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
      END IF
      IF (WHITE) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 1.000 1.000 1.000'
         WRITE(LUVRML,'(A)') ' diffuseColor 1.000 1.000 1.000'
      ELSE
         WRITE(LUVRML,'(A)') ' ambientColor 0.200 0.200 0.200'
         WRITE(LUVRML,'(A)') ' diffuseColor 0.200 0.200 0.200'
      END IF
      WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
      WRITE(LUVRML,'(A)') ' shininess 0.750 }'
      INDX = INDX + 1
      DO 10 I = 1, IATOM
C
C     Vec1 contains the position of the atom
C
         VEC1(1) = ATMARR(I,2)
         VEC1(2) = ATMARR(I,3)
         VEC1(3) = ATMARR(I,4)
C
C     Vec2 contains the eigenvector
C
         VEC2(1) = ATMARR(I,6)
         VEC2(2) = ATMARR(I,7)
         VEC2(3) = ATMARR(I,8)
         VECNRM = SQRT(DDOT(3,VEC2,1,VEC2,1))
C
C     The eigenvectors are only drawn if they have a certain length
C
         IF (VECNRM .GT. 0.1D0) THEN
            VEC3(1) = 0.0D0
            VEC3(2) = 1.0D0
            VEC3(3) = 0.0D0
C
C     We calculate a vector to rotate around, and the amount to rotate
C
            ANG = VECANG(VEC3,VEC2)
            VEC3(1) = VEC2(3)
            VEC3(2) = 0.0D0
            VEC3(3) = -VEC2(1)
C
C     The center of the vector is placed correctly
C
            VEC1(1) = VEC1(1)+0.5D0*VEC2(1)
            VEC1(2) = VEC1(2)+0.5D0*VEC2(2)
            VEC1(3) = VEC1(3)+0.5D0*VEC2(3)
            WRITE(LUVRML,'(A)') 'Separator {'
            WRITE(LUVRML,'(A,3F7.3,A)')
     &           ' Translation { translation ', 
     &           VEC1(1), VEC1(2), VEC1(3), ' }'
            WRITE(LUVRML,'(A,4F7.3,A)')
     &           ' Rotation { rotation ', 
     &           VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
            WRITE(LUVRML,'(A,F5.3,A,F5.3,A)')
     &           ' Cylinder { parts SIDES radius 0.03 height ',
     &           VECNRM,' } }'
            VEC1(1) = VEC1(1)+0.5D0*VEC2(1)
            VEC1(2) = VEC1(2)+0.5D0*VEC2(2)
            VEC1(3) = VEC1(3)+0.5D0*VEC2(3)
            WRITE(LUVRML,'(A)') 'Separator {'
            WRITE(LUVRML,'(A,3F7.3,A)')
     &           ' Translation { translation ', 
     &           VEC1(1), VEC1(2), VEC1(3), ' }'
            WRITE(LUVRML,'(A,4F7.3,A)')
     &           ' Rotation { rotation ', 
     &           VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
            WRITE(LUVRML,'(A)')
     &        ' Cone { parts ALL bottomRadius 0.060 height 0.100 } }'
         END IF
 10   CONTINUE
      RETURN
      END

C  /* Deck atmini */
      SUBROUTINE ATMINI(ATMARR,IATOM,BOHR)

C     Expand all symmetry-dependent atoms.
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <cbiwlk.h>
#include <priunit.h>
#include <codata.h>
#include <symmet.h>
      DIMENSION ATMARR(MXCENT,8)
      LOGICAL BOHR
#include <ibtfun.h>
C
C     We initialize the ATMARR array. The first index runs over all
C     atoms, the second marks the following properties:
C
C               1 - Element number
C               2 - X coordinate of atom
C               3 - Y coordinate of atom
C               4 - Z coordinate of atom
C               5 - Covalent radius
C
      FAC = XTANG
      IF (BOHR) FAC = 1.0D0
C
      IATOM = 1
      DO 10 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 20 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               ATMARR(IATOM,1) = 1.0D0*NINT(CHARGE(ICENT))
               ATMARR(IATOM,2) =
     &              FAC*PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               ATMARR(IATOM,3) =
     &              FAC*PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               ATMARR(IATOM,4) =
     &              FAC*PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               RAD = RADIUS(NINT(CHARGE(ICENT)))
               IF (RAD .LT. 0.0D0) RAD = 1.0D0
               ATMARR(IATOM,5) = RAD
               IATOM = IATOM + 1
            END IF
 20      CONTINUE
 10   CONTINUE
      IATOM = IATOM - 1
      RETURN
      END
