C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PLOT */
      SUBROUTINE PLOT(WORK,LWORK)
C=======================================================================
C     prepare {1,2,3}-dimensional plot
C     adaptation of Trond's plot routines
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "visual.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"

      LOGICAL   FILE_EXISTS, LBIT
      DIMENSION WORK(LWORK)

      CALL QENTER('PLOT')

#include "memint.h"

      CALL TITLER('{1,2,3}-dimensional plot','*',122)

      INQUIRE(FILE = 'DFCOEF',EXIST = FILE_EXISTS) !check if coefficients are on file
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'coefficients not found'
        CALL QUIT('ERROR in PLOT')
      ENDIF

!     allocate memory

      CALL MEMGET('REAL',KCMO     ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG     ,NORBT,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT_AO ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_MO ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING1 ,NESH(1),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING2 ,NESH(2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO     ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO1    ,3*NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO2    ,6*NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCNT    ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF     ,NTBAS(0)*4,
     &            WORK,KFREE,LFREE)

!     read coefficients and give info

      CALL DZERO(WORK(KCMO),N2BBASXQ)
      CALL DZERO(WORK(KEIG),NORBT)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PLOT')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),UNUSED,
     &            UNUSED,6,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')

      CALL FNDOCC(WORK(KEIG),WORK,KFREE,LFREE) !find occupation

!     the default here is to set the diagonal elements
!     of the MO density matrix to D0,
!     then to put D1 in the diagonal elements of the inactive block

!     if an orbital string is given the elements will be DZEROed again
!     and only selected elements will be set to D1
!     and later transformed to AO density matrix

!     this will not happen if a solution vector is read in
!     then MODIFIED_DMAT_AO will occupy the elements
!     and transform to AO basis

      CALL DFILL(WORK(KSTRING1),D0,NESH(1))
      CALL DFILL(WORK(KSTRING1),D1,NISH(1))
      IF(NFSYM .GT. 1) THEN
        CALL DFILL(WORK(KSTRING2),D0,NESH(2))
        CALL DFILL(WORK(KSTRING2),D1,NISH(2))
      ENDIF

      IF(LBIT(LBIT_PLOT_GENERAL,7)) THEN !use orbital string
        IF(ORBSTRING1(4:5) .EQ. '..') THEN
          READ(ORBSTRING1(1:3),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(6:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(3:4) .EQ. '..') THEN
          READ(ORBSTRING1(1:2),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(5:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(2:3) .EQ. '..') THEN
          READ(ORBSTRING1(1:1),'(I3)') ISTRING1_BEGIN
          READ(ORBSTRING1(4:8),'(I3)') ISTRING1_END
          GOTO 1
        ENDIF
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'wrong orbital string format'
        CALL QUIT('ERROR in PLOT')
 1      CONTINUE
        IF(
     &          (.NOT. (ISTRING1_BEGIN .LE. ISTRING1_END))
     &     .OR. (.NOT. (ISTRING1_BEGIN .GE. 0))
     &     .OR. (.NOT. (ISTRING1_END .LE. NESH(1)))) THEN
          WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
          WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
          CALL QUIT('ERROR in PLOT')
        ENDIF
        CALL DFILL(WORK(KSTRING1),D0,NESH(1))
        DO I = ISTRING1_BEGIN,ISTRING1_END
          IF(ISTRING1_BEGIN .NE. 0) WORK(KSTRING1-1+I) = D1
        ENDDO

        IF(NFSYM .GT. 1) THEN
          IF(ORBSTRING2(4:5) .EQ. '..') THEN
            READ(ORBSTRING2(1:3),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(6:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          IF(ORBSTRING2(3:4) .EQ. '..') THEN
            READ(ORBSTRING2(1:2),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(5:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          IF(ORBSTRING2(2:3) .EQ. '..') THEN
            READ(ORBSTRING2(1:1),'(I3)') ISTRING2_BEGIN
            READ(ORBSTRING2(4:8),'(I3)') ISTRING2_END
            GOTO 2
          ENDIF
          WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
          WRITE(LUPRI,'(A)') 'wrong orbital string format'
          CALL QUIT('ERROR in PLOT')
 2        CONTINUE
          IF(
     &            (.NOT. (ISTRING2_BEGIN .LE. ISTRING2_END))
     &       .OR. (.NOT. (ISTRING2_BEGIN .GE. 0))
     &       .OR. (.NOT. (ISTRING2_END .LE. NESH(2)))) THEN
            WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
            WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
            CALL QUIT('ERROR in PLOT')
          ENDIF
          CALL DFILL(WORK(KSTRING2),D0,NESH(2))
          DO I = ISTRING2_BEGIN,ISTRING2_END
            IF(ISTRING2_BEGIN .NE. 0) WORK(KSTRING2-1+I) = D1
          ENDDO
        ENDIF !IF(NFSYM .GT. 1) THEN
      ENDIF !IF(LBIT(LBIT_PLOT_GENERAL,7)) THEN !use orbital string

!     generate density matrix

      IF(LBIT(LBIT_PLOT_GENERAL,4)) THEN !read solution vector
        CALL MODIFIED_DMAT_AO(WORK(KCMO),
     &                        WORK(KDMAT_AO),WORK(KTMAT_AO),
     &                        WORK(KSTRING1),WORK(KSTRING2),
     &                        WORK(KFREE),LFREE)
      ELSE
        CALL OCCUPY_DMAT_MO(WORK(KDMAT_MO),
     &                      WORK(KSTRING1),WORK(KSTRING2))
        CALL DMAT_MO_2_AO(WORK(KCMO),
     &                    WORK(KDMAT_AO),
     &                    WORK(KDMAT_MO),
     &                    WORK(KFREE),LFREE)
      ENDIF
      CALL DSCAL(N2BBASXQ,D2,WORK(KDMAT_AO),1)

!     call PLOT routine

      CALL PLOT_SUB(WORK(KDMAT_AO),
     &              WORK(KGAO),
     &              WORK(KGAO1),
     &              WORK(KGAO2),
     &              WORK(KNCNT),
     &              WORK(KBUF),
     &              WORK(KFREE),LFREE)

!     release memory

      CALL MEMREL('PLOT',WORK,KCMO,KCMO,KFREE,LFREE)
      CALL QEXIT('PLOT')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PLOT_SUB */
      SUBROUTINE PLOT_SUB(DMAT_AO,
     &                    GAO,GAO1,GAO2,
     &                    NCNT,BUF,
     &                    WORK,LWORK)
C=======================================================================
C     do the {1,2,3}-dimensional plot
C     adaptation of Trond's plot routines
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "codata.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NR_QUANTITIES = 19)
      PARAMETER (MAX_NR_PAIRS = 200)

#include "dcbgen.h"

      CHARACTER QUANTITY_NAME(NR_QUANTITIES)*30,
     &          PAIR_NAME_A(MAX_NR_PAIRS)*3,
     &          PAIR_NAME_B(MAX_NR_PAIRS)*3

      CHARACTER*20 FORMAT_STRING
      CHARACTER*30 FILE_NAME

#include "visual.h"
#include "dcborb.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(LWORK),
     &          LUNIT(3,NR_QUANTITIES),
     &          QUANTITY(NR_QUANTITIES,0:3,NCUBE(3)),
     &          PAIR_COOR_A(MAX_NR_PAIRS,3),
     &          PAIR_COOR_B(MAX_NR_PAIRS,3),
     &          STEP_3D(3,0:3),
     &          ATOM_XYZ(4,NUCDEP)

#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"
#include "sdpre.h"
#include "ibtfun.h"

#include "memint.h"

      QUANTITY_NAME(1) = 'density'
      QUANTITY_NAME(2) = 'spin_density'
      QUANTITY_NAME(3) = 'current_density'
      QUANTITY_NAME(4) = 'ELF'
      QUANTITY_NAME(5) = 'r_times_j'

      QUANTITY_NAME(10) = 'opt0'
      QUANTITY_NAME(11) = 'opt1'
      QUANTITY_NAME(12) = 'opt2'
      QUANTITY_NAME(13) = 'opt3'
      QUANTITY_NAME(14) = 'opt4'
      QUANTITY_NAME(15) = 'opt5'
      QUANTITY_NAME(16) = 'opt6'
      QUANTITY_NAME(17) = 'opt7'
      QUANTITY_NAME(18) = 'opt8'
      QUANTITY_NAME(19) = 'opt9'

!     open files

      LUNIT_BASE = 13
      DO I = 1,3 !1D, 2D, 3D
        IF(LBIT(LBIT_PLOT_GENERAL,I)) THEN
          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
              WRITE(FILE_NAME,'(A5,I1,A2,A)') 'PLOT.', I, 'D.',
     &                                  QUANTITY_NAME(IQUANTITY)
              OPEN(LUNIT_BASE,
     &             FILE = FILE_NAME,
     &             STATUS = 'NEW',
     &             ACCESS = 'SEQUENTIAL',
     &             FORM = 'FORMATTED')
              LUNIT(I,IQUANTITY) = LUNIT_BASE
              LUNIT_BASE = LUNIT_BASE + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      IF(LBIT(LBIT_PLOT_GENERAL,1)) THEN !1D

        NR_PAIRS  = 0
        DO I = 1,NUCDEP
          DO J = I+1,NUCDEP
            NR_PAIRS = NR_PAIRS + 1

            WRITE(PAIR_NAME_A(NR_PAIRS),'(A2,I1)') NAMDEP(I), I
            PAIR_COOR_A(NR_PAIRS,1) = CORD(1,I)
            PAIR_COOR_A(NR_PAIRS,2) = CORD(2,I)
            PAIR_COOR_A(NR_PAIRS,3) = CORD(3,I)

            WRITE(PAIR_NAME_B(NR_PAIRS),'(A2,I1)') NAMDEP(J), J
            PAIR_COOR_B(NR_PAIRS,1) = CORD(1,J)
            PAIR_COOR_B(NR_PAIRS,2) = CORD(2,J)
            PAIR_COOR_B(NR_PAIRS,3) = CORD(3,J)
          ENDDO
        ENDDO

        DO IPAIR = 1,NR_PAIRS !loop over pairs

!         write pair label to file

          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
              WRITE(LUNIT(1,IQUANTITY),'(A5)') '*****'
              WRITE(LUNIT(1,IQUANTITY),'(A)')
              WRITE(LUNIT(1,IQUANTITY),'(A3,A1,A3)')
     &                       PAIR_NAME_A(IPAIR), '-', PAIR_NAME_B(IPAIR)
              WRITE(LUNIT(1,IQUANTITY),'(A)')
            ENDIF
          ENDDO

!         calculate distance and step

          DX = PAIR_COOR_B(IPAIR,1) - PAIR_COOR_A(IPAIR,1)
          DY = PAIR_COOR_B(IPAIR,2) - PAIR_COOR_A(IPAIR,2)
          DZ = PAIR_COOR_B(IPAIR,3) - PAIR_COOR_A(IPAIR,3)

          DISTANCE = DSQRT(DX*DX + DY*DY + DZ*DZ)

          DX = DX/DISTANCE
          DY = DY/DISTANCE
          DZ = DZ/DISTANCE

          ASTEP = STEP_1D/XTANG !-> angstrom
          NR_POINTS = INT(DISTANCE/ASTEP)
          IF(MOD(NR_POINTS,2) .EQ. 1) NR_POINTS = NR_POINTS + 1
          ASTEP = DISTANCE/SDPRE(NR_POINTS)

          IPOINT1 = -NR_POINTS/2
          IPOINT2 = NR_POINTS - IPOINT1
          DO IPOINT = IPOINT1,IPOINT2 !loop along line
            PP  = SDPRE(IPOINT)*ASTEP
            PX  = PAIR_COOR_A(IPAIR,1) + PP*DX
            PY  = PAIR_COOR_A(IPAIR,2) + PP*DY
            PZ  = PAIR_COOR_A(IPAIR,3) + PP*DZ

            CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                 GAO,GAO1,GAO2,
     &                                 PX,PY,PZ,
     &                                 NCNT,BUF,
     &                                 QUANTITY,
     &                                 WORK,KFREE,LFREE)

!           write buffer

            PP  = XTANG*PP
            DO IQUANTITY = 1,NR_QUANTITIES
              IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
                WRITE(LUNIT(1,IQUANTITY),'(F10.6,3E20.10)')
     &                PP, QUANTITY(IQUANTITY,0,1)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF !IF(LBIT(LBIT_PLOT_GENERAL,1)) THEN !1D

      IF(LBIT(LBIT_PLOT_GENERAL,3)) THEN !3D

!       calculate cube parameters

!       NCUBE(0)      number of atoms
!       NCUBE(1)      number of points in the "x"-direction
!       NCUBE(2)      number of points in the "y"-direction
!       NCUBE(3)      number of points in the "z"-direction

!       STEP_3D(*,0)  coordinates of initial point
!       STEP_3D(*,1)  step vector for the "x"-direction
!       STEP_3D(*,2)  step vector for the "y"-direction
!       STEP_3D(*,3)  step vector for the "z"-direction

!       ATOM_XYZ(1,I) x-coordinate of atom I
!       ATOM_XYZ(2,I) y-coordinate of atom I
!       ATOM_XYZ(3,I) z-coordinate of atom I
!       ATOM_XYZ(4,I) charge of atom I

        NCUBE(0) = NUCDEP
        CALL DZERO(STEP_3D,12)
        DO I = 1,3
          TEMP = CORD(I,1)
          STEP_3D(I,0) = TEMP
          STEP_3D(I,I) = TEMP
        ENDDO

        I = 0
        DO J = 1,NUCIND !loop over all centers
          K = ISTBNU(J)
          DO L = 0,MAXOPR
            IF(IBTAND(L,K) .EQ. 0) THEN
              I = I + 1
              DO M = 1,3
                ATOM_XYZ(M,I) = PT(IBTAND(ISYMAX(M,1),L))*CORD(M,J)
                STEP_3D(M,0) = MIN(STEP_3D(M,0),ATOM_XYZ(M,I))
                STEP_3D(M,M) = MAX(STEP_3D(M,M),ATOM_XYZ(M,I))
              ENDDO
              ATOM_XYZ(4,I) = CHARGE(J)
            ENDIF
          ENDDO
        ENDDO

!       add some space, calculate step sizes and convert step vectors to angstrom

        DO I = 1,3
          STEP_3D(I,0) = STEP_3D(I,0) - ADD_3D
          STEP_3D(I,I) = (STEP_3D(I,I) + (ADD_3D*D2))*XTANG
          STEP_3D(I,I) = (STEP_3D(I,I) - STEP_3D(I,0))/(NCUBE(I) - 1)
        ENDDO

!       write header

        DO IQUANTITY = 1,NR_QUANTITIES
          IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
            WRITE(LUNIT(3,IQUANTITY),'(A)') 'some'
            WRITE(LUNIT(3,IQUANTITY),'(A)') 'text'
            DO I = 0,3
              WRITE(LUNIT(3,IQUANTITY),'(I5,4F12.6)')
     &                                    NCUBE(I), (STEP_3D(J,I),J=1,3)
            ENDDO
            DO I = 1,NUCDEP
              WRITE(LUNIT(3,IQUANTITY),'(I5,4F12.6)')
     &         NINT(ATOM_XYZ(4,I)), ATOM_XYZ(4,I), (ATOM_XYZ(J,I),J=1,3)
            ENDDO
          ENDIF
        ENDDO

!       grid walk

        IF(LBIT(LBIT_PLOT_2D_OPTIONS,1) .OR.
     &     LBIT(LBIT_PLOT_2D_OPTIONS,2)) THEN
          NCUBE(3) = 1
        ENDIF
        IF(LBIT(LBIT_PLOT_2D_OPTIONS,3) .OR.
     &     LBIT(LBIT_PLOT_2D_OPTIONS,4)) THEN
          NCUBE(1) = 1
        ENDIF

        DO ICUBE = 1,NCUBE(1)
          DO JCUBE = 1,NCUBE(2)
            DO KCUBE = 1,NCUBE(3)

              PX = STEP_3D(1,0) + (ICUBE-1)*STEP_3D(1,1)
     &                          + (JCUBE-1)*STEP_3D(1,2)
     &                          + (KCUBE-1)*STEP_3D(1,3)
              PY = STEP_3D(2,0) + (ICUBE-1)*STEP_3D(2,1)
     &                          + (JCUBE-1)*STEP_3D(2,2)
     &                          + (KCUBE-1)*STEP_3D(2,3)
              PZ = STEP_3D(3,0) + (ICUBE-1)*STEP_3D(3,1)
     &                          + (JCUBE-1)*STEP_3D(3,2)
     &                          + (KCUBE-1)*STEP_3D(3,3)

              IF(LBIT(LBIT_PLOT_2D_OPTIONS,1)) THEN
                PX = -6.0D0+(ICUBE-1)*(12.0D0/(NCUBE(1)-1))
                PY = -6.0D0+(JCUBE-1)*(12.0D0/(NCUBE(2)-1))
                PZ = D0
              ENDIF
              IF(LBIT(LBIT_PLOT_2D_OPTIONS,2)) THEN
                PX = -2.0D0+(ICUBE-1)*(4.0D0/(NCUBE(1)-1))
                PY = -2.0D0+(JCUBE-1)*(4.0D0/(NCUBE(2)-1))
                PZ = D0
              ENDIF
              IF(LBIT(LBIT_PLOT_2D_OPTIONS,3)) THEN
                PY = -6.0D0+(KCUBE-1)*(12.0D0/(NCUBE(3)-1))
                PZ = -6.0D0+(JCUBE-1)*(12.0D0/(NCUBE(2)-1))
                PX = D0
              ENDIF
              IF(LBIT(LBIT_PLOT_2D_OPTIONS,4)) THEN
                PY = -2.0D0+(KCUBE-1)*(4.0D0/(NCUBE(3)-1))
                PZ = -2.0D0+(JCUBE-1)*(4.0D0/(NCUBE(2)-1))
                PX = D0
              ENDIF

              CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                   GAO,GAO1,GAO2,
     &                                   PX,PY,PZ,
     &                                   NCNT,BUF,
     &                                   QUANTITY(1,0,KCUBE),
     &                                   WORK,KFREE,LFREE)
            ENDDO

!           write buffer

            FORMAT_STRING = '(6E13.5)'
            IF(LBIT(LBIT_PLOT_2D_OPTIONS,5)) THEN
              FORMAT_STRING = '(200E13.5)'
            ENDIF

            DO IQUANTITY = 1,NR_QUANTITIES
              IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
                WRITE(LUNIT(3,IQUANTITY),FORMAT_STRING)
     &               (QUANTITY(IQUANTITY,0,M),M=1,NCUBE(3))
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        IF(LBIT(LBIT_PLOT_GENERAL,8)) THEN !integrate

!         3D integration

          LUNIT_3D_GRID = 44

          CALL DFTGRD(DMAT_AO,WORK(KFREE),LFREE,1)

          REWIND LUNIT_3D_GRID
          READ(LUNIT_3D_GRID) NR_POINTS

          CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

          CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &                WORK(KW),
     &                LUNIT_3D_GRID,NR_POINTS)

          CALL AROUND('3D quadrature')
          CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                            GAO,GAO1,GAO2,
     &                            WORK(KX),WORK(KY),WORK(KZ),
     &                            WORK(KW),
     &                            NR_POINTS,
     &                            NCNT,BUF,
     &                            WORK,KFREE,LFREE)

          CALL MEMREL('PLOT_SUB',WORK,KX,KX,KFREE,LFREE)

!         2D integration

          LUNIT_2D_GRID = 45

          N = 30
          DIMY = 10.0D0
          DIMZ = 10.0D0
          CALL LEGENDRE_2D_GRID(N,DIMY,DIMZ,LUNIT_2D_GRID)

          REWIND LUNIT_2D_GRID
          READ(LUNIT_2D_GRID) NR_POINTS

          CALL MEMGET('REAL',KX,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KY,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KZ,NR_POINTS,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KW,NR_POINTS,WORK,KFREE,LFREE)

          CALL REAQUA(WORK(KX),WORK(KY),WORK(KZ),
     &                WORK(KW),
     &                LUNIT_2D_GRID,NR_POINTS)

          CALL AROUND('2D Gauss-Legendre quadrature')
          CALL INTEGRATE_QUANTITY(DMAT_AO,
     &                            GAO,GAO1,GAO2,
     &                            WORK(KX),WORK(KY),WORK(KZ),
     &                            WORK(KW),
     &                            NR_POINTS,
     &                            NCNT,BUF,
     &                            WORK,KFREE,LFREE)

          CALL MEMREL('PLOT_SUB',WORK,KX,KX,KFREE,LFREE)

        ENDIF !IF(LBIT(LBIT_PLOT_GENERAL,8)) THEN !integrate
      ENDIF !IF(LBIT(LBIT_PLOT_GENERAL,3)) THEN !3D

!     close files

      DO I = 1,3 !1D, 2D, 3D
        IF(LBIT(LBIT_PLOT_GENERAL,I)) THEN
          DO IQUANTITY = 1,NR_QUANTITIES
            IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
              CLOSE(LUNIT(I,IQUANTITY),
     &              STATUS = 'KEEP')
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck GET_QUANTITY_IN_POINT */
      SUBROUTINE GET_QUANTITY_IN_POINT(DMAT_AO,
     &                                 GAO,GAO1,GAO2,
     &                                 PX,PY,PZ,
     &                                 NCNT,BUF,
     &                                 QUANTITY,
     &                                 WORK,KFREE,LFREE)
C=======================================================================
C     get quantity in point
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (D23 = 2.0D0/3.0D0,D53 = 5.0D0/3.0D0)
      PARAMETER (D32 = 3.0D0/2.0D0,D54 = 5.0D0/4.0D0)
      PARAMETER (NR_QUANTITIES = 19)

#include "visual.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbxrs.h"

      LOGICAL   LBIT
      DIMENSION DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          NCNT(*),
     &          BUF(*),
     &          QUANTITY(NR_QUANTITIES,0:3),
     &          VECTOR(3),
     &          WORK(*)

      CALL DZERO(QUANTITY,NR_QUANTITIES)

      IF(LBIT(LBIT_PLOT_QUANTITY,4)) THEN !ELF
        CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,0)

        CALL GETGRHO(VECTOR,IREP_PLOT,GAO,GAO1,DMAT_AO,BUF)
        GRADIENT_SQUARED = VECTOR(1)**D2 + VECTOR(2)**D2 + VECTOR(3)**D2
      ELSE
        CALL GETSOS(GAO,GAO1,UNUSED,UNUSED,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,0)
      ENDIF

!     density

      CALL GETRHO(DENSITY,IREP_PLOT,GAO,DMAT_AO,BUF)
      QUANTITY(1,0) = DENSITY

      IF(LBIT(LBIT_PLOT_QUANTITY,4)) THEN !ELF

        ELFKIN = D0
        DO I = 1,3
          DO IC = 1,MC
            TEMP = D0
            CALL GETRH1(TEMP,IC,IREP_PLOT,GAO1((I-1)*NTBAS(0)+1),
     &                  DMAT_AO,BUF)
            ELFKIN = ELFKIN + TEMP
          ENDDO
        ENDDO
        ELFKIN = ELFKIN/D2

        PI =  3.1415926536D0
        FACTOR = 0.3D0*((D3*PI*PI)**D23)

        H = FACTOR*(DENSITY**D53)
        G = ELFKIN - GRADIENT_SQUARED/(8.0D0*DENSITY)

        ELFCHI = G/H

        QUANTITY(4,0) = D1/(D1 + ELFCHI*ELFCHI)
!       QUANTITY(4,0) = ELF**D32 !MBE1
!       QUANTITY(4,0) = ELF**D54 !MBE2
      ENDIF

!     magnetization and current density

      DO IQUANTITY = 2,3
        IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)
          CALL GET_CURRENT(IQUANTITY,
     &                     VECTOR,
     &                     UNUSED,
     &                     IREP_PLOT,GAO,UNUSED,DMAT_AO,BUF)

          QUANTITY(IQUANTITY,1) = VECTOR(1)
          QUANTITY(IQUANTITY,2) = VECTOR(2)
          QUANTITY(IQUANTITY,3) = VECTOR(3)

          QUANTITY(IQUANTITY,0) = VECTOR(1)**D2
     &                          + VECTOR(2)**D2
     &                          + VECTOR(3)**D2
          QUANTITY(IQUANTITY,0) = DSQRT(QUANTITY(IQUANTITY,0))

          IF(LBIT(LBIT_PLOT_2D_OPTIONS,1) .OR.
     &       LBIT(LBIT_PLOT_2D_OPTIONS,2)) THEN
            WRITE(0,'(4F30.20)') PX, PY, (VECTOR(IXYZ),IXYZ=1,2)
          ENDIF
          IF(LBIT(LBIT_PLOT_2D_OPTIONS,3) .OR.
     &       LBIT(LBIT_PLOT_2D_OPTIONS,4)) THEN
            WRITE(0,'(4F30.20)') PX, PY, (VECTOR(IXYZ),IXYZ=2,3)
          ENDIF
#else
          CALL QUIT('not implemented')
#endif
        ENDIF
      ENDDO

!     r \times j

      DO IQUANTITY = 5,5
        IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
#if defined (MOD_CDFT)
          CALL GET_CURRENT(3,
     &                     VECTOR,
     &                     UNUSED,
     &                     IREP_PLOT,GAO,UNUSED,DMAT_AO,BUF)

          VECTOR_1 = PY*VECTOR(3) - PZ*VECTOR(2)
          VECTOR_2 = PZ*VECTOR(1) - PX*VECTOR(3)
          VECTOR_3 = PX*VECTOR(2) - PY*VECTOR(1)

          VECTOR(1) = VECTOR_1
          VECTOR(2) = VECTOR_2
          VECTOR(3) = VECTOR_3

          QUANTITY(IQUANTITY,1) = VECTOR(1)
          QUANTITY(IQUANTITY,2) = VECTOR(2)
          QUANTITY(IQUANTITY,3) = VECTOR(3)

          QUANTITY(IQUANTITY,0) = VECTOR(1)**D2
     &                          + VECTOR(2)**D2
     &                          + VECTOR(3)**D2
          QUANTITY(IQUANTITY,0) = DSQRT(QUANTITY(IQUANTITY,0))

          IF(LBIT(LBIT_PLOT_2D_OPTIONS,1) .OR.
     &       LBIT(LBIT_PLOT_2D_OPTIONS,2)) THEN
            WRITE(0,'(4F30.20)') PX, PY, (VECTOR(IXYZ),IXYZ=1,2)
          ENDIF
          IF(LBIT(LBIT_PLOT_2D_OPTIONS,3) .OR.
     &       LBIT(LBIT_PLOT_2D_OPTIONS,4)) THEN
            WRITE(0,'(4F30.20)') PX, PY, (VECTOR(IXYZ),IXYZ=2,3)
          ENDIF
#else
          CALL QUIT('not implemented')
#endif
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck MODIFIED_DMAT_AO */
      SUBROUTINE MODIFIED_DMAT_AO(CMO,DMAT_AO,TMAT_AO,
     &                            STRING1,STRING2,
     &                            WORK,LWORK)
C=======================================================================
C     fill MO density matrix
C     with components of the solution vector
C     these are read in from XVECFILE
C     then transform to AO density matrix
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "visual.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbxrs.h"

      CHARACTER LAB1*16,
     &          LAB2*16,
     &          TYP1*2,
     &          TYP2*2,
     &          LBUF*16,
     &          PBUF*1
      LOGICAL   FILE_EXISTS, LBIT
      DIMENSION CMO(*),
     &          DMAT_AO(*),
     &          TMAT_AO(*),
     &          STRING1(*),
     &          STRING2(*),
     &          WORK(LWORK)

      CALL QENTER('MODIFIED_DMAT_AO')

#include "memint.h"

      INQUIRE(FILE = XVECFILE,EXIST = FILE_EXISTS) !check if solution vector file is there and open it
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'solution vector file not found'
        CALL QUIT('ERROR in PLOT')
      ENDIF
      LUNIT_XVEC = 17
      OPEN(LUNIT_XVEC,
     &     FILE = XVECFILE,
     &     STATUS = 'OLD',
     &     ACCESS = 'SEQUENTIAL',
     &     FORM = 'UNFORMATTED')

!     get dimensions of record NXVEC

      NBDIM_EE    = 0
      NBDIM_EP    = 0
      NR_READS    = 0
      NR_READS_EE = 0
      NR_READS_EP = 0

      REWIND LUNIT_XVEC
  1   CONTINUE
      READ(LUNIT_XVEC) LAB1
      BACKSPACE(LUNIT_XVEC)
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 10
      ENDIF
      READ(LUNIT_XVEC) LAB1, LAB2, TYP1
      BACKSPACE(LUNIT_XVEC)

      IF(TYP1(1:1) .EQ. 'E') THEN
        NR_READS    = NR_READS + 1
        NR_READS_EE = NR_READS_EE + 1
        IF(TYP1(2:2) .EQ. '+') THEN
          IH_EE = 1
        ELSEIF(TYP1(2:2) .EQ. '-') THEN
          IH_EE = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUNIT_XVEC) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &                   JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &                   UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NR_READS_EE .EQ. NXVEC) THEN
          JSYMOP_EE = JSYMOP1
          NBDIM_EE  = LEN/NZ
        ENDIF
        READ(LUNIT_XVEC)

      ELSEIF(TYP1(1:1) .EQ. 'P') THEN
        NR_READS    = NR_READS + 1
        NR_READS_EP = NR_READS_EP + 1
        IF(TYP1(2:2) .EQ. '+') THEN
          IH_EP = 1
        ELSEIF(TYP1(2:2) .EQ. '-') THEN
          IH_EP = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUNIT_XVEC) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &                   JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &                   UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NR_READS_EP .EQ. NXVEC) THEN
          JSYMOP_EP = JSYMOP1
          NBDIM_EP  = LEN/NZ
        ENDIF
        READ(LUNIT_XVEC)

      ELSE
        GOTO 12
      ENDIF
      GOTO 1
 10   CONTINUE

!     check consistency

!     LBIT(LBIT_PLOT_GENERAL,5) skip +- (SKIPEP)
!     LBIT(LBIT_PLOT_GENERAL,6) skip ++ (SKIPEE)

      IF(NR_READS .EQ. 0)
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5) .AND.
     &   .NOT. LBIT(LBIT_PLOT_GENERAL,6) .AND.
     &   .NOT. NR_READS_EE .EQ. NR_READS_EP)
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5) .AND.
     &   NR_READS_EP .EQ. 0)
     &   GOTO 13

      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6) .AND.
     &   NR_READS_EE .EQ. 0)
     &   GOTO 14

      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5) .AND.
     &   NR_READS_EP .LT. NXVEC)
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6) .AND.
     &   NR_READS_EE .LT. NXVEC)
     &   GOTO 11

      IF(LBIT(LBIT_PLOT_GENERAL,5) .AND .
     &   LBIT(LBIT_PLOT_GENERAL,6))
     &   GOTO 15

!     error handling

      GOTO 20
 11   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'invalid NXVEC'
        WRITE(LUPRI,'(A)') 'or solution vector file corrupt'
        CALL QUIT('ERROR in PLOT')
 12   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'invalid operator label'
        WRITE(LUPRI,'(A)') 'solution vector file corrupt'
        CALL QUIT('ERROR in PLOT')
 13   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'EP requested but no EP rotation'
        WRITE(LUPRI,'(A)') 'on solution vector file'
        CALL QUIT('ERROR in PLOT')
 14   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'EE requested but no EE rotation'
        WRITE(LUPRI,'(A)') 'on solution vector file'
        CALL QUIT('ERROR in PLOT')
 15   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'SKIPEP = SKIPEE = .TRUE.'
        CALL QUIT('ERROR in PLOT')
 20   CONTINUE

!     allocate memory

      CALL MEMGET('REAL',KBVEC_EE ,NBDIM_EE*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC_EP ,NBDIM_EP*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT    ,N2ORBXQ,
     &            WORK,KFREE,LFREE)

!     now read the vector

      REWIND LUNIT_XVEC
      DO I = 1,NXVEC
        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6)) THEN ! "EE" part
          READ(LUNIT_XVEC)
          IF(I .EQ. NXVEC) THEN
            CALL DZERO(WORK(KBVEC_EE),NBDIM_EE*NZ)
            CALL READT(LUNIT_XVEC,NBDIM_EE*NZ,WORK(KBVEC_EE))
          ELSE
            READ(LUNIT_XVEC)
          ENDIF
        ELSEIF(NR_READS_EE .GT. 0) THEN
          READ(LUNIT_XVEC)
          READ(LUNIT_XVEC)
        ENDIF

        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5)) THEN ! "EP" part
          READ(LUNIT_XVEC)
          IF(I .EQ. NXVEC) THEN
            CALL DZERO(WORK(KBVEC_EP),NBDIM_EP*NZ)
            CALL READT(LUNIT_XVEC,NBDIM_EP*NZ,WORK(KBVEC_EP))
          ELSE
            READ(LUNIT_XVEC)
          ENDIF
        ELSEIF(NR_READS_EP .GT. 0) THEN
          READ(LUNIT_XVEC)
          READ(LUNIT_XVEC)
        ENDIF
      ENDDO

      CLOSE(LUNIT_XVEC,STATUS = 'KEEP') !close solution vector file

!     scatter vector onto MO density matrix

      IHERMITICITY = IH_EE
      ICOUNTER_EE  = 1
      ICOUNTER_EP  = 1
      CALL DZERO(WORK(KBMAT),N2ORBXQ)

      DO IZ = 1,NZ
        IF(NFSYM .EQ. 1) THEN
          CALL SCATTER_BVEC(WORK(KBVEC_EE),WORK(KBVEC_EP),
     &                      ICOUNTER_EE,ICOUNTER_EP,
     &                      WORK(KBMAT),
     &                      IZ,
     &                      1,1)
        ELSE
          IF(IHERMITICITY .EQ. 1) THEN
            CALL SCATTER_BVEC(WORK(KBVEC_EE),WORK(KBVEC_EP),
     &                        ICOUNTER_EE,ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,2)
            CALL SCATTER_BVEC(WORK(KBVEC_EE),WORK(KBVEC_EP),
     &                        ICOUNTER_EE,ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,1)
          ELSE
            CALL SCATTER_BVEC(WORK(KBVEC_EE),WORK(KBVEC_EP),
     &                        ICOUNTER_EE,ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,1)
            CALL SCATTER_BVEC(WORK(KBVEC_EE),WORK(KBVEC_EP),
     &                        ICOUNTER_EE,ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,2)
          ENDIF
        ENDIF
      ENDDO
      CALL DSCAL(N2ORBXQ,-D1,WORK(KBMAT),1)

      JSYMOP    = JSYMOP_EE
      IREP_PLOT = JSYMOP - 1
      JOPSY     = JBTOF(IREP_PLOT,1)

      DO IZ = 1,NZ
        IOFF = NORBT*NORBT*(IZ - 1)
        CALL FULMAT('S',NORBT,NORBT,WORK(KBMAT+IOFF)) !symmetrize
      ENDDO

!     get AO density matrix

      CALL DZERO(DMAT_AO,N2BBASXQ)
      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6)) THEN !if not skip ++ (SKIPEE)
        CALL DMOLI3(0,IH_EE,DMAT_AO,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
      ENDIF
      IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5)) THEN !if not skip +- (SKIPEP)
        CALL DZERO(TMAT_AO,N2BBASXQ)
        CALL DMOLI3(1,IH_EP,TMAT_AO,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
        CALL DAXPY(N2BBASXQ,D1,TMAT_AO,1,DMAT_AO,1)
      ENDIF

      IF(NZ .LT. 4) THEN
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREP_PLOT)
          IOFF = NTBAS(0)*NTBAS(0)*(IZ-1)
          CALL Q2BPHASE('D',IQ,1,DMAT_AO(1+IOFF)) !insert half-phase
        ENDDO
      ENDIF

!     release memory

      CALL MEMREL('MODIFIED_DMAT_AO',WORK,KBVEC_EE,KBVEC_EE,KFREE,LFREE)
      CALL QEXIT('MODIFIED_DMAT_AO')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck SCATTER_BVEC */
      SUBROUTINE SCATTER_BVEC(BVEC_EE,BVEC_EP,
     &                        ICOUNTER_EE,ICOUNTER_EP,
     &                        BMAT,
     &                        IZ,
     &                        IFSYM1,IFSYM2)
C=======================================================================
C     scatter solution vector onto BMAT
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "visual.h"
#include "dgroup.h"
#include "dcbbas.h"

      LOGICAL   LBIT
      DIMENSION BVEC_EE(*),
     &          BVEC_EP(*),
     &          BMAT(NORBT,NORBT,NZ),
     &          NEI(2),
     &          NES(2),
     &          NP(2)

      DO I = 1,2
        NEI(I) = NISH(I)
        NES(I) = NESH(I) - NISH(I)
        NP(I)  = NPSH(I)
      ENDDO

      II1 = 1
      II1 = II1 + (NESH(1) + NPSH(1))*(IFSYM1 - 1)
      II1 = II1 + NP(IFSYM1)
      II2 = II1 + NEI(IFSYM1) - 1
      IS1 = 1
      IS1 = IS1 + (NESH(1) + NPSH(1))*(IFSYM2 - 1)
      IS2 = IS1 + NP(IFSYM2)  - 1
      IS3 = IS2 + NEI(IFSYM2) + 1
      IS4 = IS3 + NES(IFSYM2) - 1

      DO II = II1,II2
        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5)) THEN !if not skip +- (SKIPEP)
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVEC_EP(ICOUNTER_EP)
            ICOUNTER_EP = ICOUNTER_EP + 1
          ENDDO
        ENDIF
        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6)) THEN !if not skip ++ (SKIPEE)
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVEC_EE(ICOUNTER_EE)
            ICOUNTER_EE = ICOUNTER_EE + 1
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DMAT_MO_2_AO */
      SUBROUTINE DMAT_MO_2_AO(CMO,DMAT_AO,DMAT_MO,WORK,LWORK)
C=======================================================================
C     transform MO to AO density matrix
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION CMO(*),
     &          DMAT_AO(*),
     &          DMAT_MO(*),
     &          WORK(LWORK)

#include "memint.h"

      CALL DZERO(DMAT_AO,N2BBASXQ)

      DO I = 1,NFSYM
        CALL QTRANS('MOAO','S',D0,
     &              NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &              DMAT_AO(1+I2BASX(I,I)),
     &              NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &              DMAT_MO(1+I2ORBX(I,I)),
     &              NORBT,NORBT,1,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,0)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck OCCUPY_DMAT_MO */
      SUBROUTINE OCCUPY_DMAT_MO(DMAT_MO,STRING1,STRING2)
C=======================================================================
C     occupy MO density matrix
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION DMAT_MO(*),
     &          STRING1(*),
     &          STRING2(*)

      CALL DZERO(DMAT_MO,N2ORBXQ)

      K = 0
      DO I = 1,NFSYM
        IF(NISH(I) .GT. 0) THEN
          DO J = 1,NESH(I)
            K = K + 1
            JOFF = (NORB(1)*(I-1) + NPSH(I))*NORBT
            JOFF = JOFF + NORB(1)*(I-1) + NPSH(I)
            JOFF = JOFF + (J-1)*NORBT + J - 1

            IF(I .EQ. 1) DMAT_MO(1+JOFF) = STRING1(J)
            IF(I .EQ. 2) DMAT_MO(1+JOFF) = STRING2(J)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck INTEGRATE_QUANTITY */
      SUBROUTINE INTEGRATE_QUANTITY(DMAT_AO,
     &                              GAO,GAO1,GAO2,
     &                              PX,PY,PZ,
     &                              WEIGHT,
     &                              NR_POINTS,
     &                              NCNT,BUF,
     &                              WORK,KFREE,LFREE)
C=======================================================================
C     integrate
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                            last revision: october 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NR_QUANTITIES = 19)

#include "visual.h"

      CHARACTER QUANTITY_NAME(NR_QUANTITIES)*30
      LOGICAL   LBIT,
     &          QUANTITY_IS_A_VECTOR(NR_QUANTITIES)
      DIMENSION QUANTITY_INTEGRAL(NR_QUANTITIES,0:3),
     &          QUANTITY(NR_QUANTITIES,0:3),
     &          DMAT_AO(*),
     &          GAO(*),
     &          GAO1(*),
     &          GAO2(*),
     &          PX(*),
     &          PY(*),
     &          PZ(*),
     &          WEIGHT(*),
     &          NCNT(*),
     &          BUF(*),
     &          WORK(*)

      QUANTITY_NAME(1) = 'density'
      QUANTITY_NAME(2) = 'spin_density'
      QUANTITY_NAME(3) = 'current_density'
      QUANTITY_NAME(4) = 'ELF'
      QUANTITY_NAME(5) = 'r_times_j'
      QUANTITY_IS_A_VECTOR(1) = .FALSE.
      QUANTITY_IS_A_VECTOR(2) = .TRUE.
      QUANTITY_IS_A_VECTOR(3) = .TRUE.
      QUANTITY_IS_A_VECTOR(4) = .FALSE.
      QUANTITY_IS_A_VECTOR(5) = .TRUE.

      QUANTITY_NAME(10) = 'opt0'
      QUANTITY_NAME(11) = 'opt1'
      QUANTITY_NAME(12) = 'opt2'
      QUANTITY_NAME(13) = 'opt3'
      QUANTITY_NAME(14) = 'opt4'
      QUANTITY_NAME(15) = 'opt5'
      QUANTITY_NAME(16) = 'opt6'
      QUANTITY_NAME(17) = 'opt7'
      QUANTITY_NAME(18) = 'opt8'
      QUANTITY_NAME(19) = 'opt9'
      QUANTITY_IS_A_VECTOR(10) = .FALSE.
      QUANTITY_IS_A_VECTOR(11) = .FALSE.
      QUANTITY_IS_A_VECTOR(12) = .FALSE.
      QUANTITY_IS_A_VECTOR(13) = .FALSE.
      QUANTITY_IS_A_VECTOR(14) = .FALSE.
      QUANTITY_IS_A_VECTOR(15) = .FALSE.
      QUANTITY_IS_A_VECTOR(16) = .FALSE.
      QUANTITY_IS_A_VECTOR(17) = .FALSE.
      QUANTITY_IS_A_VECTOR(18) = .FALSE.
      QUANTITY_IS_A_VECTOR(19) = .FALSE.

      CALL DZERO(QUANTITY_INTEGRAL,NR_QUANTITIES*4)

      DO IPOINT = 1,NR_POINTS
        CALL GET_QUANTITY_IN_POINT(DMAT_AO,
     &                             GAO,GAO1,GAO2,
     &                             PX(IPOINT),
     &                             PY(IPOINT),
     &                             PZ(IPOINT),
     &                             NCNT,BUF,
     &                             QUANTITY,
     &                             WORK,KFREE,LFREE)

        DO IQUANTITY = 1,NR_QUANTITIES
          IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
            QUANTITY_INTEGRAL(IQUANTITY,0)
     &                          = QUANTITY_INTEGRAL(IQUANTITY,0)
C    &                          + QUANTITY(IQUANTITY,0)*WEIGHT(IPOINT)
     &                          + WEIGHT(IPOINT)
     &                          * EXP(-(PY(IPOINT)**D2))
     &                          * EXP(-(PZ(IPOINT)**D2))

            IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
              DO IXYZ = 1,3
                QUANTITY_INTEGRAL(IQUANTITY,IXYZ)
     &                       = QUANTITY_INTEGRAL(IQUANTITY,IXYZ)
     &                       + QUANTITY(IQUANTITY,IXYZ)*WEIGHT(IPOINT)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      WRITE(LUPRI,'(A8,A29,3A20)') 'quantity',
     &                             'scalar/norm',
     &                             'x-component',
     &                             'y-component',
     &                             'z-component'

      WRITE(LUPRI,'(A)')

      DO IQUANTITY = 1,NR_QUANTITIES
        IF(LBIT(LBIT_PLOT_QUANTITY,IQUANTITY)) THEN
          IF(QUANTITY_IS_A_VECTOR(IQUANTITY)) THEN
            WRITE(LUPRI,'(A20,4E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,0),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,1),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,2),
     &                                   QUANTITY_INTEGRAL(IQUANTITY,3)
          ELSE
            WRITE(LUPRI,'(A20,E20.10)') QUANTITY_NAME(IQUANTITY),
     &                                  QUANTITY_INTEGRAL(IQUANTITY,0)
          ENDIF
        ENDIF
      ENDDO

      WRITE(LUPRI,'(A)')
      CALL PRSYMB(LUPRI,'-',100,0)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_POLYNOMIAL */
      SUBROUTINE LEGENDRE_POLYNOMIAL(P,N,X)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX)

      P(0) = 1
      P(1) = X

      DO I = 1,N-1
        P(I+1) = (2*I + 1)*X*P(I) - I*P(I-1)
        P(I+1) = P(I+1)/(I+1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_POLYNOMIAL_DERIVATIVES */
      SUBROUTINE LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      DIMENSION P(0:NMAX),
     &          DP(NMAX)

      DO I = 1,N
        DP(I) = (I*X*P(I) - I*P(I-1))/(X*X - D1)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LOBATTO_1D_GRID */
      SUBROUTINE LOBATTO_1D_GRID(N,X,W)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX),
     &          W(NMAX),
     &          P(0:NMAX),
     &          DP(NMAX)

      X_OLD = -D1
      CALL LEGENDRE_POLYNOMIAL(P,N,X_OLD)
      CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X_OLD)

      DO I = 2,N-1
        STEP  = 1.0D-4
        DO WHILE(DABS(STEP) .GT. THRESHOLD)
          Y_OLD = DP(N-1)
          X_NEW = X_OLD + STEP

          CALL LEGENDRE_POLYNOMIAL(P,N,X_NEW)
          CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X_NEW)
          Y_NEW = DP(N-1)

          IF((Y_OLD/Y_NEW) .LT. D0) THEN
            STEP = -STEP/10.0D0
          ENDIF

          X_OLD = X_NEW
        ENDDO

 1      CONTINUE
        X(I) = X_NEW

        W(I) = N*(N - 1)*P(N-1)*P(N-1)
        W(I) = D2/W(I)
      ENDDO

!     end points

      X(1) = -D1
      X(N) =  D1

      W(1) = D2/(N*(N - D1))
      W(N) = W(1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_1D_GRID */
      SUBROUTINE LEGENDRE_1D_GRID(N,X,W)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)
      PARAMETER (THRESHOLD = 1.0D-12)

      DIMENSION X(NMAX),
     &          W(NMAX),
     &          P(0:NMAX),
     &          DP(NMAX)

      X_OLD = -D1
      CALL LEGENDRE_POLYNOMIAL(P,N,X_OLD)

      DO I = 1,N
        STEP  = 1.0D-4
        DO WHILE(DABS(STEP) .GT. THRESHOLD)
          Y_OLD = P(N)
          X_NEW = X_OLD + STEP

          CALL LEGENDRE_POLYNOMIAL(P,N,X_NEW)
          Y_NEW = P(N)

          IF((Y_OLD/Y_NEW) .LT. D0) THEN
            STEP = -STEP/10.0D0
          ENDIF

          X_OLD = X_NEW
        ENDDO

 1      CONTINUE
        X(I) = X_NEW

        CALL LEGENDRE_POLYNOMIAL_DERIVATIVES(P,DP,N,X(I))
        W(I) = (D1 - X(I)*X(I))*DP(N)*DP(N)
        W(I) = D2/W(I)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck LEGENDRE_2D_GRID */
      SUBROUTINE LEGENDRE_2D_GRID(N,DIMI,DIMJ,LUNIT)
C=======================================================================
C                                                           radovan bast
C                                           last revision: november 2006
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (NMAX = 1000)

      CHARACTER*30 FILE_NAME
      DIMENSION PX(NMAX*NMAX),
     &          PY(NMAX*NMAX),
     &          PZ(NMAX*NMAX),
     &          WEIGHT(NMAX*NMAX),
     &          X(NMAX),
     &          W(NMAX)

      WRITE(FILE_NAME,'(A5,I2)') 'GRID.', LUNIT

      OPEN(LUNIT,
     &     FILE = FILE_NAME,
     &     STATUS = 'UNKNOWN',
     &     ACCESS = 'SEQUENTIAL',
     &     FORM = 'UNFORMATTED')
      REWIND LUNIT

      CALL LEGENDRE_1D_GRID(N,X,W)

      IPOINT = 0
      DO I = 1,N
        DO J = 1,N
          IPOINT = IPOINT + 1

          PX(IPOINT) = D0
          PY(IPOINT) = (X(I)+D1)*DIMI/D2
          PZ(IPOINT) = (X(J)+D1)*DIMJ/D2 - DIMJ/D2

          WEIGHT(IPOINT) = W(I)*W(J)*(DIMI/D2)*(DIMJ/D2)
        ENDDO
      ENDDO
      NR_POINTS = IPOINT

      CALL WRTQUA(PX,PY,PZ,WEIGHT,LUNIT,NR_POINTS)

      RETURN
      END
