C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck PLOT */
      SUBROUTINE PLOT(WORK,LWORK)
C=======================================================================
C     prepare {1,2,3}-dimensional plot
C     adaptation of Trond's plot routines
C-----------------------------------------------------------------------
C     radovan bast                           last revision: october 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
!-----------------------------------------------------------------------
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "visual.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"

      LOGICAL   TOBE
      DIMENSION WORK(LWORK)

      CALL QENTER('PLOT')

#include "memint.h"

      CALL TITLER('{1,2,3}-dimensional plot','*',122)

      INQUIRE(FILE = 'DFCOEF',EXIST = TOBE) !check if coefficients are on file
      IF(.NOT. TOBE) THEN
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
      CALL MEMGET('REAL',KSTRING1 ,NISH(1),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING2 ,NISH(2),
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
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),DUMMY,
     &            DUMMY,6,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')

      CALL FNDOCC(WORK(KEIG),WORK,KFREE,LFREE) !find occupation


!     zero out orbitals?

      CALL DFILL(WORK(KSTRING1),D1,NISH(1))
      IF(NFSYM .GT. 1) CALL DFILL(WORK(KSTRING2),D1,NISH(2))

      IF(USE_ORBSTRING) THEN
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
     &     .OR. (.NOT. (ISTRING1_END .LE. NISH(1)))) THEN
          WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
          WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
          CALL QUIT('ERROR in PLOT')
        ENDIF
        CALL DFILL(WORK(KSTRING1),D0,NISH(1))
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
     &       .OR. (.NOT. (ISTRING2_END .LE. NISH(1)))) THEN
            WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
            WRITE(LUPRI,'(A)') 'wrong orbital string dimension'
            CALL QUIT('ERROR in PLOT')
          ENDIF
          CALL DFILL(WORK(KSTRING2),D0,NISH(2))
          DO I = ISTRING2_BEGIN,ISTRING2_END
            IF(ISTRING2_BEGIN .NE. 0) WORK(KSTRING2-1+I) = D1
          ENDDO
        ENDIF !IF(NFSYM .GT. 1) THEN
      ENDIF !IF(USE_ORBSTRING) THEN


!     generate density matrix

      IF(READXVECFILE) THEN
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
     &              WORK,KFREE,LFREE)


!     release memory

      CALL MEMREL('PLOT',WORK,KCMO,KCMO,KFREE,LFREE)
      CALL QEXIT('PLOT')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck plot */
      SUBROUTINE PLOT_SUB(DMAT_AO,
     &                GAO,GAO1,GAO2,
     &                NCNT,BUF,
     &                WORK,KFREE,LFREE)
C  .....................................................................
C  .  description:                                                     .
C  .  do the plot                                                      .
C  .  copy-and-paste-and-modify-a-little-job                           .
C  .  used Trond's respective routines                                 .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "codata.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "dcbgen.h"
C
      CHARACTER PAIRNAMEA(100)*3,PAIRNAMEB(100)*3
      CHARACTER*20 FORMATSTRING
#include "visual.h"
#include "dcborb.h"
C
      DIMENSION DMAT_AO(*),
     &          GAO(*),GAO1(*),GAO2(*),
     &          NCNT(*),BUF(*),
     &          WORK(*)
C
      DIMENSION DENSITY(0:2,NCUBE(3)),
     &          SPINDENSITY(NCUBE(3)),
     &          CURRENTDENSITY(NCUBE(3)),
     &          ELF(NCUBE(3)),
     &          DLAPLACIAN(NCUBE(3)),
     &          RXJ(NCUBE(3)),
     &          PAIRCOORA(100,3),PAIRCOORB(100,3),
     &          D3DSTEP(3,0:3),ATOMXYZ(4,NUCDEP)
C
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"
#include "sdpre.h"
#include "ibtfun.h"
C
C
      IPRINT = IPLOT
C
C
C     open files
C     ==========
C
      LU1D = 13
      LU3D = 14
C
      IF(DO1DPLOT) THEN
        IF(PLOTDENSITY) THEN
          LU1DDENSITY = LU1D
          LU1D = LU1D+1
          OPEN(LU1DDENSITY,
     &         FILE = 'PLOT1D.DENSITY',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTELF) THEN
          LU1DELF = LU1D
          LU1D = LU1D+1
          OPEN(LU1DELF,
     &         FILE = 'PLOT1D.ELF',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTLAPLACIAN) THEN
          LU1DLAPLACIAN = LU1D
          LU1D = LU1D+1
          OPEN(LU1DLAPLACIAN,
     &         FILE = 'PLOT1D.LAPLACIAN',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
      ENDIF
      IF(DO3DPLOT) THEN
        IF(PLOTDENSITY) THEN
          LU3DDENSITY = LU3D
          LU3D = LU3D+1
          OPEN(LU3DDENSITY,
     &         FILE = 'PLOT3D.DENSITY.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTELF) THEN
          LU3DELF = LU3D
          LU3D = LU3D+1
          OPEN(LU3DELF,
     &         FILE = 'PLOT3D.ELF.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTLAPLACIAN) THEN
          LU3DLAPLACIAN = LU3D
          LU3D = LU3D+1
          OPEN(LU3DLAPLACIAN,
     &         FILE = 'PLOT3D.LAPLACIAN.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTSPINDENSITY) THEN
          LU3DSPINDENSITY = LU3D
          LU3D = LU3D+1
          OPEN(LU3DSPINDENSITY,
     &         FILE = 'PLOT3D.SPINDENSITY.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTCURRENTDENSITY) THEN
          LU3DCURRENTDENSITY = LU3D
          LU3D = LU3D+1
          OPEN(LU3DCURRENTDENSITY,
     &         FILE = 'PLOT3D.CURRENTDENSITY.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
        IF(PLOTRXJ) THEN
          LU3DRXJ = LU3D
          LU3D = LU3D+1
          OPEN(LU3DRXJ,
     &         FILE = 'PLOT3D.RXJ.cube',STATUS = 'NEW',
     &         ACCESS = 'SEQUENTIAL',FORM = 'FORMATTED')
        ENDIF
      ENDIF
C
C
C     1D ?
C     ====
C
      IF(DO1DPLOT) THEN
C
C
C       loop over pairs of nuclear centers
C       ==================================
C
        NPAIRS  = 0
        DO I = 1,NUCDEP
          DO J = I+1,NUCDEP
C
            NPAIRS = NPAIRS+1
C
            WRITE(PAIRNAMEA(NPAIRS),'(A2,I1)') NAMDEP(I), I
            PAIRCOORA(NPAIRS,1) = CORD(1,I)
            PAIRCOORA(NPAIRS,2) = CORD(2,I)
            PAIRCOORA(NPAIRS,3) = CORD(3,I)
C
            WRITE(PAIRNAMEB(NPAIRS),'(A2,I1)') NAMDEP(J), J
            PAIRCOORB(NPAIRS,1) = CORD(1,J)
            PAIRCOORB(NPAIRS,2) = CORD(2,J)
            PAIRCOORB(NPAIRS,3) = CORD(3,J)
C
          ENDDO
        ENDDO
C
C
C       loop over pairs
C       ===============
C
        DO IPAIR = 1,NPAIRS
C
C
C         write pair label to file
C         ========================
C
          IF(PLOTDENSITY) THEN
            WRITE(LU1DDENSITY,'(A5)') '*****'
            WRITE(LU1DDENSITY,'(A)')
            WRITE(LU1DDENSITY,'(A3,A1,A3)') PAIRNAMEA(IPAIR), '-',
     &                                      PAIRNAMEB(IPAIR)
            WRITE(LU1DDENSITY,'(A)')
          ENDIF
          IF(PLOTELF) THEN
            WRITE(LU1DELF,'(A5)') '*****'
            WRITE(LU1DELF,'(A)')
            WRITE(LU1DELF,'(A3,A1,A3)') PAIRNAMEA(IPAIR), '-',
     &                                  PAIRNAMEB(IPAIR)
            WRITE(LU1DELF,'(A)')
          ENDIF
          IF(PLOTLAPLACIAN) THEN
            WRITE(LU1DLAPLACIAN,'(A5)') '*****'
            WRITE(LU1DLAPLACIAN,'(A)')
            WRITE(LU1DLAPLACIAN,'(A3,A1,A3)') PAIRNAMEA(IPAIR), '-',
     &                                        PAIRNAMEB(IPAIR)
            WRITE(LU1DLAPLACIAN,'(A)')
          ENDIF
C
C
C         calculate distance and step
C         ===========================
C
          DX = PAIRCOORB(IPAIR,1)-PAIRCOORA(IPAIR,1)
          DY = PAIRCOORB(IPAIR,2)-PAIRCOORA(IPAIR,2)
          DZ = PAIRCOORB(IPAIR,3)-PAIRCOORA(IPAIR,3)
          DISTANCE = DSQRT(DX*DX + DY*DY + DZ*DZ)
          DX       = DX/DISTANCE
          DY       = DY/DISTANCE
          DZ       = DZ/DISTANCE
C
C
C         -> angstrom
C         ===========
C
          ASTEP = D1DSTEP/XTANG
          NPOINTS = INT(DISTANCE/ASTEP)
          IF(MOD(NPOINTS,2).EQ.1) NPOINTS = NPOINTS+1
          ASTEP = DISTANCE/SDPRE(NPOINTS)
C
C
C         loop along line
C         ===============
C
          IPOINT1 = -NPOINTS/2
          IPOINT2 = NPOINTS-IPOINT1
          DO IPOINT = IPOINT1,IPOINT2
            PP  = SDPRE(IPOINT)*ASTEP
            PX  = PAIRCOORA(IPAIR,1)+PP*DX
            PY  = PAIRCOORA(IPAIR,2)+PP*DY
            PZ  = PAIRCOORA(IPAIR,3)+PP*DZ
C
C
C           get necessary stuff
C           ===================
C
            CALL GETQUANTITY(DMAT_AO,
     &                       GAO,GAO1,GAO2,
     &                       PX,PY,PZ,
     &                       NCNT,BUF,IPRINT,
     &                       DENSITY,
     &                       SPINDENSITY,
     &                       CURRENTDENSITY,
     &                       ELF,
     &                       DLAPLACIAN,
     &                       RXJ,
     &                       WORK,KFREE,LFREE)
C
C
C           write to file
C           =============
C
            PP  = XTANG*PP
C
            IF(PLOTDENSITY) THEN
              WRITE(LU1DDENSITY,'(F10.6,3E20.10)')
     &                                         PP, (DENSITY(I,1),I=0,MC)
            ENDIF
            IF(PLOTELF) THEN
              WRITE(LU1DELF,'(F10.6,3E20.10)') PP, ELF(1)
            ENDIF
            IF(PLOTLAPLACIAN) THEN
              WRITE(LU1DLAPLACIAN,'(F10.6,3E20.10)') PP, DLAPLACIAN(1)
            ENDIF
C
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C
C     3D ?
C     ====
C
      IF(DO3DPLOT) THEN
C
C
C       calculate cube parameters
C       =========================
C
C       NCUBE(0)     - number of atoms
C       NCUBE(1)     - number of points in the "x"-direction
C       NCUBE(2)     - number of points in the "y"-direction
C       NCUBE(3)     - number of points in the "z"-direction
C
C       D3DSTEP(*,0) - coordinates of initial point
C       D3DSTEP(*,1) - step vector for the "x"-direction
C       D3DSTEP(*,2) - step vector for the "y"-direction
C       D3DSTEP(*,3) - step vector for the "z"-direction
C
C       ATOMXYZ(1,I) - x-coordinate of atom I
C       ATOMXYZ(2,I) - y-coordinate of atom I
C       ATOMXYZ(3,I) - z-coordinate of atom I
C       ATOMXYZ(4,I) - charge of atom I
C
        NCUBE(0) = NUCDEP
        CALL DZERO(D3DSTEP,12)
        DO I = 1,3
          TEMP = CORD(I,1)
          D3DSTEP(I,0) = TEMP
          D3DSTEP(I,I) = TEMP
        ENDDO
C
C
C       loop over all centers
C       =====================
C
        I = 0
        DO J = 1,NUCIND
          K = ISTBNU(J)
          DO L = 0,MAXOPR
            IF(IBTAND(L,K).EQ.0) THEN
              I = I+1
              DO M = 1,3
                ATOMXYZ(M,I) = PT(IBTAND(ISYMAX(M,1),L))*CORD(M,J)
                D3DSTEP(M,0) = MIN(D3DSTEP(M,0),ATOMXYZ(M,I))
                D3DSTEP(M,M) = MAX(D3DSTEP(M,M),ATOMXYZ(M,I))
              ENDDO
              ATOMXYZ(4,I) = CHARGE(J)
            ENDIF
          ENDDO
        ENDDO
C
C
C       add some space, calculate step sizes and convert step vectors to angstrom
C       =========================================================================
C
        DO I = 1,3
          D3DSTEP(I,0) = D3DSTEP(I,0)-D3DSTART
          D3DSTEP(I,I) = (D3DSTEP(I,I)+D3DNORM)*XTANG
          D3DSTEP(I,I) = (D3DSTEP(I,I)-D3DSTEP(I,0))/(NCUBE(I)-1)
        ENDDO
C
C
C       write header
C       ============
C
        IF(PLOTDENSITY) THEN
          WRITE(LU3DDENSITY,'(A)') 'some'
          WRITE(LU3DDENSITY,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DDENSITY,'(I5,4F12.6)')
     &                                    NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DDENSITY,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
        IF(PLOTELF) THEN
          WRITE(LU3DELF,'(A)') 'some'
          WRITE(LU3DELF,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DELF,'(I5,4F12.6)') NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DELF,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
        IF(PLOTLAPLACIAN) THEN
          WRITE(LU3DLAPLACIAN,'(A)') 'some'
          WRITE(LU3DLAPLACIAN,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DLAPLACIAN,'(I5,4F12.6)')
     &                                    NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DLAPLACIAN,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
        IF(PLOTSPINDENSITY) THEN
          WRITE(LU3DSPINDENSITY,'(A)') 'some'
          WRITE(LU3DSPINDENSITY,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DSPINDENSITY,'(I5,4F12.6)')
     &                                    NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DSPINDENSITY,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
        IF(PLOTCURRENTDENSITY) THEN
          WRITE(LU3DCURRENTDENSITY,'(A)') 'some'
          WRITE(LU3DCURRENTDENSITY,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DCURRENTDENSITY,'(I5,4F12.6)')
     &                                    NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DCURRENTDENSITY,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
        IF(PLOTRXJ) THEN
          WRITE(LU3DRXJ,'(A)') 'some'
          WRITE(LU3DRXJ,'(A)') 'text'
          DO I = 0,3
            WRITE(LU3DRXJ,'(I5,4F12.6)') NCUBE(I), (D3DSTEP(J,I),J=1,3)
          ENDDO
          DO I = 1,NUCDEP
            WRITE(LU3DRXJ,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
C
C
C       do the grid walk
C       ================
C
        IF(CUT_XY.OR.CUT_XY_SMALL) NCUBE(3) = 1
        IF(CUT_YZ.OR.CUT_YZ_SMALL) NCUBE(1) = 1
C
        DO ICUBE = 1,NCUBE(1)
          DO JCUBE = 1,NCUBE(2)
            DO KCUBE = 1,NCUBE(3)
C
              PX = D3DSTEP(1,0)+(ICUBE-1)*D3DSTEP(1,1)
     &                         +(JCUBE-1)*D3DSTEP(1,2)
     &                         +(KCUBE-1)*D3DSTEP(1,3)
              PY = D3DSTEP(2,0)+(ICUBE-1)*D3DSTEP(2,1)
     &                         +(JCUBE-1)*D3DSTEP(2,2)
     &                         +(KCUBE-1)*D3DSTEP(2,3)
              PZ = D3DSTEP(3,0)+(ICUBE-1)*D3DSTEP(3,1)
     &                         +(JCUBE-1)*D3DSTEP(3,2)
     &                         +(KCUBE-1)*D3DSTEP(3,3)
C
              IF(CUT_XY) THEN
                PX = -6.0D0+(ICUBE-1)*(12.0D0/(NCUBE(1)-1))
                PY = -6.0D0+(JCUBE-1)*(12.0D0/(NCUBE(2)-1))
                PZ = D0
              ENDIF
C
              IF(CUT_XY_SMALL) THEN
                PX = -2.0D0+(ICUBE-1)*(4.0D0/(NCUBE(1)-1))
                PY = -2.0D0+(JCUBE-1)*(4.0D0/(NCUBE(2)-1))
                PZ = D0
              ENDIF
C
              IF(CUT_YZ) THEN
                PY = -6.0D0+(KCUBE-1)*(12.0D0/(NCUBE(3)-1))
                PZ = -6.0D0+(JCUBE-1)*(12.0D0/(NCUBE(2)-1))
                PX = D0
              ENDIF
C
              IF(CUT_YZ_SMALL) THEN
                PY = -2.0D0+(KCUBE-1)*(4.0D0/(NCUBE(3)-1))
                PZ = -2.0D0+(JCUBE-1)*(4.0D0/(NCUBE(2)-1))
                PX = D0
              ENDIF
C
C
C             get necessary stuff
C             ===================
C
              CALL GETQUANTITY(DMAT_AO,
     &                         GAO,GAO1,GAO2,
     &                         PX,PY,PZ,
     &                         NCNT,BUF,IPRINT,
     &                         DENSITY(0,KCUBE),
     &                         SPINDENSITY(KCUBE),
     &                         CURRENTDENSITY(KCUBE),
     &                         ELF(KCUBE),
     &                         DLAPLACIAN(KCUBE),
     &                         RXJ(KCUBE),
     &                         WORK,KFREE,LFREE)
C
C
            ENDDO
C
C
C           write buffer
C           ============
C
            FORMATSTRING = '(6E13.5)'
            IF(FORMAT_200) FORMATSTRING = '(200E13.5)'
C
            IF(PLOTDENSITY) THEN
             WRITE(LU3DDENSITY,FORMATSTRING)
     &            (DENSITY(0,M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTELF) THEN
              WRITE(LU3DELF,FORMATSTRING)
     &             (ELF(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTLAPLACIAN) THEN
              WRITE(LU3DLAPLACIAN,FORMATSTRING)
     &             (DLAPLACIAN(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTSPINDENSITY) THEN
              WRITE(LU3DSPINDENSITY,FORMATSTRING)
     &             (SPINDENSITY(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTCURRENTDENSITY) THEN
              WRITE(LU3DCURRENTDENSITY,FORMATSTRING)
     &             (CURRENTDENSITY(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTRXJ) THEN
              WRITE(LU3DRXJ,FORMATSTRING)
     &             (RXJ(M),M=1,NCUBE(3))
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C
C     close files
C     ===========
C
      IF(DO1DPLOT) THEN
        IF(PLOTDENSITY) THEN
          CLOSE(LU1DDENSITY,STATUS = 'KEEP')
        ENDIF
        IF(PLOTELF) THEN
          CLOSE(LU1DELF,STATUS = 'KEEP')
        ENDIF
        IF(PLOTLAPLACIAN) THEN
          CLOSE(LU1DLAPLACIAN,STATUS = 'KEEP')
        ENDIF
      ENDIF
      IF(DO3DPLOT) THEN
        IF(PLOTDENSITY) THEN
          CLOSE(LU3DDENSITY,STATUS = 'KEEP')
        ENDIF
        IF(PLOTELF) THEN
          CLOSE(LU3DELF,STATUS = 'KEEP')
        ENDIF
        IF(PLOTLAPLACIAN) THEN
          CLOSE(LU3DLAPLACIAN,STATUS = 'KEEP')
        ENDIF
        IF(PLOTSPINDENSITY) THEN
          CLOSE(LU3DSPINDENSITY,STATUS = 'KEEP')
        ENDIF
        IF(PLOTCURRENTDENSITY) THEN
          CLOSE(LU3DCURRENTDENSITY,STATUS = 'KEEP')
        ENDIF
        IF(PLOTRXJ) THEN
          CLOSE(LU3DRXJ,STATUS = 'KEEP')
        ENDIF
      ENDIF
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck getquantity */
      SUBROUTINE GETQUANTITY(DMAT_AO,
     &                       GAO,GAO1,GAO2,
     &                       PX,PY,PZ,
     &                       NCNT,BUF,IPRINT,
     &                       DENSITY,
     &                       SPINDENSITY_NORM,
     &                       CURRENTDENSITY_NORM,
     &                       ELF,
     &                       DLAPLACIAN,
     &                       RXJ_NORM,
     &                       WORK,KFREE,LFREE)
C=======================================================================
C     description ...
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (D23 = 2.0D0/3.0D0,D53 = 5.0D0/3.0D0)
      PARAMETER (D32 = 3.0D0/2.0D0,D54 = 5.0D0/4.0D0)
C-----------------------------------------------------------------------
#include "visual.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbxrs.h"
C-----------------------------------------------------------------------
      DIMENSION DMAT_AO(*),
     &          GAO(*),GAO1(*),GAO2(*),
     &          NCNT(*),BUF(*),
     &          DENSITY(0:2),
     &          GRADIENT(3),
     &          SPINDENSITY(3),
     &          CURRENTDENSITY(3),
     &          RXJ(3),
     &          WORK(*)
C-----------------------------------------------------------------------
C
C
C     pre-stuff
C     =========
C
      IF(PLOTLAPLACIAN) THEN
        CALL GETSOS(GAO,GAO1,GAO2,TEMP,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),2,.FALSE.,IPRINT)
      ELSE
        IF(PLOTELF) THEN
          CALL GETSOS(GAO,GAO1,TEMP,TEMP,NCNT,PX,PY,PZ,
     &                BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,IPRINT)
        ELSE
          CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
     &                BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,IPRINT)
        ENDIF
      ENDIF
C
C
C     density
C     =======
C
      DENSITY(0) = D0
      DO IC = 1,MC
        DENSITY(IC) = D0
        CALL GETRH1(DENSITY(IC),IC,IREP_PLOT,GAO,DMAT_AO,BUF)
        DENSITY(0) = DENSITY(0) + DENSITY(IC)
      ENDDO
C
C
C     gradient
C     ========
C
      CALL GETGRHO(GRADIENT,IREP_PLOT,GAO,GAO1,DMAT_AO,BUF)
      GRADSQUARED = GRADIENT(1)**D2+GRADIENT(2)**D2+GRADIENT(3)**D2
      GRADNORM = DSQRT(GRADSQUARED)
C
C
C     elf
C     ===
C
      IF(PLOTELF) THEN
C
        ELFKIN = D0
C
        DO I = 1,3
          DO IC = 1,MC
            SOMETHING = D0
            CALL GETRH1(SOMETHING,IC,IREP_PLOT,GAO1((I-1)*NTBAS(0)+1),
     &                  DMAT_AO,BUF)
            ELFKIN = ELFKIN+SOMETHING
          ENDDO
        ENDDO
C
        ELFKIN = ELFKIN/D2
C
        PI =  3.1415926536D0
        PRE = 0.3D0*((D3*PI*PI)**D23)
        H = PRE*(DENSITY(0)**D53)
C
        G = ELFKIN-GRADSQUARED/(8.0D0*DENSITY(0))
C
        ELFCHI = G/H
C
C       IF(G.GT.1.0D-8) THEN
          ELF = D1/(D1+ELFCHI*ELFCHI)
C       ELSE
C         ELF = D0
C       ENDIF
C
C       MBE1
C       ELF = ELF**D32
C
C       MBE2
C       ELF = ELF**D54
C
      ENDIF
C
C
C     laplacian
C     =========
C
      IF(PLOTLAPLACIAN) THEN
        CALL QUIT('visualization of the laplacian not implemented')
C       CALL GETLRHO(XX,XY,XZ,YY,YZ,ZZ,
C    &               GAO,GAO1,GAO2,DMAT_AO,BUF)
C       DLAPLACIAN = XX+YY+ZZ
C       DLAPLACIAN = -DLAPLACIAN
      ENDIF
C
C
C     spin and current
C     ================
C
      IF(PLOTSPINDENSITY) THEN
#if defined (MOD_CDFT)
C      CALL GET_CURRENT_DENSITY(.FALSE.,DUMMY,
C    &                        .TRUE.,SPINDENSITY,
C    &                        IREP_PLOT,GAO,DMAT_AO,BUF)
C
       SPINDENSITY_NORM = SPINDENSITY(1)**D2
     &                  + SPINDENSITY(2)**D2
     &                  + SPINDENSITY(3)**D2
       SPINDENSITY_NORM = DSQRT(SPINDENSITY_NORM)
C
       IF(CUT_XY.OR.CUT_XY_SMALL) THEN
         WRITE(0,'(4F30.20)') PX, PY, (SPINDENSITY(I),I=1,2)
       ENDIF
       IF(CUT_YZ.OR.CUT_YZ_SMALL) THEN
         WRITE(0,'(4F30.20)') PY, PZ, (SPINDENSITY(I),I=2,3)
       ENDIF
#else
      CALL QUIT('perturbed magnetization density not implemented')
#endif
      ENDIF
C
      IF(PLOTCURRENTDENSITY) THEN
#if defined (MOD_CDFT)
C      CALL GET_CURRENT_DENSITY(.TRUE.,CURRENTDENSITY,
C    &                        .FALSE.,DUMMY,
C    &                        IREP_PLOT,GAO,DMAT_AO,BUF)
C
       CURRENTDENSITY_NORM = CURRENTDENSITY(1)**D2
     &                     + CURRENTDENSITY(2)**D2
     &                     + CURRENTDENSITY(3)**D2
       CURRENTDENSITY_NORM = DSQRT(CURRENTDENSITY_NORM)
C
       IF(CUT_XY.OR.CUT_XY_SMALL) THEN
         WRITE(0,'(4F30.20)') PX, PY, (CURRENTDENSITY(I),I=1,2)
       ENDIF
       IF(CUT_YZ.OR.CUT_YZ_SMALL) THEN
         WRITE(0,'(4F30.20)') PY, PZ, (CURRENTDENSITY(I),I=2,3)
       ENDIF
#else
      CALL QUIT('perturbed current density not implemented')
#endif
      ENDIF
C
      IF(PLOTRXJ) THEN
#if defined (MOD_CDFT)
C      CALL GET_CURRENT_DENSITY(.TRUE.,CURRENTDENSITY,
C    &                        .FALSE.,DUMMY,
C    &                        IREP_PLOT,GAO,DMAT_AO,BUF)
C
       RXJ(1) = PY*CURRENTDENSITY(3) - PZ*CURRENTDENSITY(2)
       RXJ(2) = PZ*CURRENTDENSITY(1) - PX*CURRENTDENSITY(3)
       RXJ(3) = PX*CURRENTDENSITY(2) - PY*CURRENTDENSITY(1)
C
       RXJ_NORM = RXJ(1)**D2
     &          + RXJ(2)**D2
     &          + RXJ(3)**D2
       RXJ_NORM = DSQRT(RXJ_NORM)
C
       IF(CUT_XY.OR.CUT_XY_SMALL) THEN
         WRITE(0,'(4F30.20)') PX, PY, (RXJ(I),I=1,2)
       ENDIF
       IF(CUT_YZ.OR.CUT_YZ_SMALL) THEN
         WRITE(0,'(4F30.20)') PY, PZ, (RXJ(I),I=2,3)
       ENDIF
#else
      CALL QUIT('.RXJ not implemented')
#endif
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck modified_dmat */
      SUBROUTINE MODIFIED_DMAT_AO(CMO,DMAT_AO,TMATBF,
     &                         STRING1,STRING2,
     &                         WORK,LWORK)
C=======================================================================
C     fill density matrix
C     with components of the solution vector
C     these are read in from XVECFILE
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcborb.h"
#include "visual.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbxrs.h"
C-----------------------------------------------------------------------
      LOGICAL TOBE
C-----------------------------------------------------------------------
      DIMENSION CMO(*),DMAT_AO(*),TMATBF(*),
     &          STRING1(*),STRING2(*),
     &          WORK(LWORK)
C-----------------------------------------------------------------------
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16,PBUF*1
C-----------------------------------------------------------------------
      CALL QENTER('MODIFIED_DMAT_AO')
#include "memint.h"
C-----------------------------------------------------------------------
C
C
C     check if solution vector file is there and open it
C     ==================================================
C
      INQUIRE(FILE = XVECFILE,EXIST = TOBE)
      IF(.NOT.TOBE) THEN
        WRITE(LUPRI,'(A)') 'ERROR in PLOT:'
        WRITE(LUPRI,'(A)') 'solution vector file not found'
        CALL QUIT('ERROR in PLOT')
      ENDIF
      LUXVEC = 17
      OPEN(LUXVEC,FILE = XVECFILE,STATUS = 'OLD',
     &            ACCESS = 'SEQUENTIAL',FORM = 'UNFORMATTED')
C
C
C     get NBDIMs of record NXVEC
C     ==========================
C
      NBDIMEE  = 0
      NBDIMEP  = 0
      NREADS   = 0
      NREADSEE = 0
      NREADSEP = 0
      REWIND LUXVEC
 1    CONTINUE
      READ(LUXVEC) LAB1
      BACKSPACE(LUXVEC)
      IF(LAB1.EQ.'END_OF_THIS_FILE') THEN
        GOTO 10
      ENDIF
      READ(LUXVEC) LAB1, LAB2, TYP1
      BACKSPACE(LUXVEC)
C
      IF(TYP1(1:1).EQ.'E') THEN
        NREADS   = NREADS+1
        NREADSEE = NREADSEE+1
        IF(TYP1(2:2).EQ.'+') THEN
          IHEE = 1
        ELSEIF(TYP1(2:2).EQ.'-') THEN
          IHEE = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUXVEC) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &               JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &               DUM, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NREADSEE.EQ.NXVEC) THEN
          JSYMOPEE = JSYMOP1
          NBDIMEE = LEN/NZ
        ENDIF
        READ(LUXVEC)
C
      ELSEIF(TYP1(1:1).EQ.'P') THEN
        NREADS   = NREADS+1
        NREADSEP = NREADSEP+1
        IF(TYP1(2:2).EQ.'+') THEN
          IHEP = 1
        ELSEIF(TYP1(2:2).EQ.'-') THEN
          IHEP = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUXVEC) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &               JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &               DUM, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NREADSEP.EQ.NXVEC) THEN
          JSYMOPEP = JSYMOP1
          NBDIMEP = LEN/NZ
        ENDIF
        READ(LUXVEC)
C
      ELSE
        GOTO 12
      ENDIF
      GOTO 1
 10   CONTINUE
C
C
C     consistency
C     ===========
C
      IF(NREADS.EQ.0) GOTO 11
      IF(.NOT.(PLOTSKIPEE).AND.(.NOT.PLOTSKIPEP)) THEN
        IF(.NOT.(NREADSEE.EQ.NREADSEP)) GOTO 11
      ENDIF
      IF((.NOT.PLOTSKIPEP).AND.(NREADSEP.EQ.0)) GOTO 13
      IF((.NOT.PLOTSKIPEE).AND.(NREADSEE.EQ.0)) GOTO 14
      IF((.NOT.PLOTSKIPEP).AND.(NREADSEP.LT.NXVEC)) GOTO 11
      IF((.NOT.PLOTSKIPEE).AND.(NREADSEE.LT.NXVEC)) GOTO 11
      IF(PLOTSKIPEE.AND.PLOTSKIPEP) GOTO 15
C
C
C     error handling
C     ==============
C
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
C
C
C     memory allocation
C     =================
C
      CALL MEMGET('REAL',KBVECEE   ,NBDIMEE*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVECEP   ,NBDIMEP*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT     ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
C
C
C     now read the vector
C     ===================
C
      REWIND LUXVEC
      DO I = 1,NXVEC
        IF(.NOT.PLOTSKIPEE) THEN
          READ(LUXVEC)
          IF(I.EQ.NXVEC) THEN
            CALL DZERO(WORK(KBVECEE),NBDIMEE*NZ)
            CALL READT(LUXVEC,NBDIMEE*NZ,WORK(KBVECEE))
          ELSE
            READ(LUXVEC)
          ENDIF
        ELSEIF(NREADSEE.GT.0) THEN
          READ(LUXVEC)
          READ(LUXVEC)
        ENDIF
        IF(.NOT.PLOTSKIPEP) THEN
          READ(LUXVEC)
          IF(I.EQ.NXVEC) THEN
            CALL DZERO(WORK(KBVECEP),NBDIMEP*NZ)
            CALL READT(LUXVEC,NBDIMEP*NZ,WORK(KBVECEP))
          ELSE
            READ(LUXVEC)
          ENDIF
        ELSEIF(NREADSEP.GT.0) THEN
          READ(LUXVEC)
          READ(LUXVEC)
        ENDIF
      ENDDO
C
C
C     close solution vector file
C     ==========================
C
      CLOSE(LUXVEC,STATUS = 'KEEP')
C
C
C     scatter vector
C     ==============
C
      IHERMITICITY = IHEE
      ICOUNTEREE   = 1
      ICOUNTEREP   = 1
      CALL DZERO(WORK(KBMAT),N2ORBXQ)
C
      DO IZ = 1,NZ
        IF(NFSYM.EQ.1) THEN
          CALL SCATTER_BVEC(WORK(KBVECEE),WORK(KBVECEP),
     &                      ICOUNTEREE,ICOUNTEREP,
     &                      WORK(KBMAT),
     &                      IZ,
     &                      1,1)
        ELSE
          IF(IHERMITICITY.EQ.1) THEN
            CALL SCATTER_BVEC(WORK(KBVECEE),WORK(KBVECEP),
     &                        ICOUNTEREE,ICOUNTEREP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,2)
            CALL SCATTER_BVEC(WORK(KBVECEE),WORK(KBVECEP),
     &                        ICOUNTEREE,ICOUNTEREP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,1)
          ELSE
            CALL SCATTER_BVEC(WORK(KBVECEE),WORK(KBVECEP),
     &                        ICOUNTEREE,ICOUNTEREP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,1)
            CALL SCATTER_BVEC(WORK(KBVECEE),WORK(KBVECEP),
     &                        ICOUNTEREE,ICOUNTEREP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,2)
          ENDIF
        ENDIF
      ENDDO
      CALL DSCAL(N2ORBXQ,-D1,WORK(KBMAT),1)
C
C
C     symmetrize
C     ==========
C
      JSYMOP    = JSYMOPEE
      IREP_PLOT = JSYMOP - 1
      JOPSY     = JBTOF(IREP_PLOT,1)
C
      DO IZ = 1,NZ
        IOFF = NORBT*NORBT*(IZ-1)
        CALL FULMAT('S',NORBT,NORBT,WORK(KBMAT+IOFF))
C       CALL FULMAT('A',NORBT,NORBT,WORK(KBMAT+IOFF))
      ENDDO
C
C
C     get density matrix
C     ==================
C
      CALL DZERO(DMAT_AO,N2BBASXQ)
      IF(.NOT.PLOTSKIPEE) THEN
        CALL DMOLI3(0,IHEE,DMAT_AO,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
      ENDIF
      IF(.NOT.PLOTSKIPEP) THEN
        CALL DZERO(TMATBF,N2BBASXQ)
        CALL DMOLI3(1,IHEP,TMATBF,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
        CALL DAXPY(N2BBASXQ,D1,TMATBF,1,DMAT_AO,1)
      ENDIF
C
      IF(NZ.LT.4) THEN
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREP_PLOT)
          IOFF = NTBAS(0)*NTBAS(0)*(IZ-1)
          CALL Q2BPHASE('D',IQ,1,DMAT_AO(1+IOFF))
        ENDDO
      ENDIF
C
C
C     release memory
C     ==============
C
      CALL MEMREL('MODIFIED_DMAT_AO',WORK,KBVECEE,KBVECEE,KFREE,LFREE)
      CALL QEXIT('MODIFIED_DMAT_AO')
C-----------------------------------------------------------------------
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck scatter_bvec */
      SUBROUTINE SCATTER_BVEC(BVECEE,BVECEP,
     &                        ICOUNTEREE,ICOUNTEREP,
     &                        BMAT,
     &                        IZ,
     &                        IFSYM1,IFSYM2)
C=======================================================================
C     scatter solution vector onto BMAT
C-----------------------------------------------------------------------
C     radovan bast                              last revision: july 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C-----------------------------------------------------------------------
#include "dcborb.h"
#include "visual.h"
#include "dgroup.h"
#include "dcbbas.h"
C-----------------------------------------------------------------------
      DIMENSION BVECEE(*),BVECEP(*),
     &          BMAT(NORBT,NORBT,NZ),
     &          NEI(2),NES(2),NP(2)
C-----------------------------------------------------------------------
      DO I = 1,2
        NEI(I) = NISH(I)
        NES(I) = NESH(I) - NISH(I)
        NP(I)  = NPSH(I)
      ENDDO
C
      II1 = 1
      II1 = II1 + (NESH(1)+NPSH(1))*(IFSYM1-1)
      II1 = II1 + NP(IFSYM1)
      II2 = II1 + NEI(IFSYM1) - 1
      IS1 = 1
      IS1 = IS1 + (NESH(1)+NPSH(1))*(IFSYM2-1)
      IS2 = IS1 + NP(IFSYM2)  - 1
      IS3 = IS2 + NEI(IFSYM2) + 1
      IS4 = IS3 + NES(IFSYM2) - 1
C
      DO II = II1,II2
        IF(.NOT.PLOTSKIPEP) THEN
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVECEP(ICOUNTEREP)
            ICOUNTEREP = ICOUNTEREP + 1
          ENDDO
        ENDIF
        IF(.NOT.PLOTSKIPEE) THEN
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVECEE(ICOUNTEREE)
            ICOUNTEREE = ICOUNTEREE + 1
          ENDDO
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck DMAT_MO_2_AO */
      SUBROUTINE DMAT_MO_2_AO(CMO,DMAT_AO,DMAT_MO,WORK,LWORK)
C=======================================================================
C     transform MO to AO density matrix
C-----------------------------------------------------------------------
C     radovan bast                           last revision: october 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
!-----------------------------------------------------------------------
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION CMO(*),DMAT_AO(*),DMAT_MO(*),WORK(LWORK)

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
C     radovan bast                           last revision: october 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
!-----------------------------------------------------------------------
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION DMAT_MO(*),STRING1(*),STRING2(*)

      CALL DZERO(DMAT_MO,N2ORBXQ)

      K = 0
      DO I = 1,NFSYM
        IF(NISH(I) .GT. 0) THEN
          DO J = 1,NISH(I)
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
