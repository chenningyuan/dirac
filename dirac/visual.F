C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck plot_mem */
      SUBROUTINE PLOT_MEM(WORK,LWORK)
C  .....................................................................
C  .  description:                                                     .
C  .  prepare {1,2,3}-dimensional plot                                 .
C  .  copy-and-paste-and-modify-a-little-job                           .
C  .  used Trond's respective routines                                 .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dcbplot.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"
C
      LOGICAL TOBE
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('PLOT_MEM')
C
#include "memint.h"
C
      CALL TITLER('{1,2,3}-dimensional plot','*',122)
C
C
C     check if coefficients are on file
C     =================================
C
      INQUIRE(FILE = 'DFCOEF',EXIST = TOBE)
      IF(.NOT.TOBE) THEN
        WRITE(LUPRI,'(/A)') 'PLOT_MEM: Coefficients not found!'
        CALL QUIT('PLOT_MEM: Coefficients not found!')
      ENDIF
C
C
C     memory allocation
C     =================
C
      CALL MEMGET('REAL',KCMO       ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATBF    ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMATORB   ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORBSTRING1,NISH(1),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KORBSTRING2,NISH(2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO       ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO1      ,3*NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO2      ,6*NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCNT      ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOOR      ,6,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF       ,NTBAS(0)*4,
     &            WORK,KFREE,LFREE)
C
C
C     read coefficients and give info
C     ===============================
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PLOT_MEM')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),TEMP,TEMP,
     &            TOTERG,3,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')
C
C
C     zero out orbitals?
C     ==================
C
      CALL DFILL(WORK(KORBSTRING1),D1,NISH(1))
      IF(NFSYM.GT.1) CALL DFILL(WORK(KORBSTRING2),D1,NISH(2))
      IF(USEORBSTRING) THEN
        IF(ORBSTRING1(4:5).EQ.'..') THEN
          READ(ORBSTRING1(1:3),'(I3)') IORBSTRING1BEGIN
          READ(ORBSTRING1(6:8),'(I3)') IORBSTRING1END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(3:4).EQ.'..') THEN
          READ(ORBSTRING1(1:2),'(I3)') IORBSTRING1BEGIN
          READ(ORBSTRING1(5:8),'(I3)') IORBSTRING1END
          GOTO 1
        ENDIF
        IF(ORBSTRING1(2:3).EQ.'..') THEN
          READ(ORBSTRING1(1:1),'(I3)') IORBSTRING1BEGIN
          READ(ORBSTRING1(4:8),'(I3)') IORBSTRING1END
          GOTO 1
        ENDIF
        WRITE(LUPRI,'(/A)') 'ERROR: wrong orbital string format!'
        GOTO 2
 1      CONTINUE
C
        IF(
     &         (.NOT.(IORBSTRING1BEGIN.LE.IORBSTRING1END))
     &     .OR.(.NOT.(IORBSTRING1BEGIN.GE.0))
     &     .OR.(.NOT.(IORBSTRING1END.LE.NISH(1)))) THEN
          WRITE(LUPRI,'(/A)') 'ERROR: wrong orbital string dimension!'
          GOTO 2
        ENDIF
        CALL DFILL(WORK(KORBSTRING1),D0,NISH(1))
        DO I = IORBSTRING1BEGIN,IORBSTRING1END
          WORK(KORBSTRING1-1+I) = D1
        ENDDO
 2      CONTINUE
        
        IF(NFSYM.GT.1) THEN
          IF(ORBSTRING2(4:5).EQ.'..') THEN
            READ(ORBSTRING2(1:3),'(I3)') IORBSTRING2BEGIN
            READ(ORBSTRING2(6:8),'(I3)') IORBSTRING2END
            GOTO 3
          ENDIF
          IF(ORBSTRING2(3:4).EQ.'..') THEN
            READ(ORBSTRING2(1:2),'(I3)') IORBSTRING2BEGIN
            READ(ORBSTRING2(5:8),'(I3)') IORBSTRING2END
            GOTO 3
          ENDIF
          IF(ORBSTRING2(2:3).EQ.'..') THEN
            READ(ORBSTRING2(1:1),'(I3)') IORBSTRING2BEGIN
            READ(ORBSTRING2(4:8),'(I3)') IORBSTRING2END
            GOTO 3
          ENDIF
          WRITE(LUPRI,'(/A)') 'ERROR: wrong orbital string format!'
          GOTO 4
 3        CONTINUE
C
          IF(
     &           (.NOT.(IORBSTRING2BEGIN.LE.IORBSTRING2END))
     &       .OR.(.NOT.(IORBSTRING2BEGIN.GE.0))
     &       .OR.(.NOT.(IORBSTRING2END.LE.NISH(2)))) THEN
            WRITE(LUPRI,'(/A)') 'ERROR: wrong orbital string dimension!'
            GOTO 4
          ENDIF
          CALL DFILL(WORK(KORBSTRING2),D0,NISH(2))
          DO I = IORBSTRING2BEGIN,IORBSTRING2END
            WORK(KORBSTRING2-1+I) = D1
          ENDDO
 4        CONTINUE
        ENDIF
      ENDIF
C
C
C     generate density matrix
C     =======================
C
      CALL GENDEN(WORK(KDMATBF),WORK(KCMO),1,IPLOT)
C
C     blocked
C     CALL PERDMAT(WORK(KCMO),WORK(KDMATBF),WORK(KFREE),LFREE)
C
      CALL OCCUPYDMATORB(WORK(KDMATORB),
     &                   WORK(KORBSTRING1),WORK(KORBSTRING2))
      CALL DMATORB2BF(WORK(KCMO),WORK(KDMATBF),WORK(KDMATORB),
     &                WORK(KFREE),LFREE)
      CALL DSCAL(N2BBASXQ,D2,WORK(KDMATBF),1)
C
C
C     call PLOT routine
C     =================
C
      CALL PLOT(WORK(KDMATBF),
     &          WORK(KGAO),WORK(KGAO1),WORK(KGAO2),
     &          WORK(KCMO),WORK(KNCNT),WORK(KCOOR),WORK(KBUF),
     &          WORK,KFREE,LFREE)
C
C
C     memory deallocation
C     ===================
C
      CALL MEMREL('PLOT_MEM',WORK,KCMO,KCMO,KFREE,LFREE)
      CALL QEXIT('PLOT_MEM')
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck plot */
      SUBROUTINE PLOT(DMATBF,
     &                GAO,GAO1,GAO2,
     &                CMO,NCNT,COOR,BUF,
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
#include "dcbplot.h"
#include "dcborb.h"
C
      DIMENSION DMATBF(N2BBASX,NZ),
     &          GAO(*),GAO1(*),GAO2(*),
     &          CMO(*),NCNT(*),COORD(3,2),BUF(*),
     &          WORK(*)
C
Cvorsicht       DENSITY(0:2) war frueher im 1D
      DIMENSION DENSITY(0:2,NCUBE(3)),
     &          SPINDENSITY(NCUBE(3)),
     &          CURRENTDENSITY(NCUBE(3)),
     &          ELF(NCUBE(3)),
     &          DLAPLACIAN(NCUBE(3)),
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
            CALL DOSTUFF(DMATBF,GAO,GAO1,GAO2,CMO,PX,PY,PZ,
     &                   NCNT,COORD,BUF,IPRINT,
     &                   DENSITY,
     &                   SPINDENSITY(1),
     &                   CURRENTDENSITY(1),
     &                   ELF(1),
     &                   DLAPLACIAN(1),
     &                   WORK,KFREE,LFREE)
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
            WRITE(LU3DLAPLACIAN,'(I5,4F12.6)') NCUBE(I),
     &                                              (D3DSTEP(J,I),J=1,3)
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
            WRITE(LU3DSPINDENSITY,'(I5,4F12.6)') NCUBE(I),
     &                                              (D3DSTEP(J,I),J=1,3)
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
            WRITE(LU3DCURRENTDENSITY,'(I5,4F12.6)') NCUBE(I),
     &                                              (D3DSTEP(J,I),J=1,3)
          ENDDO      
          DO I = 1,NUCDEP
            WRITE(LU3DCURRENTDENSITY,'(I5,4F12.6)')
     &            NINT(ATOMXYZ(4,I)), ATOMXYZ(4,I), (ATOMXYZ(J,I),J=1,3)
          ENDDO
        ENDIF
C
C
C       do the grid walk
C       ================
C
C     { 3D
        DO ICUBE = 1,NCUBE(1)
          DO JCUBE = 1,NCUBE(2)
            DO KCUBE = 1,NCUBE(3)
              PX = D3DSTEP(1,0)+(ICUBE-1)*D3DSTEP(1,1)
     &                         +(JCUBE-1)*D3DSTEP(1,2)
     &                         +(KCUBE-1)*D3DSTEP(1,3)
              PY = D3DSTEP(2,0)+(ICUBE-1)*D3DSTEP(2,1)
     &                         +(JCUBE-1)*D3DSTEP(2,2)
     &                         +(KCUBE-1)*D3DSTEP(2,3)
              PZ = D3DSTEP(3,0)+(ICUBE-1)*D3DSTEP(3,1)
     &                         +(JCUBE-1)*D3DSTEP(3,2)
     &                         +(KCUBE-1)*D3DSTEP(3,3)
C       3D }
C     { xy - plane
C       DO ICUBE = 1,NCUBE(1)
C         DO JCUBE = 1,NCUBE(2)
C           NCUBE(3) = 1
C           DO KCUBE = 1,NCUBE(3)
C             PX = -6.0D0+(ICUBE-1)*(12.0D0/(NCUBE(1)-1))
C             PY = -6.0D0+(JCUBE-1)*(12.0D0/(NCUBE(2)-1))
C             PZ = 0.2D0
C       xy - plane }
C     { yz - plane
C       NCUBE(1) = 1
C       DO ICUBE = 1,NCUBE(1)
C         DO JCUBE = 1,NCUBE(2)
C           DO KCUBE = 1,NCUBE(3)
C             PX = D0
C             PY = -3.0D0+(JCUBE-1)*(6.0D0/(NCUBE(2)-1))
C             PZ = -3.0D0+(KCUBE-1)*(6.0D0/(NCUBE(3)-1))
C       yz - plane }
C
C
C             get necessary stuff
C             ===================
C
              CALL DOSTUFF(DMATBF,GAO,GAO1,GAO2,CMO,PX,PY,PZ,
     &                     NCNT,COORD,BUF,IPRINT,
     &                     DENSITY(0,KCUBE),
     &                     SPINDENSITY(KCUBE),
     &                     CURRENTDENSITY(KCUBE),
     &                     ELF(KCUBE),
     &                     DLAPLACIAN(KCUBE),
     &                     WORK,KFREE,LFREE)
C
C
            ENDDO
C
C
C           write buffer
C           ============
C
            IF(PLOTDENSITY) THEN
              WRITE(LU3DDENSITY,'(6E13.5)')
     &                                       (DENSITY(0,M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTELF) THEN
              WRITE(LU3DELF,'(6E13.5)') (ELF(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTLAPLACIAN) THEN
              WRITE(LU3DLAPLACIAN,'(6E13.5)')
     &                                      (DLAPLACIAN(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTSPINDENSITY) THEN
              WRITE(LU3DSPINDENSITY,'(6E13.5)')
     &                                     (SPINDENSITY(M),M=1,NCUBE(3))
            ENDIF
            IF(PLOTCURRENTDENSITY) THEN
              WRITE(LU3DCURRENTDENSITY,'(6E13.5)')
     &                                  (CURRENTDENSITY(M),M=1,NCUBE(3))
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
      ENDIF
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dostuff */
      SUBROUTINE DOSTUFF(DMATBF,GAO,GAO1,GAO2,CMO,PX,PY,PZ,
     &                   NCNT,COORD,BUF,IPRINT,
     &                   DENSITY,
     &                   SPINDENSITY,
     &                   CURRENTDENSITY,
     &                   ELF,
     &                   DLAPLACIAN,
     &                   WORK,KFREE,LFREE)
C  .....................................................................
C  .  description:                                                     .
C  .                                                                   .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (D23 = 2.0D0/3.0D0,D53 = 5.0D0/3.0D0)
      PARAMETER (D32 = 3.0D0/2.0D0,D54 = 5.0D0/4.0D0)
C
#include "dcbplot.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbxrs.h"
C
      DIMENSION DMATBF(N2BBASX,NZ),
     &          GAO(*),GAO1(*),GAO2(*),CMO(*),
     &          NCNT(*),COORD(3,2),BUF(*),
     &          WORK(*)
C
      DIMENSION DENSITY(0:2),GRADIENT(3),
     &          CURRENT(4),SPIN(4)
C
C
C     pre-stuff
C     =========
C
      IF(PLOTLAPLACIAN) THEN
        CALL GETSOS(GAO,GAO1,GAO2,TEMP,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),2,.FALSE.,IPRINT)
      ELSE
        CALL GETSOS(GAO,GAO1,TEMP,TEMP,NCNT,PX,PY,PZ,
     &              BUF,WORK(KFREE),LFREE,NTBAS(0),1,.FALSE.,IPRINT)
C       CALL GETSOS(GAO,TEMP,TEMP,TEMP,TEMP,PX,PY,PZ,
C    &              BUF,WORK(KFREE),LFREE,NTBAS(0),0,.FALSE.,IPRINT)
      ENDIF
C
C
C     density
C     =======
C
      DENSITY(0) = D0
      DO IC = 1,MC
        DENSITY(IC) = D0
        CALL GETRH1(DENSITY(IC),IC,0,GAO,DMATBF,BUF)
        DENSITY(0) = DENSITY(0)+DENSITY(IC)
      ENDDO
C
C
C     gradient
C     ========
C
      CALL GETGRHO(GRADIENT,0,GAO,GAO1,DMATBF,BUF)
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
            CALL GETRH1(SOMETHING,IC,0,GAO1((I-1)*NTBAS(0)+1),
     &                  DMATBF,BUF)
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
C     blocked
C     IF(PLOTLAPLACIAN) THEN
C       CALL GETLRHO(XX,XY,XZ,YY,YZ,ZZ,
C    &               GAO,GAO1,GAO2,DMATBF,BUF)
CC      BLUB       = GRADIENT(1)*XX*GRADIENT(1)
CC   &             + GRADIENT(2)*YY*GRADIENT(2)
CC   &             + GRADIENT(3)*ZZ*GRADIENT(3)
CC   &             + GRADIENT(1)*XY*GRADIENT(2)
CC   &             + GRADIENT(1)*XZ*GRADIENT(3)
CC   &             + GRADIENT(2)*YZ*GRADIENT(3)
C       DLAPLACIAN = XX+YY+ZZ
C       DLAPLACIAN = -DLAPLACIAN
C     ENDIF
C
C
C     spin and current
C     ================
C
      IF(PLOTSPINDENSITY) THEN
        CALL GETCURRENTDENSITY(.FALSE.,DUMMY,.TRUE.,SPIN,
     &                         0,GAO,DMATBF,BUF)
        SPINDENSITY = SPIN(4)
C       SPINDENSITY = SPIN(3)
C       IF(ABS(SPIN(4)).GT.0.0001D0) THEN
C         WRITE(0,'(4F20.6)') PX, PY, (SPIN(I),I=1,2)
C         WRITE(0,'(4F20.6)') PY, PZ, (SPIN(I),I=2,3)
C       ENDIF
      ENDIF
      IF(PLOTCURRENTDENSITY) THEN
        CALL GETCURRENTDENSITY(.TRUE.,CURRENT,.FALSE.,DUMMY,
     &                         0,GAO,DMATBF,BUF)
        CURRENTDENSITY = CURRENT(4)
C       CURRENTDENSITY = CURRENT(1)
C       IF(ABS(CURRENT(4)).GT.0.0001D0) THEN
C         WRITE(0,'(4F20.6)') PX, PY, (CURRENT(I),I=1,2)
C         WRITE(0,'(4F20.6)') PY, PZ, (CURRENT(I),I=2,3)
C       ENDIF
      ENDIF
C
C
C     slice
C     =====
C
      IF(PZ.GT.D0) THEN
C       DENSITY = D0
C       SPINDENSITY = D0
C       CURRENTDENSITY = D0
C       ELF = D0
C       DLAPLACIAN = D0
      ENDIF
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck perdmat */
      SUBROUTINE PERDMAT(CMO,DMATBF,WORK,LWORK)
C  .....................................................................
C  .  description:                                                     .
C  .  fill density matrix                                              .
C  .  with components of the solution vector                           .
C  .  these are read in from XVECFILE                                  .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
#include "priunit.h"
C
#include "dcbqal.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "dcborb.h"
#include "dcbplot.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "dcbxrs.h"
C
      LOGICAL TOBE
C
      DIMENSION CMO(*),DMATBF(*),WORK(LWORK)
C
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
C
C
      CALL QENTER('PERDMAT')
C
#include "memint.h"
C
C
C     init
C     ====
C
      NIA    = NESH(1)+NESH(2)
      NI     = NISH(1)+NISH(2)
      NP     = NORBT-NIA
      NA     = NIA-NI
      NBDIM  = NI*NA
      JSYMOP = 1
      JOPSY  = 1
C
C
C     memory allocation
C     =================
C
      CALL MEMGET('REAL',KBVEC     ,NBDIM*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT     ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
C
C
C     check if solution vector file is there
C     ======================================
C
      INQUIRE(FILE = 'XVECFILE',EXIST = TOBE)
      IF(.NOT.TOBE) THEN
C       panic
      ENDIF
C
C
C     open solution vector file
C     =========================
C
      LUXVEC = 17
      OPEN(LUXVEC,FILE = 'XVECFILE',STATUS = 'OLD',
     &            ACCESS = 'SEQUENTIAL',FORM = 'UNFORMATTED')
C
C
C     read vectors
C     ============
C
      REWIND LUXVEC
C
      DO I = 1,NXVEC
C
        READ(LUXVEC) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &               JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &               DUM, LEN, INTFLG, ERGRSP, NBAS, NORBS
C
        CALL DZERO(WORK(KBVEC),NBDIM*NZ)
        CALL READT(LUXVEC,LEN,WORK(KBVEC))
C
        IREP = JSYMOP1-1
C
      ENDDO
C     WRITE(LUPRI,*) 'inside PERDMAT: BVEC'
      CALL PRQMAT(WORK(KBVEC),NBDIM,1,NBDIM,1,
     &            NZ,IPQTOQ(1,IREP),LUPRI)
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
      CALL DZERO(WORK(KBMAT),N2ORBXQ)
      DO IZ = 1,NZ
        DO II = 1,NI
          DO IA = 1,NA
            IOFF = 0
            IOFF = IOFF+NORBT*NORBT*(IZ-1)
            IOFF = IOFF+NORBT*NP
            IOFF = IOFF+NORBT*(II-1)
            IOFF = IOFF+NP
            IOFF = IOFF+NI
            IOFF = IOFF+IA
            JOFF = 0
            JOFF = JOFF+NBDIM*(IZ-1)
            JOFF = JOFF+(II-1)*NA+IA
            WORK(KBMAT+IOFF) = WORK(KBVEC+JOFF)
          ENDDO
        ENDDO
      ENDDO
C
C
C     symmetrize
C     ==========
C
      DO IZ = 1,NZ
        IOFF = NORBT*NORBT*(IZ-1)
        CALL FULMAT('A',NORBT,NORBT,WORK(KBMAT+IOFF))
      ENDDO
C     WRITE(LUPRI,*) 'inside PERDMAT: BMAT'
C     CALL PRQMAT(WORK(KBMAT),NORBT,NORBT,NORBT,NORBT,
C    &            NZ,IPQTOQ(1,IREP),LUPRI)
C
C
C     get density matrix
C     ==================
C
C     CALL DMOLI3(0*1+0*2,1,DMATBF,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
      CALL DMOLI3(0*1+0*2,-1,DMATBF,CMO,WORK(KBMAT),WORK,KFREE,LFREE)
C
C     CALL MATRIXIO(DMATBF,NTBAS(0),NTBAS(0),NZ,.FALSE.,'dmat')
C     WRITE(LUPRI,*) 'inside PERDMAT: DMATBF'
C     CALL PRQMAT(DMATBF,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
C    &            NZ,IPQTOQ(1,IREP),LUPRI)
C
C
C     release memory
C     ==============
C
      CALL MEMREL('PERDMAT',WORK,KBVEC,KBVEC,KFREE,LFREE)
      CALL QEXIT('PERDMAT')
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck dmatorb2bf */
      SUBROUTINE DMATORB2BF(CMO,DMATBF,DMATORB,WORK,LWORK)
C  .....................................................................
C  .  description:                                                     .
C  .  transform orbital density matrix to bf density matrix            .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
C
      DIMENSION CMO(*),DMATBF(*),DMATORB(*),WORK(LWORK)
C
#include "memint.h"
C
C
      CALL DZERO(DMATBF,N2BBASXQ)
C
C     CALL PRQMAT(DMATORB,NORBT,NORBT,NORBT,NORBT,
C    &            NZ,IPQTOQ(1,0),0)
C
      DO I = 1,NFSYM
        CALL QTRANS('MOAO','S',D0,
C       CALL QTRANS('MOAO','A',D0,
     &              NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &              DMATBF(1+I2BASX(I,I)),
     &              NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &              DMATORB(1+I2ORBX(I,I)),
     &              NORBT,NORBT,1,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,IPRINT)
      ENDDO
C
C     CALL PRQMAT(CMO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
C    &            NZ,IPQTOQ(1,0),0)
C     CALL PRQMAT(DMATBF,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
C    &            NZ,IPQTOQ(1,0),0)
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* deck occupydmatorb */
      SUBROUTINE OCCUPYDMATORB(DMATORB,ORBSTRING1,ORBSTRING2)
C  .....................................................................
C  .  description:                                                     .
C  .  occupy orbital density matrix                                    .
C  .....................................................................
C  .  author:                                          last revision:  .
C  .  radovan bast                                     october 2005    .
C  .....................................................................
C
#include "implicit.h"
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
C
C
      DIMENSION DMATORB(*),ORBSTRING1(*),ORBSTRING2(*)
C
C
      CALL DZERO(DMATORB,N2ORBXQ)
C
      K = 0
      DO I = 1,NFSYM
        IF(NISH(I).GT.0) THEN
          DO J = 1,NISH(I)
            K = K+1
            JOFF = (NORB(1)*(I-1)+NPSH(I))*NORBT
            JOFF = JOFF+NORB(1)*(I-1)+NPSH(I)
            JOFF = JOFF+(J-1)*NORBT+J-1
C
            IF((I.EQ.1).AND.(ORBSTRING1(J).EQ.D1)) THEN
              DMATORB(1+JOFF) = D1
            ENDIF
            IF((I.EQ.2).AND.(ORBSTRING2(J).EQ.D1)) THEN
              DMATORB(1+JOFF) = D1
            ENDIF
C
          ENDDO
        ENDIF
      ENDDO
C
C
      RETURN
      END
