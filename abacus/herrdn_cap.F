!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!dirac_copyright_end

#ifdef MOD_CAP
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck capout */
      SUBROUTINE CAPOUT(IPRINT,
     &  IQM,NBLCK,JCO,NUC,NRC,SEG,
     &  KATOM,KANG,KSETS,KBLOCK,
     &  KPRIM,CPRIMU,WORK,LWORK)
C***********************************************************************
C
C    Dumps files for the CAP integral module.
C
C    Called from BASOUT(herrdn.F)
C
C    Written by Miro Ilias, January 2007, Tel Aviv
C
C***********************************************************************
      !use cap
#include "implicit.h"
#include "priunit.h"

#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "primit.h"
#include "chrxyz.h"
#include "nuclei.h"
#include "symmet.h"
#include "dummy.h"
#include "infpar.h"

      LOGICAL SEG,HERMIT
      DIMENSION IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &     NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &     SEG(KBLOCK,KSETS),
     &     ALPHA(KPRIM,KBLOCK,KSETS),
     &     CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)

      DIMENSION WORK(LWORK)
C     CHARACTER SPDCAR
      CHARACTER*4 ATOMNAME
#include "ibtfun.h"

C     .... only for master or single process !!! 

      IF (MYTID .EQ. MPARID) THEN

      KIPCON = 1
      KLAST  = KIPCON + KPRIM
C      IF (KLAST.GT.LWORK) CALL STOPIT('BASOUT','CAPOUT',KLAST,LWORK)
      IF (KLAST.GT.LWORK) CALL QUIT('CAPOUT: KLAST.GT.LWORK !')

      IF (IPRINT.GE.2) THEN
       CALL HEADER('Output from CAPOUT ',-1)
       WRITE(LUPRI,'(2X,A)')
     & 'Preparing formatted file "CAPINT.inp" for the CAP module.'
      ENDIF
 
       ! ... for cap module
      !  napg = 0
C!
C!   Prepare another unique file with basis functions for Sajeev's code
C!
      LUTMP2 = -1
      LUTMP  = -1
      CALL GPOPEN(LUTMP2,'CAPINT.inp',
     & 'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)

C     CALL GPOPEN(LUTMP,
C    &  'CAPLBL.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)

C     ... count the total number of nuclei
      NUCTOT = 0
      DO ICENT = 1, NUCIND
       NUCTOT = NUCTOT + NUCDEG(ICENT)
      ENDDO

C      ... write the total number of nuclei
      WRITE(LUTMP2,'(1X,I3)') NUCTOT
      IF (IPRINT.GT.3) THEN
        WRITE(LUPRI,'(2X,A,I3)')
     &  '...total number of nuclei:',NUCTOT
      ENDIF

C      ... cycle over all nuclei ...
      ICRX = 1
      ICRY = 2
      ICRZ = 3

      ICENT = 0

      DO IATYP = 1, NONTYP
      DO JATYP = 1, NONT(IATYP)

        ICENT = ICENT + 1

C     DO ICENT = 1, NUCIND
C ............... basis functions ..................
         IF (IPRINT.GE.3) THEN
          write(LUPRI,'(2X,A,I2,A,I2)')
     &    'ICENT=',ICENT, '  NDEG=',NUCDEG(ICENT)
         ENDIF

         MULCNT = ISTBNU(ICENT)
         ATOMNAME   = NAMEX(3*ICENT)(1:4)
            IF (MULT(MULCNT) .EQ. 1) THEN
              IF (IPRINT.GE.3) THEN
               WRITE (LUPRI,'(I4,3X,A,5X,A,3X,F15.10)')
     *               ICRX, ATOMNAME,CHRXYZ(-1), CORD(1,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *               ICRY, CHRXYZ(-2), CORD(2,ICENT)
               WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *               ICRZ, CHRXYZ(-3), CORD(3,ICENT)
               ENDIF

               ICRX = ICRX + 3
               ICRY = ICRY + 3
               ICRZ = ICRZ + 3

               write(LUTMP2,*)
     &         ICENT,'  ',ATOMNAME,' charge:',CHARGE(ICENT)
               write(LUTMP2,*) CORD(1,ICENT),CORD(2,ICENT),CORD(3,ICENT)

              IF (IPRINT.GE.3) THEN
               write(LUPRI,*)
     &         ICENT,'  ',ATOMNAME,' charge:',CHARGE(ICENT)
               write(LUPRI,*) CORD(1,ICENT),CORD(2,ICENT),CORD(3,ICENT)
              ENDIF

C             ...   write out basis functions, exponents and contraction coefficients
             IF (IPRINT.GE.3) THEN
              WRITE(LUPRI,'(/A)')
     &        'Large component basis functions for CAP'
             ENDIF
             CALL CAPWRIBAS(IQM(1,1),NBLCK(1,1),
     &           JCO(1,1,1),NUC(1,1),NRC(1,1),SEG(1,1),WORK(KIPCON),
     &         KATOM,KANG,KBLOCK,KPRIM,CPRIMU,
     &          LUTMP2,LUTMP,IATYP,0,IPRINT)

             IF (NBASIS.GT.NLARGE) THEN
               IF (IPRINT.GE.3) THEN
                 WRITE(LUPRI,'(/A)')
     &           'Small component basis functions for CAP'
               ENDIF
               CALL CAPWRIBAS(IQM(1,2),NBLCK(1,2),
     &             JCO(1,1,2),NUC(1,2),NRC(1,2),SEG(1,2),WORK(KIPCON),
     &            KATOM,KANG,KBLOCK,KPRIM,CPRIMU(1,1,1,2),
     &            LUTMP2,LUTMP,IATYP,NPLSH,IPRINT)
             ENDIF
             write(lutmp2,*) 'E   0   0'

            ELSE

               JATOM = 0
               DO ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                     JATOM = JATOM + 1

                    IF (IPRINT.GE.3) THEN
                     WRITE (LUPRI,'(I4,3X,A,I2,3X,A,3X,F15.10)')
     *                     ICRX, ATOMNAME,JATOM,CHRXYZ(-1),
     *                     PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')
     *                     ICRY,CHRXYZ(-2),
     *                     PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                     WRITE (LUPRI,'(I4,12X,A,3X,F15.10,/)')
     *                     ICRZ,CHRXYZ(-3),
     *                     PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                     ENDIF

                     ICRX = ICRX + 3
                     ICRY = ICRY + 3
                     ICRZ = ICRZ + 3

C              write(LUTMP2,*) ICENT,'  ',ATOMNAME,' ',JATOM,
               write(LUTMP2,'(I3,A,A,A,I2,A,F12.6)')
     &         ICENT,'  ',ATOMNAME,' ',JATOM,' charge:',CHARGE(ICENT)
C               ... coordinates
               write(LUTMP2,*)
     &              PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT),
     &              PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT),
     &              PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)

               IF (IPRINT.GE.3) THEN
               write(LUPRI,*) ICENT,'  ',ATOMNAME,' ',JATOM,
     &              ' charge:',CHARGE(ICENT)
C               ... coordinates
               write(LUPRI,*)
     &              PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT),
     &              PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT),
     &              PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               ENDIF

C             ...   write out basis functions, exponents and contraction coefficients
             IF (IPRINT.GE.3) THEN
               write(LUPRI,'(/A)') '--- large component basis set !'
             ENDIF

            CALL CAPWRIBAS(IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &      NUC(1,1),NRC(1,1),SEG(1,1),WORK(KIPCON),
     &      KATOM,KANG,KBLOCK,KPRIM,CPRIMU,
     &      LUTMP2,LUTMP,IATYP,0,IPRINT)

             IF (NBASIS.GT.NLARGE) THEN
               IF (IPRINT.GE.3) THEN
                 write(LUPRI,'(/A)') '--- small component basis set !'
               ENDIF
               CALL CAPWRIBAS(IQM(1,2),NBLCK(1,2),
     &        JCO(1,1,2),NUC(1,2),NRC(1,2),SEG(1,2),WORK(KIPCON),
     &        KATOM,KANG,KBLOCK,KPRIM,CPRIMU(1,1,1,2),
     &        LUTMP2,LUTMP,IATYP,NPLSH,IPRINT)
             ENDIF
             write(lutmp2,*) 'E   0   0'
            END IF
           ENDDO    
         ENDIF

         ENDDO    
         ENDDO    

         CALL GPCLOSE(LUTMP2,'KEEP')

C        
         LUTMP3=-1
         CALL GPOPEN(LUTMP3,
     &  'CAPDIR.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)
         WRITE(LUTMP3,*) 'DIRAC unique order of atomic basis functions'

         IF (IPRINT.GE.0) THEN
          WRITE(LUPRI,'(2X,A)')
     &    'CAPOUT: preparing the "CAPDIR.inp" file with DIRAC'//
     &    ' unique order of atomic basis functions. '
         ENDIF

         IF (IPRINT.GE.2) THEN
           CALL HEADER('CAPWRIBAS2: Contracted Orbitals again',1)
         ENDIF

          CALL CAPWRIBAS2(IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &     NUC(1,1),NRC(1,1),SEG(1,1),WORK(KIPCON),
     &     KATOM,KANG,KBLOCK,KPRIM,CPRIMU,
     &     LUTMP3,0,IPRINT)

         IF (NBASIS.GT.NLARGE) THEN
C          ... small component set
          CALL CAPWRIBAS2(IQM(1,2),NBLCK(1,2),
     &      JCO(1,1,2),NUC(1,2),NRC(1,2),SEG(1,2),WORK(KIPCON),
     &      KATOM,KANG,KBLOCK,KPRIM,CPRIMU(1,1,1,2),
     &      LUTMP3,NPLSH,IPRINT)
         ENDIF
         CALL GPCLOSE(LUTMP3,'KEEP')

CMI      CALL CAPSYMOUT

      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wribas */
      SUBROUTINE CAPWRIBAS(IQM,NBLCK,JCO,NUC,NRC,SEG,IPCON,
     &     KATOM,KANG,KBLOCK,KPRIM,CPRIMU,
     &     LUTMP2,LUTMP,ICENT,IPOSEXP,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Fill the CAPINT.inp file with input for CAP program
C
C   Called from CAPOUT
C
C   Written by Miro Ilias, Tel Aviv, January 2007
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      !use cap

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "primit.h"

C#include "ccom.h"
C#include "mxcent.h"
C#include "nuclei.h"

      LOGICAL SEG, NRMPRI
      CHARACTER*1 BLNK(8)
C     CHARACTER CHRSEG*10,DEG*2

      DIMENSION IQM(KATOM),NBLCK(KATOM),
!    &     JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &     JCO(KANG,KATOM),NUC(*),NRC(*),
!    &     CPRIMU(KPRIM,KPRIM,KBLOCK),SEG(KBLOCK)
     &     CPRIMU(KPRIM,KPRIM,*),SEG(KBLOCK)
      DIMENSION IPCON(KPRIM)

      DATA BLNK /'E','S','P','D','F','G','H','I'/

      CALL QENTER('CAPWRIBAS')

      IF (IPRINT.GE.3) THEN
         write(LUPRI,*) '*** In CAPWRIBAS ***'
         write(LUPRI,*) 'ICENT=  ',ICENT
C       ... run over blocks in NBLCK(I)
         write(LUPRI,*) ' nu of blocks of shells, NBLCK:',NBLCK(ICENT) 
         write(LUPRI,*) ' nu of blocks of shells, IQM:',IQM(ICENT) 
       ENDIF
C       ... start of the block
       IEXP = IPOSEXP
       KBLOCK = 0
       IIQM = 0
       DO I = 1, ICENT-1
         KBLOCK = KBLOCK + NBLCK(I)
C        KBLOCK = KBLOCK + IQM(I)
C        DO II = 1 , IQM(I)
         DO II = 1 , NBLCK(I)
C       ... find out pointer to beginning of exponents
CMI       IF ( (IIQM+II).NE.0 ) THEN ! boundaries fix
CMI       IF ((IIQM+II).LE.KBLOCK .AND. (IIQM+II).NE.0) THEN ! boundaries fix
           IEXP = IEXP + NUC(IIQM+II)
CMI       ENDIF
         ENDDO
C        IIQM = IIQM + IQM(I)
         IIQM = IIQM + NBLCK(I)
       ENDDO
 
       INUC = 0
       IBLOCK = KBLOCK + 1
       DO J = 1, IQM(ICENT)
CC     DO J = 1, NBLCK(ICENT)
C        IBLOCK = KBLOCK + J
C        ... run over subblocks (small comp. set usually contain some)
        DO K = 1, JCO(J,ICENT)
C        IBLOCK = KBLOCK + J
         IF (IPRINT.GE.3) THEN
          write(LUPRI,*) 'IEXP pointer:',IEXP
          write(LUPRI,*) 'INUC pointer:',INUC,' IBLOCK:',IBLOCK
          write(LUPRI,*) 'uncontracted shells: ',NUC(IBLOCK),
     &    ' --->  contracted shells:',NRC(IBLOCK)
C         write(LUTMP2,*) 'TYPE:',BLNK(1+J)
          write(LUPRI,*) 'TYPE:',BLNK(1+J)
          write(LUPRI,*) NUC(IBLOCK),' exponents and contract.coeff.:'
         ENDIF

          DO L_CON = 1, NRC(IBLOCK) 
            IBASCOUNT = 0
            DO L_UNC =1 , NUC(IBLOCK)
              CCOEF = CPRIMU(L_UNC,L_CON,IBLOCK)
              IF (CCOEF.NE.0.0D0) THEN
                IBASCOUNT = IBASCOUNT + 1
              ENDIF
            ENDDO
 
            WRITE(LUTMP2,*) BLNK(1+J),' ',IBASCOUNT

            !.... increment napg (of cap module)
           ! napg = napg + IBASCOUNT

            IF (IPRINT.GE.3) THEN
             WRITE(LUPRI,*) BLNK(1+J),' ',IBASCOUNT
            ENDIF

            IBASCOUNT = 0

            DO L_UNC =1 , NUC(IBLOCK)
              CCOEF = CPRIMU(L_UNC,L_CON,IBLOCK)
              IF (CCOEF.NE.0.0D0) THEN
                IBASCOUNT = IBASCOUNT + 1
                write(LUTMP2,*) IBASCOUNT,
     &          PRIEXP(IEXP+INUC+L_UNC),CCOEF
                IF (IPRINT.GE.2) THEN
                  write(LUPRI,*) IBASCOUNT,
     &            PRIEXP(IEXP+INUC+L_UNC),CCOEF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          INUC = INUC + NUC(IBLOCK)
          IBLOCK = IBLOCK + 1
       ENDDO
      ENDDO

      CALL QEXIT('CAPWRIBAS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wribas2 */
      SUBROUTINE CAPWRIBAS2(IQM,NBLCK,JCO,NUC,NRC,SEG,IPCON,
     &     KATOM,KANG,KBLOCK,KPRIM,CPRIMU,
     &     LUTMP3,IPOSEXP,IPRINT)
C*************************************************************************
C
C    Fill the file LU=LUTMP3 with DIRAC atomic basis functions order
C
C    Called from CAPOUT two-times: for Large and for Small component 
C    basis functions.
C
C    Written by Miro Ilias (mainly by cat&paste), Tel Aviv, January 2007
C
C*************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "primit.h"

#include "ccom.h"
#include "mxcent.h"
#include "nuclei.h"

      LOGICAL SEG
      CHARACTER*1 BLNK(5)
      CHARACTER CHRSEG*10,DEG*2

      DIMENSION IQM(KATOM),NBLCK(KATOM),
CMI  &      JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &      JCO(KANG,KATOM),NUC(*),NRC(*),
     &      CPRIMU(KPRIM,KPRIM,KBLOCK),SEG(KBLOCK)
      DIMENSION IPCON(KPRIM),IUNNO(MXCENT,2)

      INTEGER FINDIUNNO

!radovan: initialize iunno
!         valgrind issues use of uninitialized values
!         in findiunno
      iunno = 0

C ... make unique numbers on centers
      ICENT   = 0
      IX = 0
      DO I = 1,NONTYP
         DO  N = 1,NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            DO LL = 1, NDEG
             IX = IX + 1
             IF (IX.LE.MXCENT) THEN
                IUNNO(IX,1) = ICENT
                IUNNO(IX,2) = LL
             ELSE
                CALL QUIT(
     &    'CAPWRIBAS2: Increase size of allocation of IUNNO(:) MXCENT!')
             ENDIF
              IF (IPRINT.GE.3) THEN
              write(lupri,*) ' ICENT=',IUNNO(IX,1),' LL=',IUNNO(IX,2),
     &        ' indx unique=',IX
              ENDIF
            ENDDO
          ENDDO
      ENDDO
      
C.   ... whether large or small component set
      IF(IPOSEXP.EQ.0) THEN
        IPRIMD  = 0
        IORBD   = 0
        IPRIM   = 0
      ELSE
        IPRIMD  = NPLRG
        IORBD   = NLARGE
        IPRIM   = NPLSH
      ENDIF

      ICENT   = 0
      IBLOCK  = 0
      DO 500 I = 1,NONTYP
         DO 510 N = 1,NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = IBLOCK
            DO 600 J = 1,IQM(I)
               NCOMP = KHK(J)
               DO 610 K = 1,JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IPSTRT = IPRIM + 1
                  IPRIM  = IPRIM + NNUC
                  DO 700 L = 1,NNRC
                     JPRIM = 0
                     JPRIMD = IPRIMD
                     DO 750 M = IPSTRT,IPRIM
                        JPRIMD = JPRIMD + 1
                        IF(PRICCF(M,L).NE.0.00D0) THEN
                          JPRIM = JPRIM + 1
                          IPCON(JPRIM) = JPRIMD
                        ENDIF
  750                CONTINUE
                     IOFF = 0
                     DEG = '  '
C                    ... run over degeneracy of nuclei...
                     DO 810 LL = 1,NDEG

                       IF(NDEG.GT.1) WRITE(DEG,'(A1,I1)') '#',LL
C                      IF(NDEG.GT.1) WRITE(LUTMP3,'(A1,I1)') '#',LL

                       ITYP = NHKOFF(J)
                       DO 800 ICOMP = 1, NCOMP
                         ITYP = ITYP + 1
                         IORBD = IORBD + 1

                         IUNNO2 = FINDIUNNO(ICENT,LL,IUNNO,IPRINT)

                         IF (IPRINT .GT. 1) THEN
                            WRITE(LUPRI,1040)
     &                      IORBD,NAMN(ICENT),DEG,GTOTYP(ITYP),
     &                      ((IPCON(M)+IOFF),M=1,JPRIM)
                         END IF

C                           WRITE(LUTMP3,'(I4,2X,A4,I2,I2,I2,2X,A4)')
C    &                   IORBD,NAMN(ICENT),ICENT,LL,IUNNO2,GTOTYP(ITYP)

                            WRITE(LUTMP3,'(I4,2X,I3,2X,A4)')
     &                   IORBD,IUNNO2,GTOTYP(ITYP)

                        IOFF = IOFF + NNUC

  800                  CONTINUE
  810                CONTINUE
  700             CONTINUE
                  IPRIMD = IPRIMD + NNUC*NDEG*NCOMP
  610          CONTINUE
  600       CONTINUE
  510    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
  500 CONTINUE
C     IF (IPREAD .GT. 0) THEN
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,'(/)')
      END IF
 1040 FORMAT(I4,2X,A4,A2,2X,A4,24I4)
 1050 FORMAT(I4,2X,A4,4X,A4,12I6)
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck FINDIUNNO */
      INTEGER FUNCTION FINDIUNNO(ICENT,LL,IUNNO,IPRINT)
C***********************************************************************
C
C Find out unique index of the nuclei.
C 
C Used byt CAPWRIBAS2.
C
C Written by Miro Ilias, Tel Aviv, February 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      DIMENSION IUNNO(MXCENT,2)

      INNO = 0
      DO I = 1 , MXCENT
        if (iunno(i, 1) == icent) then
           if (iunno(i, 2) == ll) then
              inno = i
           end if
        end if
      ENDDO

      IF (INNO.EQ.0) THEN
       call QUIT('error in INTEGER FUNCTION FINDIUNNO !')
      ENDIF
  
      FINDIUNNO = INNO

      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck capsymout */
      SUBROUTINE CAPSYMOUT
!*****************************************************************************
! 
!    Prepare identical CAPDIRSYM1.inp,CAPDIRSYM2.inp files for subsequent 
!    symmetry adaptation of CAP integrals obtained from Sajeev's code
!
!    Called from CAPOUT.
!
!    Copy&Paste mainly by Miro Ilias, Tel Aviv, 2007
!
!*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pgroup.h"
C#include "pincom.h"
#include "cbirea.h"
#include "aosotr.h"
#include "chrsgn.h"
#include "dummy.h"
      CHARACTER COMP(2)*1
C
#include "ibtfun.h"

      LUTMP1=-1
      CALL GPOPEN(LUTMP1,
     &'CAPDIRSYM1.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)

      LUTMP2=-1
      CALL GPOPEN(LUTMP2,
     &'CAPDIRSYM2.inp',' ',' ','FORMATTED',IDUMMY,.FALSE.)

      WRITE(LUTMP1,*) 'DIRAC symmetry adapted atomic basis functions'
      WRITE(LUTMP2,*) 'DIRAC symmetry adapted atomic basis functions'

C*****************************************************************************
C       S Y M M E T R Y      O R B I T A L S
C*****************************************************************************
      IF (IPREAD.GE.2) THEN
        CALL HEADER('CAPSYMOUT: Symmetry Orbitals',1)
        WRITE (LUPRI,'(A,6X,8I5)')
     &  '  Number of orbitals in each symmetry: ',(NAOS(I),I=1,MAXREP+1)
      ENDIF
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
        IF (IPREAD.GE.2) THEN
          WRITE (LUPRI,'(A,8I5)')
     &   '  Number of large orbitals in each symmetry: ',
     &    (NCOS(I,1),I=1,MAXREP+1)
          WRITE (LUPRI,'(A,8I5)')
     &   '  Number of small orbitals in each symmetry: ',
     &    (NCOS(I,2),I=1,MAXREP+1)
        ENDIF
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
      I=0
      DO 10 LA=1,MAXREP+1
        NBI=NAOS(LA)
        IF (NBI.GT.0) THEN
          IF (IPREAD.GE.2) THEN
           WRITE (LUPRI,'(//A,2X,A3,A1,I2,A1/)')
     &          '  Symmetry',REP(LA-1),'(',LA,')'
          ENDIF
          DO 20 L=1,NBI
            I=I+1
            ICENT = IPCEN(I)
            J=NUCDEG(ICENT)

            IF (IPREAD.GE.2) THEN
              WRITE (LUPRI,1060) I,NAMN(ICENT),COMP(ICLASS(I)),
     +          GTOTYP(IPTYP(I)), ITRAN(I,1),
     *         (CHRSGN(NINT(CTRAN(I,K))),ITRAN(I,K),K=2,J)

              WRITE (LUPRI,*) I,'  >>>>  ',J
              WRITE (LUPRI,*) 
     +         '           ',ITRAN(I,1),' ',
C    *         ( CHRSGN(NINT(CTRAN(I,K))),' ',
C    &            CTRAN(I,K),' ',NINT(CTRAN(I,K)),' ',
     &         ( (NINT(CTRAN(I,K))*ITRAN(I,K)), K=2, J )
            ENDIF

C                  ... write to the file
C      ==================================================
            WRITE (LUTMP1,*) I,'       ',J
            WRITE (LUTMP2,*) I,'       ',J

            WRITE (LUTMP1,*)
     +        '               ',ITRAN(I,1), ' ',
     *    ( (NINT(CTRAN(I,K))*ITRAN(I,K)), K = 2, J )

            WRITE (LUTMP2,*)
     +        '               ',ITRAN(I,1), ' ',
     *    ( (NINT(CTRAN(I,K))*ITRAN(I,K)), K = 2, J )

   20     CONTINUE
        ELSE
          IF (IPREAD.GE.2) THEN
            WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1)')
     *      'No orbitals in symmetry',REP(LA-1),'(',LA,')'
          ENDIF
        END IF
  10  CONTINUE
C     IF (IPREAD .GT. 10) THEN
C        WRITE(LUPRI,'(10X,A,/)') 'Symmetry pointer indices'
C        WRITE(LUPRI,'(4X,8I4)')
C    &      ((IPTSYM(I,J),J = 0,MAXREP),I=1,NORBS)
C        IF (NFLOAT .GT. 0) WRITE (LUPRI,'(/10X,A,I4)')
C    *       'Number of floating orbitals :',NFLOAT
C     END IF

C     WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
C    &        ' Symmetries of electric field:',
C    &        (REP(ISYMAX(I,1)),'(',ISYMAX(I,1)+1,')',I=1,3)
C     WRITE (LUPRI,'(/1X,A,3(2X,A3,A1,I1,A1))')
C    &        ' Symmetries of magnetic field:',
C    &        (REP(ISYMAX(I,2)),'(',ISYMAX(I,2)+1,')',I=1,3)

 1060 FORMAT(I5,3X,A6,1X,A1,3X,A4,5X,I5,7(2X,A,1X,I5))

      CALL GPCLOSE(LUTMP1,'KEEP')
      CALL GPCLOSE(LUTMP2,'KEEP')

      IF (IPREAD.GE.0) THEN
       WRITE(LUPRI,'(2X,A)') 
     & 'CAPSYMOUT: Files "CAPDIRSYM1.inp" and "CAPDIRSYM2.inp" for'//
     & ' CAP proper basis function order ready.'
      ENDIF

      RETURN
      END
#else
! miro: fix for OpenMPI ifort (mpif90); empty *.F file not accepted 
      SUBROUTINE DUMMY_ROUTINE
      END
#endif
