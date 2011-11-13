C 
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
C     Generates boxes for the local density-fitting scheme. The local density- 
C     fitting scheme fits the density of a box and a surrounding buffer.
C     The fitted density of the box is used and the fitted density in the buffer
C     removed.
      SUBROUTINE setfck3box(KNUMBOX,KCBOX,LFREE,KFREE,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION WORK(*)
C
      CALL QENTER('SETBOX')
C
      KNUMBOX = KFREE
      KCBOX   = KNUMBOX + NIBOX * MXCENT
      KLAST   = KCBOX + 9 * MXCENT
      LWRK    = KLAST - KFREE
      IF (LWRK .GT. LWORK) 
     &  CALL STOPIT('LIN_DF','SETBOX',LWRK,LFREE)
      LFREE   = LFREE - LWRK
      KFREE   = KLAST
C
      CALL SETGLOBALBOX(WORK(KNUMBOX),WORK(KCBOX))
      IBOX = 1
      NBOX = 1

      IF (IPRINT .GT. 10) THEN
        CALL HEADER('SETFCK3BOX',-1)
        WRITE(LUPRI,*) 'Global box given by (before any splitting):'
        WRITE(LUPRI,*) '********* GLOBAL BOX *********'
        CALL PRINTBOX(1,WORK(KNUMBOX),WORK(KCBOX),LUPRI)
      ENDIF

C     Partitioning into boxes (with #orb in box+buff < max-value)
 10   IF (IBOX.LE.NBOX) THEN
        NORBBOXBUF = NUMORBBOXBUF(IBOX,WORK(KNUMBOX))
        IF (NORBBOXBUF .LE. MAXORBBOX) THEN
          IBOX = IBOX + 1
        ELSE
          CALL SPLITBOX(IBOX,NBOX,WORK(KNUMBOX),WORK(KCBOX),IPRINT)
        ENDIF
        GOTO 10
      ENDIF
      NBOXES = NBOX

      IF (IPRINT .GT. 10) THEN
        WRITE(LUPRI,*) ' '
        WRITE(LUPRI,*) ' '
        WRITE(LUPRI,*) 'Number of boxes                      ', NBOX
        WRITE(LUPRI,*) 'Buffer range                         ', RBUFFER
        WRITE(LUPRI,*) 'Max number of orbitals in box+buffer ', 
     &                 MAXORBBOX
        IF (IPRINT .GT. 25) THEN
          WRITE(LUPRI,*) ' '
          WRITE(LUPRI,*) ' '
          DO I=1,NBOX
            WRITE(LUPRI,*) '********* BOX number ', I, ' *********'
            CALL PRINTBOX(I,WORK(KNUMBOX),WORK(KCBOX),LUPRI)
          ENDDO
        ENDIF
      ENDIF
      CALL QEXIT('SETBOX')
      RETURN
      END

      SUBROUTINE SET_BOX_DEFAULTS()
#include "implicit.h"
#include "fck3box.h"
      MAXORBBOX = 10000
      MAXSPLIT  = 100
      RBUFFER   = 15.0d0
      RETURN
      END
  
      SUBROUTINE GALPHA_FROM_3CENT(GALPHA,NA,DMAT,NB,THREE_CENT,
     &                              INDBOX,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
      DIMENSION GALPHA(NA),DMAT(NB,NB),THREE_CENT(NB,NB,NA),
     &          INDBOX(NAOBTCH), IAB(NB), IAUX(NA)
      LOGICAL IAB, IAUX, FIRSTAUX

      CALL DZERO(GALPHA,NA)

      DO I=1,NA
        IAUX(I) = .FALSE.
      ENDDO
      DO I=1,NB
        IAB(I)  = .FALSE.
      ENDDO

      FIRSTAUX = .TRUE.
      DO I=1,NAOBTCH
       ICLASS = KCLSBTCH(I)
       IF (FIRSTAUX .AND. (ICLASS .EQ. 3)) THEN
         IOFF = KAOBTCH(1,I)
         FIRSTAUX = .FALSE.
       ENDIF
       IF (INDBOX(I) .NE. 0) THEN
        NHK  = NHKTBTCH(I)
        KHK  = NHK * 2 - 1
        NCTF = NCTFBTCH(I)
        IF (ICLASS .EQ. 1) THEN
          IBASE = KAOBTCH(1,I)
          DO J=1,KHK*NCTF
            IAB(IBASE) = .TRUE.
            IBASE = IBASE + 1
          ENDDO
        ELSE
          IBASE = KAOBTCH(1,I) - IOFF
          DO J=1,KHK*NCTF
            IBASE = IBASE + 1
            IAUX(IBASE) = .TRUE.
          ENDDO
        ENDIF
       ENDIF
      ENDDO

      DO I=1,NA
       IF (IAUX(I)) THEN
        DO J=1,NB
         IF (IAB(J)) THEN
          DO K=1,NB
           IF (IAB(K)) THEN
            GALPHA(I) = GALPHA(I) + THREE_CENT(J,K,I) * DMAT(J,K)
           ENDIF
          ENDDO
         ENDIF
        ENDDO
       ENDIF
      ENDDO
        

      RETURN
      END

C     Routine for adding a projection of a vector from a global to a 
C     local region or vice versa. Project =.TRUE. means that a local 
C     vector should be projected to a global vector via ind_ab. Example let 
C       A=(1,1,1,1,1), NA=5
C       B = (3,5), NB=2 
C       IND_AB=(2,4) and project=.FALSE.
C     The result then becomes
C       A=(1,4,1,6,1)
C
      SUBROUTINE ADD_PROJECT_VEC(A,NA,B,NB,IND_AB,PROJECT,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL PROJECT
      DIMENSION A(NA),B(NB),IND_AB(NB)
      IF (PROJECT) THEN
        DO I=1,NB
          B(I) = B(I) + A(IND_AB(I))
        ENDDO
      ELSE
        DO I=1,NB
          A(IND_AB(I)) = A(IND_AB(I)) + B(I)
        ENDDO
      ENDIF
      IF (IPRINT .GT. 35) THEN
        CALL HEADER('ADD_PROJECT_VEC',-1)
        IF (PROJECT) THEN
          WRITE(LUPRI,*) 'From global vector A to local vector B'
        ELSE
          WRITE(LUPRI,*) 'From local vector B to global vector A'
        ENDIF
        WRITE(LUPRI,*) 'Via the index:'
        DO I=1,NB
          WRITE(LUPRI,*) I, IND_AB(I)
        ENDDO
        WRITE(LUPRI,*) 'Global vector A:'
        CALL OUTPUT(A,1,1,1,NA,1,NA,1,LUPRI)
        WRITE(LUPRI,*) 'Local vector B:'
        CALL OUTPUT(B,1,1,1,NB,1,NB,1,LUPRI)
      ENDIF
      RETURN
      END


C     Routine for projection a vector from a global to a local region
C     or vice versa. Project =.TRUE. means that a local vector should 
C     be projected to a global vector via ind_ab. Example let 
C       A=(1,1,1,1,1), NA=5
C       B = (3,5), NB=2 
C       IND_AB=(2,4) and project=.FALSE.
C     The result then becomes
C       A=(1,3,1,5,1)
C
      SUBROUTINE PROJECT_VEC(A,NA,B,NB,IND_AB,NDIM,PROJECT,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL PROJECT
      DIMENSION A(NA,NDIM),B(NB,NDIM),IND_AB(NB)

      DO IDIM=1,NDIM
        IF (PROJECT) THEN
          DO I=1,NB
            B(I,IDIM) = A(IND_AB(I),IDIM)
          ENDDO
        ELSE
          DO I=1,NB
            A(IND_AB(I),IDIM) = B(I,IDIM)
          ENDDO
        ENDIF
      ENDDO
      IF (IPRINT .GT. 35) THEN
        CALL HEADER('PROJECT_VEC',-1)
        IF (PROJECT) THEN
          WRITE(LUPRI,*) 'From global vector A to local vector B'
        ELSE
          WRITE(LUPRI,*) 'From local vector B to global vector A'
        ENDIF
        WRITE(LUPRI,*) 'Via the index:'
        DO I=1,NB
          WRITE(LUPRI,*) I, IND_AB(I)
        ENDDO
        DO IDIM=1,NDIM
        WRITE(LUPRI,*) 'Global vector A, IDIM = ', IDIM, ' :'
        CALL OUTPUT(A(1,IDIM),1,1,1,NA,1,NA,1,LUPRI)
        WRITE(LUPRI,*) 'Local vector B, IDIM = ', IDIM, ' :'
        CALL OUTPUT(B(1,IDIM),1,1,1,NB,1,NB,1,LUPRI)
        ENDDO
      ENDIF
      RETURN
      END


      SUBROUTINE PROJECT_MAT(A,NA,B,NB,IND_AB,PROJECT,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL PROJECT
      DIMENSION A(NA,NA),B(NB,NB),IND_AB(NB)
      IF (PROJECT) THEN
        DO I=1,NB
          DO J=1,NB
            B(I,J) = A(IND_AB(I),IND_AB(J))
          ENDDO
        ENDDO
      ELSE
        DO I=1,NB
          DO J=1,NB
            A(IND_AB(I),IND_AB(J)) = B(I,J)
          ENDDO
        ENDDO
      ENDIF
      IF (IPRINT .GT. 35) THEN
        CALL HEADER('PROJECT_MAT',-1)
        IF (PROJECT) THEN
          WRITE(LUPRI,*) 'From global matrix A to local matrix B'
        ELSE
          WRITE(LUPRI,*) 'From local matrix B to global matrix A'
        ENDIF
        WRITE(LUPRI,*) 'Via the index:'
        DO I=1,NB
          WRITE(LUPRI,*) I, IND_AB(I)
        ENDDO
        WRITE(LUPRI,*) 'Global matrix A:'
        CALL OUTPUT(A,1,NA,1,NA,NA,NA,1,LUPRI)
        WRITE(LUPRI,*) 'Local matrix B:'
        CALL OUTPUT(B,1,NB,1,NB,NB,NB,1,LUPRI)
      ENDIF
      RETURN
      END

C     Subroutine index_loc_glob creates an index that makes a projection 
C     from the global to the local auxilliary orbitals. I.e. so that the
C     density-fitting coeficients obtained by fitting the density in the 
C     box will be put in the correct places in the global density-fitting 
C     coefficients
C     INPUT:  NORB
C             IND - see find_actoa_boxbuf
C     OUTPUT: INDEX_LG - If the local orbital index is I, the global
C                        orbital index is INDEX_LG(I)
      SUBROUTINE INDEX_LOC_GLOB(INDEX_LG,NORB,IND,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
#include "priunit.h"
      DIMENSION INDEX_LG(NORB),IND(NAOBTCH)
      LOGICAL FIRSTAUX
      FIRSTAUX=.TRUE.
      DO IBATCH=1,NAOBTCH
        ICLASS = KCLSBTCH(IBATCH)
        IF (ICLASS .EQ. 3) THEN
          NINDEX = IND(IBATCH)
          IF (FIRSTAUX) THEN
            IOFF = KAOBTCH(1,IBATCH) - 1
            FIRSTAUX = .FALSE.
          ENDIF
          IF (NINDEX .NE. 0) THEN
            NHK    = NHKTBTCH(IBATCH)
            KHK    = NHK * 2 -1
            NCTF   = NCTFBTCH(IBATCH)
            IAUXORB = KAOBTCH(1,IBATCH) - IOFF
            NTYP   = NTYPBTCH(IBATCH)
            IF (NTYP.EQ.1) THEN
              DO I=NINDEX,NINDEX + KHK*NCTF - 1
                INDEX_LG(I) = IAUXORB
                IAUXORB = IAUXORB + 1
              ENDDO
            ELSE
              DO J=1,NTYP
                KHK = J*2-1
                DO I=NINDEX,NINDEX + KHK*NCTF - 1
                  INDEX_LG(I) = IAUXORB
                  IAUXORB = IAUXORB + 1
                ENDDO
                NINDEX = NINDEX + KHK*NCTF
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C     Printouts 
      IF (IPRINT.GT.25) THEN
        WRITE(LUPRI,*) 'Index for projection from global to local ',
     &                 '(or opposite) density-fitting:'
        DO I=1,NORB
          WRITE(LUPRI,*) I, INDEX_LG(I)
        ENDDO
      ENDIF
      RETURN
      END


      SUBROUTINE INDEX_GLOBAL_LOC(I_GL,NGLOBAL,I_LG,NLOCAL)
#include "implicit.h"
      DIMENSION I_GL(NGLOBAL), I_LG(NLOCAL)
      CALL IZERO(I_GL,NGLOBAL)
      DO I=1,NLOCAL
        I_GL(I_LG(I)) = I
      ENDDO
      RETURN
      END

C     Function number_of_boxes returns the number of boxes
      FUNCTION NUMBER_OF_BOXES()
#include "implicit.h"
#include "fck3box.h"
      NUMBER_OF_BOXES = NBOXES
      RETURN
      END



C     Subroutine pick_act_ind activates all batches wid ind .ne. 0, and the 
C     zero-exponent AO-batch.
      SUBROUTINE PICK_ACT_IND(IND,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION IND(NAOBTCH)
      DO I=1,NAOBTCH
        ICLASS = KCLSBTCH(I)
        ACTVAOBTCH(I) = (.NOT. (IND(I) .EQ. 0)) .OR. (ICLASS .EQ. 0 )
      ENDDO
      RETURN
      END



C     Subroutine find_actao_boxbuf finds the active batches in the box and the buffer
C     INPUT:   NUMBOX, CORDBOX, IBOX
C     OUTPUT:  INDBOX - for each auxilliary AO-batch it gives the starting 
C                       orbital-index for an AO in the box or buffer region around
C                       the box.
C                       A regular AO in the box or buffer gets a value -1
C                       All auxilliary AO's in the box get indexing from 1 to NBOX
C                       All auxilliary AO's in the buffer get indexing from NBOX+1 to NBOXBUF
C                       All AO's not in the box gets a value of 0
      SUBROUTINE FIND_ACTAO_BOXBUF(INDBOX,IBOX,NUMBOX,CORDBOX,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
#include "priunit.h"
      DIMENSION INDBOX(NAOBTCH), NUMBOX(NIBOX,*),
     &          CORDBOX(3,3,*),CORDBOXBUF(3,3)
      LOGICAL BOX, BOXBUF, AUX, ISINBOX

      SIGNC = 1.d0
      MINUS = -1.d0
      DO J=1,2
        SIGNC = SIGNC*MINUS
        DO I=1,3
          CORDBOXBUF(I,J)=CORDBOX(I,J,IBOX) + SIGNC*RBUFFER
        ENDDO
      ENDDO
      IORBBOX = 1
      IORBBUF = NUMBOX(2,IBOX) + 1
      DO IBATCH=1, NAOBTCH
        LCLASS = KCLSBTCH(IBATCH)
        AUX = LCLASS .EQ. 3
        IF ((LCLASS .EQ. 1) .OR. AUX) THEN
          NCENT  = NCNTBTCH(IBATCH)
          BOXBUF = ISINBOX(CORD(1,NCENT), CORDBOXBUF)
          IF (BOXBUF) THEN
            BOX = ISINBOX(CORD(1,NCENT), CORDBOX(1,1,IBOX))
            NHK    = NHKTBTCH(IBATCH)
            KHK    = NHK * 2 -1
            NTYP   = NTYPBTCH(IBATCH)
            IF (NTYP .NE. 1) KHK = NTYP*NTYP
            NCTF   = NCTFBTCH(IBATCH)
            IF (BOX) THEN
              IF (AUX) THEN
                INDBOX(IBATCH) = IORBBOX
                IORBBOX = IORBBOX + KHK * NCTF
              ELSE
                INDBOX(IBATCH) = -1
              ENDIF
            ELSE
              IF (AUX) THEN
                INDBOX(IBATCH) = IORBBUF
                IORBBUF = IORBBUF + KHK * NCTF
              ELSE
                INDBOX(IBATCH) = 0
              ENDIF
            ENDIF
          ELSE
            INDBOX(IBATCH) = 0
          ENDIF
        ENDIF
      ENDDO
C     Output
      IF (IPRINT .GT. 25) THEN
        CALL HEADER('Output from FIND_ACTAO_BOXBUF',-1)
        WRITE(LUPRI,*) 'Active batches in box number ', IBOX, ':'
        DO I=1,NAOBTCH - 1
          IF (.NOT. INDBOX(I) .EQ. 0) THEN
            WRITE(LUPRI,*) 'AO-batch ',I,' active, index = ',INDBOX(I)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END




      SUBROUTINE PRINTBOX(IBOX,NUMBOX,CORDBOX,LPRI)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION NUMBOX(NIBOX,*), CORDBOX(3,3,*)
      WRITE(LPRI,*) 'Number of orbitals in box and buffer', 
     &              NUMBOX(1,IBOX)
      WRITE(LPRI,*) 'Number of orbitals in box           ', 
     &              NUMBOX(2,IBOX)
      WRITE(LPRI,*) 'Number of orbitals in buffer        ', 
     &              NUMBOX(3,IBOX)
      WRITE(LPRI,*) 'Number of times box has been split  ',
     &              NUMBOX(4,IBOX)
      WRITE(LPRI,*) 'Coordinates corner 1 (x,y,z)        ',
     &              (CORDBOX(I,1,IBOX), I=1,3)
      WRITE(LPRI,*) 'Coordinates corner 2 (x,y,z)        ',
     &              (CORDBOX(I,2,IBOX), I=1,3)
      WRITE(LPRI,*) 'Dimensions of the box (x,y,z)       ',
     &              (CORDBOX(I,3,IBOX), I=1,3)
      RETURN
      END


      FUNCTION NUMORBBOX(IBOX,NUMBOX)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION NUMBOX(NIBOX,*)
      NUMORBBOX = NUMBOX(2,IBOX)
      RETURN
      END

      FUNCTION NUMORBBOXBUF(IBOX,NUMBOX)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION NUMBOX(NIBOX,*)
      NUMORBBOXBUF = NUMBOX(1,IBOX)
      RETURN
      END


      SUBROUTINE SPLITBOX(IBOX,NBOX,NUMBOX,CORDBOX,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
      PARAMETER (D5 = 0.5d0)
      DIMENSION NUMBOX(NIBOX,*), CORDBOX(3,3,*)
      LOGICAL   HALF
      NBOX   = NBOX + 1
     
      HALF = .TRUE.
      IF (HALF) THEN 
        NSPLIT = NUMBOX(4,IBOX) + 1
        IF (NSPLIT .GT. MAXSPLIT) THEN
          WRITE(LUPRI,*) 'NSPLIT = ', NSPLIT,
     &      ' greater that MAXSPLIT = ', MAXSPLIT
          CALL QUIT('NSPLIT .GT. MAXSPLIT')
        ENDIF
C       Pick up largest extent
        RMAX = 0.0d0
        DO I=1,3
          RCUR = CORDBOX(I,3,IBOX)
          IF (RCUR.GT.RMAX) THEN
            RMAX = RCUR
            IMAX = I
          ENDIF
        ENDDO
        CALL COPY_BOX(IBOX,NBOX,NUMBOX,CORDBOX)
        CORDBOX(IMAX,3,IBOX) = CORDBOX(IMAX,3,IBOX) * D5
        CORDBOX(IMAX,2,IBOX) = CORDBOX(IMAX,1,IBOX) + 
     &                           CORDBOX(IMAX,3,IBOX)
        CORDBOX(IMAX,3,NBOX) = CORDBOX(IMAX,3,IBOX)
        CORDBOX(IMAX,1,NBOX) = CORDBOX(IMAX,2,IBOX)

        CALL COUNT_ORBBOXBUF(INBOX,INBUFFER,IBOX,CORDBOX,IPRINT)
        NUMBOX(1,IBOX) = INBOX + INBUFFER
        NUMBOX(2,IBOX) = INBOX
        NUMBOX(3,IBOX) = INBUFFER
        NUMBOX(4,IBOX) = NSPLIT
        IF (INBOX .EQ. 0) THEN
          CALL COPY_BOX(NBOX,IBOX,NUMBOX,CORDBOX)
          NBOX = NBOX - 1
          CALL COUNT_ORBBOXBUF(INBOX,INBUFFER,IBOX,CORDBOX,IPRINT)
          IF (INBOX .EQ. 0) CALL QUIT('Error in Splitbox')
          NUMBOX(1,IBOX) = INBOX + INBUFFER
          NUMBOX(2,IBOX) = INBOX
          NUMBOX(3,IBOX) = INBUFFER
          NUMBOX(4,IBOX) = NSPLIT
        ELSE
          CALL COUNT_ORBBOXBUF(INBOX,INBUFFER,NBOX,CORDBOX,IPRINT)
          NUMBOX(1,NBOX) = INBOX + INBUFFER
          NUMBOX(2,NBOX) = INBOX
          NUMBOX(3,NBOX) = INBUFFER
          NUMBOX(4,NBOX) = NSPLIT
          IF (INBOX.EQ.0) NBOX = NBOX - 1
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE COPY_BOX(IBOX,IBOX_COPY,NUMBOX,CORDBOX)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
      DIMENSION NUMBOX(NIBOX,*), CORDBOX(3,3,*)
      DO I=1,NIBOX
        NUMBOX(I,IBOX_COPY) = NUMBOX(I,IBOX)
      ENDDO
      DO I=1,3
        DO J=1,3
          CORDBOX(I,J,IBOX_COPY) = CORDBOX(I,J,IBOX)
        ENDDO
      ENDDO
      RETURN
      
      END

      SUBROUTINE COUNT_ORBBOXBUF(INBOX,INBUFFER,IBOX,CORDBOX,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
      DIMENSION CORDBOX(3,3,*), CORDBOXBUF(3,2)
      LOGICAL BOXBUF, BOX, ISINBOX

      SIGNC = 1.d0
      MINUS = -1.d0
      DO J=1,2
        SIGNC = SIGNC*MINUS
        DO I=1,3
          CORDBOXBUF(I,J)=CORDBOX(I,J,IBOX) + SIGNC*RBUFFER
        ENDDO
      ENDDO
      INBOX = 0
      INBUFFER = 0
      DO IBATCH=1, NAOBTCH
        LCLASS = KCLSBTCH(IBATCH)
C       Only the auxilliary density-fitting basis is important for determining
C       the size of the boxes.
        IF (LCLASS .EQ. 3) THEN
          NCENT  = NCNTBTCH(IBATCH)
          BOXBUF = ISINBOX(CORD(1,NCENT), CORDBOXBUF)
          IF (BOXBUF) THEN
            BOX = ISINBOX(CORD(1,NCENT), CORDBOX(1,1,IBOX))
            NHK    = NHKTBTCH(IBATCH)
            KHK    = NHK * 2 -1
            NTYP   = NTYPBTCH(IBATCH)
            IF (NTYP .NE. 1) KHK = NTYP*NTYP
            NCTF   = NCTFBTCH(IBATCH)
            IF (BOX) THEN
              INBOX = INBOX + KHK * NCTF
            ELSE
              INBUFFER = INBUFFER + KHK * NCTF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C     Printouts
      IF (IPRINT .GT. 35) THEN
        CALL HEADER('COUNT_ORBBOXBUF',-1)
        WRITE(LUPRI,*) 'Corner 1 of buffer (x,y,z): ', 
     &                 (CORDBOXBUF(I,1),I=1,3)
        WRITE(LUPRI,*) 'Corner 2 of buffer (x,y,z): ', 
     &                 (CORDBOXBUF(I,2),I=1,3)
        WRITE(LUPRI,*) 'Corner 1 of box (x,y,z): ', 
     &                 (CORDBOX(I,1,IBOX),I=1,3)
        WRITE(LUPRI,*) 'Corner 2 of box (x,y,z): ', 
     &                 (CORDBOX(I,2,IBOX),I=1,3)
        WRITE(LUPRI,*) 'BOX NUMBER ', IBOX, ' INBOX=', INBOX, 
     &                 ' INBUFFER=', INBUFFER
      ENDIF
      RETURN
      END

      LOGICAL FUNCTION ISINBOX(CORD,BOX)
#include "implicit.h"
      DIMENSION CORD(3), BOX(3,2)
      ISINBOX = .TRUE.
      DO I=1,3
        IF ( (CORD(I) .GT. BOX(I,2)) .OR. (CORD(I) .LE. BOX(I,1)) ) 
     &    ISINBOX = .FALSE.
      ENDDO
      RETURN
      END

      SUBROUTINE SETGLOBALBOX(NUMBOX,CORDBOX)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
#include "fck3box.h"
#include "nuclei.h"
      PARAMETER (RANGEMAX = 1.d15)
C     CORDBOX: 1st index x,y,z, 2nd index min,max,range (max-min)
C     NUMBOX:  1  # of orbitals in box and buffer
C              2  # of orbitals in box
C              3  # of orbitals in buffer
C              4  # of times box has been split
      DIMENSION NUMBOX(NIBOX), CORDBOX(3,3)
      DO I=1,3
        CORDBOX(I,1) = RANGEMAX
        CORDBOX(I,2) = -RANGEMAX
      ENDDO
      NORBIT = 0
      DO IBATCH=1, NAOBTCH
        LCLASS = KCLSBTCH(IBATCH)
C       Only the auxilliary density-fitting basis is important for determining
C       the size of the boxes, and we use only the aux. basis here to be 
C       consequent.
        IF (LCLASS .EQ. 3) THEN
          NCENT  = NCNTBTCH(IBATCH)
          DO I=1,3
            CORDBOX(I,1) = MIN(CORDBOX(I,1),CORD(I,NCENT))
            CORDBOX(I,2) = MAX(CORDBOX(I,2),CORD(I,NCENT))
          ENDDO
          NHK    = NHKTBTCH(IBATCH)
          KHK    = NHK * 2 -1
          NTYP   = NTYPBTCH(IBATCH)
          IF (NTYP .NE. 1) KHK = NTYP*NTYP
          NCTF   = NCTFBTCH(IBATCH)
          NORBIT = NORBIT + KHK * NCTF
        ENDIF
      ENDDO
      DO I=1,3
        CORDBOX(I,1) = CORDBOX(I,1) - THRBOX
        CORDBOX(I,2) = CORDBOX(I,2) + THRBOX
        CORDBOX(I,3) = CORDBOX(I,2) - CORDBOX(I,1)
      ENDDO
      NUMBOX(1) = NORBIT
      NUMBOX(2) = NORBIT
      NUMBOX(3) = 0
      NUMBOX(4) = 0
      RETURN
      END






