C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE RTRACTL(CMO,SZTRA,WORK,LWORK)
C***********************************************************************
C
C     Input : ITRLVL - level of integrals transformation (see below)
C
C     Output: 4-index transformed integrals on disk.
C
C      
C     ITRLVL = -1: Transform all integrals (pq|rs)
C                     (which is more than we need, but for a quick start
C                      this is good)
C     ITRLVL =  0: Transform all integrals (uv|xy) for CI & MC energy
C     ITRLVL =  1: (pv|xy) for CI & MC gradients
C     ITRLVL =  2: (pq|xy) and (px|qy) for Newton-Raphson and NEO CI & MC
C                  (not implemented yet!)
C     ITRLVL =  3: (pq|ry) 
C     ITRLVL =  4: (pq|rs)
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER SZTRA(2,4)*(*)
      DIMENSION CMO(*),WORK(*)
C
#include "memint.h"
C
C
      CALL QENTER('RTRACTL')
      CALL RTKTIME(.TRUE.,7)
C
      IF ( IPROPT .GT. 3 ) THEN
         WRITE(LUPRI,'(3A,I1,A)')
     &        ' (RTRACTL) Integrals in transformation = ',OPT_INTTYP,
     &        ' (INTFLG=',IOPT_INTFLG,')'
         WRITE(LUPRI,'(A,L1)')
     &        '           No positronic integrals       = ',OPT_NOPFQ
      END IF
C
C
C     Print orbital strings
C
      IF ( IPROPT .GE. 3 ) THEN
         DO I = 1, 4
            DO J = 1,NFSYM
               WRITE(LUPRI,'(A,I1,A,I1,2A)')
     &              '* Orbital string ',I,
     &              ' for corep ',J,' = ',SZTRA(J,I)
            END DO
         END DO
      END IF
C
C     Transfer information to moltra common blocks
C
      CALL SETUP_DCBTRA(SZTRA,IOPT_INTFLG,OPT_SKIPEP,0)
C
C     Call 4-index transformation driver
C
      CALL RTRACTL1(CMO,WORK(KFREE),LFREE)
C
C
      CALL RTKTIME(.FALSE.,7)
      CALL QEXIT('RTRACTL')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck iniszt */
      SUBROUTINE INISZT(ITRLVL,SZTRA,USENOROT,NOROT)
C***********************************************************************
C
C     Set orbital strings for 4-index transformation.
C
C     USENOROT controls whether to use NOROT or not.
C
C     Written by J. Thyssen - Wed Mar 14 15:03:26 "MET 2001
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbopt.h"
C
      CHARACTER*(*) SZTRA(2,4)
      LOGICAL USENOROT
      DIMENSION NOROT(*)
C
C     ITRLVL = 0: (uu|uu) transformation
C     ITRLVL = 1: (gu|uu) transformation
C     ITRLVL = 2: (gg|uu) + (gu|gu) transformation
C     ITRLVL = 3: (gg|gu) transformation
C     ITRLVL = 4: (gg|gg) transformation
C
      IF (ITRLVL .EQ. 0) THEN
         CALL INISZ1('AAAA',SZTRA,USENOROT,NOROT)
      ELSE IF (ITRLVL .EQ. 1) THEN
         IF (OPT_SKIPEP .OR. OPT_NOPFQ) THEN
            CALL INISZ1('AAEA',SZTRA,USENOROT,NOROT)
         ELSE
            CALL INISZ1('AAGA',SZTRA,USENOROT,NOROT)
         END IF
      ELSE IF (ITRLVL .EQ. 3) THEN
         IF (OPT_SKIPEP) THEN
            CALL INISZ1('EAEE',SZTRA,USENOROT,NOROT)
         ELSE IF (OPT_NOPFQ) THEN
            CALL INISZ1('EAGE',SZTRA,USENOROT,NOROT)
         ELSE
            CALL INISZ1('GAGG',SZTRA,USENOROT,NOROT)
         END IF
      ELSE
         IF (OPT_SKIPEP) THEN
            CALL INISZ1('EEEE',SZTRA,USENOROT,NOROT)
         ELSE IF (OPT_NOPFQ) THEN
            CALL INISZ1('EEGE',SZTRA,USENOROT,NOROT)
         ELSE
            CALL INISZ1('GGGG',SZTRA,USENOROT,NOROT)
         END IF
      END IF
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inisz1 */
      SUBROUTINE INISZ1(C,SZTRA,USENOROT,NOROT)
C***********************************************************************
C
C     Set orbital string for 4-index transformation.
C
C     C is control:
C
C     'I' for inactive index,
C     'A' for active index,
C     'V' for virtual index,
C     'E' for electronic index,
C     'P' for positronic index,
C     'G' for general index.
C
C     Written by J. Thyssen - Feb 19 2001
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbidx.h"
C
      CHARACTER C*4
      CHARACTER*(*) SZTRA(2,4)
      LOGICAL USENOROT
      DIMENSION NOROT(*)
      DIMENSION ILIST(MAXORB,2), NLIST(2), IDX(MAXORB)
C
C     *******************
C     *** Index 1 & 2 ***
C     *******************
C
C     (1) 'AA'
C     (2) 'GG'
C     (3) 'EE'
C
      DO I = 1, 4
         DO IFSYM = 1, NFSYM
            SZTRA(IFSYM,I) = ' '
         END DO
      END DO
      IOFF0 = 0
      IF (NFSYM .GT. 1) THEN
         IT = 0
         DO J = 1, NORB(1)
            IF (.NOT. USENOROT .OR. NOROT(J) .EQ. 0) IT = IT + 1
         END DO
         IOFFG = IT
         IT = 0
         DO J = 1, NESH(1)
            IF (.NOT. USENOROT .OR. NOROT(IDXE2G(J)) .EQ. 0) IT = IT + 1
         END DO
         IOFFE = IT
      END IF
      print *,'hj: IOFFG, IOFFE:',IOFFG,IOFFE
C
C
      IF (C(1:2) .EQ. 'AA') THEN
C        
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,1),NIDXT(1),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,2),NIDXT(2),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX2,NIDX2)
C
      ELSE IF (C(1:2) .EQ. 'EA') THEN
C
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,1),NIDXT(1),
     &        NESH,IOFF0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,2),NIDXT(2),
     &        NASH,IOFFE,USENOROT,NOROT,IIDX2,NIDX2)
C        
      ELSE IF (C(1:2) .EQ. 'GA') THEN
C
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,1),NIDXT(1),
     &        NORB,IOFF0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,2),NIDXT(2),
     &        NASH,IOFFG,USENOROT,NOROT,IIDX2,NIDX2)
C        
      ELSE IF (C(1:2) .EQ. 'EE') THEN
C        
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,1),NIDXT(1),
     &        NESH,INT0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,2),NIDXT(2),
     &        NESH,INT0,USENOROT,NOROT,IIDX2,NIDX2)
C        
      ELSE IF (C(1:2) .EQ. 'GE') THEN
C        
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,1),NIDXT(1),
     &        NORB,INT0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,2),NIDXT(2),
     &        NESH,INTG,USENOROT,NOROT,IIDX2,NIDX2)
C        
      ELSE IF (C(1:2) .EQ. 'GG') THEN
C        
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,1),NIDXT(1),
     &        NORB,INT0,USENOROT,NOROT,IIDX1,NIDX1)
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,2),NIDXT(2),
     &        NORB,INT0,USENOROT,NOROT,IIDX2,NIDX2)
      ELSE
         CALL QUIT('*** ERROR in INISZ1 *** ' //
     &        'C(3:4) not implemented ' // C(3:4) )
      END IF
C     
chj#ifdef UNDEF
      do j = 1, 2
         write(6,*) 'index ',j, ' t --> g pointer'
         do i = 1, nidxt(j)
            write(6,*) i,idxt2g(i,j)
         end do
      end do
chj#endif
C     
      CALL INISZ3(IDXT2G(1,1),NIDXT(1),ILIST,NLIST)
      CALL INISZ4(ILIST,NLIST,SZTRA(1,1))
C
chj#ifdef UNDEF
      do j = 1, nfsym
         write(6,*) 'ilist',(ilist(i,j),i=1,nlist(j))
         write(6,*) 'sztra = ',sztra(j,1)
      end do
chj#endif
C     
      CALL INISZ3(IDXT2G(1,2),NIDXT(2),ILIST,NLIST)
      CALL INISZ4(ILIST,NLIST,SZTRA(1,2))
C
chj#ifdef UNDEF
      do j = 1, nfsym
         write(6,*) 'ilist',(ilist(i,j),i=1,nlist(j))
         write(6,*) 'sztra = ',sztra(j,2)
      end do
chj#endif
C
C
C     *******************
C     *** Index 3 & 4 ***
C     *******************
C
C
C     (same as above except that all offsets are replaced with 0)
C
C
      IF (C(3:4) .EQ. 'AA') THEN
C        
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,3),NIDXT(3),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,4),NIDXT(4),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
C
      ELSE IF (C(3:4) .EQ. 'EA') THEN
C
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,3),NIDXT(3),
     &        NESH,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,4),NIDXT(4),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
C        
      ELSE IF (C(3:4) .EQ. 'GA') THEN
C
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,3),NIDXT(3),
     &        NORB,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(IDXU2G,.TRUE.,IDXT2G(1,4),NIDXT(4),
     &        NASH,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
C        
      ELSE IF (C(3:4) .EQ. 'EE') THEN
C        
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,3),NIDXT(3),
     &        NESH,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,4),NIDXT(4),
     &        NESH,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
C        
      ELSE IF (C(3:4) .EQ. 'GE') THEN
C        
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,3),NIDXT(3),
     &        NORB,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(IDXE2G,.TRUE.,IDXT2G(1,4),NIDXT(4),
     &        NESH,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
C        
      ELSE IF (C(3:4) .EQ. 'GG') THEN
C        
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,3),NIDXT(3),
     &        NORB,IOFF0,USENOROT,NOROT,IIDX3,NIDX3)
         CALL INISZ2(DUMMY,.FALSE.,IDXT2G(1,4),NIDXT(4),
     &        NORB,IOFF0,USENOROT,NOROT,IIDX4,NIDX4)
      ELSE
         CALL QUIT('*** ERROR in INISZ1 *** ' //
     &        'C(3:4) not implemented ' // C(3:4) )
      END IF
C     
#ifdef UNDEF
      do j = 3, 4
         write(6,*) 'index ',j, ' t --> g pointer'
         do i = 1, nidxt(j)
            write(6,*) i,idxt2g(i,j)
         end do
      end do
#endif
C     
      CALL INISZ3(IDXT2G(1,3),NIDXT(3),ILIST,NLIST)
      CALL INISZ4(ILIST,NLIST,SZTRA(1,3))
C
#ifdef UNDEF
      do j = 1, nfsym
         write(6,*) 'ilist',(ilist(i,j),i=1,nlist(j))
         write(6,*) 'sztra = ',sztra(j,3)
      end do
#endif
C     
      CALL INISZ3(IDXT2G(1,4),NIDXT(4),ILIST,NLIST)
      CALL INISZ4(ILIST,NLIST,SZTRA(1,4))
C
#ifdef UNDEF
      do j = 1, nfsym
         write(6,*) 'ilist',(ilist(i,j),i=1,nlist(j))
         write(6,*) 'sztra = ',sztra(j,4)
      end do
C     
      write(6,*) 'iidx',iidx1,iidx2,iidx3,iidx4
      write(6,*) 'nidx',nidx1,nidx2,nidx3,nidx4
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inisz2 */
      SUBROUTINE INISZ2(IDX,IDXFLG,IDXT2G,NIDXT,
     &     NORB,IOFF,USENOROT,NOROT,IIDX,NIDX)
C***********************************************************************
C
C     Set orbital string for 4-index transformation.
C
C     Written by J. Thyssen - Wed Mar 14 17:20:32 "MET 2001
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
C
      DIMENSION IDX(*)
      LOGICAL IDXFLG
      DIMENSION IDXT2G(*)
      DIMENSION NORB(*)
      LOGICAL USENOROT
      DIMENSION NOROT(*)
      DIMENSION IIDX(2), NIDX(2)
C
      IA = 0
      IT = 0
      DO IFSYM = 1, NFSYM
         IIDX(IFSYM) = IT
         DO I = 1, NORB(IFSYM)
            IA = IA + 1
            IF (IDXFLG) THEN
               IG = IDX(IA)
            ELSE
               IG = IA
            END IF
            IF (.NOT. USENOROT .OR. NOROT(IG) .EQ. 0) THEN
               IT = IT + 1
               IDXT2G(IT) = IG
            END IF
         END DO
         IF (NFSYM .GT. 1) THEN
            DO I = NORB(IFSYM) + 1, IOFF
               IT = IT + 1
               IDXT2G(IT) = -1
            END DO
         END IF
         NIDX(IFSYM) = IT - IIDX(IFSYM)
      END DO
      NIDXT = IT
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inisz3 */
      SUBROUTINE INISZ3(IDX,NIDX,ILIST,NLIST)
C***********************************************************************
C
C     Transform IDX to ``normalized list''
C
C     Written by J. Thyssen - Wed Mar 14 17:20:32 "MET 2001
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "dcbidx.h"
#include "dcborb.h"
C
      DIMENSION IDX(*)
      DIMENSION ILIST(MAXORB,2), NLIST(2)
C
      NLIST(1) = 0 
      NLIST(2) = 0 
      DO I = 1, NIDX
         ID = IDX(I)
         IF (ID .GT. 0) THEN
            ITYP = IOBTYP(ID)
            ISYM = IFSMO(ID)
            NLIST(ISYM) = NLIST(ISYM) + 1
            IF (ITYP .EQ. JTPOSI) THEN
              ILIST(NLIST(ISYM),ISYM) = ID - IORB(ISYM) - NPSH(ISYM) - 1
            ELSE
              ILIST(NLIST(ISYM),ISYM) = ID - IORB(ISYM) - NPSH(ISYM)
            END IF
         END IF
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck inisz4 */
      SUBROUTINE INISZ4(ILIST,NLIST,SZTRA)
C***********************************************************************
C
C     Transform ILIST to orbital string.
C
C     Written by J. Thyssen - Wed Mar 14 17:20:32 "MET 2001
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
C
#include "dgroup.h"
C
      DIMENSION ILIST(MAXORB,2), NLIST(2)
      CHARACTER*(*) SZTRA(2)
C
      LOGICAL FIRST
      CHARACTER TEGN*1
C
      N = 1
      M = 12
      DO 10 I = 1, NFSYM
         FIRST = .TRUE.
         IF (NLIST(I) .EQ. 0) GOTO 10
         J = 1
         INTST = 0
 20      CONTINUE
            IF (J .GT. NLIST(I)) GOTO 30
            IF (INTST .EQ. 0) THEN
               INTST = ILIST(J,I)
               J = J + 1
            ELSE
               IF (ILIST(J-1,I) + 1 .NE. ILIST(J,I)) THEN
C                 ... interval finished
                  IF (FIRST) THEN
C                    ... no comma
                     IF (ILIST(J-1,I) .EQ. INTST) THEN
                        WRITE(SZTRA(I)(N:M), '(I5,7X)')
     &                       INTST
                     ELSE
                        WRITE(SZTRA(I)(N:M), '(I5,A2,I5)')
     &                       INTST, '..', ILIST(J-1,I)
                     END IF
                     N = M + 1
                     M = M + 13
                     FIRST = .FALSE.
                  ELSE
C                    ... with comma
                     IF (ILIST(J-1,I) .EQ. INTST) THEN
                        WRITE(SZTRA(I)(N:M), '(A1,I5,7X)')
     &                       ',',INTST
                     ELSE
                        WRITE(SZTRA(I)(N:M), '(A1,I5,A2,I5)')
     &                       ',', INTST, '..', ILIST(J-1,I)
                     END IF
                     N = M + 1
                     M = M + 13
                  END IF
                  INTST = 0
               ELSE
                  J = J + 1
               END IF
            END IF
            GOTO 20
 30      CONTINUE
C
C        last interval
C
         IF (INTST .NE. 0) THEN
            IF (FIRST) THEN
C              ... no comma
               IF (ILIST(J-1,I) .EQ. INTST) THEN
                  WRITE(SZTRA(I)(N:M), '(I5,7X)')
     &                 INTST
               ELSE
                  WRITE(SZTRA(I)(N:M), '(I5,A2,I5)')
     &                 INTST, '..', ILIST(NLIST(I),I)
               END IF
            ELSE
C              ... with comma
               IF (ILIST(J-1,I) .EQ. INTST) THEN
                  WRITE(SZTRA(I)(N:M), '(A1,I5,7X)')
     &                 ',',INTST
               ELSE
                  WRITE(SZTRA(I)(N:M), '(A1,I5,A2,I5)')
     &                 ',', INTST, '..', ILIST(NLIST(I),I)
               END IF
            END IF
         END IF
C
C        Check for overflow of SZTRA
C
         IF ( M .GT. LEN(SZTRA(I)) ) THEN
            WRITE(LUPRI,*)
     &         'INISZ4: length of SZTRA too short for char string:',
     &         SZTRA(I)
            CALL QUIT(
     &         'INISZ4: length of SZTRA too short for char string')
         END IF
C
C        Remove blanks
C
         NB = 0
         DO 50 J = 1,M
            TEGN = SZTRA(I)(J:J)
            IF (TEGN.EQ.' ') GOTO 50
            NB = NB + 1
            SZTRA(I)(NB:NB) = TEGN
 50      CONTINUE
         DO J = NB+1,M
            SZTRA(I)(J:J) = ' '
         END DO
C
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setup_dcbtra */
      SUBROUTINE SETUP_DCBTRA(SZ,INTFLG,LNOPAIR,IPRINT)
C***********************************************************************
C
C     Input : SZ     - orbitals strings
C             IPRINT - print level
C             INTFLG - integral classes to include in transformation
C
C     Output: none
C
C     Sets the common blocks used in the 4-index transformation.
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( DM1 = -1.0D00 )
C
      CHARACTER SZ(2,4)*(*)
      LOGICAL LNOPAIR
C
#include "dcbtra.h"
C     
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            TRA4_INDSTR(I,J) = SZ(J,I)
         END DO
      END DO
      DO I = 1,2
         TRA_INDSTR(I) = ' '
         TRA_CORSTR(I) = ' '
         DO J = 1,2
            TRA2_INDSTR(I,J) = ' '
         END DO
      END DO
C
C     /DCBTRI/
C
      ITRA_INTFLG = INTFLG
      ITRA_INTFL2 = INTFLG
      ITRA_INTFL4 = INTFLG
c     IPRTRA = IPRINT - 2
C TODO debug : extra print because no integrals for Levy-Leblond /hjaaj aug 2002
      IPRTRA = 4
      ISTRAT = 4
C
C     /DCBTRL/
C
      TRA_ANTIS = .FALSE.
      TRA_MSOUT = .TRUE.
      NO2IND = .TRUE.
      NO4IND = .FALSE.
      NOPAIR = LNOPAIR
      PRPTRA = .FALSE.
      RCORBS = .FALSE.
      FKRMC  = .TRUE.
      NOMDCINT = .TRUE.
C
C     Don't scatter 4INDX.. to all nodes, we only need on the master.
C
      MDCSCAT = .FALSE.
C      
C     /DCBTRR/
C
      THROUT = 1.0D-16
      SCRTRA = DM1
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rtractl1 */
      SUBROUTINE RTRACTL1(CMO,WORK,LWORK)
C***********************************************************************
C
C     Driver routine for KRMCSCF 4-index transformation
C     Based on PAMTRA from trapam.F.
C
C     Written by J. Thyssen - Dec 11 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "mxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
      CHARACTER*24 DAYTID
      LOGICAL TRASAM
      DIMENSION CMO(*),WORK(*),
     &          NSTR(2,0:2,4),KVEC(2,4),KQ(2,4),KIBE(2,4),
     &          IR(2,2),NR(2),NQ(2,4),NQT(4),KE(2,4),
     &          NSPC(2,0:2),KQC(2),NSTRT(4)
C
      CALL QENTER('RTRACTL1')
#include "memint.h"
      KFRSAV = KFREE
      DUM = DUMMY
C
      CALL GETTIM(CPU1,WALL1)
      WRITE(LUPRI,'(/A)')
     &     ' (RTRACTL1) Starting 4-index integral transformation.'
      IF (IPRTRA .GE. 3) THEN
         CALL TRAHI(IPRTRA,4)
      END IF
C
C     Set up index arrays for integral transformation
C     -----------------------------------------------
C
      IF(NOPAIR) THEN
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
      DO I = 1,4
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET('INTE',KVEC(IFRP,I),NR(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = -1
C         ... -1 tells NUMLST to _not_ reorder orbital list
          CALL  NUMLST(TRA4_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 IFRP,NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('PAMTRA',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO      
C
C     Print section
C     =============
C
      IF ( IPRTRA .GE. 3 ) THEN
         WRITE(LUPRI,'(/A)')
     &        '(RTRACTL1) Orbital ranges for 4-index transformation:'
         DO I = 1, NFSYM
            CALL TRAPRI(4,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),
     &           WORK(KVEC(I,3)),WORK(KVEC(I,4)),NSTR)
         END DO
      END IF
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 4 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 4
         NQT(I) = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQT(I) = NQT(I) + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM) 
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
      ENDDO
C
C     Allocate the memory for the coefficients, check their relations
C
      DO I = 1, 4
C
         DO J = I, 1, -1
            K = 2
            IF (NFSYM.EQ.1) K=1
            IF (TRASAM(WORK(KVEC(1,I)),WORK(KVEC(K,I)),
     &                 WORK(KVEC(1,J)),WORK(KVEC(K,J)),
     &          NSTR(1,0,I),NSTR(1,0,J))) ISAME(I) = J
         ENDDO
C
         IF (ISAME(I).NE.I) THEN
            KQ(1,I) = KQ(1,ISAME(I))
         ELSE
c           IF (I.NE.1) CALL QUIT ('All ranges should be equal')
            CALL MEMGET('REAL',KQ(1,I),NQT(I),WORK,KFREE,LFREE)
         ENDIF
         KQ(2,I) = KQ(1,I) + NQ(1,I)
C
      ENDDO
C
C     Selected eigenvalues
C
      DO I = 1,4
         CALL MEMGET('REAL',KE(1,I),0,WORK,KFREE,LFREE)
         KE(2,I) = KE(1,I)
      ENDDO
C
C     Allocate memory for the integer array with information about the
C     spinors.
C
      DO I = 1, 4
         DO IFRP = 1, NFSYM
            CALL MEMGET('INTE',KIBE(IFRP,I),NSTR(IFRP,0,I),
     &           WORK,KFREE,LFREE)
         ENDDO
      ENDDO
C
C     Select the set that we need
C
C     TODO: selection IBEIG not implemented yet
C
      CALL MEMGET('INTE',KIBEIG,NORBT,WORK,KFREE,LFREE)
      CALL IZERO(WORK(KIBEIG),NORBT)
C
      DO I = 1, 4 
        DO IFRP = 1, NFSYM
          IF(NSTR(IFRP,0,I).GT.0) THEN
             CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   NSTR(IFRP,0,I),WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
             CALL SELIBEIG(WORK(KIBEIG),IORB(IFRP),IFRP,
     &                     WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                     NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          ENDIF
        ENDDO
      ENDDO
C
C     Throw away IBEIG again 
C
      CALL MEMREL('RTRACTL1.CMO',WORK,KFRSAV,KIBEIG,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      CALL PAMTR1(WORK,KFREE,LFREE,IPRTRA,KQ,KE,KIBE,
     &     NDMOQR,ICMOQR,NSTR,.FALSE.,DUM,
     &     TRA_MSOUT,TRA_ANTIS,ITRA_INTFL4)
C     
C     Print timing information
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
C
      IMINS  = INT(WALL)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(WALL) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(/A,I5.2,A,I2.2,A,I2.2)')
     &      ' (RTRACTL1) Total wall time used :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CPU    = CPU2 - CPU1
      IMINS  = INT(CPU)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(CPU) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &      ' (RTRACTL1) Total CPU  time used :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)')
     &     ' (RTRACTL1) Transformation ended at : ',DAYTID
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RTRACTL1')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2 */
      SUBROUTINE RGETH2(FQ,H2AC,PV,WORK,LWORK)
C***********************************************************************
C
C     Calculate FQ and get H2AC (two-electron integrals
C     with active indices).
C
C     Input:
C        PV   - two-electron density matrix
C
C     Output:
C        FQ   -   aux. FQ fock matrix
C        H2AC -   two-electron integrals with active indices.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      DIMENSION FQ(*)
      DIMENSION H2AC(*)
      DIMENSION PV(*)
      DIMENSION WORK(*)
C
C
#include "memint.h"
C
      CALL QENTER('RGETH2')
      CALL RTKTIME(.TRUE.,3)
C
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' .OR.
     &     OPT_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
C
C        Two electron integrals
C        ----------------------
C
C        Get FQ and active two-electron integrals.
C
         CALL RGETH21(.TRUE.,.TRUE.,PV,FQ,H2AC,
     &        WORK(KFREE),LFREE,IPROPT)
C
      ELSE
         WRITE(LUPRI,'(/,A,/,2A)')
     $        '*** ERROR in RGETH2 ***',
     $        'FQ and H2AC not implemented for CI program ',
     $        OPT_CIPROGRAM
         CALL QUIT('*** ERROR in RGETH2 ***')
      END IF
      CALL RTKTIME(.FALSE.,3)
      CALL QEXIT('RGETH2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth21 */
      SUBROUTINE RGETH21(DOFQ,DOH2AC,PV,FQ,H2AC,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Calculate FQ and H2AC.
C
C     Input:
C        DOFQ  - calculate the auxiliary Fock matrix FQ
C        DOH2AC - collect active two-electron integrals.
C        PV - two-electron density matrix
C
C     Output:
C        FQ - auxiliary Fock matrix
C        (H2AC written to disk)
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : L. Visscher - Mar 15 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbidx.h"
#include "dcbham.h"
C
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION WORK(*)
      LOGICAL DOFQ, DOH2AC
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION IKLB(2)
C     strings
      CHARACTER CDATE*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
      logical fopen
C
#include "memint.h"
C
      CALL QENTER('RGETH21')
C
C     Read information from integral file
C     -----------------------------------
C
      LUINT = 1 
      OPEN (LUINT, FILE = 'MDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      IF ( IPRINT .GE. 5 ) WRITE(LUPRI,'(4A)')
     &     ' (RGETH2) Time and date read from MDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,2*NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
C     Boson irrep for IJ indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXB12))
C     Boson irrep for KL indices
      CALL READI(LUINT,NFPCK34T*2,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C     --------------------
C
      IIKLB = 0
      LUINT = 56
C     LUINT = 56 + IREPIJ, IREPIJ = 1,NFSYM are actually used
      DO IREPIJ = 1, NFSYM
         IKLB(IREPIJ) = IIKLB
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
         IIKLB = IIKLB + NFPCK34 (IREPIJ)
      END DO
C
C     Allocate memory for integrals (H2XY)
C     
      CALL MEMGET('REAL',KHMAT,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XY,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2YX,N2ORBX*NCLASS,WORK,KFREE,LFREE)
c
#ifdef UNDEF
      write(6,*) '(RGETH2) IDXKR12'
      do i = 1, nfpck12t
         write(6,*) I,
     &        IGETINT(WORK(KDXKR12),1+(I-1)*2),
     &        IGETINT(WORK(KDXKR12),I*2)
      end do
#endif
c
C
C     Read integrals (x>y)
C     --------------------
C
C     zero arrays
C
      IF (DOFQ) CALL DZERO(FQ,NASHT*NORBT*NZ)
      IF (DOH2AC) CALL DZERO(H2AC,NASHT*NASHT*NNASHX*NZ*3)
C
C     Loop over records in integrals files.
C
      NXY = 0
      DO 10 IREPIJ = 1, NFSYM
         IREC = 0
         DO 20 IXY = 1, NFPCK12(IREPIJ)
            NXY = NXY + 1
            IREC = IREC + 1
            CALL GTXY(NXY,IX,IY,IXU,IYU,WORK(KDXKR12))
            IREPIJB = IGETINT(WORK(KDXB12),NXY*2-1)
#ifdef UNDEF
            write(6,*) 'x & y: ',ix,iy
            write(6,*) 'ative : ',ixu,iyu
            write(6,*) ' '
#endif
C           We only need active-active distributions.
C
            IF (IXU .LT. 0 .OR. IYU .LT. 0) GOTO 20
C
C           ...also skip if X<Y
C
            IF (IXU .LT. IYU) GOTO 20
C
            IUU = MAX(IXU,IYU) * ( MAX(IXU,IYU) - 1 ) / 2
     &           + MIN(IXU,IYU)
C
C           Read integrals
C           
            NKL = NFPCK34(IREPIJ)
C
            READ (LUINT + IREPIJ, REC=IREC,
     &           ERR=8000, IOSTAT=IOSVAL)
     &           (WORK(KH2XY+I-1),I=1,NKL*NCLASS) 
#ifdef UNDEF
            ii = 0
            do i = 1, 4 * nz
               write(6,*) 'class ',i,' before transformation'
               call output(work(kh2xy+ii),1,norbt,1,norbt,
     &              norbt,norbt,1,lupri)
               ii = ii + nkl
            end do
#endif
C
C           Reorder integrals to full matrix
C           --------------------------------
C
            CALL DZERO(WORK(KH2YX),N2ORBX*NCLASS)
C
            IREPABR = 0
            DO IREPAB = 0, NBSYM - 1
               IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
C            
C                    Transform symmetry packed
C                    integrals to full matrix.
C                    (take care of diff. integral transformations etc.)
C
                     CALL INTB2X(WORK(KH2XY),IZ2,IREPABR,IREPAB,
     &                    WORK(KH2YX),IREPIJ,NKL)
C                    CALL INTB2X(HT,IZ2,IREPABR,IREPAB,HX,IREPIJ,NKL)
C
C                    Transform boson-symmetry array to
C                    full
C
C                    (FIXME: write me)
C                    
                     IF (SPINFR) THEN
                        write(6,*) 'SPINFR not implemented in RGETH2'
                        CALL QUIT('SPINFR not implemented in RGETH2')
                     END IF
C
                  END DO
               END IF
            END DO
C
C           Transform to Dirac (NZ,3) format
C           --------------------------------
C
            CALL CL2NZ3(WORK(KH2YX),WORK(KHMAT),N2ORBX,IREPIJ,
     &           IREPIJB,WORK(KDXB34),IKLB)
C           CALL CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C
C           Get (gu) and (ug) blocks + H2AC
C           -------------------------------
C
            CALL X2GU(WORK(KHMAT),WORK(KH2XY),WORK(KH2YX),
     &           H2AC,IUU,IXU .EQ. IYU, DOFQ, DOH2AC)
C           CALL X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C
#ifdef UNDEF
            do i = 1, 3 * nz
               write(6,*) 'xy matrix no ',i
               call output(work(kh2xy + (i-1)*norbt*nasht),
     &              1,norbt,1,nasht,norbt,nasht,1,lupri)
            end do
            if ( ix .ne. iy) then
               do i = 1, 3 * nz
                  write(6,*) 'yx matrix no ',i
                  call output(work(kh2yx + (i-1)*norbt*nasht),
     &                 1,norbt,1,nasht,norbt,nasht,1,lupri)
               end do
            end if
#endif
C
C           Do the gazillion matrix multiplies to construct
C           the auxililary Fock matrix FQ.
C           
            IF (DOFQ) THEN
               CALL RADDFQ(FQ,WORK(KH2XY),WORK(KH2YX),PV,
     &              IUU,IXU,IYU,IREPIJ,IPRINT)
            END IF
C
 20      CONTINUE
 10   CONTINUE
C
C     Close integral files
C     --------------------
C
      DO IFSYM = 1, NFSYM
         CLOSE (LUINT + IFSYM, STATUS = 'KEEP')
      END DO
C
      IF (IPRINT .GE. 10) THEN
         IF (DOH2AC) THEN
            WRITE(LUPRI,'(A)') ' (RGETH21) H2AC matrix'
            CALL PRDNZ3(H2AC,NASHT,NNASHX,NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IF (DOFQ) THEN
            WRITE(LUPRI,'(A)') ' (RGETH21) FQ matrix'
            CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END IF
C
      CALL QEXIT('RGETH21')
      RETURN
C
C     Read errors
C
 8000 CONTINUE
      WRITE(LUPRI,'(/A)')
     &     '*** ERROR in RGETH2 *** Read error'
 8100 CONTINUE
      WRITE(LUPRI,'(1X,A,I4)')
     &     'IREPIJ = ',IREPIJ,
     &     'IREC   = ',IREC,
     &     'IOSTAT = ',IOSVAL,
     &     'I      = ',I,
     &     'IX     = ',IX,
     &     'IY     = ',IY,
     &     'NXY    = ',NXY
      CALL QUIT('*** ERROR in RGTETH2 *** Error reading integral file')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck raddfq */
      SUBROUTINE RADDFQ(FQ,H2XY,H2YX,PV,NXY,IX,IY,IREPIJ,IPRINT)
C***********************************************************************
C
C     Add up the infinite number of terms that constitutes the FQ fock
C     matrix.
C
C     For proper documentation see Ph.D. thesis of Joern Thyssen where
C     all the formulae below should be documented :-)
C
C     Input:
C        IREPIJ: symmetry of integrals (currently unused)
C        H2XY, H2YX: (p,v) and (v,p)^T integrals for compound index XY
C        NXY   : the XY compound index
C        IX,IY : x and y
C 
C     Output:
C        FQ    : the aux. Fock matrix
C
C     Written by J. Thyssen - Nov 16 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D1 = 1.0D00 )
C
#include "dgroup.h"     
#include "dcborb.h" 
#include "dcbopt.h" 
C
      DIMENSION H2XY(NORBT,NASHT,NZ,3)
      DIMENSION H2YX(NORBT,NASHT,NZ,3)
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Output from RADDFQ',-1)
         WRITE(LUPRI,'(A,I3)') 'NXY       = ',NXY
         WRITE(LUPRI,'(A,I3)') 'IX        = ',IX
         WRITE(LUPRI,'(A,I3)') 'IY        = ',IY
         WRITE(LUPRI,'(A,I3)') 'IREPIJ    = ',IREPIJ
         WRITE(LUPRI,'(/A)') ' FQ on entry'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
c        IF (IX .EQ. 2 .AND. IY .EQ. 2) THEN
            NPVCOL = NASHT * NNASHX
            DO I = 1, 3
               WRITE(LUPRI,'(A,I2,A)') ' PV density matrix',I,' of 3'
               CALL PRQMAT(PV(1,1,NXY,1,I),NASHT,NASHT,NASHT,NPVCOL,
     &              NZ,IPQTOQ(1,0),LUPRI)
            END DO
            DO I = 1, 3
               WRITE(LUPRI,'(A,I2,A)') ' H2XY matrix',I,' of 3'
               CALL PRQMAT(H2XY(1,1,1,I),NORBT,NASHT,NORBT,NASHT,
     &              NZ,IPQTOQ(1,0),LUPRI)
            END DO
            IF (IX .NE. IY) THEN
               DO I = 1, 3
                  WRITE(LUPRI,'(A,I2,A)') ' H2YX matrix',I,' of 3'
                  CALL PRQMAT(H2YX(1,1,1,I),NORBT,NASHT,NORBT,NASHT,
     &                 NZ,IPQTOQ(1,0),LUPRI)
               END DO
            END IF
c        END IF
      END IF
C
C
C     ***********************
C     *** Real part of FQ ***
C     ***********************
C
C
C     Term 1 (real part)
C     ------------------
C
C
C     \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;1,1}
C      
      CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &     H2XY(1,1,1,1),NORBT,
     &     PV(1,1,NXY,1,1),NASHT,
     &     D1,FQ(1,1,1),NORBT)
C
      IF (IX .NE. IY) THEN
C
C        \sum_{x>y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;1,1} 
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,D1,
     &        H2YX(1,1,1,1),NORBT,
     &        PV(1,1,NXY,1,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
C
C     Term 2 (real part)
C     ------------------
C
C
      IF (NZ .GE. 2) THEN
C
C        - \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;2,1}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,2,1),NORBT,
     &        PV(1,1,NXY,2,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
         IF (IX .NE. IY) THEN
C
C           - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;2,1}
C           
            CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &           H2YX(1,1,2,1),NORBT,
     &           PV(1,1,NXY,2,1),NASHT,
     &           D1,FQ(1,1,1),NORBT)
C
         END IF
C
      END IF
C
C
C     Term 3 and 4 (real part)
C     ------------------------
C
C
      IF (NZ .GE. 4 .AND. IX .NE. IY) THEN
C
C        \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;3,1}
C        
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,3,1),NORBT,
     &        PV(1,1,NXY,3,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C        
C
C        - \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;4,1}
C        
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,4,1),NORBT,
     &        PV(1,1,NXY,4,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C        
      END IF
C
C
C     Term 5 (real part)
C     ------------------
C
C
      IF (IX .NE. IY) THEN
C        
C        \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;1,2}
C         
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,1,2),NORBT,
     &        PV(1,1,NXY,1,2),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
C
C     Term 6 (real part)
C     ------------------
C
C
      IF (NZ .GE. 2 .AND. IX .NE. IY) THEN
C        
C        - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;2,2}
C         
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,2,2),NORBT,
     &        PV(1,1,NXY,2,2),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
      IF (NZ .GE. 4) THEN
C
C
C        Term 7 (real part)
C        ------------------
C
C        
C        \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;3,2}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,3,2),NORBT,
     &        PV(1,1,NXY,3,2),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C        
         IF (IX .NE. IY) THEN
C           
C           \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;3,3} 
C           
            CALL DGEMM('N','N',NORBT,NASHT,NASHT,D1,
     &           H2YX(1,1,3,3),NORBT,
     &           PV(1,1,NXY,3,3),NASHT,
     &           D1,FQ(1,1,1),NORBT)
C           
         END IF
C
C
C        Term 8 (real part)
C        ------------------
C
C        
C        - \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;4,2}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,4,2),NORBT,
     &        PV(1,1,NXY,4,2),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C        
         IF (IX .NE. IY) THEN
C           
C           - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;4,3} 
C           
            CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &           H2YX(1,1,4,3),NORBT,
     &           PV(1,1,NXY,4,3),NASHT,
     &           D1,FQ(1,1,1),NORBT)
C           
         END IF
C
      END IF
C
C
C     Term 9 (real part)
C     ------------------
C
C
      IF (IX .NE. IY) THEN
C
C        \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;1,3}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,1,3),NORBT,
     &        PV(1,1,NXY,1,3),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
C
C     Term 10 (real part)
C     ------------------
C
C
      IF (NZ .GE. 2 .AND. IX .NE. IY) THEN
C
C        - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;2,3}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,2,3),NORBT,
     &        PV(1,1,NXY,2,3),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
C
      IF (NZ .GE. 4 .AND. IX .NE. IY) THEN
C
C
C        Term 11 (real part)
C        ------------------
C
C
C        \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;3,1}
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,D1,
     &        H2YX(1,1,3,1),NORBT,
     &        PV(1,1,NXY,3,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
C
C
C        Term 12 (real part)
C        ------------------
C
C
C        - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;4,1}
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &        H2YX(1,1,4,1),NORBT,
     &        PV(1,1,NXY,4,1),NASHT,
     &        D1,FQ(1,1,1),NORBT)
C
      END IF
C
C     Skip remaining 36 terms if NZ = 1.
C
      IF (NZ .EQ. 1) GOTO 100
C
C
C     *************************
C     *** i-imag part of FQ ***
C     *************************
C
C
C     Term 1 (i-imag part)
C     --------------------
C
C
C     \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{1,1} P_{qv,xy;2,1}
C      
      CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &     H2XY(1,1,1,1),NORBT,
     &     PV(1,1,NXY,2,1),NASHT,
     &     D1,FQ(1,1,2),NORBT)
C
      IF (IX .NE. IY) THEN
C
C        - \sum_{x>y} \sum_v \Twoint{vp|xy}_{1,1} P_{vq,xy;2,1} 
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &        H2YX(1,1,1,1),NORBT,
     &        PV(1,1,NXY,2,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
      END IF
C
C
C     Term 2 (i-imag part)
C     --------------------
C
C
C     + \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{2,1} P_{qv,xy;1,1}
C
      CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &     H2XY(1,1,2,1),NORBT,
     &     PV(1,1,NXY,1,1),NASHT,
     &     D1,FQ(1,1,2),NORBT)
C     
      IF (IX .NE. IY) THEN
C        
C        - \sum_{x > y} \sum_v \Twoint{vp|xy}_{2,1} P_{vq,xy;1,1}
C        
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &        H2YX(1,1,2,1),NORBT,
     &        PV(1,1,NXY,1,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C        
      END IF
C
C
C     Term 3 and 4 (i-imag part)
C     --------------------------
C
C
      IF (NZ .GE. 4 .AND. IX .NE. IY) THEN
C
C        \sum_{x > y} \sum_v \Twoint{pv|xy}_{3,1} P_{qv,xy;4,1}
C        
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,3,1),NORBT,
     &        PV(1,1,NXY,4,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C        
C
C        \sum_{x > y} \sum_v \Twoint{pv|xy}_{4,1} P_{qv,xy;3,1}
C        
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,4,1),NORBT,
     &        PV(1,1,NXY,3,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C        
      END IF
C
C
      IF (IX .NE. IY) THEN
C        
C
C        Term 5 (i-imag part)
C        ------------------
C
C
C        \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,2} P_{qv,xy;2,2}
C         
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,1,2),NORBT,
     &        PV(1,1,NXY,2,2),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
C
C        Term 6 (i-imag part)
C        --------------------
C
C
C        \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,2} P_{qv,xy;1,2}
C         
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,2,2),NORBT,
     &        PV(1,1,NXY,1,2),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
      END IF
C
      IF (NZ .GE. 4) THEN
C
C
C        Term 7 (i-imag part)
C        --------------------
C
C        
C        \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{3,2} P_{qv,xy;4,2}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,3,2),NORBT,
     &        PV(1,1,NXY,4,2),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C        
         IF (IX .NE. IY) THEN
C           
C           - \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,3} P_{vq,xy;4,3} 
C           
            CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &           H2YX(1,1,3,3),NORBT,
     &           PV(1,1,NXY,4,3),NASHT,
     &           D1,FQ(1,1,2),NORBT)
C           
         END IF
C
C
C        Term 8 (i-imag part)
C        --------------------
C
C        
C        \sum_{x \geq y} \sum_v \Twoint{pv|xy}_{4,2} P_{qv,xy;3,2}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,D1,
     &        H2XY(1,1,4,2),NORBT,
     &        PV(1,1,NXY,3,2),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C        
         IF (IX .NE. IY) THEN
C           
C           - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,3} P_{vq,xy;3,3} 
C           
            CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &           H2YX(1,1,4,3),NORBT,
     &           PV(1,1,NXY,3,3),NASHT,
     &           D1,FQ(1,1,2),NORBT)
C           
         END IF
C
      END IF
C
C
      IF (IX .NE. IY) THEN
C
C
C        Term 9 (i-imag part)
C        --------------------
C
C
C        - \sum_{x>y} \sum_v \Twoint{pv|xy}_{1,3} P_{qv,xy;2,3}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,1,3),NORBT,
     &        PV(1,1,NXY,2,3),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
C        
C        Term 10 (i-imag part)
C        ---------------------
C
C
C        - \sum_{x>y} \sum_v \Twoint{pv|xy}_{2,3} P_{qv,xy;1,3}
C
         CALL DGEMM('N','T',NORBT,NASHT,NASHT,DM1,
     &        H2XY(1,1,2,3),NORBT,
     &        PV(1,1,NXY,1,3),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
      END IF
C
C
      IF (NZ .GE. 4 .AND. IX .NE. IY) THEN
C
C
C        Term 11 (i-imag part)
C        ---------------------
C
C
C        \sum_{x>y} \sum_v \Twoint{vp|xy}_{3,1} P_{vq,xy;4,1}
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,D1,
     &        H2YX(1,1,3,1),NORBT,
     &        PV(1,1,NXY,4,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
C
C
C        Term 12 (i-imag part)
C        ---------------------
C
C
C        - \sum_{x>y} \sum_v \Twoint{vp|xy}_{4,1} P_{vq,xy;3,1}
C
         CALL DGEMM('N','N',NORBT,NASHT,NASHT,DM1,
     &        H2YX(1,1,4,1),NORBT,
     &        PV(1,1,NXY,3,1),NASHT,
     &        D1,FQ(1,1,2),NORBT)
C
      END IF
C
C     Skip remaining 24 terms if NZ = 2.
C
      IF (NZ .EQ. 2) GOTO 100
C
C
      call quit('*** ERROR in RADDFQ *** NZ = 4 not implemented')
C
C     Output section
C     --------------
C
C
 100  CONTINUE
      IF (IPRINT .GE. 20) THEN
         WRITE(LUPRI,'(A)') 'FQ on exit'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tknanobl */
      SUBROUTINE TKNONABL(H2FULL,NKL,H2XY,H2YX,H2AC,NXY,
     &     IREPIJ,NSTR,ONLYXY,DOFQ,DOH2AC)
C***********************************************************************
C
C     Transform H2FULL from packed format to H2(NORBT,NASHT,*).
C
C     Input:
C        IREPIJ: symmetry of H2FULL
C        H2FULL: symmetry packed:
C                   IREPIJ = gerade:
C                       index3(gerade) * index4(gerade) + 
C                       index3(ungerade) * index4(ungerade)
C                   IREPIJ = ungerade:
C                       index3(gerade) * index4(ungerade) + 
C                       index3(ungerade) * index4(gerade)
C        NKL   : leading dimension of H2FULL
C        ONLYXY: calculate only H2XY (not H2YX)
C 
C     Output:
C        H2XY   : integrals in (NORBT,NASHT,*) format.
C        H2YX   : integrals in (NORBT,NASHT,*) format.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dgroup.h"     
#include "dcborb.h" 
#include "dcbopt.h" 
#include "dcbidx.h"
C
      DIMENSION H2FULL(NKL,NZ,3)
      DIMENSION H2XY(NORBT,NASHT,NZ,3)
      DIMENSION H2YX(NORBT,NASHT,NZ,3)
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3)
      LOGICAL   ONLYXY, DOH2AC, DOFQ
      DIMENSION NSTR(2,0:2,4)
C
      IF (DOFQ) THEN
         CALL DZERO(H2XY, NORBT * NASHT * NZ * 3)
         IF (.NOT. ONLYXY) CALL DZERO(H2YX, NORBT * NASHT * NZ * 3)
      END IF
C
      IF (JTRLVL .EQ. 0 .AND. DOFQ) THEN
         WRITE(LUPRI,'(A,/,A,I3)')
     &        '*** ERROR in TKNANOBL ***',
     &        'Illegal transformation level: ',JTRLVL
         CALL QUIT('*** ERROR in TKNANOBL ***')
      ELSE IF (JTRLVL .EQ. 2 ) THEN
         WRITE(LUPRI,'(A,/,A,I3)')
     &        '*** ERROR in TKNANOBL ***',
     &        'Not implemented; transformation level: ',JTRLVL
         CALL QUIT('*** ERROR in TKNANOBL ***')
      ELSE 
C
C        (uu|uu), (gg|gu) or (gg|gg) integrals transformation.
C
         IPV = 0
         IVV = 0
         DO IREPV = 1, NFSYM
            IREPP = MOD (IREPV + IREPIJ, 2)  + 1
            DO IV = 1, NSTR(IREPV,0,4)
               IVV = IVV + 1
               DO IP = 1, NSTR(IREPP,0,3)
                  IF (IREPP .EQ. 1) THEN
                     IPP = IP
                  ELSE
                     IPP = IP + NSTR(1,0,3)
                  END IF
                  IPV = IPV + 1
C
C                 Get active indices from general/active indices
C
                  IF (.NOT. OPT_NOPFQ) THEN
                     IPVG = IPV
                     IVVG = IVV
                     IPPG = IPP
                  ELSE
                     IPVG = IDXE2G(IPV)
                     IVVG = IDXE2G(IVV)
                     IPPG = IDXE2G(IPP)
                  END IF
C
                  IF (JTRLVL .EQ. 0) THEN
                     IACP = IPP
                     IACV = IVV
                  ELSE IF (JTRLVL .EQ. 1) THEN
                     IACP = IDXG2U(IPPG)
                     IACV = IVV
                  ELSE 
                     IACP = IDXG2U(IPPG)
                     IACV = IDXG2U(IVVG)
                  END IF
C
C                 Elements for H2AC
C
                  IF (DOH2AC) THEN
                     IF (IACP .GT. 0 .AND. IACV .GT. 0) THEN
                        DO I3 = 1, 3
                           DO IZ = 1, NZ
                              H2AC(IACP,IACV,NXY,IZ,I3) = 
     &                             H2FULL(IPV,IZ,I3)
                           END DO
                        END DO
                     END IF
                  END IF
C
C                 Elements for FQ
C
                  IF (DOFQ) THEN
C
C                    (p,v) elements
C                     
                     IF (IACV .GT. 0) THEN
                        DO I3 = 1, 3
                           DO IZ = 1, NZ
                              H2XY(IPPG,IACV,IZ,I3) = 
     &                             H2FULL(IPV,IZ,I3)
                           END DO
                        END DO
                     END IF
C
C                    (v,p) elements
C                    
                     IF (IACP .GT. 0 .AND. .NOT. ONLYXY) THEN
                        DO I3 = 1, 3
                           DO IZ = 1, NZ
                              H2YX(IVVG,IACP,IZ,I3) = 
     &                             H2FULL(IPV,IZ,I3)
                           END DO
                        END DO
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END IF
C
      IF (IPROPT .GE. 30 .AND. DOH2AC) THEN
         WRITE(LUPRI,'(A,I3,A)') ' (TKNANOBL) H2AC(*,*,',NXY,'): '
         DO I3 = 1, 3
            DO IZ = 1, NZ
               WRITE(LUPRI,'(A,I3,A,I3,A)')
     &              ' Matrix no. (',IZ,',',I3,'): '
               CALL OUTPUT(H2AC(1,1,NXY,IZ,I3),1,NASHT,1,NASHT,
     &              NASHT,NASHT,1,LUPRI)
            END DO
         END DO
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndrecno */
      SUBROUTINE FNDRECNO(IX,IY,IDXKR12,NFPCK12,NFPCK12T,
     &     IREPIJ,IREC)
C***********************************************************************
C
C     Find record number for distribution (*,*,IX,IY) integrals.
C
C
C     Input:
C        IX, IY  : active orbital indices
C        IDXKR12 : gives orbitals for compound index NXY
C        NFPCK12 : number of compund indices for each symmetry
C        NFPCK12T: sum over NFPCK12
C       
C     Output:
C        IREPIJ  : irrep for compound (IX,IY)
C        IREC    : record number
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbidx.h"
C
      DIMENSION NFPCK12(2)
      DIMENSION IDXKR12(2,*)
C
      IF (JTRLVL .GE. 0 .AND. JTRLVL .LE. 2) THEN
C        (uu|uu), (gu|uu), or (gg|uu) integrals
         IGX = IX
         IGY = IY
      ELSE IF (JTRLVL .EQ. 3) THEN
C        (gg|gu) integrals
         IGX = IDXU2G(IX)
         IGY = IDXU2G(IY)
      ELSE
C        (gg|gg) integrals
         IGX = IDXU2G(IX)
         IGY = IDXU2G(IY)
      END IF
C
      NIJ = 0
      DO IFSYM = 1, NFSYM
         IRECNO = 0
         DO IJ = 1, NFPCK12(IFSYM)
            NIJ = NIJ + 1
            IRECNO = IRECNO + 1
            I = IDXKR12(1,NIJ)
            J = IDXKR12(2,NIJ)
            IF ( I .EQ. IGX .AND. J .EQ. IGY) THEN
               IREPIJ = IFSYM
               IREC = IRECNO
               RETURN
            END IF
         END DO
      END DO
      WRITE(LUPRI,'(/A/A)')
     &     '*** ERROR in FNDRECNO ***',
     &     'Internal error...'
      WRITE(LUPRI,'(4(A,I4))')
     &     ' IX = ',IX,'; IY = ',IY,'; IGX = ',IGX,'; IGY = ',IGY
      WRITE(LUPRI,'(/A)') ' Dump of IDXKR12: '
      WRITE(LUPRI,'(3I5)') (I,IDXKR12(1,I),IDXKR12(2,I),I=1,NFPCK12T)
      CALL QUIT('*** ERROR in FNDRECNO *** Internal error')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igetint */
      INTEGER FUNCTION IGETINT(IX,I)
C***********************************************************************
C
C     Return element I of array IX
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION IX(*)
C
      IGETINT = IX(I)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2tx */
      SUBROUTINE RGETH2TX(FQX,FQT,H2ACX,PV,PVT,
     &     BEVEC,BPVEC,JXOPE,JXOPP,
     &     IBEVC,IBPVC,IBTYP,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Calculate FQX and get H2ACX, i.e. one-index transformed
C     FQ and H2AC. Also calculate FQT.
C     
C
C     Input:
C        PV   - two-electron density matrix
C        PVT  - symm. two-electron transition density matrix
C
C     Output:
C        FQX  -   one-index transformed FQ fock matrix
C        FQT  -   transition density FQ fock matrix
C        H2ACX-   two-electron integrals with active indices.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
C
      DIMENSION FQX(*),FQT(*), H2ACX(*)
      DIMENSION PV(*),PVT(*)
      DIMENSION BEVEC(*),BPVEC(*)
      DIMENSION IBEVC(*),IBPVC(*),IBTYP(*)
      DIMENSION JXOPE(*),JXOPP(*)
      DIMENSION WORK(*)
C
C
#include "memint.h"
C
      CALL QENTER('RGETH2TX')
      CALL RTKTIME(.TRUE.,4)
C
      IF ( XRS_CIPROGRAM .EQ. 'GOSCI' .OR.
     &     XRS_CIPROGRAM .EQ. 'LUCIAREL' ) THEN
C
C        Two electron integrals
C        ----------------------
C
C        Get FQ and active two-electron integrals.
C
         CALL RGETH2TX1(.TRUE.,.TRUE.,PV,PVT,FQX,FQT,H2ACX,
     &        BEVEC,BPVEC,JXOPE,JXOPP,IBEVC,IBPVC,IBTYP,IPRINT,
     &        WORK(KFREE),LFREE)
C
      ELSE
         WRITE(LUPRI,'(/,A,/,2A)')
     $        '*** ERROR in RGETH2TX ***',
     $        'FQX, FQT, and H2ACX not implemented for CI program ',
     $        XRS_CIPROGRAM
         CALL QUIT('*** ERROR in RGETH2TX ***')
      END IF
      CALL RTKTIME(.FALSE.,4)
      CALL QEXIT('RGETH2TX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2tx1 */
      SUBROUTINE RGETH2TX1(DOFQ,DOH2AC,PV,PVT,FQX,FQT,H2ACX,
     &     BEVEC,BPVEC,JXOPE,JXOPP,
     &     IBEVC,IBPVC,IBTYP,IPRINT,WORK,LWORK)
C***********************************************************************
C
C
C
C     Input:
C        DOFQ  - calculate FQ matrices
C        DOH2AC - collect active two-electron integrals.
C        PV - two-electron density matrix
C        PVT - symm. transition two-electron density matrix
C
C     Output:
C        FQX - auxiliary Fock matrix from one-index 
C              transformed integrals.
C        FQT - auxiliary Fock matrix from symm. transition 
C              density matrix.
C        H2ACX - one-index transformed integrals.
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbxrs.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbidx.h"
#include "dcbham.h"
C
      DIMENSION FQX(NORBT,NASHT,NZ,*)
      DIMENSION FQT(NORBT,NASHT,NZ,*)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION PVT(NASHT,NASHT,NNASHX,NZ,3,*)
      DIMENSION H2ACX(NASHT,NASHT,NNASHX,NZ,3,*)
      DIMENSION JXOPE(*),JXOPP(*)
      DIMENSION BEVEC(*),BPVEC(*)
      DIMENSION IBEVC(*),IBPVC(*),IBTYP(2,*)
      DIMENSION WORK(*)
      LOGICAL DOFQ, DOH2AC
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION IKLB(2)
C     strings
      CHARACTER CDATE*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
C     logicals
      LOGICAL DOFQT, DOFQX, DOH2ACX
      logical fopen
C
#include "memint.h"
C
      CALL QENTER('RGETH2TX1')
C
      DOFQT = DOFQ .AND. NCSIM .GT. 0
      DOFQX = DOFQ .AND. NOSIM .GT. 0 .AND. (.NOT. XRS_NOFQX)
      DOH2ACX = DOH2AC .AND. NOSIM .GT. 0
C
C
      IF ( IPRINT .GE. 3 ) THEN
         WRITE(LUPRI,'(/2(A,L1))')
     &        ' (RGETH2TX1) DOFQ = ',DOFQ,' and DOH2AC = ',DOH2AC
         WRITE(LUPRI,'(3(A,L1)/)')
     &        '         ==> DOFQX = ',DOFQX,', DOFQT = ',DOFQT,
     &        ', and DOH2ACX = ',DOH2ACX
      END IF
C
C
C     debug debug
c     calculate FQT with RGET2 module
c     CALL RGETH2(FQT,H2ACX,PVT,WORK(KFREE),LFREE)
c     call quit('ralle ralle')
C
C     *******************************************
C     *** Read information from integral file ***
C     *******************************************
C
      LUINT = 1 
      OPEN (LUINT, FILE = 'MDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      IF ( IPRINT .GE. 3 ) WRITE(LUPRI,'(/4A/)')
     &     ' (RGETH2TX1) time and date read from MDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,2*NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
C     Boson irrep for IJ indices
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXB12))
C     Boson irrep for KL indices
      CALL READI(LUINT,NFPCK34T*2,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C     --------------------
C
      LUINT = 56
C     LUINT = 56 + IREPIJ, IREPIJ = 1,NFSYM are actually used
      IIKLB = 0
      DO IREPIJ = 1, NFSYM
         IKLB(IREPIJ) = IIKLB
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
         IIKLB = IIKLB + NFPCK34 (IREPIJ)
      END DO
C
C     Allocate memory for integrals (H2XY)
C     
      CALL MEMGET('REAL',KHMAT,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHMA2,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XY,N2ORBX*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2YX,N2ORBX*NCLASS,WORK,KFREE,LFREE)
C
C     ****************************
C     *** Read integrals (x>y) ***
C     ****************************
C
C     zero arrays
C     -----------
C
      IF (DOFQ) THEN
         CALL DZERO(FQX,NASHT*NORBT*NZ*NOSIM)
         CALL DZERO(FQT,NASHT*NORBT*NZ*NCSIM)
      END IF
      IF (DOH2AC) CALL DZERO(H2ACX,NASHT*NASHT*NNASHX*NZ*3*NOSIM)
C
C     Loop over records on file
C     -------------------------
C
#ifdef UNDEF
      write(6,*) 'nfpck12',nfpck12,nfpck12t
      write(6,*) 'nfpck34',nfpck34,nfpck34t
      call flshfo(6)
#endif
      NXY = 0
      IUU = 0
      DO 10 IREPXY = 1, NFSYM
         IREC = 0
         DO 20 IXY = 1, NFPCK12(IREPXY)
            IREC = IREC + 1
            NXY = NXY + 1
            CALL GTXY(NXY,IX,IY,IXU,IYU,WORK(KDXKR12))
            IREPIJB = IGETINT(WORK(KDXB12),NXY*2-1)
C
C           We assume it is a (gu|gg) or (gg|gg) integral
C           transformation.
C
C           case (xy) distribution:
C
C           active-active: 
C              (1) calculate FQT
C              (2) one-index transform (pv)
C              (3) calculate (xy|p~v~) P_{qv,xy} contribution
C              (4) calculate H2ACX integrals
C              (5) calculate (pv|x~y~) P_{qv,xy} contribution
C                  (if we have active-active rotations, e.g GAS)
C
C           general-active:
C              (1) calculate (pv|x~y~) P_{qv,xy} contribution
C
C           else:
C              nothing (jump to 20)
C
#ifdef UNDEF
            write(6,*) 'ix,..',ix,iy,ixu,iyu
            call flshfo(6)
#endif
            IF (IXU .LT. 0 .AND. IYU .LT. 0) GO TO 20
C
C           If there are no orbital trial vectors we only need FQT.
C           
            IF (NOSIM .EQ. 0 .AND. (IYU .LT. 0 .OR. IXU .LT. 0))
     &           GO TO 20
C
C
C           Read integrals
C           --------------
C
            NKL = NFPCK34(IREPXY)
C
            READ (LUINT + IREPXY, REC=IREC,
     &           ERR=8000, IOSTAT=IOSVAL)
     &           (WORK(KHMA2+I-1),I=1,NKL*NCLASS) 
#ifdef UNDEF
            write(6,*) 'integrals...'
            call output(work(khma2),1,nkl,1,nclass,nkl,nclass,1,lupri)
#endif
C
C           Reorder integrals to full matrix
C           --------------------------------
C
            CALL DZERO(WORK(KHMAT),N2ORBX*NCLASS)
C
            IREPABR = 0
            DO IREPAB = 0, NBSYM - 1
               IF (IREPXY .EQ. JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
C            
C                    Transform symmetry packed
C                    integrals to full matrix.
C                    (take care of diff. integral transformations etc.)
C
                     CALL INTB2X(WORK(KHMA2),IZ2,IREPABR,IREPAB,
     &                    WORK(KHMAT),IREPXY,NKL)
C                    CALL INTB2X(HT,IZ2,IREPABR,IREPAB,HX,IREPIJ,NKL)
C
C                    Transform boson-symmetry array to
C                    full
C
C                    (FIXME: write me)
C                    
                     IF (SPINFR) THEN
                        write(6,*) 'SPINFR not implemented in RGETH2TX1'
                        CALL QUIT('SPINFR not implemented in RGETH2TX1')
                     END IF
C
                  END DO
               END IF
            END DO
C
            IF (IXU .GT. 0 .AND. IYU .GT. 0) THEN
C
C              Somehow the 4-index program doesn't have index 1 & 2 
C              triangular for (gu|gg) transformations, so to avoid
C              doing both (21) and (12) distributions:
               IF (IXU .LT. IYU) GOTO 20
C
C
C              
C              **********************************
C              *** Active-active distribution ***
C              **********************************
C
C
C              Increment active-active distribution counter
C              --------------------------------------------
C
               IUU = IUU + 1
C              nope, calculate IUU:
               IUU = MAX(IXU,IYU) * ( MAX(IXU,IYU) - 1 ) / 2
     &              + MIN(IXU,IYU)
C              
C              (1) Make FQT
C              ------------
C
               IF ( DOFQT ) THEN
C                 
C
C                 Transform to Dirac (NZ,3) format
C                 --------------------------------
C
                  CALL CL2NZ3(WORK(KHMAT),WORK(KHMA2),N2ORBX,IREPXY,
     &                 IREPIJB,WORK(KDXB34),IKLB)
C                 CALL CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C
C                 Get (gu) and (ug) blocks
C                 -------------------------------
C
                  CALL X2GU(WORK(KHMA2),WORK(KH2XY),WORK(KH2YX),
     &                 DUMMY,DUMMY,IXU .EQ. IYU, DOFQT, .FALSE.)
C                 CALL X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C
C                 Calculate FQT:
C                 --------------
C                 
C                 Loop over all conf. trial vectors and for each
C                 symm. two-electron transition density matrix do
C                 the gazillion matrix multiplies to construct
C                 the auxiliary Fock matrix FQT.
C                 
                  DO ICSIM = 1, NCSIM
                     CALL RADDFQ(FQT(1,1,1,ICSIM),
     &                    WORK(KH2XY),WORK(KH2YX),PVT(1,1,1,1,1,ICSIM),
     &                    IUU,IXU,IYU,IREPXY,IPRINT)
                  END DO
               END IF
C
               IF ( DOFQX .OR. DOH2ACX ) THEN
C                 
C                 (2) One-index transform (pv)
C                 ----------------------------
C                 
C                 The integrals are stored in the double-quaternion form.
C                 We loop over quaternion 2, and one-index transform
c                 quaternion 1.
C                 
C                 Note, that index 2 is always (gg).
C                 
                  call flshfo(6)
                  IREP = JSYMOP - 1
C                 
C                 e-e and e-p trial vectors
C                 -------------------------
C                 
                  DO I = 1, NOSIM
C                    
C                    One-index transform integrals
C                    -----------------------------
C
                     CALL RGETH2TX2(I,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                    WORK(KHMAT),WORK(KH2XY),JXOPE,JXOPP,
     &                    IREPXY,NKL,.TRUE.,IXU .EQ. IYU,
     &                    WORK(KFREE),LFREE)
C
C                    
C                    Transform to Dirac (NZ,3) format
C                    --------------------------------
C
                     CALL CL2NZ3(WORK(KH2XY),WORK(KHMA2),N2ORBX,IREPXY,
     &                    IREPIJB,WORK(KDXB34),IKLB)
C                    CALL CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C                    
C
C                    Get (gu) and (ug) blocks + H2ACX
C                    --------------------------------
C                    
                     CALL X2GU(WORK(KHMA2),WORK(KH2XY),WORK(KH2YX),
     &                    H2ACX(1,1,1,1,1,I),
     &                    IUU,IXU .EQ. IYU, DOFQX, DOH2ACX)
C                    CALL X2GU(HX,H2XY,H2YX,H2AC,IUU,ONLYXY,DOFQ,DOH2AC)
C
C                    
C                    Calculate ``matrix-multiply'' part of FQX:
C                    ------------------------------------------
C                    
                     IF ( DOFQX ) THEN
                        CALL RADDFQ(FQX(1,1,1,I),
     &                       WORK(KH2XY),WORK(KH2YX),PV,
     &                       IUU,IXU,IYU,IREPXY,IPRINT)
                     END IF
C
                  END DO
C
               END IF
            END IF
C
            IF ((IXU .GT. 0 .OR. IYU .GT. 0).AND. DOFQX) THEN
C
               IF (IXU .GT. 0 .AND. IYU .GT. 0) THEN
C
C                 We only need the active-active ddot terms
C                 for GASSCF wave functions.
C
C
                  IF (MCTYPE .NE. JGAS) GOTO 20
C
C                 Somehow the 4-index program doesn't have index 1 & 2 
C                 triangular for (gu|gg) transformations, so to avoid
C                 doing both (21) and (12) distributions:
C
                  IF (IXU .LT. IYU) GOTO 20
C
               END IF
C
C              (3) calculate (pv|x~y~) P_{qv,xy} contribution
C              ----------------------------------------------
C              
C              If MC type is GAS/RAS then we also need the 
C              pv = active-active. This code must be copied
C              up to the RADDFQ call just above.
C
C              (3.1) One-index transform integrals.
C
C              Note, that index 2 is always (gg).
C                 
               IREP = JSYMOP - 1
C                 
C              Loop over e-e and e-p trial vectors
C              -----------------------------------
C                 
               DO I = 1, NOSIM
C                 
C                 One index transform integrals
C                 -----------------------------
C
                  CALL RGETH2TX2(I,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                 WORK(KHMAT),WORK(KH2XY),JXOPE,JXOPP,
     &                 IREPXY,NKL,.FALSE.,.FALSE.,
     &                 WORK(KFREE),LFREE)
C
C                    
C                 Transform to Dirac (NZ,3) format
C                 --------------------------------
C
                  CALL CL2NZ3(WORK(KH2XY),WORK(KHMA2),N2ORBX,IREPXY,
     &                 IREPIJB,WORK(KDXB34),IKLB)
C                 CALL CL2NZ3(VCL,VNZ3,NELEM,IREPIJ,IREPIJB,IREPKLB,IKLB)
C                    
C
C                 Get (uu) block
C                 --------------
C
C                 (use the H2AC-branch in X2GU)
C                    
                  KH2VP = KH2XY
                  CALL X2UU(WORK(KHMA2),WORK(KH2VP))
C                 CALL X2UU(HX,H2XY)
C
C
C                 Add contribution to FQX
C                 -----------------------
C
                  CALL RADDFQX(FQX(1,1,1,I),
     &                 WORK(KH2VP),PV,
     &                 IX,IY,IREPXY,IPRINT)
C
               END DO
C
            END IF
C
 20      CONTINUE
 10   CONTINUE
C
C     Close integral files
C     ---------------------
C
      DO IREPIJ = 1, NFSYM
         CLOSE (LUINT+IREPIJ, STATUS='KEEP')
      END DO
C     
      CALL MEMREL('RGETH2TX1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Symmetrize H2ACX matrices
C     -------------------------
C     
C     Add (u~v~|xy) to (uv|x~y~)
C
      IF (IPRINT .GE. 20) THEN
         DO IOSIM = 1, NOSIM
            WRITE(LUPRI,'(/A,I3)')
     $           ' H2ACX matrix before SYMMH2AC no. ', IOSIM
            CALL PRDNZ3(H2ACX(1,1,1,1,1,IOSIM),NASHT,NNASHX,NZ,
     &           IPQTOQ(1,0),LUPRI)
         END DO
      END IF
      CALL SYMMH2AC(H2ACX,NOSIM,WORK(KFREE),LFREE)
C
C     Print section
C     -------------
C
      IF (IPRINT .GE. 20) THEN
         DO IOSIM = 1, NOSIM
            WRITE(LUPRI,'(/A,I3)') ' (RGETH2TX1) H2ACX matrix no. ',
     &           IOSIM
            CALL PRDNZ3(H2ACX(1,1,1,1,1,IOSIM),NASHT,NNASHX,NZ,
     &           IPQTOQ(1,0),LUPRI)
         END DO
C
         DO IOSIM = 1, NOSIM
            WRITE(LUPRI,'(/A,I3)') ' (RGETH2TX1) FQX matrix no. ',
     &           IOSIM
            CALL PRQMAT(FQX(1,1,1,IOSIM),NORBT,NASHT,NORBT,NASHT,
     &           NZ,IPQTOQ(1,0),LUPRI)
         END DO
C
         DO ICSIM = 1, NCSIM
            WRITE(LUPRI,'(/A,I3)') ' (RGETH2TX1) FQT matrix no. ',
     &           ICSIM
            CALL PRQMAT(FQT(1,1,1,ICSIM),NORBT,NASHT,NORBT,NASHT,
     &           NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
      CALL QEXIT('RGETH2TX1')
      RETURN
C
C     Read errors
C
 8000 CONTINUE
      WRITE(LUPRI,'(/A)')
     &     '*** ERROR in RGETH2 *** Read error'
 8100 CONTINUE
      WRITE(LUPRI,'(1X,A,I4)')
     &     'IREPIJ = ',IREPIJ,
     &     'IREC   = ',IREC,
     &     'IOSTAT = ',IOSVAL,
     &     'I      = ',I,
     &     'IX     = ',IX,
     &     'IY     = ',IY,
     &     'NXY    = ',NXY
      CALL QUIT('*** ERROR in RGTETH2 *** Error reading integral file')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtxy */
      SUBROUTINE GTXY(IXY,IX,IY,IXU,IYU,IDXKR12)
C***********************************************************************
C
C     Return orbital no. from compound index IXY.
C     Also, return active orbital no.s
C
C     Written by J. Thyssen - Dec 3 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
C     
      DIMENSION IDXKR12(2,*)
C
C     Get X and Y indices:
C     
      IX = IDXKR12(1,IXY)
      IY = IDXKR12(2,IXY)
C
C
      IX = IDXT2G(IX,1)
      IY = IDXT2G(IY,2)
C
C
      IXU = IDXG2U(IX)
      IYU = IDXG2U(IY)
C
c     write(6,*) 'gtxy(2): ',ixu,iyu
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck symmh2ac */
      SUBROUTINE SYMMH2AC(H2AC,N,WORK,LWORK)
C***********************************************************************
C
C     Symmetrize H2AC.
C
C     H2AC(I,J,NKL) = H2AC(I,J,NKL) + H2AC(K,L,NIJ)
C     H2AC(K,L,NIJ) = H2AC(I,J,NKL)
C
C     Written by J. Thyssen - Dec 5 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION H2AC(*)
      DIMENSION WORK(*)
C
#include "memint.h"
C
      CALL QENTER('SYMMH2AC')
      CALL MEMGET('REAL',KH2ACO,NASHT*NASHT*NNASHX*3*NZ,
     &     WORK,KFREE,LFREE)
      CALL SYMMH2AC1(H2AC,N,WORK(KH2ACO))
      CALL MEMREL('SYMMH2AC',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('SYMMH2AC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck symmh2ac1 */
      SUBROUTINE SYMMH2AC1(H2AC,N,H2ACO)
C***********************************************************************
C
C     Symmetrize H2AC.
C
C     H2AC(I,J,NKL) = H2AC(I,J,NKL) + H2AC(K,L,NIJ)
C     H2AC(K,L,NIJ) = H2AC(I,J,NKL)
C
C     Written by J. Thyssen - Dec 5 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3,*)
      DIMENSION H2ACO(NASHT,NASHT,NNASHX,NZ,3)
C
      LOGICAL SWAPIJ
C
      if (nz .eq. 4) call quit('nz=4 not implemented')
      DO IN = 1, N
         CALL DCOPY(NASHT*NASHT*NNASHX*3*NZ,H2AC(1,1,1,1,1,IN),1,
     &        H2ACO,1)
C
         NIJ = 0
         DO I = 1, NASHT
            DO J = 1, I
               NIJ = NIJ + 1
               DO L = 1, NASHT
C                 
C                 The K<L part:
C                 -------------
C                 
                  DO K = 1, L
                     NKL = L * ( L - 1 ) / 2 + K
C                    
C                    NZ = 1:
C                    -------
C                    
C                    (kl|ij) + (ij|kl) = (kl|ij) + (ji|lk)^*
C                    
                     H2AC(K,L,NIJ,1,1,IN) =
     &                    H2AC(K,L,NIJ,1,1,IN)
     &                    + H2ACO(J,I,NKL,1,1)
C                    
C                    (kL|iJ) + (iJ|kL) = (kL|iJ) - (iJ|lK)
C                    
                     H2AC(K,L,NIJ,1,2,IN) =
     &                    H2AC(K,L,NIJ,1,2,IN)
     &                    - H2ACO(I,J,NKL,1,2)
C                    
C                    (Kl|iJ) + (iJ|Kl) = (Kl|iJ) + (Ji|lK)^*
C                    
                     H2AC(K,L,NIJ,1,3,IN) =
     &                    H2AC(K,L,NIJ,1,3,IN)
     &                    + H2ACO(J,I,NKL,1,3)
C                    
C                    NZ = 2:
C                    -------
C                    
                     IF (NZ .GE. 2) THEN
C                       
C                       (kl|ij) + (ij|kl) = (kl|ij) + (ji|lk)^*
C                       
                        H2AC(K,L,NIJ,2,1,IN) =
     &                       H2AC(K,L,NIJ,2,1,IN)
     &                       - H2ACO(J,I,NKL,2,1)
C                       
C                       (kL|iJ) + (iJ|kL) = (kL|iJ) - (iJ|lK)
C                       
                        H2AC(K,L,NIJ,2,2,IN) =
     &                       H2AC(K,L,NIJ,2,2,IN)
     &                       - H2ACO(I,J,NKL,2,2)
C                       
C                       (Kl|iJ) + (iJ|Kl) = (Kl|iJ) + (Ji|lK)^*
C                       
                        H2AC(K,L,NIJ,2,3,IN) =
     &                       H2AC(K,L,NIJ,2,3,IN)
     &                       - H2ACO(J,I,NKL,2,3)
C                       
                     END IF
C                    
C                    NZ = 4:
C                    -------
C                    
                     IF (NZ .GE. 4) THEN
C                       
C                       (kl|iJ) + (iJ|kl) = (kl|iJ) + (Ji|lk)^*
C                       
                        H2AC(K,L,NIJ,3,1,IN) =
     &                       H2AC(K,L,NIJ,3,1,IN)
     &                       + H2ACO(J,I,NKL,3,3)
C                       
                        H2AC(K,L,NIJ,4,1,IN) =
     &                       H2AC(K,L,NIJ,4,1,IN)
     &                       - H2ACO(J,I,NKL,4,3)
C                       
C                       (kL|ij) + (ij|kL) = (kL|ij) - (ij|lK)
C                       
                        H2AC(K,L,NIJ,3,2,IN) =
     &                       H2AC(K,L,NIJ,3,2,IN)
     &                       - H2ACO(I,J,NKL,3,1)
C                       
                        H2AC(K,L,NIJ,4,2,IN) =
     &                       H2AC(K,L,NIJ,4,2,IN)
     &                       - H2ACO(I,J,NKL,4,1)
C                       
C                       (Kl|ij) + (ij|Kl) = (Kl|ij) - (ij|Lk)
C                       
                        H2AC(K,L,NIJ,3,3,IN) =
     &                       H2AC(K,L,NIJ,3,3,IN)
     &                       - H2ACO(I,J,NKL,3,2)
C                       
                        H2AC(K,L,NIJ,4,3,IN) =
     &                       H2AC(K,L,NIJ,4,3,IN)
     &                       - H2ACO(I,J,NKL,4,2)
C                       
                     END IF
                  END DO
C                 
C                 The K>L part:
C                 -------------
C                 
                  DO K = L + 1, NASHT
                     NKL = K * ( K - 1 ) / 2 + L
C
C                    NZ = 1:
C                    -------
C
C                    (kl|ij) + (ij|kl)
C
                     H2AC(K,L,NIJ,1,1,IN) =
     &                    H2AC(K,L,NIJ,1,1,IN) +
     &                    H2ACO(I,J,NKL,1,1)
C
C                    (kL|iJ) + (iJ|kL)
C
                     H2AC(K,L,NIJ,1,2,IN) =
     &                    H2AC(K,L,NIJ,1,2,IN) +
     &                    H2ACO(I,J,NKL,1,2)
C
C                    (Kl|iJ) + (iJ|Kl) = (Kl|iJ) + (Ij|kL)^*
C
                     H2AC(K,L,NIJ,1,3,IN) =
     &                    H2AC(K,L,NIJ,1,3,IN) +
     &                    H2ACO(I,J,NKL,1,3)
C
                     IF (NZ .GE. 2) THEN
C
C                       NZ = 2:
C                       -------
C                       
C                       (kl|ij) + (ij|kl)
C                       
                        H2AC(K,L,NIJ,2,1,IN) =
     &                       H2AC(K,L,NIJ,2,1,IN) +
     &                       H2ACO(I,J,NKL,2,1)
C                       
C                       (kL|iJ) + (iJ|kL)
C                       
                        H2AC(K,L,NIJ,2,2,IN) =
     &                       H2AC(K,L,NIJ,2,2,IN) +
     &                       H2ACO(I,J,NKL,2,2)
C                       
C                       (Kl|iJ) + (iJ|Kl) = (Kl|iJ) + (Ij|kL)^*
C                       
                        H2AC(K,L,NIJ,2,3,IN) =
     &                       H2AC(K,L,NIJ,2,3,IN) -
     &                       H2ACO(I,J,NKL,2,3)
C
                     END IF
C
                     IF (NZ .GE. 4) THEN
C
C                       NZ = 4:
C                       -------
C                       
C                       (kl|iJ) + (iJ|kl)
C                       
                        H2AC(K,L,NIJ,3,1,IN) =
     &                       H2AC(K,L,NIJ,3,1,IN) +
     &                       H2ACO(I,J,NKL,3,2)
C
                        H2AC(K,L,NIJ,4,1,IN) =
     &                       H2AC(K,L,NIJ,4,1,IN) +
     &                       H2ACO(I,J,NKL,4,2)
C                       
C                       (kL|ij) + (ij|kL)
C                       
                        H2AC(K,L,NIJ,3,2,IN) =
     &                       H2AC(K,L,NIJ,3,2,IN) +
     &                       H2ACO(I,J,NKL,3,1)
C                       
                        H2AC(K,L,NIJ,4,2,IN) =
     &                       H2AC(K,L,NIJ,4,2,IN) +
     &                       H2ACO(I,J,NKL,4,1)
C                       
C                       (Kl|ij) + (ij|Kl) = (Kl|ij) - (ji|kL)^*
C                       
                        H2AC(K,L,NIJ,3,3,IN) =
     &                       H2AC(K,L,NIJ,3,3,IN) -
     &                       H2ACO(J,I,NKL,3,1)
C
                        H2AC(K,L,NIJ,4,3,IN) =
     &                       H2AC(K,L,NIJ,4,3,IN) +
     &                       H2ACO(J,I,NKL,4,1)
C
                     END IF
C
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2tx2 */
      SUBROUTINE RGETH2TX2(ITRVEC,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &     HMAT,HMATX,JXOPE,JXOPP,IREPXY,NKL,TFLAG,ONLYXY,
     &     WORK,LWORK)
C***********************************************************************
C
C     Calculate one-index transform of the integrals in HMAT.
C
C     Input:
C        I: orbital trial vector no.
C        BEVEC, BPVEC: orbital trial vectors
C        IBEVC, IBPVC, IBTYP: usual index vectors
C        HMAT: integrals in double quaternion format
C        JXOPE, JXOPP: orbital rotations
C        TFLAG: controls which part of HMATX to generate
C               - general-active and active-general (true)
C               - active-active (false)
C
C     Output: 
C        HMATX: one-index transformed integrals in double quaternion 
C               format.
C
C     Written by J. Thyssen - Dec 6 2000
C     Last revision : 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbxrs.h"
C
      DIMENSION BEVEC(*), BPVEC(*)
      DIMENSION IBEVC(*), IBPVC(*), IBTYP(2,*)
      DIMENSION HMAT(N2ORBX,NZ,NZ,*), HMATX(N2ORBX,NZ,NZ,*)
      DIMENSION JXOPE(*), JXOPP(*)
      LOGICAL   TFLAG, ONLYXY
      DIMENSION WORK(*)
C
      DIMENSION JGENMO(2), JACTMO(2), NGENMO(2), NACTMO(2)
C
#include "ibtfun.h"
#include "memint.h"
C
      CALL QENTER('RGETH2TX2')
C
C     Initialize JGENMO etc.
C     ----------------------
C
      DO I = 1, NFSYM
         JACTMO(I) = IORB(I) + NPSH(I) + NISH(I) + 1
         NACTMO(I) = NASH(I)
         IF ( NZXOPP .GT. 0) THEN
C           ...we have positronic rotations
            IF ( NZXOPE .GT. 0) THEN
               JGENMO(I) = IORB(I) + 1
               NGENMO(I) = NORB(I)
            ELSE
               JGENMO(I) = IORB(I) + 1
               NGENMO(I) = NPSH(I) + NOCC(I)
            END IF
         ELSE
C           ...no positronic rotations
            IF ( NZXOPE .GT. 0) THEN
               JGENMO(I) = IORB(I) + NPSH(I) + 1
               NGENMO(I) = NESH(I)
            ELSE
               WRITE(LUPRI,'(/A/A/)')
     &              '*** ERROR in RGETH2TX2 ***',
     &              'no positronic nor electronic rotations...'
               CALL QUIT('*** ERROR in RGETH2TX2 ***')
            END IF
         END IF
      END DO
C
      CALL MEMGET('REAL',KBMAT,N2ORBXQ,WORK,KFREE,LFREE)
C
C     Find Hermicity of trial vector
C     ------------------------------
C
      IF (ITRVEC .LE. NESIM) THEN
         IH   = IBEVC(NERED+ITRVEC)
         IBOFF = (ITRVEC-1) * NZXOPEQ + 1
      ELSE
         IH   = IBPVC(NPRED+ITRVEC-NESIM)
         IBOFF = (ITRVEC-NESIM-1) * NZXOPPQ + 1
      END IF
      IH = IBTYP(2,IH)
      IREP = JSYMOP - 1
      NCLASS = NZ * NZ * NBSYM / NFSYM
C                    
C     Scatter BVEC into BMAT
C     ----------------------
C     
      CALL DZERO(WORK(KBMAT),N2ORBXQ)
      IF (ITRVEC .LE. NESIM) THEN
         CALL XOPSCT(BEVEC(IBOFF),JXOPE,NZXOPE,
     &        WORK(KBMAT),NORBT,NZ)
      ELSE
         CALL XOPSCT(BPVEC(IBOFF),JXOPP,NZXOPP,
     &        WORK(KBMAT),NORBT,NZ)
      END IF
      IF(IPRXRS.GE.6) THEN
         WRITE(LUPRI,'(A,I5)') 
     &        ' (RGETH2TX2): Scattered'//
     &        ' orbital trial vector no.:',ITRVEC
         CALL PRQMAT(WORK(KBMAT),NORBT,NORBT,NORBT,NORBT,
     &        NZ,IPQTOQ(1,IREP),LUPRI)
      ENDIF
C                    
C     Construct full BMAT
C     -------------------
C     
      IOFF = 0
      DO IZ = 1,NZ
         IQ = IPQTOQ(IZ,IREP)
         IS = IHQMAT(IQ,IH)
         IF(IS.EQ.1) THEN
            CALL FULMAT('A',NORBT,NORBT,WORK(KBMAT+IOFF))
         ELSEIF(IS.EQ.2) THEN
            CALL FULMAT('S',NORBT,NORBT,WORK(KBMAT+IOFF))
         ENDIF
         IOFF = IOFF + N2ORBX
      ENDDO
C                    
      IF(IPRXRS.GE.6) THEN
         WRITE(LUPRI,'(A,I5)') 
     &        ' (RGETH2TX2): Symm. scattered orbital'//
     &        ' trial vector no.:',ITRVEC
         CALL PRQMAT(WORK(KBMAT),NORBT,NORBT,NORBT,NORBT,
     &        NZ,IPQTOQ(1,IREP),LUPRI)
      ENDIF
C
C     Loop over quaternion 2
C     ----------------------
C                    
      CALL DZERO(HMATX,N2ORBX * NCLASS)
      IREPABR = 0
      DO IREPAB = 0, NBSYM - 1
         IF (IREPXY .EQ. JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1
            DO IZ2 = 1, NZ
               IQ2 = IPQTOQ(IZ2,IREPAB)
C              
C              Calculate one-index transform:
C              
C              \sum_r b_{pr} F_{rv} -  F_{pr} b_{rv}
C              
C              IREPM is the product of IREPAB and the symmetry of the
C              B operator:
C
               IREPM = IBTXOR(IREPAB,JSYMOP-1)
C
               IF (TFLAG) THEN
C                 
C                 We need both general-active block:
C                 
                  CALL RTR1H1(JGENMO,NGENMO,JACTMO,NACTMO,
     &                 HMATX(1,1,IZ2,IREPABR),IPQTOQ(1,IREPM),
     &                 WORK(KBMAT),IPQTOQ(1,JSYMOP-1),
     &                 HMAT(1,1,IZ2,IREPABR),IPQTOQ(1,IREPAB),IREPXY,
     &                 WORK,KFREE,LFREE)
C                 
C                 ...and active-general block for X .ne. Y:
C                 (zero out active-active block to avoid double-counting)
C                 
                  IF (.NOT. ONLYXY) THEN
                     CALL BLKZERO(JACTMO,NACTMO,JACTMO,NACTMO,
     &                    HMATX(1,1,IZ2,IREPABR),NORBT,NORBT,
     &                    NZ,NFSYM,IREPXY)
                     CALL RTR1H1(JACTMO,NACTMO,JGENMO,NGENMO,
     &                    HMATX(1,1,IZ2,IREPABR),IPQTOQ(1,IREPM),
     &                    WORK(KBMAT),IPQTOQ(1,JSYMOP-1),
     &                    HMAT(1,1,IZ2,IREPABR),IPQTOQ(1,IREPAB),
     &                    IREPXY,
     &                    WORK,KFREE,LFREE)
                  END IF
               ELSE
                  CALL RTR1H1(JACTMO,NACTMO,JACTMO,NACTMO,
     &                 HMATX(1,1,IZ2,IREPABR),IPQTOQ(1,IREPM),
     &                 WORK(KBMAT),IPQTOQ(1,JSYMOP-1),
     &                 HMAT(1,1,IZ2,IREPABR),IPQTOQ(1,IREPAB),IREPXY,
     &                 WORK,KFREE,LFREE)
               END IF
            END DO
         END IF
      END DO
      CALL MEMREL('RGETH2TX2',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('RGETH2TX2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtuubl */
      SUBROUTINE GTUUBL(H2GG,NKL,H2VP,IREPVP,NSTR)
C***********************************************************************
C
C     Transform H2GG from packed format to H2UU
C
C     Input:
C        IREPVP: symmetry of H2FULL
C        H2GG: symmetry packed:
C                   IREPVP = gerade:
C                       index3(gerade) * index4(gerade) + 
C                       index3(ungerade) * index4(ungerade)
C                   IREPVP = ungerade:
C                       index3(gerade) * index4(ungerade) + 
C                       index3(ungerade) * index4(gerade)
C        NKL   : leading dimension of H2GG
C 
C     Output:
C        H2VP   : integrals in (NASHT,NASHT,*) format.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "dgroup.h"     
#include "dcborb.h" 
#include "dcbopt.h" 
#include "dcbidx.h"
C
      DIMENSION H2GG(NKL,NZ,3)
      DIMENSION H2VP(NASHT,NASHT,NZ,3)
      DIMENSION NSTR(2,0:2,4)
C
      CALL DZERO(H2VP, NASHT * NASHT * NZ * 3)
C
      IF (JTRLVL .EQ. 0) THEN
         WRITE(LUPRI,'(A,/,A,I3)')
     &        '*** ERROR in GTUUBL ***',
     &        'Illegal transformation level: ',JTRLVL
         CALL QUIT('*** ERROR in GTUUBL ***')
      ELSE IF (JTRLVL .EQ. 2 ) THEN
         WRITE(LUPRI,'(A,/,A,I3)')
     &        '*** ERROR in GTUUBL ***',
     &        'Not implemented; transformation level: ',JTRLVL
         CALL QUIT('*** ERROR in GTUUBL ***')
      ELSE 
C
C        (gg|gu) or (gg|gg) integrals transformation.
C
         IXY = 0
         IYY = 0
         DO IREPY = 1, NFSYM
            IREPX = MOD (IREPY + IREPVP, 2) + 1
            DO IY = 1, NSTR(IREPY,0,4)
               IYY = IYY + 1
               DO IX = 1, NSTR(IREPX,0,3)
                  IF (IREPX .EQ. 1) THEN
                     IXX = IX
                  ELSE
                     IXX = IX + NSTR(1,0,3)
                  END IF
                  IXY = IXY + 1
C
C                 Get active indices
C
                  IF (JTRLVL .EQ. 1) THEN
                     IF (.NOT. OPT_NOPFQ) THEN
                        IXU = IDXG2U(IXX)
                     ELSE
                        IXU = IDXG2U(IDXE2G(IXX))
                     END IF
                     IYU = IYY
                  ELSE
                     IF (.NOT. OPT_NOPFQ) THEN
                        IXU = IDXG2U(IXX)
                        IYU = IDXG2U(IYY)
                     ELSE
                        IXU = IDXG2U(IDXE2G(IXX))
                        IYU = IDXG2U(IDXE2G(IYY))
                     END IF
                  END IF
C
C                 (v,p) elements
C                     
                  IF (IXU .GT. 0 .AND. IYU .GT. 0) THEN
                     DO I3 = 1, 3
                        DO IZ = 1, NZ
                           H2VP(IXU,IYU,IZ,I3) = 
     &                          H2GG(IXY,IZ,I3)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck raddfqx */
      SUBROUTINE RADDFQX(FQX,H2VPX,PV,IV,IP,IREPVP,IPRINT)
C***********************************************************************
C
C
C     Add the FQX_{pq} += \sum_{vxy} (pv|x~y~) P_{qv,xy}
C     terms to the FQX matrix
C
C     For proper documentation see Ph.D. thesis of Joern Thyssen where
C     all the formulae below should be documented :-)
C
C     Input:
C        H2VP: the (vp|x~y~) integrals
C        PV  : the two-electron density matrix
C        IV, IP: (v,p) index
C        IREPVP: the symmetry of the integrals (currently unused)
C 
C     Output:
C        FQX    : the \sum_{vxy} (pv|x~y~) P_{qv,xy} term added to
C                 existing FQX.
C
C     Written by J. Thyssen - Dec 6 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "consts.h"
C
#include "dgroup.h"     
#include "dcborb.h" 
#include "dcbopt.h" 
#include "dcbidx.h" 
C
      DIMENSION H2VPX(NASHT,NASHT,NZ,3)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION FQX(NORBT,NASHT,NZ)
C
C     For each part of FQ (1, i, j, or k) we give
C     the type 'PV' or 'VP' for the integrals, and
C     'QV' or 'VQ' for the density matrix.
C
C     First index is over NZ for PV, H2VPX
C     Second index is 1..3
C     Third index is NZ for FQ.
C
      CHARACTER*2 szPV(4,3,4), szQV(4,3,4)
      CHARACTER*2 szPV1, szQV1
      DIMENSION dFAC(4,3,4)
      DIMENSION iPVINT(2,4,3,4)
      DIMENSION iQVINT(2,4,3,4)
C
      DATA szPV /
     &     'PV', 'PV', 'PV', 'PV', ! real FQ
     &     'PV', 'PV', 'PV', 'PV',
     &     'PV', 'PV', 'VP', 'VP',
     &     'PV', 'PV', 'PV', 'PV', ! i-imag FQ
     &     'PV', 'PV', 'PV', 'PV',
     &     'PV', 'PV', 'VP', 'VP',
     &     'PV', 'PV', 'PV', 'PV', ! j-imag FQ
     &     'PV', 'PV', 'PV', 'PV',
     &     'VP', 'VP', 'PV', 'PV',
     &     'PV', 'PV', 'PV', 'PV', ! k-imag FQ
     &     'PV', 'PV', 'PV', 'PV',
     &     'VP', 'VP', 'PV', 'PV' /
      DATA szQV /
     &     'QV', 'QV', 'QV', 'QV', ! real FQ
     &     'QV', 'QV', 'QV', 'QV',
     &     'QV', 'QV', 'VQ', 'VQ',
     &     'QV', 'QV', 'QV', 'QV', ! i-imag FQ
     &     'QV', 'QV', 'QV', 'QV',
     &     'QV', 'QV', 'VQ', 'VQ',
     &     'QV', 'QV', 'VQ', 'VQ', ! j-imag FQ
     &     'QV', 'QV', 'VQ', 'VQ',
     &     'QV', 'QV', 'QV', 'QV',
     &     'QV', 'QV', 'VQ', 'VQ', ! k-imag FQ
     &     'QV', 'QV', 'VQ', 'VQ',
     &     'QV', 'QV', 'QV', 'QV' /
      DATA dFAC /
     &     D1  , DM1 , DP5 , DMP5,
     &     DP5 , DMP5, D1  , DM1,
     &     DP5 , DMP5, DP5 , DMP5,
     &     D1  , D1  , DP5 , DP5,
     &     DP5 , DP5 , D1  , D1,
     &     DMP5, DMP5, DMP5, DMP5,
     &     D1  , DM1 , D1  , DM1,
     &     DMP5, DP5 , DP5 , DMP5,
     &     DMP5, DP5 , DP5 , DMP5,
     &     D1  , D1  , D1  , D1,
     &     DP5 , DP5 , DP5 , DP5,
     &     DP5 , DP5 , DP5 , DP5 /
      DATA iPVINT /
     &     1,1, 2,1, 3,1, 4,3,
     &     1,2, 2,2, 3,2, 4,2,
     &     1,3, 2,3, 3,1, 4,1,
     &     1,1, 2,1, 3,1, 4,1,
     &     1,2, 2,2, 3,2, 4,2,
     &     1,3, 2,3, 3,1, 4,1,
     &     1,1, 2,1, 3,2, 4,2,
     &     1,3, 2,3, 1,2, 2,2,
     &     3,1, 4,1, 3,1, 4,1,
     &     1,1, 2,1, 3,2, 4,2,
     &     1,3, 2,3, 1,2, 2,2,
     &     3,1, 4,1, 3,1, 4,1/
      DATA iQVINT /
     &     1,1, 2,1, 3,1, 4,3,
     &     1,2, 2,2, 3,2, 4,2,
     &     1,3, 2,3, 3,1, 4,1,
     &     2,1, 1,1, 4,1, 3,1,
     &     2,2, 1,2, 4,2, 3,2,
     &     2,3, 1,3, 4,1, 3,1,
     &     3,3, 4,3, 1,1, 2,1,
     &     3,1, 4,1, 3,1, 4,1,
     &     1,2, 2,2, 1,3, 2,3,
     &     4,3, 3,3, 2,1, 1,1,
     &     4,1, 3,1, 4,1, 3,1,
     &     2,2, 1,2, 2,3, 1,3/
C
C     So out of this we can establish the first term as:
C
C     D1 * \Twoint{pv,xy;1,1} P_{qv,xy;1,1}
C
C     and that the 11th term is
C
C     DP5 * \Twoint{vq,xy;3,1} P_{vq,xy;3,1}
C
C     etc...
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Output from RADDFQX',-1)
         WRITE(LUPRI,'(A,I3)') 'IV        = ',IV
         WRITE(LUPRI,'(A,I3)') 'IP        = ',IP
         WRITE(LUPRI,'(A,I3)') 'IREPVP    = ',IREPVP
         WRITE(LUPRI,'(/A/)') ' (RADDFQX) FQX on entry'
         CALL PRQMAT(FQX,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         DO I = 1, 3
            WRITE(6,*) '(RADDFQX) integral matrix no. ',I
            CALL PRQMAT(H2VPX(1,1,1,I),NASHT,NASHT,NASHT,NASHT,NZ,
     &           IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
C     (vp) is active-general
C
      IVU = IDXG2U(IV)
      IF (IVU .GT. 0) THEN
c        write(6,*) '(vp active general)'
         DO IQU = 1, NASHT
            NVQ = MAX(IVU,IQU) * ( MAX(IVU,IQU) - 1 ) / 2 + MIN(IQU,IVU)
            DO IZFQ = 1, NZ
               DO I3 = 1, 3
                  DO IZ = 1, NZ
C                    
#ifdef UNDEF
                     WRITE(LUPRI,*)
     &                    szPV(IZ,I3,IZFQ),szQV(IZ,I3,IZFQ),
     &                    dFAC(IZ,I3,IZFQ),
     &                    iPVINT(1,IZ,I3,IZFQ),iPVINT(2,IZ,I3,IZFQ),
     &                    iQVINT(1,IZ,I3,IZFQ),iQVINT(2,IZ,I3,IZFQ),
     &                    IVU,IQU,NVQ,IZFQ,I3,IZ
#endif
C
                     FAC = dFAC(IZ,I3,IZFQ)
                     szPV1 = szPV(IZ,I3,IZFQ)
                     szQV1 = szQV(IZ,I3,IZFQ)
                     iPVINT1 = iPVINT(1,IZ,I3,IZFQ)
                     iPVINT2 = iPVINT(2,IZ,I3,IZFQ)
                     iQVINT1 = iQVINT(1,IZ,I3,IZFQ)
                     iQVINT2 = iQVINT(2,IZ,I3,IZFQ)
                     NVQ1 = NVQ
                     CALL RADDFQX1(
     &                    szPV1,szQV1,FAC,
     &                    iPVINT1,iPVINT2,
     &                    iQVINT1,iQVINT2,
     &                    H2VPX,PV,NVQ1,IQU,IVU,
     &                    FQX(IP,IQU,IZFQ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C     (pv) is active-general
C
      IPU = IDXG2U(IP)
      IF (IPU .GT. 0 .AND. IPU .NE. IVU) THEN
c        write(6,*) '(pv active general)'
         DO IQU = 1, NASHT
            NPQ = MAX(IPU,IQU) * ( MAX(IPU,IQU) - 1 ) / 2 + MIN(IPU,IQU)
            DO IZFQ = 1, NZ
               DO I3 = 1, 3
                  DO IZ = 1, NZ
C                    
#ifdef UNDEF
                     WRITE(LUPRI,*)
     &                    szPV(IZ,I3,IZFQ),szQV(IZ,I3,IZFQ),
     &                    dFAC(IZ,I3,IZFQ),
     &                    iPVINT(1,IZ,I3,IZFQ),iPVINT(2,IZ,I3,IZFQ),
     &                    iQVINT(1,IZ,I3,IZFQ),iQVINT(2,IZ,I3,IZFQ),
     &                    IPU,IQU,-NPQ,IZFQ,I3,IZ
#endif
C
                     FAC = dFAC(IZ,I3,IZFQ)
                     szPV1 = szPV(IZ,I3,IZFQ)
                     szQV1 = szQV(IZ,I3,IZFQ)
                     iPVINT1 = iPVINT(1,IZ,I3,IZFQ)
                     iPVINT2 = iPVINT(2,IZ,I3,IZFQ)
                     iQVINT1 = iQVINT(1,IZ,I3,IZFQ)
                     iQVINT2 = iQVINT(2,IZ,I3,IZFQ)
                     NPQ1 = -NPQ
                     CALL RADDFQX1(
     &                    szPV1,szQV1,FAC,
     &                    iPVINT1,iPVINT2,
     &                    iQVINT1,iQVINT2,
     &                    H2VPX,PV,NPQ1,IQU,IPU,
     &                    FQX(IV,IQU,IZFQ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C     Output section
C
      IF (IPRINT .GE. 20) THEN
         WRITE(LUPRI,'(/A/)') ' (RADDFQX) FQX on exit'
         CALL PRQMAT(FQX,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
      RETURN
      END
C
      SUBROUTINE RADDFQX1(szPV,szQV,FAC,IZPV,I3PV,IZQV,I3QV,
     &     H2VPX,PV,NVQ,IQU,IVU,FQX)
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION H2VPX(NASHT,NASHT,NZ,3)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      CHARACTER szPV*2, szQV*2

C
      CHARACTER szPVXY*2, szQVXY*2
C
C     Calculate \sum_{pv|xy}_{i,j} P_{qv,xy;k,l}
C
C     from integrals {xy|vp}
C
      IF (NVQ .LT. 0) THEN
         IF (szPV .EQ. 'PV') THEN
            szPV = 'VP'
         ELSE
            szPV = 'PV'
         END IF
      END IF
c     write(6,*) 'before RADDFQX2',szPV,FAC,IZPV,I3PV
      CALL RADDFQX2(szPV,FAC,IZPV,I3PV,szPVXY)
c     write(6,*) 'after RADDFQX2',szPV,FAC,IZPV,I3PV,szPVXY
C
      IF (NVQ .LT. 0) THEN
         IF (szQV .EQ. 'QV') THEN
            szQV = 'VQ'
         ELSE
            szQV = 'QV'
         END IF
      END IF
c     write(6,*) 'before RADDFQX2',szQV,FAC,IZQV,I3QV
      CALL RADDFQX2(szQV,FAC,IZQV,I3QV,szQVXY)
c     write(6,*) 'after RADDFQX2',szQV,FAC,IZQV,I3QV,szQVXY
C
C     If NVQ > 0:
C
C     We have density matrices: (**|vq),
C     so if v<q then get (**|qv) instead.
C
C     If NVQ < 0:
C
C     We have density matrices: (**|qp
C
      IF ( (IVU .LT. IQU) .AND. NVQ .GT. 0) THEN
         CALL RADDFQX3(szQV,FAC,IZQV,I3QV,szQVXY)
      END IF
      IF ( (IQU .LT. IVU) .AND. NVQ .LT. 0) THEN
         CALL RADDFQX3(szQV,FAC,IZQV,I3QV,szQVXY)
      END IF
      NVQ = ABS(NVQ)
c     write(6,*) 'after RADDFQX3',szQV,FAC,IZQV,I3QV,szQVXY
c     write(6,*) 'NVQ = ',NVQ
C
C     Do the DDOT:
C
c     write(6,*) 'FQX before ',FQX
      CALL RADDFQX4(szPVXY,szQVXY,FAC,
     &     H2VPX(1,1,IZPV,I3PV),PV(1,1,NVQ,IZQV,I3QV),
     &     FQX)
c     write(6,*) 'FQX before ',FQX
C
      RETURN
      END
C
      SUBROUTINE RADDFQX4(szPVXY,szQVXY,FAC,H2VPX,PV,FQX)
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
#include "dcborb.h"
C
      DIMENSION H2VPX(NASHT,NASHT), PV(NASHT,NASHT)
      CHARACTER szPVXY*2, szQVXY*2
C
      IF (szPVXY .EQ. szQVXY) THEN
C
C        \sum_{xy} H2VPX(x,y) PV(x,y)
C        or
C        \sum_{xy} H2VPX(y,x) PV(y,x)
C
         FQX = FQX + FAC * DDOT(N2ASHX,H2VPX,1,PV,1)
C
      ELSE
C
C        \sum_{xy} H2VPX(x,y) PV(y,x)
C        or
C        \sum_{xy} H2VPX(y,x) PV(x,y)
C
         DO I = 1, NASHT
            DO J = 1, NASHT
               FQX = FQX + FAC * H2VPX(I,J) * PV(J,I)
            END DO
         END DO
      END IF
      RETURN
      END
C
C
      SUBROUTINE RADDFQX3(sz,FAC,IZ,I3,szXY)
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
      CHARACTER sz*2, szXY*2
C
      LOGICAL INTCH
C
      INTCH = .FALSE.
C
C     (Eq. A36 in my Ph.D. thesis)
C
      IF (I3 .EQ. 1) THEN
         IF (IZ .EQ. 1) THEN
            INTCH = .TRUE.
         ELSE IF (IZ .EQ. 2) THEN
            INTCH = .TRUE.
            FAC = FAC * DM1
         ELSE
            FAC = FAC * DM1
         ENDIF
      ELSE IF (I3 .EQ. .2) THEN
         IF (IZ .EQ. 1 .OR. IZ .EQ. 2) THEN
            FAC = FAC * DM1
         ELSE IF (IZ .EQ. 3) THEN
            INTCH = .TRUE.
            I3 = 3
         ELSE
            INTCH = .TRUE.
            FAC = FAC * DM1
            I3 = 3
         END IF
      ELSE
         IF (IZ .EQ. 1 .OR. IZ .EQ. 2) THEN
            FAC = FAC * DM1
         ELSE IF (IZ .EQ. 3) THEN
            INTCH = .TRUE.
            I3 = 2
         ELSE
            INTCH = .TRUE.
            FAC = FAC * DM1
            I3 = 2
         END IF
      END IF
C
      IF (INTCH) THEN
         IF (szXY .EQ. 'XY') THEN
            szXY = 'YX'
         ELSE
            szXY = 'XY'
         END IF
      END IF
C
      RETURN
      END 
C
C
      SUBROUTINE RADDFQX2(sz,FAC,IZ,I3,szXY)
#include "priunit.h"
#include "implicit.h"
C
#include "consts.h"
C
      CHARACTER SZ*2, szXY*2
C
      szXY = 'XY'
      IF (I3 .EQ. 1) THEN
C
C
         IF (IZ .EQ. 1) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (pv|xy) = (xy|pv) = (yx|vp)^*
               szXY = 'YX'
C           ELSE
C              We want (vp|xy) = (xy|vp)
C              szXY = 'XY'
            END IF
C
C
         ELSE IF (IZ .EQ. 2) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
               szXY = 'YX'
               FAC = FAC * DM1
            END IF
C
C
         ELSE IF (IZ .EQ. 3) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (pv|xY) = (xY|pv) = (Yx|vp)^* (class 3,2 and 4,2)
               szXY = 'YX'
               I3 = 3
            ELSE
C              We want (vp|xY) = (xY|vp) (class 3,3 and 4,3)
               IZ= 3
               I3= 3
            END IF
C
C
         ELSE
C
C
            IF (sz(2:2) .EQ. 'V') THEN
               szXY = 'YX'
               I3 = 4
               FAC = FAC * DM1
            ELSE
               IZ = 4
               I3 = 3
            END IF
C
C
         END IF
         RETURN
      ELSE IF (I3 .EQ. 2) THEN
C
C
         IF (IZ .EQ. 1 .OR. IZ .EQ. 2) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (pV|xY) = (xY|pV) = (yX|vP)
               szXY = 'YX'
C           ELSE
C              We want (vP|xY) = (xY|vP)
C              szXY = 'XY'
            END IF
C
C
         ELSE
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (pV|xy) = (xy|pV) = - (xy|vP)
               I3 = 1
               FAC = FAC * DM1
            ELSE
C              We want (vP|xy) = (xy|vP)
               I3 = 1
            END IF
C
C
         END IF
C
C
      ELSE
C
C
         IF (IZ .EQ. 1) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (Pv|xY) = (xY|Pv) = (Yx|vP)^*
               szXY = 'YX'
            ELSE
C              We want (Vp|xY) = (xY|Vp) = - (xY|Pv) = - (Yx|vP)^*
               szXY = 'YX'
               FAC = FAC * DM1
            END IF
C
C
         ELSE IF (IZ .EQ. 2) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
               szXY = 'YX'
               FAC = FAC * DM1
            ELSE
               szXY = 'YX'
            END IF
C
C
         ELSE IF (IZ .EQ. 3) THEN
C
C
            IF (sz(2:2) .EQ. 'V') THEN
C              We want (Pv|xy) = (xy|Pv) = (yx|vP)^*
               szXY = 'YX'
               I3 = 1
            ELSE
C              We want (Vp|xy) = (xy|Vp) = - (xy|Pv) = - (yx|vP)^*
               szXY = 'YX'
               I3 = 1
               FAC = FAC * DM1
            END IF
C
C
         ELSE
C
C
            IF (sz(2:2) .EQ. 'V') THEN
               szXY = 'YX'
               I3 = 1
               FAC = FAC * DM1
            ELSE
               szXY = 'YX'
               I3 = 1
            END IF
C
C
         END IF
      END IF
      RETURN
      END
C
C

#ifdef UNDEF
C     Loop over q
C
C     add to FQX_{pq}:
C
      IVU = IDXG2U(IV)
      IF ( IVU .GT. 0) THEN
C        
C        Loop over Q <= V:
C
         DO IQU = 1, IVU
C
            NVQ = IVU * ( IVU - 1 ) / 2 + IQU
C
C           ************************
C           *** Real part of FQX ***
C           ************************
C
C           Real terms:
C           -----------
C
C           Below the following terms:
C              +     (pv|xy)_{1,1} P_{qv,xy;1,1}
C              + 1/2 (pv|xy)_{1,2} P_{qv,xy;1,2}
C              + 1/2 (pv|xy)_{1,3} P_{qv,xy;1,3}
C
C           
C           FQX_{pq} +=
C             \sum_{xy} \Twoint{vp|yx}_{1,1} P_{vq,yx;1,1}
C           + 1/2 \sum_{xy}
C                \Twoint{v\bar{p}|\bar{x}y}_{1,2} P_{v\bar{q},\bar{x}y;1,2}
C           + 1/2 \sum_{xy}
C                \Twoint{v\bar{p}|x\bar{y}}_{1,3} P_{v\bar{q},x\bar{y};1,3}
C           
            FQX(IP,IQU,1) = FQX(IP,IQU,1) 
     &           +       DDOT(N2ASHX,H2VPX(1,1,1,1),1,
     &                               PV(1,1,NVQ,1,1),1)
     &           + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &                               PV(1,1,NVQ,1,2),1) 
     &           + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,
     &                               PV(1,1,NVQ,1,3),1)
C
            IF (NZ .GE. 2) THEN
C
C
C              ************************
C              *** Real part of FQX ***
C              ************************
C
C
C              Complex terms:
C              --------------
C
C              Below the following terms:
C                 - 1/2 (pv|xy)_{2,1} P_{qv,xy;2,1}
C                 - 1/2 (pv|xy)_{2,2} P_{qv,xy;2,2}
C                 - 1/2 (pv|xy)_{2,3} P_{qv,xy;2,3}
C
C       - \sum_{xy} [ - \Twoint{yx|vp}_{2,1} ] [ - P_{yx,vq;2,1} ]
C       - 1/2 \sum_{xy} [ - \Twoint{yx|vp}_{2,2} ] [ - P_{yx,vq;2,2} ]
C       - 1/2 \sum_{xy} [ - \Twoint{yx|vp}_{2,3} ] [ - P_{xy,vq;2,3} ]
C
               FQX(IP,IQU,1) = FQX(IP,IQU,1)
     &              -       DDOT(N2ASHX,H2VPX(1,1,2,1),1,
     &                                  PV(1,1,NVQ,2,1),1)
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,2),1,
     &                                  PV(1,1,NVQ,2,2),1)
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,3),1,
     &                                  PV(1,1,NVQ,2,3),1)
C
C
C              **************************
C              *** i-imag part of FQX ***
C              **************************
C
C              Below are the following terms:
C                 +     (pv|xy)_{1,1} P_{qv,xy;2,1}
C                 +     (pv|xy)_{2,1} P_{qv,xy;1,1}
C                 + 1/2 (pv|xy)_{1,2} P_{qv,xy;2,2}
C                 + 1/2 (pv|xy)_{2,2} P_{qv,xy;1,2}
C                 - 1/2 (pv|xy)_{1,3} P_{qv,xy;2,3}
C                 - 1/2 (pv|xy)_{2,3} P_{qv,xy;1,3}
C
C       +     \sum{xy} [ + \Twoint{yx|vp}_{1,1} ] [ - P_{yx,vq;2,1} ]
C       +     \sum{xy} [ - \Twoint{yx|vp}_{2,1} ] [ + P_{yx,vq;1,1} ]
C       + 1/2 \sum{xy} [ + \Twoint{yx|vp}_{1,2} ] [ - P_{yx,vq;2,2} ]
C       + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{2,2} ] [ + P_{yx,vq;1,2} ]
C       - 1/2 \sum{xy} [ + \Twoint{yx|vp}_{1,3} ] [ - P_{yx,vq;2,3} ]
C       - 1/2 \sum{xy} [ - \Twoint{yx|vp}_{2,3} ] [ + P_{yx,vq;1,3} ]
C
               FQX(IP,IQU,2) = FQX(IP,IQU,2)
     &              -       DDOT(N2ASHX,H2VPX(1,1,1,1),1,
     &                                  PV(1,1,NVQ,2,1),1) 
     &              -       DDOT(N2ASHX,H2VPX(1,1,2,1),1,
     &                                  PV(1,1,NVQ,1,1),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &                                  PV(1,1,NVQ,2,2),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,2),1,
     &                                  PV(1,1,NVQ,1,3),1) 
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &                                  PV(1,1,NVQ,2,3),1) 
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,2,3),1,
     &                                  PV(1,1,NVQ,1,3),1) 
C
C
            END IF
C
            IF (NZ .GE. 4) THEN
C
C
C              ************************
C              *** Real part of FQX ***
C              ************************
C
C
C              Below the following terms:
C                 + 1/2 (pv|xy)_{3,1} P_{qv,xy;3,1}
C                 - 1/2 (pv|xy)_{4,1} P_{qv,xy;4,1}
C                 +     (pv|xy)_{3,2} P_{qv,xy;3,2}
C                 -     (pv|xy)_{4,2} P_{qv,xy;4,2}
C                 + 1/2 (vp|xy)_{3,1} P_{vq,xy;3,1}
C                 - 1/2 (vp|xy)_{4,1} P_{vq,xy;4,1}
C              
C
C       + 1/2 \sum_{xy} [ + \Twoint{yx|vp}_{3,3} ] [ + P_{yx,vq;3,3} ]
C       - 1/2 \sum_{xy} [ - \Twoint{yx|vp}_{3,4} ] [ - P_{yx,vq;3,4} ]
C       + \sum_{xy} [ - \Twoint{xy|vp}_{3,1} ] [ - P_{xy,vq;3,1} ]
C       - \sum_{xy} [ - \Twoint{xy|vp}_{4,1} ] [ - P_{xy,vq;4,1} ]
C       + 1/2 \sum_{xy} [ + \Twoint{xy|vp}_{3,2} ] [ + P_{xy,vq;3,2} ]
C       - 1/2 \sum_{xy} [ - \Twoint{xy|vp}_{3,4} ] [ - P_{xy,vq;3,4} ]
C
               FQX(IP,IQU,1) = FQX(IP,IQU,1)
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,3,3),1,
     &                                  PV(1,1,NVQ,3,3),1)
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,4),1,
     &                                  PV(1,1,NVQ,3,4),1)
     &              +       DDOT(N2ASHX,H2VPX(1,1,3,1),1,
     &                                  PV(1,1,NVQ,3,1),1)
     &              -       DDOT(N2ASHX,H2VPX(1,1,4,1),1,
     &                                  PV(1,1,NVQ,4,1),1)
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,3,2),1,
     &                                  PV(1,1,NVQ,3,2),1)
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,2),1,
     &                                  PV(1,1,NVQ,3,2),1)
C
C
C              **************************
C              *** i-imag part of FQX ***
C              **************************
C
C
C              Below are the following terms:
C                 + 1/2 (pv|xy)_{3,1} P_{qv,xy;4,1}
C                 + 1/2 (pv|xy)_{4,1} P_{qv,xy;3,1}
C                 +     (pv|xy)_{3,2} P_{qv,xy;4,2}
C                 +     (pv|xy)_{4,2} P_{qv,xy;3,2}
C                 - 1/2 (vp|xy)_{3,1} P_{vq,xy;4,1}
C                 - 1/2 (vp|xy)_{4,1} P_{vq,xy;3,1}
C
C       + 1/2 \sum{xy} [ + \Twoint{yx|vp}_{3,3} ] [ - P_{yx,vq;4,3} ]
C       + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{4,3} ] [ + P_{yx,vq;3,3} ]
C       +     \sum{xy} [ - \Twoint{xy|vp}_{3,1} ] [ - P_{xy,vq;4,1} ]
C       +     \sum{xy} [ - \Twoint{xy|vp}_{4,1} ] [ - P_{xy,vq;3,1} ]
C       - 1/2 \sum{xy} [ + \Twoint{xy|vp}_{3,2} ] [ + P_{xy,vq;4,2} ]
C       - 1/2 \sum{xy} [ + \Twoint{xy|vp}_{4,2} ] [ + P_{xy,vq;3,2} ]
C
               FQX(IP,IQU,2) = FQX(IP,IQU,2)
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,3),1,
     &                                  PV(1,1,NVQ,4,3),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,4,3),1,
     &                                  PV(1,1,NVQ,3,3),1) 
     &              +       DDOT(N2ASHX,H2VPX(1,1,3,1),1,
     &                                  PV(1,1,NVQ,4,1),1) 
     &              +       DDOT(N2ASHX,H2VPX(1,1,4,1),1,
     &                                  PV(1,1,NVQ,3,1),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,2),1,
     &                                  PV(1,1,NVQ,4,2),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,4,2),1,
     &                                  PV(1,1,NVQ,3,2),1) 
C
C
C              **************************
C              *** j-imag part of FQX ***
C              **************************
C
C
C              Below are the following terms:
C                 +     (pv|xy)_{1,1} P_{qv,xy;3,3}
C                 -     (pv|xy)_{2,1} P_{qv,xy;4,3}
C                 +     (pv|xy)_{3,2} P_{vq,xy;1,1}
C                 -     (pv|xy)_{4,2} P_{vq,xy;2,1}
C
C                 - 1/2 (pv|xy)_{1,3} P_{qv,xy;3,1}
C                 + 1/2 (pv|xy)_{2,3} P_{qv,xy;4,1}
C                 + 1/2 (pv|xy)_{1,2} P_{vq,xy;3,1}
C                 - 1/2 (pv|xy)_{2,2} P_{vq,xy;4,1}
C
C                 - 1/2 (vp|xy)_{3,1} P_{qv,xy;1,2}
C                 + 1/2 (vp|xy)_{4,1} P_{qv,xy;2,2}
C                 + 1/2 (pv|xy)_{3,1} P_{qv,xy;1,3}
C                 - 1/2 (pv|xy)_{4,1} P_{qv,xy;2,3}
C
C        +     \sum{xy} [ + \Twoint{yx|vp}_{1,1} ] [ - P_{yx,vq;3,1} ]
C        -     \sum{xy} [ - \Twoint{yx|vp}_{2,1} ] [ - P_{yx,vq;4,1} ]
C        +     \sum{xy} [ - \Twoint{xy|vp}_{3,1} ] [ + P_{xy,vq;1,1} ]
C             \sum{xy} [ - \Twoint{xy|vp}_{4,1} ] [ + P_{xy,vq;2,1} ]
C
C        - 1/2 \sum{xy} [ + \Twoint{yx|vp}_{1,3} ] [ + P_{yx,vq;3,3} ]
C        + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{2,3} ] [ - P_{yx,vq;4,3} ]
C        + 1/2 \sum{xy} [ + \Twoint{xy|vp}_{1,2} ] [ - P_{xy,vq;3,2} ]
C        - 1/2 \sum{xy} [ - \Twoint{xy|vp}_{2,2} ] [ - P_{xy,vq;4,2} ]
C
C        - 1/2 \sum{xy} [ - \Twoint{yx|vp}_{3,2} ] [ + P_{yx,vq;1,2} ]
C        + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{4,2} ] [ - P_{yx,vq;2,2} ]
C        + 1/2 \sum{xy} [ + \Twoint{xy|vp}_{3,3} ] [ - P_{xy,vq;1,3} ]
C        - 1/2 \sum{xy} [ - \Twoint{xy|vp}_{4,3} ] [ + P_{xy,vq;2,3} ]
C
               FQX(IP,IQU,3) = FQX(IP,IQU,3)
     &            -       DDOT(N2ASHX,H2VPX(1,1,1,1),1,
     &                                PV(1,1,NVQ,3,1),1) 
     &            -       DDOT(N2ASHX,H2VPX(1,1,2,1),1,
     &                                PV(1,1,NVQ,4,1),1) 
     &            -       DDOT(N2ASHX,H2VPX(1,1,3,1),1,
     &                                PV(1,1,NVQ,1,1),1) 
     &            +       DDOT(N2ASHX,H2VPX(1,1,4,1),1,
     &                                PV(1,1,NVQ,2,1),1) 
     &            - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,
     &                                PV(1,1,NVQ,3,3),1) 
     &            + DP5 * DDOT(N2ASHX,H2VPX(1,1,2,3),1,
     &                                PV(1,1,NVQ,4,3),1) 
     &            - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &                                PV(1,1,NVQ,3,2),1) 
     &            - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,2),1, !FIXME: sign???
     &                                PV(1,1,NVQ,4,2),1) 
     &            + DP5 * DDOT(N2ASHX,H2VPX(1,1,3,2),1,
     &                                PV(1,1,NVQ,1,2),1) 
     &            + DP5 * DDOT(N2ASHX,H2VPX(1,1,4,2),1,
     &                                PV(1,1,NVQ,2,2),1) 
     &            - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,3),1,
     &                                PV(1,1,NVQ,1,3),1) 
     &            + DP5 * DDOT(N2ASHX,H2VPX(1,1,4,3),1, 
     &                                PV(1,1,NVQ,2,3),1) 
C
C
C              **************************
C              *** j-imag part of FQX ***
C              **************************
C
C
C              Below are the following terms:
C                 +     (pv|xy)_{1,1} P_{qv,xy;4,3}
C                 +     (pv|xy)_{2,1} P_{qv,xy;3,3}
C                 +     (pv|xy)_{3,2} P_{vq,xy;2,1}
C                 +     (pv|xy)_{4,2} P_{vq,xy;1,1}
C
C                 + 1/2 (pv|xy)_{1,3} P_{qv,xy;4,1}
C                 + 1/2 (pv|xy)_{2,3} P_{qv,xy;3,1}
C                 + 1/2 (pv|xy)_{1,2} P_{vq,xy;4,1}
C                 + 1/2 (pv|xy)_{2,2} P_{vq,xy;3,1}
C
C                 + 1/2 (vp|xy)_{3,1} P_{qv,xy;2,2}
C                 + 1/2 (vp|xy)_{4,1} P_{qv,xy;1,2}
C                 + 1/2 (pv|xy)_{3,1} P_{qv,xy;2,3}
C                 + 1/2 (pv|xy)_{4,1} P_{qv,xy;1,3}
C
C        +     \sum{xy} [ + \Twoint{yx|vp}_{1,1} ] [ - P_{yx,vq;4,1} ]
C        -     \sum{xy} [ - \Twoint{yx|vp}_{2,1} ] [ - P_{yx,vq;3,1} ]
C        +     \sum{xy} [ - \Twoint{xy|vp}_{3,1} ] [ - P_{xy,vq;2,1} ]
C        -     \sum{xy} [ - \Twoint{xy|vp}_{4,1} ] [ + P_{xy,vq;1,1} ]
C
C        - 1/2 \sum{xy} [ + \Twoint{yx|vp}_{1,3} ] [ - P_{yx,vq;4,3} ]
C        + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{2,3} ] [ + P_{yx,vq;3,3} ]
C        + 1/2 \sum{xy} [ + \Twoint{xy|vp}_{1,2} ] [ - P_{xy,vq;4,2} ]
C        - 1/2 \sum{xy} [ - \Twoint{xy|vp}_{2,2} ] [ - P_{xy,vq;3,2} ]
C
C        - 1/2 \sum{xy} [ - \Twoint{yx|vp}_{3,2} ] [ - P_{yx,vq;2,2} ]
C        + 1/2 \sum{xy} [ - \Twoint{yx|vp}_{4,2} ] [ + P_{yx,vq;1,2} ]
C        + 1/2 \sum{xy} [ + \Twoint{xy|vp}_{3,3} ] [ + P_{xy,vq;2,3} ]
C        - 1/2 \sum{xy} [ - \Twoint{xy|vp}_{4,3} ] [ - P_{xy,vq;1,3} ]
C
               FQX(IP,IQU,4) = FQX(IP,IQU,4)
     &              -       DDOT(N2ASHX,H2VPX(1,1,1,1),1,
     &                                  PV(1,1,NVQ,4,1),1) 
     &              -       DDOT(N2ASHX,H2VPX(1,1,2,1),1,
     &                                  PV(1,1,NVQ,3,1),1) 
     &              +       DDOT(N2ASHX,H2VPX(1,1,3,1),1,
     &                                  PV(1,1,NVQ,2,1),1) 
     &              +       DDOT(N2ASHX,H2VPX(1,1,4,1),1,
     &                                  PV(1,1,NVQ,1,1),1) 
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,
     &                                  PV(1,1,NVQ,4,3),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,3),1,
     &                                  PV(1,1,NVQ,3,3),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &                                  PV(1,1,NVQ,4,2),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,2,2),1,   
     &                                  PV(1,1,NVQ,3,2),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,3,2),1,
     &                                  PV(1,1,NVQ,2,2),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,4,2),1,
     &                                  PV(1,1,NVQ,1,2),1) 
     &              + DP5 * DDOT(N2ASHX,H2VPX(1,1,3,3),1,
     &                                  PV(1,1,NVQ,2,3),1) 
     &              - DP5 * DDOT(N2ASHX,H2VPX(1,1,4,3),1, 
     &                                  PV(1,1,NVQ,1,3),1) 
C
C
C
            END IF
         END DO
         DO IQU = IVU+1, NASHT
C
C           Loop over Q > V:
C           
C           FQX_{pq} +=
C             \sum_{xy} \Twoint{vp|yx}_{1,1} P_{qv,xy;1,1}
C       + 1/2 \sum_{xy}
C         \Twoint{v\bar{p}|\bar{x}y}_{1,2} [-P_{v\bar{q},\bar{x}y;1,2}]
C       + 1/2 \sum_{xy}
C         \Twoint{v\bar{p}|x\bar{y}}_{1,3} [-P_{v\bar{q},x\bar{y};1,3}]
C           
            NVQ = IQU * ( IQU - 1 ) / 2 + IVU
            DTERM = D0
            DO IX = 1, NASHT
               DO IY = 1, NASHT
                  DTERM = DTERM +
     &                 H2VPX(IX,IY,1,1) * PV(IY,IX,NVQ,1,1)
               END DO
            END DO
            FQX(IP,IQU,1) = FQX(IP,IQU,1) 
     &           + DTERM
     &           - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,
     &           PV(1,1,NVQ,1,2),1) 
     &           - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,
     &           PV(1,1,NVQ,1,3),1)
         END DO
      END IF
C
C     If (vp) is general-active we need to consider the
C     (pv) distribution instead.
C      
      IPU = IDXG2U(IP)
      IF ( IPU .GT. 0 .AND. IPU .NE. IVU ) THEN
         DO IQU = 1, IPU-1
C           
C           FQX_{vq} += \sum_{xy} \Twoint{vp|xy}_{1,1} P_{xy,qp;1,1}
C           = \sum_{xy} \Twoint{vp|xy}_{1,1} P_{yx,pq;1,1}
C           
            NQP = IPU * ( IPU - 1 ) / 2 + IQU
            DTERM = D0
            DO IX = 1, NASHT
               DO IY = 1, NASHT
                  DTERM = DTERM +
     &                 H2VPX(IX,IY,1,1) * PV(IY,IX,NQP,1,1)
               END DO
            END DO
            FQX(IV,IQU,1) = FQX(IV,IQU,1) 
     &           + DTERM
     &           - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,PV(1,1,NQP,1,2),1)
     &           - DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,PV(1,1,NQP,1,3),1)
         END DO
         DO IQU = IPU, NASHT
C           
C           FQX_{vq} += \sum_{xy} \Twoint{vp|xy}_{1,1} P_{xy,qp;1,1}
C           
            NQP = IQU * ( IQU - 1 ) / 2 + IPU
            FQX(IV,IQU,1) = FQX(IV,IQU,1) 
     &           + DDOT(N2ASHX,H2VPX(1,1,1,1),1,PV(1,1,NQP,1,1),1)
     &           + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,2),1,PV(1,1,NQP,1,2),1)
     &           + DP5 * DDOT(N2ASHX,H2VPX(1,1,1,3),1,PV(1,1,NQP,1,3),1)
         END DO
      END IF
      RETURN
      END 

#endif
