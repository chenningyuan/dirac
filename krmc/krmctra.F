C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE RTRACTL(ITRLVL,WORK,LWORK)
C*****************************************************************************
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
C     ITRLVL =  2: (pq|xy) for Newton-Raphson and NEO CI & MC
C
C
C     Written by J. Thyssen - Jun 27 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcborb.h>
C
      CHARACTER SZTRA(4,2)*80
      DIMENSION WORK(*)
C
#include <memint.h>
C
C
      CALL QENTER('RTRACTL')
C
      IF ( IPROPT .GT. 5 ) THEN
         CALL HEADER('Output from RTRACTL',-1)
         WRITE(LUPRI,'(A,I2)')
     $        '* Transformation level (ITRLVL) = ',ITRLVL
         WRITE(LUPRI,'(3A,I1,A)')
     $        '* Integrals in transformation = ',OPT_INTTYP,
     $        ' (INTFLG=',IOPT_INTFLG,')'
      END IF
C
      IF ( ITRLVL .LT. -1 .OR. ITRLVL .GT. 4 .OR. ITRLVL .EQ. 2 ) THEN
         WRITE(LUPRI,'(2(/,A),I2)')
     $        '*** ERROR in RTRACTL ***',
     $        'Illegal transformation level (ITRLVL) = ',ITRLVL
         CALL QUIT('*** ERROR in RTRACTL ***')
      END IF
C
C     Initialize all indeces with active orbitals.
C
      IF (ITRLVL .EQ. 1 .OR. ITRLVL .EQ. 0) THEN
         DO I = 1,4
            DO J = 1, NFSYM
               IF ( NASH (J) .GT. 0 ) THEN
                  WRITE(SZTRA(I,J),'(A,I4)')
     $                 '1..', NASH(J)
               ELSE
                  SZTRA(I,J) = ' '
               END IF
            END DO
         END DO
         IF (ITRLVL .EQ. 1 ) THEN
C           (uu|gu) integrals
            DO J = 1, NFSYM
               IF ( NESH (J) .GT. 0 ) THEN
                  IF ( NPSH (J) .GT. 0 ) THEN
                     WRITE(SZTRA(3,J),'(I4,A,I4)')
     $                 -NPSH(J), '..-1, 1..', NESH(J)
                  ELSE
                     WRITE(SZTRA(3,J),'(A,I4)')
     $                    '1..', NESH(J)
                  END IF
               ELSE
                  IF ( NPSH ( J ) .GT. 0 ) THEN
                     WRITE(SZTRA(3,J),'(I4,A)')
     $                    -NPSH(J), '..-1'
                  ELSE
                     SZTRA(3,J) = ' '
                  END IF
               END IF
            END DO
         END IF
      ELSE IF (ITRLVL .EQ. 2 ) THEN
C        cannot be implemented with the current structure...
      ELSE 
C        initialize all indeces with general orbitals
         DO I = 1, 4
            DO J = 1, NFSYM
               IF ( NESH (J) .GT. 0 ) THEN
                  IF ( NPSH (J) .GT. 0 ) THEN
                     WRITE(SZTRA(I,J),'(I4,A,I4)')
     $                    -NPSH(J), '..-1, 1..', NESH(J)
                  ELSE
                     WRITE(SZTRA(I,J),'(A,I4)')
     $                    '1..', NESH(J)
                  END IF
               ELSE
                  IF ( NPSH ( J ) .GT. 0 ) THEN
                     WRITE(SZTRA(I,J),'(I4,A)')
     $                    -NPSH(J), '..-1'
                  ELSE
                     SZTRA(I,J) = ' '
                  END IF
               END IF
            END DO
         END DO
         IF (ITRLVL .EQ. 3) THEN
C           (gu|gg) integrals
            DO J = 1, NFSYM
               IF ( NASH (J) .GT. 0 ) THEN
                  WRITE(SZTRA(2,J),'(A,I4)')
     $                 '1..', NASH(J)
               ELSE
                  SZTRA(2,J) = ' '
               END IF
            END DO
         END IF
      END IF
C
C     Print orbital strings
C
      IF ( IPROPT .GE. 5 ) THEN
         DO I = 1, 4
            DO J = 1,NFSYM
               WRITE(LUPRI,'(A,I1,A,I1,2A)')
     $              '* Orbital string ',I,
     $              ' for irrep ',J,' = ',SZTRA(I,J)
            END DO
         END DO
      END IF
      write(6,*) 'nopqfq ',opt_nopfq
C
C     Transfer information to moltra common blocks
C
      CALL SETUP_DCBTRA(SZTRA,IOPT_INTFLG,OPT_NOPFQ,0)
C
C     Call 4-index transformation driver
C
C
      CALL PAMTRA(WORK(KFREE),LFREE)
C
C
      CALL QEXIT('RTRACTL')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rtractl */
      SUBROUTINE SETUP_DCBTRA(SZ,INTFLG,LNOPAIR,IPRINT)
C*****************************************************************************
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
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( DM1 = -1.0D00 )
C
      CHARACTER SZ(4,2)*(*)
      LOGICAL LNOPAIR
C
#include <dcbtra.h>
C     
C     /DCBTC/
C
      DO I = 1, 4
         DO J = 1, 2
            TRA4_INDSTR(I,J) = SZ(I,J)
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
      IPRTRA = IPRINT
      ISTRAT = 4
C
C     /DCBTRL/
C
      TRA_ANTIS = .FALSE.
      TRA_MSOUT  = .FALSE.
      NO2IND = .TRUE.
      NO4IND = .FALSE.
      NOPAIR = LNOPAIR
      PRPTRA = .FALSE.
      RCORBS = .FALSE.
      FKRMC = .TRUE.
      NOMDCINT = .TRUE.
C      
C     /DCBTRR/
C
      THROUT = 1.0D-16
      SCRTRA = DM1
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth2 */
      SUBROUTINE RGETH2(FQ,H2AC,WORK,LWORK)
C*****************************************************************************
C
C     Calculate FQ and get H2AC (two-electron integrals
C     with active indeces).
C
C     Input:
C
C     Output:
C        FQ   -   aux. FQ fock matrix
C        H2AC -   two-electron integrals with active indeces.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dcbgen.h>
#include <dgroup.h>
C
      DIMENSION FQ(*)
      DIMENSION H2AC(*)
      DIMENSION WORK(*)
C
C
#include <memint.h>
C
      CALL QENTER('RGETH2')
C
      IF ( OPT_CIPROGRAM .EQ. 'GOSCI' ) THEN
C
C        Get density matrix
C        ------------------
C
C        Read density matrix.
C         
         KLEN = NASHT * NASHT * NNASHX * NZ * 3
         CALL MEMGET('REAL',KPV,KLEN,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KPVM,(2*NASHT)**4 * MIN0(NZ,2),
     &        WORK,KFREE,LFREE)
         CALL REAKRMC(LUKRM3,'PVMO    ',
     &        WORK(KPVM),(2*NASHT)**4 * MIN0(NZ,2))
C
C        Transform from Molfdir basis to Dirac (NZ,3) basis.
C        
         CALL M2DNZ3(WORK(KPVM),WORK(KPV),IPROPT)
C
         CALL MEMREL('RGETH2.PV',WORK,KWORK,KPVM,KFREE,LFREE)
C
C
C        Two electron integrals
C        ----------------------
C
C        Get FQ and active two-electron integrals.
C
         CALL RGETH21(.TRUE.,.TRUE.,WORK(KPV),FQ,H2AC,
     &        WORK(KFREE),LFREE,IPROPT)
C
C        Write FQ and H2AC to disk.
C
         CALL WRTKRMC(LUKRM3,'FQMO    ',FQ,NASHT * NORBT * NZ)
         CALL WRTKRMC(LUKRM3,'H2AC    ',H2AC,NASHT*NASHT*NNASHX*NZ*3)
         CALL MEMREL('RGETH2.after FQ',WORK,KWORK,KWORK,KFREE,LFREE)
      ELSE
         WRITE(LUPRI,'(/,A,/,2A)')
     $        '*** ERROR in RGETH2 ***',
     $        'FQ and H2AC not implemented for CI program ',
     $        OPT_CIPROGRAM
         CALL QUIT('*** ERROR in RGETH2 ***')
      END IF
      CALL QEXIT('RGETH2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgeth21 */
      SUBROUTINE RGETH21(DOFQ,DOH2AC,PV,FQ,H2AC,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C
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
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
C
#include <dcbopt.h>
#include <dcborb.h>
#include <dgroup.h>
#include <dcbidx.h>
C
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION H2AC(NASHT,NASHT,NNASHX,NZ,3)
      DIMENSION WORK(*)
      LOGICAL DOFQ, DOH2AC
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
C     strings
      CHARACTER CDATE*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
C
#include <memint.h>
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
      IF ( IPRINT .GE. 3 ) WRITE(LUPRI,'(1X,4A)')
     &     '- time and date read from MDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indeces
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
C     Boson irrep for IJ indeces
      CALL READI(LUINT,NFPCK12T,WORK(KDXB12))
C     Boson irrep for KL indeces
      CALL READI(LUINT,NFPCK34T,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C     --------------------
C
      LUINT = 0
      DO IREPIJ = 1, NFSYM
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
      END DO
C
C     Allocate memory for integrals (H2XY)
C     
      NKL = MAX0( NFPCK34(1), NFPCK34(2) )
      CALL MEMGET('REAL',KHMAT,NKL*NCLASS,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XY,NKL*NCLASS,WORK,KFREE,LFREE)
C     (allocate KH2XY with NCLASS so we can use it as scratch)
C     Use KH2XY as two different matrices: H2XY and H2YX
      KH2YX = KH2XY + 3 * NZ *
     &     MAX0( NORB(1) * NASH(1),
     &     NORB(2) * NASH(2) )
C
C     Read integrals (x>y)
C     --------------------
C
C     zero arrays
C
      IF (DOFQ) CALL DZERO(FQ,NASHT*NORBT*NZ)
      IF (DOH2AC) CALL DZERO(H2AC,NASHT*NASHT*NNASHX*NZ*3)
C
C     Loop over distributions
C
      NXY = 0
      DO IX = 1, NASHT
         DO IY = 1, IX
            NXY = NXY + 1
C
C           Find record number and irep
C
            CALL FNDRECNO(IX,IY,WORK(KDXKR12),NFPCK12,NFPCK12T,
     &           IREPIJ,IREC)
            write(6,*) 'ix iy',IX,IY,IREPIJ,IREC
C
C           Read integrals
C           
            NKL = NFPCK34(IREPIJ)
C
            READ (LUINT + IREPIJ, REC=IREC)
     &           (WORK(KH2XY+I-1),I=1,NKL*NCLASS) 
#ifdef UNDEF
            ii = 0
            do i = 1, 3 * nz
               write(6,*) 'class ',i,' before transformation'
               call output(work(kh2xy+ii),1,norbt,1,norbt,
     &              norbt,norbt,1,lupri)
               ii = ii + nkl
            end do
#endif
C
C
C           Transform to Dirac (NZ,3) format.
C
            CALL DZERO(WORK(KHMAT),NKL * 3 * NZ )
            CALL CL2NZ3(WORK(KH2XY),WORK(KHMAT),NKL,IREPIJ)
            ii = 0
            do i = 1, 3 * nz
               write(6,*) 'class ',i,' after transformation'
               call output(work(khmat+ii),1,norbt,1,norbt,
     &              norbt,norbt,1,lupri)
               ii = ii + nkl
            end do
C
C           - take out NORBT * NASHT  XY- and YX-blocks.
C           - make H2AC.
C
            CALL TKNONABL(WORK(KHMAT),NKL,
     &           WORK(KH2XY),WORK(KH2YX),H2AC,NXY,
     &           IREPIJ,NSTR,IX .EQ. IY, DOFQ, DOH2AC)
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
C
C           Do the gazillion matrix multiplies to construct
C           the auxililary Fock matrix FQ.
C           
            IF (DOFQ) THEN
               CALL RADDFQ(FQ,WORK(KH2XY),WORK(KH2YX),PV,
     &              NXY,IX,IY,IREPIJ,IPRINT)
            END IF
C
         END DO
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,3(D14.7,1X),3(/,16X,3(D14.7,1X)))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('H2AC matrix in Dirac (NZ,3) format',-1)
         IJ = 0
         DO I = 1,NASHT
            DO J = 1,I
               IJ = IJ + 1
               DO K = 1,NASHT
                  DO L = 1,NASHT
                     WRITE(LUPRI,9000)
     &                    K,L,I,J,
     &                    ((H2AC(K,L,IJ,IZ,ICL),
     &                    IZ = 1,NZ),ICL=1,3)
                  END DO
               END DO
            END DO
         END DO
C
         CALL HEADER('FQ matrix',-1)
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL QEXIT('RGETH21')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck raddfq */
      SUBROUTINE RADDFQ(FQ,H2XY,H2YX,PV,NXY,IX,IY,IREPIJ,IPRINT)
C*****************************************************************************
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
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER ( D1 = 1.0D00 )
C
#include <dgroup.h>     
#include <dcborb.h> 
#include <dcbopt.h> 
C
      DIMENSION H2XY(NORBT,NASHT,NZ,3)
      DIMENSION H2YX(NORBT,NASHT,NZ,3)
      DIMENSION FQ(NORBT,NASHT,NZ)
      DIMENSION PV(NASHT,NASHT,NNASHX,NZ,3)
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('Output from RADDFQ',-1)
         WRITE(LUPRI,'(A,I3)') 'NXY       = ',NXY
         WRITE(LUPRI,'(A,I3)') 'IX        = ',IX
         WRITE(LUPRI,'(A,I3)') 'IY        = ',IY
         WRITE(LUPRI,'(A,I3)') 'IREPIJ    = ',IREPIJ
         WRITE(LUPRI,'(A)') 'FQ on entry'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
C     FIXME! Only implemented for real groups...
C
      IF (NZ .NE. 1) THEN
         CALL QUIT('*** AARRGGH, some lazy person did'
     &        //' not implement NZ>1')
      END IF
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
C     Output section
C
      IF (IPRINT .GE. 0) THEN
         WRITE(LUPRI,'(A)') 'FQ on exit'
         CALL PRQMAT(FQ,NORBT,NASHT,NORBT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
      END IF
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tknanobl */
      SUBROUTINE TKNONABL(H2FULL,NKL,H2XY,H2YX,H2AC,NXY,
     &     IREPIJ,NSTR,ONLYXY,DOFQ,DOH2AC)
C*****************************************************************************
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
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
C
#include <dgroup.h>     
#include <dcborb.h> 
#include <dcbopt.h> 
#include <dcbidx.h>
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
      IF (JTRLVL .EQ. 0) THEN
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
C        (gg|gu) or (gg|gg) integrals transformation.
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
C                 Get active indeces from general/active indeces
C
                  IF (JTRLVL .EQ. 1) THEN
                     IACP = IDXG2U(IPP)
                     IACV = IVV
                  ELSE 
                     IACP = IDXG2U(IPP)
                     IACV = IDXG2U(IVV)
                  END IF
C
C                 Elements for H2AC
C
                  IF (DOH2AC) THEN
                     IF (IACP .GT. 0 .AND. IACV .GT. 0) THEN
                        write(6,*) 'iacp iacv nxy',iacp,iacv,nxy,
     &                       h2full(ipv,1,1)
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
                              H2XY(IPP,IACV,IZ,I3) = 
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
                              H2YX(IVV,IACP,IZ,I3) = 
     &                             H2FULL(IPV,IZ,I3)
                           END DO
                        END DO
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck cl2nz3 */
      SUBROUTINE CL2NZ3(VCL,VNZ3,NELEM,IREPIJ)
C*****************************************************************************
C
C     Transform from 4INDXXXX0-integral-file format to (NZ,3) format.
C
C     Input:
C        VCL   : elements in (*,NCLASS) format
C        NELEM : # of elements
C        IREPIJ: irep of VCL
C       
C     Output:
C        VNZ3  : elements in (*,NZ,3) format
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>     
C
      DIMENSION VCL(NELEM,NZ,NZ,*)
      DIMENSION VNZ3(NELEM,NZ,3)
      DIMENSION DINT(2,8)
C
      CALL DZERO(DINT,16)
C
      DO I = 1, NELEM
         DO ICLASS = 1, 3
            IREPABR = 0
            DO IREPAB = 0, NBSYM - 1
               IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ1 = 1, NZ
                     IQ1 = IPQTOQ(IZ1,IREPAB)
                     DO IZ2 = 1,NZ
                        IQ2 = IPQTOQ(IZ2,IREPAB)
                        IFAC = IFACNZ3(ICLASS,1,IQ1,IQ2)
                        IF (IFAC .NE. 0) THEN
                           VNZ3(I,1,ICLASS) = VNZ3(I,1,ICLASS) + 
     &                          IFAC * VCL(I,IZ2,IZ1,IREPABR)
                        END IF
                        IF (NZ .GT. 1) THEN
                           IFAC = IFACNZ3(ICLASS,2,IQ1,IQ2)
                           IF (IFAC .NE. 0) THEN
                              VNZ3(I,2,ICLASS) = VNZ3(I,2,ICLASS) + 
     &                             IFAC * VCL(I,IZ2,IZ1,IREPABR)
                           END IF
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END DO
         IF (NZ .EQ. 4) THEN
            DO ICLASS = 1, 3
               IREPABR = 0
               DO IREPAB = 0, NBSYM - 1
                  IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                     IREPABR = IREPABR + 1
                     DO IZ1 = 1, NZ
                        IQ1 = IPQTOQ(IZ1,IREPAB)
                        DO IZ2 = 1,NZ
                           IQ2 = IPQTOQ(IZ2,IREPAB)
                           IFAC = IFACNZ3(ICLASS,3,IQ1,IQ2)
                           IF (IFAC .NE. 0) THEN
                              VNZ3(I,3,ICLASS) = VNZ3(I,3,ICLASS) + 
     &                             IFAC * VCL(I,IZ2,IZ1,IREPABR)
                           END IF
                           IFAC = IFACNZ3(ICLASS,4,IQ1,IQ2)
                           IF (IFAC .NE. 0) THEN
                              VNZ3(I,4,ICLASS) = VNZ3(I,4,ICLASS) + 
     &                             IFAC * VCL(I,IZ2,IZ1,IREPABR)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF
C
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ifacnz3 */
      FUNCTION IFACNZ3(I3,IZ,IQ1,IQ2)
C*****************************************************************************
C
C     Factors for integrals.
C
C     Written by Joern Thyssen - Nov 15 2000
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
C
      COMMON /SMOELF/ IFAC(4,3,4,4)
C
      LOGICAL FIRST
      DATA FIRST / .TRUE. /
      SAVE FIRST
C
      IF (FIRST) THEN
C        Initialize IFAC array
         FIRST = .FALSE.
         CALL IZERO(IFAC,3*4*4*4)
C
C        IZ = 1
C
         IFAC(1,1,1,1) = +1
         IFAC(1,1,2,2) = -1
C
         IFAC(1,2,3,3) = +1
         IFAC(1,2,4,4) = -1
C
         IFAC(1,3,3,3) = -1
         IFAC(1,3,4,4) = -1
C
         IF (NZ .GE. 2) THEN
            IFAC(2,1,1,2) = +1
            IFAC(2,1,2,1) = +1
C
            IFAC(2,2,3,4) = +1
            IFAC(2,2,4,3) = +1
C
            IFAC(2,3,3,4) = -1
            IFAC(2,3,4,3) = +1
         END IF
         IF (NZ .EQ. 4) THEN
            IFAC(3,1,1,3) = +1
            IFAC(3,1,2,4) = -1
C
            IFAC(3,2,3,1) = +1
            IFAC(3,2,4,2) = -1
C
            IFAC(3,3,3,1) = -1
            IFAC(3,3,4,2) = -1
C
            IFAC(4,1,1,4) = +1
            IFAC(4,1,2,3) = +1
C
            IFAC(4,2,3,2) = +1
            IFAC(4,2,4,1) = +1
C
            IFAC(4,3,3,2) = -1
            IFAC(4,3,4,1) = +1
         END IF
      END IF
C
      IFACNZ3 = IFAC(IZ,I3,IQ1,IQ2)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndrecno */
      SUBROUTINE FNDRECNO(IX,IY,IDXKR12,NFPCK12,NFPCK12T,
     &     IREPIJ,IREC)
C*****************************************************************************
C
C     Find record number for distribution (*,*,IX,IY) integrals.
C
C
C     Input:
C        IX, IY  : active orbital indeces
C        IDXKR12 : gives orbitals for compound index NXY
C        NFPCK12 : number of compund indeces for each symmetry
C        NFPCK12T: sum over NFPCK12
C       
C     Output:
C        IREPIJ  : irrep for compound (IX,IY)
C        IREC    : record number
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcbidx.h>
C
      DIMENSION NFPCK12(2)
      DIMENSION IDXKR12(2,*)
C
      write(6,*) 'active orbitals',IX,IY
      IF (JTRLVL .GE. 0 .AND. JTRLVL .LE. 2) THEN
C        (uu|uu), (gu|uu), or (gg|uu) integrals
         IGX = IX
         IGY = IY
      ELSE IF (JTRLVL .EQ. 3) THEN
C        (gg|gu) integrals
         IGX = IDXU2G(IX)
         IGY = IY
      ELSE
C        (gg|gg) integrals
         IGX = IDXU2G(IX)
         IGY = IDXU2G(IY)
      END IF
C
      write(6,*) '"general" orbitals',IGX,IGY
C
      NIJ = 0
      DO IFSYM = 1, NFSYM
         DO IJ = 1, NFPCK12(IFSYM)
            NIJ = NIJ + 1
            I = IDXKR12(1,NIJ)
            J = IDXKR12(2,NIJ)
            IF ( I .EQ. IGX .AND. J .EQ. IGY) THEN
               IREPIJ = IFSYM
               IREC = NIJ
               RETURN
            END IF
         END DO
      END DO
      WRITE(LUPRI,'(A,/,A)')
     &     '*** ERROR in FNDRECNO ***',
     &     'Internal error...'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gosci_rgeth2 */
      SUBROUTINE REAINTS(VUUUU,VGUUU,VGGUU,VGUGU,NDIMU,NDIMG,
     &     IFLAG,ITRLVL,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Read integrals and distribute in VUUUU, VGUUU, VGGUU and VGUGU.
C
C     IFLAG: bit 1: generate VUUUU
C            bit 2: generate VGUUU
C            bit 3: generate VGGUU
C            bit 4: generate VGUGU
C
C     Input:
C
C     Output:
C        VUUUU - integrals of type (uu|uu) 
C        VGUUU - integrals of type (gu|uu) 
C        VGGUU - integrals of type (gg|uu) 
C        VGUGU - integrals of type (gu|gu)
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
C
#include <dcbopt.h>
#include <dgroup.h>
#include <dcbidx.h>
C
      DIMENSION VUUUU(NDIMU,NDIMU,NDIMU,NDIMU,*)
      DIMENSION VGUUU(NDIMG,NDIMU,NDIMU,NDIMU,*)
      DIMENSION VGGUU(NDIMG,NDIMG,NDIMU,NDIMU,*)
      DIMENSION VGUGU(NDIMG,NDIMU,NDIMG,NDIMU,*)
      DIMENSION WORK(*)
C     index arrays
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
C     logicals
      LOGICAL LBIT
      LOGICAL LUUUU,LGUUU,LGGUU,LGUGU
      LOGICAL LDZERO
C     strings
      CHARACTER CDATE*8,TIMEX*8,FNAM*5,FNODE*10,LFNAM*10
C
      DIMENSION ICLFAC(16), ICLFAC4(4)
C
      DATA ICLFAC4 / +1,-1,-1,+1 /
C
#include <sdpre.h>
#include <memint.h>
C
      CALL QENTER('REAINTS')
C
      write(6,*) 'iprnt',iprnt
      IPRINT = 4
C
C     Initialize ICLFAC
C
      CALL IZERO(ICLFAC,16)
      IF (NZ .EQ. 1) THEN
         DO I = 1, 4
            ICLFAC(I) = ICLFAC4(I)
         END DO
      ELSE IF (NZ .EQ. 2) THEN
         DO I = 0,3
            ICLFAC(2*I+1) = ICLFAC4(I+1)
            ICLFAC(2*I+2) = ICLFAC4(I+1)
         END DO
      ELSE
         DO I = 0,3
            ICLFAC(4*I+1) = ICLFAC4(I+1)
            ICLFAC(4*I+2) = ICLFAC4(I+1)
            ICLFAC(4*I+3) = ICLFAC4(I+1)
            ICLFAC(4*I+4) = ICLFAC4(I+1)
         END DO
      END IF
C
      IF (IFLAG .LT. 0) THEN
         LDZERO = .FALSE.
         IFLAG = -IFLAG
      ELSE
         LDZERO = .TRUE.
      END IF
C
      IF ( NDIMG .GE. 1000 ) 
     &     CALL QUIT('Internal error in REAINTS (NDIMG>1000)')
C
      LUUUU = LBIT(IFLAG,1)
      LGUUU = LBIT(IFLAG,2)
      LGGUU = LBIT(IFLAG,3)
      LGUGU = LBIT(IFLAG,4)
C
      IF (IPRINT .GE. 3) THEN
         CALL HEADER('Output from REAINTS',-1)
         WRITE(LUPRI,'(1X,A,I2)') '- integral flag (IFLAG): ',IFLAG
         WRITE(LUPRI,'(1X,A,I2)')
     &        '- integral transformation level (ITRLVL): ',ITRLVL
         IF (LUUUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (uu|uu) integrals into VUUUU'
         IF (LGUUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gu|uu) integrals into VGUUU'
         IF (LGGUU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gg|uu) integrals into VGGUU'
         IF (LGUGU)
     &        WRITE(LUPRI,'(1X,A)')
     &        '- distribute (gu|gu) integrals into VGUGU'
      END IF
C
C
C     Zero out integral arrays and check for consistency
C
      IF (LUUUU .AND. LDZERO) CALL DZERO(VUUUU, NDIMU**4 * NZ**2 )
      IF (LGUUU) THEN
         IF ( ITRLVL .GT. -1 .AND. ITRLVL .LT. 1 ) THEN
            WRITE(LUPRI,'(A,/,A,I3,A)')
     &           '*** ERROR in REAINTS ***',
     &           'Integral level = ',ITRLVL,
     &           ' but (gu|uu) integrals requested!'
            CALL QUIT('*** ERROR in REAINTS ***')
         END IF
         IF (LDZERO) CALL DZERO(VGUUU, NDIMG * NDIMU**3 * NZ**2 )
      END IF
      IF (LGGUU .OR. LGUGU) THEN
         IF ( ITRLVL .GT. -1 .AND. ITRLVL .LT. 2 ) THEN
            WRITE(LUPRI,'(A,/,A,I3,A)')
     &           '*** ERROR in REAINTS ***',
     &           'Integral level = ',ITRLVL,
     &           ' but (gg|uu) and/or (gu|gu) integrals requested!'
            CALL QUIT('*** ERROR in REAINTS ***')
         END IF
         IF (LDZERO) THEN
            CALL DZERO(VGGUU, NDIMG**2 * NDIMU**2 * NZ**2 )
            CALL DZERO(VGUGU, NDIMG**2 * NDIMU**2 * NZ**2 )
         END IF
      END IF
C
C     Make index vector general index --> active index
C
      IF ( ITRLVL .GT. 0 ) THEN
         NGENERAL = ITRLVL
      ELSE
         NGENERAL = 4
      END IF
#ifdef UNDEF
      write(6,*) 'idxg2u = ',(IDXG2U(I),I=1,NDIMG)
#endif
C
C     Read integral info from MDINFO
C
      LUINT = 1 
      OPEN (LUINT, FILE = 'MDINFO', STATUS = 'OLD',
     &     FORM = 'UNFORMATTED' )
C
C     Read time & date for creation
C
      READ (LUINT) DATEX, TIMEX
C
      IF ( IPRINT .GE. 3 ) WRITE(LUPRI,'(1X,4A)')
     &     '- time and date read from MDINFO: ',
     &     TIMEX, ' ',DATEX
C
      READ (LUINT) NFPCK12, NFPCK12T
      READ (LUINT) NFPCK34, NFPCK34T
      READ (LUINT) NSTR
C     
      CALL MEMGET('INTE',KDXKR12,2*NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB12,NFPCK12T,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KDXB34,NFPCK34T,WORK,KFREE,LFREE)
C
C     I and J indeces
      CALL READI(LUINT,NFPCK12T*2,WORK(KDXKR12))
#ifdef UNDEF      
      do i = 1, 2 * nfpck12t
         write(6,*) i,IGETINT(WORK(KDXKR12),I)
      end do
#endif
C     Boson irrep for IJ indeces
      CALL READI(LUINT,NFPCK12T,WORK(KDXB12))
C     Boson irrep for KL indeces
      CALL READI(LUINT,NFPCK34T,WORK(KDXB34))
C
      CLOSE (LUINT, STATUS = 'KEEP')
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Open integral files.
C
      LUINT = 0
      DO IREPIJ = 1, NFSYM
         IRECL = IRECLEN(NCLASS * NFPCK34 (IREPIJ), 0, 0)
         IRECL = MAX0(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)
         OPEN (LUINT+IREPIJ,FILE=FNODE,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=IRECL)
      END DO
C
C     Read and distribute integrals.
C
C     - IJ is triangular.
C     - KL is not triangular.
C
      NBSYMP = NBSYM / NFSYM
C
      NIJOFF = 0
      NKLOFF = 0
      DO IREPIJ = 1, NFSYM
#ifdef UNDEF      
         write(6,*) 'REAINT IREPIJ = ',IREPIJ
#endif
         NIJ = NFPCK12(IREPIJ)
         NKL = NFPCK34(IREPIJ)
#ifdef UNDEF      
         write(6,*) 'REAINT NIJ = ',NIJ
         write(6,*) 'REAINT NKL = ',NKL
#endif
         CALL MEMGET('REAL',KHMAT,NKL*NZ*NZ*NBSYMP,WORK,KFREE,LFREE)
         DO IJ = 1, NIJ
#ifdef UNDEF      
            write(6,*) '  REAINT IJ = ', IJ
#endif
            IKR = IGETINT(WORK(KDXKR12),2*(IJ-1)+1 + NIJOFF)
            JKR = IGETINT(WORK(KDXKR12),2*IJ + NIJOFF)
            IREPIJB = IGETINT(WORK(KDXB12),IJ + NIJOFF)
#ifdef UNDEF      
            write(6,*) '  REAINT IKR = ', IKR
            write(6,*) '  REAINT JKR = ', JKR
            write(6,*) '  REAINT IREPIJB = ', IREPIJB
            call flshfo(6)
#endif
            READ (LUINT+IREPIJ,REC=IJ) (WORK(KHMAT+I-1),I=1,NKL*NCLASS)
c#ifdef UNDEF      
            write(6,*) 'int ',(work(khmat+i-1),i=1,nkl*nclass)
c#endif
C
            NZKL = 0
            IKR1 = IKR
            JKR1 = JKR
            KL = 0
            LL = 0
            DO IREPL = 1, NFSYM
               IREPK = MOD(IREPL + IREPIJ, 2 ) + 1
               DO L = 1, NSTR(IREPL,0,4)
                  LL = LL + 1
                  DO K = 1, NSTR(IREPK,0,3)
                     IF (IREPK .EQ. 1) THEN
                        KK = K
                     ELSE
                        KK = K + NSTR(1,0,3)
                     END IF
                     KL = KL + 1
C
#ifdef UNDEF
                     WRITE(6,'(A,I2,A,I2,A,I2,A,I2,A,4(3X,D14.8))')
     &                    '    INT(',IKR,',',JKR,',',KK,',',LL,') = ',
     &                    (WORK(KHMAT + KL + (I-1)*NKL - 1),I=1,NCLASS)
#endif
C
C                    Get active indeces from general/active indeces.
C
                     IF (ITRLVL .EQ. 0) THEN
                        IAC = IKR
                        JAC = JKR
                        KAC = KK
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 1) THEN
                        IAC = IKR
                        JAC = JKR
                        KAC = IDXG2U(KK)
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 2) THEN
C                       ...not implemented...
                        call quit('ralle ralle')
                        IAC = IDXG2U(IKR)
                        JAC = IDXG2U(JKR)
                        KAC = KK
                        LAC = LL
                     ELSE IF (ITRLVL .EQ. 3) THEN
                        IAC = IDXG2U(IKR)
                        JAC = JKR
                        KAC = IDXG2U(KK)
                        LAC = IDXG2U(LL)
                     ELSE 
                        IAC = IDXG2U(IKR)
                        JAC = IDXG2U(JKR)
                        KAC = IDXG2U(KK)
                        LAC = IDXG2U(LL)
                     END IF
C
C                    Distribute integrals
C                    
                     DO I = 1, NCLASS
                        IF (LUUUU .AND.
     &                       IAC .GT. 0 .AND. JAC .GT. 0 .AND.
     &                       KAC .GT. 0 .AND. LAC .GT. 0) THEN
                           VUUUU(KAC,LAC,IAC,JAC,I) =
     &                          WORK(KHMAT + KL + (I-1)*NKL - 1)
                           VUUUU(KAC,LAC,JAC,IAC,I) =
     &                          ICLFAC(I) *
     &                          WORK(KHMAT + KL + (I-1)*NKL - 1)
                        END IF
                        IF (LGUUU) THEN
                           IF (JAC .GT. 0 .AND. IAC .GT. 0
     &                          .AND. LAC .GT. 0) THEN
                              VGUUU(KK,LAC,IAC,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                           IF (IAC .GT. 0 .AND. LAC .GT. 0
     &                          .AND. JAC .GT. 0) THEN
                              VGUUU(KK,LAC,JAC,IAC,I) =
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                        IF (LGGUU) THEN
                           IF (IAC .GT. 0 .AND. IAC .GT. 0) THEN
                              VGGUU(KK,LL,IAC,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                              VGGUU(KK,LL,JAC,IAC,I) = 
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                        IF (LGUGU) THEN
                           IF (JAC .GT. 0 .AND. LAC .GT. 0) THEN
                              VGUGU(KK,LAC,IKR,JAC,I) = 
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                           IF (IAC .GT. 0 .AND. LAC .GT. 0) THEN
                              VGUGU(KK,LAC,JKR,IAC,I) =
     &                             ICLFAC(I) *
     &                             WORK(KHMAT + KL + (I-1)*NKL - 1)
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Close files
C
      DO IREPIJ = 1,NFSYM
         CLOSE (LUINT+IREPIJ, STATUS = 'KEEP' )
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,4(D14.7,1X),3(/,16X,4(D14.7,1X)))
C
      CALL MEMREL('REAINTS',WORK,KWORK,KWORK,KFREE,LFREE)
C
      IF (IPRINT .GE. 3) THEN
         IF (LUUUU) THEN
            CALL HEADER('VUUUU matrix',-1)
            DO I = 1,NDIMU
               DO J = 1,NDIMU
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VUUUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGUUU) THEN
            CALL HEADER('VGUUUU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMU
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGUUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGGUU) THEN
            CALL HEADER('VGGUU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMG
                  DO K = 1,NDIMU
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGGUU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (LGUGU) THEN
            CALL HEADER('VGUGU matrix',-1)
            DO I = 1,NDIMG
               DO J = 1,NDIMU
                  DO K = 1,NDIMG
                     DO L = 1,NDIMU
                        WRITE(LUPRI,9000)
     &                       I,J,K,L,(VGUGU(I,J,K,L,M),M=1,NCLASS)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END IF
C
      CALL QEXIT('REAINTS')
C
      RETURN
10000 STOP 'END READING MDCINT'
10010 STOP 'ERROR READING MDCINT'
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck igetint */
      INTEGER FUNCTION IGETINT(IX,I)
C*****************************************************************************
C
C     Return element I of array IX
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION IX(*)
C
      IGETINT = IX(I)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck D2MUUUU */
      SUBROUTINE D2MUUUU(VDUUUU,VMUUUU,IPRINT)
C*****************************************************************************
C
C     Transform a UUUU matrix in Dirac format to a UUUU matrix
C     in Molfdir format.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION VDUUUU(NASHT,NASHT,NASHT,NASHT,NZ,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NFPCK12(2), NFPCK34(2)
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
C
#include <ibtfun.h>
#include <memint.h>
C
      CALL DZERO(VMUUUU,(2*NASHT)**4 * MIN0(2,NZ))
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DEFCL(NCLASS)
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
      DO IREPI = 1, NFSYM
      DO II = 1+IASH(IREPI),IASH(IREPI)+NASH(IREPI)
         DO IREPJ = 1, NFSYM
         IREPIJ = IBTXOR(IREPI-1,IREPJ-1) + 1
         DO JJ = 1+IASH(IREPJ),IASH(IREPJ)+NASH(IREPJ)
            DO KK = 1, NASHT
               DO LL = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C                 
                  CALL DZERO(DINT,16)
                  DO ICLASS = 1, 4
                     IREPABR = 0
                     DO IREPAB = 0, NBSYM - 1
                        IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                           IREPABR = IREPABR + 1
                           DO IZ1 = 1, NZ
                              IQ1 = IPQTOQ(IZ1,IREPAB)
                              DO IZ2 = 1,NZ
                                 IQ2 = IPQTOQ(IZ2,IREPAB)
                                 IFAC = INTCL(ICLASS,IQ1,IQ2)
#ifdef UNDEF
                                 if (ii .eq. 1 .and. jj .eq. 2 .and.
     &                                kk .eq. 1 .and. ll .eq. 2 ) then
                                    write(6,*) 'iz1,iz2,iq1,iq2',
     &                                   iz1,iz2,iq1,iq2,iclass,ifac
                                    write(6,*) 'dint before',
     &                                   dint(iclass,1)
                                    write(6,*) 'vd',vduuuu(ii,jj,kk,ll,
     &                                   iz2,iz1,irepabr)
                                 end if
#endif
                                 IF (IFAC .NE. 0) THEN
                                    DINT(ICLASS,1) = DINT(ICLASS,1) +
     &                                   IFAC * VDUUUU(II,JJ,KK,LL,
     &                                   IZ2,IZ1,IREPABR)
                                 END IF
#ifdef UNDEF
                                 if (ii .eq. 1 .and. jj .eq. 2 .and.
     &                                kk .eq. 1 .and. ll .eq. 2 ) then
                                    write(6,*) 'dint after',
     &                                   dint(iclass,1)
                                 end if
#endif
                                 IF (NZ .GT. 1) THEN
                                    IFAC = INTCL(ICLASS+4,IQ1,IQ2)
                                    IF (IFAC .NE. 0) THEN
                                       DINT(ICLASS,2) = DINT(ICLASS,2) +
     &                                      IFAC * VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR)
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
                  IF (NZ .EQ. 4) THEN
                     DO ICLASS = 1, 4
                        IREPABR = 0
                        DO IREPAB = 0, NBSYM - 1
                           IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                              IREPABR = IREPABR + 1
                              DO IZ1 = 1, NZ
                                 IQ1 = IPQTOQ(IZ1,IREPAB)
                                 DO IZ2 = 1,NZ
                                    IQ2 = IPQTOQ(IZ2,IREPAB)
                                    IFAC = INTCL(ICLASS+8,IQ1,IQ2)
                                    IF (IFAC .NE. 0) THEN
                                       DINT(ICLASS+4,1) =
     &                                      DINT(ICLASS+4,1) +
     &                                      IFAC * VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR)
                                    END IF
                                    IFAC = INTCL(ICLASS+12,IQ1,IQ2)
                                    IF (IFAC .NE. 0) THEN
                                       DINT(ICLASS+4,2) =
     &                                      DINT(ICLASS+4,2) +
     &                                      IFAC * VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR)
                                    END IF
                                 END DO
                              END DO
                           END IF
                        END DO
                     END DO
                  END IF
C
C                 Distribute integrals
C                 
#ifdef UNDEF
                  write(6,'(1P,A,8I3,/,3(4D15.8/),4D15.8)')
     &                 'd2u: ',i,j,k,l,itr,jtr,ktr,ltr,dint
#endif
                  DO IZ = 1,2
                     IS = -2*IZ+3
                     VMUUUU(I,J,K,L,IZ) = DINT(1,IZ)
                     VMUUUU(ITR,J,K,L,IZ) = DINT(7,IZ)
                     VMUUUU(I,JTR,K,L,IZ) = DINT(8,IZ)
                     VMUUUU(ITR,JTR,K,L,IZ) = IS * DINT(2,IZ)
                     VMUUUU(I,J,KTR,L,IZ) = DINT(5,IZ)
                     VMUUUU(ITR,J,KTR,L,IZ) = IS * DINT(3,IZ)
                     VMUUUU(I,JTR,KTR,L,IZ) = DINT(4,IZ)
                     VMUUUU(ITR,JTR,KTR,L,IZ) = -IS * DINT(6,IZ)
c
                     VMUUUU(I,J,K,LTR,IZ) = DINT(6,IZ)
                     VMUUUU(ITR,J,K,LTR,IZ) = IS * DINT(4,IZ)
                     VMUUUU(I,JTR,K,LTR,IZ) = DINT(3,IZ)
                     VMUUUU(ITR,JTR,K,LTR,IZ) = -IS * DINT(5,IZ)
                     VMUUUU(I,J,KTR,LTR,IZ) = DINT(2,IZ)
                     VMUUUU(ITR,J,KTR,LTR,IZ) = -IS * DINT(8,IZ)
                     VMUUUU(I,JTR,KTR,LTR,IZ) = -IS * DINT(7,IZ)
                     VMUUUU(ITR,JTR,KTR,LTR,IZ) = IS * DINT(1,IZ)
                  END DO
               END DO
            END DO
            END DO
         END DO
      END DO
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('VUUUU matrix in molfdir format',-1)
         DO I = 1,2*NASHT
            DO J = 1,2*NASHT
               DO K = 1,2*NASHT
                  DO L = 1,2*NASHT
                     IF ( ABS(VMUUUU(I,J,K,L,1)) .GT. 0 .OR.
     &                    ABS(VMUUUU(I,J,K,L,2)) .GT. 0) 
     &                    WRITE(LUPRI,9000)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,MIN0(2,NZ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck m2duuuu */
      SUBROUTINE M2DUUUU(VMUUUU,VDUUUU,IPRINT)
C*****************************************************************************
C
C     Transform a UUUU matrix in Molfdir format to a UUUU matrix
C     in Dirac format.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION VDUUUU(NASHT,NASHT,NASHT,NASHT,NZ,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
C
#include <ibtfun.h>
#include <memint.h>
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      write(6,*) 'nclass',nclass,nasht,nasht**4 * nclass
      call flshfo(6)
      CALL DZERO(VDUUUU,NASHT**4 * NCLASS)
      write(6,*) 'after dzero'
      call flshfo(6)
C
      CALL DEFINVCL(NCLASS)
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
      DO IREPI = 1, NFSYM
      DO II = 1+IASH(IREPI),IASH(IREPI)+NASH(IREPI)
         DO IREPJ = 1, NFSYM
         IREPIJ = IBTXOR(IREPI-1,IREPJ-1) + 1
         DO JJ = 1+IASH(IREPJ),IASH(IREPJ)+NASH(IREPJ)
            DO KK = 1, NASHT
               DO LL = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C                 
                  CALL DZERO(DINT,16)
C
C                 Calculate integral classes.
C                 
                  DO IZ = 1,2
                     DINT(1,IZ) = VMUUUU(I,J,K,L,IZ)
                     DINT(2,IZ) = VMUUUU(I,J,KTR,LTR,IZ)
                     DINT(3,IZ) = VMUUUU(I,JTR,K,LTR,IZ)
                     DINT(4,IZ) = VMUUUU(I,JTR,KTR,L,IZ)
                     DINT(5,IZ) = VMUUUU(I,J,KTR,L,IZ)
                     DINT(6,IZ) = VMUUUU(I,J,K,LTR,IZ)
                     DINT(7,IZ) = VMUUUU(ITR,J,K,L,IZ)
                     DINT(8,IZ) = VMUUUU(I,JTR,K,L,IZ)
                  END DO
                  DO ICLASS = 1, 4
                     IREPABR = 0
                     DO IREPAB = 0, NBSYM - 1
                        IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                           IREPABR = IREPABR + 1
                           DO IZ1 = 1, NZ
                              IQ1 = IPQTOQ(IZ1,IREPAB)
                              DO IZ2 = 1,NZ
                                 IQ2 = IPQTOQ(IZ2,IREPAB)
                                 RFAC = RINVINTCL(ICLASS,IQ1,IQ2)
#ifdef UNDEF
                                 if (ii .eq. 1 .and. jj .eq. 2 .and.
     &                                kk .eq. 1 .and. ll .eq. 2 ) then
                                    write(6,*) 'iz1,iz2,iq1,iq2',
     &                                   iz1,iz2,iq1,iq2,iclass,ifac
                                    write(6,*) 'dint before',
     &                                   dint(iclass,1)
                                    write(6,*) 'vd',vduuuu(ii,jj,kk,ll,
     &                                   iz2,iz1,irepabr)
                                 end if
#endif
                                 IF (RFAC .NE. 0) THEN
                                    VDUUUU(II,JJ,KK,LL,
     &                                   IZ2,IZ1,IREPABR) =
     &                                   VDUUUU(II,JJ,KK,LL,
     &                                   IZ2,IZ1,IREPABR) +
     &                                   RFAC * DINT(ICLASS,1)
                                 END IF
#ifdef UNDEF
                                 if (ii .eq. 1 .and. jj .eq. 2 .and.
     &                                kk .eq. 1 .and. ll .eq. 2 ) then
                                    write(6,*) 'dint after',
     &                                   dint(iclass,1)
                                 end if
#endif
                                 IF (NZ .GT. 1) THEN
                                    RFAC = RINVINTCL(ICLASS+4,IQ1,IQ2)
                                    IF (RFAC .NE. 0) THEN
                                       VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) =
     &                                      VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) +
     &                                      RFAC * DINT(ICLASS,2)
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
                  IF (NZ .EQ. 4) THEN
                     DO ICLASS = 1, 4
                        IREPABR = 0
                        DO IREPAB = 0, NBSYM - 1
                           IF (IREPIJ .EQ. JBTOF(IREPAB,1)) THEN
                              IREPABR = IREPABR + 1
                              DO IZ1 = 1, NZ
                                 IQ1 = IPQTOQ(IZ1,IREPAB)
                                 DO IZ2 = 1,NZ
                                    IQ2 = IPQTOQ(IZ2,IREPAB)
                                    RFAC = RINVINTCL(ICLASS+8,IQ1,IQ2)
                                    IF (RFAC .NE. 0) THEN
                                       VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) =
     &                                      VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) +
     &                                      RFAC * DINT(ICLASS+4,1)
                                    END IF
                                    RFAC = RINVINTCL(ICLASS+12,IQ1,IQ2)
                                    IF (RFAC .NE. 0) THEN
                                       VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) =
     &                                      VDUUUU(II,JJ,KK,LL,
     &                                      IZ2,IZ1,IREPABR) +
     &                                      RFAC * DINT(ICLASS+4,2)
                                    END IF
                                 END DO
                              END DO
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            END DO
         END DO
      END DO
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,4(D14.7,1X),3(/,16X,4(D14.7,1X)))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('VUUUU matrix in Dirac format',-1)
         DO I = 1,NASHT
            DO J = 1,NASHT
               DO K = 1,NASHT
                  DO L = 1,NASHT
                     WRITE(LUPRI,9000)
     &                    I,J,K,L,
     &                    (((VDUUUU(I,J,K,L,IZ1,IZ2,ICL),
     &                    IZ1 = 1,NZ),IZ2=1,NZ),ICL=1,NCLASS)
                  END DO
               END DO
            END DO
         END DO
      END IF
      call quit('ralle ralle')
C
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck m2dnz3 */
      SUBROUTINE M2DNZ3(VMUUUU,VDUUUU,IPRINT)
C*****************************************************************************
C
C     Transform a UUUU matrix in Molfdir format to a UUUU matrix
C     in Dirac (NZ,3) format.
C
C     Written by J. Thyssen - Jul 21 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
C
#include <ibtfun.h>
#include <memint.h>
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VDUUUU,NASHT * NASHT * NNASHX * NZ * 3)
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
C
      IJ = 0
      DO II = 1, NASHT
         DO JJ = 1, II
            IJ = IJ + 1
            DO KK = 1, NASHT
               DO LL = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C
C                 Calculate integral classes.
C
                  VDUUUU(KK,LL,IJ,1,1) = VMUUUU(K,L,I,J,1)
                  VDUUUU(KK,LL,IJ,1,2) = VMUUUU(K,LTR,I,JTR,1)
                  VDUUUU(KK,LL,IJ,1,3) = VMUUUU(KTR,L,I,JTR,1)
C
                  IF (NZ .GE. 2) THEN
                     VDUUUU(KK,LL,IJ,2,1) = VMUUUU(K,L,I,J,2)
                     VDUUUU(KK,LL,IJ,2,2) = VMUUUU(K,LTR,I,JTR,2)
                     VDUUUU(KK,LL,IJ,2,3) = VMUUUU(KTR,L,I,JTR,2)
                  END IF
C
                  IF (NZ .EQ. 4) THEN
                     VDUUUU(KK,LL,IJ,3,1) = VMUUUU(I,J,I,JTR,1)
                     VDUUUU(KK,LL,IJ,4,1) = VMUUUU(I,J,I,JTR,2)
                     VDUUUU(KK,LL,IJ,3,2) = VMUUUU(K,LTR,I,J,1)
                     VDUUUU(KK,LL,IJ,4,2) = VMUUUU(K,LTR,I,J,2)
                     VDUUUU(KK,LL,IJ,3,3) = VMUUUU(KTR,L,I,J,1)
                     VDUUUU(KK,LL,IJ,4,3) = VMUUUU(KTR,L,I,J,2)
                  END IF
               END DO
            END DO
         END DO
      END DO
C
 9000 FORMAT(1P,1X,4I3,4X,3(D14.7,1X),3(/,16X,3(D14.7,1X)))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('VUUUU matrix in Dirac (NZ,3) format',-1)
         IJ = 0
         DO I = 1,NASHT
            DO J = 1,I
               IJ = IJ + 1
               DO K = 1,NASHT
                  DO L = 1,NASHT
                     WRITE(LUPRI,9000)
     &                    K,L,I,J,
     &                    ((VDUUUU(K,L,IJ,IZ,ICL),
     &                    IZ = 1,NZ),ICL=1,3)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dnz32m */
      SUBROUTINE DNZ32M(VDUUUU,VMUUUU,IPRINT)
C*****************************************************************************
C
C     Transform a UUUU matrix in Dirac (NZ,3) format to a UUUU matrix
C     in Molfdir format.
C
C     Written by J. Thyssen - Nov 20 2000
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION VDUUUU(NASHT,NASHT,NNASHX,NZ,*)
      DIMENSION VMUUUU(2*NASHT,2*NASHT,2*NASHT,2*NASHT,*)
C
      DIMENSION NSTR(2,0:2,4)
      DIMENSION KR(-300:300)
      DIMENSION DINT(8,2)
C
#include <ibtfun.h>
#include <memint.h>
C
      NCLASS = NZ * NZ * NBSYM / NFSYM
      CALL DZERO(VMUUUU, NASHT**4 * MIN0(NZ,2))
C
      II = 0
      JJ = 0
      DO IFSYM = 1, NFSYM
         DO I = 1, NASH(IFSYM)
            II = II + 1
            JJ = JJ + 1
            KR(II) = JJ
            KR(-II) = JJ + NASH(IFSYM)
         ENDDO
         JJ = JJ + NASH(IFSYM)
      ENDDO
 9000 FORMAT(1P,1X,4I3,4X,3(D14.7,1X),3(/,16X,3(D14.7,1X)))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('Output from DNZ32M',-1)
         CALL HEADER('Input: VUUUU matrix in Dirac (NZ,3) format',-1)
         IJ = 0
         DO I = 1,NASHT
            DO J = 1,I
               IJ = IJ + 1
               DO K = 1,NASHT
                  DO L = 1,NASHT
                     WRITE(LUPRI,9000)
     &                    K,L,I,J,
     &                    ((VDUUUU(K,L,IJ,IZ,ICL),
     &                    IZ = 1,NZ),ICL=1,3)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      DO II = 1, NASHT
         DO JJ = 1, NASHT
            IJ = MAX(II,JJ) * ( MAX(II,JJ) - 1 ) / 2 + MIN(II,JJ)
            write(6,*) 'II,JJ,IJ',II,JJ,IJ
            DO LL = 1, NASHT
               DO KK = 1, NASHT
                  I = KR(II)
                  J = KR(JJ)
                  K = KR(KK)
                  L = KR(LL)
                  ITR = KR(-II)
                  JTR = KR(-JJ)
                  KTR = KR(-KK)
                  LTR = KR(-LL)
C
C                 Canonical integrals.
C                 (kl|ij), (kL|iJ), (Kl|iJ)
C
                  VMUUUU(K,L,I,J,1) = VDUUUU(KK,LL,IJ,1,1)
                  VMUUUU(K,LTR,I,JTR,1) = VDUUUU(KK,LL,IJ,1,2)
                  VMUUUU(KTR,L,I,JTR,1) = VDUUUU(KK,LL,IJ,1,3)
C
C                 (LK|ij), (kl|JI), (LK|JI)
C
                  VMUUUU(LTR,KTR,I,J,1) = VDUUUU(KK,LL,IJ,1,1)
                  VMUUUU(K,L,JTR,ITR,1) = VDUUUU(KK,LL,IJ,1,1)
                  VMUUUU(LTR,KTR,JTR,ITR,1) = VDUUUU(KK,LL,IJ,1,1)
C
C                 (Lk|Ji)^*, (lK|Ji)^*
C
                  VMUUUU(LTR,K,JTR,I,1) = VDUUUU(KK,LL,IJ,1,2)
                  VMUUUU(L,KTR,JTR,I,1) = VDUUUU(KK,LL,IJ,1,3)
C
                  IF (NZ .GE. 2) THEN
C
C                    Canonical integrals.
C                    (kl|ij), (kL|iJ), (Kl|iJ)
C
                     VMUUUU(K,L,I,J,2) = VDUUUU(KK,LL,IJ,2,1)
                     VMUUUU(K,LTR,I,JTR,2) = VDUUUU(KK,LL,IJ,2,2)
                     VMUUUU(KTR,L,I,JTR,2) = VDUUUU(KK,LL,IJ,2,3)
C
C                    (LK|ij), (kl|JI), (LK|JI)
C
                     VMUUUU(LTR,KTR,I,J,2) = VDUUUU(KK,LL,IJ,2,1)
                     VMUUUU(K,L,JTR,ITR,2) = VDUUUU(KK,LL,IJ,2,1)
                     VMUUUU(LTR,KTR,JTR,ITR,2) = VDUUUU(KK,LL,IJ,2,1)
C
C                    (Lk|Ji)^*, (lK|Ji)^*
C
                     VMUUUU(LTR,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,2)
                     VMUUUU(L,KTR,JTR,I,2) = - VDUUUU(KK,LL,IJ,2,3)
C
                  END IF
C
                  IF (NZ .EQ. 4) THEN
C
C                    Canonical integrals
C                    (kl|iJ)
C
                     VMUUUU(K,L,I,JTR,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(K,L,I,JTR,2) = VDUUUU(KK,LL,IJ,4,1)
C
C                    (kL|ij) 
C
                     VMUUUU(K,LTR,I,J,1) = VDUUUU(KK,LL,IJ,3,2)
                     VMUUUU(K,LTR,I,J,2) = VDUUUU(KK,LL,IJ,4,2)
C
C                    (Kl|ij)
C
                     VMUUUU(KTR,L,I,J,1) = VDUUUU(KK,LL,IJ,3,3)
                     VMUUUU(KTR,L,I,J,2) = VDUUUU(KK,LL,IJ,4,3)
C
C                    The rest...
C
C                    (Kl|JI)
C
                     VMUUUU(KTR,L,JTR,ITR,1) = VDUUUU(KK,LL,IJ,3,3)
                     VMUUUU(KTR,L,JTR,ITR,2) = VDUUUU(KK,LL,IJ,4,3)
C
C                    (kL|JI) 
C
                     VMUUUU(K,LTR,ITR,JTR,1) = VDUUUU(KK,LL,IJ,3,2)
                     VMUUUU(K,LTR,ITR,JTR,2) = VDUUUU(KK,LL,IJ,4,2)
C
C                    (lk|Ji)^*, (KL|Ji)^*
C
                     VMUUUU(L,K,JTR,I,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(L,K,JTR,I,2) = - VDUUUU(KK,LL,IJ,4,1)
                     VMUUUU(KTR,LTR,JTR,I,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(KTR,LTR,JTR,I,2) = - VDUUUU(KK,LL,IJ,4,1)
C
C                    (LK|iJ)
C
                     VMUUUU(LTR,KTR,I,JTR,1) = VDUUUU(KK,LL,IJ,3,1)
                     VMUUUU(LTR,KTR,I,JTR,2) = VDUUUU(KK,LL,IJ,4,1)
C
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C
C
 9001 FORMAT(1P,1X,4I3,4X,2(D14.7,1X))
C
      IF (IPRINT .GE. 0) THEN
         CALL HEADER('Output: VUUUU matrix in molfdir format',-1)
         DO I = 1,2*NASHT
            DO J = 1,2*NASHT
               DO K = 1,2*NASHT
                  DO L = 1,2*NASHT
                     IF ( ABS(VMUUUU(I,J,K,L,1)) .GT. 0 .OR.
     &                    ABS(VMUUUU(I,J,K,L,2)) .GT. 0) 
     &                    WRITE(LUPRI,9001)
     &                    I,J,K,L,(VMUUUU(I,J,K,L,IZ),IZ=1,MIN0(2,NZ))
                  END DO
               END DO
            END DO
         END DO
      END IF
C
      RETURN
      END 




