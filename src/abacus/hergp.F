!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

C
C FILE    : abacus/hergp.F
C
C  /* Deck readsq */
      SUBROUTINE READSQ(IT,N,IX)
      DIMENSION IX(N)
      READ (IT,END=10) IX
      RETURN
   10 CONTINUE
      IX(N)=-1
      RETURN
      END
C  /* Deck writsq */
      SUBROUTINE WRITSQ(IT,N,IX)
      DIMENSION IX(N)
      WRITE (IT) IX
      RETURN
      END
C  /* Deck readda */
      SUBROUTINE READDA(IT,IU,N,IX)
      DIMENSION IX(N)
      READ(IT,REC=IU) IX
      RETURN
      END
C  /* Deck writda */
      SUBROUTINE WRITDA(IT,IU,N,IX)
      DIMENSION IX(N)
      WRITE(IT,REC=IU) IX
      RETURN
      END
C  /* Deck zergrd */
      SUBROUTINE ZERGRD
      use mxcent_mod, only: MXCOOR
#include "implicit.h"
#include "trkoor.h"
#include "taymol.h"
!      use trkoor_mod, only: NCOOR
!      use taymol_mod, only: GRDMOL
      GRDMOL(1:NCOOR) = 0.0D0
      RETURN
      END
C  /* Deck zerhes */
      SUBROUTINE ZERHES
      use mxcent_mod, only: MXCOOR
#include "implicit.h"
#include "trkoor.h"
#include "taymol.h"
!      use trkoor_mod, only: NCOOR
!      use taymol_mod, only: HESMOL
      HESMOL(1:NCOOR,1:NCOOR) = 0.0D0
      RETURN
      END
C  /* Deck addgrd */
      SUBROUTINE ADDGRD(GRAD)
      use mxcent_mod, only: MXCOOR
#include "implicit.h"
!      use trkoor_mod, only: NCOOR
!      use taymol_mod, only: GRDMOL
      INTEGER :: I
#include "trkoor.h"
      REAL*8 GRAD(NCOOR)
#include "taymol.h"
      DO I = 1, NCOOR
         GRDMOL(I) = GRDMOL(I) + GRAD(I)
      END DO
      RETURN
      END
C  /* Deck addhes */
      SUBROUTINE ADDHES(HESS)
      use mxcent_mod, only: MXCOOR
#include "implicit.h"
!      use taymol_mod, only: HESMOL
!      use trkoor_mod, only: NCOOR
      INTEGER :: I,J
#include "trkoor.h"
      DIMENSION HESS(MXCOOR,MXCOOR)
#include "taymol.h"
C
      DO 100 I = 1, NCOOR
         DO 200 J = 1, NCOOR
            HESMOL(J,I) = HESMOL(J,I) + HESS(J,I)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck prigrd */
      SUBROUTINE PRIGRD(GRAD,CSTRA,SCTRA)
C
C     tuh Mar 1985
C     tuh Jun 28 1988 - modified for symmetry
C
      use maxaqn_mod, only: MXAQN, MXQN
!      use nuclei_mod, only: NAMEX, NUCDEP
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use priunit_mod, only: LUPRI
!      use abainf_mod, only: DOSYM
      use symmet_mod, only: IPTCOR, MAXREP, NCRREP
#include "implicit.h"
      INTEGER :: I,J
C
#include "abainf.h"
#include "nuclei.h"
      DIMENSION GRAD(MXCOOR), CGRAD(MXCOOR), CSTRA(*), SCTRA(*)
      IF (.NOT.DOSYM(1)) THEN
         WRITE (LUPRI,'(2X,A)')
     &  ' Gradient not calculated - '//
     &  ' totally symmetric distortions not requested.'
      ELSE IF (MAXREP .EQ. 0) THEN
         IOFF = 0
         DO 100 IATOM = 1, NUCDEP
            WRITE (LUPRI,1000) NAMDEP(IATOM), (GRAD(IOFF+J),J=1,3)
            IOFF = IOFF + 3
  100    CONTINUE
      ELSE
         DO 200 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(25X,A6,F17.10)') NAMEX(IPTCOR(I,1)), GRAD(I)
  200    CONTINUE
         WRITE (LUPRI,'(//)')
         CALL TRAGRD(GRAD,CGRAD,CSTRA,SCTRA,NCRREP(0,1),3*NUCDEP)
         IOFF = 0
         DO 300 IATOM = 1, NUCDEP
            WRITE (LUPRI,1000) NAMDEP(IATOM),(CGRAD(IOFF+J),J=1,3)
            IOFF = IOFF + 3
  300    CONTINUE
      END IF
      WRITE (LUPRI,'(//)')
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT (1X,A6,F17.10,2F24.10)
      END
C  /* Deck prihes */
      SUBROUTINE PRIHES (SHESS,KEY,CSTRA,SCTRA)
C
C     Adapted from Nelson's OUTPAK for Hessians March 1985 tuh
C     Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
C     Jun 28 1988 tuh - modified for symmetry
C
!      use nuclei_mod, only: NAMEX, NATOMS, NUCDEP
      use symmet_mod, only: IPTCOR, MAXREP, NCRREP
      use priunit_mod, only: LUPRI
!      use abainf_mod, only: WRITE
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
#include "implicit.h"
      INTEGER :: I,J,K
      PARAMETER (D0=0.D00, KCOL=6)
      CHARACTER*(*) KEY
      DIMENSION SHESS(*), CSTRA(*), SCTRA(*)
#include "abainf.h"
#include "nuclei.h"
      INTEGER BEGIN, NROW, NCOL, NROW_DIM
      real(8), allocatable :: CHESS(:,:)
C
      IF (KEY .EQ. 'HESMOL') THEN
         NROW = 3*NUCDEP
         NROW_DIM = NROW
         ITYPE = 1
      ELSE IF (KEY .EQ. 'CENTERS') THEN
         NROW = 3*NUCDEP
         NROW_DIM = MXCOOR
         ITYPE = 1
      ELSE IF (KEY .EQ. 'ATOMS') THEN
         NROW = 3*NATOMS
         NROW_DIM = MXCOOR
         ITYPE = 1
      ELSE IF (KEY .EQ. 'SPNSPN') THEN
         NROW = 3*NUCDEP
         NROW_DIM = MXCOOR
         ITYPE = 2
      ELSE
         WRITE (LUPRI,'(//A/2A//)') ' ERROR in PRIHES ',
     &      '     - incorrect specification of keyword: ', KEY
         CALL QUIT('Incorrect keyword in PRIHES.')
      END IF
      IF (MAXREP .EQ. 0) THEN
         LAST = MIN(NROW,KCOL)
         BEGIN= 1
   50    CONTINUE
         WRITE (LUPRI,1000) (NAMEX(I),I = BEGIN,LAST)
         WRITE (LUPRI,'()')
         NCOL = 1
         DO 100 K = BEGIN,NROW
            DO 200 I = 1,NCOL
               IF (SHESS(K+NROW_DIM*(BEGIN-2+I)) .NE. D0) GO TO 400
  200       CONTINUE
            GO TO 300
  400       WRITE (LUPRI,2000) NAMEX(K),
     &            (SHESS(K+NROW_DIM*(BEGIN-2+J)),J=1,NCOL)
            IF (MOD(K,3) .EQ. 0) WRITE (LUPRI,'()')
  300       IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
  100    CONTINUE
         WRITE (LUPRI,'()')
         LAST = MIN(LAST+KCOL,NROW)
         BEGIN= BEGIN+NCOL
         IF (BEGIN.LE.NROW) GO TO 50
      ELSE
         NOFF = 0
         DO 500 IREP = 0, MAXREP
            IF (DOSYM(IREP+1) .AND. (NCRREP(IREP,ITYPE) .GT. 0)) THEN
               WRITE (LUPRI,'(/11X,A,I1/)') 'Symmetry ',IREP + 1
               NROW  = NOFF + NCRREP(IREP,ITYPE)
               LAST  = MIN(NROW,NOFF + KCOL)
               BEGIN = NOFF + 1
C
  550          CONTINUE
                  WRITE(LUPRI,1000)(NAMEX(IPTCOR(I,ITYPE)),I=BEGIN,LAST)
                  WRITE(LUPRI,'()')
                  NCOL = 1
                  DO 40 K = BEGIN, NROW
                     WRITE (LUPRI,2000) NAMEX(IPTCOR(K,ITYPE)),
     &                     (SHESS(K+NROW_DIM*(BEGIN-2+J)),J=1,NCOL)
                     IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40             CONTINUE
                  WRITE (LUPRI,'()')
                  LAST  = MIN(LAST + KCOL,NROW)
                  BEGIN = BEGIN + NCOL
               IF (BEGIN.LE.NROW) GO TO 550
            END IF
            NOFF = NOFF + NCRREP(IREP,ITYPE)
  500    CONTINUE
C
C        Print in non-symmetry basis
C
         WRITE (LUPRI,'(//)')
         NCOOR = 3*NUCDEP
         allocate( CHESS(NCOOR,NCOOR) )
         CALL DZERO(CHESS,NCOOR*NCOOR)
         CALL TRAHES(SHESS,NROW_DIM,CHESS,CSTRA,SCTRA,NCOOR,NCOOR,ITYPE)
         CALL PR2DER(CHESS,NCOOR,NCOOR,LUPRI)
         deallocate( CHESS )
      END IF
      WRITE (LUPRI,'()')
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT (8X,6(3X,A6,3X),(3X,A6,3X))
 2000 FORMAT (1X,A6,6F12.6)
      END
C  /* Deck tragrd */
      SUBROUTINE TRAGRD(SGRAD,CGRAD,CSTRA,SCTRA,NCR0,NCOOR)
C
C     Transform molecular gradient from symmetry nuclear coordinates
C     to cartesian nuclear coordinates
C
      use mxcent_mod, only: MXCOOR
#include "implicit.h"
      DIMENSION SGRAD(MXCOOR), CGRAD(MXCOOR)
      DIMENSION CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR)
      CALL TRACOR(CSTRA,SCTRA,1,NCOOR,0)
      CALL DGEMM('N','N',NCOOR,1,NCR0,1.D0,
     &           SCTRA(1,1),NCOOR,
     &           SGRAD,MXCOOR,0.D0,
     &           CGRAD,MXCOOR)
      RETURN
      END
C  /* Deck trahes */
      SUBROUTINE TRAHES(SHESS,NDIMS,CHESS,CSTRA,SCTRA,NDIMC,NCOOR,ITYPE)
C
C     Transform molecular Hessian from symmetry nuclear coordinates
C     to cartesian nuclear coordinates
C
#include "implicit.h"
      INTEGER :: I,J
      DIMENSION SHESS(NDIMS,NDIMS), CHESS(NDIMC,NDIMC),
     &          CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR)
      real(8), allocatable :: AMAT(:,:)
C
      CALL TRACOR(CSTRA,SCTRA,ITYPE,NCOOR,0)
      DO 100 I = 1, NCOOR
         DO 200 J = 1, I - 1
            SHESS(J,I) = SHESS(I,J)
  200    CONTINUE
  100 CONTINUE
      allocate( AMAT(NCOOR,NCOOR) )
      CALL DGEMM('N','N',NCOOR,NCOOR,NCOOR,1.D0,
     &           SCTRA(1,1),NCOOR,
     &           SHESS,NDIMS,0.D0,
     &           AMAT,NCOOR)
      CALL DGEMM('N','T',NCOOR,NCOOR,NCOOR,1.D0,
     &           AMAT,NCOOR,
     &           SCTRA(1,1),NCOOR,0.D0,
     &           CHESS,NDIMC)
      deallocate( AMAT )
      RETURN
      END

C  /* Deck prigeolu */
      SUBROUTINE PRIGEOLU(LU,COORD)
C
C     OC cloned PRIGEO (april 2009) of Jun 29 1988 tuh
C
!      use nuclei_mod, only: ISTBNU, NAMDEP, NUCIND
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use symmet_mod, only: ISYMAX, MAXOPR, PT
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
#include "implicit.h"
      INTEGER :: I
C
      DIMENSION COORD(3,*)
C
#include "nuclei.h"

C
      CHARACTER*8 ATNAME 
      LOGICAL FNONB
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 200 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               ATNAME='        '
               FNONB = .FALSE.
               INOTBLANK=0
               DO I=1,6
                  IF (NAMDEP(IATOM)(I:I).NE. ' ') THEN 
                     INOTBLANK=INOTBLANK+1
                     ATNAME(INOTBLANK:INOTBLANK)=
     *                  NAMDEP(IATOM)(I:I)
                     FNONB=.TRUE.
                  ENDIF 
                  IF (FNONB .AND.(NAMDEP(IATOM)(I:I).EQ. ' ')) THEN 
                    GOTO 17 
                  ENDIF 
               ENDDO 
  17           CONTINUE 
               WRITE (LU,'(A,3(2X,F22.16))') ATNAME, 
     *           (PT(IAND(ISYMAX(I,1),ISYMOP))*COORD(I,ICENT),I=1,3)
            END IF
  200    CONTINUE
  100 CONTINUE
      IATOM = 0
      WRITE(LU,*) " Atomic labels "
      DO 300 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 400 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               WRITE (LU,'(A)') NAMDEP(IATOM)
            END IF
  400    CONTINUE
  300 CONTINUE

      RETURN
      END
C  /* Deck prigeo */
      SUBROUTINE PRIGEO(COORD)
C
C     Jun 29 1988 tuh
C
      use symmet_mod, only: ISYMAX, MAXOPR, PT
      use maxaqn_mod, only: MXAQN, MXQN
!      use nuclei_mod, only: ISTBNU, NUCIND
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
#include "implicit.h"
      INTEGER :: I
C
      DIMENSION COORD(3,*)
C
#include "nuclei.h"
C
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 200 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               WRITE (LUPRI,'(1X,A,F17.10,2F24.10)') NAMDEP(IATOM),
     *           (PT(IAND(ISYMAX(I,1),ISYMOP))*COORD(I,ICENT),I=1,3)
            END IF
  200    CONTINUE
  100 CONTINUE
      WRITE (LUPRI,'(//)')
      RETURN
      END
C  /* Deck symhes */
      SUBROUTINE SYMHES(HESS,NROW_HESS)
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use nuclei_mod, only: NUCDEP
      use kinds_mod, only: dp
#include "implicit.h"
      DIMENSION HESS(NROW_HESS,NROW_HESS)
      REAL(dp), PARAMETER :: DP5 = 0.5D0
      INTEGER :: I,J
#include "nuclei.h"
C
      NCOORD = 3*NUCDEP
      DO 100 I = 1, NCOORD
         DO 200 J = 1, I
            AVERAG = DP5*(HESS(I,J) + HESS(J,I))
            HESS(I,J) = AVERAG
            HESS(J,I) = AVERAG
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck molchr */
      SUBROUTINE MOLCHR(ICHRG)
C
C     Calculates total charge of molecule
C
      use maxaqn_mod, only: MXAQN, MXQN
!      use inforb_mod, only: NISHT
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use infinp_mod, only: NACTEL
!      use nuclei_mod, only: NUCIND
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: FMULT
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D00
C
#include "nuclei.h"
#include "infinp.h"
#include "inforb.h"
C
C     Nuclear charge
C
      CHRGNU = D0
      DO 100 IATOM = 1, NUCIND
         CHRGNU = CHRGNU + FMULT(ISTBNU(IATOM))*CHARGE(IATOM)
  100 CONTINUE
C
C     Subtract number of electrons
C
      ICHRG = NINT(CHRGNU) - (2*NISHT + NACTEL)
      RETURN
      END
C  /* Deck rmolchr */
      SUBROUTINE RMOLCHR(ICHRG)
C
C     Calculates sum of all nuclear charges
C
      use maxaqn_mod, only: MXAQN, MXQN
!      use nuclei_mod, only: NUCIND
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: FMULT
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D00
C
#include "nuclei.h"
C
C     Nuclear charge
C
      CHRGNU = D0
      DO 100 IATOM = 1, NUCIND
         CHRGNU = CHRGNU + FMULT(ISTBNU(IATOM))*CHARGE(IATOM)
  100 CONTINUE
      ICHRG = NINT(CHRGNU)
      RETURN
      END
C  /* Deck disotp */
      FUNCTION DISOTP(IATOM,ISOTOP,TYPE)
C
C     NOTE: Isotopes are sorted according to abundance,
C     i.e. DISOTP(IATOM,1,TYPE) will return the most abundant
C     isotope etc.
C
C      DATNUC(1,ISOTOP,IATOM) - atomic mass
C      DATNUC(2,ISOTOP,IATOM) - abundance
C      DATNUC(3,ISOTOP,IATOM) - spin
C      DATNUC(4,ISOTOP,IATOM) - magnetic moment
C      DATNUC(5,ISOTOP,IATOM) - quadrupole moment
C
C     Extended to Z=54 on march 10, 1994, K.Ruud
C     Extended to Z=86 feb. 1996, S. Kirpekar
C     Extended to Z=109 1997/06/20 Joern Thyssen
C     Corrected Z=55-109 for correct units 1999/04/23 W.A. de Jong
C     Extended to Z=118 1997/05/28 Trond Saue
C     Extended to Z=138 2023/01/25  Kirk Peterson
C       using data from Z.Zhou et al., Atomic Data Nuc. Data Tables 114, 262 (2017)
C     Corrected elements 108 and 109, data from 
C       Handbook of Chemistry and Physics, 83rd Edition, by Thierfelder and R. Bast
C     Updated magnetic and quadrupole moments,  2003/03/06 Joost van Stralen
C     Updated 209Po magnetic moment (2016/02/03 Andre Gomes)
C       NJ Stone, TABLE OF NUCLEAR MAGNETIC DIPOLE AND ELECTRIC QUADRUPOLE MOMENTS,
C       IAEA  International Nuclear Data Committee, INDC(NDS)-0658, February 2014
C
C     Proton mass and electron charge:
C        1986 CODATA Recommended Values
C
C     Nuclear masses:
C        A. H. Wapstra and K. Bos, Atomic Data and Nuclear Tables 19 (1977) 177
C
C     Abundancies:
C        Handbook of Chemistry and Physics, 73rd Edition
C        
C     Nuclear moments and spins:
C        P. Raghavan, Atomic Data and Nuclear Data Tables 42 (1989) 189
C
C        2014-03-03 Radovan Bast: updated Ag value as suggested by Michal Jaszunski in private communication
C                                 corrected values for Nd (lines interchanged) and Tl (typo)
C                                 as suggested by Michal Jaszunski in private communication
C
C     Quadrupole moments:
C        P.Pykkoe
C        Molecular Physics 99 (2001) 1617
C 
C              or the same
C
C        R. K. Harris et al
C        Pure Applied Chemistry 77 (2001) 1795
C        (IUPAC Recommendations 2001)
C
C     Nuclear masses, Abundancies, nuclear moments, spins 
C     and quadrupole moments for Z= 55 to Z = 86:
C       I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu
C       Quantities, Units and Symbols in Physical Chemistry
C       (IUPAC, Blackwell Scientific, Oxford, 1988)
C
C 07 Febr 2016 Miro Ilias: added nuclear mass for the Tc atom to enable automatic symmetry detection
C (otherwise was getting NaNs), based on
C https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
C   G. Audi, A. H. Wapstra Nucl. Phys A. 1993, 565, 1-65 
C   G. Audi, A. H. Wapstra Nucl. Phys A. 1995, 595, 409-480.
C
C 14 Sep 2020 - Agustin Aucar added data for 223Ra and 225Ra.
C     Masses and nuclear spins from
C       Haynes, Handbook of Chemistry and Physics, 97th Edition, 2017
C     Magnetic moments from
C       NJ Stone, TABLE OF NUCLEAR MAGNETIC DIPOLE AND ELECTRIC QUADRUPOLE MOMENTS,
C       IAEA  International Nuclear Data Committee, INDC(NDS)-0658, February 2014
C
      use priunit_mod, only: LUPRI
!      use codata_mod, only: ECHARGE, PMASS, UMASS, XFAMU
!      use pi_mod, only: PI
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,J,K
#include "codata.h"
#include "pi.h"
      INTEGER, PARAMETER :: MAXISO = 6, MAXCHR = 138
      REAL(dp), PARAMETER :: D0 = 0.0D0, D4= 4.0D0
      REAL(dp), PARAMETER :: THRESH = 1.0D-10
C
      CHARACTER*(*) TYPE
      DIMENSION DATNUC(5,MAXISO,0:MAXCHR)
C
C     H - Ne
C     ======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=0,10) /
C
C
C     Dummy:
C
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
     &  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,  0.000000D0,
C
C     H:
C
     &  1.007825D0, 99.985000D0,   .500000D0,  2.792847D0,   .000000D0,
     &  2.014102D0,   .015000D0,  1.000000D0,   .857438D0,   .002860D0,
     &  3.016049D0,   .000000D0,   .500000D0,  2.978962D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     He:
C
     &  4.002603D0, 99.999870D0,   .000000D0,   .000000D0,   .000000D0,
     &  3.016029D0,   .000130D0,   .500000D0, -2.127625D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Li:
C
     &  7.016005D0, 92.500000D0,  1.500000D0,  3.256427D0,  -.040100D0,
     &  6.015123D0,  7.500000D0,  1.000000D0,   .822047D0,  -.000808D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Be:
C
     &  9.012183D0,100.000000D0,  1.500000D0, -1.177800D0,   .052880D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     B:
C
     & 11.009305D0, 80.100000D0,  1.500000D0,  2.688649D0,   .040590D0,
     & 10.012938D0, 19.900000D0,  3.000000D0,  1.800645D0,   .084590D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     C:
C
     & 12.000000D0, 98.900000D0,   .000000D0,   .000000D0,   .000000D0,
     & 13.003355D0,  1.100000D0,   .500000D0,   .702412D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     N:
C
     & 14.003074D0, 99.630000D0,  1.000000D0,   .403761D0,   .020440D0,
     & 15.000109D0,   .370000D0,   .500000D0,  -.283189D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     O:
C
     & 15.994915D0, 99.760000D0,   .000000D0,   .000000D0,   .000000D0,
     & 17.999159D0,   .200000D0,   .000000D0,   .000000D0,   .000000D0,
     & 16.999131D0,   .040000D0,  2.500000D0, -1.893790D0,  -.025580D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     F:
C
     & 18.998403D0,100.000000D0,   .500000D0,  2.628868D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ne:
C
     & 19.992439D0, 90.480000D0,   .000000D0,   .000000D0,   .000000D0,
     & 21.991384D0,  9.250000D0,   .000000D0,   .000000D0,   .000000D0,
     & 20.993845D0,   .270000D0,  1.500000D0,  -.661797D0,   .101550D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
C
C    Na - Ar
C    =======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=11,18) /
C
C     Na:
C
     & 22.989770D0,100.000000D0,  1.500000D0,  2.217656D0,   .104000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Mg:
C
     & 23.985045D0, 78.990000D0,   .000000D0,   .000000D0,   .000000D0,
     & 25.982595D0, 11.010000D0,   .000000D0,   .000000D0,   .000000D0,
     & 24.985839D0, 10.000000D0,  2.500000D0,  -.855450D0,   .199400D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Al:
C
     & 26.981541D0,100.000000D0,  2.500000D0,  3.641507D0,   .146600D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Si:
C
     & 27.976928D0, 92.230000D0,   .000000D0,   .000000D0,   .000000D0,
     & 28.976496D0,  4.670000D0,   .500000D0,  -.555290D0,   .000000D0,
     & 29.973772D0,  3.100000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     P:
C
     & 30.973763D0,100.000000D0,   .500000D0,  1.131600D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     S:
C
     & 31.972072D0, 95.020000D0,   .000000D0,   .000000D0,   .000000D0,
     & 33.967868D0,  4.210000D0,   .000000D0,   .000000D0,   .000000D0,
     & 32.971459D0,   .750000D0,  1.500000D0,   .643821D0,  -.067800D0,
     & 35.967079D0,   .020000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Cl:
C
     & 34.968853D0, 75.770000D0,  1.500000D0,   .821874D0,  -.081650D0,
     & 36.965903D0, 24.230000D0,  1.500000D0,   .684124D0,  -.064350D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ar:
C
     & 39.962383D0, 99.600000D0,   .000000D0,   .000000D0,   .000000D0,
     & 35.967546D0,   .337000D0,   .000000D0,   .000000D0,   .000000D0,
     & 37.962732D0,   .063000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
C
C     K - Ca
C     ======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=19,20) /
C
C     K:
C
     & 38.963708D0, 93.258100D0,  1.500000D0,   .391507D0,   .058500D0,
     & 40.961825D0,  6.730200D0,  1.500000D0,   .214893D0,   .071100D0,
     & 39.963999D0,   .011700D0,  4.000000D0, -1.298100D0,  -.073000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ca:
C
     & 39.962591D0, 96.941000D0,   .000000D0,   .000000D0,   .000000D0,
     & 43.955485D0,  2.086000D0,   .000000D0,   .000000D0,   .000000D0,
     & 41.958622D0,   .647000D0,   .000000D0,   .000000D0,   .000000D0,
     & 47.952532D0,   .187000D0,   .000000D0,   .000000D0,   .000000D0,
     & 42.958770D0,   .135000D0,  3.500000D0, -1.317643D0,  -.040800D0,
     & 45.953689D0,   .004000D0,   .000000D0,   .000000D0,   .000000D0/
C
C     Sc - Zn
C     =======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=21,30) /
C
C     Sc:
C
     & 44.955914D0,100.000000D0,  3.500000D0,  4.756487D0,  -.220000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ti:
C
     & 47.947947D0, 73.800000D0,   .000000D0,   .000000D0,   .000000D0,
     & 45.952633D0,  8.000000D0,   .000000D0,   .000000D0,   .000000D0,
     & 46.951765D0,  7.300000D0,  2.500000D0,  -.788480D0,   .302000D0,
     & 48.947871D0,  5.500000D0,  3.500000D0, -1.104170D0,   .247000D0,
     & 49.944786D0,  5.400000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     V:
C
     & 50.943963D0, 99.750000D0,  3.500000D0,  5.148706D0,  -.052000D0,
     & 49.947161D0,   .250000D0,  6.000000D0,  3.345689D0,   .210000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Cr:
C
     & 51.940510D0, 83.790000D0,   .000000D0,   .000000D0,   .000000D0,
     & 52.940651D0,  9.500000D0,  1.500000D0,  -.474540D0,  -.150000D0,
     & 49.946046D0,  4.345000D0,   .000000D0,   .000000D0,   .000000D0,
     & 53.938882D0,  2.365000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Mn:
C
     & 54.938046D0,100.000000D0,  2.500000D0,  3.468719D0,   .330000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Fe:
C
     & 55.934939D0, 91.720000D0,   .000000D0,   .000000D0,   .000000D0,
     & 53.939612D0,  5.900000D0,   .000000D0,   .000000D0,   .000000D0,
     & 56.935396D0,  2.100000D0,   .500000D0,   .090623D0,   .000000D0,
     & 57.933278D0,   .280000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Co:
C
     & 58.933198D0,100.000000D0,  3.500000D0,  4.627000D0,   .420000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ni:
C
     & 57.935347D0, 68.077000D0,   .000000D0,   .000000D0,   .000000D0,
     & 59.930789D0, 26.223000D0,   .000000D0,   .000000D0,   .000000D0,
     & 61.928346D0,  3.634000D0,   .000000D0,   .000000D0,   .000000D0,
     & 60.931059D0,  1.140000D0,  1.500000D0,  -.750020D0,   .162000D0,
     & 63.927968D0,  0.926000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Cu:
C
     & 62.929599D0, 69.170000D0,  1.500000D0,  2.227206D0,  -.220000D0,
     & 64.927792D0, 30.830000D0,  1.500000D0,  2.381610D0,  -.204000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Zn:
C
     & 63.929145D0, 48.600000D0,   .000000D0,   .000000D0,   .000000D0,
     & 65.926035D0, 27.900000D0,   .000000D0,   .000000D0,   .000000D0,
     & 67.924846D0, 18.800000D0,   .000000D0,   .000000D0,   .000000D0,
     & 66.927129D0,  4.100000D0,  2.500000D0,   .875479D0,   .150000D0,
     & 69.925325D0,   .600000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
C
C     Ga - Kr
C     =======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=31,36) /
C
C     Ga:
C
     & 68.925581D0, 60.108000D0,  1.500000D0,  2.016589D0,   .171000D0,
     & 70.924701D0, 39.892000D0,  1.500000D0,  2.562266D0,   .107000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ge:
C
     & 73.921179D0, 35.940000D0,   .000000D0,   .000000D0,   .000000D0,
     & 71.922080D0, 27.660000D0,   .000000D0,   .000000D0,   .000000D0,
     & 69.924250D0, 21.240000D0,   .000000D0,   .000000D0,   .000000D0,
     & 72.923464D0,  7.720000D0,  4.500000D0,  -.879468D0,  -.196000D0,
     & 75.921403D0,  7.440000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     As:
C
     & 74.921596D0,100.000000D0,  1.500000D0,  1.439475D0,   .314000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Se:
C
     & 79.916521D0, 49.610000D0,   .000000D0,   .000000D0,   .000000D0,
     & 77.917304D0, 23.770000D0,   .000000D0,   .000000D0,   .000000D0,
     & 75.919207D0,  9.360000D0,   .000000D0,   .000000D0,   .000000D0,
     & 81.916709D0,  8.740000D0,   .000000D0,   .000000D0,   .000000D0,
     & 76.919908D0,  7.630000D0,   .500000D0,   .535042D0,   .000000D0,
     & 73.922477D0,  0.890000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Br:
C
     & 78.918336D0, 50.690000D0,  1.500000D0,  2.106400D0,   .313000D0,
     & 80.916290D0, 49.310000D0,  1.500000D0,  2.270562D0,   .261500D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Kr:
C
     & 83.911506D0, 57.000000D0,   .000000D0,   .000000D0,   .000000D0,
     & 85.910614D0, 17.300000D0,   .000000D0,   .000000D0,   .000000D0,
     & 81.913483D0, 11.600000D0,   .000000D0,   .000000D0,   .000000D0,
     & 82.914134D0, 11.500000D0,  4.500000D0,  -.970669D0,   .259000D0,
     & 79.916375D0,  2.250000D0,   .000000D0,   .000000D0,   .000000D0,
     & 77.920397D0,  0.350000D0,   .000000D0,   .000000D0,   .000000D0/
C
C     Rb:
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=37,45) /
C
     & 84.911800D0, 72.170000D0,  2.500000D0,  1.353352D0,   .276000D0,
     & 86.909184D0, 27.830000D0,  1.500000D0,  2.751818D0,   .133500D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Sr:
C
     & 87.905625D0, 82.580000D0,   .000000D0,   .000000D0,   .000000D0,
     & 85.909273D0,  9.860000D0,   .000000D0,   .000000D0,   .000000D0,
     & 86.908890D0,  7.000000D0,  4.500000D0, -1.093603D0,   .335000D0,
     & 83.913428D0,   .560000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Y:
C
     & 88.905856D0,100.000000D0,   .500000D0,  -.137415D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Zr:
C
     & 89.904708D0, 51.450000D0,   .000000D0,   .000000D0,   .000000D0,
     & 93.906319D0, 17.380000D0,   .000000D0,   .000000D0,   .000000D0,
     & 91.905039D0, 17.150000D0,   .000000D0,   .000000D0,   .000000D0,
     & 90.905644D0, 11.220000D0,  2.500000D0, -1.303620D0,  -.176000D0,
     & 95.908272D0,  2.800000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Nb:
C
     & 92.906378D0,100.000000D0,  4.500000D0,  6.170500D0,  -.320000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Mo:
C
     & 97.905405D0, 24.130000D0,   .000000D0,   .000000D0,   .000000D0,
     & 95.904676D0, 16.680000D0,   .000000D0,   .000000D0,   .000000D0,
     & 94.905838D0, 15.920000D0,  2.500000D0,  -.914200D0,  -.022000D0,
     & 93.905086D0, 14.840000D0,   .000000D0,   .000000D0,   .000000D0,
     & 99.907473D0,  9.630000D0,   .000000D0,   .000000D0,   .000000D0,
     & 96.906018D0,  9.550000D0,  2.500000D0,  -.933500D0,  0.255000D0,
C
C     Tc:
C
     & 97.907216D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ru:
C
     &101.904348D0, 31.600000D0,   .000000D0,   .000000D0,   .000000D0,
     &103.905422D0, 18.600000D0,   .000000D0,   .000000D0,   .000000D0,
     &100.905581D0, 17.100000D0,  2.500000D0,  -.718800D0,   .457000D0,
     & 98.905937D0, 12.700000D0,  2.500000D0,  -.641300D0,   .079000D0,
     & 99.904218D0, 12.600000D0,   .000000D0,   .000000D0,   .000000D0,
     & 95.907596D0,  5.540000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Rh:
C
     &102.905503D0,100.000000D0,   .500000D0,  -.088400D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0/
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=46,54) /
C
C     Pd:
C
     &105.903475D0, 27.330000D0,   .000000D0,   .000000D0,   .000000D0,
     &107.903894D0, 26.460000D0,   .000000D0,   .000000D0,   .000000D0,
     &104.905075D0, 22.330000D0,  2.500000D0,  -.642000D0,   .660000D0,
     &109.905169D0, 11.720000D0,   .000000D0,   .000000D0,   .000000D0,
     &103.904026D0, 11.140000D0,   .000000D0,   .000000D0,   .000000D0,
     &101.905609D0,  1.020000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Ag:
C
     &106.905095D0, 51.839000D0,   .500000D0,  -.113570D0,   .000000D0,
     &108.904754D0, 48.161000D0,   .500000D0,  -.130691D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Cd:
C
     &113.903361D0, 28.730000D0,   .000000D0,   .000000D0,   .000000D0,
     &111.902761D0, 24.130000D0,   .000000D0,   .000000D0,   .000000D0,
     &110.904182D0, 12.800000D0,   .500000D0,  -.594886D0,  -.850000D0,
     &109.903007D0, 12.490000D0,   .000000D0,   .000000D0,   .000000D0,
     &112.904401D0, 12.220000D0,   .500000D0,  -.622301D0,   .000000D0,
     &115.904758D0,  7.490000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     In:
C
     &114.903875D0, 95.700000D0,  4.500000D0,  5.540800D0,   .810000D0,
     &112.904056D0,  4.300000D0,  4.500000D0,  5.528900D0,   .799000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Sn:
C
     &119.902199D0, 32.590000D0,   .000000D0,   .000000D0,   .000000D0,
     &117.901607D0, 24.220000D0,   .000000D0,   .000000D0,   .000000D0,
     &115.901744D0, 14.530000D0,   .000000D0,   .000000D0,   .000000D0,
     &118.903310D0,  8.580000D0,   .500000D0, -1.047280D0,   .000000D0,
     &116.902954D0,  7.680000D0,   .500000D0, -1.001040D0,   .000000D0,
     &123.905271D0,  5.790000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Sb:
C
     &120.903824D0, 57.360000D0,  2.500000D0,  3.363400D0,  -.360000D0,
     &122.904222D0, 42.640000D0,  3.500000D0,  2.549800D0,  -.490000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Te:
C
     &129.906229D0, 33.870000D0,   .000000D0,   .000000D0,   .000000D0,
     &127.904464D0, 31.700000D0,   .000000D0,   .000000D0,   .000000D0,
     &125.903310D0, 18.930000D0,   .000000D0,   .000000D0,   .000000D0,
     &124.904435D0,  7.120000D0,   .500000D0,  -.888505D0,   .000000D0,
     &123.902825D0,  4.790000D0,   .000000D0,   .000000D0,   .000000D0,
     &121.903055D0,  2.590000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     I:
C
     &126.904477D0,100.000000D0,  2.500000D0,  2.813273D0,  -.710000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
     &   .000000D0,   .000000D0,   .000000D0,   .000000D0,   .000000D0,
C
C     Xe:
C
     &131.904148D0, 26.900000D0,   .000000D0,   .000000D0,   .000000D0,
     &128.904780D0, 26.400000D0,   .500000D0,  -.777976D0,   .000000D0,
     &130.905076D0, 21.200000D0,  1.500000D0,   .691862D0,  -.114000D0,
     &133.905395D0, 10.400000D0,   .000000D0,   .000000D0,   .000000D0,
     &135.907219D0,  8.900000D0,   .000000D0,   .000000D0,   .000000D0,
     &129.903510D0,  4.100000D0,   .000000D0,   .000000D0,   .000000D0/
C
C
C
C
C    Cs - Rn*
C    =======
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=55,64) /
C
C     Cs:
C
     &132.905429D0,100.000000D0,  3.50000D0,   2.582025D0,  -0.00343D0,
     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
     &  0.00000D0,   0.000000D0,  0.00000D0,   0.000000D0,   0.00000D0,
C
C     Ba:
C   
     &137.905232D0, 71.70000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &136.905812D0, 11.23000D0,   1.50000D0,  0.937365D0,   0.245000D0,
     &135.904553D0,  7.85400D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &134.905665D0,  6.59200D0,   1.50000D0,  0.837943D0,   0.160000D0,
     &133.904486D0,  2.41700D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &131.905042D0,  0.10100D0,   0.00000D0,  0.000000D0,   0.000000D0,
C
C     La:
C
     &138.906347D0, 99.90980D0,   3.50000D0,  2.7830455D0,  0.200000D0,
     &137.907105D0,  0.09020D0,   5.00000D0,  3.7136460D0,  0.450000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
C
C     Ce
C
     &139.905433D0, 88.48000D0,   0.00000D0,  0.000000D0,    0.00000D0,
     &141.909241D0, 11.08000D0,   0.00000D0,  0.000000D0,    0.00000D0,
     &137.905985D0,  0.25000D0,   0.00000D0,  0.000000D0,    0.00000D0,
     &135.907140D0,  0.19000D0,   0.00000D0,  0.000000D0,    0.00000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,   0.00000D0,  0.000000D0,   0.000000D0,
C
C     Pr:
C
     &140.907647D0,100.0000D0,    2.50000D0,  4.275400D0,  -0.058900D0,
     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     & 0.000000D0,  0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
C
C     Nd:
C
     &141.907719D0, 27.130000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &143.910083D0, 23.800000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &145.913113D0, 17.190000D0,   0.00000D0,  0.000000D0,   0.000000D0,
     &142.909810D0, 12.180000D0,   3.50000D0, -1.065000D0,  -0.630000D0,
     &144.912570D0,  8.300000D0,   3.50000D0, -0.650000D0,  -0.330000D0,
     &147.916889D0,  5.760000D0,   0.00000D0,  0.000000D0,   0.000000D0,
C
C     Pm:
C
     &144.912743D0,100.000000D0,   2.50000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
C
C     Sm:
C
     &151.919728D0, 26.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &153.922205D0, 22.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &146.914894D0, 15.000000D0,   3.500000D0,-0.814800D0,  -0.259000D0,
     &148.917180D0, 13.800000D0,   3.500000D0,-0.671700D0,   0.075000D0,
     &147.914819D0, 11.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &149.917273D0,  7.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Eu:
C
     &152.921225D0, 52.200000D0,   2.500000D0, 1.533000D0,   2.412000D0,
     &150.919702D0, 47.800000D0,   2.500000D0, 3.471700D0,   0.903000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
     &  0.00000D0,   0.00000D0,    0.00000D0,  0.000000D0,   0.000000D0,
C
C     Gd:
C
     &157.924019D0, 24.840000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &159.927049D0, 21.860000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &155.922118D0, 20.470000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &156.923956D0, 15.650000D0,   1.500000D0,-0.337260D0,   1.350000D0,
     &154.922618D0, 14.800000D0,   1.500000D0,-0.257230D0,   1.270000D0,
     &153.920861D0,  2.180000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C     Tb:
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=65,74) /
C
     &158.925342D0,100.000000D0,   1.500000D0, 2.014000D0,   1.432000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Dy:
C
     &163.929171D0, 28.200000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &161.926795D0, 25.500000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &162.928728D0, 24.900000D0,   2.500000D0, 0.672600D0,   2.648000D0,
     &160.926930D0, 18.900000D0,   2.500000D0,-0.480300D0,   2.507000D0,
     &159.925193D0,  2.340000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &157.924277D0,  0.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ho:
C
     &164.930319D0,100.000000D0,   3.500000D0, 4.173000D0,   3.580000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Er:
C
     &165.930290D0, 33.600000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &167.932368D0, 26.800000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &166.932368D0, 22.950000D0,   3.500000D0,-0.563850D0,   3.565000D0,
     &169.935461D0, 14.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &163.929198D0,  1.610000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &161.928775D0,  0.140000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Tm:
C
     &168.934212D0,100.000000D0,   0.500000D0,-0.231600D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Yb:
C
     &173.938859D0, 31.800000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &171.936378D0, 21.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &172.938208D0, 16.120000D0,   2.500000D0,-0.679890D0,   2.800000D0,
     &170.936323D0, 14.300000D0,   0.500000D0, 0.493670D0,   0.000000D0,
     &175.942564D0, 12.700000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &169.934759D0,  3.050000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Lu:
C
     &174.940770D0, 97.410000D0,   3.500000D0, 2.232700D0,   3.490000D0,
     &175.942679D0,  2.590000D0,   7.000000D0, 3.169200D0,   4.970000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Hf:
C
     &179.9465457D0,35.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &177.943696D0, 27.297000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &176.943217D0, 18.606000D0,   3.500000D0, 0.793500D0,   3.365000D0,
     &178.9458122D0,13.629000D0,   4.500000D0,-0.640900D0,   3.793000D0,
     &175.941406D0,  5.206000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &173.940044D0,  0.162000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ta:
C
     &180.947992D0, 99.988000D0,   3.500000D0, 2.370500D0,   3.170000D0,
     &179.947462D0,  0.012000D0,   8.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     W:
C
     &183.950928D0, 30.670000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &185.954357D0, 28.600000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &181.948202D0, 26.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &182.950928D0, 14.300000D0,   0.500000D0, 0.11778476,   0.000000D0,
     &179.947462D0,  0.162000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C     Re:
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=75,80) /
C
     &186.955744D0, 62.600000D0,   2.500000D0, 3.219700D0,   2.070000D0,
     &184.952951D0, 37.400000D0,   2.500000D0, 3.187100D0,   2.180000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Os:
C
     &191.961467D0, 41.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &189.958436D0, 26.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &188.958436D0, 16.100000D0,   1.500000D0, 0.659933D0,   0.856000D0,
     &187.955830D0, 13.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &186.955741D0,  1.600000D0,   0.500000D0, 0.06465189,   0.000000D0,
     &185.953830D0,  1.580000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ir:
C
     &192.962917D0, 62.600000D0,   1.500000D0, 0.163700D0,   0.751000D0,
     &190.960584D0, 37.400000D0,   1.500000D0, 0.150700D0,   0.816000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Pt:
C
     &194.964766D0, 33.800000D0,   0.500000D0, 0.609520D0,   0.000000D0,
     &193.962655D0, 32.900000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &195.964926D0, 25.300000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &197.967869D0,  7.200000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &191.961019D0,  0.790000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &189.959917D0,  0.010000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Au:
C   
     &196.966543D0,100.000000D0,   1.500000D0, 0.148158D0,   0.547000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Hg:
C
     &201.970617D0, 29.860000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &199.968300D0, 23.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &198.968254D0, 16.870000D0,   0.500000D0, 0.50588549D0, 0.850000D0,
     &200.970277D0, 13.180000D0,   1.500000D0,-0.5602257D0,  0.385000D0,
     &197.966743D0,  9.970000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &203.973467D0,  6.870000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C     Tl:
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=81,86) /
C
     &204.974401D0, 70.476000D0,   0.500000D0, 1.63821461D0, 0.000000D0,
     &202.972320D0, 29.524000D0,   0.500000D0, 1.62225787D0, 0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Pb:
C   
     &207.976627D0, 52.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &205.975872D0, 24.100000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &206.975872D0, 22.100000D0,   0.500000D0, 0.592583D0,   0.000000D0,
     &203.973020D0,  1.400000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Bi:
C
     &208.980374D0,100.000000D0,   4.500000D0, 4.110600D0,  -0.516000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Po:
C
     &208.982404D0,100.000000D0,   0.500000D0, 0.688000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     At:
C
     &209.987126D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Rn:
C
     &222.017571D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C
C     Fr - Am
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=87,95) /
C
C     Fr:
C
     &223.019733D0,100.000000D0,   1.500000D0, 1.170000D0,   1.170000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ra:
C
     &226.025403D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &223.018500D0,  0.000000D0,   1.500000D0, 0.271000D0,   0.000000D0,
     &225.023610D0,  0.000000D0,   0.500000D0,-0.733800D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ac:
C
     &227.027750D0,100.000000D0,   1.500000D0, 1.100000D0,   1.700000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Th:
C
     &232.0380508D0,100.00000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Pa:
C
     &231.035880D0,100.000000D0,   1.500000D0, 2.010000D0,  -1.720000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     U:
C
     &238.0507847D0,99.274500D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &235.0439242D0, 0.720000D0,   3.500000D0,-0.380000D0,   4.936000D0,
     &234.0409468D0, 0.005500D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &233.039628D0,  0.000000D0,   2.500000D0, 0.590000D0,   3.663000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Np:
C
     &237.0481678D0,100.00000D0,   2.500000D0, 3.140000D0,   3.886000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Pu:
C
     &244.064199D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Am:
C
     &243.061373D0,100.000000D0,   2.500000D0, 1.610000D0,   4.210000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C
C     Cm - Mt(109)
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=96,109) /
C
C     Cm:
C
     &247.070347D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Bk:
C
     &247.070300D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Cf:
C
     &251.079580D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Es:
C
     &252.082944D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Fm:
C
     &257.095099D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Md:
C
     &258.098570D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     No:
C
     &259.100931D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Lr(103):
C
     &260.105320D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Rf(104):
C
     &261.108690D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Db(105):
C
     &262.113760D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Sg(106):
C
     &263.118220D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Bh(107):
C
     &262.122930D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Hs(108):
C
     &269.134100D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &263.128700D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &264.128400D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &265.130000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &267.137100D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Mt(109):
C
     &267.138000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &266.137900D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &268.138800D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C     110 - 118
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=110,118) /
C
C     Ds(110):
C
     &268.143500D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Rg(111):
C
     &272.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Cn(112):
C
     &285.177120D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Nh(113):
C     -- no data --- 
     &283.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Fl(114):
C
     &289.190420D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Mc(115):
C     -- no data --
     &294.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Lv(116):
C                   
     &293.204490D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ts(117):
C     -- no data --
     &306.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Og(118):
C     -- no data --
     &310.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=119,138) /
C
C     Uue(119):
C     -- no data --
     &314.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubn(120):
C     -- no data --
     &318.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubu(121):
C     -- no data --
     &330.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubb(122):
C     -- no data --
     &334.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubt(123):
C     -- no data --
     &337.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubq(124):
C     -- no data --
     &340.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubp(125):
C     -- no data --
     &344.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubh(126):
C     -- no data --
     &347.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubs(127):
C     -- no data --
     &350.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Ubo(128):
C     -- no data --
     &354.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C    Ube (129):
C     -- no data --
     &357.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utn(130):
C     -- no data --
     &361.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utu(131):
C     -- no data --
     &364.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utb(132):
C     -- no data --
     &367.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utt(133):
C     -- no data --
     &371.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utq(134):
C     -- no data --
     &374.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Utp(135):
C     -- no data --
     &378.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Uth(136):
C     -- no data --
     &381.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Uts(137):
C     -- no data --
     &385.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
C
C     Uto(138):
C     -- no data --
     &388.000000D0,100.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0,
     &  0.000000D0,  0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0/
C
C
      IF (ISOTOP .GT. MAXISO) THEN
         WRITE (LUPRI,'(//,A,2(/,A,I5),A)')
     &         ' ISOTOP too large in DISOTP.',
     &         ' Input value:  ',ISOTOP,
     &         ' Maximum value:',MAXISO,
     &         ' Program cannot continue.'
         CALL QUIT('MAXISO exceeded in DISOTP')
      END IF
      IF (IATOM .GT. MAXCHR) THEN
         WRITE (LUPRI,'(//,A,2(/,A,I5),A)')
     &         ' IATOM too large in DISOTP.',
     &         ' Input value:  ',IATOM,
     &         ' Maximum value:',MAXCHR,
     &         ' Program cannot continue.'
         CALL QUIT('MAXCHR exceeded in DISOTP')
      END IF
C
      IF (IATOM .LE. 0) THEN
C        This is a floating orbital, a point charge,
C        or an auxiliary basis set /Mar 2004 hjaaj
         DISOTP = D0
      ELSE IF (TYPE .EQ. 'MASS') THEN
         DISOTP = DATNUC(1,ISOTOP,IATOM)
      ELSE IF (TYPE .EQ. 'MASS_in_AU') THEN
         DISOTP = DATNUC(1,ISOTOP,IATOM) * XFAMU
      ELSE IF (TYPE .EQ. 'A') THEN
         DISOTP = NINT(DATNUC(1,ISOTOP,IATOM))
      ELSE IF (TYPE .EQ. 'ABUNDANCE') THEN
         DISOTP = DATNUC(2,ISOTOP,IATOM)
      ELSE IF (TYPE .EQ. 'SPIN') THEN
         DISOTP = DATNUC(3,ISOTOP,IATOM)
      ELSE IF (TYPE .EQ. 'MMOM') THEN
         DISOTP = DATNUC(4,ISOTOP,IATOM)
      ELSE IF (TYPE .EQ. 'GVAL') THEN
         SPIN = DATNUC(3,ISOTOP,IATOM)
         IF (SPIN .GT. THRESH) THEN
            DISOTP = DATNUC(4,ISOTOP,IATOM)/SPIN
         ELSE
            DISOTP = D0
         END IF
      ELSE IF (TYPE .EQ. 'LARMOR') THEN
         SPIN = DATNUC(3,ISOTOP,IATOM)
         IF (SPIN .GT. THRESH) THEN
            DISOTP = ABS(ECHARGE*DATNUC(4,ISOTOP,IATOM)
     &                   /(D4*PI*SPIN*UMASS*PMASS))
         ELSE
            DISOTP = D0
         END IF
      ELSE IF (TYPE .EQ. 'QMOM') THEN
         DISOTP = DATNUC(5,ISOTOP,IATOM)
      ELSE IF (TYPE .EQ. 'NEUTRONS') THEN
         No_of_NEUTRONS = NINT(DATNUC(1,ISOTOP,IATOM)) - IATOM
         DISOTP = No_of_NEUTRONS
      ELSE
         WRITE (LUPRI,'(//,3A,/,A)')
     &         ' Keyword ',TYPE,' unknown in DISOTP.',
     &         ' Program cannot continue.'
         CALL QUIT('Illegal keyword in DISOTP')
      END IF
      RETURN
      END
C  /* Deck unpkij */
      SUBROUTINE UNPKIJ(IJ,I,J)
C
C     This routine unpacks triangular index, inspired by P. Taylor
C     T. Helgaker & H. Koch
C     Wed Nov 27 10:38:04  1991
C
#include "implicit.h"
      INTEGER :: I,J
C
      RIJ = 2*IJ
      I = INT(SQRT(RIJ + 0.25D00) + 0.4999999D00)
      J = IJ - I*(I - 1)/2
      RETURN
C
      END
C  /* Deck pr2der */
      SUBROUTINE PR2DER(CHESS,NROW,NDIM,LUPRI)
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use nuclei_mod, only: NAMDPX
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,J,K
      INTEGER, PARAMETER :: KCOL = 6
      REAL(dp), PARAMETER :: D0 = 0.0D0
      INTEGER BEGIN
      DIMENSION CHESS(NDIM,NDIM)
#include "nuclei.h"
      LAST = MIN(NROW,KCOL)
      BEGIN= 1
 51   CONTINUE
      WRITE (LUPRI,1000) (NAMDPX(I),I = BEGIN,LAST)
      WRITE (LUPRI,'()')
      NCOL = 1
      DO 101 K = BEGIN,NROW
         DO 201 I = 1,NCOL
            IF (CHESS(K,(BEGIN-1)+I) .NE. D0) GO TO 401
 201     CONTINUE
         GO TO 301
 401     WRITE (LUPRI,2000) NAMDPX(K),
     *        (CHESS(K,(BEGIN-1)+J),J=1,NCOL)
         IF (MOD(K,3) .EQ. 0) WRITE (LUPRI,'()')
 301     IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
 101  CONTINUE
      WRITE (LUPRI,'()')
      LAST = MIN(LAST+KCOL,NROW)
      BEGIN= BEGIN+NCOL
      IF (BEGIN.LE.NROW) GO TO 51
      RETURN
 1000 FORMAT (8X,6(3X,A6,3X),(3X,A6,3X))
 2000 FORMAT (1X,A6,6F12.6)
      END
C  /* Deck fcpri */
      SUBROUTINE FCPRI(AMAT,KEY,CSTRA,SCTRA)
C
C     tuh 131288
C
      use symmet_mod, only: IPTCNT, IPTCOR, IPTXYZ, MAXREP, NAXREP,
     &  NCRREP
      use priunit_mod, only: LUPRI
!      use nuclei_mod, only: NUCDEP, NUCIND
!      use codata_mod, only: ALPHA2
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxaqn_mod, only: MXAQN, MXQN
      use maxorb_mod, only: MXCORB
!      use chrxy_mod, only: CHRXYZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,K
#include "codata.h"
#include "chrxyz.h"
      REAL(dp), PARAMETER :: D1 = 1.D0, D10E6 = 1.D6
      LOGICAL DONE, CARTRA
      CHARACTER FIELD*1, KEY*(*)
      DIMENSION AMAT(3,MXCOOR), CMAT(3,MXCOOR), CSTRA(*), SCTRA(*)
!#include "abainf.h"
#include "nuclei.h"
C
      CARTRA = .TRUE.
C
C     Atomic polar tensors (au)
C
      IF (KEY .EQ. 'APT') THEN
         ITYPEF = 1
         ITYPEC = 1
         FIELD  = 'E'
         FACTOR = D1
C
C     Atomic axial tensors (au)
C
      ELSE IF (KEY .EQ. 'AAT') THEN
         ITYPEF = 2
         ITYPEC = 1
         FIELD  = 'B'
         FACTOR = D1
C
C     Nuclear shieldings (ppm)
C
      ELSE IF (KEY .EQ. 'SIGMA') THEN
         ITYPEF = 2
         ITYPEC = 2
         FIELD  = 'B'
         FACTOR = D10E6*ALPHA2
      ELSE IF (KEY .EQ. 'SIGMAS') THEN
         CARTRA = .FALSE.
         ITYPEF = 2
         ITYPEC = 2
         FIELD  = 'B'
         FACTOR = D10E6*ALPHA2
cLig added instructions how to print shieldings with CTOCD
C
C     Nuclear shieldings (ctocd)
C
      ELSE IF (KEY .EQ. 'CTOCD') THEN
         CARTRA = .FALSE.
         ITYPEF = 2
         ITYPEC = 2
         FIELD  = 'B'
         FACTOR = D1
cLig
C
C     Nuclear shieldings (unscaled)
C
      ELSE IF (KEY .EQ. 'SIGMANO') THEN
         ITYPEF = 2
         ITYPEC = 2
         FIELD  = 'B'
         FACTOR = D1
      ELSE
         WRITE (LUPRI,'(//A/2A/)') ' ERROR in FCPRI ',
     &     '    - Illegal keyword ',KEY
         CALL QUIT('Illegal keyword in FCPRI')
      END IF
      IOFFAX = 0
      DONE = .FALSE.
      DO 100 IREP = 0, MAXREP
         NAXIS = NAXREP(IREP,ITYPEF)
#if defined (PRG_DIRAC)
         IF (NCRREP(IREP,ITYPEC) .GT. 0) THEN
#else
         IF (DOSYM(IREP + 1) .AND. (NAXIS .GT. 0)
     &                       .AND. (NCRREP(IREP,ITYPEC) .GT. 0)) THEN
#endif
            DONE = .TRUE.
            IF (MAXREP .GT. 0) THEN
               WRITE (LUPRI,'(/4X,A,I1)') 'Symmetry ', IREP + 1
            END IF
            WRITE (LUPRI,'(23X,3(2A,13X))')
     &         (FIELD,CHRXYZ(-IPTXYZ(I,IREP,ITYPEF)),I=1,NAXIS)
            WRITE (LUPRI,'()')
            DO 200 IATOM = 1, NUCIND
               DO 300 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,ITYPEC)
                  IF (ISCOOR .GT. 0) THEN
                     WRITE (LUPRI,'(8X,A6,3F15.8)')
     &                  NAMEX(IPTCOR(ISCOOR,ITYPEC)),
     &                  (FACTOR*AMAT(IOFFAX + I,ISCOOR), I=1,NAXIS)
                  END IF
  300          CONTINUE
               IF (MAXREP .EQ. 0) WRITE (LUPRI,'()')
  200       CONTINUE
         END IF
         IOFFAX = IOFFAX + NAXIS
  100 CONTINUE
      IF (.NOT.DONE) THEN
         WRITE (LUPRI,'(2X,A)') ' No elements calculated - '//
     &                         ' appropriate symmetries not requested.'
      END IF
      WRITE (LUPRI,'(/)')
      IF ((MAXREP.GT.0) .AND. DONE .AND. CARTRA) THEN
         CALL TRADIP(AMAT,CMAT,CSTRA,SCTRA,3*NUCDEP,ITYPEF,ITYPEC)
         WRITE (LUPRI,'(/,23X,3(A,13X),/)')
     &               FIELD//'x', FIELD//'y', FIELD//'z'
         DO 400 ICOOR = 1, 3*NUCDEP
            WRITE (LUPRI,'(8X,A6,3F15.8)') NAMDPX(ICOOR),
     &                              (FACTOR*CMAT(K,ICOOR),K=1,3)
            IF (MOD(ICOOR,3) .EQ. 0) WRITE (LUPRI,'()')
 400     CONTINUE
      END IF
      RETURN
      END
C  /* Deck tradip */
      SUBROUTINE TRADIP(SDIPDR,CDIPDR,CSTRA,SCTRA,NCOOR,ITYPEF,ITYPEC)
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTAX
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D1 = 1.0D0
      DIMENSION SDIPDR(3,MXCOOR), CDIPDR(3,MXCOOR),
     *          XMATRI(3,MXCOOR), PERDIP(3,3),
     *          CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR)
C
      CALL DZERO(PERDIP,9)
      DO 100 ICOOR = 1, 3
         PERDIP(ICOOR,IPTAX(ICOOR,ITYPEF)) = D1
  100 CONTINUE
C.....get transformation matrix SCTRA from symmetry to Cartesian coordinates
      CALL TRACOR(CSTRA,SCTRA,ITYPEC,NCOOR,0)
C.....calculate XMATRI = PERDIP * SDIPDR
      CALL DGEMM('N','N',3,NCOOR,3,1.D0,
     &           PERDIP,3,
     &           SDIPDR,3,0.D0,
     &           XMATRI,3)
C.....calculate CDIPDR = XMATRI*SCTRA(t)
      CALL DGEMM('N','T',3,NCOOR,NCOOR,1.D0,
     &           XMATRI,3,
     &           SCTRA(1,1),NCOOR,0.D0,
     &           CDIPDR,3)
      RETURN
      END
C  /* Deck dp0pri */
      SUBROUTINE DP0PRI(DIPOLE)
!
! prints the nuclear dipole moment
!
      use symmet_mod, only: ISYMAX, NAXREP
!      use codata_mod, only: DEBYE
      use mxcent_mod, only: MXCENT
      use maxaqn_mod, only: MXAQN, MXQN
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB
!      use chrxyz_mod, only: CHRXYZ
#include "implicit.h"
      INTEGER :: I
#include "chrxyz.h"
#include "codata.h"
!#include "cbidip.h"
      DIMENSION DIPOLE(3)
      IF (NAXREP(0,1) .EQ. 0) THEN
         WRITE (LUPRI,'(20X,A)')
     &      'All dipole components are zero by symmetry'
      ELSE
         WRITE (LUPRI,'(31X,A/)') 'au             Debye'
         DO 100 I = 1, 3
            IF (ISYMAX(I,1) .EQ. 0) THEN
               WRITE (LUPRI,'(20X,A,2F16.8)') CHRXYZ(-I),
     *                DIPOLE(I), DEBYE*DIPOLE(I)
            END IF
  100    CONTINUE
         WRITE (LUPRI,'(/,23X,A,F12.8,A)') ' 1 a.u. = ',DEBYE,' Debye'
      END IF
      !WRITE (LUPRI,'(//)')
      RETURN
      END
C  /* Deck cmmass */
      SUBROUTINE CMMASS(GEOM,AMASS,NATTYP,NUMIS,IPRINT)
C
C     Determine center of mass in orgcom.h(CMXYZ) from information in common blocks
C     K.Ruud, June-94
C
C
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
!      use codata_mod, only: XFAMU
!      use orgcom_mod, only: CMXYZ, ORIGIN
!      use nuclei_mod, only: IZATOM, NATOMS, NUCIND
      use symmet_mod, only: ISYMAX, MAXOPR, PT
#include "implicit.h"
#include "codata.h"
      REAL*8    GEOM(3*NATOMS), AMASS(NATOMS)
      INTEGER   NATTYP(NATOMS), NUMIS(NATOMS)
#include "nuclei.h"
#include "orgcom.h"
C
      ICHTOT = 0
      JATOM  = 0
C
      DO 10 IATOM = 1, NUCIND
         DO 20 ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = IZATOM(IATOM)
               IF (ICHARG .GT. 0) THEN
                  ICHTOT = ICHTOT + ICHARG
                  DO 30 KCOOR = 1, 3
                     GEOM(3*JATOM + KCOOR) =
     &                    PT(IAND(ISYMAX(KCOOR,1),ISYMOP))
     &                    *CORD(KCOOR,IATOM)
 30               CONTINUE
                  JATOM = JATOM + 1
                  NATTYP(JATOM) = ICHARG
                  NUMIS(JATOM)  = ISOTOP(IATOM)
               END IF
            END IF
 20      CONTINUE
 10   CONTINUE
C
      CALL VIBMAS(AMASS,TOTMAS,NUMIS,NATTYP,JATOM,GEOM,CMXYZ,IPRINT-1)


      TELCMAS = ICHTOT/XFAMU
      TNUCMAS = TOTMAS - TELCMAS
C     ... Total NUClear MASs = TOTal MASs - Total ELeCtron MASs
      RETURN
      END
C  /* Deck vibmas */
      SUBROUTINE VIBMAS(AMASS,TOTMAS,NUMIS,NATTYP,NNATOM,
     &                  GEOM,ORIGIN,IPRINT)
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use cbisol_mod, only: SOLVNT
      use priunit_mod, only: LUPRI
!      use codata_mod, only: XTANG
!      use nuclei_mod, only: NAMDEP
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: K
#include "codata.h"
C
C     Sets up masses and calculates total mass for an isotopic
C     molecule
C
      REAL(dp), PARAMETER :: D0 = 0.0D0, D100 = 100.0D0, DP01 = 0.01D0,
     &                       BIGMAS = 1.0D20
C
      DIMENSION NUMIS(*), NATTYP(*), AMASS(*), GEOM(3,NNATOM), ORIGIN(3)
C
#include "cbisol.h"
#include "nuclei.h"
C
      IF (IPRINT .GT. 0) CALL HEADER('Isotopic Masses',-1)
      TOTMAS = D0
      ABUND  = D100
      CALL DZERO(ORIGIN,3)
      DO 100 IATOM = 1, NNATOM
         IF (SOLVNT .AND. IATOM .EQ. NNATOM) THEN
            RMASS = BIGMAS
            AMASS(IATOM) = RMASS
         ELSE
            RMASS = DISOTP(NATTYP(IATOM),NUMIS(IATOM),'MASS')
            AMASS(IATOM) = RMASS
            TOTMAS       = TOTMAS + RMASS
            QMASS = DISOTP(NATTYP(IATOM),NUMIS(IATOM),'ABUNDANCE')
            ABUND        = ABUND*DP01*QMASS
            ORIGIN(1) = ORIGIN(1) + RMASS*GEOM(1,IATOM)
            ORIGIN(2) = ORIGIN(2) + RMASS*GEOM(2,IATOM)
            ORIGIN(3) = ORIGIN(3) + RMASS*GEOM(3,IATOM)
         END IF
         IF (IPRINT .GT. 0)
     &      WRITE (LUPRI,'(27X,A6,2X,F12.6)')  NAMDEP(IATOM), RMASS
 100  CONTINUE
      IF(TOTMAS.NE.D0) THEN
        ORIGIN(1) = ORIGIN(1)/TOTMAS
        ORIGIN(2) = ORIGIN(2)/TOTMAS
        ORIGIN(3) = ORIGIN(3)/TOTMAS
      ENDIF
      IF (IPRINT .GE. 0) THEN
        WRITE(LUPRI,'(/23X,A,F12.6,A)')'* Total mass: ',TOTMAS,' amu'
        WRITE(LUPRI,'( 23X,A,F8.3,A)') '* Natural abundance:',ABUND,' %'
        WRITE (LUPRI,'(//,(A,4X,3F20.15)))')
     &    '* Center-of-mass coordinates (a.u.):',(ORIGIN(K),K=1,3),
     &    '* Center-of-mass coordinates (A)   :',(ORIGIN(K)*XTANG,K=1,3)
      END IF
      RETURN
      END

C  /* Deck chargebc */
      SUBROUTINE CHARGEBC
C     Determine nuclei charge barycenter in orgcom.h(NCHBC) from information in common blocks
C     J. J. Aucar, 2024

      use symmet_mod, only: ISYMAX, MAXOPR, PT
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use maxaqn_mod, only: MXAQN, MXQN
      use priunit_mod, only: LUPRI
      use memory_allocator
      IMPLICIT NONE

#include "nuclei.h"
!#include "codata.h"
#include "orgcom.h"

C     Variables Declaration
      INTEGER :: IATOM, ISYMOP
      REAL*8  :: CHTOT

C     Initialize array
      NCHBC(1:3) = 0.0d0

C     Initialize total charge to Zero
      CHTOT = 0

C     Nuclei charge barycenter calculation
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               !If charge of IATOM is greater than 0 (i.e. if not point charge)
               IF (IZATOM(IATOM) .GT. 0) THEN
                  ! use effective nuclear charge in CHARGE(:) (less than real charge for ECPs)
!      use orgcom_mod, only: NCHBC, ORIGIN
!      use nuclei_mod, only: CHARGE, CORD, ISTBNU, IZATOM, NUCIND
                  ! instead of real nuclear charge in IZATOM(:)
                  CHTOT = CHTOT + CHARGE(IATOM)
                  NCHBC(1) = NCHBC(1) + CHARGE(IATOM)*
     &                    PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  NCHBC(2) = NCHBC(2) + CHARGE(IATOM)*
     &                    PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  NCHBC(3) = NCHBC(3) + CHARGE(IATOM)*
     &                    PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
               END IF
            END IF
         END DO
      END DO

      IF (CHTOT.NE.0.0D0) THEN
         NCHBC = NCHBC / CHTOT
      ELSE
         WRITE (LUPRI,'(//A/2A//)') ' ERROR in CHARGEBC ',
     &          '     - Total nuclei charge equals zero.'
         CALL QUIT('Total nuclei charge equals zero.')
      END IF

      RETURN
      END

      SUBROUTINE TRACOR(CSTRA,SCTRA,ITYPE,NCOOR,IPRINT)
C
C     Sets up transformation matrices between Cartesian (ITYPE=1)
C     or rotational (ITYPE=2) to symmetry nuclear coordinates
C
C     CSTRA - transformation matrix Cartesian to symmetry coordinates
C     SCTRA - inverse matrix
C
C     tuh Jun 29 1988
C
!      use ibtfun_mod, only: IBTXOR
!      use nuclei_mod, only: ISTBNU, NUCIND, NUCNUM
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: IPTCNT, ISYMAX, MAXOPR, MAXREP, NCRREP, PT
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,J
      REAL(dp), PARAMETER :: D0 = 0.0D0
C
#include "nuclei.h"
C
      DIMENSION CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR)
C
#include "ibtfun.h"
C
      CALL DZERO(CSTRA,NCOOR*NCOOR)
      CALL DZERO(SCTRA,NCOOR*NCOOR)
C
C     *****************************************************************
C     **** Transformation matrix Cartesian to symmetry coordinates ****
C     *****************************************************************
C
      DO 100 IREP = 0, MAXREP
         NCRREP(IREP,ITYPE) = 0
         DO 110 ICENT = 1, NUCIND
            MULCNT = ISTBNU(ICENT)
            DO 120 ICOOR = 1, 3
               IF (ITYPE .EQ. 1) THEN
                  ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     NCRREP(IREP,1) = NCRREP(IREP,1) + 1
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,1))
                     DO 130 ISYMOP = 0, MAXOPR
                        IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                           ICCOOR = 3*(NUCNUM(ICENT,ISYMOP+1) - 1)+ICOOR
                           CSTRA(ISCOOR,ICCOOR)=PT(IBTAND(ISYMOP,IVAR))
                        END IF
 130                 CONTINUE
                  END IF
               ELSE
                  LSCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,2)
                  IF (LSCOOR .GT. 0) THEN
                     NCRREP(IREP,2) = NCRREP(IREP,2) + 1
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,2))
                     DO 140 ISYMOP = 0, MAXOPR
                        IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                           ICCOOR = 3*(NUCNUM(ICENT,ISYMOP+1) - 1)+ICOOR
                           CSTRA(LSCOOR,ICCOOR)=PT(IBTAND(ISYMOP,IVAR))
                        END IF
 140                 CONTINUE
                  END IF
               END IF
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF (IPRINT .GE. 20) THEN
         IF (ITYPE .EQ. 1) THEN
            WRITE (LUPRI,'(/,1X,A,8I3,/)')
     *           ' Number of coordinates in each symmetry: ',
     *           (NCRREP(I,1),I=0,MAXREP)
            CALL HEADER('Transformation from Cartesian to symmetry '
     *                //'coordinates',-1)
         ELSE
            CALL HEADER('Transformation from rotational cartesian '//
     &                  'to symmetry coordinates',-1)
         END IF
         DO 200 ISYM = 1, NCOOR
            WRITE (LUPRI,'(9F6.2)') (CSTRA(ISYM,J),J=1,NCOOR)
 200     CONTINUE
      END IF
C
C     Inverted matrix
C
      DO 210 I = 1, NCOOR
         DO 220 J = 1, NCOOR
           SCTRA(I,J) = CSTRA(J,I)
  220    CONTINUE
  210 CONTINUE
      DO 230 ICOL = 1, NCOOR
         SUM1 = D0
         SUM2 = D0
         DO 240 IROW = 1, NCOOR
            SUM1 = SUM1 + ABS(SCTRA(IROW,ICOL))
  240    CONTINUE
         DO 250 IROW = 1, NCOOR
            SCTRA(IROW,ICOL) = SCTRA(IROW,ICOL)/SUM1
  250    CONTINUE
  230 CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETMSS(IZ,DMASS,ABUND,INFO)
C***********************************************************************
C
C     Get mass of the most abundant isotope of charge IZ
C
C***********************************************************************
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0=0.0D0
      PARAMETER (ISO=1, MAXCHR=138)
C
      INFO = -1
      IF(IZ.GT.MAXCHR) THEN
        RETURN
      ENDIF
C     DISOTP returns the most abundant isotope
C     for ISO = 1.
      ABUND = DISOTP(IZ,ISO,'ABUNDANCE')
      IF (ABUND .GT. D0) THEN
         INFO = 0
         DMASS = DISOTP(IZ,ISO,'MASS')
      END IF
C
      RETURN
      END
