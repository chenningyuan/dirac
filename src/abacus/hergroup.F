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

C FILE : hergroup.F
C
C Purpose: set up point group symmetry information
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck syminp */
      SUBROUTINE SYMINP(NSYMOP,KASYM,IFXYZ,CLASS)
C*****************************************************************************
C
C     This subroutine sets up point group symmetry and looks at
C     the behaviour of principal axes and rotations under point group
C     symmetry
C
C     tsaue - 940825 - major revision: tagged on SYMGRP + polish
C
C*****************************************************************************
!      use cbirea_mod, only: DIRAC, IPREAD
      use symmet_mod, only: FMULT, IPTAX, IPTXYZ, ISYMAX, MAXREP, MULT,
     &  NAXREP, PT
      use maxaqn_mod, only: MXAQN, MXQN
!      use ibtfun_mod, only: IBTAND, IBTSHR, IBTXOR
      use maxorb_mod, only: MXCORB
      use priunit_mod, only: LUPRI
      use mxcent_mod, only: MXCENT
#include "implicit.h"
      INTEGER :: I,J,K
C
#include "cbirea.h"
      CHARACTER*1 KASYM(3,3)
      CHARACTER*(*) CLASS
      DIMENSION IFXYZ(3),IS(0:7),IGEN(3)
      DATA IS/0,1,1,2,1,2,2,3/
#include "ibtfun.h"
C
C     Initialization
C     ==============
C
C     PT is parity of a bitstring:
C       1 for an even number of ones: 000,011,110,101
C      -1 for an odd  number of ones: 001,010,100,111
C
      PT(0) =  1.0D0
      PT(1) = -1.0D0
      PT(2) = -1.0D0
      PT(3) =  1.0D0
      PT(4) = -1.0D0
      PT(5) =  1.0D0
      PT(6) =  1.0D0
      PT(7) = -1.0D0
      DO 5 I = 1,3
         IFXYZ(I)    = 0
         ISYMAX(I,1) = 0
         IGEN(I)     = 0
    5 CONTINUE
C
C     Determine:
C     IGEN(I)   - basic operations
C     ISYMAX(I) - behavior of principal axes under basic operations
C     ===============================================================
C
      IAXIS = 0
      MAXREP = 2**NSYMOP - 1
      IF (NSYMOP.GT.0.AND.IPREAD .GT. 0) THEN
         CALL HEADER('Symmetry Operations',1)
         WRITE (LUPRI,'(1X,A,I2,/)') ' Symmetry operations:',NSYMOP
      END IF
      DO 100 J = 1,NSYMOP
        DO 110 I = 1,3
        IF (KASYM(I,J).NE.' ') THEN
          K = ICHAR(KASYM(I,J)) - ICHAR('W')
          IGEN(J)     = IGEN(J)     + 2**(K-1)
          ISYMAX(K,1) = ISYMAX(K,1) + 2**(J-1)
        END IF
  110   CONTINUE
        IAXIS = IBTOR(IAXIS,IGEN(J))
  100 CONTINUE
C
C     Determine IFXYZ
C     ===============
C     Do we really need it ????
C
      IND = 0
      DO 140 I = 1,NSYMOP
        IND = IBTOR(IND,IGEN(I))
  140 CONTINUE
      DO 141 I = 1,3
        IFXYZ(I) = IBTAND(IBTSHR(IND,(I-1)),1)
  141 CONTINUE
C
C     Determine:
C     ISYMAX(I,2) - behaviour of principal rotations under basic operations
C     =====================================================================
C
      ISYMAX(1,2) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
      ISYMAX(2,2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
      ISYMAX(3,2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
C
C     Determine:
C     IPTAX   - coordinate axis: pointer analogous to IPTSYM
C     NAXREP  - number of coordinate axis in each symmetry
C     ========================================================
C
      CALL IZERO(IPTAX,6)
      CALL IZERO(IPTXYZ(1,0,1),48)
      DO 200 ITYPE = 1,2
        IPTAXI = 0
        DO 205 IREP = 0, MAXREP
           NAXIS = 0
           DO 210 ICOOR = 1, 3
              IF (IBTXOR(IREP,ISYMAX(ICOOR,ITYPE)) .EQ. 0) THEN
                 NAXIS  = NAXIS + 1
                 IPTAXI = IPTAXI + 1
                 IPTAX (ICOOR,ITYPE)      = IPTAXI
                 IPTXYZ(NAXIS,IREP,ITYPE) = ICOOR
              END IF
  210      CONTINUE
           NAXREP(IREP,ITYPE) = NAXIS
  205   CONTINUE
  200 CONTINUE
      IF (IPREAD.GT. 20) THEN
         WRITE (LUPRI,'(1X,A,3I5)') ' IPTAX(*,1)  ',
     &                               (IPTAX(I,1),I=1,3)
         WRITE (LUPRI,'(1X,A,3I5)') ' IPTAX(*,2)  ',
     &                               (IPTAX(I,2),I=1,3)
         WRITE (LUPRI,'(1X,A,8I5)') ' NAXREP(*,1) ',
     &                               (NAXREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(1X,A,8I5)') ' NAXREP(*,2) ',
     &                               (NAXREP(I,2),I=0,MAXREP)
      END IF
C
C     Determine group & properties
C
      CALL SYMGRP(IGEN,NSYMOP,CLASS)
C
C     For DIRAC: double group symmetry
C
      IF(DIRAC) CALL DBLGRP
C
C     Determine MULT(I) and FMULT(I) - multiplicity of center
C     =======================================================
C
      DO 101 I = 0,7
         MULT(I)  = 2**MAX(0,NSYMOP-IS(I))
         FMULT(I) = dble(MULT(I))
  101 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck symgrp */
      SUBROUTINE SYMGRP(IGEN,NSYMOP,CLASS)
C*****************************************************************************
C
C     Given the group generators, this subroutine will identify
C     Abelian subgroup and set up group multiplication table,
C     character table and direct product table. Irreps are identified.
C
C     tsaue - august 1994
C     Sep 24 1996 - tsaue : Included GROUPS,JSOP and IPAR in COMMON SYMMET.
C                           Modified activation of LSYMOP
C                           GROUPS and SYMOP are initialized in BLOCK DATA
C
C*****************************************************************************
      use mxcent_mod, only: MXCENT
      use maxorb_mod, only: MXCORB
!      use argoscomch_mod, only: ITYPAG
!      use ibtfun_mod, only: IBTAND, IBTSHR, IBTXOR
      use priunit_mod, only: LUPRI
      use maxaqn_mod, only: MXAQN, MXQN
!      use argoscom_mod, only: IGRPAG, NDPTAG, NDPTAGS
!      use cbirea_mod, only: IPREAD
      use symmet_mod, only: IPAR, IXVAL, JSOP, MAXREP, NINVC, NREFL,
     &  NROTS, PT
!      use pgroup_mod, only: GROUP, GROUPS, REP, SYMOP
      use xmlout
#include "implicit.h"
      INTEGER :: I,J,K
C
#include "pgroup.h"
#include "cbirea.h"
      LOGICAL LSYMOP(0:7)
      DIMENSION IGEN(3),IROTS(3),IREFL(3),JPAR(0:7)
      CHARACTER*(*) CLASS
      DATA IROTS/3,5,6/
      DATA IREFL/4,2,1/
      DATA JPAR/ 1,-1,-1, 1,-1, 1, 1,-1/
#include "ibtfun.h"
#include "argoscom.h"
#include "argoscomch.h"
      DO 5 I = 0,7
        LSYMOP(I) = .FALSE.
    5 CONTINUE
C
C     Activate all symmetry operations of the group
C     =============================================
C
      LSYMOP(0) = .TRUE.
      JSOP(0)   = 0
      IPAR(0)   = 1
      DO 10 I   = 1,MAXREP
        I0 = IBTAND(1,I)          *IGEN(1)
        I1 = IBTAND(1,IBTSHR(I,1))*IGEN(2)
        I2 = IBTAND(1,IBTSHR(I,2))*IGEN(3)
        IND = IBTXOR(IBTXOR(I0,I1),I2)
        LSYMOP(IND) = .TRUE.
        IPAR(I)     = JPAR(IND)
   10 CONTINUE
C
C     List group operations in preferred order,
C       that is construct pointer LSOP
C     =========================================
C
C     Identity
C
      IND = 0
      JSOP(IND) = 0
C
C     Rotations
C
      NROTS = 0
      DO 40 I = 1,3
      IF(LSYMOP(IROTS(I))) THEN
        IND         = IND + 1
        JSOP(IND) = IROTS(I)
        NROTS       = NROTS + 1
      ENDIF
   40 CONTINUE
C
C     Inversion
C
      NINVC = 0
      IF(LSYMOP(7)) THEN
        IND         = IND + 1
        JSOP(IND) = 7
        NINVC       = 1
      ENDIF
C
C     Reflections
C
      NREFL = 0
      DO 50 I = 1,3
      IF(LSYMOP(IREFL(I))) THEN
        IND         = IND + 1
        JSOP(IND) = IREFL(I)
        NREFL       = NREFL + 1
      ENDIF
   50 CONTINUE
      IF(IND.NE.MAXREP) CALL QUIT('SYMGRP:IND.NE.MAXREP!')
C
C     Classify group
C     ==============
C     tsaue - Here I have devised a highly empirical formula, but
C             it works !!!
C
      IGROUP = MIN(7,NINT((4*NROTS+8*NINVC+6*NREFL)/3.0))
      GROUP  = GROUPS(IGROUP)
C
C     export it to ECP module
C     =======================
C
      ndptag = ndptags(IGROUP)
      IGRPAG = IGROUP
C
C     Generate character table
C     ========================
C
      DO 60 I = 0,MAXREP
        IXVAL(0,I) = 1
        DO 70 J = 1,NSYMOP
          IXVAL(IGEN(J),I) = PT(IBTAND(IBTSHR(I,(J-1)),1))
          DO 80 K = 1,(J-1)
            IND      = IBTXOR(IGEN(J),IGEN(K))
            IXVAL(IND,I)  = IXVAL(IGEN(J),I)*IXVAL(IGEN(K),I)
            DO 90 L = 1,(K-1)
              IXVAL(IBTXOR(IND,IGEN(L)),I)
     &           = IXVAL(IND,I)*IXVAL(IGEN(L),I)
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Classify irrep
C     ==============
C
      DO 100 I = 0,MAXREP
        REP(I) = 'A  '
        IPOS = 2
C
C       Rotational symmetry
C
        IF(NROTS.EQ.3) THEN ! D2h, D2
          IND = (1-IXVAL(JSOP(1),I))+(1-IXVAL(JSOP(2),I))/2
          IF(IND.NE.0) THEN
            REP(I)(1:1) = 'B'
            REP(I)(2:2) = CHAR(ICHAR('0')+IND)
            IPOS = 3
          ENDIF
        ELSEIF(NROTS.EQ.1) THEN ! C2, C2v, C2h
          IF(IXVAL(JSOP(1),I).EQ.-1) REP(I)(1:1) = 'B'
          IF(NREFL.EQ.2) THEN ! C2v
            IF(IBTAND(IBTSHR(JSOP(1),1),1).EQ.1) THEN
              IND = 2
            ELSE
              IND = 3
            ENDIF
            IF(IXVAL(JSOP(IND),I).EQ.1) THEN
              REP(I)(2:2) = '1'
            ELSE
              REP(I)(2:2) = '2'
            ENDIF
          ENDIF
        ELSEIF(NREFL.EQ.1) THEN ! Cs
C
C       Mirror symmetry
C
          IF(IXVAL(JSOP(1),I).EQ.1) THEN
            REP(I)(2:2) = ''''
          ELSEIF(IXVAL(JSOP(1),I).EQ.-1) THEN
            REP(I)(2:2) = '"'
          ENDIF
        ENDIF
C
C       Inversion symmetry
C
        IF(NINVC.EQ.1) THEN
          IND = NROTS+1
          IF(IXVAL(JSOP(IND),I).EQ.1) THEN
            REP(I)(IPOS:IPOS) = 'g'
          ELSE
            REP(I)(IPOS:IPOS) = 'u'
          ENDIF
        ENDIF
  100 CONTINUE
C
C     Output section
C     ==============
C
      IF (doxml) THEN
         call xml_begin('point-group')
         write(luxml,*) GROUP
         call xml_end('point-group')
      ENDIF

      IF(IPREAD.GT.0) THEN
        CALL HEADER('SYMGRP:Point group information',-1)
C
C       Group name
C
        IF (CLASS(1:3) .EQ. 'N/A') THEN
           WRITE(LUPRI,'(A,A3)') 'Point group: ',GROUP
        ELSE
           WRITE(LUPRI,'(A,A)') 'Full group is:  ',CLASS
           WRITE(LUPRI,'(A,A3)') 'Represented as: ',GROUP
        END IF
C
C       Group generators
C
        IF(NSYMOP.GT.0) THEN
          WRITE(LUPRI,'(/3X,A/)')
     &        '* The point group was generated by:'
          DO 110 I = 1,NSYMOP
          IF    (SYMOP(IGEN(I))(1:1).EQ.'C') THEN
            WRITE(LUPRI,'(6X,3A)')
     &        'Rotation about the ',SYMOP(IGEN(I))(3:3),'-axis'
          ELSEIF(SYMOP(IGEN(I))(1:1).EQ.'O') THEN
            WRITE(LUPRI,'(6X,3A)')
     &        'Reflection in the ',SYMOP(IGEN(I))(2:3),'-plane'
          ELSE
            WRITE(LUPRI,'(6X,A)')
     &        'Inversion centre'
          ENDIF
  110     CONTINUE
C
C         Group multiplication table
C
          WRITE(LUPRI,'(/3X,A/)') '* Group multiplication table'
          WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(SYMOP(JSOP(I)),
     &        I = 0,MAXREP)
          WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
          DO 120 I = 0,MAXREP
            WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,A3,1X))')
     &        SYMOP(JSOP(I)),'|',
     &        (SYMOP(IBTXOR(JSOP(I),JSOP(J))),J = 0,MAXREP)
  120     CONTINUE
        ENDIF
C
C       Character table
C
        WRITE(LUPRI,'(/3X,A/)') '* Character table'
        WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))') '|',(SYMOP(JSOP(J)),
     &               J = 0,MAXREP)
        WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
        DO 130 I = 0,MAXREP
          WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,I3,1X))')
     &      REP(I),'|',(IXVAL(JSOP(J),I),J=0,MAXREP)
c
c     For ECP module
c     ==============
            itypag(i) = REP(I)
c
  130   CONTINUE
C
C       Direct product table
C
        WRITE(LUPRI,'(/3X,A/)') '* Direct product table'
        WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(REP(I),I = 0,MAXREP)
        WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
        DO 140 I = 0,MAXREP
          WRITE(LUPRI,'(3X,1X,A3,1X,A1,8(1X,A3,1X))')
     &      REP(I),'|',(REP(IBTXOR(I,J)),J = 0,MAXREP)
  140   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dblgrp */
      SUBROUTINE DBLGRP
C*****************************************************************************
C
C     This routine will analyze the fermion irreps in the molecular symmetry
C     and, if possible, set up a transformation within AO - basis to a
C     symmetry-adapted spinor basis. details are given below.
C
C     Written by T.Saue - october 1994 - odense
C     Sep 24 1996 - tsaue: Updated version
C     Febr.2007, M.Ilias - added operator type 20
C
C*****************************************************************************
      use maxaqn_mod, only: MXAQN, MXQN
!      use cbirea_mod, only: IPREAD
      use mxcent_mod, only: MXCENT
!      use dgroup_mod, only: FREP, IBTOSF, IFDIRP, IHQMAT, IPQTOQ, IQDEF, IQMAP, IQPH, IQTOPQ, ITOT, JMROI, JQBAS, JSPINR, NBSYM, NFSYM, NZ, NZ_IN_CI
      use symmet_mod, only: ISYMAX, MAXREP, NINVC
!      use ibtfun_mod, only: IBTAND, IBTSHR, IBTXOR
      use maxorb_mod, only: MXCORB
      use priunit_mod, only: LUPRI
!      use pgroup_mod, only: REP
      use quaternion_algebra
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I,J
      REAL(dp), PARAMETER :: DM1 = -1.0D0, D1 = 1.0D0
C
#include "cbirea.h"
#include "pgroup.h"
#include "dgroup.h"
C
      LOGICAL LBUF(8)
C
#include "ibtfun.h"
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL IZERO(IQTOPQ,32)
      CALL IZERO(IPQTOQ,32)
C
C     Determine number of boson and fermion irreps
C     ============================================
C
      NBSYM = MAXREP + 1
      NFSYM = 1
      IF(NINVC.EQ.1) NFSYM = 2
      IF(NFSYM.EQ.1) THEN
        FREP(1) = 'E1 '
      ELSE
        FREP(1) = 'E1g'
        FREP(2) = 'E1u'
      ENDIF
C
C     Determine NZ - parameter
C     ========================
C       NZ = 4  quaternionic group        no    totally symmetric rotations
C       NZ = 2  complex group             one   totally symmetric rotation
C       NZ = 1  real group                three totally symmetric rotations
C
C       NZ_in_CI is 1 for real CI coefficients, 2 for complex CI coefficients
C
      NZ = 1
      DO 30 I = 1,3
        IF(ISYMAX(I,2).EQ.0) THEN
          NZ   = NZ + 1
          ITOT = I
        ENDIF
   30 CONTINUE
      NZ_in_CI = MIN(2,NZ)
C
C     Align rotations and quaternion irreps,
C     for complex groups this is important....
C
C     The three quaternion units are equivalent and can
C     thus be interchanged. Normally one makes the mapping
C       (x,y,z) --> (k,j,i)
C     and this is the default, including complex groups with
C     a totally symmetric z-rotation
C
C     Variables IXQ, IYQ, and IZQ gives the mapping from coordinates
C     to quaternion units 1, 2, and 3, whereas IQ1, IQ2, and IQ3 gives
C     the reverse mapping.
C
      IQ1 = 3
      IQ2 = 2
      IQ3 = 1
      IZQ = 2
      IYQ = 3
      IXQ = 4
C.....Complex groups
C.....has only one totally symmetric rotation. Since the triple product of the three rotations is
C.....totally symmetric G0, this implies that the other two rotations have the same symmetry G1
C.....Let X0 be a totally symmetric basis function, and X1 a basis function of this second symmetry
      IF(NZ.EQ.2) THEN
C.........totally symmetric z-rotation uses default mapping:
C         (x,y,z) --> (k,j,i)
C         Non-zero coefficients of X0: ( C0 + iCz)
C         Non-zero coefficients of X1: (jCy + kCx) = k(Cx + iCy) (k shifted over to X1, as indicated by JMROI)
        IF(ITOT.EQ.1) THEN
C.........totally symmetric x-rotation
C         (x,y,z) --> (i,k,j)
C         Non-zero coefficients of X0: ( C0 + iCx)
C         Non-zero coefficients of X1: (jCz + kCy) = k(Cy + iCz) (k shifted over to X1, as indicated by JMROI)
          IQ1 = 1
          IQ2 = 3
          IQ3 = 2
          IXQ = 2
          IZQ = 3
          IYQ = 4
        ELSEIF(ITOT.EQ.2) THEN
C.........totally symmetric y-rotation
C         (x,y,z) --> (j,i,k)
C         Non-zero coefficients of X0: ( C0 + iCy)
C         Non-zero coefficients of X1: (jCz + kCx) = k(Cx + iCz) (k shifted over to X1, as indicated by JMROI)
          IQ1 = 2
          IQ2 = 1
          IQ3 = 3
          IYQ = 2
          IXQ = 3
          IZQ = 4
        ENDIF
      ENDIF
C
C     Set up an array that points back to the default assignment of quaternion units
C
      IQMAP(1)   = 1
      IQMAP(IZQ) = 2
      IQMAP(IYQ) = 3
      IQMAP(IXQ) = 4
C
C
C     Set up direct product table for fermion ircops
C     ==============================================
C
C     Consider a Kramers basis of spinors {p,P} where p and P referred to
C     the unbarred and barred spinors of a Kramers pair, respectively.
C     Limiting attention to the binary groups, that is, D2h and subgroups,
C     there are at most two fermion ircops, depending on whether inversion is
C     present or not. We shall use subscripts 1 and 2 to denote the gerade
C     and ungerade representations, respectively.
C
C     Assuming the spinor structure to be given below the boson irrep of the real and imaginary
C     parts of spinor products pq and pQ are given below.
C     For complex groups having a single totally symmetric rotation,
C     the rotation will have the symmetry of Im(p_1,q_1)......
C
C     Gerade product:
C
C     Re(p_1,q_1)
      IFDIRP(1,1) = 0
C     Im(p_1,q_1)
      IFDIRP(2,1) = ISYMAX(IQ1,2)
C     Re(p_1,Q_1)
      IFDIRP(3,1) = ISYMAX(IQ2,2)
C     Im(p_1,Q_1)
      IFDIRP(4,1) = ISYMAX(IQ3,2)
C
C     Ungerade product:
C
C     Re(p_1,q_2)
      IFDIRP(1,2) = IXYZ
C     Im(p_1,q_2)
      IFDIRP(2,2) = ISYMAX(IQ1,1)
C     Re(p_1,Q_2)
      IFDIRP(3,2) = ISYMAX(IQ2,1)
C     Im(p_1,Q_2)
      IFDIRP(4,2) = ISYMAX(IQ3,1)
C
C
C     Determine distribution of boson irreps in spinors
C     =================================================
C     All 4-component operators span boson irreps.
C     JBTOF shows for given operator symmetry IBREP
C     connects spinors of fermion ircops IFSYM and JFSYM,
C     that is, JFSYM = JBTOF(IBREP,IFSYM)
C
C     A 2- (or 4-) spinor spans fermion irreps.
C     However, the real and imaginary parts of the alpha- and beta-
C     components are scalar functions and span boson irreps.
C
C     The array JSPINR(IM,IC,IFSYM) gives the boson irrep
C     of the real and imaginary parts of the alpha/beta-components
C     of the large (IC=1) and small (IC=2) components of
C     2-spinors of fermion ircop IFSYM.
C
C     The default structure of large component spinors are:
C
C       Re[La]: JSPINR(1,1,1) = G0 (totally symmetric)
C       Im[La]: JSPINR(2,1,1) = Rz (z-rotation)       i
C       Re[Lb]: JSPINR(3,1,1) = Ry (y-rotation)       j
C       Im[Lb]: JSPINR(4,1,1) = Rx (x-rotation)       k
C
C     The default structure of small component spinors are:
C
C       Re[Sa]: JSPINR(1,2,1) = xyz
C       Im[Sa]: JSPINR(2,2,1) = z
C       Re[Sb]: JSPINR(3,2,1) = y
C       Im[Sb]: JSPINR(4,2,1) = x
C
C     The choice of Re[La] being totally symmetric corresponds to a choice of phase.
C
C     For complex groups with a totally symmetric rotation other than z, the spinor
C     structure is modified such that the totally symmetric rotation enters Im[La]
C     Note that this corresponds to choosing a different axis than z for spin quantization
C     and so the (Pauli) spin matrices are reordered accordingly.
C
C     Specifically:
C     Totally symmetric x-rotation:
C       Re[La]: JSPINR(1,1,1) = G0 (totally symmetric)
C       Im[La]: JSPINR(2,1,1) = Rx i
C       Re[Lb]: JSPINR(3,1,1) = Rz j
C       Im[Lb]: JSPINR(4,1,1) = Ry k
C     Totally symmetric y-rotation:
C       Re[La]: JSPINR(1,1,1) = G0 (totally symmetric)
C       Im[La]: JSPINR(2,1,1) = Ry i
C       Re[Lb]: JSPINR(3,1,1) = Rx j
C       Im[Lb]: JSPINR(4,1,1) = Rz k
C
C
      DO IM = 1,4
        JBTOF(IFDIRP(IM,1),1) = 1
        JSPINR(IM,1,1) = IFDIRP(IM,1)
        JBTOF(IFDIRP(IM,2),2) = 1
        JSPINR(IM,2,1) = IFDIRP(IM,2)
      ENDDO
C
C     If there is inversion symmetry, the structures of the
C     large and small components are reversed for the ungerade fermion ircop:
      IF(NFSYM.EQ.2) THEN
        DO IM = 1,4
          JBTOF(IFDIRP(IM,2),1) = 2
          JSPINR(IM,1,2) = JSPINR(IM,2,1)
          JBTOF(IFDIRP(IM,1),2) = 2
          JSPINR(IM,2,2) = JSPINR(IM,1,1)
        ENDDO
      ENDIF
C
C     Determine boson symmetries connected by a given
C     fermion irrep
C
      CALL ICOPY(16,-1,0,IBTOSF(0,1),1)
      IF(NFSYM.EQ.1) THEN
        DO ISYM = 1,NBSYM
          JFSYM(ISYM,1) = ISYM
          JFSYM(ISYM,2) = ISYM
          IBTOSF(ISYM-1,1) = ISYM
          IBTOSF(ISYM-1,2) = ISYM
        ENDDO
      ELSEIF(NFSYM.EQ.2) THEN
        ISYMA = 0
        ISYMB = 0
        DO ISYM = 1,NBSYM
          IF(JBTOF(ISYM-1,1).EQ.1) THEN
            ISYMA = ISYMA + 1
            JFSYM(ISYMA,1) = ISYM
            IBTOSF(ISYM-1,1) = ISYMA
          ELSE
            ISYMB = ISYMB + 1
            JFSYM(ISYMB,2) = ISYM
            IBTOSF(ISYM-1,2) = ISYMB
          ENDIF
        ENDDO
      ENDIF
C
C     Matrix packing
C     ==============
C
      JMROI(1) = 1
      IQDEF(1) = 1
      IQDEF(2) = 2
      IQDEF(3) = 3
      IQDEF(4) = 4
C
C     Quaternion groups
C
      IF(NZ.EQ.4) THEN
        JMROI (2) = 1
        JMROI (3) = 1
        JMROI (4) = 1
        IF(NFSYM.EQ.1) THEN
C       C1:
          DO I = 1,4
            IPQTOQ(I,   0) = I
          ENDDO
        ELSE
C       Ci:
          DO I = 1,4
            IPQTOQ(I,   0) = I
            IPQTOQ(I,IXYZ) = I
          ENDDO
        ENDIF
C
C     Complex groups:
C
      ELSEIF(NZ.EQ.2) THEN
C                        ! Z Y X
          JMROI (2) = 1  ! z y x
          JMROI (3) = 4  ! y z z
          JMROI (4) = 4  ! x x y
C         ..non-totally symmetric rotation
          IRNT           = ISYMAX(IQ2,2)
          IPQTOQ(1,0   ) = 1
          IPQTOQ(2,0   ) = 2
          IPQTOQ(1,IRNT) = 4
          IPQTOQ(2,IRNT) = 3
        IF(NFSYM.EQ.2) THEN
          IRNU = IBTXOR(IRNT,IXYZ)
          DO IZ = 1,NZ
            IPQTOQ(IZ,IXYZ) = IPQTOQ(IZ,0)
            IPQTOQ(IZ,IRNU) = IPQTOQ(IZ,IRNT)
          ENDDO
        ENDIF
      ELSE
        JMROI(2)  = 2
        JMROI(3)  = 3
        JMROI(4)  = 4
        IF(NFSYM.EQ.1) THEN
          IPQTOQ(1,0) = 1
          DO I = 2,4
            JMROI(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
          ENDDO
        ELSE
          IPQTOQ(1,   0) = 1
          IPQTOQ(1,IXYZ) = 1
          DO I = 2,4
            JMROI(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
            IPQTOQ(1,ISYMAX(5-I,1)) = I
          ENDDO
        ENDIF
      ENDIF
C
C     Assign IQTOPQ
C
      DO IREP = 0,MAXREP
        DO IPQ = 1,NZ
          IQ = IPQTOQ(IPQ,IREP)
          IQTOPQ(IQ,IREP) = IPQ
        ENDDO
      ENDDO
C
C     Assign JQBAS - gives quaternion vector for given
C     component and bosonirrep
C     ================================================
C
      DO 60 IFSYM = 1,NFSYM
        DO 70 IC = 1,2
          DO 80 IM = 1,4
            JQBAS(JSPINR(IM,IC,IFSYM),IC) = JMROI(IM)
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Give signs of quaternion phases: 1,i,-j,k
C
      IQPH(1,1) =  1
      IQPH(2,1) =  1
      IQPH(3,1) =  1
      IQPH(4,1) =  1
C
C     ...and their Hermitian conjugates: 1,-i, j,-k
C
      IQPH(1,2) =  1
      IQPH(2,2) = -1
      IQPH(3,2) = -1
      IQPH(4,2) = -1
C
C     Matrix symmetry
C     ===============
C
C     The matrix of an operator that is symmetric or
C     antisymmetric under time reversal has the following
C     structure:
C
C       A   B          A^{dagger} =   hA
C     -tB* tA*         B^T        = -thB
C
C     The matrix symmetry of matrices A and B
C     can be summarized as follows:
C
C     IH   ITIM   AR  AI  BR  BI
C      1     1     1   2   2   2
C     -1     1     2   1   1   1
C      1    -1     1   2   1   1
C     -1    -1     2   1   2   2
C
C     Only ITIM=1 is used (ITIM=-1 is transferred to IH=-1
C     by extracting an imaginary i), and the matrix symmetry
C     of matrices A and B is stored in IHQMAT(4,IH):
C
C     H-
      IHQMAT(1,-1) = 2
      IHQMAT(2,-1) = 1
      IHQMAT(3,-1) = 1
      IHQMAT(4,-1) = 1
C     H0
      IHQMAT(1, 0) = 0
      IHQMAT(2, 0) = 0
      IHQMAT(3, 0) = 0
      IHQMAT(4, 0) = 0
C     H+
      IHQMAT(1, 1) = 1
      IHQMAT(2, 1) = 2
      IHQMAT(3, 1) = 2
      IHQMAT(4, 1) = 2
C
C     IRQMAT gives the irrep of a given component of
C     quaternion matrix:
C
      DO IZ = 1,4
        DO IBRP = 0,MAXREP
          IRQMAT(IZ,IBRP) = IBTXOR(IFDIRP(IZ,1),IBRP)
        ENDDO
      ENDDO
C
C     Define generic 4-component operators
C
      CALL MOPDEF(IXYZ,IXQ,IYQ,IZQ)
C
C     Output section
C     ==============
C
      IF(IPREAD.GE.1) THEN
        CALL TITLER('Output from DBLGRP','*',103)
        IF    (NFSYM.EQ.1) THEN
          WRITE(LUPRI,'(3X,A,2X,A3)')
     &      '* One fermion irrep: ',FREP(1)
        ELSEIF(NFSYM.EQ.2) THEN
          WRITE(LUPRI,'(3X,A,2(2X,A3))')
     &      '* Two fermion irreps:',FREP(1),FREP(2)
        ENDIF
        IF    (NZ.EQ.4) THEN
          WRITE(LUPRI,'(3X,A)') '* Quaternionic group. NZ = 4'
        ELSEIF(NZ.EQ.2) THEN
          WRITE(LUPRI,'(3X,A)') '* Complex group. NZ = 2'
        ELSEIF(NZ.EQ.1) THEN
          WRITE(LUPRI,'(3X,A)') '* Real group. NZ = 1'
        ENDIF
        WRITE(LUPRI,'(3X,A)') '* Direct product decomposition:'
        IND = 1
        DO 110 I1 = 1,NFSYM
          DO 120 I2 = I1,NFSYM
            IND = IND + 1
            J = MOD(IND,2)+1
            WRITE(LUPRI,'(10X,2(A3,A),3(A3,A),A3)')
     &        FREP(I2),' x ',FREP(I1),' : ',
     &        REP(IFDIRP(1,J)),' + ',
     &        REP(IFDIRP(2,J)),' + ',
     &        REP(IFDIRP(3,J)),' + ',
     &        REP(IFDIRP(4,J))
  120     CONTINUE
  110   CONTINUE
        CALL HEADER('Spinor structure',-1)
        WRITE(LUPRI,'(/2(3X,A,I2,9X)/)')
     &     ('* Fermion irrep no.:',I,I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('La','|',REP(JSPINR(1,1,I)),
     &            '(',JQBAS(JSPINR(1,1,I),1),')',
     &               REP(JSPINR(2,1,I)),
     &            '(',JQBAS(JSPINR(2,1,I),1),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('Sa','|',REP(JSPINR(1,2,I)),
     &            '(',JQBAS(JSPINR(1,2,I),2),')',
     &               REP(JSPINR(2,2,I)),
     &            '(',JQBAS(JSPINR(2,2,I),2),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('Lb','|',REP(JSPINR(3,1,I)),
     &            '(',JQBAS(JSPINR(3,1,I),1),')',
     &               REP(JSPINR(4,1,I)),
     &            '(',JQBAS(JSPINR(4,1,I),1),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &   ('Sb','|',REP(JSPINR(3,2,I)),
     &            '(',JQBAS(JSPINR(3,2,I),2),')',
     &             REP(JSPINR(4,2,I)),
     &            '(',JQBAS(JSPINR(4,2,I),2),')','|',I=1,NFSYM)
        CALL HEADER('Quaternion symmetries',-1)
        WRITE(LUPRI,'(4X,A3,2X,A4)') 'Rep','T(+)'
        CALL PRSYMB(LUPRI,'-',29,4)
        DO IBRP = 0,NBSYM-1
          WRITE(LUPRI,'(4X,A3,4(2X,A1))')
     &       REP(IBRP),(QUNIT(IPQTOQ(IZ,IBRP)),IZ=1,NZ)
        ENDDO
        CALL FLSHFO(LUPRI)
      ENDIF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mopdef */
      SUBROUTINE MOPDEF(IXYZ,IXQ,IYQ,IZQ)
C***********************************************************************
C
C     Define generic one-electron operators
C
C***********************************************************************
      use maxorb_mod, only: MXCORB
      use symmet_mod, only: ISYMAX
      use mxcent_mod, only: MXCENT
      use maxaqn_mod, only: MXAQN, MXQN
!      use dgroup_mod, only: JCOM, JM4, JM4POS, JM4REP, JM4TIM, MCMP
#include "implicit.h"
C
C
#include "dgroup.h"
C
C     ********************************************
C     *****   Define 4-component matrices   ******
C     ********************************************
C
C     Time symmetric matrices:
C        I_4, gamma_5, beta, beta x gamma_5
C     Time anti-symmetric matrices:
C        Sigma, alpha, beta x Sigma, beta x alpha
C
C     NOTE: An imaginary phase is inserted if necessary
C           to make all matrices time symmetric to accomodate
C           the quaternion symmetry scheme of DIRAC.
C           JM4TIM reports the time reversal symmetry without imaginary phase
C
C     The matrices are given in terms of their quaternion forms.
C     The forms are obtained by a reordering of the 4-component matrices
C     corresponding to
C
C     | La |     | La |
C     | Lb | --> | Sa |
C     | Sa |     | Lb |
C     | Sb |     | Sb |
C
C     Time symmetric matrices then attain the structure:
C
C     |   A    B  |
C     |  -B*   A* |
C
C     with the quaternion form given by
C
C     A + Bj
C
C     The end result is:
C
C     -------------------------------------------------------------------------------------------
C     |            | 1              | i                 | j                 | k                 |
C     -------------------------------------------------------------------------------------------
C     | I_2        | I_4            | i(Sigma_z)        | i(Sigma_y)        | i(Sigma_x)        |
C     | sigma_z    | beta           | i(beta x Sigma_z) | i(beta x Sigma_y) | i(beta x Sigma_x) |
C     | i(sigma_y) | beta x gamma_5 | i(beta x alpha_z) | i(beta x alpha_y) | i(beta x alpha_x) |
C     | sigma_x    | gamma_5        | i(alpha_z)        | i(alpha_y)        | i(alpha_x)        |
C     -------------------------------------------------------------------------------------------
C
C     The spatial symmetry of these basis matrices is given in terms of the totally symmetric irrep Gamma_0,
C     the rotations R_q, the coordinates q and the function xyz:
C
C     Gamma_0: I_4       , beta
C     R_z    : i(Sigma_z), i(beta x Sigma_z)
C     R_y    : i(Sigma_y), i(beta x Sigma_y)
C     R_x    : i(Sigma_x), i(beta x Sigma_x)
C     xyz    : gamma_5   , beta x gamma_5
C     z      : i(alpha_z), i(beta x alpha_z)
C     y      : i(alpha_y), i(beta x alpha_y)
C     x      : i(alpha_x), i(beta x alpha_x)
C
C     This information is contained in JM4REP.
C
C     Commutation and anticommutation relations:
C     beta commutes with Sigma_i and Sigma_i commutes with gamma_5
C     beta anticommutes with gamma_5
C     alpha_i = gamma_5 x Sigma_i
C
C     For each basis matrix we define
C       JM4REP: spatial symmetry
C       JM4POS: position in quaternion matrix (1 -> 1, 2 -> i, 3 -> j, 4 -> k)
C       M4COMB: character string (LL,SL,LS,SS) giving the overall phase for each block
C
C     0. I_4 --> (1)x(I_2)
C     ====================
C
      JM4REP(0)  = 0
      JM4POS(0)  = 1
      JM4TIM(0)  = 1
C     M4COMB(0)  = '+00+'
C
C     1. i(alpha_z) --> (i)x(sigma_x)
C     ===============================
C
      JM4REP(1)  = ISYMAX(3,1)
      JM4POS(1)  = IZQ
      JM4TIM(1)  = -1
C     M4COMB(1)  = '0++0'
C
C     2. i(alpha_y) --> (j)x(sigma_x)
C     ===============================
C
      JM4REP(2)  = ISYMAX(2,1)
      JM4POS(2)  = IYQ
      JM4TIM(2)  = -1
C     M4COMB(2)  = '0++0'
C
C     3. i(alpha_x) --> (k)x(sigma_x)
C     ===============================
C
      JM4REP(3)  = ISYMAX(1,1)
      JM4POS(3)  = IXQ
      JM4TIM(3)  = -1
C     M4COMB(3)  = '0++0'
C
C     4. gamma_5 --> (1)x(sigma_x)
C     ============================
C
      JM4REP(4)  = IXYZ
      JM4POS(4)  = 1
      JM4TIM(4)  = 1
C     M4COMB(4)  = '0++0'
C
C     5. i(Sigma_z) --> (i)x(I_2)
C     ===========================
C
      JM4REP(5)  = ISYMAX(3,2)
      JM4POS(5)  = IZQ
      JM4TIM(5)  = -1
C     M4COMB(5)  = '+00+'
C
C     6. i(Sigma_y) --> (j)x(I_2)
C     ===========================
C
      JM4REP(6)  = ISYMAX(2,2)
      JM4POS(6)  = IYQ
      JM4TIM(6)  = -1
C     M4COMB(6)  = '+00+'
C
C     7. i(Sigma_x) --> (k)x(I_2)
C     ===========================
C
      JM4REP(7)  = ISYMAX(1,2)
      JM4POS(7)  = IXQ
      JM4TIM(7)  = -1
C     M4COMB(7)  = '+00+'
C
C     8. beta --> (1)x(sigma_z)
C     =========================
C
      JM4REP(8)  = 0
      JM4POS(8)  = 1
      JM4TIM(8)  = 1
C     M4COMB(8)  = '+00-'
C
C     9. beta x gamma_5 --> (1)x(i.sigma_y)
C     =====================================
C
      JM4REP(9)  = IXYZ
      JM4POS(9)  = 1
      JM4TIM(9)  = 1
C     M4COMB(9)  = '0-+0'
C
C     10. i(beta x Sigma_z) --> (i)x(sigma_z)
C     =======================================
C
      JM4REP(10) = ISYMAX(3,2)
      JM4POS(10) = IZQ
      JM4TIM(10) = -1
C     M4COMB(10) = '+00-'
C
C     11. i(beta x Sigma_y) --> (j)x(sigma_z)
C     =======================================
C
      JM4REP(11) = ISYMAX(2,2)
      JM4POS(11) = IYQ
      JM4TIM(11) = -1
C     M4COMB(11) = '+00-'
C
C     12. i(beta x Sigma_x) --> (k)x(sigma_z)
C     =======================================
C
      JM4REP(12) = ISYMAX(1,2)
      JM4POS(12) = IXQ
      JM4TIM(12) = -1
C     M4COMB(12) = '+00-'
C
C     13. i(beta x alpha_z) --> (i)x(i.sigma_y)
C     =========================================
C
      JM4REP(13) = ISYMAX(3,1)
      JM4POS(13) = IZQ
      JM4TIM(13) = -1
C     M4COMB(13) = '0-+0'
C
C     14. i(beta x alpha_y) --> (j)x(i.sigma_y)
C     =========================================
C
      JM4REP(14) = ISYMAX(2,1)
      JM4POS(14) = IYQ
      JM4TIM(14) = -1
C     M4COMB(14) = '0-+0'
C
C     15. i(beta x alpha_x) --> (k)x(i.sigma_y)
C     =========================================
C
      JM4REP(15) = ISYMAX(1,1)
      JM4POS(15) = IXQ
      JM4TIM(15) = -1
C     M4COMB(15) = '0-+0'
C
C     16. (beta - I_4)/2
C     ==================
C
      JM4REP(16) = 0
      JM4POS(16) = 1
      JM4TIM(16) = 1
C     M4COMB(16) = '000-'
C
C     Operators 1 to 16 are all time-reversal symmetric operators.
C     Only operator 17 is a time-reversal anti-symmetric operator.
C     Operator 17 is equal to operator 9, but taken as T antisymmetric.
C
C     17. beta x gamma_5 (T anti-symmetric)
C     =====================================
C
      JM4REP(17) = IXYZ
      JM4POS(17) = 1
      JM4TIM(17) = -1
C     M4COMB(17) = '0-+0'
C
C
C     ***************************************
C     ***** Definition of full operator *****
C     ***************************************
C
C     MCMP : number of components
C     JM4  : pointer to matrix operator for a given component
C     JCOM : coefficient of component
C
C     1. P                 * diagonal times scalar operator
C     =====================================================
C
      MCMP(1)    =  1
      JM4 (1,1)  =  0
      JCOM(1,1)  =  1
C
C     2. i[alpha_x]P       * x-component of alpha times scalar operator
C     =================================================================
C
      MCMP(2)    =  1
      JM4 (1,2)  =  3
      JCOM(1,2)  =  1
C
C     3. i[alpha_y]P       * y-component of alpha times scalar operator
C     =================================================================
C
      MCMP(3)    =  1
      JM4 (1,3)  =  2
      JCOM(1,3)  =  1
C
C     4. i[alpha_z]P       * z-component of alpha times scalar operator
C     =================================================================
C
      MCMP(4)    =  1
      JM4 (1,4)  =  1
      JCOM(1,4)  =  1
C
C     5. i[alpha x P]_x    * x-component of the vector product of alpha and vector operator
C     =====================================================================================
C
      MCMP(5)    =  2
      JM4 (1,5)  =  2
      JM4 (2,5)  =  1
      JCOM(1,5)  =  1
      JCOM(2,5)  = -1
C
C     6. i[alpha x P]_y    * y-component of the vector product of alpha and vector operator
C     =====================================================================================
C
      MCMP(6)    =  2
      JM4 (1,6)  =  1
      JM4 (2,6)  =  3
      JCOM(1,6)  =  1
      JCOM(2,6)  = -1
C
C     7. i[alpha x P]_z    * z-component of the vector product of alpha and vector operator
C     =====================================================================================
C
      MCMP(7)    =  2
      JM4 (1,7)  =  3
      JM4 (2,7)  =  2
      JCOM(1,7)  =  1
      JCOM(2,7)  = -1
C
C     8. i[alpha.P]        * dot-product of alpha and vector operator
C     ===============================================================
C
      MCMP(8)    =  3
      JM4 (1,8)  =  3
      JM4 (2,8)  =  2
      JM4 (3,8)  =  1
      JCOM(1,8)  =  1
      JCOM(2,8)  =  1
      JCOM(3,8)  =  1
C
C     9. [gamma5]P         * gamma5 times scalar operator
C     ===================================================
C
      MCMP(9)    =  1
      JM4 (1,9)  =  4
      JCOM(1,9)  =  1
C
C     10. i[Sigma_x]P      * x-component of Sigma times scalar operator
C     =================================================================
C
      MCMP(10)   =  1
      JM4 (1,10) =  7
      JCOM(1,10) =  1
C
C     11. i[Sigma_y]P      * y-component of Sigma times scalar operator
C     =================================================================
C
      MCMP(11)   =  1
      JM4 (1,11) =  6
      JCOM(1,11) =  1
C
C     12. i[Sigma_z]P      * z-component of Sigma times scalar operator
C     =================================================================
C
      MCMP(12)   =  1
      JM4 (1,12) =  5
      JCOM(1,12) =  1
C
C     13. i[beta Sigma_x]P * x-component of beta Sigma times scalar operator
C     ======================================================================
C
      MCMP(13)   =  1
      JM4 (1,13) =  12
      JCOM(1,13) =  1
C
C     14. i[beta Sigma_y]P * y-component of beta Sigma times scalar operator
C     ======================================================================
C
      MCMP(14)   =  1
      JM4 (1,14) =  11
      JCOM(1,14) =  1
C
C     15. i[beta Sigma_z]P * z-component of beta Sigma times scalar operator
C     ======================================================================
C
      MCMP(15)   =  1
      JM4 (1,15) =  10
      JCOM(1,15) =  1
C
C     16. i[beta alpha_x]P * x-component of beta alpha times scalar operator
C     ======================================================================
C
      MCMP(16)   =  1
      JM4 (1,16) =  15
      JCOM(1,16) =  1
C
C     17. i[beta alpha_y]P * y-component of beta alpha times scalar operator
C     ======================================================================
C
      MCMP(17)   =  1
      JM4 (1,17) =  14
      JCOM(1,17) =  1
C
C     18. i[beta alpha_z]P * z-component of beta alpha times scalar operator
C     ======================================================================
C
      MCMP(18)   =  1
      JM4 (1,18) =  13
      JCOM(1,18) =  1
C
C     19. [beta]P          * beta times scalar operator
C     =================================================
C
      MCMP(19)   =  1
      JM4 (1,19) =  8
      JCOM(1,19) =  1
C
C     20. i[Sigma.P]       * dot-product of Sigma and vector operator
C     ===============================================================
C
      MCMP(20)   =  3
      JM4 (1,20) =  7
      JM4 (2,20) =  6
      JM4 (3,20) =  5
      JCOM(1,20) =  1
      JCOM(2,20) =  1
      JCOM(3,20) =  1
C
C     21. [beta gamma5]P   * beta gamma5 times scalar operator
C     ========================================================
C
      MCMP(21)   =  1
      JM4 (1,21) =  9
      JCOM(1,21) =  1
C
C     22. i[beta alpha.P]  * dot-product of beta alpha and vector operator
C     ====================================================================
C
      MCMP(22)   =  3
      JM4 (1,22) =  15
      JM4 (2,22) =  14
      JM4 (3,22) =  13
      JCOM(1,22) =  1
      JCOM(2,22) =  1
      JCOM(3,22) =  1
C
C     23. i[beta gamma5]P  * i beta gamma5 times scalar operator
C     ==========================================================
C
      MCMP(23)   =  1
      JM4 (1,23) =  17
      JCOM(1,23) =  1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION IQFROMPQ(IBSYM,IFSYM,IC,IZ)
C***********************************************************************
C
C     Give quaternion unit corresponding to symmetry-packed one
C     IBSYM is index of boson irreps of fermion ircop IFSYM
C     Written by T. Saue April 10 2013
C
C***********************************************************************
      use symmet_mod, only: ISYMAX
      use maxaqn_mod, only: MXAQN, MXQN
!      use dcbibt_mod, only: IBTAND, IBTSHR
      use maxorb_mod, only: MXCORB
      use mxcent_mod, only: MXCENT
!      use dgroup_mod, only: PT
#include "implicit.h"
C
      INTEGER :: I,III,J,K
#include "dgroup.h"
      DIMENSION III(2)
#include "dcbibt.h"
C.....III is the totally symmetric irrep for the large components
      III(1) = 0
C.....III is the symmety of XYZ for the small components
      III(2) = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C.....find the boson irrep that corresponds to index IBSYM of IFSYM for component IC
      ITMP     = MOD(IFSYM+IC,2) + 1
      ITMP     = JFSYM(IBSYM,ITMP) - 1
      ITMP     = IBTXOR(ITMP,III(IC))
      ITMP     = IPQTOQ(IZ,ITMP)
      IQFROMPQ = IQMAP(ITMP)
      END
C --- end of hergroup.F ---
