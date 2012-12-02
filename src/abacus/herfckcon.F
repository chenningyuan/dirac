!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

#if defined (VAR_SINGLE)
#include "single.h"
#endif
C
#ifdef MOD_ERI
C  /* Deck fckcon */
      SUBROUTINE FCKCON(FMAT,DMAT,NDMAT,I,AOINT,IND,
     &        NINT,NIND,NBUF,IX,IY1,FAC)
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
C  This routine tries to replace INTFC1(her2fck.F), FOKDI1(eri2fck.F)
C  and FCKDS1, FCKDS2,FCKDT3 (sirfck.F)
C
C We keep to outer loop over NDMAT(in calling routine), but we vectorize with
C IVECLN / NDMAT temoprary matrices. This will use the same 
C amount of scratch memory we use in fokdi1 in eri. 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dcbham.h"
      PARAMETER (D0=0.0D0,DP25=0.25D0,DP5=0.5D0)
      INTEGER A,B,C,D
      PARAMETER (K = 1)
      DIMENSION FMAT(NBASIS,NBASIS,NDMAT,K)
      DIMENSION DMAT(NBASIS,NBASIS,NDMAT), AOINT(NINT),IND(NIND,4)
C
C
      IF (IX.NE.2) THEN
        IY = IY1
      ELSE
C       ... No Coulomb contribution for
C           antisymmetric density matrix
        IY = IY1 - MOD(IY1,2)
      END IF
C
      IF (HFXFAC.EQ.D0 .AND. HFXMU.EQ.D0) THEN
        IY = MOD(IY,2)
C       ... only Coulomb contribution
      ELSE
         DP5X  =  DP5*HFXFAC
         DP25X = DP25*HFXFAC
      END IF
      IF (HFXMU.NE.D0) IY = 2 
      IF (IY .EQ. 0) GO TO 9999
C     ... nothing to do!
C
        IF    (IX.EQ.1.OR.IX.EQ.2) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO J = 1,NBUF
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC   
                EINT = -DP25X*DINT       
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I) 
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I) 
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I) 
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I) 
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I) 
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I) 
              END DO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO J = 1,NBUF
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                EINT = -DP25X*DINT
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              END DO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO J = 1,NBUF
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I)
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I)
              END DO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'FCKCON ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IX*10+IY,
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in FCKCON: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO J = 1,NBUF
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
             ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO J = 1,NBUF
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
            END DO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO J = 1,NBUF
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
            END DO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'FCKCON ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IX*10+IY1,
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in FCKCON: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      ' ERROR, specified IFCTYP not implemented yet',
     &      '        specified IFCTYP was',IX*10+IY1,
     &      '        for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in FCKCON: specified IFCTYP not implemented.')
        ENDIF
C
 9999 CONTINUE
      RETURN
C
      END
C
#ifdef MOD_DNF
C  /* Deck df_fckcon */
      SUBROUTINE DF_FCKCON(FMAT,DMAT,DF_FMAT,DF_DMAT,NDMAT,I,AOINT,IND,
     &                     NINT,NIND,NBUF,IX,IY1,FAC,MODE)
C
C
C       MODE 1 : Contract integrals with a fitted density matrix and produce
C                ordinary Fock matrix
C
C       MODE 2 : Contract integrals with regular density matrix and produce
C                fit function weight vector
C
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
C  This routine is a trimmed version of FCKCON
C
C We keep to outer loop over NDMAT(in calling routine), but we vectorize with
C IVECLN / NDMAT temoprary matrices. This will use the same 
C amount of scratch memory we use in fokdi1 in eri. 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dcbham.h"
      PARAMETER (D0=0.0D0,DP25=0.25D0,DP5=0.5D0)
      INTEGER A,B,C,D
      DIMENSION FMAT(NBASIS,NBASIS,NDMAT)
      DIMENSION DMAT(NBASIS,NBASIS,NDMAT), AOINT(NINT),IND(NIND,4)
      DIMENSION DF_DMAT(NBASISAUX,NDMAT),DF_FMAT(NBASISAUX,NDMAT)
      DIMENSION ICLASSES(4)
C
C
      IF (IX.NE.2) THEN
        IY = IY1
      ELSE
C       ... No Coulomb contribution for
C           antisymmetric density matrix
        IY = IY1 - MOD(IY1,2)
      END IF
C
C     Check whether this is a pure DFT calculation (no exchange)
C
      IF (HFXFAC.EQ.D0 .AND. HFXMU.EQ.D0) THEN
        IY = MOD(IY,2)
C       ... only Coulomb contribution
      ELSE
         DP5X  =  DP5*HFXFAC
         DP25X = DP25*HFXFAC
      END IF
      IF (HFXMU.NE.D0) IY = 2 
C
      iy = 1
      IF (IY .EQ. 0) GO TO 9999
C     ... nothing to do!
C
      IF (MODE.EQ.1) THEN
C       We will process integrals of the type (ab|nf)
C       where |a> and |b> are ordinary functions,
C       |f> is a fit function and |n> the null function
        ICLASSES(1) = 1
        ICLASSES(2) = 1
        ICLASSES(3) = 0
        ICLASSES(4) = 3
        CALL INDEX_CANON(IND(1,3),NIND,NBUF)
        CALL DF_INDEX_SHIFT(IND,NIND,IND,NIND,NBUF,ICLASSES,4,1)
        IF  (IY.NE.1) THEN
C
C         Coulomb contributions only
C         ==========================
C
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DF_FCKCON ERROR, specified IFCTYP not possible here',
     &      '                 specified IFCTYP was',IX*10+IY,
     &      '                 for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DF_FCKCON: can only do Coulomb type.')
        ELSEIF  (IX.EQ.1.OR.IX.EQ.2) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
            DO J = 1,NBUF
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                FMAT(A,B,I) = FMAT(A,B,I) + DINT*DF_DMAT(D,I)
              END DO
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
            DO J = 1,NBUF
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DF_DMAT(D,I)
              FMAT(B,A,I) = FMAT(B,A,I) + DINT*DF_DMAT(D,I)
            END DO
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      ' ERROR, specified IFCTYP not implemented yet',
     &      '        specified IFCTYP was',IX*10+IY1,
     &      '        for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DF_FCKCON: specified IFCTYP not implemented.')
        ENDIF
      ELSEIF (MODE.EQ.2) THEN
C       We will process integrals of the type (nf|ij)
C       where |c> and |d> are ordinary functions,
C       |f> is a fit function and |n> the null function
        ICLASSES(1) = 0
        ICLASSES(2) = 3
        ICLASSES(3) = 1
        ICLASSES(4) = 1
        CALL INDEX_CANON(IND(1,1),NIND,NBUF)
        CALL DF_INDEX_SHIFT(IND,NIND,IND,NIND,NBUF,ICLASSES,4,1)
        IF  (IY.NE.1) THEN
C
C         Coulomb contributions only
C         ==========================
C
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DF_FCKCON ERROR, specified IFCTYP not possible here',
     &      '                 specified IFCTYP was',IX*10+IY,
     &      '                 for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DF_FCKCON: can only do Coulomb type.')
        ELSEIF  (IX.EQ.1.OR.IX.EQ.2) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
            DO J = 1,NBUF
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                DF_FMAT(B,I) = DF_FMAT(B,I) + DINT*DMAT(C,D,I)
              END DO
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
            DO J = 1,NBUF
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              DF_FMAT(B,I) = DF_FMAT(B,I) + GCD
            END DO
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      ' ERROR, specified IFCTYP not implemented yet',
     &      '        specified IFCTYP was',IX*10+IY1,
     &      '        for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DF_FCKCON: specified IFCTYP not implemented.')
        ENDIF
      ENDIF
C
 9999 CONTINUE
      RETURN
C
      END
#endif /* MOD_DNF */
#else /* MOD_ERI */
C dummy routine to keep some compilers happy 
      SUBROUTINE FCKCON(FMAT,DMAT,NDMAT,I,AOINT,IND,NINT,NIND,NBUF,
     &                  IX,IY1,FAC)    
      CALL QUIT('*** FCKCON needs ERI to be activated ***')
      RETURN 
      END
#endif /* MOD_ERI */
