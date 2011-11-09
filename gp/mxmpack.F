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

C
C
C  /* Deck mxm */
#if !defined (SYS_CRAY) && !defined (VAR_MOTECC)
      SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
#include "implicit.h"
      DIMENSION A(NAR,NAC), B(NAC,NBC), C(NAR,NBC)
#if defined (VAR_F90)
C
C     tuh Mar 1988 - Alliant spends a lot of time setting up this call,
C     it might therefore be best to avoid calling MATMUL
C
      C=MATMUL(A,B)
#else
#if defined (VAR_ESSL)
C
C ESSL DGEMUL uses same cpu time as DGEMM /930330-hjaaj
C     but has no problem with NaN /931022-hjaaj
      CALL DGEMUL(A,NAR,'N',B,NAC,'N',C,NAR,NAR,NAC,NBC)
#endif
#endif
#if !defined (VAR_ESSL) && !defined (VAR_F90)
#if defined (SYS_CONVEX) || defined (SYS_IRIX) || defined (SYS_AIX) || defined (VAR_BLAS3) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX) || defined (SYS_SX) || defined (SYS_DARWIN)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C     Necessary to zero out because of NaN, tuh
C
      CALL DZERO(C,NAR*NBC)
      CALL DGEMM('N','N',NAR,NBC,NAC,D1,A,NAR,B,NAC,D0,C,NAR)
#else
C
C     C(I,J) = A(I,K) * B(K,J)
C
      DO 100 J = 1,NBC
         DO 100 I = 1,NAR
  100       C(I,J) = A(I,1)*B(1,J)
C
      DO 300 K = 2,NAC
         DO 300 J = 1,NBC
            BKJ = B(K,J)
            DO 300 I = 1,NAR
  300          C(I,J) = C(I,J) + A(I,K) * BKJ
C
#endif
#endif
      RETURN
      END
#endif
C  /* Deck mxmd */
#if defined (VAR_NOTUSED  )
C930522-hjaaj: MXMD not used in Hermit/Abacus any more
#if !defined (VAR_MOTECC)
#if defined (VAR_F90)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C
C     tuh Mar 1988 - Alliant spends a lot of time setting up this call,
C     it might therefore be best to avoid calling MATMUL
C
#include "implicit.h"
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
      C(1:NARA,:)=MATMUL(A,B(1:NACA,:))
      RETURN
      END
#endif
#if defined (SYS_CONVEX) || defined (VAR_ESSL) || defined (VAR_BLAS3)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C     890905 tuh
#include "implicit.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
      CALL DZERO(C,NAR*NBC)
      CALL DGEMM('N','N',NARA,NBC,NACA,D1,A,NARA,B,NAC,D0,C,NAR)
      RETURN
      END
#endif
#if !defined (SYS_CONVEX) && !defined (VAR_ESSL) && !defined (VAR_F90) && !defined (VAR_BLAS3)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C
C     Modified MXM by Trygve U. Helgaker
C
#include "implicit.h"
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
C
C     C(1:NARA,:)=MATMUL(A,B(1:NACA,:))
C     C(I,J) = A(I,K) * B(K,J)
C
      DO 100 J = 1,NBC
         DO 100 I = 1,NARA
  100       C(I,J) = A(I,1)*B(1,J)
C
      DO 300 K = 2,NACA
         DO 300 J = 1,NBC
            BKJ = B(K,J)
            DO 300 I = 1,NARA
  300          C(I,J) = C(I,J) + A(I,K) * BKJ
C
      RETURN
      END
#endif
#endif
#else
C
C     This is only to avoid the object file being empty, causing problems
C     on Cray C90./jth 980716
C
      SUBROUTINE GNYF
      END
#endif
