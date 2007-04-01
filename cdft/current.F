C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_CURRENT */
      SUBROUTINE GET_CURRENT(DO_CURRENT,
     &                       DO_DERVS,
     &                       CURRENT,
     &                       CURRENT_DERVS,
     &                       F,
     &                       G,
     &                       IU,
     &                       IV,
     &                       IX,
     &                       IY,
     &                       ITRANSFORM,
     &                       IREPOP,
     &                       GAO,
     &                       GAO1,
     &                       DMAT,
     &                       BUF)
C=======================================================================
C
C     CURRENT_DERVS(1,1) = derivative of x-component wrt x
C     CURRENT_DERVS(2,1) = derivative of x-component wrt y
C     CURRENT_DERVS(3,1) = derivative of x-component wrt z
C
C     CURRENT_DERVS(1,2) = derivative of y-component wrt x
C     CURRENT_DERVS(2,2) = derivative of y-component wrt y
C     CURRENT_DERVS(3,2) = derivative of y-component wrt z
C
C     CURRENT_DERVS(1,3) = derivative of z-component wrt x
C     CURRENT_DERVS(2,3) = derivative of z-component wrt y
C     CURRENT_DERVS(3,3) = derivative of z-component wrt z
C
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "cdft.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
      LOGICAL   DO_CURRENT,
     &          DO_DERVS
      DIMENSION CURRENT(3),
     &          CURRENT_DERVS(3,3),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#include "ibtfun.h"

      IF(DO_CURRENT) THEN
        CALL DZERO(CURRENT,3)
      ENDIF
      IF(DO_DERVS) THEN
        CALL DZERO(CURRENT_DERVS,9)
      ENDIF

      DO IXYZ = 1,3
#if defined (VAR_G77) || defined (SYS_SX)
        IOFF = NBAST*NBAST
        JOFF = NBAST
#else
        IOFF = NTBAS(0)*NTBAS(0)
        JOFF = NTBAS(0)
#endif
        IQ = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))
        IOFF = IOFF*(IQ-1)

        INTERMEDIATE = IBTXOR(ISYMAX(IXYZ,ITRANSFORM),IREPOP)

        IF(DO_CURRENT) THEN
          CALL GET_CURRENT_SUB(CURRENT(IXYZ),
     &                         F,IU,IV,
     &                         GAO,GAO,
     &                         INTERMEDIATE,DMAT(1+IOFF),BUF)
          CALL GET_CURRENT_SUB(CURRENT(IXYZ),
     &                         G,IX,IY,
     &                         GAO,GAO,
     &                         INTERMEDIATE,DMAT(1+IOFF),BUF)
        ENDIF
        IF(DO_DERVS) THEN
          DO JXYZ = 1,3
            CALL GET_CURRENT_SUB(CURRENT_DERVS(JXYZ,IXYZ),
     &                           D2*F,IU,IV,
     &                           GAO,
     &                           GAO1(1+JOFF*(JXYZ-1)),
     &                           INTERMEDIATE,DMAT(1+IOFF),BUF)
            CALL GET_CURRENT_SUB(CURRENT_DERVS(JXYZ,IXYZ),
     &                           D2*G,IX,IY,
     &                           GAO,
     &                           GAO1(1+JOFF*(JXYZ-1)),
     &                           INTERMEDIATE,DMAT(1+IOFF),BUF)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* DIV_VECTOR */
      SUBROUTINE DIV_VECTOR(S,
     &                      A_DERVS)
C=======================================================================
C     s = \nabla \cdot a
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION A_DERVS(3,3)

      S = 0.0D0
      DO IXYZ = 1,3
        S = S + A_DERVS(IXYZ,IXYZ)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* ROT_VECTOR */
      SUBROUTINE ROT_VECTOR(V,
     &                      A_DERVS)
C=======================================================================
C     v = \nabla \times a
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION V(3),
     &          A_DERVS(3,3)

      V(1) = A_DERVS(2,3) - A_DERVS(3,2)
      V(2) = A_DERVS(3,1) - A_DERVS(1,3)
      V(3) = A_DERVS(1,2) - A_DERVS(2,1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VECTOR_CROSS_VECTOR */
      SUBROUTINE VECTOR_CROSS_VECTOR(V,
     &                               A,
     &                               B)
C=======================================================================
C     v = a \times b
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION V(3),
     &          A(3),
     &          B(3)

      V(1) = A(2)*B(3) - A(3)*B(2)
      V(2) = A(3)*B(1) - A(1)*B(3)
      V(3) = A(1)*B(2) - A(2)*B(1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_V */
      SUBROUTINE GET_V(VECTOR_V,
     &                 VECTOR_J,
     &                 VECTOR_J_DERVS,
     &                 DENSITY,
     &                 DENSITY_GRADIENT)
C=======================================================================
C
C                                            j(r)
C     returns vorticity v(r) = \nabla \times ----
C                                            n(r)
C
C     needs j, derivatives of j, n, and grad(n)
C
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: march 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION VECTOR_V(3),
     &          VECTOR_TEMP(3),
     &          VECTOR_J(3),
     &          VECTOR_J_DERVS(3,3),
     &          DENSITY_GRADIENT(3)

      DO IXYZ = 1,3
        VECTOR_TEMP(IXYZ) = DENSITY_GRADIENT(IXYZ)/DENSITY
      ENDDO

      CALL VECTOR_CROSS_VECTOR(VECTOR_V,VECTOR_TEMP,VECTOR_J)

      CALL ROT_VECTOR(VECTOR_TEMP,VECTOR_J_DERVS)

      DO IXYZ = 1,3
        VECTOR_V(IXYZ) = (VECTOR_TEMP(IXYZ) - VECTOR_V(IXYZ))/DENSITY
      ENDDO

      RETURN
      END

















C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DENSITY */
      SUBROUTINE GET_DENSITY(IGET,
     &                       DENSITY,
     &                       IOPERATOR_LEFT,
     &                       ICOMBINE,
     &                       IOPERATOR_RIGHT,
     &                       IREPOP,
     &                       GAO,
     &                       GAO1,
     &                       DMAT,
     &                       BUF)
C=======================================================================
C
C     IOPERATOR_LEFT  = 1      "nothing"
C                       2      \nabla
C
C     ICOMBINE        = 1      "nothing"
C                       2      \cdot
C                       3      \times
C
C     IOPERATOR_RIGHT = 1      1_{4x4}
C                       2      c \alpha
C                       3      \Sigma
C                       4      \beta\Sigma
C                       5      \gamma_5
C
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"

      DIMENSION DENSITY(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      JOFF = NBAST
#else
      JOFF = NTBAS(0)
#endif

      K = 1
      L = 1
      M = 2
      N = 2
      A = D1
      B = D1
      ITRANSFORM = 0

      IF((IOPERATOR_RIGHT .EQ. 2) .OR. (IOPERATOR_RIGHT .EQ. 5)) THEN
        K = 1
        L = 2
        M = 2
        N = 1
      ENDIF

      IF(IOPERATOR_RIGHT .EQ. 4) THEN
        A = -137.0359998D0
        B = -137.0359998D0
      ENDIF
      IF(IOPERATOR_RIGHT .EQ. 4) B = -D1

      IF((IOPERATOR_RIGHT .EQ. 3) .OR. (IOPERATOR_RIGHT .EQ. 4)) THEN
        ITRANSFORM = 2
      ENDIF
      IF(IOPERATOR_RIGHT .EQ. 2) ITRANSFORM = 1

      IF(ICOMBINE .EQ. 1) THEN

        IEND = 3
        IF(IOPERATOR_RIGHT .EQ. 1) IEND = 1

        DO IXYZ = 1,IEND

          IF(IOPERATOR_LEFT .EQ. 1) THEN

            IF(IGET .EQ. 1) DENSITY(IXYZ) = D0
            CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                           A,K,L,
     &                           IXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)
            CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                           B,M,N,
     &                           IXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)

          ELSEIF(IOPERATOR_LEFT .EQ. 2) THEN

            DO JXYZ = 1,3

              IF(IGET .EQ. 1) DENSITY(JXYZ+3*(IXYZ-1)) = D0
              CALL GET_DENSITY_SUB(IGET,DENSITY(JXYZ+3*(IXYZ-1)),
     &                             D2*A,K,L,
     &                             IXYZ,ITRANSFORM,IREPOP,
     &                             GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                             DMAT,BUF)
              CALL GET_DENSITY_SUB(IGET,DENSITY(JXYZ+3*(IXYZ-1)),
     &                             D2*B,M,N,
     &                             IXYZ,ITRANSFORM,IREPOP,
     &                             GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                             DMAT,BUF)
            ENDDO
          ENDIF
        ENDDO

      ELSEIF(ICOMBINE .EQ. 2) THEN

        IF(IGET .EQ. 1) DENSITY(1) = D0
        DO JXYZ = 1,3

          CALL GET_DENSITY_SUB(IGET,DENSITY(1),
     &                         D2*A,K,L,
     &                         JXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                         DMAT,BUF)
          CALL GET_DENSITY_SUB(IGET,DENSITY(1),
     &                         D2*B,M,N,
     &                         JXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                         DMAT,BUF)
        ENDDO

      ELSEIF(ICOMBINE .EQ. 3) THEN

        DO IXYZ = 1,3

          IF(IGET .EQ. 1) DENSITY(IXYZ) = D0
   
          JXYZ = IXYZ + 1
          KXYZ = JXYZ + 1

          IF(JXYZ .GT. 3) JXYZ = JXYZ - 3
          IF(KXYZ .GT. 3) KXYZ = KXYZ - 3

          CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                         D2*A,K,L,
     &                         JXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(KXYZ-1)),
     &                         DMAT,BUF)
          CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                         D2*B,M,N,
     &                         JXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(KXYZ-1)),
     &                         DMAT,BUF)
          CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                         -D2*A,K,L,
     &                         KXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                         DMAT,BUF)
          CALL GET_DENSITY_SUB(IGET,DENSITY(IXYZ),
     &                         -D2*B,M,N,
     &                         KXYZ,ITRANSFORM,IREPOP,
     &                         GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                         DMAT,BUF)

        ENDDO

      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DENSITY_SUB */
      SUBROUTINE GET_DENSITY_SUB(IGET,
     &                           DENSITY,
     &                           F,
     &                           IC1,
     &                           IC2,
     &                           IXYZ,
     &                           ITRANSFORM,
     &                           IREPOP,
     &                           GAO_LEFT,
     &                           GAO_RIGHT,
     &                           DMAT,
     &                           BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcbgen.h"

      DIMENSION GAO_LEFT(*),
     &          GAO_RIGHT(*),
     &          DMAT(*),
     &          BUF(*)

#include "ibtfun.h"

      ICHECK = IC1 + IC2
      IF((MC .EQ. 1) .AND. (ICHECK .GT. 2)) GOTO 999

      IF(ITRANSFORM .EQ. 0) THEN

        CALL GET_DENSITY_SUBSUB(IGET,
     &                          DENSITY,
     &                          F,
     &                          IC1,
     &                          IC2,
     &                          IREPOP,
     &                          GAO_LEFT,
     &                          GAO_RIGHT,
     &                          DMAT,
     &                          BUF)

      ELSE

#if defined (VAR_G77) || defined (SYS_SX)
        IOFF = NBAST*NBAST
#else
        IOFF = NTBAS(0)*NTBAS(0)
#endif
        IQ   = IQTOPQ(5-IXYZ,ISYMAX(IXYZ,2))

        CALL GET_DENSITY_SUBSUB(IGET,
     &                          DENSITY,
     &                          F,
     &                          IC1,
     &                          IC2,
     &                          IBTXOR(ISYMAX(IXYZ,ITRANSFORM),IREPOP),
     &                          GAO_LEFT,
     &                          GAO_RIGHT,
     &                          DMAT(1+IOFF*(IQ-1)),
     &                          BUF)

      ENDIF

999   CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DENSITY_SUBSUB */
      SUBROUTINE GET_DENSITY_SUBSUB(IGET,
     &                              DENSITY,
     &                              F,
     &                              IC1,
     &                              IC2,
     &                              IREPOP,
     &                              GAO_LEFT,
     &                              GAO_RIGHT,
     &                              DMAT,
     &                              BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION DMAT(NBAST,NBAST)
#else
      DIMENSION DMAT(NTBAS(0),NTBAS(0))
#endif
      DIMENSION GAO_LEFT(*),
     &          GAO_RIGHT(*),
     &          BUF(*)

#include "ibtfun.h"

      DO IREP1 = 0,NBSYM-1

        N1 = NBBAS(IREP1,IC1)
        IF(N1 .GT. 0) THEN

          IREP2 = IBTXOR(IREP1,IREPOP)
          N2    = NBBAS(IREP2,IC2)
          IF(N2 .GT. 0) THEN

            I1 = IBBAS(IREP1,IC1) + 1
            I2 = IBBAS(IREP2,IC2) + 1

            IF(IGET .EQ. 1) THEN !get

              CALL DGEMV('N',NBBAS(IREP1,IC1),NBBAS(IREP2,IC2),D1,
     &                   DMAT(I1,I2),NTBAS(0),GAO_RIGHT(I2),1,D0,BUF,1)
              DENSITY = DENSITY
     &                + F*DDOT(NBBAS(IREP1,IC1),GAO_LEFT(I1),1,BUF,1)

            ELSE !put (DMAT is then EXCMAT)

              DO J = 1,N2
                FACTOR = F*DENSITY*GAO_RIGHT(I2+J)
                CALL DAXPY(N1,FACTOR,GAO_LEFT(I1+1),1,DMAT(I1+1,I2+J),1)
              ENDDO

            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
