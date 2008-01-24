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
C     /* GET_DENSITY */
      SUBROUTINE GET_DENSITY(IOPTION,
     &                       SOME_VECTOR,
     &                       DENSITY,
     &                       IREPOP,
     &                       GAO,
     &                       GAO1,
     &                       DMAT,
     &                       BUF)
C=======================================================================
C
C     IOPTION is a four digit integer:
C
C     1 get   1 "nothing"     1 "nothing"   1 1_{4x4}
C     2 put   2 \nabla        2 \cdot       2 c \alpha
C             3 some vector   3 \times      3 \Sigma
C                                           4 \beta\Sigma
C                                           5 \gamma_5
C
C     in several places (grep for GAO1)
C     hermitian or antihermitian DMAT is assumed
C     these places will fail in case of a general DMAT
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

      LOGICAL   LDEC, LGET
      DIMENSION SOME_VECTOR(*),
     &          DENSITY(*),
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

      LGET = LDEC(IOPTION,1,1)

      K = 1
      L = 1
      M = 2
      N = 2
      A = D1
      B = D1
      ITRANSFORM = 0

      IF(LDEC(IOPTION,4,2) .OR.
     &   LDEC(IOPTION,4,5)) THEN
        K = 1
        L = 2
        M = 2
        N = 1
      ENDIF

      IF(LDEC(IOPTION,4,2)) THEN
C.......using current value of c
        A = -CVAL
        B = -CVAL
      ENDIF
      IF(LDEC(IOPTION,4,4)) B = -D1

      IF(LDEC(IOPTION,4,3) .OR.
     &   LDEC(IOPTION,4,4)) THEN
        ITRANSFORM = 2
      ENDIF
      IF(LDEC(IOPTION,4,2)) ITRANSFORM = 1

!     no funny stuff

      IF(LDEC(IOPTION,3,1)) THEN

        IEND = 3
        IF(LDEC(IOPTION,4,1)) IEND = 1

        DO IXYZ = 1,IEND
          IF(LDEC(IOPTION,2,1)) THEN
            CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                           DENSITY(IXYZ),
     &                           A,B,K,L,M,N,
     &                           IXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)
          ELSEIF(LDEC(IOPTION,2,2)) THEN
            DO JXYZ = 1,3
              CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                             DENSITY(JXYZ+3*(IXYZ-1)),
     &                             D2*A,D2*B,K,L,M,N,
     &                             IXYZ,ITRANSFORM,IREPOP,
     &                             GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                             DMAT,BUF)
            ENDDO
          ENDIF
        ENDDO

!     \cdot

      ELSEIF(LDEC(IOPTION,3,2)) THEN

        IF(LGET) DENSITY(1) = D0
        DO JXYZ = 1,3
          IF(LDEC(IOPTION,2,2)) THEN
            CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                           DENSITY(1),
     &                           D2*A,D2*B,K,L,M,N,
     &                           JXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                           DMAT,BUF)
          ELSEIF(LDEC(IOPTION,2,3)) THEN
            CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                           DENSITY(1),
     &                           SOME_VECTOR(JXYZ)*A,
     &                           SOME_VECTOR(JXYZ)*B,
     &                           K,L,M,N,
     &                           JXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)
          ENDIF
        ENDDO

!     \times

      ELSEIF(LDEC(IOPTION,3,3)) THEN

        DO IXYZ = 1,3

          JXYZ = IXYZ + 1
          KXYZ = JXYZ + 1

          IF(JXYZ .GT. 3) JXYZ = JXYZ - 3
          IF(KXYZ .GT. 3) KXYZ = KXYZ - 3

          IF(LDEC(IOPTION,2,2)) THEN
            CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                           DENSITY(IXYZ),
     &                           D2*A,D2*B,K,L,M,N,
     &                           JXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO1(1+JOFF*(KXYZ-1)),
     &                           DMAT,BUF)
            CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                           DENSITY(IXYZ),
     &                           -D2*A,-D2*B,K,L,M,N,
     &                           KXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO1(1+JOFF*(JXYZ-1)),
     &                           DMAT,BUF)
          ELSEIF(LDEC(IOPTION,2,3)) THEN
            CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                           DENSITY(IXYZ),
     &                           SOME_VECTOR(KXYZ)*A,
     &                           SOME_VECTOR(KXYZ)*B,
     &                           K,L,M,N,
     &                           JXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)
            CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                           DENSITY(IXYZ),
     &                           -SOME_VECTOR(JXYZ)*A,
     &                           -SOME_VECTOR(JXYZ)*B,
     &                           K,L,M,N,
     &                           KXYZ,ITRANSFORM,IREPOP,
     &                           GAO,GAO,
     &                           DMAT,BUF)
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DENSITY_SUB */
      SUBROUTINE GET_DENSITY_SUB(LGET,
     &                           RESET,
     &                           DENSITY,
     &                           A,
     &                           B,
     &                           K,
     &                           L,
     &                           M,
     &                           N,
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

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
#endif
      LOGICAL   LGET,
     &          RESET
      DIMENSION GAO_LEFT(*),
     &          GAO_RIGHT(*),
     &          DMAT(*),
     &          BUF(*)

#include "ibtfun.h"

#if defined (VAR_G77) || defined (SYS_SX)
      ioff     = nbast*nbast
      joff_max = nbast*nbast*3
#else
      ioff     = ntbas(0)*ntbas(0)
      joff_max = ntbas(0)*ntbas(0)*3
#endif

      if(lget .and. reset) density = d0

      intermediate = irepop
      joff     = 0
      joff_min = 0
      if(.not. itransform .eq. 0) then
        intermediate = ibtxor(isymax(ixyz,itransform),irepop)
        iq = JM4POS(4-ixyz)
        iq = iqtopq(iq,isymax(ixyz,2))
C        iq = iqtopq(5-ixyz,isymax(ixyz,2))
        joff = ioff*(iq-1)
      endif
      if(joff .lt. joff_min) then
        WRITE(6,*) 'INFO..',irepop,ixyz,jm4pos(5-ixyz),iq,joff
        call quit('joff .lt. joff_min')
      endif
      if(joff .gt. joff_max) call quit('joff .gt. joff_max')

      ICHECK = K + L
      IF((MC .EQ. 1) .AND. (ICHECK .GT. 2)) GOTO 1

      CALL GET_DENSITY_SUBSUB(LGET,
     &                        DENSITY,
     &                        A,K,L,
     &                        INTERMEDIATE,
     &                        GAO_LEFT,
     &                        GAO_RIGHT,
     &                        DMAT(1+JOFF),
     &                        BUF)

1     CONTINUE

      ICHECK = M + N
      IF((MC .EQ. 1) .AND. (ICHECK .GT. 2)) GOTO 2

      CALL GET_DENSITY_SUBSUB(LGET,
     &                        DENSITY,
     &                        B,M,N,
     &                        INTERMEDIATE,
     &                        GAO_LEFT,
     &                        GAO_RIGHT,
     &                        DMAT(1+JOFF),
     &                        BUF)

2     CONTINUE

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DENSITY_SUBSUB */
      SUBROUTINE GET_DENSITY_SUBSUB(LGET,
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

      LOGICAL   LGET
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

            IF(LGET) THEN

              I1 = IBBAS(IREP1,IC1) + 1
              I2 = IBBAS(IREP2,IC2) + 1

              CALL DGEMV('N',NBBAS(IREP1,IC1),NBBAS(IREP2,IC2),D1,
     &                   DMAT(I1,I2),NTBAS(0),GAO_RIGHT(I2),1,D0,BUF,1)
              DENSITY = DENSITY
     &                + F*DDOT(NBBAS(IREP1,IC1),GAO_LEFT(I1),1,BUF,1)

            ELSE !DMAT is then EXCMAT

              I1 = IBBAS(IREP1,IC1)
              I2 = IBBAS(IREP2,IC2)

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
