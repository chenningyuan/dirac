!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.sdu.dk
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_ROTM */
      SUBROUTINE GET_ROTM(LGET,
     &                    VECTOR,
     &                    IREPOP,
     &                    GAO,
     &                    GAO1,
     &                    DMAT,
     &                    BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2233
      IF(LGET) IOPTION = 1233

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 VECTOR,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_ROTJ */
      SUBROUTINE GET_ROTJ(LGET,
     &                    VECTOR,
     &                    IREPOP,
     &                    GAO,
     &                    GAO1,
     &                    DMAT,
     &                    BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2232
      IF(LGET) IOPTION = 1232

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 VECTOR,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_M_DERVS */
      SUBROUTINE GET_M_DERVS(LGET,
     &                       VECTOR_DERVS,
     &                       IREPOP,
     &                       GAO,
     &                       GAO1,
     &                       DMAT,
     &                       BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION VECTOR_DERVS(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2213
      IF(LGET) IOPTION = 1213

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 VECTOR_DERVS,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_J_DERVS */
      SUBROUTINE GET_J_DERVS(LGET,
     &                       VECTOR_DERVS,
     &                       IREPOP,
     &                       GAO,
     &                       GAO1,
     &                       DMAT,
     &                       BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION VECTOR_DERVS(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2212
      IF(LGET) IOPTION = 1212

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 VECTOR_DERVS,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DIVM */
      SUBROUTINE GET_DIVM(LGET,
     &                    SCALAR,
     &                    IREPOP,
     &                    GAO,
     &                    GAO1,
     &                    DMAT,
     &                    BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2223
      IF(LGET) IOPTION = 1223

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 SCALAR,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_DIVJ */
      SUBROUTINE GET_DIVJ(LGET,
     &                    SCALAR,
     &                    IREPOP,
     &                    GAO,
     &                    GAO1,
     &                    DMAT,
     &                    BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2222
      IF(LGET) IOPTION = 1222

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 SCALAR,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_ROTJPARA */
      SUBROUTINE GET_ROTJPARA(LGET,
     &                        VECTOR,
     &                        IREPOP,
     &                        GAO,
     &                        GAO1,
     &                        DMAT,
     &                        BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#include "dcbbas.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      joff = NBAST
#else
      joff = NTBAS(0)
#endif

      K = 1
      L = 1
      M = 2
      N = 2
      A = D1
      B = D1
      ITRANSFORM = 0

      DO ixyz = 1,3

        jxyz = ixyz + 1
        kxyz = jxyz + 1

        IF(jxyz .GT. 3) jxyz = jxyz - 3
        IF(kxyz .GT. 3) kxyz = kxyz - 3

        CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                       VECTOR(ixyz),
     &                       A,B,K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+joff*(jxyz-1)),
     &                       GAO1(1+joff*(kxyz-1)),
     &                       DMAT,BUF)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(ixyz),
     &                       -A,-B,K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+joff*(kxyz-1)),
     &                       GAO1(1+joff*(jxyz-1)),
     &                       DMAT,BUF)

      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_GR_NH */
      SUBROUTINE GET_GR_NH(LGET,
     &                     VECTOR,
     &                     IREPOP,
     &                     GAO,
     &                     GAO1,
     &                     DMAT,
     &                     BUF)
C=======================================================================
C                                                           radovan bast
C                                               last revision: july 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#include "dcbbas.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      joff = NBAST
#else
      joff = NTBAS(0)
#endif

      K = 1
      L = 1
      M = 2
      N = 2
      A = D1
      B = D1
      ITRANSFORM = 0

      DO ixyz = 1,3

        CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                       VECTOR(ixyz),
     &                       A,B,K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO,
     &                       GAO1(1+joff*(ixyz-1)),
     &                       DMAT,BUF)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(ixyz),
     &                       A,B,K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+joff*(ixyz-1)),
     &                       GAO,
     &                       DMAT,BUF)

      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* GET_M_DERVS_NH */
!     SUBROUTINE GET_M_DERVS_NH(LGET,
!    &                          VECTOR_DERVS,
!    &                          IREPOP,
!    &                          GAO,
!    &                          GAO1,
!    &                          DMAT,
!    &                          BUF)
!=======================================================================
!                                                           radovan bast
!                                               last revision: july 2007
!=======================================================================

!include "implicit.h"
!include "priunit.h"

!     PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

!     LOGICAL   LGET
!     DIMENSION VECTOR_DERVS(*),
!    &          GAO(*),
!    &          GAO1(*),
!    &          DMAT(*),
!    &          BUF(*)

!include "dcbbas.h"

!if defined (VAR_G77) || defined (SYS_SX)
!     EQUIVALENCE (NBAST,NTBAS(0))
!     joff = NBAST
!else
!     joff = NTBAS(0)
!endif

!     K = 1
!     L = 1
!     M = 2
!     N = 2
!     A = D1
!     B = D1
!     ITRANSFORM = 2

!     DO ixyz = 1,3
!       DO jxyz = 1,3

!         CALL GET_DENSITY_SUB(LGET,.TRUE.,
!    &                         VECTOR_DERVS(jxyz+3*(ixyz-1)),
!    &                         A,B,K,L,M,N,
!    &                         ixyz,ITRANSFORM,IREPOP,
!    &                         GAO,
!    &                         GAO1(1+joff*(jxyz-1)),
!    &                         DMAT,BUF)
!         CALL GET_DENSITY_SUB(LGET,.FALSE.,
!    &                         VECTOR_DERVS(jxyz+3*(ixyz-1)),
!    &                         A,B,K,L,M,N,
!    &                         ixyz,ITRANSFORM,IREPOP,
!    &                         GAO1(1+joff*(jxyz-1)),
!    &                         GAO,
!    &                         DMAT,BUF)

!       ENDDO
!     ENDDO

!     RETURN
!     END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_JPARA */
      SUBROUTINE GET_JPARA(LGET,
     &                     VECTOR,
     &                     IREPOP,
     &                     GAO,
     &                     GAO1,
     &                     DMAT,
     &                     BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)
      PARAMETER (DP12 = 1.0D0/2.0D0)

#include "dummy.h"

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

      IOPTION = 2211
      IF(LGET) IOPTION = 1211

      IF(.NOT. LGET) CALL DSCAL(3,-DP12,VECTOR(1),1)

      CALL GET_DENSITY(IOPTION,
     &                 DUMMY,
     &                 VECTOR,
     &                 IREPOP,
     &                 GAO,
     &                 GAO1,
     &                 DMAT,
     &                 BUF)

      IF(LGET) CALL DSCAL(3,-DP12,VECTOR(1),1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* GET_ROTVPARA */
      SUBROUTINE GET_ROTVPARA(LGET,
     &                        VECTOR,
     &                        DENSITY,
     &                        DENSITY_GRADIENT,
     &                        IREPOP,
     &                        GAO,
     &                        GAO1,
     &                        DMAT,
     &                        BUF)
C=======================================================================
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      LOGICAL   LGET
      DIMENSION VECTOR(*),
     &          DENSITY_GRADIENT(*),
     &          GAO(*),
     &          GAO1(*),
     &          DMAT(*),
     &          BUF(*)

#include "dcbbas.h"

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
      A = D1/DENSITY
      B = D1/DENSITY
      ITRANSFORM = 0

      DO IXYZ = 1,3

        JXYZ = IXYZ + 1
        KXYZ = JXYZ + 1

        IF(JXYZ .GT. 3) JXYZ = JXYZ - 3
        IF(KXYZ .GT. 3) KXYZ = KXYZ - 3

        CALL GET_DENSITY_SUB(LGET,.TRUE.,
     &                       VECTOR(IXYZ),
     &                       A,
     &                       B,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+JOFF*(JXYZ-1)),
     &                       GAO1(1+JOFF*(KXYZ-1)),
     &                       DMAT,BUF)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(IXYZ),
     &                       -A,
     &                       -B,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+JOFF*(KXYZ-1)),
     &                       GAO1(1+JOFF*(JXYZ-1)),
     &                       DMAT,BUF)

        C = -DENSITY_GRADIENT(JXYZ)/(DENSITY*D2)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(IXYZ),
     &                       A*C,
     &                       B*C,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO,
     &                       GAO1(1+JOFF*(KXYZ-1)),
     &                       DMAT,BUF)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(IXYZ),
     &                       -A*C,
     &                       -B*C,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+JOFF*(KXYZ-1)),
     &                       GAO,
     &                       DMAT,BUF)

        C = DENSITY_GRADIENT(KXYZ)/(DENSITY*D2)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(IXYZ),
     &                       A*C,
     &                       B*C,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO,
     &                       GAO1(1+JOFF*(JXYZ-1)),
     &                       DMAT,BUF)
        CALL GET_DENSITY_SUB(LGET,.FALSE.,
     &                       VECTOR(IXYZ),
     &                       -A*C,
     &                       -B*C,
     &                       K,L,M,N,
     &                       1,ITRANSFORM,IREPOP,
     &                       GAO1(1+JOFF*(JXYZ-1)),
     &                       GAO,
     &                       DMAT,BUF)

      ENDDO

      RETURN
      END
