#include <iratdef.h>
#include <lrecao.h>
#if defined (VAR_OLDMOLECUL)
#if !defined (SYS_CRAY) && !defined (SYS_T3D)
      PARAMETER (LRINT = LRECAO, LIT = 500)
C     ... LRINT = MAX(2*901,LRECAO)
C     ... 901 is sufficient for LUINTA buffer length 600 when irat=2
#else
      PARAMETER (LRINT = LRECAO, LIT = 500)
C     ... LRINT = MAX(2*1201,LRECAO)
C     ... 1201 is sufficient for LUINTA buffer length 600 when irat=1
#endif
C         lbuf = 600, lrint = (IRAT*lbuf+lbuf+1 - 1)/irat + 1
C         LRECAO is needed for LUSORT in ORDER
      COMMON/N_SORT/ RINT(LRINT),IT(LIT),IBATCH(666),LASTAD(106),
     *             KEEP(8),NBAS(8),IS(LIT),IBAS(8),NREC(106),NSYM,
     *             IDATA(18),NSOINT(106),LINT,LBATCH,LUSORT,LUORD,
     *             NSOBAT(106),INTSYM,IPRINT,IPRFIO,IW
      INTEGER IINT(IRAT*LRINT)
      EQUIVALENCE (RINT, IINT)
#else
      PARAMETER (LRINT = 2*LRECAO, LIT = 500)
      COMMON/N_SORT/ RINT(LRINT,IT(LIT),IBATCH(666),LASTAD(106),
     *             KEEP(8),NBAS(8),IS(LIT),IBAS(8),NREC(106),NSYM,
     *             IDATA(18),NSOINT(106),LINT,LBATCH,LUSORT,LUORD,
     *             NSOBAT(106),INTSYM,IPRINT
#endif
