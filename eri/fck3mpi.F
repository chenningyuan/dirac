C 
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
C aspg, 2006-06-19
C
C since the eri/fck3/densfit code in not really available in parallel form
C at this time, we apply some tricks to get the program to compile in 
C parallel mode. These are to have dummy funcion in place of MPI3INI
C and have FCK3COLLECT, which sits in the old fck3driver and shouldn't
C need to be compiled, available here. 

#if defined (VAR_MPI)
#define USE_FCK3COLLECT
#define USE_PARALLEL_ERI_DUMMY
#endif

#ifdef USE_PARALLEL_ERI_DUMMY
C aspg, 2006-06-19
C for the time being put a dummy call to the parallel routines 
C
      SUBROUTINE MPI3INI(DMAT, NDMT, GABAB, GABAL, DGABAL, WORK, LWORK)
         CALL QUIT('Parallel Implementation of FCK3 not available yet.')
      RETURN
      END
#else 
      SUBROUTINE MPI3INI(DMAT, NDMT, GABAB, GABAL, DGABAL, WORK, LWORK)
#include "implicit.h"
#if defined (VAR_MPI)
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
#include "cbifck3.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3thr.h"
#include "densfit.h"
#include "pbc.h"
#include "infpar.h"
#include "maxmom.h"
#include "xyzpow.h"
#include "maxaqn.h"
#include "gamcom.h"
#include "nuclei.h"
#include "inforb.h"
      DIMENSION DMAT(*), WORK(*), GABAB(*), GABAL(*), DGABAL(*)
      IPRINT = 0
      IPRTYP = 8
      IF (MYNUM.EQ.0) THEN
        CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
        CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
        CALL DF_DIMS()
        CALL FCK3_DIMS(NDMT)
        IF (DOFMM) CALL QUIT('Parallel code not yet implemented for mm')
        IF (DOPBC)
     &    CALL QUIT('Parallel code not yet implemented for pbc')
c     ELSE
c       CALL FCK3INI
      ENDIF

C     ERI initialization
      CALL MPIXBCAST(ISTEP,MXAQNM,'INTEGER',MASTER)
      CALL MPIXBCAST(MVAL,MXAQNM,'INTEGER',MASTER)
      CALL MPIXBCAST(NVAL,MXAQNM,'INTEGER',MASTER)
c     CALL MPIXBCAST(MAXJ,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TABFJW,121*(MAXJ+7),'DOUBLE',MASTER)
      CALL MPIXBCAST(WVAL,MAXJ+1,'DOUBLE',MASTER)


#ifdef MOD_PBC
C     PBC initializations
      CALL MPIXBCAST(DOPBC,1,'LOGICA',MASTER)
      CALL MPIXBCAST(PBCACTIVE,1,'LOGICA',MASTER)
#endif


C     FCK3 initializations
      CALL MPIXBCAST(DOFCK3,1,'LOGICA',MASTER)
      CALL MPIXBCAST(FIXTHR,1,'LOGICA',MASTER)
      CALL MPIXBCAST(CHCKNAI,1,'LOGICA',MASTER)
      CALL MPIXBCAST(SPEXP,1,'LOGICA',MASTER)
      CALL MPIXBCAST(NOSCRPRIM,1,'LOGICA',MASTER)
      CALL MPIXBCAST(NOSKIPCL,1,'LOGICA',MASTER)
      CALL MPIXBCAST(SCRFBTCH,1,'LOGICA',MASTER)
      CALL MPIXBCAST(NCLFCK3,1,'LOGICA',MASTER)
      CALL MPIXBCAST(PACKFTUV,1,'LOGICA',MASTER)
      CALL MPIXBCAST(SKALLCL,1,'LOGICA',MASTER)
      CALL MPIXBCAST(OLDMOM,1,'LOGICA',MASTER)
      CALL MPIXBCAST(FULL_J,1,'LOGICA',MASTER)
      CALL MPIXBCAST(CONTRACTED,1,'LOGICA',MASTER)
      CALL MPIXBCAST(BRFREE,1,'LOGICA',MASTER)
      CALL MPIXBCAST(NCSODB,1,'LOGICA',MASTER)
      CALL MPIXBCAST(NCSTUV,1,'LOGICA',MASTER)
c     CALL MPIXBCAST(,1,'LOGICA',MASTER)
c     CALL MPIXBCAST(,1,'LOGICA',MASTER)
      CALL MPIXBCAST(IPRFCK3,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IFTHRS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(THRCLS,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(THRSPEXP,1,'DOUBLE',MASTER)
      
      
C     General initializations
      CALL MPIXBCAST(DMAT,N2BASX*NDMT,'DOUBLE',MASTER)
      IF (DOFIT) THEN
        CALL MPIXBCAST(GABAB,N2BASX*NDMT,'DOUBLE',MASTER)
        CALL MPIXBCAST(GABAL,NBASISAUX*NDMT,'DOUBLE',MASTER)
        IF (RHSFIT)
     &    CALL MPIXBCAST(DGABAL,NBASISAUX*NDMT,'DOUBLE',MASTER)
      ENDIF

C     Readin initializations
C     /NUCLEI/
      CALL MPIXSCOM(CHARGE,NPAUX,4)
C     /NUCLEC/
      CALL MPIXSCOM(NAMEX,NAMN(MXCENT),4)
C     /MULBAS/
      CALL MPIXBCAST(MULBSI,MXCENT,'INTEGER',MASTER)
#endif      

      RETURN
      END

#if defined (VAR_MPI)

      SUBROUTINE FCK3DRV_SLAVE(WORK,MWORK,IPRINT)
#include "implicit.h"
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "densfit.h"
      DIMENSION WORK(MWORK)
      CALL DF_DIMS()
      IF (.NOT. DOFIT) NBASISAUX=0
      CALL FCK3_DIMS(NDMT)
      KFMAT   = 2
      KDMAT   = KFMAT   + N2BASX*NDMT
      KGABAB  = KDMAT   + N2BASX*NDMT
      KGABAL  = KGABAB  + N2BASX*NDMT
      KDGABAL = KGABAL  + NBASISAUX*NDMT
      KLAST   = KDGABAL + NBASISAUX*NDMT
      IF (KLAST .GT. MWORK) CALL STOPIT('DF_SLAVE','NODSTR',KLAST,MWORK)
      LFREE = MWORK - (KLAST - 2)
      KFREE = KLAST
      CALL DZERO(WORK(KFMAT),N2BASX*NDMT)
      CALL FCK3DRV (WORK(KFMAT),WORK(KDMAT),NDMT,WORK(KGABAB),
     &              WORK(KGABAL),WORK(KDGABAL),WORK(KFREE),LFREE)

      RETURN
      END



      SUBROUTINE DF_DIMS()
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mxcent.h"
#include "nuclei.h"
#include "densfit.h"
      INCLUDE "mpif.h"
      CALL MPI_BCAST(NBASISAUX,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOFIT,1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,IERR)
      RETURN
      END                  
      
#endif
#endif /* USE_PARALLEL_ERI_DUMMY */



#if defined (USE_FCK3COLLECT)
      SUBROUTINE fck3collect(FMAT,NDMAT,tm,WORK,LFREE,IPRINT)
#include "implicit.h"
#include "inforb.h"
      DIMENSION FMAT(N2BASX,NDMAT), WORK(*)
      INCLUDE 'mpif.h'
#include "maxorb.h"
#include "infpar.h"
#include "priunit.h"
c
      DIMENSION tms(NODTOT+1)
      NSIZE = N2BASX*NDMAT
      IF (LFREE .LT. NSIZE) THEN
         WRITE(LUPRI,"(I9,'*',I2,'=',I9,' needed, ',I9,' available.')")
     &        N2BASX,NDMAT,N2BASX*NDMAT,LFREE
        Call QUIT('Insufficient memory in fck3collect')
      END IF
      CALL DCOPY(NSIZE,FMAT,1,WORK,1)
      CALL MPI_REDUCE(WORK,FMAT,NSIZE,MPI_DOUBLE_PRECISION,MPI_SUM,
     &     MASTER,MPI_COMM_WORLD,IERR)
      iF(iprint.ge.0) THEN
c     for a good measure, collect times to to provide more info about
c     load balancing.
         CALL MPI_GATHER(tm,1,MPI_DOUBLE_PRECISION,
     &        tms,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,IERR)
         if(IERR.NE.0) WRITE(LUPRI,*) 'error collecting times:',IERR
         WRITE(LUPRI,"('Node times:',/,20(10F7.2,/))")
     &        (tms(i),i=1,NODTOT+1)
      END IF
      END
#endif /* USE_FCK3COLLECT */
