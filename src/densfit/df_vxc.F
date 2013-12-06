!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
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

C /* Deck DF_AuxVXC_UpdateCFit */
      SUBROUTINE DF_AuxVXC_UpdateCFit(DMAT,CFIT,WORK,LWORK)
C
C Description
C -----------
C This routine is an auxiliary routine for the evaluation of the 
C exchange-correlation potential using density C fitting, 
C following the scheme of Birkenheuer et at (IJQC 102, 743, 2005)
C
C In it we call numerical integration integration routines to set
C the quantity <V_{XC}[\rho_{fit}]|\alpha> (alpha=fit function), 
C and then with it and C with Gjk (the inverse of Vab = [(a|b)]), 
C which was previously computed. With these quantities we calculate 
C
C      b_k = \sum_j <V_{XC}[\rho_{fit}]|\alpha_j> G_{jk} 
C
C in order to update the fit coefficients and therefore include the
C exchange-correlation contribution.
C
C Variables
C ---------
C CFIT : fit coefficients.
C DMAT : density matrix for the orbital basis
C
C Results
C -------
C Upon completion, the fit coefficients will be updated as
C
C    c_k (with XC) = c_k (original) + b_k 
C
C History
C -------
C 2006-09-21: pilot implementation by Andre Gomes and Radovan Bast, Strasbourg
C
C Remarks
C -------
C at this point, both Gjk and the integral list are C read from 
C files "FITXCG" and "FITVXC", respectively, in order to avoid 
C much change to the current function calls.
C 
        use dft_cfg

#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#else
#include "gnrinf.h"
#endif
#include "densfit.h"
      DOUBLE PRECISION CFIT(*), DMAT(*), WORK(*)
      integer :: response_order
#include "memint.h"
C
C first step is to numerically integrate <V_{XC}[\rho_{fit}]|\alpha>
C

#ifdef PRG_DIRAC

      if (using_diffdens) then
       call quit('unsupported option')
      endif

      response_order = 0
      if (dodflr) then
        response_order = 1
      end if

!fixme call XC integrator

#endif /* ifdef PRG_DIRAC */


C
C after that, we get gkj and the integrals above from file...
C
      CALL MEMGET('REAL',KGJK,NBASISAUX*NBASISAUX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVXCI,NBASISAUX,WORK,KFREE,LFREE)

      CALL READ_MATRIX(WORK(KGJK),NBASISAUX,NBASISAUX,'FITXCG')

      if (DODFLR) then 
         CALL READ_MATRIX(WORK(KVXCI),1,NBASISAUX,'LRDFVX')
      else
         CALL READ_MATRIX(WORK(KVXCI),1,NBASISAUX,'FITVXC')
      endif

C .. and create b_k = \sum_j [<integral above>]Gjk 
C    which is then used to update cfit...
C
      DO K = 1, NBASISAUX
         BK = DDOT(NBASISAUX,WORK(KVXCI),1,
     &             WORK(KGJK + NBASISAUX*(K-1)),1)
         CFIT(K)  = CFIT(K) + BK 
      ENDDO
C
C and we're done!
C
      CALL MEMREL('DF_Set_AuxVXC_CFit-2',WORK,KVXCI,KVXCI,KFREE,LFREE)
      CALL MEMREL('DF_Set_AuxVXC_CFit-1',WORK,KGJK,KGJK,KFREE,LFREE)
      END

C /* Deck DF_FitXC_Build_Gjk */
      SUBROUTINE DF_FitXC_Build_GKJ(G,IDimK,IdimJ,Vab,ConsQInt,CFitQ)
C
C Description 
C -----------
C This subroutine calculates the matrix Gkj, used in the procedure 
C outlined by Birkenheuer et at (IJQC 102, 743, 2005) to use quantities 
C from density fitting in the construction of the exchange-correlation 
C potential. 
C 
C Variables
C ---------
C Gjk         : "projected inverse" of the charge fit overlap matrix Vab
C Vab         : copy of the (alpha|beta) matrix
C IDimK,IDimJ : dimensions of Gkj
C ConsQInt    : auxiliary integrals to define the charge conservation 
C               constraint
C CFitQ       : coefficients for the charge conservation constraint
C
C                 \sum_alpha (alpha|beta) CFitQ_alpha = ConsQInt_alpha
C
C Results
C -------
C Upon return, the matrix Gkj will hold the projected inverse,
C
C G_{kj} = \sum_m 
C           [ 
C            \delta_{km} 
C             - \frac{CFitQ_k  ConsQInt_m}{\sum_l CFitQ_l ConsQInt_l}
C              ]
C               (Vab)^{-1}_{mj}
C
C and Vab will be destroyed and replaced by its inverse
C
C History
C -------
C 2006-09-21: pilot implementation by Andre Gomes and Radovan Bast, Strasbourg 
C
C Remarks
C ------
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "densfit.h"
#include "priunit.h"
Caspg these array sizes have to be changed, dynamic allocation needed...
         LOGICAL DEBUG
         INTEGER IPIV(MXPRIM), LWORK
         DOUBLE PRECISION Vab(*),WORK(MXPRIM),ConsQInt(*),CFitQ(*),G(*)
C#include "memint.h"
C

C first step: invert Vab. 
C
C ... We begin by performing a LU factorization of it...
      LDA   = MAX(1,IDimK)
      LWORK = 2*LDA 

      CALL DGETRF(IDimK,IDimJ,Vab,LDA,IPIV,INFO)
      IF ( INFO.NE.0 ) THEN
         WRITE (LUPRI,*) 'DGETRF returns ',INFO
         CALL QUIT('DF_FitXC_Build_Gjk: Error in factoring Vab')
      ENDIF

C ... which is then used to compute the inverse, and place it in Vab
      CALL DGETRI(LDA,Vab,LDA,IPIV,WORK,LWORK,INFO)
      IF ( INFO.NE.0 ) THEN
         WRITE (LUPRI,*) 'DGETRI returns ',INFO
         CALL QUIT('DF_FitXC_Build_Gjk: Error in inverting Vab')
      ENDIF

C
C second step: get G
C
C ... first the denominator in bracket in the formula above is calculated,
C     provided we have the charge conservation stuff ... 
      IF (DF_CONSTR) THEN
         DENOM = DDOT(NBASISAUX,ConsQInt,1,CFitQ,1)
C        SAFE  = DLAMCH('SFMIN')
         SAFE  = 10E-16 
         IF (DENOM.LT.SAFE) THEN
            WRITE (LUPRI,*) 'Denominator is too small!'
            CALL QUIT('DF_FitXC_Build_Gjk: Error in setting up Gkj')
         ENDIF
      ENDIF ! df_constr

C ... then we loop over the indexes of G, forming the value in brackets
C     and multiplying it with the  ...
      DO K = 1, IDimK
         DO J = 1, IDimJ
            IF (DF_CONSTR) THEN
               dsum = 0
               DO M = 1, IDimJ
                  VINV = Vab((M-1)*IDimK + J)
                  A    = CFitQ(K) * ConsQInt(M)  
                  term = (IDLTA(K,M)-(A/DENOM))*VINV
                  dsum = dsum + term
               ENDDO
               G((K-1)*IDimK+J) = dsum 
            ELSE
               VINV = Vab((K-1)*IDimK + J)
               G((K-1)*IDimK+J) = VINV 
            ENDIF
         ENDDO
      ENDDO
C ... and now we are done! 
      END
C /* Deck IDLTA */
      FUNCTION IDLTA(I,J)
C
C aspg, 2006-09-19
C calculates \delta_{ij}
C
      INTEGER I,J
      IF (I.EQ.J) THEN
         IDLTA = 1
      ELSE
         IDLTA = 0 
      ENDIF
      RETURN
      END
C /* Deck SAVE_MATRIX */
      SUBROUTINE SAVE_MATRIX(MATRIX,IDIMA,IDIMB,LABEL)
C
C aspg, 2006-09-19
C functions to save matrices/vectors in the density fitting procedure
C
#include "implicit.h"
      CHARACTER*6 LABEL
      CHARACTER*20 FMT
      LOGICAL DEBUG
      DOUBLE PRECISION MATRIX(*)

C     DEBUG = .TRUE.
      DEBUG = .FALSE.
C
      IF (DEBUG) THEN
         write (6,*) 'debug, in save_matrix:',idima,idimb,label
      ENDIF
    
      IUNIT = -1
      CALL GPOPEN(IUNIT,LABEL,'UNKNOWN','SEQUENTIAL','FORMATTED',
     &     IDUM,LDUM)
      WRITE(IUNIT,'(2I8)') IDIMA,IDIMB  
      DO I=1,IDIMA
        WRITE(IUNIT,'(66(D37.28))')
     &          (MATRIX((I-1)*IDIMA+J), J=1,IDIMB)
        IF (DEBUG) THEN
           FMT='(2(I8)1(D37.28))'
           write(6,*) FMT
           WRITE(6,FMT)
     &          (I,J,MATRIX((I-1)*IDIMA+J), J=1,IDIMB)
        ENDIF
      ENDDO
      CALL GPCLOSE(IUNIT,'KEEP')
      RETURN
      END

C /* Deck READ_MATRIX */
      SUBROUTINE READ_MATRIX(MATRIX,IDIMA,IDIMB,LABEL)
C
C aspg, 2006-09-19
C functions to save matrices/vectors in the density fitting procedure
C
#include "implicit.h"
      CHARACTER*6 LABEL
      CHARACTER*20 FMT
      LOGICAL DEBUG
      DOUBLE PRECISION MATRIX(*)

C     DEBUG = .TRUE.
      DEBUG = .FALSE.
      IUNIT = -1
      CALL GPOPEN(IUNIT,LABEL,'UNKNOWN','SEQUENTIAL','FORMATTED',
     &            IDUM,LDUM)
      READ(IUNIT,'(2I8)') IDIMC,IDIMD 
      IF ((IDIMA.NE.IDIMC).OR.(IDIMB.NE.IDIMD)) THEN
        CALL QUIT('Dimensions do not match in READ_MATRIX')
      ELSE
          DO I=1,IDIMA
            READ(IUNIT,'(66(D37.28))')
     &              (MATRIX((I-1)*IDIMA+J), J=1,IDIMB)
            IF (DEBUG) THEN
               FMT='(66(D37.28))'
               write(6,*) FMT
               WRITE(6,FMT)
     &              (MATRIX((I-1)*IDIMA+J), J=1,IDIMB)
            ENDIF
          ENDDO
          CALL GPCLOSE(IUNIT,'KEEP')
      ENDIF
      RETURN
      END
C /* Deck INIT_DD_FITVXC_FILES */
C
      SUBROUTINE DF_FitXC_InitFilesDD(NBASIS,WORK,LWORK)
C
C Description
C -----------
C This subroutine initializes the files which will hold the information
C regarding fit coefficients and auxiliary coefficients for the evaluation
C of the exchange-correlation contribution, when differential densities
C are used in the SCF procedure... This is so because we cannot feed
C "differential" fit coefficients to obtail the density in the DFT part,
C and also cannot use the B_k coefficients in full, since for the SCF
C the C_k's are differential...
C
C Arguments
C ---------
C NBASIS: number of auxiliary basis functions
C WORK  : work array
C LWORK : current positon of memory 
C
C Results
C -------
C After it is called, the files "FITDDC" and "FITDDB" will be created
C and its elements set to zero.
C
C "FITDDC" will contain the "non-differential" fit coefficients,
C that is, the sum of the differential fit coefficients accumulated during 
C the scf proc. 
C
C "FITDDB" will containing the coefficients B_k which multiply the 
C exchance-correlation contribution from a previous scf iteration.
C
C Remarks
C -------
C Should be placed before routines that perform the scf iterations.
C If it is somewhere else then there's something wrong...
C
C History
C -------
C 2006-10-15: Initial implementation, Andre Gomes
C
#include "implicit.h"
      DOUBLE PRECISION WORK(*)
#include "memint.h"
C
      CALL MEMGET('REAL',KPREV,NBASIS,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KPREV),NBASIS)
C
      CALL SAVE_MATRIX(WORK(KPREV),1,NBASIS,'FITDDB')
      CALL SAVE_MATRIX(WORK(KPREV),1,NBASIS,'FITDDC')
C
      CALL MEMREL('Init_DD_FitVXC_Files',WORK,KPREV,KPREV,KFREE,LFREE)
      RETURN
      END
