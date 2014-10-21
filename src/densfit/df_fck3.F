!dirac_copyright_start
!      Copyright (c) 2014 by the authors of DIRAC.
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

C  /* Deck lin_df_coulomb */
      SUBROUTINE LIN_DF_Coulomb (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                           WORK,LWORK)
C*****************************************************************************
C    
C     Driver routine for linear-scaling density-fitting Coulomb
C        
C     Written by Simen Reine, february 2006
C     
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "iteration.h"
      DIMENSION DMAT(NBASIS*NBASIS,NDMT), FMAT(NBASIS*NBASIS,NDMT), 
     &          WORK(LWORK)
      LOGICAL IS_DIAGZERO
#include "memint.h" 
C
      CALL QENTER('LIN_DF')
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C     Simplified test to chech if FMAT is zero
C
      IF (.NOT.IS_DIAGZERO(FMAT,NBASIS))
     & CALL QUIT('FMAT must be dzeroed in LIN_DF_Coulomb')
     
      IPRINT = IPRFCK

      SYMDENS = (IFCTYP/10).EQ.1

      IF (.NOT. LIN_DF) THEN
cSimen
        DO IDMT = 1, NDMT
cSimen
c       CALL DENS_FIT(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,WORK,LWORK)
        CALL DENS_FIT(FMAT(1,IDMT),DMAT(1,IDMT),1,ISYMDM,IFCTYP,
     &                IPRFCK,WORK,LWORK)
cSimen
        ENDDO
cSimen
      ELSE
        IF (.NOT.DOFCK3) THEN
          WRITE(LUPRI,*) 'LIN_DF_Coulomb driver implemented ',
     &                   'with FCK3-scheme only (Activate with',
     &                   'keyword .DOFCK3 under **INTEG)'
          CALL QUIT('Ended in LIN_DF_Coulomb')
        ENDIF

cTODO   Create projected (alpha|beta) matrix


C       CFIT   Full fitting coefficients
C       GABAB  GAB (g_abab) matrix
C       GABAL  Density fitting GALPHA
C              (alpha|alpha) vector
C       AUXINT Charge integrals
           CALL MEMGET ('REAL',KCFITSAV,NBASISAUX*NDMT,WORK,KFREE,LFREE)
        CALL MEMGET ('REAL',KCFIT,NBASISAUX*NDMT,WORK,KFREE,LFREE)
        CALL MEMGET ('REAL',KGABAB,NBASIS*NBASIS,WORK,KFREE,LFREE)
        CALL MEMGET ('REAL',KGABAL,NBASISAUX,WORK,KFREE,LFREE)
        CALL MEMGET ('REAL',KAUXINT,NBASISAUX,WORK,KFREE,LFREE)

C       We will have no symmetry between density 1 and 2
        ODTR12 = .FALSE.

C       Determines the dimensions of the density matrix.
        CALL DMATDIM

C       Generate the 2-index Cauchy-Schwarz integrals, and normalize the 
C       auxilliary basis with the 1-index auilliary Cauchy-Schartz integrals.
        CALL SCRINT(WORK(KGABAB),WORK(KGABAL),WORK(KFREE),LFREE,IPRINT)

C       Calculate attenuated integrals, i.e. over the oparator:
C       erfc(omega*r12)/r12 instead of 1/r12
        CALL SET_ATTGAM(DF_ATTENU,DF_ATTPAR)

        CALL SETFCK3BOX(KNUMBOX,KCORDBOX,LFREE,KFREE,WORK,LFREE,IPRINT)

        NUMBOXES = NUMBER_OF_BOXES()
        WRITE(LUPRI,*) 'NUMBOXES = ', NUMBOXES 

C       Call to grenerate the one-electron integrals over the fit set
        CALL AUXCHA(WORK(KAUXINT),WORK(KFREE),LFREE,IPRINT)
        
cSimen
        CALL DZERO(WORK(KCFIT),NBASISAUX*NDMT)
cSimen
        DO IBOX=1,NUMBOXES
          CALL DF_BOX(WORK(KCFIT),IBOX,WORK(KNUMBOX),
     &                WORK(KCORDBOX),FMAT,DMAT,NDMT,WORK(KGABAB),
     &                WORK(KGABAL),ISYMDM,IFCTYP,WORK(KFREE),LFREE,
     &                IPRINT)
        ENDDO ! Loop over IBOX

        IF (IPRINT.GT.25) THEN
          CALL HEADER('CFIT',-1)
          CALL OUTPUT(WORK(KCFIT),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
        ENDIF

        SU = 0.0d0
        DO I=1,NBASISAUX
          SU = SU + WORK(KCFIT+I-1) * WORK(KAUXINT+I-1)
        ENDDO
        Write(lupri,*) ' '
        Write(lupri,*) 'Total density-fitted charge ', SU

        IPRINT = 1

        CALL SET_ATTGAM(.FALSE.,0.0d0)

        CALL TIMER('ITSOLV',TIMSTR,TIMEND)

Caspg+rbast
        IF (DF_FITVXC) THEN
           CALL TIMER('START ',TIMSTR,TIMEND)
           CALL DCOPY  (NBASISAUX*NDMT,WORK(KCFIT),1,WORK(KCFITSAV),1)

           CALL DF_AuxVXC_UpdateCFit(DMAT,WORK(KCFIT),WORK(KFREE),LFREE)
           CALL TIMER('VxcUpFit',TIMSTR,TIMEND)
        ENDIF
Caspg+rbast

        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL SETAOBTCH(SPEXP,CONTRACTED,WORK(KFREE),LFREE,IPRINT)
        CALL ACT_ALL_BATCHES
        CALL FITFCK3(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,
     &       WORK(KGABAB),WORK(KGABAL),WORK(KCFIT),IPRFCK,
     &       WORK(KFREE),LFREE)
        CALL TIMER('JLINBU',TIMSTR,TIMEND)

Caspg+rbast
        IF (DF_FITVXC) THEN
           CALL TIMER('START ',TIMSTR,TIMEND)
           CALL DF_CAL_ApproxERRKS(DMAT,FMAT,WORK(KCFITSAV),NDMT,ISYMDM,
     &                       IFCTYP,WORK(KGABAL),WORK(KFREE),LFREE,
     &                       IPRINT)
           CALL TIMER('AERRKS',TIMSTR,TIMEND)
        ENDIF
Caspg+rbast

        IF (ROBUST_E) THEN
          IF (DF_FITVXC) THEN
             WRITE (LUPRI,*) 'WARNING! Robust energy correction not '//
     &                       'yet implemented for fitted XC runs.' 
          ELSE
             CALL DF_CAL_ROB_E(DMAT,FMAT,WORK(KCFIT),NDMT,ISYMDM,IFCTYP,
     &                         WORK(KGABAL),WORK(KFREE),LFREE,IPRINT)
          ENDIF
        ENDIF

        IF (IPRFCK.GT.25) THEN
          CALL HEADER('Fitted Fock matrix',-1)
          CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        ENDIF
C
C       Release the memory allocated for the Fock matrix build. This
C       will also check for errors due to memory bound failures
C
        CALL MEMREL('LIN_DF',WORK,KAUXINT,KAUXINT,KFREE,LFREE)
        CALL MEMREL('LIN_DF',WORK,KGABAL,KGABAL,KFREE,LFREE)
        CALL MEMREL('LIN_DF',WORK,KGABAB,KGABAB,KFREE,LFREE)
        CALL MEMREL('LIN_DF',WORK,KCFIT,KCFIT,KFREE,LFREE)
        CALL MEMREL('LIN_DF',WORK,KCFITSAV,KCFITSAV,KFREE,LFREE)
      ENDIF
      CALL QEXIT('LIN_DF')
      RETURN
      END
C  /* Deck dens_fit */
      SUBROUTINE DENS_FIT (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                     WORK,LWORK)
C*****************************************************************************
C    
C     DENS_FIT: Entree routine for density fitting in Dirac & Dalton.
C        
C     Written by Luuk Visscher, october 2004
C     Mdified by Simen Reine, april 2005
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "iteration.h"
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
      LOGICAL ROBUST_3I
#include "memint.h" 
C
      CALL QENTER('DENS_FIT')
      CALL TIMER('START ',TIMSTR,TIMEND)

      WRITE(LUPRI,*) 'Entering obsolete DENS_FIT driver'

      IPRINT = IPRFCK
C     Roubust driver is needed if the fitting coefficients have a 
C     significant error. Currently turned off, but might prove useful 
C     for larger systems. /SR
      ROBUST_3I = .FALSE.
C
C     DENS_F contains the fitted density of an earlier run. Use this 
C     fitted density if it exists for the first iteration. 
C     Might be useful for initial guess; projecting a smaller basis 
C     to a larger one. Currently turned off. /SR
      INQUIRE(file='DENS_F',EXIST=F_DENS)
c      F_DENS = F_DENS .AND. (ITER_SCF.EQ.0)
      F_DENS = .FALSE.

      SYMDENS = (IFCTYP/10).EQ.1

      IF (.NOT.DOFCK3) THEN
        WRITE(LUPRI,*) 'New DENS_FIT driver implemented ',
     &       'with FCK3-scheme. Must either use ATOMDF'
        WRITE(LUPRI,*) 'and DOFCK3, or ATOMDF and OLD_DF'
        CALL QUIT('Ended in DENS_FIT')
      ENDIF

C
C     Allocate memory for fit coefficients, work arrays,
C     contraction coefficients and index array
C
      CALL MEMGET ('REAL',KCFIT,NBASISAUX*NDMT,WORK,KFREE,LFREE)
C     Density fitting (alpha|beta) matrix
      CALL MEMGET ('REAL',KALBE,NBASISAUX*NBASISAUX,WORK,KFREE,LFREE)
C     GAB (g_abab) matrix
      CALL MEMGET ('REAL',KGABAB,NBASIS*NBASIS,WORK,KFREE,LFREE)
C     Density fitting GALPHA (alpha|alpha) vector
      CALL MEMGET ('REAL',KGABAL,NBASISAUX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KWVEC,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KGABAL_IT,NBASISAUX*NDMT,WORK,KFREE,LFREE)

C     For ROBUST_3I fitting scheme only
      IF (ROBUST_3I) THEN
        CALL MEMGET ('REAL',KALAB,NBASISAUX*NBASIS*NBASIS*NDMT,
     &       WORK,KFREE,LFREE)
      ENDIF

      WRITE(LUPRI,*) 'Memory requirement in DENS_FIT ', KFREE - 1

C     We will have no symmetry between density 1 and 2
      ODTR12 = .FALSE.

C     Determines the dimensions of the density matrix.
      CALL DMATDIM

C     Set AO- and OD-batches (and classes)
      RHSFIT = .FALSE.
      LHSFIT = .FALSE.
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,
     &     IPRFCK3)
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &     KFREE,LFREE,WORK,LFREE,DUMMY,DUMMY,DUMMY,DUMMY,
     &     NDMT,IPRFCK3,.false.)

C     Activates all batches for generation of (ab|ab), 
C     (alpha|alpha) and (alpha|beta) integrals.
      CALL ACT_ALL_BATCHES

C     Generate the 2-index Cauchy-Schwarz integrals
      CALL GETINT3(WORK(KGABAB),WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.FALSE.,.FALSE.,
     &     IDUMMY,IDUMMY,.TRUE.)
C      CALL TIMER('GABAB ',TIMSTR,TIMEND)

C     Generate the 1-index Cauchy-Schwarz integrals over the fit set 
      CALL GETINT3(WORK(KGABAL),WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     IDUMMY,IDUMMY,.TRUE.)



C     Normalize the fit set.
      CALL NRMFIT(WORK(KGABAL),WORK(KFREE),LFREE,IPRINT)
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,
     &     IPRFCK3)
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &     KFREE,LFREE,WORK,LFREE,DUMMY,DUMMY,DUMMY,DUMMY,
     &     NDMT,IPRFCK3,.false.)

C     Generate the normalized 1-index Cauchy-Schwarz integrals over the 
C     fit set (should produce ones)
C     FIXME: these are now obsolete for screening and should be removed.
      CALL GETINT3(WORK(KGABAL),WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     IDUMMY,IDUMMY,.TRUE.)
C      CALL TIMER('GABALN',TIMSTR,TIMEND)

      IF (.FALSE.) THEN
C       Tests the normalization of the density-fitting basis
        DEVIATION = 0.0d0
        DO I=1,NBASISAUX
          DEVIATION = MAX(DEVIATION,ABS(WORK(KGABAL+I-1)-1.0d0))
        ENDDO
        WRITE(LUPRI,*) 'Maximum df-normilization deviation ', DEVIATION
        THRESHOLD = 1.0d-14
        IF (DEVIATION .GT. THRESHOLD) 
     &       CALL QUIT('Deviation greater than threshold')

c       Writes the (alpha|beta) matrix to file.
        IUNIT = -1
        CALL GPOPEN(IUNIT,'ALBEMA','UNKNOWN','SEQUENTIAL','FORMATTED',
     &       IDUM,LDUM)
        DO I=1,NBASISAUX
          WRITE(IUNIT,'(66(D27.18))') 
     &         (WORK(KALBE+(I-1)*NBASISAUX+J-1), J=1,NBASISAUX)
        ENDDO
        CALL GPCLOSE(IUNIT,'KEEP')
      ENDIF

C     Calculate attenuated integrals, i.e. over the oparator:
C     erfc(omega*r12)/r12 instead of 1/r12
      CALL SET_ATTGAM(DF_ATTENU,DF_ATTPAR)

C     Generate the 2-index repulsion integrals over the fit set
      CALL GETINT3(WORK(KALBE),WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     IDUMMY,IDUMMY,.FALSE.)
c     CALL TIMER('ALBETA',TIMSTR,TIMEND)


      IF (IPRINT.GT.25) THEN
        CALL HEADER('(ALBE)',-1)
        CALL OUTPUT(WORK(KALBE),1,NBASISAUX,1,NBASISAUX,
     &              NBASISAUX,NBASISAUX,1,LUPRI)
        CALL HEADER('GABAB',-1)
        CALL OUTPUT(WORK(KGABAB),1,NBASIS,1,NBASIS,
     &              NBASIS,NBASIS,1,LUPRI)
        CALL HEADER('GABALP',-1)
        CALL OUTPUT(WORK(KGABAL),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF

C     Roubust driver is needed if the fitting coefficients have a 
C     significant error. Might prove useful. \SR
      IF (ROBUST_3I) THEN
C       Generate the 3-index repulsion integrals
        CALL GETINT3(WORK(KALAB),WORK(KIODBCLS),WORK(KIODBTCH),
     &       WORK(KRODBTCH),
     &       WORK(KFREE),LFREE,IPRINT,.FALSE.,.TRUE.,
     &       IDUMMY,IDUMMY,.FALSE.)
c       CALL TIMER('ABAL  ',TIMSTR,TIMEND)

        CALL ROBUST_3I_DRV(FMAT,DMAT,WORK(KALBE),NBASISAUX,WORK(KALAB),
     &       NBASIS,NBASIS,WORK(KFREE),LFREE,IPRINT)
c       CALL TIMER('ROBUST_3I',TIMSTR,TIMEND)
      ELSE
        IF (F_DENS) THEN
          write(LUPRI,*) '*** RESTART FROM FITTED DENSITY ON DISK ',
     &                   '- READ FROM DENS_F  ***'
          CALL READ_FIT_DENS(WORK(KCFIT),NBASISAUX,ITER_SCF,.TRUE.)
        ELSE
          CALL DRIVE_FIT (FMAT,DMAT,WORK(KALBE),WORK(KWVEC),NDMT,
     &         ISYMDM,IFCTYP,WORK(KIODBCLS),WORK(KIODBTCH),
     &         WORK(KRODBTCH),
     &         WORK(KGABAB),WORK(KGABAL),WORK(KFREE),LFREE,
     &         WORK(KCFIT), IPRFCK)
        ENDIF

C       Saving of fitted densities turned off. /SR
        IF (.FALSE.) THEN
          CALL SAVE_FIT_DENS(WORK(KCFIT),NBASISAUX,ITER_SCF,.FALSE.)
c          CALL TIMER('FITING',TIMSTR,TIMEND)
        ENDIF

        IF (IPRINT.GT.25) THEN
          CALL HEADER('CFIT',-1)
          CALL OUTPUT(WORK(KCFIT),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
        ENDIF

        CALL SET_ATTGAM(.FALSE.,0.0d0)

c        CALL TIMER('ITSOLV',TIMSTR,TIMEND)

        CALL FITFCK3(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,
     &       WORK(KGABAB),WORK(KGABAL),WORK(KCFIT),IPRFCK,
     &       WORK(KFREE),LFREE)

        CALL TIMER('JBUILD',TIMSTR,TIMEND)
C
#if defined (DEBUG_DIRAC_FCKBUILD_DF) 
        IF (IPRFCK.GE.0) THEN
#else
        IF (IPRFCK.GT.2) THEN
#endif
          CALL HEADER('Fitted Fock matrix',-1)
          CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        ENDIF
      ENDIF
C
C     Release the memory allocated for the Fock matrix build. This
C     will also check for errors due to memory bound failures
C

C     For ROBUST_3I fitting scheme only
      IF (ROBUST_3I) THEN
        CALL MEMREL('DENS_FIT',WORK,KALAB,KALAB,KFREE,LFREE)
      ENDIF
      CALL MEMREL('DENS_FIT',WORK,KGABAL_IT,KGABAL_IT,KFREE,LFREE)
      CALL MEMREL('DENS_FIT',WORK,KWVEC,KWVEC,KFREE,LFREE)
      CALL MEMREL('DENS_FIT',WORK,KGABAL,KGABAL,KFREE,LFREE)
      CALL MEMREL('DENS_FIT',WORK,KGABAB,KGABAB,KFREE,LFREE)
      CALL MEMREL('DENS_FIT',WORK,KALBE,KALBE,KFREE,LFREE)
      CALL MEMREL('DENS_FIT',WORK,KCFIT,KCFIT,KFREE,LFREE)
      CALL QEXIT('DENS_FIT')
      RETURN
      END
C
C
C  /* Deck fitfck3_old */
      SUBROUTINE FITFCK3(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,
     &                  GABAB,GABAL,DGABAL,IPRFCK,
     &                  WORK,LWORK)
C
C   PURPOSE :   Similar to HERFCK, but adjusted for density fitting.
C               Driver routine for direct calculation of the
C               two-electron part of the fock matrices by means
C               of density-fitting.
C               We assume the densities are fitted into a three indeces 
C               vector DMAT(NBASISAUX,1,NDMAT) and that the fock matrix 
C               is a full square matrix. Both matrices without symmetry 
C               reduction.
C
C               SR and FP 15/2-2005.
C
#include "implicit.h"
#include "dummy.h"
      LOGICAL KCLERI
C
      DIMENSION FMAT(*), DMAT(*), ISYMDM(*), IFCTYP(*), 
     &          WORK(*), GABAB(*), GABAL(*), DGABAL(*)
C
C
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
#include "cbifck3.h"
#include "priunit.h"
#include "densfit.h"
C
#include "memint.h"
C
      CALL QENTER('FITFCK3')
c      
#ifdef MOD_FMM
      IF (DOFMM) THEN
        KCLERI = NCLERI
        NCLERI = .TRUE.
        NCLHER = .TRUE.
      END IF
#endif

c     CALL TIMER('START ',TIMSTR,TIMEND)
      DO IDMT=1,NDMT
        ISTR = N2BASX*(IDMT - 1) + 1
        ISTRAUX = NBASISAUX*(IDMT - 1) + 1
C
C aspg, 2006-08-14:
C
C due to the differences in the way dirac and dalton's dft code 
C work at the time of writing (in dirac's setfck_ ifctyp is 
C changed to "coulomb"(=1) whereas in dalton apparently ifcktyp 
C not changed and remains as "coulomb+exchange"(=3) ), the
C if statement below for dirac tests for both the pure coulomb
C and coulomb+exchange cases 
C
#if defined (PRG_DIRAC)
        IF ((MOD(IFCTYP(IDMT),10) .EQ. 1).OR.
     &      (MOD(IFCTYP(IDMT),10) .EQ. 3)) THEN
#else
        IF ( MOD(IFCTYP(IDMT),10) .EQ. 3 ) THEN
#endif
          RHSFIT = .TRUE.
          LHSFIT = .FALSE.
          CALL FCK3DRV(FMAT(ISTR),DMAT(ISTR),1,GABAB,GABAL,
     &                 DGABAL(ISTRAUX),WORK,LWORK)
#ifdef MOD_FMM
          IF (DOFMM) THEN
           CALL FMMFCK(FMAT(ISTR),DMAT(ISTR),WORK,LWORK)
          END IF
#else
C         IF (DOFMM) CALL QUIT('FMM requested but not enabled')
#endif
        END IF
      END DO

#ifdef MOD_FMM
      IF (DOFMM) NCLERI = KCLERI
#endif

      CALL QEXIT('FITFCK3')

c     CALL TIMER('FOCKAB',TIMSTR,TIMEND)
 
      RETURN
      END
C
C  /* Deck drive_fit */
       SUBROUTINE DRIVE_FIT (FMAT,DMAT,WMAT,WVEC,NDMT,ISYMDM,
     &                       IFCTYP,IODBCLS,IODBTCH,RODBTCH,
     &                       GABAB,GABAL,WORK,LWORK,CFIT,IPRINT)
C*****************************************************************************
C    
C     DRIVE_FIT: Driver for density fitting.
C        
C     Written by Luuk Visscher, november 2003
C     Mdified by Simen Reine, september 2005
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "densfit.h" 
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK),
     &          IFCTYP(NDMT), ISYMDM(NDMT),
     &          IODBCLS(*),IODBTCH(*),RODBTCH(*), 
     &          GABAB(*),GABAL(*), CFIT(*)
#include "nuclei.h"
#include "symmet.h"
C
#include "memint.h" 
      CALL QENTER('DRIVE_FIT')
C
C     **************************************************
C     ************* Carry out the fitting ************** 
C     **************************************************
C
      !If requested, loop over symmetry unique center pairs (diatom fit)
      IF (DIATOM) THEN
        CALL QUIT('Not working for DIATOM scheme')
      ELSE
        CALL FIT_DENSITY (FMAT,DMAT,WMAT,WVEC,NDMT,ISYMDM,IFCTYP,
     &       GABAB,GABAL,WORK(KFREE),LFREE,CFIT,
     &       IPRINT)
      END IF
C
      CALL QEXIT('DRIVE_FIT')
      RETURN
      END
C  /* Deck fit_density */
      SUBROUTINE FIT_DENSITY (FMAT,DMAT,WMAT,WVEC,NDMT,ISYMDM,IFCTYP,
     &                        GABAB,GABAL,WORK,LWORK,CFIT,IPRINT)
C*****************************************************************************
C    
C     FIT_DENSITY: Global density fitting procedure (not robust).
C        
C     Written by Luuk Visscher, november 2003
C     Mdified by Simen Reine, april 2005
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK),
     &          IFCTYP(NDMT), ISYMDM(NDMT), 
     &          GABAB(*),GABAL(*), CFIT(*)
#include "nuclei.h"
#include "symmet.h"
#include "inforb.h"
#include "densfit.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
#include "memint.h" 
C
C     We take as first index the fit set functions on either center A or B,
C     and as second index the null function. The remaining two are ordinary
C     basis functions on either center A or B.
C     
      CALL QENTER('FIT_DENSITY')
c     CALL TIMER('START ', TIMSTR, TIMEND)
C
C     Initialize the weight vector and subsequently fill it by
C     contracting the active integrals with density matrix
C
      CALL DZERO (WVEC,NBASISAUX*NDMT)

      CALL GALFCK3(WVEC,DMAT,NDMT,ISYMDM,IFCTYP,IPRINT,
     &     GABAB,GABAL,WORK(KFREE),LFREE)
c     CALL TIMER('GALPHA', TIMSTR, TIMEND)

C     
C     Saving the LHS-vector. Currently turned off \SR
C
      IF (.FALSE.) THEN
        IUNIT = -1
        CALL GPOPEN(IUNIT,'GALPHA','UNKNOWN','SEQUENTIAL','FORMATTED',
     &       IDUM,LDUM)
        WRITE(IUNIT,'(66D27.18)') (WVEC(J), J=1,NBASISAUX)
        CALL GPCLOSE(IUNIT,'KEEP')
      ENDIF

C     
C     Do some printing in case we're really desparate
C
      IF (IPRINT.GT.10) THEN
         CALL HEADER('Full weight vector',-1)
         CALL OUTPUT(WVEC,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
C
C     Solve the equation for the fit coefficients
C
C     Direct Linear equation solver (with Cholesky-decomposition)
      IF (.TRUE.) THEN

      IF (NDMT .GT. 1) THEN
        CALL MEMGET ('REAL',KWMATCP,NBASISAUX*NBASISAUX,WORK,
     &               KFREE,LFREE)
        CALL MEMGET ('REAL',KWVECCP,NBASISAUX*NDMT,WORK,
     &               KFREE,LFREE)
        CALL DCOPY(NBASISAUX*NDMT,WVEC,1,WORK(KWVECCP),1)
        DO IDMT=1,NDMT
          ISTRAUX = NBASISAUX*(IDMT - 1) + 1
c         Need a copy because DPOSV manipulates the matrix
          CALL DCOPY(NBASISAUX*NBASISAUX,WMAT,1,WORK(KWMATCP),1)
          CALL DPOSV ( 'L', NBASISAUX, 1, WORK(KWMATCP), NBASISAUX,
     &                 WORK(KWVECCP+ISTRAUX-1), NBASISAUX, INFO )
        ENDDO
        CALL DCOPY(NBASISAUX*NDMT,WORK(KWVECCP),1,CFIT,1)
        CALL MEMREL('FIT_DENSITY',WORK,KWVECCP,KWVECCP,KFREE,LFREE)
        CALL MEMREL('FIT_DENSITY',WORK,KWMATCP,KWMATCP,KFREE,LFREE)
      ELSE
        CALL MEMGET ('REAL',KWVECCP,NBASISAUX,WORK,
     &               KFREE,LFREE)
        CALL DCOPY(NBASISAUX,WVEC,1,WORK(KWVECCP),1)
        CALL DPOSV ( 'L', NBASISAUX, 1, WMAT, NBASISAUX,
     &               WORK(KWVECCP), NBASISAUX, INFO )
        CALL DCOPY(NBASISAUX,WORK(KWVECCP),1,CFIT,1)
        CALL MEMREL('FIT_DENSITY',WORK,KWVECCP,KWVECCP,KFREE,LFREE)
      ENDIF
c       CALL TIMER('LINSOL', TIMSTR, TIMEND)
C
C     Iterative solver. Tested but slow convergence due to large
C     condition-number. Could still prove useful with a more local
C     constraint of the density residual, i.e. with a more sprase/
C     better precondition (alpha|beta) matrix. Turned off. \SR
      ELSE
C
        CALL MEMGET ('REAL',KWMATCP,NBASISAUX*NBASISAUX,WORK,
     &       KFREE,LFREE)
        CALL DCOPY(NBASISAUX*NBASISAUX,WMAT,1,WORK(KWMATCP),1)

        IF (NDMT .GT. 1) THEN
          CALL QUIT('Not implemented!')
        ELSE
          TOLERANCE = 1.0d-11
          MAXITER   = 100
          CALL MEMGET ('REAL',KWVECCP,NBASISAUX,WORK,
     &               KFREE,LFREE)
C         INITIAL GUESS
          call asolve(GABAL,NBASISAUX,WVEC,WORK(KWVECCP))

          CALL linbcg(NBASISAUX,WVEC,WORK(KWVECCP),2,TOLERANCE,
     &         MAXITER,ITER,ERROR,
     &         NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
C     With the full (al|be) matrix
     &         WORK(KWMATCP),LFREE)
C     With the direct product (al|be)c_be via the fck3drv
c     &         WORK(KFREE),LFREE)

          CALL MEMREL('FIT_DENSITY',WORK,KWVECCP,KWVECCP,KFREE,LFREE)
        ENDIF
        WRITE(LUPRI,*) 'Iterative solver for fitting coefficients ',
     &                 'used ',ITER,' iterations; error is ', ERROR
        CALL MEMREL('FIT_DENSITY',WORK,KWMATCP,KWMATCP,KFREE,LFREE)
c       CALL TIMER('ITLES ', TIMSTR, TIMEND)
      ENDIF
C

C
C     We again have some information that may be useful
C
      IF (IPRINT.GT.4) THEN
         WRITE (LUPRI,'(//8X,A,2I6//)') 'INFO code from DPOSV',INFO
         CALL HEADER('Fit coefficients',-1)
         CALL OUTPUT(CFIT,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
C
      CALL QEXIT('FIT_DENSITY')
      RETURN
      END
C
C     Gets galpha dependent on the scheme chosen
C
      SUBROUTINE DF_GET_GALPHA(GALPHA,NDIMAUX,DMAT,NDMT,ISYMDM,IFCTYP,
     &                         GABAB,GABAL,I_LG,IPRINT,WORK,LWORK)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "densfit.h" 
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
      DIMENSION GALPHA(NDIMAUX*NDMT),DMAT(*),WORK(LWORK),ISYMDM(*),
     &          IFCTYP(NDMT),GABAB(*),GABAL(*), I_LG(*)
#include "nuclei.h"
#include "symmet.h"
      LOGICAL SPEXP_SAVE,DOFMM_SAVE
C
#include "memint.h" 
C
      CALL QENTER('GET_GALP')
C
      IF (.NOT.DOBOX) THEN
        CALL ACT_ALL_BATCHES
        CALL DZERO (GALPHA,NDIMAUX*NDMT)
        CALL GALFCK3(GALPHA,DMAT,NDMT,ISYMDM,IFCTYP,IPRINT,
     &               GABAB,GABAL,WORK(KFREE),LFREE)
      ELSE
        CALL MEMGET('REAL',KGALFULL,NBASISAUX*NDMT,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KGALFULL),NBASISAUX*NDMT)
        SPEXP_SAVE = SPEXP
        DOFMM_SAVE = DOFMM
        SPEXP = .FALSE.
        DOFMM = .FALSE.   
        CALL GALFCK3(WORK(KGALFULL),DMAT,NDMT,ISYMDM,IFCTYP,IPRINT,
     &               GABAB,GABAL,WORK(KFREE),LFREE)
        SPEXP = SPEXP_SAVE
        DOFMM = DOFMM_SAVE
        CALL PROJECT_VEC(WORK(KGALFULL),NBASISAUX,GALPHA,
     &                   NDIMAUX,I_LG,NDMT,.TRUE.,IPRINT)
        CALL MEMREL('GET_GALP',WORK,KGALFULL,KGALFULL,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT('GET_GALP')
C
      RETURN
      END
C
C
C
      FUNCTION DF_GET_CHARGE (DMAT,WORK,LWORK)
#include "implicit.h"
#include "inforb.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "nuclei.h"
#include "fck3aob.h"
#include "densfit.h"
#include "cbifck3.h"
      DIMENSION DMAT(*), KCCF(2,2), WORK(*)
      LOGICAL LFOUND, ACTIVE1, ACTIVE2
      PARAMETER (DP0 = 0.0d0, DP5 = 0.5d0, DP1 = 1.0d0 )
#include "memint.h"
C
      CALL QENTER('GET_CHAR')
C     Get overlap-matrix
c     CALL SETAOBTCH(SPEXP,CONTRACTED,WORK(KFREE),LFREE,0)
      CALL MEMGET('REAL',KSMAT,NBASIS*NBASIS,WORK,KFREE,LFREE)
#if defined PRG_DIRAC
      CALL GTOVLX(WORK(KSMAT),DP1)
#else
      CALL MEMGET('REAL',KSMTMP,NBASIS*(NBASIS+1)/2,WORK,KFREE,LFREE)
      CALL RDONEL('OVERLAP ',LFOUND,WORK(KSMTMP),NBASIS*(NBASIS+1)/2)
      IF (.NOT.LFOUND) CALL QUIT('OVERLAP integrals not found on disk'//
     &                          ' in routine DF_GET_CHARGE')
      CALL DSPTSI(NBASIS,WORK(KSMTMP),WORK(KSMAT))
      CALL MEMREL('GET_CHAR',WORK,KSMTMP,KSMTMP,KFREE,LFREE)
#endif
      IF (DOBOX) THEN
        DCHARGE = DP0
        DO IAOB1=1,NAOBTCH
          LCLASS = KCLSBTCH(IAOB1)
          ACTIVE1 = ACTVAOBTCH(IAOB1) 
          IF (LCLASS.EQ.1) THEN
            CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,J1,KHKT1,KCKT1,
     &                  KCSQ1,KEXP1,KCCF(1,1),IAOB1)
            DO IAOB2=1,NAOBTCH
              LCLASS = KCLSBTCH(IAOB2)
              ACTIVE2 = ACTVAOBTCH(IAOB2) 
              IF ((LCLASS.EQ.1).AND.(ACTIVE1.OR.ACTIVE2)) THEN
                CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,J2,KHKT2,KCKT2,
     &                      KCSQ2,KEXP2,KCCF(1,2),IAOB2)
                WEIGHT = DP5
                IF (ACTIVE1.AND.ACTIVE2) WEIGHT = DP1
                DO ICOMP1=1,NTYP1
                  IF ( NTYP1 .EQ. 2 ) 
     &               CALL SETJKKK(J1,KCSQ1,KCKT1,KHKT1,ICOMP1)
                  DO ICOMP2=1,NTYP2
                    IF ( NTYP2 .EQ. 2 ) 
     &               CALL SETJKKK(J2,KCSQ2,KCKT2,KHKT2,ICOMP2)
                     DCHARGE = DCHARGE + DF_CHARGE_BLOCK(DMAT,
     &                                   WORK(KSMAT),WEIGHT,
     &                                   KHKT1,NCTF1,ICOMP1,IAOB1,
     &                                   KHKT2,NCTF2,ICOMP2,IAOB2)
                  ENDDO   !ICOMP2
                ENDDO   !ICOMP1
              ENDIF   !LCLASS&ACTIVE
            ENDDO   !IAOB2
          ENDIF   !LCLASS
        ENDDO   !IAOB1
        DF_GET_CHARGE = DCHARGE
      ELSE
        DF_GET_CHARGE = TRACE_MAT_PRO(DMAT,WORK(KSMAT),NBASIS,NBASIS)
      ENDIF
      CALL MEMREL('GER_CHAR',WORK,KSMAT,KSMAT,KFREE,LFREE)
C
      CALL QEXIT('GET_CHAR')
C
      RETURN
      END

      FUNCTION TRACE_MAT_PRO(A,B,N1,N2)
#include "implicit.h"
      DIMENSION A(N1,N2),B(N2,N1)
      TRACE = 0.0d0
      DO I=1,N1
        DO J=1,N2
          TRACE = TRACE + A(I,J)*B(J,I)
        ENDDO
      ENDDO
      TRACE_MAT_PRO = TRACE
      RETURN
      END


      FUNCTION DF_CHARGE_BLOCK(DMAT,SMAT,WEIGHT,
     &                         KHKT1,NCTF1,ISPCOMP1,IAOB1,
     &                         KHKT2,NCTF2,ISPCOMP2,IAOB2)
#include "implicit.h"
#include "inforb.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "nuclei.h"
#include "fck3aob.h"
      DIMENSION DMAT(NBASIS,NBASIS), SMAT(NBASIS,NBASIS) 
      Q = 0.0d0
      DO ICTF1 = 1, NCTF1
        ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1) + ICTF1 - 1, 1 )
        DO ICTF2 = 1, NCTF2
          ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1) + ICTF2 - 1, 1 )
          DO ISP1=1, KHKT1
           IAO1 = ISTAO1+ISP1-1
           DO ISP2=1, KHKT2
            IAO2 = ISTAO2 + ISP2 - 1
            Q = Q + DMAT(IAO1,IAO2) * SMAT(IAO2,IAO1)
           ENDDO
          ENDDO
        ENDDO
      ENDDO
      DF_CHARGE_BLOCK = Q*WEIGHT
      RETURN
      END

      SUBROUTINE DF_GET_VMAT(VMAT,I_LG,NDIMAUX,WORK,LWORK,IPRINT,
C     Only temporary
     &                       IBOX,NUMBOX,CORDBOX)
C     Only temporary
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "iteration.h"
      DIMENSION VMAT(*), I_LG(*), WORK(LWORK),
C     Only temporary
     &          NUMBOX(*), CORDBOX(*)
C     Only temporary
#include "memint.h"
      CALL QENTER('GET_VMAT')
C
c     IF (DOBOX) THEN
C       Only temporary
c       CALL MEMGET('INTEGER',KINDBOX,NBASISAUX,WORK,KFREE,LFREE)
C       Only temporary
c       CALL MEMGET('REAL',KVMATBG,NBASISAUX*NBASISAUX,WORK,KFREE,LFREE)
c       CALL GTALBE(WORK(KVMATBG),I_LG,NDIMAUX,WORK(KFREE),LFREE,IPRINT)
c       CALL PROJECT_MAT(WORK(KVMATBG),NBASISAUX,VMAT,
c    &                   NDIMAUX,I_LG,.TRUE.,IPRINT)

c       CALL MEMREL('GET_VMAT',WORK,KVMATBG,KVMATBG,KFREE,LFREE)
c       CALL MEMREL('GET_VMAT',WORK,KINDBOX,KINDBOX,KFREE,LFREE)
c     ELSE
        CALL GTALBE(VMAT,I_LG,NDIMAUX,WORK(KFREE),LFREE,IPRINT)
c     ENDIF
C
      CALL QEXIT('GET_VMAT')
C
      RETURN
      END
C  /* Deck df_gendrv */
      SUBROUTINE DF_GENDRV (CFIT,DMAT,NDIMAUX,NDMT,I_LG,ISYMDM,IFCTYP,
     &                       GABAB,GABAL,IBOX,
     &                       WORK,LWORK,IPRINT)
C*****************************************************************************
C    
C     DF_GENDRV: General driver for density fitting.
C        
C     Written by Erik Tellgren and Simen Reine april 2006
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "densfit.h" 
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
      DIMENSION DMAT(*),WORK(LWORK),ISYMDM(*),IFCTYP(NDMT), 
     &          GABAB(*),GABAL(*), CFIT(NDIMAUX,NDMT), I_LG(*)
#include "nuclei.h"
#include "symmet.h"
C
#include "memint.h" 
C
      CALL QENTER('DFGENDRV')
      CALL MEMGET('REAL',KGALPHA,NDIMAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAUXINT,NDIMAUX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMAT,NDIMAUX*NDIMAUX,WORK,KFREE,LFREE)

Caspg+rbast
      IF (DF_FITVXC) THEN
         CALL MEMGET('REAL',KGJK,NDIMAUX*NDIMAUX,WORK,KFREE,LFREE)
      ENDIF
Caspg+rbast

C
      CALL DZERO (CFIT,NDIMAUX*NDMT)
C     GALPHA = sum_cd (alpha|cd)Dcd
      CALL DF_GET_GALPHA(WORK(KGALPHA),NDIMAUX,DMAT,NDMT,ISYMDM,IFCTYP,
     &                   GABAB,GABAL,I_LG,IPRINT,WORK(KFREE),LFREE)
C
C     Get the auxilliary 2-center matrix
      CALL DF_GET_VMAT(WORK(KVMAT),I_LG,NDIMAUX,
     &                 WORK(KFREE),LFREE,IPRINT,
C     Only temporary
     &                 IBOX,NUMBOX,CORDBOX)
C
c     CALL HEADER('VMAT in DF_GENDRV',-1)
c     CALL OUTPUT(WORK(KVMAT),1,NDIMAUX,1,NDIMAUX,NDIMAUX,NDIMAUX,1,2)
c     CALL HEADER('GALPHA in DF_GENDRV',-1)
c     CALL OUTPUT(WORK(KGALPHA),1,NDIMAUX,1,1,NDIMAUX,1,1,2)
C

Caspg+rbast, 2006-09-21 
      IF (SAVE_VAB) THEN
         CALL SAVE_MATRIX(WORK(KVMAT),NBASISAUX,NBASISAUX,'ALBEMA')
         IF (.false.) THEN
            CALL HEADER('ccVMAT in DF_GENDRV',-1)
            CALL OUTPUT(WORK(KVMAT),1,NBASISAUX,1,NBASISAUX,
     &               NBASISAUX,NBASISAUX,1,2)
         ENDIF
      ENDIF
      IF (DF_FITVXC) THEN
C to avoid having to setup another KVMAT for the XC-Fit stuff at
C the end, we now always set the "keep matrix" variable in linsol
C set to true
         CALL LINSOL(WORK(KVMAT),CFIT,WORK(KGALPHA),NDIMAUX,NDMT,
     &               .TRUE.,WORK(KFREE),LFREE)
      ELSE
         CALL LINSOL(WORK(KVMAT),CFIT,WORK(KGALPHA),NDIMAUX,NDMT,
     &               DF_CONSTR,WORK(KFREE),LFREE)
      ENDIF
      IF (SAVE_CFIT) THEN
          if (DODFLR) then
          write (lupri,*) 'debug: saving intermediate fit coefs for lr'
          CALL SAVE_MATRIX(CFIT,1,NBASISAUX,'LRDFC1')
          else
          CALL SAVE_MATRIX(CFIT,1,NBASISAUX,'FITCO1')
          endif
      ENDIF
Caspg+rbast, 2006-09-21

C
      IF (DF_CONSTR) THEN
C       Grenerate the one-electron integrals over the fitting basis
        CALL DF_GET_AUXCHA_LOC(WORK(KAUXINT),I_LG,NDIMAUX,
     &                         WORK(KFREE),LFREE,IPRINT)
        CALL MEMGET('REAL',KCORR,NDIMAUX,WORK,KFREE,LFREE)
C

Caspg+rbast, 2006-09-21
C
         IF (DF_FITVXC) THEN
            CALL LINSOL(WORK(KVMAT),WORK(KCORR),WORK(KAUXINT),NDIMAUX,1,
     &                 .TRUE.,WORK(KFREE),LFREE)
         ELSE
            CALL LINSOL(WORK(KVMAT),WORK(KCORR),WORK(KAUXINT),NDIMAUX,1,
     &                 .FALSE.,WORK(KFREE),LFREE)
         ENDIF 
C if we are asked to save the fit coefficients, those for the charge
C conservation constraint are also saved, but in a different file
         IF (SAVE_CFIT) THEN
            if (DODFLR) then
            write (*,*) 'debug: saving second fit coefs for lr'
            CALL SAVE_MATRIX(CFIT,1,NBASISAUX,'LRDFCH')
            else
            CALL SAVE_MATRIX(WORK(KCORR),1,NBASISAUX,'FITCHC')
            endif
         ENDIF
Cdebug
         if (.false.) then
          CALL HEADER('Constrained fit coefs ',-1)
          CALL OUTPUT(WORK(KCORR),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
        ENDIF
Cdebug
Caspg+rbast, 2006-09-21

C
        WRITE(LUPRI,'(A,I4)') 'Charge of density/fitted '//
     &                         'density of IBOX', IBOX
        WRITE(LUPRI,'(A)') ' IDMT      Normal           Free fit    '//
     &                     '  Constrained fit    Max cfit change'
        DO IDMT=1,NDMT
          BOX_CHARGE = DF_GET_CHARGE(DMAT,WORK(KFREE),LFREE)
          CHARGE_FIT =
     &         DDOT(NDIMAUX,WORK(KAUXINT),1,CFIT(1,IDMT),1) 
          DLAMBDA = (BOX_CHARGE - 
     &         DDOT(NDIMAUX,WORK(KAUXINT),1,CFIT(1,IDMT),1)) /
     &         DDOT(NDIMAUX,WORK(KAUXINT),1,WORK(KCORR),1)
          CALL DAXPY(NDIMAUX,DLAMBDA,WORK(KCORR),1,CFIT(1,IDMT),1)
          CORRECTED_CHARGE =
     &         DDOT(NDIMAUX,WORK(KAUXINT),1,CFIT(1,IDMT),1) 
          COEFF_MAX = 0.0d0
          DO I=1,NDIMAUX
           COEFF_MAX = MAX(COEFF_MAX,ABS(WORK(KCORR+I-1)))
          ENDDO
          COEFF_MAX = COEFF_MAX * ABS(DLAMBDA)
          WRITE(LUPRI,'(I3,4F18.12)') IDMT,BOX_CHARGE,CHARGE_FIT,
     &                                CORRECTED_CHARGE,COEFF_MAX
        ENDDO
Caspg+rbast, 2006-09-21
C setting up Gjk matrix for the XC contribution in DFT
        IF (DF_FITVXC) THEN
C          write (6,*) 'debug: ook'
           CALL DF_FitXC_Build_GKJ(WORK(KGJK),NBASISAUX,NBASISAUX,
     &                           WORK(KVMAT),WORK(KAUXINT),WORK(KCORR))
           CALL SAVE_MATRIX(WORK(KGJK),NBASISAUX,NBASISAUX,'FITXCG')
        ENDIF

Caspg+rbast, 2006-09-21
        CALL MEMREL('DFGENDRV',WORK,KCORR,KCORR,KFREE,LFREE)
      ELSE
        CALL DF_GET_AUXCHA_LOC(WORK(KAUXINT),I_LG,NDIMAUX,
     &                         WORK(KFREE),LFREE,IPRINT)
        WRITE(LUPRI,'(A,I4)') 'Charge of density/fitted '//
     &                         'density of IBOX', IBOX
        WRITE(LUPRI,'(A)') ' IDMT      Normal           Free fit    '
        DO IDMT=1,NDMT
          BOX_CHARGE = DF_GET_CHARGE(DMAT,WORK(KFREE),LFREE)
          CHARGE_FIT =
     &         DDOT(NDIMAUX,WORK(KAUXINT),1,CFIT(1,IDMT),1) 
          WRITE(LUPRI,'(I3,4F18.12)') IDMT,BOX_CHARGE,CHARGE_FIT
        ENDDO

Caspg+rbast, 2006-09-21
C setting up Gjk matrix for the XC contribution in DFT
        IF (DF_FITVXC) THEN
           CALL DF_FitXC_Build_GKJ(WORK(KGJK),NBASISAUX,NBASISAUX,
     &                           WORK(KVMAT),DUMMY,DUMMY)
           CALL SAVE_MATRIX(WORK(KGJK),NBASISAUX,NBASISAUX,'FITXCG')
        ENDIF
Caspg+rbast, 2006-09-21

      ENDIF
Caspg, 2006-09-21
C fit coefficients for this iteration are saved to file
      write (6,*) 'goo'
      IF (SAVE_CFIT) THEN
          write (6,*) 'foo'
          if (DODFLR) then
          write (*,*) 'debug: saving final fit coefs for lr'
          CALL SAVE_MATRIX(CFIT,1,NBASISAUX,'LRDFCE')
          else
          CALL SAVE_MATRIX(CFIT,1,NBASISAUX,'FITCOE')
          endif
      ENDIF
      IF (DF_FITVXC) THEN
         CALL MEMREL('dfgendrv-kgjk',WORK,KGJK,KGJK,KFREE,LFREE)
      ENDIF
Caspg, 2006-09-21
      CALL MEMREL('DFGENDRV',WORK,KVMAT,KVMAT,KFREE,LFREE)
      CALL MEMREL('DFGENDRV',WORK,KAUXINT,KAUXINT,KFREE,LFREE)
      CALL MEMREL('DFGENDRV',WORK,KGALPHA,KGALPHA,KFREE,LFREE)
      CALL QEXIT('DFGENDRV')
      RETURN
      END
C  /* Deck galfck3 */
      SUBROUTINE GALFCK3(CVEC,DMAT,NDMT,ISYMDM,IFCTYP,IPRINT,
     &                   GABAB,GABAL,WORK,LWORK)
C
C   PURPOSE :   Driver routine for direct calculation of the
C               gamma_alpha (GAL) vector by means of FCK3. The GAL-
C               vector it the right-hand side of the equation for 
C               the density-fitting coefficients:
C
C                    SUM_beta{(alpha|beta) c_beta} = gamma_alpha,
C
C              where
C
C                    gamma_alpha = SUM_ab{(ab|alpha) Dab}.
C
C               DMAT is the density matrix (in MO-representation), 
C               CVEC returns the vector with elements gamma_alpha.
C
C               SR and FP 18/2-2005.
C
#include "implicit.h"
      LOGICAL KCLERI
C
      DIMENSION CVEC(*), DMAT(*), ISYMDM(*), IFCTYP(*), 
     &          GABAB(*),GABAL(*),WORK(*)
C
C
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
#include "cbifck3.h"
#include "priunit.h"
#include "densfit.h"
C
#include "memint.h"
C
      CALL QENTER('GALFCK3')
      
      DO IDMT=1,NDMT
        ISTR = N2BASX*(IDMT - 1) + 1
        ISTRAUX = NBASISAUX*(IDMT - 1) + 1
C
C aspg, 2006-08-14:
C
C due to the differences in the way dirac and dalton's dft code 
C work at the time of writing (in dirac's setfck_ ifctyp is 
C changed to "coulomb"(=1) whereas in dalton apparently ifcktyp 
C not changed and remains as "coulomb+exchange"(=3) ), the
C if statement below for dirac tests for both the pure coulomb
C and coulomb+exchange cases 
C
#if defined (PRG_DIRAC)
        IF ((MOD(IFCTYP(IDMT),10) .EQ. 1).OR.
     &      (MOD(IFCTYP(IDMT),10) .EQ. 3)) THEN
#else
        IF ( MOD(IFCTYP(IDMT),10) .EQ. 3 ) THEN
#endif
Caspg
C         RHSFIT build of Fab (Fitted 2nd electron part of Fock-matrix)
C         LHSFIT build of gamma_alpha (GAL) vector. 
          RHSFIT = .FALSE.
          LHSFIT = .TRUE.
          CALL FCK3DRV(CVEC(ISTRAUX),DMAT(ISTR),1,GABAB,GABAL,
     &                DUMMY, WORK,LWORK)
#ifdef MOD_FMM
          IF (DOFMM) THEN
            CALL FMMFCK(CVEC(ISTRAUX),DMAT(ISTR),WORK,LWORK)
          END IF
#else
          IF (DOFMM) CALL QUIT('FMM requested but not enabled')
#endif
        END IF
      END DO
C     
C     Do some printing in case we're really desparate
C
      IF (IPRINT.GT.20) THEN
         CALL HEADER('Full G_alpha vector',-1)
         CALL OUTPUT(CVEC,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF

      CALL QEXIT('GALFCK3')
      RETURN
      END



      SUBROUTINE ROBUST_3I_DRV(FMAT,DMAT,ALBE,NAUX,ABAL,NA,NB,
     &                      WORK,LWORK,IPRINT)
C
C     Driver for ROBUST_3I fitting (only Coulomb integrals)
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMAT(*), DMAT(*), ALBE(*), ABAL(*),WORK(*)
C
#include "memint.h"
C
      IF (IPRINT.GT.10) THEN
        CALL HEADER('ROBUST_3I_DRV',-1)
        WRITE(LUPRI,*) 'NA, NB, NAUX ', NA, NB, NAUX
      ENDIF
      IF (IPRINT.GT.30) THEN
        CALL HEADER('FMAT into ROBUST_3I_DRV',-1)
        CALL OUTPUT(FMAT,1,NA,1,NB,NA,NB,1,LUPRI)
      ENDIF

      CALL MEMGET ('REAL',KCABAL,NAUX*NA*NB,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KCALPH,NAUX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KALBEI,NAUX*NAUX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KIDUMY,NAUX*NAUX*2+NAUX,WORK,KFREE,LFREE)

      CALL DCOPY(NAUX*NAUX,ALBE,1,WORK(KALBEI),1)
#ifdef MOD_DNF
      CALL INVERT_BY_DIAG2(WORK(KALBEI), WORK(KIDUMY),
     &                     WORK(KIDUMY+NAUX*NAUX),
     &                     WORK(KIDUMY+2*NAUX*NAUX),NAUX)
#endif

      IF (IPRINT.GT.30) THEN
        CALL HEADER('(alpha|beta)',-1)
        CALL OUTPUT(ALBE,1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
        CALL HEADER('(alpha|beta) inverse',-1)
        CALL OUTPUT(WORK(KALBEI),1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
      ENDIF

      CALL FAB_ROBUST_3I(FMAT,DMAT,ALBE,WORK(KALBEI),NAUX,ABAL,NA,NB,
     &               WORK(KCABAL),WORK(KCALPH),
     &               WORK(KIDUMY),WORK(KFREE),LFREE,IPRINT)

      CALL MEMREL('ROBUST_3I_DRV',WORK,KIDUMY,KIDUMY,KFREE,LFREE)
      CALL MEMREL('ROBUST_3I_DRV',WORK,KALBEI,KALBEI,KFREE,LFREE)
      CALL MEMREL('ROBUST_3I_DRV',WORK,KCALPH,KCALPH,KFREE,LFREE)
      CALL MEMREL('ROBUST_3I_DRV',WORK,KCABAL,KCABAL,KFREE,LFREE)
      RETURN
      END


      SUBROUTINE FAB_ROBUST_3I(FMAT,DMAT,ALBE,ALBEIN,NAUX,ABAL,NA,NB,
     &                     CABAL,CALPHA,RHOALP,WORK,LWORK,IPRINT)
C
C     Builds Fock-matrix by ROBUST_3I fitting
C
C        F_ab^ROBUST_3I = sum_alpha (ab|alpha) c_alpha
C                        + sum_alpha c_alpha^ab rho_alpha
C
C     with
C
C        c_alpha_cd = sum_beta (alpha|beta)^-1 (beta|cd),
C
C        c_alpha = sum_cd c_alpha^cd D_cd
C
C     and
C
C        rho_alpha = sum_cd (alpha|cd) D_cd
C                      - sum_beta (alpha|beta) c_alpha.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMAT(NA,NB), DMAT(NA,NB), ALBE(NAUX,NAUX), 
     &          ALBEIN(NAUX,NAUX), ABAL(NA,NB,NAUX), 
     &          CABAL(NA,NB,NAUX), CALPHA(NAUX), 
     &          RHOALP(NAUX), WORK(*)
C
#include "memint.h"
      IF (IPRINT.GT.20) THEN
        CALL HEADER('FMAT when entering FAB_ROBUST_3I',-1)
        CALL OUTPUT(FMAT,1,NA,1,NB,NA,NB,1,LUPRI)
      ENDIF
      DO IAUX1 = 1, NAUX
        CALPHA(IAUX1) = 0.d0
        RHOALP(IAUX1) = 0.d0
        DO IA = 1, NA
          DO IB = 1, NB
            CABAL(IA,IB,IAUX1) = ALBEIN(IAUX1,1)*ABAL(IA,IB,1)
            DO IAUX2 = 2, NAUX
              CABAL(IA,IB,IAUX1) = CABAL(IA,IB,IAUX1) + 
     &             ALBEIN(IAUX1,IAUX2)*ABAL(IA,IB,IAUX2)
            ENDDO
            CALPHA(IAUX1) = CALPHA(IAUX1) + 
     &           CABAL(IA,IB,IAUX1) * DMAT(IA,IB)
            RHOALP(IAUX1) = RHOALP(IAUX1) + 
     &           ABAL(IA,IB,IAUX1) * DMAT(IA,IB)
          ENDDO
        ENDDO
      ENDDO
      DO IAUX1 = 1, NAUX
        DO IAUX2 = 1, NAUX
          RHOALP(IAUX1) = RHOALP(IAUX1) - 
     &         ALBE(IAUX1,IAUX2) * CALPHA(IAUX2)
        ENDDO
      ENDDO
      DO IA = 1, NA
        DO IB = 1, NB
          DO IAUX = 1, NAUX
            FMAT(IA,IB) = FMAT(IA,IB) +
     &           ABAL(IA,IB,IAUX) * CALPHA(IAUX) +
     &           CABAL(IA,IB,IAUX) * RHOALP(IAUX)
          ENDDO
        ENDDO
      ENDDO

      IF (IPRINT.GT.30) THEN
        CALL HEADER('Fit coefficients',-1)
        WRITE(LUPRI,*) 'NAUX, NA, NB ', NAUX, NA, NB
        CALL OUTPUT(CABAL,1,NAUX,1,NA*NB,NAUX,NA*NB,1,LUPRI)
        CALL HEADER('Contracted fit coefficients',-1)
        CALL OUTPUT(CALPHA,1,NAUX,1,1,NAUX,1,1,LUPRI)
        CALL HEADER('Rho_alpha',-1)
        CALL OUTPUT(RHOALP,1,NAUX,1,1,NAUX,1,1,LUPRI)
      ENDIF
      IF (IPRINT.GT.20) THEN
        CALL HEADER('FMAT when leaving FAB_ROBUST_3I',-1)
        CALL OUTPUT(FMAT,1,NA,1,NB,NA,NB,1,LUPRI)
      ENDIF
      RETURN
      END

C
C The following four subroutines are not used in the current
C version of code
C
C  /* Deck pick_act_batch */
      SUBROUTINE PICK_ACT_BATCH (NCENTS,ICENTS,ICLASSES,IPRINT)
C*****************************************************************************
C    
C     PICK_AOAB: Define active batches in density fitting call to ERI.
C        
C     Written by Luuk Visscher, october 2004
C     
C     Arrays needed to pick the active batches:
C     - NCENTS :   the number of active centers for a given index (maximally 2)
C     - ICENTS :   the corresponding center indices 
C     - ICLASSES : the active class (can only be one at a time).
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
      DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)
C
      DO I = 1, NAOBTCH
        ICENT = NCNTBTCH(I)
        ICLASS = KCLSBTCH(I)
        ACTVAOBTCH(I) = .FALSE.
        DO J = 1, 4
C         The shell must be of the right class
          IF (ICLASS.EQ.ICLASSES(J)) THEN
C         The center must be active (ncents=0 means all centers active)
            IF (NCENTS(J).EQ.0) THEN
              ACTVBTH(I,J) = .TRUE.
            ELSEIF (NCENTS(J).EQ.1.OR.NCENTS(J).EQ.2) THEN
              ACTVBTH(I,J) = .FALSE.
              DO K = 1, NCENTS(J)
                IF (ICENTS(K,J).EQ.ICENT) ACTVBTH(I,J) = .TRUE.
              ENDDO
            ELSE
              CALL QUIT ("Wrong value for NCENTS in PICK_ACT_BATCH")
            ENDIF
          ELSE
            ACTVBTH(I,J) = .FALSE.
          ENDIF
          IF (ACTVBTH(I,J)) ACTVAOBTCH(I) = .TRUE.
        END DO
      END DO
C     
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from PICK_ACT_BATCH',-1)
         DO J = 1, 4
             IF (NCENTS(J).EQ.0) THEN
                WRITE (LUPRI,'(1X,A,I1,A)')
     &         ' Index ',J,', all centers active'
             ELSE
                WRITE (LUPRI,'(1X,A,I1,A,2I5)')
     &         ' Index ',J,', active centers :',
     &           (ICENTS(K,J),K=1,NCENTS(J))
             ENDIF
         ENDDO
         WRITE (LUPRI,'(1X,A,4I5,/)')
     &      ' Active classes for each index:',
     &      (ICLASSES(I),I=1,4)
         WRITE (LUPRI,'(1X,A,/)') ' ACTVBTH in PICKAO:'
         DO I = 1, NAOBTCH
            WRITE (LUPRI,'(10X,I5,5X,4L5)') I, (ACTVBTH(I,J),J=1,4)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck nfun_act_batch */
      SUBROUTINE NFUN_ACT_BATCH (IND,N_FIT,INDEX)
C*****************************************************************************
C    
C     Counts number of functions in active batches for index IND
C     Sets up index array that can be used in gather/scatter
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "fck3aob.h"
C
      DIMENSION INDEX(*)
C
      N_FIT = 0
      J = 0
      DO ISH = 1, NAOBTCH
        IF (ACTVBTH(ISH,IND)) THEN
          NHK =  NHKTBTCH(ISH)
          KHK = NHK * 2 - 1
          NORBBTH = KHK * NCTFBTCH(ISH)
          INDX1   = KAOBTCH(1,ISH) - 1
          DO I = 1, NORBBTH
            N_FIT = N_FIT + 1
            INDEX(N_FIT) = INDX1 + I
          ENDDO
          NTYP = NTYPBTCH(ISH)
          IF ( NTYP .EQ. 2 ) THEN
            INDX2   = KAOBTCH(2,ISH) - 1
            DO I = 1, NORBBTH
              N_FIT = N_FIT + 1
              INDEX(N_FIT) = INDX2 + I
            ENDDO
          ENDIF
        ENDIF
c        DO I = 1, NORBBTH
c           J = J + 1
c           IF (ACTVBTH(ISH,IND)) THEN
c              N_FIT = N_FIT + 1
c              INDEX(N_FIT) = J
c           ENDIF
c        ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck extr_act_batch */
      SUBROUTINE EXTR_ACT_BATCH (FULL_MAT,N_FULL,RED_MAT,N_RED,
     &                           INDA,NA,INDB,NB,MODE)
C*****************************************************************************
C    
C     EXTR_ACTB: Extract/insert the active parts out of /in to a full matrix.
C     MODE 1 : Extraction
C          2 : Insertion
C          3 : Addition of reduced matrix to full matrix
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION FULL_MAT(N_FULL,*),RED_MAT(N_RED,*)
      DIMENSION INDA(*),INDB(*)
C
      IF (MODE.EQ.1) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               RED_MAT(I,J) = FULL_MAT(IF,JF)
            ENDDO
         ENDDO
      ELSEIF (MODE.EQ.2) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               FULL_MAT(IF,JF) = RED_MAT(I,J)
            ENDDO
         ENDDO
      ELSEIF (MODE.EQ.3) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               FULL_MAT(IF,JF) = FULL_MAT(IF,JF) + RED_MAT(I,J)
            ENDDO
         ENDDO
      ELSE
         CALL QUIT ('Wrong MODE in EXTR_ACT_BATCH')
      ENDIF
C     
      RETURN
      END


      SUBROUTINE SAVE_CONV_FIT_DENS
#include "implicit.h"
#include "iteration.h"
#include "densfit.h"
      IF (DOFIT) THEN
#ifdef MOD_LINEARS
        CALL di_SAVE_CONV_FIT_DENS(NFITDIM,ITER_SCF)
#else
        write (*,*) 'a replacement of di_save_conv_fit ',
     &              'needs to be implemented'
#endif
      ENDIF
      RETURN
      END

      SUBROUTINE SAVE_FIT_DENS(CFIT,NFIT,ITER,LAST)
#include "implicit.h"
#include "iteration.h"
      DIMENSION CFIT(NFIT)
      CHARACTER*6 FILENAME
      LOGICAL LAST

      NFITDIM = NFIT
      WRITE(FILENAME,'(A4,I2)') 'DENS', ITER
      IF (LAST) THEN
        FILENAME='DENS_F'
      ENDIF
      IUNIT = -1
      CALL GPOPEN(IUNIT,FILENAME,'UNKNOWN','SEQUENTIAL','FORMATTED',
     &            IDUM,LDUM)
      WRITE(IUNIT,'(I8)') NFIT
      DO I=1,NFIT
        WRITE(IUNIT,'(F18.6)') CFIT(I)
      ENDDO
      CALL GPCLOSE(IUNIT,'KEEP')

      RETURN
      END

      SUBROUTINE READ_FIT_DENS(CFIT,NFIT,ITER,FIRST)
#include "implicit.h"
      DIMENSION CFIT(NFIT)
      LOGICAL FIRST
      CHARACTER*6 FILENAME

      WRITE(FILENAME,'(A4,I2)') 'DENS', ITER
      IF (FIRST) THEN
        FILENAME='DENS_F'
      ENDIF
      IUNIT = -1
      CALL GPOPEN(IUNIT,FILENAME,'UNKNOWN','SEQUENTIAL','FORMATTED',
     &            IDUM,LDUM)
      READ(IUNIT,'(I8)') IFIT
      IF (IFIT .NE. NFIT) THEN
        CALL QUIT('IFIT and NFIT do not match in call READ_FIT_DENS')
      ELSE
        DO I=1,NFIT
          READ(IUNIT,'(F18.6)') CFIT(I)
        ENDDO
        CALL GPCLOSE(IUNIT,'KEEP')
      ENDIF
      RETURN
      END

      SUBROUTINE SET_ITER(ITER)
#include "implicit.h"
#include "iteration.h"
      ITER_SCF = ITER
      RETURN
      END

      SUBROUTINE SET_ATTGAM(ATTSET,OMEGA)
#include "implicit.h"
#include "cbieri.h"
#include "ericom.h"
      LOGICAL ATTSET
      ATTGAM = ATTSET
      OMEGA2 = OMEGA*OMEGA
      RETURN
      END

      LOGICAL FUNCTION GET_F_DENS()
#include "iteration.h"
      GET_F_DENS = F_DENS
      RETURN
      END

      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      DOUBLE PRECISION sx(n),snrm
      if (itol.le.3)then
        snrm=0.
        do 11 i=1,n
          snrm=snrm+sx(i)**2
11      continue
        snrm=sqrt(snrm)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END

C     Subroutines for iterative solver of the fitting equations.
C     Convergence is poor. If it should be used, a good 
C     preconditioner must be implemented! /SR
      SUBROUTINE asolve(diag,n,b,x)
#include "implicit.h"
      DIMENSION diag(n),x(n),b(n)
      do i=1,n
c        x(i)=b(i)
        x(i)=b(i)/diag(i)/diag(i)
      enddo
      return
      END

      SUBROUTINE atimes(n,x,r,itrnsp,
     &                  NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
     &                  WORK,LWORK)
c     For direct product via FCK3DRV
c      DOUBLE PRECISION GABAL(*),WORK(*)
c     For (alpha|beta) in work
      DOUBLE PRECISION GABAL(*),WORK(n,n)
      INTEGER          NDMT,ISYMDM(*),IFCTYP(*),IPRINT,LWORK
cSimen
      INTEGER n,itrnsp,ija,NMAX
      DOUBLE PRECISION x(n),r(n),sa

      INTEGER i,j

c     For direct product via FCK3DRV

c      CALL DZERO(r,n)
c      CALL ITPROD(r,NDMT,ISYMDM,IFCTYP,
c     &       GABAL,x,IPRINT,WORK,LWORK)

c     For (alpha|beta) in work
      DO i=1,n
        r(i) = 0.0d0
        DO j=1,n
          r(i) = r(i) + WORK(i,j) * x(j)
        ENDDO
      ENDDO

      return
      END
C
      SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err,
     &                  NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
     &                  WORK,LWORK)
      DOUBLE PRECISION GABAL(*),WORK(*)
      INTEGER          NDMT,ISYMDM(*),IFCTYP(*),IPRINT,LWORK
      INTEGER iter,itmax,itol,n,NMAX
      DOUBLE PRECISION err,tol,b(*),x(*),EPS
      PARAMETER (NMAX=50000,EPS=1.d-14)
CU    USES atimes,asolve,snrm
      INTEGER j
      DOUBLE PRECISION ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,
     *znrm,p(NMAX),pp(NMAX),r(NMAX),rr(NMAX),z(NMAX),zz(NMAX),snrm
      iter=0
C
C aspg, 2006-06-23
C initializing bkden to 1
C
      bkden=1
C
      call atimes(n,x,r,0,
     &                  NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
     &                  WORK,LWORK)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
        write(2,*) 'Res ', j, r(j)
11    continue
C     call atimes(n,r,rr,0)
      znrm=1.d0
      if(itol.eq.1) then
        bnrm=snrm(n,b,itol)
        call quit('Error: itol=1 not allowed in this implementation!')
      else if (itol.eq.2) then
        call asolve(GABAL,n,b,z)
        bnrm=snrm(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
        call asolve(GABAL,n,b,z)
        bnrm=snrm(n,z,itol)
        call asolve(GABAL,n,r,z)
        znrm=snrm(n,z,itol)
      else
        call quit('illegal itol in linbcg')
      endif
      call asolve(GABAL,n,r,z)
100   if (iter.le.itmax) then
        iter=iter+1
        zm1nrm=znrm
        call asolve(GABAL,n,rr,zz)
        bknum=0.d0
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call atimes(n,p,z,0,
     &                  NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
     &                  WORK,LWORK)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call atimes(n,pp,zz,1,
     &                  NDMT,ISYMDM,IFCTYP,GABAL,IPRINT,
     &                  WORK,LWORK)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
        call asolve(GABAL,n,r,z)
        if(itol.eq.1.or.itol.eq.2)then
          znrm=1.d0
          err=snrm(n,r,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          znrm=snrm(n,z,itol)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*snrm(n,p,itol)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=snrm(n,x,itol)
          if(err.le.0.5d0*xnrm) then
            err=err/xnrm
          else
            err=znrm/bnrm
            goto 100
          endif
        endif
        write (*,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%12'%*Sim+).
C  /* Deck itprod */
      SUBROUTINE ITPROD(G_ALPHA,NDMT,ISYMDM,IFCTYP,
     &                  GABAL,C_ALPHA,IPRINT,
cSimen Added temporary to remove segmentation fault in setmaxcsbtch
     &                  GABAB,DMAT,
cSimen 
     &                  WORK,LWORK)
C
C   PURPOSE :   Driver routine for direct calculation of 
C
C                   g_alpha = sum_beta (alpha|beta) c_beta
C
C               where alpha is the fitting-basis and c are 
C               densities (multiple density vectors possible
C               through NDMT; generates multiple G_ALPHA's).
C
C               Multiple densities not debuged.
C
C               SR 19/12-2005.
C
#include "implicit.h"
#include "dummy.h"
      LOGICAL KCLERI
C
      DIMENSION G_ALPHA(*), ISYMDM(*), IFCTYP(*), 
cSimen Added temporary to remove segmentation fault in setmaxcsbtch
     &          GABAB(*),DMAT(*),
cSimen 
     &          WORK(*), GABAL(*), C_ALPHA(*)
C
C
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbihr2.h"
#include "cbieri.h"
#include "cbifck3.h"
#include "priunit.h"
#include "densfit.h"
#include "fck3thr.h"
C
      LOGICAL DOFMM_SAVE, NOSCRE_SAVE
#include "memint.h"
C
      CALL QENTER('ITPROD')
c     Turn off fmm
      DOFMM_SAVE = DOFMM
      DOFMM = .FALSE.
   
      IF (DOFMM) THEN
        KCLERI = NCLERI
        NCLERI = .TRUE.
        NCLHER = .TRUE.
      END IF

      DO IDMT=1,NDMT
        ISTR = N2BASX*(IDMT - 1) + 1
        ISTRAUX = NBASISAUX*(IDMT - 1) + 1
C
C aspg, 2006-08-14:
C
C due to the differences in the way dirac and dalton's dft code 
C work at the time of writing (in dirac's setfck_ ifctyp is 
C changed to "coulomb"(=1) whereas in dalton apparently ifcktyp 
C not changed and remains as "coulomb+exchange"(=3) ), the
C if statement below for dirac tests for both the pure coulomb
C and coulomb+exchange cases 
C
#if defined (PRG_DIRAC)
        IF ((MOD(IFCTYP(IDMT),10) .EQ. 1).OR. 
     &      (MOD(IFCTYP(IDMT),10) .EQ. 3)) THEN 
#else
        IF ( MOD(IFCTYP(IDMT),10) .EQ. 3 ) THEN
#endif
Caspg
          RHSFIT = .TRUE.
          LHSFIT = .TRUE.
          NOSCRE_SAVE = NOSCRE
          NOSCRE = .FALSE.  
          THRFCK3_SAVE = THRFCK3
          THRFCK3 = 1.0d-25
cSimen Added temporary to remove segmentation fault in setmaxcsbtch
          CALL FCK3DRV(G_ALPHA(ISTRAUX),DMAT,1,GABAB,GABAL,
c         CALL FCK3DRV(G_ALPHA(ISTRAUX),dummy,1,dummy,GABAL,
cSimen
     &                 C_ALPHA(ISTRAUX),WORK,LWORK)
          THRFCK3 = THRFCK3_SAVE
          NOSCRE = NOSCRE_SAVE
#ifdef MOD_FMM
          IF (DOFMM) THEN
cToDo       implemend LHS-fitting moments.
            CALL FMMFCK(G_ALPHA(ISTRAUX),C_ALPHA(ISTRAUX),WORK,LWORK)
          END IF
#else
          IF (DOFMM) CALL QUIT('FMM requested but not enabled')
#endif
        END IF
      END DO

      IF (DOFMM) NCLERI = KCLERI
      DOFMM = DOFMM_SAVE
      
      CALL QEXIT('ITPROD')
      RETURN
      END


      SUBROUTINE SCRINT(GABAB,GABAL,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "iteration.h"
      DIMENSION GABAB(*), GABAL(*), WORK(LWORK)
#include "memint.h" 
      CALL QENTER('SCRINT')

C     Set AO- and OD-batches (and classes)
      RHSFIT = .FALSE.
      LHSFIT = .FALSE.
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,
     &     IPRINT)
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &     KFREE,LFREE,WORK,LFREE,DUMMY,DUMMY,DUMMY,DUMMY,
     &     1,IPRINT,.false.)
      CALL ACT_ALL_BATCHES

C     Generate the 2-index Cauchy-Schwarz integrals
      CALL GETINT3(GABAB,WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.FALSE.,.FALSE.,
     &     IDUMMY,IDUMMY,.TRUE.)

      IF (IPRINT.GT.10) THEN
        CALL HEADER('GAB in SCRINT',-1)
        CALL OUTPUT(GABAB,1,NBASIS,1,NBASIS,NBASIS,NBASIS,1,LUPRI)
      ENDIF

C     Normalize the fit set.
C     1. Generate the 1-index Cauchy-Schwarz integrals over the fit set 
      CALL GETINT3(GABAL,WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     IDUMMY,IDUMMY,.TRUE.)
C     2. Normalize
      CALL NRMFIT(GABAL,WORK(KFREE),LFREE,IPRINT)
C     3. Set up batches and recalculate 1-index Cauchy-Schwarz integrals
C        over the fit set (should produce ones)
C        FIXME: these are now obsolete for screening and should be removed.
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,
     &     IPRINT)
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &     KFREE,LFREE,WORK,LFREE,DUMMY,DUMMY,DUMMY,DUMMY,
     &     1,IPRINT,.false.)
      CALL GETINT3(GABAL,WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     IDUMMY,IDUMMY,.TRUE.)
      IF (IPRINT.GT.10) THEN
        CALL HEADER('GABAL in SCRINT',-1)
        CALL OUTPUT(GABAL,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
      CALL QEXIT('SCRINT')
      RETURN
      END


      SUBROUTINE GTALBE(ALBE,I_LG,NDIMAUX,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "iteration.h"
      DIMENSION ALBE(*), WORK(LWORK), I_LG(*)
#include "memint.h" 
      CALL QENTER('GTALBE')

C     Set AO- and OD-batches (and classes)
      RHSFIT = .FALSE.
      LHSFIT = .FALSE.
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,
     &     IPRINT)
      CALL SETODBTCH(KIODBCLS,KRODBCLS,KIODBTCH,KRODBTCH,
     &     KFREE,LFREE,WORK,LFREE,DUMMY,DUMMY,DUMMY,DUMMY,
     &     1,IPRINT,.false.)
c     CALL ACT_ALL_BATCHES

C     Calculate the (al|be) matrix
      CALL GETINT3(ALBE,WORK(KIODBCLS),WORK(KIODBTCH),
     &     WORK(KRODBTCH),
     &     WORK(KFREE),LFREE,IPRINT,.TRUE.,.TRUE.,
     &     I_LG,NDIMAUX,.FALSE.)
      CALL QEXIT('GTALBE')
      RETURN
      END


      SUBROUTINE DF_BOX(CFIT,IBOX,NUMBOX,CORDBOX,FMAT,DMAT,NDMT,
     &                  GABAB,GABAL,ISYMDM,IFCTYP,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "odclss.h"
#include "inforb.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbifck3.h"
#include "densfit.h"
#include "fck3aob.h"
#include "iteration.h"
      DIMENSION CFIT(NBASISAUX,NDMT), NUMBOX(*), CORDBOX(*), 
     &          FMAT(*), DMAT(*), 
     &          GABAB(*), GABAL(*), ISYMDM(*), IFCTYP(*), WORK(LWORK)
#include "memint.h" 
C
      CALL QENTER('DF_BOX')
C
      NORBBB  = NUMORBBOXBUF(IBOX,NUMBOX)
      NORBBOX = NUMORBBOX(IBOX,NUMBOX)
C
      CALL MEMGET('REAL',KWVEC,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET('INTEGER',KINDBOX,NAOBTCH,WORK,KFREE,LFREE)
      CALL MEMGET('INTEGER',KI_LG,NORBBB,WORK,KFREE,LFREE)
C
      CALL SETAOBTCH(.FALSE.,CONTRACTED,WORK(KFREE),LFREE,IPRINT)
      CALL FIND_ACTAO_BOXBUF(WORK(KINDBOX),IBOX,NUMBOX,
     &                       CORDBOX,IPRINT)
      CALL PICK_ACT_IND(WORK(KINDBOX),IPRINT)
      CALL INDEX_LOC_GLOB(WORK(KI_LG),NORBBB,WORK(KINDBOX),IPRINT)
C
      IF (DOBOX) THEN
        CALL MEMGET('REAL',KC_LOC,NORBBB*NDMT,WORK,KFREE,LFREE)
        CALL DF_GENDRV (WORK(KC_LOC),DMAT,NORBBB,NDMT,WORK(KI_LG),
     &                  ISYMDM,IFCTYP,
     &                  GABAB,GABAL,IBOX,
     &                  WORK(KFREE),LFREE,IPRINT)
        DO IDMT=1,NDMT
          CALL ADD_PROJECT_VEC(CFIT(1,IDMT),NBASISAUX,
     &                         WORK(KC_LOC+(IDMT-1)*NORBBB),NORBBB,
     &                         WORK(KI_LG),.FALSE.,IPRINT)
        ENDDO
        CALL MEMREL('DF_BOX',WORK,KC_LOC,KC_LOC,KFREE,LFREE)
      ELSE
        CALL DF_GENDRV (CFIT,DMAT,NORBBB,NDMT,WORK(KI_LG),ISYMDM,IFCTYP,
     &                  GABAB,GABAL,IBOX,
     &                  WORK(KFREE),LFREE,IPRINT)
      ENDIF
C
      CALL MEMREL('DF_BOX',WORK,KI_LG,KI_LG,KFREE,LFREE)
      CALL MEMREL('DF_BOX',WORK,KINDBOX,KINDBOX,KFREE,LFREE)
      CALL MEMREL('DF_BOX',WORK,KWVEC,KWVEC,KFREE,LFREE)
C
      CALL QEXIT('DF_BOX')
C
      RETURN
      END

C     LINSOL solves the linear equations Ax=b, where x and b can have 
C     several columns (NRHS), and A is symmetric positive definite matrix.
C     MAKECOPY indicates wether we would like to make a copy of A (so that
C     it remains unchanged by this call) or not.
      SUBROUTINE LINSOL(AMAT,XVEC,BVEC,N,NRHS,MAKECOPY,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMAT(N,N), XVEC(N,NRHS), BVEC(N,NRHS), WORK(LWORK)
      LOGICAL MAKECOPY
      CALL QENTER('LINSOL')
      CALL DCOPY(N*NRHS,BVEC,1,XVEC,1)
      IF (.NOT. MAKECOPY) THEN
C       Overwrites the matrix AMAT
        CALL DPOSV ('L',N,NRHS,AMAT,N,XVEC,N,INFO)
      ELSE IF (N*N.LE.LWORK) THEN
C       Makes a copy of the AMAT, so that AMAT is unchanged
        CALL DCOPY(N*N,AMAT,1,WORK(1),1)
        CALL DPOSV ('L',N,NRHS,WORK(1),N,XVEC,N,INFO)
      ELSE
        CALL STOPIT('LINSOL','DPOSV ',N*N,LWORK)
      ENDIF
      IF (INFO.NE.0) THEN
        WRITE(LUPRI,*) 'No convergence in DPOSV from DF_DPOSV, INFO=',
     &                 INFO
        CALL QUIT('No convergence in DPOSV (from DF_DPOSV)')
      ENDIF
      CALL QEXIT('LINSOL')
      RETURN
      END


      SUBROUTINE DF_GET_AUXCHA_LOC(AUXINT,I_LG,NDIMAUX,
     &                             WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "densfit.h"
      DIMENSION AUXINT(NDIMAUX),I_LG(NDIMAUX),WORK(*)
#include "memint.h"
      CALL QENTER('GET_AUXC')
      IF (DOBOX) THEN
        CALL MEMGET('REAL',KAUXGLO,NBASISAUX,WORK,KFREE,LFREE)
        CALL AUXCHA(WORK(KAUXGLO),WORK(KFREE),LFREE,IPRINT)
        CALL PROJECT_VEC(WORK(KAUXGLO),NBASISAUX,AUXINT,
     &                   NDIMAUX,I_LG,1,.TRUE.,IPRINT)
        CALL MEMREL('GET_AUXC',WORK,KAUXGLO,KAUXGLO,KFREE,LFREE)
      ELSE
        IF (NDIMAUX.NE.NBASISAUX) 
     &     CALL QUIT('NDIMAUX.NE.NBASISAUX in DF_GET_AUXCHA_LOC')
        CALL AUXCHA(AUXINT,WORK(KFREE),LFREE,IPRINT)
      ENDIF
      CALL QEXIT('GET_AUXC')
      RETURN
      END

      LOGICAL FUNCTION IS_DIAGZERO(AMAT,N)
#include "implicit.h"
      DIMENSION AMAT(N,N)
      AMATMAX = 0.0d0
      DO I=1,N
        AMATMAX = MAX(AMATMAX,ABS(AMAT(I,I)))
      ENDDO
      IF (AMATMAX.GE.1.0d-14) THEN
        IS_DIAGZERO = .FALSE.
      ELSE 
        IS_DIAGZERO = .TRUE.
      ENDIF
      RETURN
      END

      FUNCTION DF_E_ROB_CORR()
#include "implicit.h"
#include "densfit.h"
#include "priunit.h"
      IF (DOFIT.AND.ROBUST_E) THEN
        DF_E_ROB_CORR = ROB_E_CORR
      ELSE
        DF_E_ROB_CORR = 0.0d0
      ENDIF
      WRITE(LUPRI,'(A,F16.10)') 'Robust density-fitting contribution '//
     &               'added to the total energy', DF_E_ROB_CORR
      RETURN
      END

      SUBROUTINE DF_CAL_ROB_E(DMAT,FMAT,CFIT,NDMT,ISYMDM,IFCTYP,GABAL,
     &                        WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "densfit.h"
      DIMENSION DMAT(*),FMAT(*),CFIT(*),ISYMDM(*),IFCTYP(*),GABAL(*),
     &          WORK(*)
#include "memint.h"
      CALL QENTER('CAL_ROB_')
      IF (NDMT.NE.1) 
     &  CALL QUIT('Robust correction not implenented for NDMT .GT. 1')
      ROB_E_CORR = TRACE_MAT_PRO(DMAT,FMAT,NBASIS,NBASIS)
      CALL MEMGET('REAL',KJALPHA,NBASISAUX*NDMT,WORK,KFREE,LFREE)
c     CALL ITPROD(WORK(KJALPHA),NDMT,ISYMDM,IFCTYP,GABAL,CFIT,IPRINT,
CAREFUL Added temporary to remove segmentation fault in setmaxcsbtch
#if 1
      CALL ITPROD(WORK(KJALPHA),NDMT,ISYMDM,IFCTYP,GABAL,CFIT,DMAT,DMAT,
     &            IPRINT,WORK(KFREE),LFREE)
#else
      CALL ITPROD(WORK(KJALPHA),NDMT,ISYMDM,IFCTYP,GABAL,CFIT,
     &            IPRINT,WORK(KFREE),LFREE)
#endif
      ROB_E_CORR = ROB_E_CORR - DDOT(NBASISAUX,CFIT,1,WORK(KJALPHA),1)
      ROB_E_CORR = ROB_E_CORR * 0.5d0
      CALL MEMREL('CAL_ROB_',WORK,KJALPHA,KJALPHA,KFREE,LFREE)
      CALL QEXIT('CAL_ROB_')
      RETURN
      END


C /* deck DF_CAL_ApproxERRKS */
C aspg, 2007-06-02:
C calculating approximated ERRKS necessary for DFT (for df-vxc), 
C as the difference between the enerfy for the (df+dfvxc) fock
C matrix and the ddot(cfit,j_alpha) matrix. my understanding of
C a suggestion from simen... 
      SUBROUTINE DF_CAL_ApproxERRKS(DMAT,FMAT,CFIT,NDMT,ISYMDM,
     &                              IFCTYP,GABAL,WORK,LWORK,IPRINT)

#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "densfit.h"
      DIMENSION DMAT(*),FMAT(*),CFIT(*),ISYMDM(*),IFCTYP(*),GABAL(*),
     &          WORK(*)
#include "memint.h"

      CALL QENTER('AppERRKS')

      IF (NDMT.NE.1)
     &  CALL QUIT('Robust correction not implenented for NDMT .GT. 1')

      dft_temp_xc_mat_energy = TRACE_MAT_PRO(DMAT,FMAT,NBASIS,NBASIS)

      CALL MEMGET('REAL',KJALPHA,NBASISAUX*NDMT,WORK,KFREE,LFREE)

      CALL ITPROD(WORK(KJALPHA),NDMT,ISYMDM,IFCTYP,GABAL,CFIT,DMAT,DMAT,
     &            IPRINT,WORK(KFREE),LFREE)
      dft_temp_xc_mat_energy = dft_temp_xc_mat_energy 
     &                       - DDOT(NBASISAUX,CFIT,1,WORK(KJALPHA),1)
      dft_temp_xc_mat_energy = dft_temp_xc_mat_energy * 0.5d0

      CALL MEMREL('AppERRKS',WORK,KJALPHA,KJALPHA,KFREE,LFREE)
      CALL QEXIT('AppERRKS')
      END



C  /* Deck dfiinp */
      SUBROUTINE DFIINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "densfit.h"
#include "fck3box.h"
      PARAMETER ( NTABLE = 13 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DATA TABLE /'.DIATOM','.DOBOX ','.BOXORB','.BUFSIZ','.CONSTR',
     &            '.NOROBU','.ATTENU','.GAUDAM','.OLD_DF','.SCREEN',
     &            '.SAVE_C','.SAVE_V','.DF_VXC'/
C
      WORD1 = WORD
C       
 100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
        
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13), I
                  END IF
 200          CONTINUE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FCK3INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FCK3INP.')
 1    CONTINUE
        !.DIATOM - fit density for pairs of atoms
        DIATOM     = .TRUE.
      GO TO 100
 2    CONTINUE
        DOBOX  = .TRUE.
        LIN_DF = .TRUE.
      GO TO 100
 3    CONTINUE
        READ (LUCMD,*) MAXORBBOX
        IF ((MAXORBBOX .LE. 1) .OR. (MAXORBBOX .GT. 100000)) THEN
C
C aspg, 2006-06-27 making ftnchek happy
C
          WRITE(LUPRI,*) 'Maxorbox input error, maxorbox=', MAXORBBOX
          CALL QUIT('Maxorbox input error')
        ENDIF
        DOBOX = .TRUE.
        LIN_DF = .TRUE.
      GO TO 100
 4    CONTINUE
        READ (LUCMD,*) RBUFFER
        IF (RBUFFER.LT.0.0d0) THEN
          WRITE(LUPRI,*) 'Rbuffer input error, Rbuffer=', RBUFFER
          CALL QUIT('Rbuffer input error')
        ENDIF
        DOBOX = .TRUE.
        LIN_DF = .TRUE.
      GO TO 100
 5    CONTINUE
        DF_CONSTR = .TRUE.
      GO TO 100
 6    CONTINUE
        ROBUST_E = .FALSE.
      GO TO 100
 7    CONTINUE
        DF_ATTENU = .TRUE.
        READ (LUCMD,*) DF_ATTPAR
        IF (DF_ATTPAR.LT.0.0d0) THEN
          WRITE(LUPRI,*) 'DF_ATTPAR input error, DF_ATTPAR=', DF_ATTPAR
          CALL QUIT('DF_ATTPAR input error')
        ENDIF
        IF (DF_GAUDMP) THEN
          WRITE(LUPRI,*)
     &     'Cannot have both .ATTENU and .GAUDAM simultaneously'
          CALL QUIT('Input error: .ATTENU and .GAUDAM simult.')
        ENDIF
      GO TO 100
 8    CONTINUE
        DF_GAUDMP = .TRUE.
        READ (LUCMD,*) DF_GAUPAR
        IF (DF_GAUPAR.LT.0.0d0) THEN
          WRITE(LUPRI,*) 'DF_GAUPAR input error, DF_GAUPAR=', DF_GAUPAR
          CALL QUIT('DF_GAUPAR input error')
        ENDIF
        IF (DF_ATTENU) THEN
          WRITE(LUPRI,*)
     &     'Cannot have both .ATTENU and .GAUDAM simultaneously'
          CALL QUIT('Input error: .ATTENU and .GAUDAM simult.')
        ENDIF
      GO TO 100
 9    CONTINUE
        LIN_DF = .FALSE.
      GO TO 100
 10   CONTINUE
        READ (LUCMD,*) DF_SCRTHR
        WRITE(LUPRI,*) 'Setting screening threshold in DF to',DF_SCRTHR
      GO TO 100
 11   CONTINUE
        WRITE(LUPRI,*) 'Saving fit coefficients at each iteration.'
        SAVE_CFIT = .TRUE.
      GO TO 100
 12   CONTINUE
        WRITE(LUPRI,*) 'Saving (alpha|beta) matrix at each iteration.'
        SAVE_VAB  = .TRUE.
      GO TO 100
 13   CONTINUE
        WRITE(LUPRI,'(A)') 'Density Fitting will be used to '//
     &   'build the Exchange-Correlation contribution',
     &   '   according to the scheme of Birkenheuer and coworkers'//
     &   '(IJQC 102, 743, 2005)'
        DF_FITVXC = .TRUE.
        SAVE_CFIT = .TRUE.
        SAVE_VAB  = .TRUE.
      GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DFIINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 CALL QUIT('Illegal keyword in DFIINP.')
            END IF
 300  CONTINUE

*     CALL HEADER('Settings for FCK3 calculation:',0)

*     IF (IPRINT .NE. IPRFDEF) THEN
*         WRITE (LUPRI,'(A,I5)') ' Print level in FCK3:',IPRFCK3
*     END IF

*     IF ( FIXTHR ) THEN
*       WRITE (LUPRI,'(A,I5)') 'Threshold is fixed.'
*     ELSE
*       WRITE (LUPRI,'(A,I5)') 'Threshold set dynamically.'
*     END IF

*     IF ( NOSKIPCL ) THEN
*       WRITE (LUPRI,'(A,I5)') 
*    &  'No skipping of the nonclassical integrals (only set to 0).'
*     END IF

      RETURN
      END
C  /* Deck dfiini */
      SUBROUTINE DFIINI
#include "implicit.h"
#include "densfit.h"
Caspg, 2006-10-15
C still thinking about the best way to communicate whether or
C not we will be using differential densities to de df routines.. 
C if not done correctly, fitvxc will not work properly
C one way would be to get the value of dodscf from
C
C#include "dcbdhf.h"
C
C but this is initalized after the integrals...
C
      !.DIATOM - fit density for pairs of atoms (default is the global fit)
      DOFIT     = .FALSE.
      DIATOM    = .FALSE.
      LHSFIT    = .FALSE.
      RHSFIT    = .FALSE.
      ROBUST_E  = .TRUE.
      DOBOX     = .FALSE.
      DF_CONSTR = .FALSE.
      DF_GAUDMP = .FALSE.
      DF_ATTENU = .FALSE.
      DF_FITVXC = .FALSE.
      LIN_DF    = .TRUE.
      DF_SCRTHR = 1.0E-10
      SAVE_CFIT = .FALSE.
      SAVE_VAB  = .FALSE.
      DODFLR    = .FALSE.
C
C     using_diffdens = .false.
      using_diffdens = .true.
C     if (dodscf.eq..false.) using_diffdens = .false.
C
C     LIN_DF    = .FALSE.
      CALL SET_BOX_DEFAULTS()
      RETURN
      END


