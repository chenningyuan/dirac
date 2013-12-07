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

C  /* Deck dens_fit_old */
      SUBROUTINE DENS_FIT_OLD (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                     WORK,LWORK)
C*****************************************************************************
C    
C     DENS_FIT_OLD: Entree routine for density fitting in Dirac & Dalton.
C        
C     Written by Luuk Visscher, october 2004
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
#include "symmet.h"
#include "cbieri.h"
#include "odclss.h"
#include "inforb.h"
#include "inftap.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
      DIMENSION KOD(12)
CC    DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)

      LOGICAL FIRST_CALL
      SAVE FIRST_CALL
      DATA FIRST_CALL /.TRUE./
#include "memint.h" 
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C     Transfer information about the number of density matrices
C     to the eri common block
C
      NDMAT = NDMT
C
C     Allocate memory for fit coefficients, work arrays,
C     contraction coefficients and index array
C     TODO : check whether allocation of FTMP is really necessary
C
      CALL MEMGET ('REAL',KCFIT,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KWMAT,NBASISAUX*NBASISAUX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KWVEC,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KFTMP,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KCCFBT,MXPRIM*MXCONT,WORK,KFREE,LFREE)
      CALL MEMGET ('INTE',KINDXB,MXSHEL*MXCONT*8,WORK,KFREE,LFREE)
C
C     Set the CCFBT and INDBX arrays and get the ODC pointers
C
      CALL ERI_POINT (WORK(KFREE),LFREE,WORK(KCCFBT),WORK(KINDXB),
     &                KOD,IPRFCK)
C
C     Update KFREE and LFREE since we want to keep the data generated in
C     ERIPOINT. We will call that part of WORK WKOD from now on.
C
      KWKOD = KFREE
      KFREE = KFREE + KOD(11) - 1
      LFREE = KOD(12)
C
C     Generate the 2-index repulsion integrals over the fit set
C
      CALL TIMER('MEMORY',TIMSTR,TIMEND)
      IF (FIRST_CALL) THEN
C        The 2-index integrals are written to disk
         WRTINT = .TRUE.
         FCKINT = .FALSE.
C
C        We have symmetry between density 1 and 2
C
         ODTR12 = .TRUE.
C
C        Set the active classes for the ERI calculation.
C        For each electron we take as first index the fit set functions and 
C        as second the null function.
C
C        Index 1
         ICLASSES(1) = 3
         NCENTS(1)   = 0
         ICENTS(1,1) = 0
C        Index 2
         ICLASSES(2) = 0
         NCENTS(2)   = 1
         ICENTS(1,2) = 0
C        Index 3
         ICLASSES(3) = 3
         NCENTS(3)   = 0
         ICENTS(1,3) = 0
C        Index 4
         ICLASSES(4) = 0
         NCENTS(4)   = 1
         ICENTS(1,4) = 0
C
         CALL PICK_ACT_BATCH_OLD (NCENTS,ICENTS,ICLASSES,IPRFCK)
C
C        Open the file to which these integrals are to be written
C
         LUINTA = -1
         CALL GPOPEN(LUINTA,'FIT_INTS','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
C        Generate and write the integrals
C
         CALL ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,0,IPRFCK,
     &                    WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                    WORK(KWKOD),WORK(KFREE),LFREE)
C
         FIRST_CALL = .FALSE.
      ENDIF
      CALL TIMER('F_CALL',TIMSTR,TIMEND)
C
C     From now on we will produce (modified) Fock matrices in ERI
C
      WRTINT = .FALSE.
      FCKINT = .TRUE.

C
C     We will have no symmetry between density 1 and 2
C
      ODTR12 = .FALSE.
C
C     Get the fit coefficients by looping over atom pairs
C
      CALL DRIVE_FIT_OLD (FMAT,DMAT,WORK(KWMAT),WORK(KWVEC),NDMT,
     &                ISYMDM,IFCTYP,
     &                WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                WORK(KWKOD),WORK(KFREE),LFREE,WORK(KCFIT),
     &                IPRFCK)
      CALL TIMER('FITING',TIMSTR,TIMEND)

cfp
c     print at most first 5000 elements of fitted density to see if denfit works
      write(lupri,*)'C_al (fitted density) vector check for DOFCK3 ',
     & DOFCK3
      write(lupri,*)'NBASISAUX ',NBASISAUX
      NRANGEEEE = NBASISAUX
      if (NBASISAUX.gt.5000) then
         NRANGEEEE = 5000
         write(lupri,*)'only first 5000 elemnts will be printed'
      end if
      do ii = 1,NRANGEEEE
         write(lupri,*) ii, WORK(KCFIT+ii-1)
      end do
cfp


      IF (DOFCK3) THEN
        CALL FITFCK3_OLD(FMAT,WORK(KCFIT),NDMT,ISYMDM,IFCTYP,IPRFCK,
     &               WORK(KFREE),LFREE)
      ELSE

C
C       Build the fock matrix with the fitted density
C       Now the first two indices are regular and the last two relate
C       to the fit set and the null function respectively.
C
C       Index 1
        ICLASSES(1) = 1
        NCENTS(1)   = 0
        ICENTS(1,1) = 0
C       Index 2
        ICLASSES(2) = 1
        NCENTS(2)   = 0
        ICENTS(1,2) = 0
C       Index 3
        ICLASSES(3) = 3
        NCENTS(3)   = 0
        ICENTS(1,3) = 0
C       Index 4
        ICLASSES(4) = 0
        NCENTS(4)   = 1
        ICENTS(1,4) = 0
C     
        CALL PICK_ACT_BATCH_OLD (NCENTS,ICENTS,ICLASSES,IPRFCK)
C     
C       Compute the integrals and build the Fock matrix
C     
        IF (IPRFCK.GT.2) THEN
          CALL HEADER('Fit coefficients ',-1)
          CALL OUTPUT(WORK(KCFIT),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
        ENDIF
C     
        CALL ERI_MKINTS (FMAT,WORK(KCFIT),NDMT,ISYMDM,IFCTYP,1,IPRFCK,
     &                   WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                   WORK(KWKOD),WORK(KFREE),LFREE)

      ENDIF

      CALL TIMER('FabBUI',TIMSTR,TIMEND)
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Fitted Fock matrix',-1)
         CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C---->The code below restores the original Fock matrix. Should be removed
C     once the implementation is tested and full debugged !
C
*     odtr12 = .true.
*     call pickao(0)
*     call dzero(fmat,ndmt*n2basx)
*     CALL ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,0,IPRFCK,
*    &                 WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
*    &                 WORK(KWKOD),WORK(KFREE),LFREE)
*     IF (IPRFCK.GT.2) THEN
*        CALL HEADER('Exact Fock matrix',-1)
*        CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
*     ENDIF
C----< End of test code.
C
C     Release the memory allocated for the Fock matrix build. This
C     will also check for errors due to memory bound failures
C
      CALL MEMREL('DENS_FIT_OLD',WORK,KCFIT,KCFIT,KFREE,LFREE)
C
      RETURN
      END
C


C  /* Deck fitfck3 */
      SUBROUTINE FITFCK3_OLD(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
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
      LOGICAL KCLERI
C
      DIMENSION FMAT(*), DMAT(*), ISYMDM(*), IFCTYP(*), 
     &          WORK(*)
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
#include "dummy.h"
#include "memint.h"
C
      CALL QENTER('FITFCK3_OLD')

c      USEERI = RUNERI
c      DO I = 1, NDMT
c        USEERI = USEERI .AND. ((IFCTYP(I).EQ.13).OR.(IFCTYP(I).EQ.3))
c      END DO
c      
#ifdef MOD_FMM
      IF (DOFMM) THEN
        KCLERI = NCLERI
        NCLERI = .TRUE.
        NCLHER = .TRUE.
      END IF
#endif
      
      DO IDMT=1,NDMT
        ISTR = N2BASX*(IDMT - 1) + 1
        ISTRAUX = NBASISAUX*(IDMT - 1) + 1
cSimen There might be a problem here!
        IF (NBASISAUX.GT.N2BASX) CALL QUIT('Error: NBASISAUX > N2BASX.')
cSimen
        IF ( MOD(IFCTYP(IDMT),10) .EQ. 3 ) THEN
          RHSFIT = .TRUE.
          LHSFIT = .FALSE.
          CALL ACT_ALL_BATCHES
          CALL FCK3DRV(FMAT(ISTR),DUMMY,1,DUMMY,DUMMY,DMAT(ISTRAUX),
     &                 WORK,LWORK)
#ifdef MOD_FMM
          IF (DOFMM) THEN
cSimen      CALL FMMFCK(FMAT(ISTR),DMAT(ISTR),WORK,LWORK)
            WRITE(LUPRI,*) 'FMM not yet implemented for densfit FCK3'
            CALL QUIT('FMM not enabled')
          END IF
#else
Caspg     IF (DOFMM) CALL QUIT('FMM requested but not enabled')
#endif
        END IF
      END DO

#ifdef MOD_FMM
      IF (DOFMM) NCLERI = KCLERI
#endif

      CALL QEXIT('FITFCK3_OLD')
      RETURN
      END
C  /* Deck eri_point */
      SUBROUTINE ERI_POINT (WORK,LWORK,CCFBT,INDXBT,KOD,IPRFCK)
C*****************************************************************************
C    
C     ERI_POINT: Initialization of ERI pointers for Fock matrix calculation.
C        
C     Taken from HRFCK2 by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "ccom.h"
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "erimem.h"
#include "aobtch.h"
#include "odbtch.h"
#include "shells.h"
#include "symmet.h"
#include "infpar.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
      DIMENSION  WORK(*), CCFBT(*), INDXBT(*)
      DIMENSION KOD(12)
#include "memint.h" 
C
C     Initialization of ERI (if we are the master)
C
      IF (SLAVE) THEN
         IPRINT = IPRFCK
      ELSE
C        
C        Initialization in ER2INI
C
         CALL ER2INI
C
         IPRINT = MAX(IPRERI,IPRFCK)
      END IF
C
C     Give permutational symmetries
C
      PMS12 = .TRUE.
      PMSAB = .FALSE.
      PMSCD = .FALSE.
C
      THRSH  = MAX(THRS,1.00D-15)
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
      EXCINT = .FALSE.
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C     
      CALL ODBCHS(KOD(1),KOD(2),KOD(3),KOD(4),KOD(5),
     &            KOD(6),KOD(7),KOD(8),KOD(9),KOD(10),
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
C
C     Save the current values of KFREE and LFREE, these are called 
C     KLAST and LWRK, resp. in HR2FCK
C
      KOD(11) = KFREE
      KOD(12) = LFREE
C     
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C     
      ICALL = 0
      CALL GETDST(ICALL,ICALL,IPRINT)
C
      RETURN
      END
C
C  /* Deck eri_mkints */
      SUBROUTINE ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IFIT_DMAT,
     &                       IPRFCK,CCFBT,INDXBT,KOD,FCKTMP,WKOD,
     &                       WORK,LWORK)
C*****************************************************************************
C    
C     ERI_MKINTS: Make integrals by calling ERI. 
C     The active classes should be set outside of this routine.
C     Control of writing versus Fock matrix construction is also done
C     on the outside.
C        
C     Second part of HRFCK2, split by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "ccom.h"
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "erimem.h"
#include "aobtch.h"
#include "odbtch.h"
#include "symmet.h"
#include "infpar.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*)
      DIMENSION KOD(12)
C
C     Transfer information about the type of Fock matrix to ERI common
C     Set print level for ERI
C
      IFITDM = IFIT_DMAT
      IPRINT = MAX(IPRERI,IPRFCK)
C     
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)
C
      KODCL1 = KOD(1)
      KODCL2 = KOD(2)
      KODBC1 = KOD(3)
      KODBC2 = KOD(4)
      KRDBC1 = KOD(5)
      KRDBC2 = KOD(6)
      KODPP1 = KOD(7)
      KODPP2 = KOD(8)
      KRDPP1 = KOD(9)
      KRDPP2 = KOD(10)
C
C     Calculate integrals
C     ===================
C
      IF (SLAVE) THEN
         CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &               WKOD(KODBC1),WKOD(KODBC2),
     &               WKOD(KRDBC1),WKOD(KRDBC2),
     &               WKOD(KODPP1),WKOD(KODPP2),
     &               WKOD(KRDPP1),WKOD(KRDPP2),
     &               DUMMY,DUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &               FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK,LWORK,IPRINT)
      ELSE
         IF (.NOT.INTSKP) THEN
            CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &                  WKOD(KODBC1),WKOD(KODBC2),
     &                  WKOD(KRDBC1),WKOD(KRDBC2),
     &                  WKOD(KODPP1),WKOD(KODPP2),
     &                  WKOD(KRDPP1),WKOD(KRDPP2),
     &                  DUMMY,DUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &                  FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK,LWORK,IPRINT)
C
C           Error message in case of insufficient memory
C
            IF (.NOT.MEMOK) THEN
               WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &            ' Not enough memory for this run of ERIFCK.',
     &            ' Available memory in ERIFCK:',LWORK,
     &            ' Required memory for ERIFCK:',LWORK + MEMADD,
     &            ' Increase memory (LWORK) by:',MEMADD
               WRITE (LUPRI,'(/,1X,A,2I5)')
     &            ' Memory requirements largest for OD classes :',
     &              MODAB,MODCD
               CALL QUIT('Insufficient memory in ERIFCK.')
            END IF
         END IF
C
C        Copy and/or print of Fock matrix if we are in that branch
C        This code will only work with regular 2-index Fock matrices,
C        not with the modified 1-index Fock matrix that is made if
C        IFITDM .EQ. 2. The IF statement below should catch this.
C
         IF (FCKINT.AND.ABS(IFITDM).LE.1) THEN
C
C           Print densities and Fock matrix
C           ===============================
C
            IF (IPRINT.GT.4) THEN
               CALL HEADER('Density and Fock matrices in ERIFCK',-1)
               KSTR = 1
               DO I = 1, NDMT
                  WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
                  CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                        NBASE,1,LUPRI)
                  WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
                  CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                        NBASE,1,LUPRI)
                  KSTR = KSTR + NBASE*NBASE
               END DO
            END IF
C
C        End of IF Block for FMAT construction
C
         END IF
C
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END
C
C  /* Deck drive_fit */
       SUBROUTINE DRIVE_FIT_OLD (FMAT,DMAT,WMAT,WVEC,NDMT,ISYMDM,
     &                       IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                       WKOD,WORK,LWORK,CFIT,IPRINT)
C*****************************************************************************
C    
C     DRIVE_FIT: Driver for density fitting.
C        
C     Written by Luuk Visscher, november 2003
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
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*), CFIT(*)
      DIMENSION KOD(12)
#include "nuclei.h"
#include "symmet.h"
C
#include "memint.h" 
C
C     Read the 2-index integrals from file
C
      CALL READ_FIT2INT(WMAT,WORK(KFREE),LFREE)
cfp
c     print at most first 5000 elements of (al|be) matrix to see if denfit works
      write(lupri,*)'(al|be) matrix check for DOFCK3 ',
     & DOFCK3
      write(lupri,*)'NBASISAUX*NBASISAUX ', NBASISAUX*NBASISAUX
      NRANGEEEE = NBASISAUX*NBASISAUX
      if (NBASISAUX*NBASISAUX.gt.5000) then
         NRANGEEEE = 5000
         write(lupri,*)'only first 5000 elements will be printed'
      end if
      do ii = 1,NRANGEEEE
         write(lupri,*) ii, WMAT(ii)
      end do
cfp

C
C     Allocate memory for the indexing of fit functions
C
      CALL MEMGET ('INTE',KINDAB,NBASISAUX,WORK,KFREE,LFREE)
C
C     Initialize the fit coefficients
C
      CALL DZERO (CFIT,NBASISAUX*NDMT)
C
C     **************************************************
C     ************* Carry out the fitting ************** 
C     **************************************************
C
      !If requested, loop over symmetry unique center pairs (diatom fit)
      IF (DIATOM) THEN
         NUM_OF_CENTS = NUCIND
      !else do global fitting
      ELSE
         NUM_OF_CENTS = 1
      END IF
C
      !Loop over symmetry unique centers
      DO ICENTB = 1, NUM_OF_CENTS
         !Loop over symmetry unique centers
         DO ICENTA = 1, ICENTB
            CALL FIT_DENSITY_OLD(ICENTA,ICENTB,FMAT,DMAT,WMAT,WVEC,NDMT,
     &                        ISYMDM,IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                        WKOD,WORK(KFREE),LFREE,CFIT,
     &                        WORK(KINDAB),IPRINT)
         END DO 
      END DO 
C
      RETURN
      END
C  /* Deck fit_density */
      SUBROUTINE FIT_DENSITY_OLD(ICENTA,ICENTB,FMAT,DMAT,WMAT,WVEC,NDMT,
     &                        ISYMDM,IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                        WKOD,WORK,LWORK,CFIT,INDAB,IPRINT)
C*****************************************************************************
C    
C     FIT_DENSITY: Fit an atom pair density.
C        
C     Written by Luuk Visscher, november 2003
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*), CFIT(*), INDAB(*)
      DIMENSION KOD(12)
#include "nuclei.h"
#include "symmet.h"
#include "inforb.h"
#include "densfit.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbieri.h"
C
CC    DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)
C
#include "memint.h" 
C
C     We take as first index the fit set functions on either center A or B,
C     and as second index the null function. The remaining two are ordinary
C     basis functions on either center A or B.
C     
      IF (DIATOM) THEN 
        !Index 1
        ICLASSES(1) = 3
        NCENTS(1)   = 2
        ICENTS(1,1) = ICENTA
        ICENTS(2,1) = ICENTB
        !Index 2
        ICLASSES(2) = 0
        NCENTS(2)   = 1
        ICENTS(1,2) = 0
        !Index 3
        ICLASSES(3) = 1
        NCENTS(3)   = 2
        ICENTS(1,3) = ICENTA
        ICENTS(2,3) = ICENTB
        !Index 4
        ICLASSES(4) = 1
        NCENTS(4)   = 2
        ICENTS(1,4) = ICENTA
        ICENTS(2,4) = ICENTB
      ELSE
        !Index 1
        ICLASSES(1) = 3
        NCENTS(1)   = 0
        ICENTS(1,1) = 0
        !Index 2
        ICLASSES(2) = 0
        NCENTS(2)   = 1
        ICENTS(1,2) = 0
        !Index 3
        ICLASSES(3) = 1
        NCENTS(3)   = 0
        ICENTS(1,3) = 0
        !Index 4
        ICLASSES(4) = 1
        NCENTS(4)   = 0
        ICENTS(1,4) = 0
      END IF
      CALL PICK_ACT_BATCH_OLD (NCENTS,ICENTS,ICLASSES,IPRINT)
C     Initialize the weight vector and subsequently fill it by
C     contracting the active integrals with density matrix
C
      CALL DZERO (WVEC,NBASISAUX*NDMT)

      CALL TIMER('START ', TIMSTR, TIMEND)
      IF (DOFCK3) THEN
        CALL GALFCK3_OLD(WVEC,DMAT,NDMT,ISYMDM,IFCTYP,IPRINT,
     &               WORK(KFREE),LFREE)
      ELSE
C
        IOPT = 2
        !To prevent many fittings for a single atom in diatom scheme
        IF (DIATOM.AND.(ICENTA.EQ.ICENTB)) IOPT = -2

        CALL ERI_MKINTS (WVEC,DMAT,NDMT,ISYMDM,IFCTYP,IOPT,IPRINT,
     &       CCFBT,INDXBT,KOD,FCKTMP,
     &       WKOD,WORK,LWORK)
      ENDIF
      CALL TIMER('GALPHA', TIMSTR, TIMEND)

cfp
c     print at most first 5000 elements of galpha vector to see if densfit works
      write(lupri,*)'galpha vector check for DOFCK3 ',DOFCK3
      write(lupri,*)'NBASISAUX ',NBASISAUX
      NRANGEEEE = NBASISAUX
      if (NBASISAUX.gt.5000) then
         NRANGEEEE = 5000
         write(lupri,*)'only first 5000 elemnts will be printed'
      end if
      do ii = 1,NRANGEEEE
         write(lupri,*) ii, WVEC(ii)
      end do
cfp

C     
C     Do some printing in case we're really desparate
C
      IF (IPRINT.GT.10) THEN
         IF (DIATOM) THEN
          WRITE (LUPRI,'(//8X,A,2I6//)') 'Fit atom pair:',ICENTA,ICENTB
         END IF
         CALL HEADER('Full weight matrix',-1)
         CALL OUTPUT(WMAT,1,NBASISAUX,1,NBASISAUX,NBASISAUX,NBASISAUX,
     &            1,LUPRI)
         CALL HEADER('Full weight vector',-1)
         CALL OUTPUT(WVEC,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
C
C     Determine how many fit functions we have (should sit on either
C     A or B) and make an index array to do gather/scatter. The indices
C     need to be shifted let them start at 1.
C
      CALL NFUN_ACT_BATCH_OLD (1,N_FIT,INDAB)
      CALL DF_INDEX_SHIFT(INDAB,N_FIT,INDAB,N_FIT,N_FIT,ICLASSES,1,1)
C
C     Allocate memory for temporary storage of A,B subblocks
C
      CALL MEMGET ('REAL',KJAB,N_FIT*N_FIT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KCFTAB,N_FIT,WORK,KFREE,LFREE)
C
C     Extract the weight vector for the pairs
C
      INDNULL = 1
      CALL EXTR_ACT_BATCH_OLD (WVEC,NBASISAUX,WORK(KCFTAB),N_FIT,
     &                     INDAB,N_FIT,INDNULL,1,1)
C
C     Extract the weight matrix from the full weight matrix.
C
      CALL EXTR_ACT_BATCH_OLD (WMAT,NBASISAUX,WORK(KJAB),N_FIT,
     &                     INDAB,N_FIT,INDAB,N_FIT,1)
C
C     Printing these matrices may sometimes be helpful
C
      IF (IPRINT.GT.4) THEN
         IF (DIATOM) THEN
          WRITE (LUPRI,'(//8X,A,2I6/)')'Fitting atom pair',ICENTA,ICENTB
         END IF
         CALL HEADER('Weight matrix',-1)
         CALL OUTPUT(WORK(KJAB),1,N_FIT,1,N_FIT,N_FIT,N_FIT,1,LUPRI)
         CALL HEADER('Weight vector',-1)
         CALL OUTPUT(WORK(KCFTAB),1,N_FIT,1,1,N_FIT,1,1,LUPRI)
      ENDIF
C
C     Solve the equation for the fit coefficients
C
      CALL TIMER('START ', TIMSTR, TIMEND)
      CALL DPOSV ( 'L', N_FIT, 1, WORK(KJAB), max(1,N_FIT),
     &             WORK(KCFTAB), max(1,N_FIT), INFO )
      CALL TIMER('LINSOL', TIMSTR, TIMEND)
C
C     We again have some information that may be useful
C
      IF (IPRINT.GT.4) THEN
         WRITE (LUPRI,'(//8X,A,2I6//)') 'INFO code from DPOSV',INFO
         CALL HEADER('Fit coefficients',-1)
         CALL OUTPUT(WORK(KCFTAB),1,N_FIT,1,1,N_FIT,1,1,LUPRI)
      ENDIF
C
C     Add the fit coefficients to the full array
C
      CALL EXTR_ACT_BATCH_OLD (CFIT,NBASISAUX,WORK(KCFTAB),N_FIT,
     &                     INDAB,N_FIT,INDNULL,1,3)
C
      CALL MEMREL('FIT_DENSITY',WORK,KJAB,KJAB,KFREE,LFREE)
C
      RETURN
      END
C  /* Deck galfck3 */
      SUBROUTINE GALFCK3_OLD(CVEC,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                  WORK,LWORK)
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
CAREFUL
      LOGICAL FIXTHRSAVE
CAREFUL
C
      DIMENSION CVEC(*), DMAT(*), ISYMDM(*), IFCTYP(*), 
     &          WORK(*)
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
CAREFUL
#include "fck3thr.h"
CAREFUL
#include "priunit.h"
#include "densfit.h"
C
#include "dummy.h"
#include "memint.h"
C
      CALL QENTER('GALFCK3_OLD')

c      USEERI = RUNERI
c      DO I = 1, NDMT
c        USEERI = USEERI .AND. ((IFCTYP(I).EQ.13).OR.(IFCTYP(I).EQ.3))
c      END DO
c      
      
      DO IDMT=1,NDMT
        ISTR = N2BASX*(IDMT - 1) + 1
        ISTRAUX = NBASISAUX*(IDMT - 1) + 1

cSimen There be a problem here
        IF (NBASISAUX.GT.N2BASX) CALL QUIT('Error: NBASISAUX > N2BASX.')
cSimen

        IF ( MOD(IFCTYP(IDMT),10) .EQ. 3 ) THEN
C         RHSFIT build of Fab (Fitted 2nd electron part of Fock-matrix)
C         LHSFIT build of gamma_alpha (GAL) vector. 
          RHSFIT = .FALSE.
          LHSFIT = .TRUE.
CAREFUL
          FIXTHRSAVE = FIXTHR
          THRFCK3SAVE = THRFCK3
          FIXTHR = .TRUE.
          THRFCK3 = 1.0D-11
CAREFUL
          CALL ACT_ALL_BATCHES
          CALL FCK3DRV(CVEC(ISTRAUX),DMAT(ISTR),1,DUMMY,DUMMY,DUMMY,
     &                 WORK,LWORK)
CAREFUL
          FIXTHR = FIXTHRSAVE
          THRFCK3 = THRFCK3SAVE
CAREFUL

        END IF
      END DO

      CALL QEXIT('GALFCK3_OLD')
      RETURN
      END
C  /* Deck pick_act_batch */
      SUBROUTINE PICK_ACT_BATCH_OLD (NCENTS,ICENTS,ICLASSES,IPRINT)
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
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "shells.h"
#include "erisel.h"
#include "aobtch.h"
      DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)
C
      DO I = 1, NAOBCH
         ICENT = NCNTBT(I)
         ICLASS = KCLSBT(I)
         DO J = 1, 4
C           The shell must be of the right class
            IF (ICLASS.EQ.ICLASSES(J)) THEN
C              The center must be active (ncents=0 means all centers active)
               IF (NCENTS(J).EQ.0) THEN
                  ACTVBT(I,J) = .TRUE.
               ELSEIF (NCENTS(J).EQ.1.OR.NCENTS(J).EQ.2) THEN
                  ACTVBT(I,J) = .FALSE.
                  DO K = 1, NCENTS(J)
                     IF (ICENTS(K,J).EQ.ICENT) ACTVBT(I,J) = .TRUE.
                  ENDDO
               ELSE
                  CALL QUIT ("Wrong value for NCENTS in PICK_ACT_BATCH")
               ENDIF
            ELSE
               ACTVBT(I,J) = .FALSE.
            ENDIF
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
         WRITE (LUPRI,'(1X,A,/)') ' ACTVBT in PICKAO:'
         DO I = 1, NAOBCH
            WRITE (LUPRI,'(10X,I5,5X,4L5)') I, (ACTVBT(I,J),J=1,4)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck nfun_act_batch */
      SUBROUTINE NFUN_ACT_BATCH_OLD (IND,N_FIT,INDEX)
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
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "shells.h"
#include "erisel.h"
#include "aobtch.h"
C
      DIMENSION INDEX(*)
C
      N_FIT = 0
      J = 0
      DO ISH = 1, NAOBCH
        DO I = 1, NORBBT(ISH)
           J = J + 1
           IF (ACTVBT(ISH,IND)) THEN
              N_FIT = N_FIT + 1
              INDEX(N_FIT) = J
           ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck extr_act_batch */
      SUBROUTINE EXTR_ACT_BATCH_OLD (FULL_MAT,N_FULL,RED_MAT,N_RED,
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
C  /* Deck read_fit2int */
      SUBROUTINE READ_FIT2INT(WMAT,WORK,LWORK)
C
C     Wrapper routine for read of the full set of 2-index fit integrals
C
#include "implicit.h"
#include "iratdef.h"
#include "inforb.h"
#include "cbieri.h"
      DIMENSION WMAT(*),WORK(LWORK)
      LBUF = LBFINP
      KIBUF = 1
      KTOP  = KIBUF + 4*LBUF
      KLEFT = LWORK - KTOP + 1
      CALL READ_FIT2INT_2(WMAT,WORK(KIBUF),LBUF,
     &                    WORK(KTOP),KLEFT)
      RETURN
      END
C  /* Deck read_fit2int_2 */
      SUBROUTINE READ_FIT2INT_2(WMAT,IBUF,LBUF,WRK,LWORK)
C
C     Read the full set of 2-index fit integrals to memory.
C     The indices are shifted so that the first fit function is placed at position 1
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
      INTEGER A, B, C, D
      DIMENSION WMAT(NBASISAUX,NBASISAUX),
     &          IBUF(4,LBUF), MBAS(8), WRK(LWORK)
      DIMENSION ICLASSES(4)
#include "ibtfun.h"
C
      KFREE = 1
      LFREE = LWORK
C
C     We will read integrals of the type (nf|nf)
C     where |f> is a fit function and |n> the null function
      ICLASSES(1) = 0
      ICLASSES(2) = 3
      ICLASSES(3) = 0
      ICLASSES(4) = 3
C
C     Two-electron part
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'FIT_INTS','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM, (MBAS(I),I=1,8), KBUF, NIBUF, NBITS
      JBUF = 0
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('INTE',KINT,LENINT,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
C
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
      CALL DZERO(WMAT,NBASISAUX*NBASISAUX)
      JBUF = 608
  200 CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUINTA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LUINTA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,IBUF,NINT)
         IF (NINT .EQ. 0) GO TO 200
         IF (NINT .LT. 0) GO TO 400
         IF (NINT .GT. 0) THEN
C           Indices are relative to the fulllist, shift them to the aux. list
            CALL DF_INDEX_SHIFT(IBUF,LBUF,IBUF,LBUF,NINT,ICLASSES,4,2)
            DO 300 I = 1, NINT
               GINT  = WRK(KINT + I - 1)
               A = IBUF(1,I)
               B = IBUF(2,I)
               C = IBUF(3,I)
               D = IBUF(4,I)
               WMAT(B,D) = GINT
               WMAT(D,B) = GINT
  300       CONTINUE
         END IF
      GO TO 200
  400 CONTINUE
      CALL GPCLOSE(LUINTA,'KEEP')
      RETURN
      END
C  /* Deck df_index_shift */
      SUBROUTINE DF_INDEX_SHIFT(IBIN,LBIN,IBUF,LBUF,NINTS,ICLASSES,
     &                          NIND,MODE)
C*****************************************************************************
C    
C     Shifts the index for the auxilliary fit function such that the first fit
C     functions get index 1 and the last nbasisaux
C
C     MODE 1 : To be used in the Fock matrix construction where the labels are
C              organized as in IBIN
C     MODE 2 : To be used in the read of integrals were the labels are 
C              organized as in IBUF
C
C     ICLASSES : class type for ech of the four indices
C     IBIN/IBUF: integral labels that are to be shifted
C     LBIN/LBUF: length of the buffer
C     NIND:      number of indices that are to be shifted
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "shells.h"
#include "nuclei.h"
#include "erisel.h"
#include "aobtch.h"
C
      DIMENSION IBIN(LBIN,NIND),IBUF(NIND,LBUF),IND_SHIFT(NIND)
      DIMENSION ICLASSES(NIND)
C
C     Determine the shifts for each index
C
      DO J = 1, NIND
         IF (ICLASSES(J).EQ.1) THEN
C           Nonrelativistic or large component functions
            IND_SHIFT(J) = 0
         ELSEIF (ICLASSES(J).EQ.2) THEN
C           Huckel or small component functions
            IND_SHIFT(J) = NLARGE
         ELSEIF (ICLASSES(J).EQ.3) THEN
C           Auxilliary functions for density fitting
            IND_SHIFT(J) = NLARGE + NSMALL
         ELSEIF (ICLASSES(J).EQ.0) THEN
C           The null function comes after the normal aux. functions
            IND_SHIFT(J) = NLARGE + NSMALL + NBASISAUX
         ENDIF
      END DO
C
      IF (MODE.EQ.1) THEN
         DO J = 1, NIND
            IF (IND_SHIFT(J).NE.0) THEN
               DO I = 1, NINTS
                  IBIN(I,J) = IBIN(I,J) - IND_SHIFT(J)
               END DO
            ENDIF
         END DO
      ELSE
         DO I = 1, NINTS
            DO J = 1, NIND
               IBUF(J,I) = IBUF(J,I) - IND_SHIFT(J)
            END DO
         END DO
      ENDIF
C
      RETURN
      END
C  /* Deck index_canon */
      SUBROUTINE INDEX_CANON(IBIN,LBIN,NINTS)
C*****************************************************************************
C    
C     Make pair of indices canonical (I>J)
C
C     IBIN: integral labels that are to be made canonical
C     LBIN: length of the buffer
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
C
      DIMENSION IBIN(LBIN,2)
C
      DO INT = 1, NINTS
         I = IBIN(INT,1)
         J = IBIN(INT,2)
         IBIN(INT,1) = MAX(I,J)
         IBIN(INT,2) = MIN(I,J)
      ENDDO
C
      RETURN
      END
cC  /* Deck dfiinp */
c      SUBROUTINE DFIINP(WORD)
c#include "implicit.h"
c#include "priunit.h"
c#include "densfit.h"
c#include "fck3box.h"
c      PARAMETER ( NTABLE = 13 )
c      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
c      DATA TABLE /'.DIATOM','.DOBOX ','.BOXORB','.BUFSIZ','.CONSTR',
c     &            '.NOROBU','.ATTENU','.GAUDAM','.XXXXXX','.XXXXXX',
c     &            '.XXXXXX','.XXXXXX','.XXXXXX'/
cC
c      WORD1 = WORD
cC
c 100    CONTINUE
c            READ (LUCMD, '(A7)') WORD
c            PROMPT = WORD(1:1)
c
c            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
c               GO TO 100
c            ELSE IF (PROMPT .EQ. '.') THEN
c               DO 200 I = 1, NTABLE
c                  IF (TABLE(I) .EQ. WORD) THEN
c                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13), I
c                  END IF
c 200          CONTINUE
c              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
c     *            '" not recognized in FCK3INP.'
c               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
c               CALL QUIT('Illegal keyword in FCK3INP.')
c 1    CONTINUE
c        !.DIATOM - fit density for pairs of atoms
c        DIATOM     = .TRUE.
c      GO TO 100
c 2    CONTINUE
c        DOBOX = .TRUE.
c      GO TO 100
c 3    CONTINUE
c        READ (LUCMD,*) MAXORBBOX
c        IF ((MAXORBBOX .LE. 1) .OR. (MAXORBBOX .GT. 100000)) THEN
c          WRITE(LUPRI,*) 'Maxorbox input error, maxorbox=', MAXORBOX
c          CALL QUIT('Maxorbox input error')
c        ENDIF
c        DOBOX = .TRUE.
c      GO TO 100
c 4    CONTINUE
c        READ (LUCMD,*) RBUFFER
c        IF (RBUFFER.LT.0.0d0) THEN
c          WRITE(LUPRI,*) 'Rbuffer input error, Rbuffer=', RBUFFER
c          CALL QUIT('Rbuffer input error')
c        ENDIF
c        DOBOX = .TRUE.
c      GO TO 100
c 5    CONTINUE 
c        DF_CONSTR = .TRUE.
c      GO TO 100
c 6    CONTINUE
c        ROBUST_E = .FALSE.
c      GO TO 100
c 7    CONTINUE
c        DF_ATTENU = .TRUE.
c        READ (LUCMD,*) DF_ATTPAR
c        IF (DF_ATTPAR.LT.0.0d0) THEN
c          WRITE(LUPRI,*) 'DF_ATTPAR input error, DF_ATTPAR=', DF_ATTPAR
c          CALL QUIT('DF_ATTPAR input error')
c        ENDIF
c        IF (DF_GAUDMP) THEN
c          WRITE(LUPRI,*) 
c     &     'Cannot have both .ATTENU and .GAUDAM simultaneously'
c          CALL QUIT('Input error: .ATTENU and .GAUDAM simult.')
c        ENDIF
c      GO TO 100
c 8    CONTINUE
c        DF_GAUDMP = .TRUE.
c        READ (LUCMD,*) DF_GAUPAR
c        IF (DF_GAUPAR.LT.0.0d0) THEN
c          WRITE(LUPRI,*) 'DF_GAUPAR input error, DF_GAUPAR=', DF_GAUPAR
c          CALL QUIT('DF_GAUPAR input error')
c        ENDIF
c        IF (DF_ATTENU) THEN
c          WRITE(LUPRI,*) 
c     &     'Cannot have both .ATTENU and .GAUDAM simultaneously'
c          CALL QUIT('Input error: .ATTENU and .GAUDAM simult.')
c        ENDIF
c      GO TO 100
c 9    CONTINUE
c      GO TO 100  
c 10   CONTINUE
c      GO TO 100  
c 11   CONTINUE
c      GO TO 100  
c 12   CONTINUE
c      GO TO 100  
c 13   CONTINUE
c      GO TO 100  
c            ELSE IF (PROMPT .EQ. '*') THEN
c               GO TO 300
c            ELSE
c              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
c     *            '" not recognized in DFIINP.'
c               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
c                 CALL QUIT('Illegal keyword in DFIINP.')
c            END IF
c 300  CONTINUE
c
c*     CALL HEADER('Settings for FCK3 calculation:',0)
c
c*     IF (IPRINT .NE. IPRFDEF) THEN
c*         WRITE (LUPRI,'(A,I5)') ' Print level in FCK3:',IPRFCK3
c*     END IF
c
c*     IF ( FIXTHR ) THEN
c*       WRITE (LUPRI,'(A,I5)') 'Threshold is fixed.'
c*     ELSE
c*       WRITE (LUPRI,'(A,I5)') 'Threshold set dynamically.'
c*     END IF
c
c*     IF ( NOSKIPCL ) THEN
c*       WRITE (LUPRI,'(A,I5)') 
c*    &  'No skipping of the nonclassical integrals (only set to 0).'
c*     END IF
c
c      RETURN
c      END
cC  /* Deck dfiini */
c      SUBROUTINE DFIINI
c#include "implicit.h"
c#include "densfit.h"
c      !.DIATOM - fit density for pairs of atoms (default is the global fit)
c      DOFIT     = .FALSE.
c      DIATOM    = .FALSE.
c      LHSFIT    = .FALSE.
c      RHSFIT    = .FALSE.
c      ROBUST_E  = .TRUE.
c      DOBOX     = .FALSE.
c      DF_CONSTR = .FALSE.
c      DF_GAUDMP = .FALSE.
c      DF_ATTENU = .FALSE.
c      CALL SET_BOX_DEFAULTS()
c      RETURN
c      END
c
