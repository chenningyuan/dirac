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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

MODULE Emb_Frozen_Orbitals

  IMPLICIT NONE
  PRIVATE

  Public :: init_Emb_Frozen_Orbitals, dealloc_Emb_Frozen_Orbitals 
  Public :: Get_Rho_Frozen_Orbitals, Get_Grad_Rho_Frozen_Orbitals

  REAL*8, Dimension(:,:)  , Allocatable, Save :: Frozen_Density_Mat
  REAL*8, Dimension(:)    , Allocatable, Save :: BUF

CONTAINS

  SUBROUTINE init_Emb_Frozen_Orbitals ()
    IMPLICIT NONE

#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"

    INTEGER :: ierr

    Allocate (Frozen_Density_Mat(NTBAS(0), NTBAS(0)), STAT=ierr )
    IF (ierr /= 0) THEN
        CALL QUIT ("init_Emb_Frozen_Orbitals: Cannot allocate memory for frozen density matrix")
    END IF

    Allocate ( BUF( NTBAS(0) ), STAT=ierr )
    IF (ierr /= 0) THEN
        CALL QUIT ("init_Emb_Frozen_Orbitals: Cannot allocate memory for BUF")
    END IF

    CALL init_Frozen_DMat ()

  END SUBROUTINE init_Emb_Frozen_Orbitals




  SUBROUTINE init_Frozen_DMat ()
    IMPLICIT NONE

#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    INTEGER :: irrep, i
    INTEGER :: n_read                                         ! number of orbitals, that have been read in

    REAL*8,  DIMENSION ( N2BBASXQ ) :: temp_frozen_dmat       ! one dimensional, to make old parts of the code happy
        
    REAL*8,  DIMENSION ( NTBAS(0), NTBAS(0), NZ ) :: temp_frozen_dmat_2

    INTEGER, DIMENSION ( NFROT )                     :: frozen_orb_nums
    REAL*8,  DIMENSION ( NTBAS(0), NZ )              :: orbital_buf
    REAL*8,  DIMENSION ( NTBAS(0), NFROT, NZ )       :: frozen_CMO

    temp_frozen_dmat = 0

    WRITE (LUPRI, *) "CJDEBUG: Entering init_Frozen_DMat"
    WRITE (LUPRI, *) "         NFSYM = ", NFSYM, " NFROT = ", NFROT

    CALL OPNFIL (LUCOEF, 'DFFROZ', 'OLD', 'embed')

    DO irrep = 1, NFSYM 

        WRITE (LUPRI, *) "  irrep = ", irrep, " NFRO(irrep) = ", NFRO(irrep)
        WRITE (LUPRI, *) "    NPSH = ", NPSH(irrep), " NISH = ", NISH(irrep)

        ! Read the numbers of the frozen orbitals (to frozen_orbital_numbers) 
        CALL READI (LUCOEF, NFRO(irrep), frozen_orb_nums)

        WRITE (LUPRI, *) "frozen_orb_nums: ", frozen_orb_nums(1:NFRO(irrep))

        n_read = 0    
        DO i = 1, NFRO(irrep)
            WRITE (LUPRI, *) " i = ", i
            WRITE (LUPRI, *) " frozen_orb_nums(i) = ", frozen_orb_nums(i)

            ! decide if we need this orbital (is it occupied?)
            IF ( (frozen_orb_nums(i) > NPSH(irrep)) .AND. (frozen_orb_nums(i) <= (NPSH(irrep)+NISH(irrep))) ) THEN
                ! count it 
                WRITE (LUPRI, *) "Read orbital! "
                n_read = n_read+1
                ! read it
                CALL READT (LUCOEF, NFBAS(irrep,0)*NZ, orbital_buf)
                frozen_CMO(:, n_read, :) = orbital_buf(:, :)
            END IF
        END DO

        WRITE (LUPRI, *) "  n_read = ", n_read

        IF (n_read > 0) THEN
            ! Construct (part) of the density matrix
            CALL DENST1 ( temp_frozen_dmat( I2BASX(irrep, irrep)+1),                  &
                                               NTBAS(0),        NTBAS(0),             &
                          frozen_CMO,          NFBAS(irrep, 0), NFRO(irrep),          &
                          1, n_read, NFBAS(irrep, 0) ) 
        END IF

!        CALL TITLER('CJDEBUG: temp_frozen_dmat','*',103)
!        DO I = 1,NFSYM
!            WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
!            CALL PRQMAT(temp_frozen_dmat(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0), &
!                            NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
!        END DO

    END DO

    CLOSE (LUCOEF, STATUS='KEEP')

    ! now get it into the correct shape
    temp_frozen_dmat_2 = reshape(temp_frozen_dmat,(/NTBAS(0),NTBAS(0),NZ/)) 

    ! only save the real part, everything else should be 0
    ! and scale it with a factor of 2 (two electrons per frozen orbital)
    Frozen_Density_Mat = 2 * temp_frozen_dmat_2 (:,:,1)

!    CALL TITLER('CJDEBUG: Frozen_Density_Matt','*',103)
!    DO I = 1,NFSYM
!        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
!        CALL OUTPUT(Frozen_Density_Mat(I2BASX(I,I)+1), 1, NFBAS(I,0), 1, NFBAS(I,0), NTBAS(0), NTBAS(0), 1, LUPRI)
!    END DO

  END SUBROUTINE init_Frozen_DMat



  SUBROUTINE dealloc_Emb_Frozen_Orbitals ()
    IMPLICIT NONE

    INTEGER :: ierr

    IF (allocated(Frozen_Density_Mat)) Deallocate (Frozen_Density_Mat, STAT=ierr)
    IF (allocated(BUF)) Deallocate(BUF, STAT=ierr)

  END SUBROUTINE dealloc_Emb_Frozen_Orbitals



  SUBROUTINE Get_Rho_Frozen_Orbitals (Frozen_Rho, GAO)
    IMPLICIT NONE

#include "dcbbas.h"
#include "dgroup.h"
    
    REAL*8,                      Intent(OUT) :: Frozen_Rho
    REAL*8, Dimension(NTBAS(0)), Intent(IN)  :: GAO

    CALL GETRHO (Frozen_Rho, 0, GAO, Frozen_Density_Mat, BUF)

  END SUBROUTINE Get_Rho_Frozen_Orbitals



  SUBROUTINE Get_Grad_Rho_Frozen_Orbitals (Frozen_Grad_Rho, GAO, GAO1)
    IMPLICIT NONE

#include "dcbbas.h"
#include "dgroup.h"

    REAL*8, Dimension(3),           Intent(OUT)   :: Frozen_Grad_Rho
    REAL*8, Dimension(NTBAS(0)),    Intent(IN)    :: GAO
    REAL*8, Dimension(NTBAS(0), 3), Intent(IN)    :: GAO1    

    CALL GETGRHO (Frozen_Grad_Rho, 0, GAO, GAO1, Frozen_Density_Mat, BUF)

  END SUBROUTINE Get_Grad_Rho_Frozen_Orbitals
  


END MODULE Emb_Frozen_Orbitals
