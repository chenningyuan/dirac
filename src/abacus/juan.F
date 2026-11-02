C/* Deck Prdrao */
      SUBROUTINE PRDRAO(DRAO,NDMAT,INDEX)
C*****************************************************************************
C
C     Print reduced matrix DRAO
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
!      use nuclei_mod, only: NAMN, NUCDEG
      use priunit_mod, only: LUPRI
      use maxorb_mod, only: MXCORB, MXSHEL
!      use blocks_mod, only: LCLASH, MAXSHL, NCNTSH, NHKTSH, NSYMBL
!      use aovec_mod, only: MXAOVC
      implicit none
!#include "implicit.h"
      INTEGER :: I,II,J,JJ,K
#include "aovec.h"
C
!#include "shells.h"
#include "blocks.h"
#include "nuclei.h"
      DIMENSION DRAO(NSYMBL,NSYMBL,NDMAT),INDEX(NSYMBL)
      CHARACTER SPDCAR*1,COMP(2)*1
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     First make pointer from symmetry dependent blocks to
C     symmetry independent blocks
C
      NBL = 0
      DO I = 1,MAXSHL
        NCENTA = NCNTSH(I)
        NDEG = NUCDEG(NCENTA)
        DO IDEG = 1,NDEG
          NBL = NBL + 1
          INDEX(NBL) = I
        ENDDO
      ENDDO
C
C     Print section
C
      CALL HEADER('Output from MKDRAO',-1)
      DO K = 1,NDMAT
        WRITE(LUPRI,'(A,I5)') '* Reduced AO-matrix no.',K
        JEND = MIN(NSYMBL,8)
        DO J = 1,NSYMBL,8
          WRITE(LUPRI,'(12X,8(I4,A1,A3,1X,3A1))')
     &         (JJ,'(',NAMN(NCNTSH(INDEX(JJ)))(1:3),
     &          COMP(LCLASH(INDEX(JJ))),
     &          SPDCAR(NHKTSH(INDEX(JJ))-1),')',JJ = J,JEND)
          DO I = 1,NSYMBL
             WRITE(LUPRI,'(I4,A1,A3,1X,3A1,1P,8D12.3)')
     &         I,'(',NAMN(NCNTSH(INDEX(I)))(1:3),
     &         COMP(LCLASH(INDEX(I))),
     &         SPDCAR(NHKTSH(INDEX(I))-1),')',
     &        (DRAO(I,II,K),II = J,JEND)
          ENDDO
          JEND = MIN(JEND+8,NSYMBL)
        ENDDO
      ENDDO
C
      RETURN
      END