      subroutine spitra2(I12O,I34O,I1234O,
     &                   ICCSMO,ISPITRA,ISPITRA_RC,
     &                   LSPINOR,IRIO,ICCSM,
     &                   NLISTO,I1234I,I12I,I34I,IOPSMIN,IOPSMOUT)
*
* Input
* ======
*
* LBATCHI : Length of batch of input integrals
* LBATCHO : Length of batch of output integrals
*
* I1234I : Permutational symmetry between particle and and two
* I12I  : Permutational symmetry between first two orbital indeces
* I34I  : Permutational symmetry between last two orbital indeces
* NIRR_PN : Number of irreps in input symmmetry group
*
* IOPSMIN : Symmetry of operator in input group
* IOPSMOUT : Symmetry of operator in output group
* NLISTO : Number of output lists
* I1234O : Permutational symmetry between particle one and two
*          for each output list
* I12O  : Permutational symmetry between the first two indeces
*         for each output list
* I34O : Permutational symmetry between the last two indeces
*        for each output list
*
*
* Main routine for the transformartion of
* two-electron integrals to from a.o. basis
* to spinors in double group symmetry
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "symm.inc"
#include "input.inc"
*. Input
      INTEGER I1234O(NLISTO),I12O(NLISTO),I34O(NLISTO)
      INTEGER ICCSMO(NLISTO)
      INTEGER ISPITRA(4,NLISTO)
*
      INTEGER IBLOCKO(4,MXBLKO),LISTO(MXBLKO)
      INTEGER INDEXO(4,MXBLKO),ISIGNO(MXBLKO),IRIO(MXBLKO)
      INTEGER IBL1(6,MXBLKO),IBL2(6,MXBLKO),IBL3(6,MXBLKO)
      INTEGER IBL4(6,MXBLKO)
*
      dimension ISPITRA_RC(NLISTO)
      dimension LSPINOR(N_SPINR_EXP)

C     CALL QENTER('DBGT2')
      IDUM = 0
C     CALL MEMMAN(IDUM,IDUM,'MARK',IDUM,'SPITR2')
*
*
      WRITE(6,*) ' Control has been transferred to SPITRA2'
      WRITE(6,*)
      WRITE(6,*) ' I12I, I34I,I1234I',I12I,I34I,I1234I
      WRITE(6,*) 'ISPITRA : '
      CALL IWRTMA(ISPITRA,4,NLISTO,4,NLISTO)
*
*. Loop over symmetry-unique blocks of input ao integrals
*
      DO IPSM = 1, NIRR_PN
*. Symmetry Range of second index
       IF(I12I.EQ.1.OR.I12I.EQ.-1) THEN
         IMXQSM = IPSM
       ELSE
         IMXQSM = NIRR_PN
       END IF
       DO IQSM = 1, IMXQSM
*. Symmetry range of third index
        IF(I1234I.EQ.1 .OR. I1234I.EQ.-1) THEN
          IMXRSM = IPSM
        ELSE
          IMXRSM = NIRR_PN
        END IF
        DO IRSM = 1, IMXRSM
*. Symmetry range of fourth index
         IMXSSM = NIRR_PN
         IF (IPSM.EQ.IRSM.AND.(I1234I.EQ.1.OR.I1234I.EQ.-1))
     &   IMXSSM = IQSM
         IF(I34I.EQ.1.OR.I34I.EQ.-1) IMXSSM = MIN(IMXSSM,IRSM)
*. Required symmetry of index 4 to produce total operator of
*. correct symmetry
         IPQSM =  IPGMULT(IPSM,IQSM)
         IPQRSM = IPGMULT(IPQSM,IRSM)
         ISSM = IPGMULT(IPQRSM,IOPSMIN)
         WRITE(6,*) ' IPQSM IPQRSM ISSM',IPQSM,IPQRSM,ISSM
         IF(1.LE.ISSM.AND.ISSM.LE.IMXSSM) THEN
*
            IF(IPRINT.GE.2) THEN
              WRITE(6,*) ' Next unique input symmetry block'
              WRITE(6,*) ' IPSM,IQSM,IRSM,ISSM',
     &        IPSM,IQSM,IRSM,ISSM
            END IF
*. Which output  blocks corresponds to this input block
            CALL INTBLK_IN_TO_OUT(IPSM,IQSM,IRSM,ISSM,
     &           I1234I,I12I,I34I,NLISTO,I1234O,I12O,I34O,
     &           ICCSMO,IMOSF_SP,ISPSF_MO,NIRR_PN,NIRR_DG,
     &           ISPITRA,ISPITRA_RC,
     &           NBLOCKO,IBLOCKO(1,1),LISTO(1),INDEXO(1,1),
     &           ISIGNO,IRIO,NCOMP)
*
*     INPUT
*    -------
*           IPSM,IQSM,IRSM,ISSM	: spatial index symmetries
*           I1234I,I12I,I34I	: inlist particle and complexity symmetries
*           NLISTO		: number of output lists
*           ISPITRA		: spinor expansion of index and output list
*           ISPITRA_RC		: r/i/c type of transformation
*
*     OUTPUT
*    --------
*           I1234O,I12O,I34O	: outlist particle and complexity symmetries
*           NBLOCKO		: number of output blocks
*           IBLOCKO		: dg symmetry for index in part. block
*           LISTO		: output list of part. block
*           INDEXO		: permutation of index in this block
*           IRIO		: real or imaginary part of this block
*
*. Intermediate blocks
*. Should not be here but...
*
*
            NRORC = NCOMP
            CALL INTERM_BLK(NBLOCKO,IBLOCKO,INDEXO,LISTO,
     &           IRIO,ISPITRA,ISPITRA_RC,LSPINOR,NRORC,
     &           NLISTO,N_SPINR_EXP,
     &           NBL1,IBL1,NBL2,IBL2,NBL3,IBL3,NBL4,IBL4,
     &           MXBLKO)
*
*     INPUT
*    -------
*           NBLOCKO,IBLOCKO,INDEXO,
*           LISTO,IRIO			: output from intblk_in_to_out
*
*           LSPINOR			: size of spinor expansion sets
*           NRORC			: number of r/i components
*
*     OUTPUT
*    --------
*           NBL1,NBL2,NBL3,NBL4	: number of blocks with X indeces transformed
*           IBL1,IBL2,IBL3,IBL4	: specification of one-index transformed
*				  block as
*				   (1,*) spinor transformation
*				   (2,*) dg symmetry
*				   (3,*) real or imaginary part
*				   (4,*) transformation of index
*					   3 or 4  if X=1,2
*					   1 or 2  if X=3,4
*				   (5,*) # of X+1 index transformed blocks
*				         this X index transf. bl. will
*				         contribute to
*				   (6,*) first X+1 index transf. bl. this
*				         bl. will contribute to
*
         END IF
        END DO
       END DO
      END DO
*
      RETURN
      END
