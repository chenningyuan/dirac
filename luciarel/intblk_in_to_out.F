      SUBROUTINE INTBLK_IN_TO_OUT(IPSM,IQSM,IRSM,ISSM,I1234I,
     &           I12I,I34I,NLISTO,I1234O,I12O,I34O,
     &           ICCSMO,IMOSF_SP,ISPSF_MO,NIRRIN,NIRROUT,
     &           ISPITRA,ISPITRA_RC,
     &           NBLOCKO,IBLOCKO,LISTO,INDEXO,ISIGNO,IRIO,NCOMP)
*
* A given symmetry block of two-electron integrals
* is given by IPSM, IQSM, IRSM, ISSM.
* The permutational symmetry of the input integrals
* are defined by I1234I,I12I,I34I,ICCSMO
*
* NLISTO lists of output integrals are to be generated.
* These are determined by I1234O,I12O,I34O.
*
* Obtain the output blocks that this given input block
* will contribute to
*
* NBLOCKO : Number of output blocks
* IBLOCKO : Symmetry labels of output blocks
*           (in double group symmetry)
* ISIGNO  : Sign from permuting input symmetry block
*           (unused!)
* LISTO :  Output list of integral to which each block belongs
* INDEXO : Mapping of index in output block to
*          index in original input symmetry block
* IRIO   : Real or imaginary block
*
*
* The largest allowed number of output blocks is set by
* MXBLKO
*
* The blocks are - if required seperated into real and imaginary
* blocks
*
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL(A-H,O-Z)
*
#include "mxpdim.inc"
#include "input.inc"
*
*. Input
      INTEGER I1234O(NLISTO),I12O(NLISTO),I34O(NLISTO)
      INTEGER ICCSMO(NLISTO)
      INTEGER IMOSF_SP(8,2),ISPSF_MO(8,2)
      INTEGER ISPITRA(4,NLISTO),ISPITRA_RC(NLISTO)
*
*. Note ISPITRA_RC(I) = 1 => This transformation is real
*                     = 2 => This transformation is imaginary
*                     = 3 => Contains both an imaginary and real part
*. Output
      INTEGER IBLOCKO(4,MXBLKO),ISIGNO(MXBLKO)
      INTEGER INDEXO (4,MXBLKO),LISTO(MXBLKO)
      INTEGER IRIO(MXBLKO)
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) ' Welcome to INTBLK_IN_TO_OUT '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' IPSM IQSM IRSM ISSM for input block :',
     &               IPSM,IQSM,IRSM,ISSM
        WRITE(6,'(A,3I3)') ' I1234I,I12I,I34I :',I1234I,I12I,I34I
*
*
        WRITE(6,*) ' Number of output lists ', NLISTO
        WRITE(6,*) ' Output integral lists '
        WRITE(6,*)
        WRITE(6,*)
     &  ' List  I12  I34  I1234 ICC  SP1  SP2  SP3  SP4'
        WRITE(6,*)
     &  '=============================================='
        DO JLISTO = 1, NLISTO
        WRITE(6,'(1H , 9I5)')
     &   JLISTO,I12O(JLISTO),I34O(JLISTO),I1234O(JLISTO),
     &   ICCSMO(JLISTO),
     &   ISPITRA(1,JLISTO),ISPITRA(2,JLISTO),ISPITRA(3,JLISTO),
     &   ISPITRA(4,JLISTO)
        END DO
*
      END IF
*
      NBLOCKO = 0
*. Each input symmetry can in general correspond to two
*. output symmetries ( pointgroupXspinfunction=>double group).
      DO IPSF = 1,2
      DO IQSF = 1,2
*. Well here we are going to insert some resrtrictions.
* We will - for the moment assume that the basic integrals
* must have same spin components for particle 1 and 2.
*. This can be removed or modified later
Ctimo      IF(IPSF.EQ.IQSF) THEN
*
      DO IRSF = 1,2
      DO ISSF = 1,2
Ctimo      IF(IRSF.EQ.ISSF) THEN
        IPSMO = IMOSF_SP(IPSM,IPSF)
        IQSMO = IMOSF_SP(IQSM,IQSF)
        IRSMO = IMOSF_SP(IRSM,IRSF)
        ISSMO = IMOSF_SP(ISSM,ISSF)
        WRITE(6,'(1H ,A,4I3)')
     &  ' IPSMO IQSMO IRSMO ISSMO ',
     &    IPSMO,IQSMO,IRSMO,ISSMO
        IF(IPSMO.NE.0.AND.IQSMO.NE.0.AND.
     &  IRSMO.NE.0.AND.ISSMO.NE.0) THEN
*. Loop over all permutations of input block
          IF(I12I.EQ.0.OR.IPSM.EQ.IQSM) THEN
            J12MX = 1
          ELSE
            J12MX = 2
          END IF
*
          if(I34I.EQ.0.OR.IRSM.EQ.ISSM) then
            J34MX = 1
          ELSE
            J34MX = 2
          END IF
*
          IF(I1234I.EQ.1.OR.(IPSM.EQ.IRSM.AND.IQSM.EQ.ISSM)) THEN
            J1234MX = 1
          ELSE
            J1234MX = 2
          END IF
*
          DO J12 = 1, J12MX
            IF(J12.EQ.1) THEN
              JPSMO = IPSMO
              JQSMO = IQSMO
              J1 = 1
              J2 = 2
            ELSE
              JPSMO = IQSMO
              JQSMO = IPSMO
              J1 = 2
              J2 = 1
            END IF
            DO J34 = 1, J34MX
              IF(J34.EQ.1) THEN
               JRSMO = IRSMO
               JSSMO = ISSMO
               J3 = 3
               J4 = 4
              ELSE
               JRSMO = ISSMO
               JSSMO = IRSMO
               J3 = 4
               J4 = 3
              END IF
              DO J1234 = 1, J1234MX
                IF(J1234.EQ.1) THEN
                  KPSMO = JPSMO
                  KQSMO = JQSMO
                  KRSMO = JRSMO
                  KSSMO = JSSMO
*
                  K1 = J1
                  K2 = J2
                  K3 = J3
                  K4 = J4
*
                ELSE
                  KPSMO = JRSMO
                  KRSMO = JPSMO
                  KQSMO = JSSMO
                  KSSMO = JQSMO
*
                  K1 = J3
                  K2 = J4
                  K3 = J1
                  K4 = J2
                END IF
*
                IF(IPRINT.GE.2) THEN
                  WRITE(6,*) ' Permutation of input block'
                  WRITE(6,*) ' K1 K2 K3 K4',K1,K2,K3,K4
                  WRITE(6,*) ' KPSMO KQSMO KRSMO KSSMO ',
     &            KPSMO,KQSMO,KRSMO,KSSMO
                END IF
*
                DO JLISTO = 1, NLISTO
*. Does symmetries satisfy restrictions imposed on
* output list JLISTO
                  I_AM_OKAY = 1
                  IF(IPSMO.LT.IQSMO.AND.I12O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
*. complex conjugation symmetry can be used to reduce
*. output lists to IPSMO .ge. IQSMO
                  IF(KPSMO.LT.KQSMO.AND.ICCSMO(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  IF(KRSMO.LT.KSSMO.AND.I34O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  KPQSMO = (KPSMO-1)*NIRRIN+KQSMO
                  KRSSMO = (KRSMO-1)*NIRRIN+KSSMO
                  IF(KPQSMO.LT.KRSSMO.AND.I1234O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  IF(I_AM_OKAY.EQ.1) THEN
*. Congratulations, you have passed all tests,
*  and you have now the right to be computed
*.
*
* One or two blocks due to real/imaginary separation
*
Ctimo                    IF(ISPITRA_RC(ISPITRA(1,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(2,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(3,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(4,JLISTO)).EQ.3) THEN
*
*  Somethings wrong here. The following seems to be appropriate:
                    if (ISPITRA_RC(JLISTO).eq.3) then
                        IMIX = 1
                    ELSE
                        IMIX = 0
                    END IF
                    IF(IMIX.EQ.0) THEN
Ctimo                      NIM = 0
Ctimo                      IF(ISPITRA_RC(ISPITRA(1,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(2,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(3,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(4,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
                       NIM = 0
                       if (ISPITRA(1,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(2,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(3,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(4,JLISTO).eq.2) NIM = NIM + 1
                    END IF
                    IF(IMIX.EQ.1) THEN
*  This transformation is complex
                      NCOMP = 2
                    ELSE
*  This transformation is real or imaginary
                      NCOMP = 1
                      IF(NIM.EQ.1.OR.NIM.EQ.3) THEN
*  This transformation is imaginary
                        IRI = 2
                      ELSE
*  This transformation is real
                        IRI = 1
                      END IF
                    END IF
                    DO JCOMP = 1, NCOMP
                      NBLOCKO = NBLOCKO + 1
                      IF(NBLOCKO.GT.MXBLKO) THEN
                        WRITE(6,*)
     &                  '  TROUBLE in INTBLK_IN_TO_OUT'
                        WRITE(6,*)
     &                  ' Too many output blocks '
                        WRITE(6,*)
     &                ' Raise MXBLKO from ', MXBLKO
                        Call Abend2( '  TROUBLE in INTBLK_IN_TO_OUT' )
                      END IF
                      IF(NCOMP.EQ.1) THEN
                        IRIO(NBLOCKO) = IRI
                      ELSE
                        IRIO(NBLOCKO) = JCOMP
                      END IF
                      LISTO(NBLOCKO) = JLISTO
*. symmetries of output block
                      IBLOCKO(1,NBLOCKO)=IPSMO
                      IBLOCKO(2,NBLOCKO)=IQSMO
                      IBLOCKO(3,NBLOCKO)=IRSMO
                      IBLOCKO(4,NBLOCKO)=ISSMO
*. Index of input block to index of output block
                      INDEXO(1,NBLOCKO) = K1
                      INDEXO(2,NBLOCKO) = K2
                      INDEXO(3,NBLOCKO) = K3
                      INDEXO(4,NBLOCKO) = K4
                    END DO
*                   /\ end  of loop over real/imag components
                  END IF
*                      (end i_am_okay check)
                END DO
*               /\ End of loop over output lists of integrals
*              given spatial symmetry
              END DO
*             /\ End of J1234 loop
            END DO
*           /\ End  of J34 Loop
          END DO
*         /\ End of J12 Loop
*
        END IF
*       /\ End if for nonvanishing symmetries
Ctimo      END IF
*     /\ End if for identity of spincomponents for particle 2
      END DO
      END DO
Ctimo      END IF
*     /\ End if for identity of spincomponents for particle 1
      END DO
      END DO
*     /\ End of loop over the spinor symmetries associated with
*     given spatial symmetry
*
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) ' Output from INTBLK_IN_TO_OUT'
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Number of output blocks', NBLOCKO
        WRITE(6,*)
        WRITE(6,*)
     &  ' Block  list   symmetries     indeces      real/imag  '
        WRITE(6,*)
     &  ' ======================================================='
        DO JBLOCKO = 1, NBLOCKO
          WRITE(6,'(1H ,I4,5x,I2,5X,4I2,5x,4I2,9X,I2)')
     &    JBLOCKO,LISTO(JBLOCKO),(IBLOCKO(I,JBLOCKO),I=1,4),
     &    (INDEXO(I,JBLOCKO),I=1,4),IRIO(JBLOCKO)
        END DO
      END IF
*
      RETURN
      END
