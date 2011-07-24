      SUBROUTINE CONTRACTION_SOLVER(WORK,KFREE,LFREE)
*
* 1. Find all possible contractions between the cluster operator and
* Hamiltonian and the Intermediates.
*
*
      DIMENSION WORK(*)
*
* 1. Find individual contractions
*
      CALL FIND_CONTRACTIONS(WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*
* 2. Find correcpondence between H and M. This is to add them up later
*
      CALL H_AND_M(WORK,KFREE,LFREE)
      CALL MEMCHK_KRCC(WORK)
*
* 3. Setup contraction order (algorithm dependant)
*
      CALL CONTRACTION_ORDER(WORK,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE FIND_CONTRACTIONS(WORK,KFREE,LFREE)
*
* Will find all possible contractions between T and H and M
*
* Now for all possible contractions in the (p,h) scheme
* Shown as Operator --Type of contarction with T --> New operator
* H for Hamiltonian
* M for intermediate (easier to see)
* T for Cluster operator (no index needed)
* Indices are free of the ph definition
*
* This a long and not very pretty and repetitative routine
*
* All possible contraction with the current scheme needs to be
* considered. These are:
*
* H22 --> M12,M20,M11,M02,M10,M01,T
*
* H21 --> M11,M10,M01,T
*
* H12 --> M02,M10,M01,T
*
* H20 --> M10,T
*
* H11 --> M01,T
*
* H02 --> M01,T
*
* H10 --> T
*
* H01 --> T
*
* M12 --> M02
*
* M20 --> T
*
* M11 --> M01,T
*
* M02 --> M01,T
*
* M10 --> T
*
* M01 --> T
*
*
* The various contractions will be listed as listed above
*
* Perhaps the longest subroutine in the world. Lasse 2010
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "interm.inc"
#include "ipoist8.inc"
*
      DIMENSION WORK(*)
*
      INTEGER NUMCOM,IFACTOR(4),INUMDUM(4),INUMCON(4,4)
*
      NTEST = 100
*
      IONE = 1
* All Hamiltonian contractions
      IHM = 1
*
* All H22 Contractions start
* 
*
* H22 --10--> M12
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 1
      NHOLE = 0
*
      NH22TOM12 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM12),
     &                           WORK(KINT12_IDXS),WORK(KINT12_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM12 = NH22TOM12 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM12,NH22TOM12,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T12,NH22TOM12,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM12,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM12,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM12,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM12,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM12,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM12,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM12 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM12),
     &                           WORK(KINT12_IDXS),WORK(KINT12_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM12+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM12+NH22TOM12),
     &                           WORK(KNH22T12+NH22TOM12),
     &                           WORK(KCOMBFACNH22TOM12+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM12+IH-1),NH22TOM12+1,IONE)
        NH22TOM12 = NH22TOM12 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM12+IH-1),NH22TOM12,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM12+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM12+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M12'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM12,WORK(KNH22TOM12),
     &                          WORK(KNH22T12),
     &                          WORK(KCOMBNH22TOM12),
     &                          WORK(KPERMNH22TOM12),
     &                          WORK(KFACNH22TOM12),
     &                          WORK(KCOMBFACNH22TOM12),
     &                          WORK(KSNH22TOM12),WORK(KFNH22TOM12),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM12),NINTER12,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --02--> M20
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 0
      NHOLE = 2
*
      NH22TOM20 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM20),
     &                           WORK(KINT20_IDXS),WORK(KINT20_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM20 = NH22TOM20 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM20,NH22TOM20,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T20,NH22TOM20,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM20,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM20,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM20,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM20,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM20,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM20,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM20 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM20),
     &                           WORK(KINT20_IDXS),WORK(KINT20_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM20+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM20+NH22TOM20),
     &                           WORK(KNH22T20+NH22TOM20),
     &                           WORK(KCOMBFACNH22TOM20+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM20+IH-1),NH22TOM20+1,IONE)
        NH22TOM20 = NH22TOM20 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM20+IH-1),NH22TOM20,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM20+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM20+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M20'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM20,WORK(KNH22TOM20),
     &                          WORK(KNH22T20),
     &                          WORK(KCOMBNH22TOM20),
     &                          WORK(KPERMNH22TOM20),
     &                          WORK(KFACNH22TOM20),
     &                          WORK(KCOMBFACNH22TOM20),
     &                          WORK(KSNH22TOM20),WORK(KFNH22TOM20),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM20),NINTER20,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --11--> M11
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 1
      NHOLE = 1
*
      NH22TOM11 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM11 = NH22TOM11 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM11,NH22TOM11,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T11,NH22TOM11,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM11,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM11,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM11,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM11,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM11,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM11,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM11 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM11+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM11+NH22TOM11),
     &                           WORK(KNH22T11+NH22TOM11),
     &                           WORK(KCOMBFACNH22TOM11+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM11+IH-1),NH22TOM11+1,IONE)
        NH22TOM11 = NH22TOM11 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM11+IH-1),NH22TOM11,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM11+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM11+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M11'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM11,WORK(KNH22TOM11),
     &                          WORK(KNH22T11),
     &                          WORK(KCOMBNH22TOM11),
     &                          WORK(KPERMNH22TOM11),
     &                          WORK(KFACNH22TOM11),
     &                          WORK(KCOMBFACNH22TOM11),
     &                          WORK(KSNH22TOM11),WORK(KFNH22TOM11),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM11),NINTER11,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --20--> M02
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 2
      NHOLE = 0
*
      NH22TOM02 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM02 = NH22TOM02 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM02,NH22TOM02,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T02,NH22TOM02,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM02,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM02,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM02,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM02,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM02,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM02,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM02 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM02+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM02+NH22TOM02),
     &                           WORK(KNH22T02+NH22TOM02),
     &                           WORK(KCOMBFACNH22TOM02+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM02+IH-1),NH22TOM02+1,IONE)
        NH22TOM02 = NH22TOM02 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM02+IH-1),NH22TOM02,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM02+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM02+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M02'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM02,WORK(KNH22TOM02),
     &                          WORK(KNH22T02),
     &                          WORK(KCOMBNH22TOM02),
     &                          WORK(KPERMNH22TOM02),
     &                          WORK(KFACNH22TOM02),
     &                          WORK(KCOMBFACNH22TOM02),
     &                          WORK(KSNH22TOM02),WORK(KFNH22TOM02),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM02),NINTER02,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --12--> M10
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 1
      NHOLE = 2
*
      NH22TOM10 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM10 = NH22TOM10 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM10,NH22TOM10,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T10,NH22TOM10,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM10,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM10,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM10,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM10,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM10,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM10,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM10 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM10+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM10+NH22TOM10),
     &                           WORK(KNH22T10+NH22TOM10),
     &                           WORK(KCOMBFACNH22TOM10+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM10+IH-1),NH22TOM10+1,IONE)
        NH22TOM10 = NH22TOM10 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM10+IH-1),NH22TOM10,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM10+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM10+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M10'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM10,WORK(KNH22TOM10),
     &                          WORK(KNH22T10),
     &                          WORK(KCOMBNH22TOM10),
     &                          WORK(KPERMNH22TOM10),
     &                          WORK(KFACNH22TOM10),
     &                          WORK(KCOMBFACNH22TOM10),
     &                          WORK(KSNH22TOM10),WORK(KFNH22TOM10),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM10),NINTER10,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --21--> M01
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 2
      NHOLE = 1
*
      NH22TOM01 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOM01 = NH22TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOM01,NH22TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22T01,NH22TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOM01,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOM01,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOM01,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOM01,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOM01,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOM01,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOM01 = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNH22TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOM01+NH22TOM01),
     &                           WORK(KNH22T01+NH22TOM01),
     &                           WORK(KCOMBFACNH22TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOM01+IH-1),NH22TOM01+1,IONE)
        NH22TOM01 = NH22TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOM01+IH-1),NH22TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO M01'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOM01,WORK(KNH22TOM01),
     &                          WORK(KNH22T01),
     &                          WORK(KCOMBNH22TOM01),
     &                          WORK(KPERMNH22TOM01),
     &                          WORK(KFACNH22TOM01),
     &                          WORK(KCOMBFACNH22TOM01),
     &                          WORK(KSNH22TOM01),WORK(KFNH22TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H22 --22--> T
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 2
      NHOLE = 2
*
      NH22TOT = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH22TOT = NH22TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH22TOT,NH22TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH22TT,NH22TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH22TOT,NH22,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH22TOT,NH22,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH22TOT,4*NH22,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH22TOT,4*NH22,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH22TOT,NH22,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH22TOT,NH22,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH22TOT = 0
*
      DO IH = 1,NH22
        CALL GET_OP_STUPID(WORK(KNH22),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH22TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH22TOT+NH22TOT),
     &                           WORK(KNH22TT+NH22TOT),
     &                           WORK(KCOMBFACNH22TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH22TOT+IH-1),NH22TOT+1,IONE)
        NH22TOT = NH22TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH22TOT+IH-1),NH22TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH22TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH22TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH22 TO T'
        CALL CONTRACTION_TESTER(NH22,WORK(KNH22),
     &                          NH22TOT,WORK(KNH22TOT),
     &                          WORK(KNH22TT),
     &                          WORK(KCOMBNH22TOT),
     &                          WORK(KPERMNH22TOT),
     &                          WORK(KFACNH22TOT),
     &                          WORK(KCOMBFACNH22TOT),
     &                          WORK(KSNH22TOT),WORK(KFNH22TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H22 contractions
*
*
* All H21 Contractions start
* 
*
* H21 --10--> M11
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 1
      NPART = 1
      NHOLE = 0
*
      NH21TOM11 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH21TOM11 = NH21TOM11 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH21TOM11,NH21TOM11,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH21T11,NH21TOM11,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH21TOM11,NH21,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH21TOM11,NH21,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH21TOM11,4*NH21,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH21TOM11,4*NH21,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH21TOM11,NH21,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH21TOM11,NH21,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH21TOM11 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           NUMCON,WORK(KFACNH21TOM11+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH21TOM11+NH21TOM11),
     &                           WORK(KNH21T11+NH21TOM11),
     &                           WORK(KCOMBFACNH21TOM11+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH21TOM11+IH-1),NH21TOM11+1,IONE)
        NH21TOM11 = NH21TOM11 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH21TOM11+IH-1),NH21TOM11,IONE)
        CALL ISETVC(WORK(KCOMBNH21TOM11+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH21TOM11+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH21 TO M11'
        CALL CONTRACTION_TESTER(NH21,WORK(KNH21),
     &                          NH21TOM11,WORK(KNH21TOM11),
     &                          WORK(KNH21T11),
     &                          WORK(KCOMBNH21TOM11),
     &                          WORK(KPERMNH21TOM11),
     &                          WORK(KFACNH21TOM11),
     &                          WORK(KCOMBFACNH21TOM11),
     &                          WORK(KSNH21TOM11),WORK(KFNH21TOM11),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM11),NINTER11,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H21 --11--> M10
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 1
      NPART = 1
      NHOLE = 1
*
      NH21TOM10 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH21TOM10 = NH21TOM10 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH21TOM10,NH21TOM10,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH21T10,NH21TOM10,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH21TOM10,NH21,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH21TOM10,NH21,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH21TOM10,4*NH21,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH21TOM10,4*NH21,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH21TOM10,NH21,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH21TOM10,NH21,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH21TOM10 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,WORK(KFACNH21TOM10+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH21TOM10+NH21TOM10),
     &                           WORK(KNH21T10+NH21TOM10),
     &                           WORK(KCOMBFACNH21TOM10+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH21TOM10+IH-1),NH21TOM10+1,IONE)
        NH21TOM10 = NH21TOM10 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH21TOM10+IH-1),NH21TOM10,IONE)
        CALL ISETVC(WORK(KCOMBNH21TOM10+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH21TOM10+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH21 TO M10'
        CALL CONTRACTION_TESTER(NH21,WORK(KNH21),
     &                          NH21TOM10,WORK(KNH21TOM10),
     &                          WORK(KNH21T10),
     &                          WORK(KCOMBNH21TOM10),
     &                          WORK(KPERMNH21TOM10),
     &                          WORK(KFACNH21TOM10),
     &                          WORK(KCOMBFACNH21TOM10),
     &                          WORK(KSNH21TOM10),WORK(KFNH21TOM10),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM10),NINTER10,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H21 --20--> M01
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 1
      NPART = 2
      NHOLE = 0
*
      NH21TOM01 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH21TOM01 = NH21TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH21TOM01,NH21TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH21T01,NH21TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH21TOM01,NH21,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH21TOM01,NH21,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH21TOM01,4*NH21,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH21TOM01,4*NH21,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH21TOM01,NH21,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH21TOM01,NH21,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH21TOM01 = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNH21TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH21TOM01+NH21TOM01),
     &                           WORK(KNH21T01+NH21TOM01),
     &                           WORK(KCOMBFACNH21TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH21TOM01+IH-1),NH21TOM01+1,IONE)
        NH21TOM01 = NH21TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH21TOM01+IH-1),NH21TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNH21TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH21TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH21 TO M01'
        CALL CONTRACTION_TESTER(NH21,WORK(KNH21),
     &                          NH21TOM01,WORK(KNH21TOM01),
     &                          WORK(KNH21T01),
     &                          WORK(KCOMBNH21TOM01),
     &                          WORK(KPERMNH21TOM01),
     &                          WORK(KFACNH21TOM01),
     &                          WORK(KCOMBFACNH21TOM01),
     &                          WORK(KSNH21TOM01),WORK(KFNH21TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H21 --21--> T
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 1
      NPART = 2
      NHOLE = 1
*
      NH21TOT = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH21TOT = NH21TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH21TOT,NH21TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH21TT,NH21TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH21TOT,NH21,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH21TOT,NH21,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH21TOT,4*NH21,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH21TOT,4*NH21,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH21TOT,NH21,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH21TOT,NH21,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH21TOT = 0
*
      DO IH = 1,NH21
        CALL GET_OP_STUPID(WORK(KNH21),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH21TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH21TOT+NH21TOT),
     &                           WORK(KNH21TT+NH21TOT),
     &                           WORK(KCOMBFACNH21TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH21TOT+IH-1),NH21TOT+1,IONE)
        NH21TOT = NH21TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH21TOT+IH-1),NH21TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH21TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH21TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH21 TO T'
        CALL CONTRACTION_TESTER(NH21,WORK(KNH21),
     &                          NH21TOT,WORK(KNH21TOT),
     &                          WORK(KNH21TT),
     &                          WORK(KCOMBNH21TOT),
     &                          WORK(KPERMNH21TOT),
     &                          WORK(KFACNH21TOT),
     &                          WORK(KCOMBFACNH21TOT),
     &                          WORK(KSNH21TOT),WORK(KFNH21TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H21 contractions
*
*
* All H12 Contractions start
*
*
* H12 --10--> M02
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 2
      NPART = 1
      NHOLE = 0
*
      NH12TOM02 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH12TOM02 = NH12TOM02 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH12TOM02,NH12TOM02,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH12T02,NH12TOM02,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH12TOM02,NH12,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH12TOM02,NH12,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH12TOM02,4*NH12,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH12TOM02,4*NH12,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH12TOM02,NH12,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH12TOM02,NH12,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH12TOM02 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,WORK(KFACNH12TOM02+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH12TOM02+NH12TOM02),
     &                           WORK(KNH12T02+NH12TOM02),
     &                           WORK(KCOMBFACNH12TOM02+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH12TOM02+IH-1),NH12TOM02+1,IONE)
        NH12TOM02 = NH12TOM02 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH12TOM02+IH-1),NH12TOM02,IONE)
        CALL ISETVC(WORK(KCOMBNH12TOM02+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH12TOM02+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH12 TO M02'
        CALL CONTRACTION_TESTER(NH12,WORK(KNH12),
     &                          NH12TOM02,WORK(KNH12TOM02),
     &                          WORK(KNH12T02),
     &                          WORK(KCOMBNH12TOM02),
     &                          WORK(KPERMNH12TOM02),
     &                          WORK(KFACNH12TOM02),
     &                          WORK(KCOMBFACNH12TOM02),
     &                          WORK(KSNH12TOM02),WORK(KFNH12TOM02),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM02),NINTER02,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H12 --02--> M10
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 2
      NPART = 0
      NHOLE = 2
*
      NH12TOM10 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH12TOM10 = NH12TOM10 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH12TOM10,NH12TOM10,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH12T10,NH12TOM10,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH12TOM10,NH12,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH12TOM10,NH12,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH12TOM10,4*NH12,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH12TOM10,4*NH12,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH12TOM10,NH12,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH12TOM10,NH12,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH12TOM10 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,WORK(KFACNH12TOM10+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH12TOM10+NH12TOM10),
     &                           WORK(KNH12T10+NH12TOM10),
     &                           WORK(KCOMBFACNH12TOM10+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH12TOM10+IH-1),NH12TOM10+1,IONE)
        NH12TOM10 = NH12TOM10 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH12TOM10+IH-1),NH12TOM10,IONE)
        CALL ISETVC(WORK(KCOMBNH12TOM10+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH12TOM10+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH12 TO M10'
        CALL CONTRACTION_TESTER(NH12,WORK(KNH12),
     &                          NH12TOM10,WORK(KNH12TOM10),
     &                          WORK(KNH12T10),
     &                          WORK(KCOMBNH12TOM10),
     &                          WORK(KPERMNH12TOM10),
     &                          WORK(KFACNH12TOM10),
     &                          WORK(KCOMBFACNH12TOM10),
     &                          WORK(KSNH12TOM10),WORK(KFNH12TOM10),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM10),NINTER10,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H12 --11--> M01
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 2
      NPART = 1
      NHOLE = 1
*
      NH12TOM01 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH12TOM01 = NH12TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH12TOM01,NH12TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH12T01,NH12TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH12TOM01,NH12,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH12TOM01,NH12,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH12TOM01,4*NH12,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH12TOM01,4*NH12,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH12TOM01,NH12,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH12TOM01,NH12,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH12TOM01 = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNH12TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH12TOM01+NH12TOM01),
     &                           WORK(KNH12T01+NH12TOM01),
     &                           WORK(KCOMBFACNH12TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH12TOM01+IH-1),NH12TOM01+1,IONE)
        NH12TOM01 = NH12TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH12TOM01+IH-1),NH12TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNH12TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH12TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH12 TO M01'
        CALL CONTRACTION_TESTER(NH12,WORK(KNH12),
     &                          NH12TOM01,WORK(KNH12TOM01),
     &                          WORK(KNH12T01),
     &                          WORK(KCOMBNH12TOM01),
     &                          WORK(KPERMNH12TOM01),
     &                          WORK(KFACNH12TOM01),
     &                          WORK(KCOMBFACNH12TOM01),
     &                          WORK(KSNH12TOM01),WORK(KFNH12TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H12 --12--> T
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 2
      NPART = 1
      NHOLE = 2
*
      NH12TOT = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH12TOT = NH12TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH12TOT,NH12TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH12TT,NH12TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH12TOT,NH12,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH12TOT,NH12,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH12TOT,4*NH12,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH12TOT,4*NH12,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH12TOT,NH12,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH12TOT,NH12,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH12TOT = 0
*
      DO IH = 1,NH12
        CALL GET_OP_STUPID(WORK(KNH12),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH12TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH12TOT+NH12TOT),
     &                           WORK(KNH12TT+NH12TOT),
     &                           WORK(KCOMBFACNH12TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH12TOT+IH-1),NH12TOT+1,IONE)
        NH12TOT = NH12TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH12TOT+IH-1),NH12TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH12TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH12TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH12 TO T'
        CALL CONTRACTION_TESTER(NH12,WORK(KNH12),
     &                          NH12TOT,WORK(KNH12TOT),
     &                          WORK(KNH12TT),
     &                          WORK(KCOMBNH12TOT),
     &                          WORK(KPERMNH12TOT),
     &                          WORK(KFACNH12TOT),
     &                          WORK(KCOMBFACNH12TOT),
     &                          WORK(KSNH12TOT),WORK(KFNH12TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H21 contractions
*
*
* All H20 Contractions start
*
*
* H20 --10--> M10
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 0
      NPART = 1
      NHOLE = 0
*
      NH20TOM10 = 0
*
      DO IH = 1,NH20
        CALL GET_OP_STUPID(WORK(KNH20),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH20TOM10 = NH20TOM10 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH20TOM10,NH20TOM10,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH20T10,NH20TOM10,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH20TOM10,NH20,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH20TOM10,NH20,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH20TOM10,4*NH20,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH20TOM10,4*NH20,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH20TOM10,NH20,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH20TOM10,NH20,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH20TOM10 = 0
*
      DO IH = 1,NH20
        CALL GET_OP_STUPID(WORK(KNH20),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM10),
     &                           WORK(KINT10_IDXS),WORK(KINT10_IDXF),
     &                           NUMCON,WORK(KFACNH20TOM10+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH20TOM10+NH20TOM10),
     &                           WORK(KNH20T10+NH20TOM10),
     &                           WORK(KCOMBFACNH20TOM10+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH20TOM10+IH-1),NH20TOM10+1,IONE)
        NH20TOM10 = NH20TOM10 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH20TOM10+IH-1),NH20TOM10,IONE)
        CALL ISETVC(WORK(KCOMBNH20TOM10+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH20TOM10+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH20 TO M10'
        CALL CONTRACTION_TESTER(NH20,WORK(KNH20),
     &                          NH20TOM10,WORK(KNH20TOM10),
     &                          WORK(KNH20T10),
     &                          WORK(KCOMBNH20TOM10),
     &                          WORK(KPERMNH20TOM10),
     &                          WORK(KFACNH20TOM10),
     &                          WORK(KCOMBFACNH20TOM10),
     &                          WORK(KSNH20TOM10),WORK(KFNH20TOM10),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM10),NINTER10,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H20 --20--> T
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 0
      NPART = 2
      NHOLE = 0
*
      NH20TOT = 0
*
      DO IH = 1,NH20
        CALL GET_OP_STUPID(WORK(KNH20),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH20TOT = NH20TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH20TOT,NH20TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH20TT,NH20TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH20TOT,NH20,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH20TOT,NH20,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH20TOT,4*NH20,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH20TOT,4*NH20,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH20TOT,NH20,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH20TOT,NH20,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH20TOT = 0
*
      DO IH = 1,NH20
        CALL GET_OP_STUPID(WORK(KNH20),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH20TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH20TOT+NH20TOT),
     &                           WORK(KNH20TT+NH20TOT),
     &                           WORK(KCOMBFACNH20TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH20TOT+IH-1),NH20TOT+1,IONE)
        NH20TOT = NH20TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH20TOT+IH-1),NH20TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH20TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH20TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH20 TO T'
        CALL CONTRACTION_TESTER(NH20,WORK(KNH20),
     &                          NH20TOT,WORK(KNH20TOT),
     &                          WORK(KNH20TT),
     &                          WORK(KCOMBNH20TOT),
     &                          WORK(KPERMNH20TOT),
     &                          WORK(KFACNH20TOT),
     &                          WORK(KCOMBFACNH20TOT),
     &                          WORK(KSNH20TOT),WORK(KFNH20TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H20 contractions
*
*
* All H11 Contractions start
*
*
* H11 --10--> M01
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 1
      NPART = 1
      NHOLE = 0
*
      NH11TOM01 = 0
*
      DO IH = 1,NH11
        CALL GET_OP_STUPID(WORK(KNH11),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH11TOM01 = NH11TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH11TOM01,NH11TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH11T01,NH11TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH11TOM01,NH11,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH11TOM01,NH11,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH11TOM01,4*NH11,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH11TOM01,4*NH11,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH11TOM01,NH11,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH11TOM01,NH11,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH11TOM01 = 0
*
      DO IH = 1,NH11
        CALL GET_OP_STUPID(WORK(KNH11),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNH11TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH11TOM01+NH11TOM01),
     &                           WORK(KNH11T01+NH11TOM01),
     &                           WORK(KCOMBFACNH11TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH11TOM01+IH-1),NH11TOM01+1,IONE)
        NH11TOM01 = NH11TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH11TOM01+IH-1),NH11TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNH11TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH11TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH11 TO M01'
        CALL CONTRACTION_TESTER(NH11,WORK(KNH11),
     &                          NH11TOM01,WORK(KNH11TOM01),
     &                          WORK(KNH11T01),
     &                          WORK(KCOMBNH11TOM01),
     &                          WORK(KPERMNH11TOM01),
     &                          WORK(KFACNH11TOM01),
     &                          WORK(KCOMBFACNH11TOM01),
     &                          WORK(KSNH11TOM01),WORK(KFNH11TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H11 --11--> T
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 1
      NPART = 1
      NHOLE = 1
*
      NH11TOT = 0
*
      DO IH = 1,NH11
        CALL GET_OP_STUPID(WORK(KNH11),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH11TOT = NH11TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH11TOT,NH11TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH11TT,NH11TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH11TOT,NH11,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH11TOT,NH11,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH11TOT,4*NH11,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH11TOT,4*NH11,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH11TOT,NH11,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH11TOT,NH11,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH11TOT = 0
*
      DO IH = 1,NH11
        CALL GET_OP_STUPID(WORK(KNH11),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH11TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH11TOT+NH11TOT),
     &                           WORK(KNH11TT+NH11TOT),
     &                           WORK(KCOMBFACNH11TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH11TOT+IH-1),NH11TOT+1,IONE)
        NH11TOT = NH11TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH11TOT+IH-1),NH11TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH11TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH11TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH11 TO T'
        CALL CONTRACTION_TESTER(NH11,WORK(KNH11),
     &                          NH11TOT,WORK(KNH11TOT),
     &                          WORK(KNH11TT),
     &                          WORK(KCOMBNH11TOT),
     &                          WORK(KPERMNH11TOT),
     &                          WORK(KFACNH11TOT),
     &                          WORK(KCOMBFACNH11TOT),
     &                          WORK(KSNH11TOT),WORK(KFNH11TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H11 contractions
*
*
* All H02 Contractions start
*
*
* H02 --01--> M01
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 2
      NPART = 0
      NHOLE = 1
*
      NH02TOM01 = 0
*
      DO IH = 1,NH02
        CALL GET_OP_STUPID(WORK(KNH02),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH02TOM01 = NH02TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH02TOM01,NH02TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH02T01,NH02TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH02TOM01,NH02,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH02TOM01,NH02,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH02TOM01,4*NH02,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH02TOM01,4*NH02,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH02TOM01,NH02,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH02TOM01,NH02,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH02TOM01 = 0
*
      DO IH = 1,NH02
        CALL GET_OP_STUPID(WORK(KNH02),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNH02TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH02TOM01+NH02TOM01),
     &                           WORK(KNH02T01+NH02TOM01),
     &                           WORK(KCOMBFACNH02TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH02TOM01+IH-1),NH02TOM01+1,IONE)
        NH02TOM01 = NH02TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNH02TOM01+IH-1),NH02TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNH02TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH02TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH02 TO M01'
        CALL CONTRACTION_TESTER(NH02,WORK(KNH02),
     &                          NH02TOM01,WORK(KNH02TOM01),
     &                          WORK(KNH02T01),
     &                          WORK(KCOMBNH02TOM01),
     &                          WORK(KPERMNH02TOM01),
     &                          WORK(KFACNH02TOM01),
     &                          WORK(KCOMBFACNH02TOM01),
     &                          WORK(KSNH02TOM01),WORK(KFNH02TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* H02 --02--> T
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 2
      NPART = 0
      NHOLE = 2
*
      NH02TOT = 0
*
      DO IH = 1,NH02
        CALL GET_OP_STUPID(WORK(KNH02),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH02TOT = NH02TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH02TOT,NH02TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH02TT,NH02TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH02TOT,NH02,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH02TOT,NH02,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH02TOT,4*NH02,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH02TOT,4*NH02,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH02TOT,NH02,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH02TOT,NH02,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH02TOT = 0
*
      DO IH = 1,NH02
        CALL GET_OP_STUPID(WORK(KNH02),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH02TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH02TOT+NH02TOT),
     &                           WORK(KNH02TT+NH02TOT),
     &                           WORK(KCOMBFACNH02TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH02TOT+IH-1),NH02TOT+1,IONE)
        NH02TOT = NH02TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH02TOT+IH-1),NH02TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH02TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH02TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH02 TO T'
        CALL CONTRACTION_TESTER(NH02,WORK(KNH02),
     &                          NH02TOT,WORK(KNH02TOT),
     &                          WORK(KNH02TT),
     &                          WORK(KCOMBNH02TOT),
     &                          WORK(KPERMNH02TOT),
     &                          WORK(KFACNH02TOT),
     &                          WORK(KCOMBFACNH02TOT),
     &                          WORK(KSNH02TOT),WORK(KFNH02TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H02 contractions
*
*
* All H10 Contractions start
*
*
* H10 --10--> T
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 0
      NPART = 1
      NHOLE = 0
*
      NH10TOT = 0
*
      DO IH = 1,NH10
        CALL GET_OP_STUPID(WORK(KNH10),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH10TOT = NH10TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH10TOT,NH10TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH10TT,NH10TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH10TOT,NH10,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH10TOT,NH10,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH10TOT,4*NH10,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH10TOT,4*NH10,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH10TOT,NH10,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH10TOT,NH10,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH10TOT = 0
*
      DO IH = 1,NH10
        CALL GET_OP_STUPID(WORK(KNH10),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH10TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH10TOT+NH10TOT),
     &                           WORK(KNH10TT+NH10TOT),
     &                           WORK(KCOMBFACNH10TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH10TOT+IH-1),NH10TOT+1,IONE)
        NH10TOT = NH10TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH10TOT+IH-1),NH10TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH10TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH10TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH10 TO T'
        CALL CONTRACTION_TESTER(NH10,WORK(KNH10),
     &                          NH10TOT,WORK(KNH10TOT),
     &                          WORK(KNH10TT),
     &                          WORK(KCOMBNH10TOT),
     &                          WORK(KPERMNH10TOT),
     &                          WORK(KFACNH10TOT),
     &                          WORK(KCOMBFACNH10TOT),
     &                          WORK(KSNH10TOT),WORK(KFNH10TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H10 contractions
*
*
* All H01 Contractions start
*
*
* H01 --01--> T
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 1
      NPART = 0
      NHOLE = 1
*
      NH01TOT = 0
*
      DO IH = 1,NH01
        CALL GET_OP_STUPID(WORK(KNH01),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NH01TOT = NH01TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNH01TOT,NH01TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNH01TT,NH01TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNH01TOT,NH01,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNH01TOT,NH01,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNH01TOT,4*NH01,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNH01TOT,4*NH01,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNH01TOT,NH01,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNH01TOT,NH01,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NH01TOT = 0
*
      DO IH = 1,NH01
        CALL GET_OP_STUPID(WORK(KNH01),IH,IHTP)
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXS),
     &                           NUMCON,WORK(KFACNH01TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNH01TOT+NH01TOT),
     &                           WORK(KNH01TT+NH01TOT),
     &                           WORK(KCOMBFACNH01TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNH01TOT+IH-1),NH01TOT+1,IONE)
        NH01TOT = NH01TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNH01TOT+IH-1),NH01TOT,IONE)
        CALL ISETVC(WORK(KCOMBNH01TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNH01TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NH01 TO T'
        CALL CONTRACTION_TESTER(NH01,WORK(KNH01),
     &                          NH01TOT,WORK(KNH01TOT),
     &                          WORK(KNH01TT),
     &                          WORK(KCOMBNH01TOT),
     &                          WORK(KPERMNH01TOT),
     &                          WORK(KFACNH01TOT),
     &                          WORK(KCOMBFACNH01TOT),
     &                          WORK(KSNH01TOT),WORK(KFNH01TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX),IHM)
      END IF
*
* End all H01 contractions
*
* All intermediates are now contracted
      IHM = 2
*
* All M12 Contractions start
*
*
* M12 --10--> M02
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 2
      NPART = 1
      NHOLE = 0
*
      NM12TOM02 = 0
*
      DO IH = 1,NINTER12
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM12),IH,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM12TOM02 = NM12TOM02 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM12TOM02,NM12TOM02,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM12T02,NM12TOM02,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM12TOM02,NINTER12,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM12TOM02,NINTER12,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM12TOM02,4*NINTER12,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM12TOM02,4*NINTER12,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM12TOM02,NINTER12,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM12TOM02,NINTER12,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM12TOM02 = 0
*
      DO IH = 1,NINTER12
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM12),IH,WORK(KINTM02),
     &                           WORK(KINT02_IDXS),WORK(KINT02_IDXF),
     &                           NUMCON,WORK(KFACNM12TOM02+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM12TOM02+NM12TOM02),
     &                           WORK(KNM12T02+NM12TOM02),
     &                           WORK(KCOMBFACNM12TOM02+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM12TOM02+IH-1),NM12TOM02+1,IONE)
        NM12TOM02 = NM12TOM02 + NCONNECTTOT
        CALL ISETVC(WORK(KFNM12TOM02+IH-1),NM12TOM02,IONE)
        CALL ISETVC(WORK(KCOMBNM12TOM02+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM12TOM02+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM12 TO M02'
        CALL CONTRACTION_TESTER(NINTER12,IDUMMY,
     &                          NM12TOM02,WORK(KNM12TOM02),
     &                          WORK(KNM12T02),
     &                          WORK(KCOMBNM12TOM02),
     &                          WORK(KPERMNM12TOM02),
     &                          WORK(KFACNM12TOM02),
     &                          WORK(KCOMBFACNM12TOM02),
     &                          WORK(KSNM12TOM02),WORK(KFNM12TOM02),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM02),NINTER02,
     &                          WORK(KINTM12),IHM)
      END IF
*
* End all M12 contractions
*
*
* All M20 Contractions start
*
*
* M20 --20--> T
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 0
      NPART = 2
      NHOLE = 0
*
      NM20TOT = 0
*
      DO IH = 1,NINTER20
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM20),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM20TOT = NM20TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM20TOT,NM20TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM20TT,NM20TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM20TOT,NINTER20,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM20TOT,NINTER20,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM20TOT,4*NINTER20,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM20TOT,4*NINTER20,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM20TOT,NINTER20,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM20TOT,NINTER20,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM20TOT = 0
*
      DO IH = 1,NINTER20
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM20),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,WORK(KFACNM20TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM20TOT+NM20TOT),
     &                           WORK(KNM20TT+NM20TOT),
     &                           WORK(KCOMBFACNM20TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM20TOT+IH-1),NM20TOT+1,IONE)
        NM20TOT = NM20TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNM20TOT+IH-1),NM20TOT,IONE)
        CALL ISETVC(WORK(KCOMBNM20TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM20TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM20 TO T'
        CALL CONTRACTION_TESTER(NINTER20,IDUMMY,
     &                          NM20TOT,WORK(KNM20TOT),
     &                          WORK(KNM20TT),
     &                          WORK(KCOMBNM20TOT),
     &                          WORK(KPERMNM20TOT),
     &                          WORK(KFACNM20TOT),
     &                          WORK(KCOMBFACNM20TOT),
     &                          WORK(KSNM20TOT),WORK(KFNM20TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM20),IHM)
      END IF
*
* End all M20 contractions
*
*
* All M11 Contractions start
*
*
* M11 --10--> M01
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 1
      NPART = 1
      NHOLE = 0
*
      NM11TOM01 = 0
*
      DO IH = 1,NINTER11
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM11),IH,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM11TOM01 = NM11TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM11TOM01,NM11TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM11T01,NM11TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM11TOM01,NINTER11,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM11TOM01,NINTER11,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM11TOM01,4*NINTER11,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM11TOM01,4*NINTER11,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM11TOM01,NINTER11,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM11TOM01,NINTER11,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM11TOM01 = 0
*
      DO IH = 1,NINTER11
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM11),IH,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNM11TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM11TOM01+NM11TOM01),
     &                           WORK(KNM11T01+NM11TOM01),
     &                           WORK(KCOMBFACNM11TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM11TOM01+IH-1),NM11TOM01+1,IONE)
        NM11TOM01 = NM11TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNM11TOM01+IH-1),NM11TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNM11TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM11TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM11 TO M01'
        CALL CONTRACTION_TESTER(NINTER11,IDUMMY,
     &                          NM11TOM01,WORK(KNM11TOM01),
     &                          WORK(KNM11T01),
     &                          WORK(KCOMBNM11TOM01),
     &                          WORK(KPERMNM11TOM01),
     &                          WORK(KFACNM11TOM01),
     &                          WORK(KCOMBFACNM11TOM01),
     &                          WORK(KSNM11TOM01),WORK(KFNM11TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KINTM11),IHM)
      END IF
*
* M11 --11--> T
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 1
      NPART = 1
      NHOLE = 1
*
      NM11TOT = 0
*
      DO IH = 1,NINTER11
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM11),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM11TOT = NM11TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM11TOT,NM11TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM11TT,NM11TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM11TOT,NINTER11,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM11TOT,NINTER11,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM11TOT,4*NINTER11,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM11TOT,4*NINTER11,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM11TOT,NINTER11,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM11TOT,NINTER11,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM11TOT = 0
*
      DO IH = 1,NINTER11
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM11),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,WORK(KFACNM11TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM11TOT+NM11TOT),
     &                           WORK(KNM11TT+NM11TOT),
     &                           WORK(KCOMBFACNM11TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM11TOT+IH-1),NM11TOT+1,IONE)
        NM11TOT = NM11TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNM11TOT+IH-1),NM11TOT,IONE)
        CALL ISETVC(WORK(KCOMBNM11TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM11TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM11 TO T'
        CALL CONTRACTION_TESTER(NINTER11,IDUMMY,
     &                          NM11TOT,WORK(KNM11TOT),
     &                          WORK(KNM11TT),
     &                          WORK(KCOMBNM11TOT),
     &                          WORK(KPERMNM11TOT),
     &                          WORK(KFACNM11TOT),
     &                          WORK(KCOMBFACNM11TOT),
     &                          WORK(KSNM11TOT),WORK(KFNM11TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM11),IHM)
      END IF
*
* End all M11 contractions
*
*
* All M02 Contractions start
*
*
* M02 --01--> M01
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 2
      NPART = 0
      NHOLE = 1
*
      NM02TOM01 = 0
*
      DO IH = 1,NINTER02
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM02),IH,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM02TOM01 = NM02TOM01 + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM02TOM01,NM02TOM01,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM02T01,NM02TOM01,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM02TOM01,NINTER02,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM02TOM01,NINTER02,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM02TOM01,4*NINTER02,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM02TOM01,4*NINTER02,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM02TOM01,NINTER02,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM02TOM01,NINTER02,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM02TOM01 = 0
*
      DO IH = 1,NINTER02
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM02),IH,WORK(KINTM01),
     &                           WORK(KINT01_IDXS),WORK(KINT01_IDXF),
     &                           NUMCON,WORK(KFACNM02TOM01+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM02TOM01+NM02TOM01),
     &                           WORK(KNM02T01+NM02TOM01),
     &                           WORK(KCOMBFACNM02TOM01+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM02TOM01+IH-1),NM02TOM01+1,IONE)
        NM02TOM01 = NM02TOM01 + NCONNECTTOT
        CALL ISETVC(WORK(KFNM02TOM01+IH-1),NM02TOM01,IONE)
        CALL ISETVC(WORK(KCOMBNM02TOM01+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM02TOM01+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM02 TO M01'
        CALL CONTRACTION_TESTER(NINTER02,IDUMMY,
     &                          NM02TOM01,WORK(KNM02TOM01),
     &                          WORK(KNM02T01),
     &                          WORK(KCOMBNM02TOM01),
     &                          WORK(KPERMNM02TOM01),
     &                          WORK(KFACNM02TOM01),
     &                          WORK(KCOMBFACNM02TOM01),
     &                          WORK(KSNM02TOM01),WORK(KFNM02TOM01),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),NINTER01,
     &                          WORK(KINTM02),IHM)
      END IF
*
* M02 --02--> T
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 2
      NPART = 0
      NHOLE = 2
*
      NM02TOT = 0
*
      DO IH = 1,NINTER02
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM02),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM02TOT = NM02TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM02TOT,NM02TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM02TT,NM02TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM02TOT,NINTER02,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM02TOT,NINTER02,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM02TOT,4*NINTER02,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM02TOT,4*NINTER02,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM02TOT,NINTER02,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM02TOT,NINTER02,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM02TOT = 0
*
      DO IH = 1,NINTER02
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM02),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,WORK(KFACNM02TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM02TOT+NM02TOT),
     &                           WORK(KNM02TT+NM02TOT),
     &                           WORK(KCOMBFACNM02TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM02TOT+IH-1),NM02TOT+1,IONE)
        NM02TOT = NM02TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNM02TOT+IH-1),NM02TOT,IONE)
        CALL ISETVC(WORK(KCOMBNM02TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM02TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM02 TO T'
        CALL CONTRACTION_TESTER(NINTER02,IDUMMY,
     &                          NM02TOT,WORK(KNM02TOT),
     &                          WORK(KNM02TT),
     &                          WORK(KCOMBNM02TOT),
     &                          WORK(KPERMNM02TOT),
     &                          WORK(KFACNM02TOT),
     &                          WORK(KCOMBFACNM02TOT),
     &                          WORK(KSNM02TOT),WORK(KFNM02TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM02),IHM)
      END IF
*
* End all M02 contractions
*
*
* All M10 Contractions start
*
*
* M10 --10--> T
*
      INI = 0
      NPARTTOT = 1
      NHOLETOT = 0
      NPART = 1
      NHOLE = 0
*
      NM10TOT = 0
*
      DO IH = 1,NINTER10
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM10),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM10TOT = NM10TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM10TOT,NM10TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM10TT,NM10TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM10TOT,NINTER10,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM10TOT,NINTER10,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM10TOT,4*NINTER10,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM10TOT,4*NINTER10,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM10TOT,NINTER10,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM10TOT,NINTER10,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM10TOT = 0
*
      DO IH = 1,NINTER10
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM10),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,WORK(KFACNM10TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM10TOT+NM10TOT),
     &                           WORK(KNM10TT+NM10TOT),
     &                           WORK(KCOMBFACNM10TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM10TOT+IH-1),NM10TOT+1,IONE)
        NM10TOT = NM10TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNM10TOT+IH-1),NM10TOT,IONE)
        CALL ISETVC(WORK(KCOMBNM10TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM10TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM10 TO T'
        CALL CONTRACTION_TESTER(NINTER10,IDUMMY,
     &                          NM10TOT,WORK(KNM10TOT),
     &                          WORK(KNM10TT),
     &                          WORK(KCOMBNM10TOT),
     &                          WORK(KPERMNM10TOT),
     &                          WORK(KFACNM10TOT),
     &                          WORK(KCOMBFACNM10TOT),
     &                          WORK(KSNM10TOT),WORK(KFNM10TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM10),IHM)
      END IF
*
* End all M10 contractions
*
*
* All M01 Contractions start
*
*
* M01 --01--> T
*
      INI = 0
      NPARTTOT = 0
      NHOLETOT = 1
      NPART = 0
      NHOLE = 1
*
      NM01TOT = 0
*
      DO IH = 1,NINTER01
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM01),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,IFACTOR,
     &                           NCONNECTTOT,
     &                           IDUMMY,
     &                           IDUMMY,
     &                           INUMDUM,
     &                           WORK,KFREE,LFREE)
        NM01TOT = NM01TOT + NCONNECTTOT 
      END DO
*
* Now allocate arrays for storage of info
*
* The resultant operator
      CALL MEMGET('INTE',KNM01TOT,NM01TOT,WORK,KFREE,LFREE)
* The cluster operator contracted with
      CALL MEMGET('INTE',KNM01TT,NM01TOT,WORK,KFREE,LFREE)
* The number of connections for a given HTYPE
      CALL MEMGET('INTE',KCOMBNM01TOT,NINTER01,WORK,KFREE,LFREE)
* The number of permutations. 
* How can NPART and NHOLE be contracted for HTYPE
      CALL MEMGET('INTE',KPERMNM01TOT,NINTER01,WORK,KFREE,LFREE)
* The permutation factor for the contractions
      CALL MEMGET('INTE',KFACNM01TOT,4*NINTER01,WORK,KFREE,LFREE)
* The number of ways with a given permutation factor      
      CALL MEMGET('INTE',KCOMBFACNM01TOT,4*NINTER01,WORK,KFREE,LFREE)
* The start for an operator (Only used for internal checking)
      CALL MEMGET('INTE',KSNM01TOT,NINTER01,WORK,KFREE,LFREE)
* The finish (Only used for internal checking)
      CALL MEMGET('INTE',KFNM01TOT,NINTER01,WORK,KFREE,LFREE)
*
* Now for the actual connections
*
      INI = 1
      NM01TOT = 0
*
      DO IH = 1,NINTER01
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KINTM01),IH,WORK(KLSOBEX_CC),
     &                           WORK(KT_IDXS),WORK(KT_IDXF),
     &                           NUMCON,WORK(KFACNM01TOT+4*(IH-1)),
     &                           NCONNECTTOT,
     &                           WORK(KNM01TOT+NM01TOT),
     &                           WORK(KNM01TT+NM01TOT),
     &                           WORK(KCOMBFACNM01TOT+4*(IH-1)),
     &                           WORK,KFREE,LFREE)
        CALL ISETVC(WORK(KSNM01TOT+IH-1),NM01TOT+1,IONE)
        NM01TOT = NM01TOT + NCONNECTTOT
        CALL ISETVC(WORK(KFNM01TOT+IH-1),NM01TOT,IONE)
        CALL ISETVC(WORK(KCOMBNM01TOT+IH-1),NCONNECTTOT,IONE)
        CALL ISETVC(WORK(KPERMNM01TOT+IH-1),NUMCON,IONE)
      END DO
* Let us see what we got out
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Contractions from NM01 TO T'
        CALL CONTRACTION_TESTER(NINTER01,IDUMMY,
     &                          NM01TOT,WORK(KNM01TOT),
     &                          WORK(KNM01TT),
     &                          WORK(KCOMBNM01TOT),
     &                          WORK(KPERMNM01TOT),
     &                          WORK(KFACNM01TOT),
     &                          WORK(KCOMBFACNM01TOT),
     &                          WORK(KSNM01TOT),WORK(KFNM01TOT),
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                          WORK(KINTM01),IHM)
      END IF
*
* End all M01 contractions
*
*
* End of the line. Was this not long and tedious?
* Should you find or know a longer subroutine then please write to me at
* lks@chem.au.dk and i will try to make this one longer.
*
      RETURN
      END
*
      SUBROUTINE GET_OP_STUPID(IOP,IH,IHTP)
* Gettting the operator the stupid way
#include "implicit.inc"
*
      INTEGER IOP(*)
*
      IHTP = IOP(IH)
*
      RETURN
      END
*
      SUBROUTINE H_AND_M(WORK,KFREE,LFREE)
*
* Overall aim is to reduce the number of contractions by only doing
* those that are connected!
*
* List all operators from intermediates which can be connected to T
* 
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "ctcc.inc"
#include "interm.inc"
*
      DIMENSION WORK(*)
*
      WRITE(6,*) ' Welcome to H and M '
*
* First find for every H-type the same M-type so they can be added later
*
* H12 + M12
*
* How do I want to store it??? Assume for now every H is included in M
* and then check if that is true!
* 
      print*,' matching 12'
      ISIMHT12 = NH12 !MIN(NINTER12,NH12)
      CALL MEMGET('INTE',KM12TOH12,ISIMHT12,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH12),NH12,WORK(KINTM12),NINTER12,
     &     WORK(KINT12_IDXS),WORK(KINT12_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM12TOH12),WORK(KLSOBEX))
* H11 + M11
      print*,' matching 11 '
      ISIMHT11 = NH11 !MIN(NINTER11,NH11)
      CALL MEMGET('INTE',KM11TOH11,ISIMHT11,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH11),NH11,WORK(KINTM11),NINTER11,
     &     WORK(KINT11_IDXS),WORK(KINT11_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM11TOH11),WORK(KLSOBEX))
* H02 + M02
      print*,' matching 02'
      ISIMHT02 = NH02 !MIN(NINTER02,NH02)
      CALL MEMGET('INTE',KM02TOH02,ISIMHT02,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH02),NH02,WORK(KINTM02),NINTER02,
     &     WORK(KINT02_IDXS),WORK(KINT02_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM02TOH02),WORK(KLSOBEX))
* H10 + M10
      print*,' matching 10'
      ISIMHT10 = NH10 !MIN(NINTER10,NH10)
      CALL MEMGET('INTE',KM10TOH10,ISIMHT10,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH10),NH10,WORK(KINTM10),NINTER10,
     &     WORK(KINT10_IDXS),WORK(KINT10_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM10TOH10),WORK(KLSOBEX))
* H01 + M01
      print*,' matching 01'
      ISIMHT01 = NH01 !MIN(NINTER01,NH01)
      CALL MEMGET('INTE',KM01TOH01,ISIMHT01,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH01),NH01,WORK(KINTM01),NINTER01,
     &     WORK(KINT01_IDXS),WORK(KINT01_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM01TOH01),WORK(KLSOBEX))
* H00 + T
      print*,' matching 00'
      ISIMHT00 = NH00 !MIN(NINTER00,NH00)
      CALL MEMGET('INTE',KM00TOH00,ISIMHT00,WORK,KFREE,LFREE)
      CALL OP_COMPARE(WORK(KNH00),NH00,WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &     WORK(KT_IDXS),WORK(KT_IDXF),MX_EXC_LEVEL,IMAXKRFLIP,
     &     WORK(KM00TOH00),WORK(KLSOBEX))
      RETURN
      END
*
      SUBROUTINE CONTRACTION_ORDER(WORK,KFREE,LFREE)
*
* Will combine the possible contractions from FIND_CONTRACTIONS.
*
* This routine is algorithm dependant
*
* Algorithm in the simplified version goes as following
*
* H20 --T-->
*
* H22 --T--> M02+H02 --T-->
*
* H22 --T--> M11+H11 --T-->
*
* H22 --T-->
*
* H21 --T-->
*
* H12 --T-->
*
* H00 
*
* H22 + H21 + H12 + H20 --T--> M10+H10 --T-->
*
* H22 + H21 + H12 + H11 + H02 --T--> M01+H01 --T-->
*
* H22 + H12 --T--> M02 --T--> M01 --T-->
*
* H22 + H21 --T--> M11 --T--> M01 --T-->
*
* H22 --T--> M12 --T--> M02 --T-->
*                           \
*                            \T--> M01 --T-->
*
* No intersection of two or more sets is needed in finding the
* global permutation parameter for the intermediates since the only time
* a global parameter is needed is when two contractions of the same
* ph-type follows each other. Here the first contraction is multiplied
* by 1/2. This simple global permutation factor stems from the
* insistance on only going "one way" in the multiple possible ways a
* given number of operators can be contracted.
*
* Will create a list of all contractions which can be looped over or??
* List already exist since all contractions are known! Will here in the
* first run not start to eliminate intermediate or final contractions*
* List already exist since all contractions are known! Will here in the
* first run not start to eliminate intermediate or final contractions.
* Therefore this will first be an empty subroutine!!!
*
*
#include "implicit.inc"
*
*
*
*
*
      RETURN
      END
