      SUBROUTINE CONTRACTION_SOLVER(WORK,KFREE,LFREE)
*
* 1. Find all possible contractions between the cluster operator and
* Hamiltonian and the Intermediates.
*
*
      DIMENSION WORK(*)
*
* 1. Find contractions
*
      CALL FIND_CONTRACTIONS(WORK,KFREE,LFREE)
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
#include "implicit.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "ipoist8.inc"
      DIMENSION WORK(*)
*
* H22 --11--> M11
*
      INI = 0
      NPARTTOT = 2
      NHOLETOT = 2
      NPART = 1
      NHOLE = 1
*
      
      DO IH = 1,NH22
        IHTP = WORK(KNH22 +IH) !Not so pretty
        CALL CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                           WORK(KLSOBEX_CC),NSPOBEX_TPE,
     &                           WORK(KLSOBEX),IHTP,WORK(KINTM11),
     &                           WORK(KINT11_IDXS),WORK(KINT11_IDXF),
     &                           WORK,KFREE,LFREE)
      END DO
*
* M11' --11--> T    (O11' now goes to M11)
*

*
* H21 --10--> M11
*

*
* M11 --10--> M01 
*
* Deallocate for M11

*
* H22 --10--> M12
*
* Allocate for M02 (Not yet in use)
* Allocate for M12

*
* M12 --10--> M02
*
* Deallocate M12

*
* M02 --02--> T
*

*
* M02 --01--> M01
*

*
* H22 --21--> M01
*

*
* H21 --20--> M01
*

*
* H12 --11--> M01
*

*
* H11 --10--> M01
*

*
* H02 --01--> M01
*

*
* M01 --01--> T
*
* Deallocate M01

*
* H22 --12--> M10
*
* Allocate M10

*
* H21 --11--> M10
*

*
* H20 --10--> M10
*

*
* M10 --10--> T
*
* Deallocate M10

*
* H22 --22--> T
*

*
* H21 --21--> T
*

*
* H12 --12--> T
*

*
* H20 --20--> T
*

*
* H11 --11--> T
*

*
* H02 --02--> T
*

*
* H10 --10--> T
*

*
* H01 --01--> T
*

*
* H00 --00--> T
*

* End of the line
*
      RETURN
      END























