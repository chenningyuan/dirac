      SUBROUTINE GSTTBLD(CTT,IATP,IASM,IBTP,IBSM,
     &                  NSASO,NSBSO,PSSIGN,IDC,
     &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
*. If ISCALE = 1, the routine scales and returns the block
*  in determinant normalization, and SCLFAC = 1.0D0
*
* If ISCALE = 0, the routine does not perform any overall
* scaling, and a scale factor is returned in SCLFAC
*
* IF ISCALE = 0, zero blocks are not set explicitly to zero,
* instead  zero is returned in SCLFAC
*
* ISCALE, SCLFAC added May 97
*
* Simplified version working only for vectors on disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION CTT(*),NSASO(NSMST, *),NSBSO(NSMST, *)
      DIMENSION SCR(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN 
        write(6,*) ' GSTTBLD, IATP,IASM,IBTP,IBSM,ISCALE'
        write(6,*)            IATP,IASM,IBTP,IBSM,ISCALE     
        WRITE(6,*) ' LUC = ', LUC
      END IF
* =================
* Read in from disc
* =================
      CALL IFRMDS(LBL,1,-1,LUC)
      IF(NTEST.GE.100) write(6,*) ' LBL = ', LBL
      IF(ISCALE.EQ.1) THEN
        CALL FRMDSC(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
      ELSE
        NO_ZEROING = 1
        CALL FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
      END IF
*
      IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
        SCLFAC = 0.0D0
      ELSE
        NAST = NSASO(IASM,IATP)
        NBST = NSBSO(IBSM,IBTP)
        IF(LBL.NE.0) THEN
          CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &         IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
        ELSE
          SCLFAC = 0.0D0
        END IF
      END IF
*
C?    WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
C?   &ISCALE,SCLFAC
*
      RETURN
      END
