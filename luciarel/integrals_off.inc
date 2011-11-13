!
! VERSION : $Revision: 6375 $
! DATE    : $Date: 2007-11-23 13:33:36 +0100 (Fri, 23 Nov 2007) $
! FILE    : integrals_off.inc
!
!     =========================================
!      Offset arrays and integral array length
!      --> for large systems, we have to avoid
!      INTEGER*4 overflow...          
!      MXINT_TP and MXINT_DIA are defined in
!      mxpdim.inc
!     =========================================
!
      INTEGER*8 IBTSOSO_TP, NTERMS, NINT12, NINT_2E, 
     &          I_TOT_DIM_T, IPOSDIA
      COMMON /INT_OFFSET/IBTSOSO_TP(MXINT_TP), NTERMS,
     &                   NINT12, NINT_2E, I_TOT_DIM_T,
     &                   IPOSDIA(MXINT_DIA)
