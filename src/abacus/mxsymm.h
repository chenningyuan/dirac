! FILE mxsymm.h

!     mxsymm.h - for symmetry detection (values set and used in hersym.F)

      REAL*8  TOLLRN, TOLLRN_2, ZERTOL
      INTEGER MAXAXS, MAXMIR,   IPRSYM

      COMMON /MXSYMM/ TOLLRN, TOLLRN_2, ZERTOL, MAXAXS, MAXMIR, IPRSYM

!     Layout of the atom arrays ATMARR(6,0:*) in hersym.F:
!       (1:3,I) coordinates, (4,I) charge CHARGE(I) (may be non-integer),
!       (5,I) original atom index, (6,I) key = IKEYZ*IZATOM(I) + ISOTOP(I).
!     Two centres can only be symmetry-related if their keys are equal
!     and their charges agree to within TOLCHG.
!     (IZATOM is 0 for point charges, so a point charge is never
!      related to a nucleus with the same charge.)
!     A warning is printed for charges that are almost equal:
!     0 < |dq| < TOLCHG (treated as equal) or TOLCHG <= |dq| < TOLCHW
!     (treated as different).
      REAL*8  TOLCHG, TOLCHW
      INTEGER IKEYZ
      PARAMETER (TOLCHG = 1.0D-6, TOLCHW = 1.0D-3, IKEYZ = 1000)

! -- end of mxsymm.h --
