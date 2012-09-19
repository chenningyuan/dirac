!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!
!
! this module contains all the functionality required to add one- and two-electron
! spin-orbit corrections to the hamiltonian in DIRAC-sorted SA-AO basis. 
! the spin-orbit corrections will be obtained as the spin-dependent parts of the 
! (converged) 4c-Fock operator in an atomic Dirac run.
!
! routines written by tsaue march 2011, collected by sknecht sep 2012
!
module atom2molecule_utils

  implicit none

  public labcount
  public ipmat
  public revindex

  private

contains

!***********************************************************************

      SUBROUTINE IPMAT(NR,NC,A,LDA,B,LDB,INDR,INDC)
!-----------------------------------------------------------------------
!     This routine will transfer elements from matrix B to matrix A
!     according to row and column pointers INDR and INDC
!
!	Written by T.Saue, March 2011
!
!-----------------------------------------------------------------------
  real(8), intent(out) :: a(lda,*)
  real(8), intent(in)  :: b(ldb,nc)
  integer, intent(in)  :: indr(nr)
  integer, intent(in)  :: indc(nc)
  integer, intent(in)  :: nc
  integer, intent(in)  :: nr
  integer, intent(in)  :: lda
  integer, intent(in)  :: ldb
!-----------------------------------------------------------------------
  integer              :: i, j, ii, jj
!-----------------------------------------------------------------------

      DO J = 1,NC
        JJ = INDC(J)
        DO I = 1,NR
          II = INDR(I)
          A(II,JJ)=B(I,J)
        END DO
      END DO

      END
!***********************************************************************

      SUBROUTINE REVINDEX(ILIST,NI,JLIST,NJ)
!-----------------------------------------------------------------------
!     Make reverse indexing for elements of ILIST.GE.0
!     Written by T. Saue March 31 2011
!-----------------------------------------------------------------------
  integer, intent(out) :: jlist(nj)
  integer, intent(in)  :: ilist(ni)
  integer, intent(in)  :: ni
  integer, intent(in)  :: nj
!-----------------------------------------------------------------------
  integer              :: i, j
!-----------------------------------------------------------------------

      J = 0
      DO I = 1, NI
        IF(ILIST(I) > 0)THEN
          J = J + 1
          IF(J > NJ) CALL QUIT('REVINDEX: Out of bounds !')
          JLIST(J) = I
        ENDIF
      END DO

      END
!***********************************************************************

      SUBROUTINE LABCOUNT(N,LIST,NBAS,ILAB,I1,I2,I3,I4)
!***********************************************************************
!
!     Do statistics or make pointers using basis function labels
!
!     ILAB = 1: A)AO-labels
!     ---------------------------------------------------------------------
!     I1   ICLS  - Basis set class: large(1), small(2) ..
!     I2   ICENT - index of center in list of symmetry-independent centers
!     I3   IDEG  - for given symmetry-independent center, index of degeneracy
!     I4   ITYP  - basis function type, e.g. s, px, dxy etc.
!
!     ILAB = 2: SO-labels
!     ---------------------------------------------------------------------
!     I1   IC    - Basis set class: large(1), small(2) ....
!     I2   ICENT - index of center in list of symmetry-independent centers
!     I3   IRP   - index of boson irrep
!     I4   ITYP  - basis function type, e.g. s, px, dxy
!
!     T. Saue Mar 31 2011
!
!***********************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "dcblab.h"
#include "dcbibt.h"
      DIMENSION LIST(NBAS)

!.....Initialize index array
      CALL ICOPY(NBAS,1,0,LIST,1)
!.....If active, search first value
      IF(I1.GT.0) THEN
        N=0
        DO I = 1,NBAS
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = IGET(ITEMP)
          IF(IVAL.NE.I1) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDDO
      ENDIF
!.....If active, search second value
      IF(I2.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = JGET(ITEMP)
          IF(IVAL.NE.I2) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
!.....If active, search third value
      IF(I3.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = KGET(ITEMP)
          IF(IVAL.NE.I3) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
!.....If active, search third value
      IF(I4.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = LGET(ITEMP)
          IF(IVAL.NE.I4) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
 
      END
!***********************************************************************

end module
