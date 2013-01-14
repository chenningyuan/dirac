!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end
!
module x2c_fragment_x2c_interface

! stefan: this module contains all the functionality to carry out an ensuing 
!         transformation
!         ========================
!         h_2c^++ = U^+ * h_4c * U
!         ========================
!         where U was obtained by means of an atom-by-atom X2C approach:
!         U =~ SUM_{all unique atomic center/fragments i} U_i = U^{@}
!
!         for more details see Peng and Reiher, JCP 2012
!
! written by sknecht jan 2013
!

  use fragment_x2c_cfg

  implicit none

  public x2c_fragment_x2c_driver

  private


contains

!***********************************************************************
  subroutine x2c_fragment_x2c_driver(nz,             &
                                     nfsym,          &
                                     naosh_all_cent)

!**********************************************************************
!
!    purpose: read pct-matrices U_fragment for atomic types and
!             generate the approximate molecular U
!
!----------------------------------------------------------------------
  integer, intent(in)  :: nz
  integer, intent(in)  :: nfsym
  integer, intent(in)  :: naosh_all_cent(nfsym)
!----------------------------------------------------------------------
  real(8), allocatable :: umol(:)
!----------------------------------------------------------------------


!     safety check
      if(.not.fragment_x2c_info%fragment_approach_enabled) return

!     default for the following procedure is to use the atomic fragments
!     ----------------------------------------------------------------------------
!     possible extensions are: unit matrices for light atoms and fragment matrices 
!     for heavy centers being close or even sharing a chemical bond.
!     ----------------------------------------------------------------------------

!     call construct_umol_from_atoms

  end subroutine x2c_fragment_x2c_driver
!***********************************************************************

  subroutine construct_umol_from_atoms(nz,             &
                                       nfsym,          &
                                       naosh_all_cent)

!**********************************************************************
!
!    purpose: read pct-matrices U_fragment for atomic types and
!             generate the approximate molecular U
!
!----------------------------------------------------------------------
  real(8), intent(inout) :: umol
  integer, intent(in)    :: nz
  integer, intent(in)    :: nfsym
  integer, intent(in)    :: naosh_all_cent(nfsym)
!----------------------------------------------------------------------
  integer, allocatable   :: rbuf(:)
  integer, allocatable   :: cbuf(:)
  integer, allocatable   :: ibuf(:)
  integer, allocatable   :: jbuf(:)
  integer                :: icent
  integer                :: iatom
!----------------------------------------------------------------------
      icent    = 1
      do iatom = 1, ncont

!       get # of basis functions (L and L+S) for this atom type
!       read atomic pctmat from file

        call put_aoblock(umol,
                         naosh_all_cent(1),
                         naosh_L(1),
                         icent,
                         NONT(IATOM),
                         NUCDEG(ICENT),
                         NZC1,
                         WORK(KADE),
                         NBAS,
                         NBAS,
                         -1,
                         NBAS,
                         NBAS,
                         -1,
                         rbuf,
                         cbuf,
                         ibuf,
                         jbuf)
        ICENT = ICENT+NONT(IATOM)
      end do

  end subroutine construct_umol_from_atoms
!***********************************************************************
end module x2c_fragment_x2c_interface

