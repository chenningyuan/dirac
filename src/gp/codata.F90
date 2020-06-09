!dirac_copyright_start
!      Copyright (c) 2019 by the authors of DIRAC.
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

module codata

  implicit none
#include "pi.h"
!  revised 1999/04/29 and combined with units.h
!  Radovan: converted to f90 2009--09--23
!  Miro: towards extension by older and newer constant sets, separated with preprocessor command; 2009--11--11
!  Agustin: updated 2020/05/18

! -----------------  Index  -----------------

! xtang   : a_0 (bohr radius) in angstrom
! echarge : e in Coulomb
! hbar    : hbar in J.s
! xfmol   : Avogadro constant
! umass   : m_u in kg (unified atomic mass unit)
! pmass   : m_p in amu
! emass   : m_e in kg
! ccm     : c (speed of light) in m/s
! autk    : Hartree-Kelvin relationship E_h/k

#if   defined (CODATA86)
!  From
!  "The 1986 adjustment of the fundamental physical constants"
!                     E. Richard Cohen and Barry N. Taylor
!  Reviews of Modern Physics, Vol. 59, No. 4, 1987

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.529177249d0
  real(8), parameter, public :: echarge  = 1.60217733d-19
  real(8), parameter, public :: hbar     = 1.05457266d-34
  real(8), parameter, public :: xfmol    = 6.0221367d23
  real(8), parameter, public :: umass    = 1.6605402d-27
  real(8), parameter, public :: pmass    = 1.007276470d0
  real(8), parameter, public :: emass    = 9.1093897d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.157733d5

#elif defined (CODATA98)
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34
  real(8), parameter, public :: xfmol    = 6.02214199d23
  real(8), parameter, public :: umass    = 1.66053873d-27
  real(8), parameter, public :: pmass    = 1.007276470d0
  real(8), parameter, public :: emass    = 9.10938188d-31
  real(8), parameter, public :: ccm      = 2.997924580d8
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA02)
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 2002"
!                     Peter J. Mohr and Barry N. Taylor
!  Reviews of Modern Physics, Vol. 77, No. 1, 2005

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.5291772108d0
  real(8), parameter, public :: echarge  = 1.60217653d-19
  real(8), parameter, public :: hbar     = 1.05457168d-34
  real(8), parameter, public :: xfmol    = 6.0221415d23
  real(8), parameter, public :: umass    = 1.66053886d-27
  real(8), parameter, public :: pmass    = 1.00727646688d0
  real(8), parameter, public :: emass    = 9.1093826d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA06)
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 2006"
!                     Peter J. Mohr, Barry N. Taylor and David B. Newell
!  Journal of Physical and Chemical Reference Data, Vol. 37, No. 3, 2008

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.52917720859d0
  real(8), parameter, public :: echarge  = 1.602176487d-19
  real(8), parameter, public :: hbar     = 1.054571628d-34
  real(8), parameter, public :: xfmol    = 6.02214179d23
  real(8), parameter, public :: umass    = 1.660538782d-27
  real(8), parameter, public :: pmass    = 1.00727646677d0
  real(8), parameter, public :: emass    = 9.10938215d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA10)
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 2010"
!                     Peter J. Mohr, Barry N. Taylor and David B. Newell
!  Journal of Physical and Chemical Reference Data, Vol. 41, No. 4, 2012

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.52917721092d0
  real(8), parameter, public :: echarge  = 1.602176565d-19
  real(8), parameter, public :: hbar     = 1.054571726d-34
  real(8), parameter, public :: xfmol    = 6.02214129d23
  real(8), parameter, public :: umass    = 1.660538921d-27
  real(8), parameter, public :: pmass    = 1.007276466812d0
  real(8), parameter, public :: emass    = 9.10938291d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577504d5

#elif defined (CODATA14)
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 2014"
!                     Peter J. Mohr, David B. Newell and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 45, No. 4, 2016

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.52917721067d0
  real(8), parameter, public :: echarge  = 1.6021766208d-19
  real(8), parameter, public :: hbar     = 1.054571800d-34
  real(8), parameter, public :: xfmol    = 6.022140857d23
  real(8), parameter, public :: umass    = 1.660539040d-27
  real(8), parameter, public :: pmass    = 1.007276466879d0
  real(8), parameter, public :: emass    = 9.10938356d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577513d5

#elif defined (CODATA18)
!  From
!  "The 2018 CODATA Recommended Values of the Fundamental Physical Constants (Web Version 8.1)"
!                     Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor
!                     Database developed by J. Baker, M. Douma, and S. Kotochigova
!  Available at http://physics.nist.gov/constants, National Institute of Standards and Technology, Gaithersburg, MD 20899

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.529177210903d0
  real(8), parameter, public :: echarge  = 1.602176634d-19
  real(8), parameter, public :: hbar     = 1.054571817d-34
  real(8), parameter, public :: xfmol    = 6.02214076d23
  real(8), parameter, public :: umass    = 1.66053906660d-27
  real(8), parameter, public :: pmass    = 1.007276466621d0
  real(8), parameter, public :: emass    = 9.1093837015d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577502480407d5

#else
! Agustin: these are the default values

!  From
!  "The 2018 CODATA Recommended Values of the Fundamental Physical Constants (Web Version 8.1)"
!                     Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor
!                     Database developed by J. Baker, M. Douma, and S. Kotochigova
!  Available at http://physics.nist.gov/constants, National Institute of Standards and Technology, Gaithersburg, MD 20899

! ------------ pure numerical constants ---------------------
  real(8), parameter, public :: xtang    = 0.529177210903d0
  real(8), parameter, public :: echarge  = 1.602176634d-19
  real(8), parameter, public :: hbar     = 1.054571817d-34
  real(8), parameter, public :: xfmol    = 6.02214076d23
  real(8), parameter, public :: umass    = 1.66053906660d-27
  real(8), parameter, public :: pmass    = 1.007276466621d0
  real(8), parameter, public :: emass    = 9.1093837015d-31
  real(8), parameter, public :: ccm      = 2.99792458d8
  real(8), parameter, public :: autk     = 3.1577502480407d5

#endif
! -------------------- common derived constants --------------------------
!                               planck   = 6.6260...d-34
  real(8), parameter, public :: planck   = hbar*2*pi
!                               xtangm10 = 0.5291...d-10
  real(8), parameter, public :: xtangm10 = xtang*1.0d-10
!                               cvel     = 137.0359...d0
  real(8), parameter, public :: cvel     = ccm*xtangm10*emass/hbar
!                               alphac   = 0.7297...d-2
  real(8), parameter, public :: alphac   = 1.0d0/cvel
!                               alpha2   = 0.5325...d-4
  real(8), parameter, public :: alpha2   = alphac*alphac
!                               xtj      = 0.4848...d-17
  real(8), parameter, public :: xtj      = hbar**2/(xtangm10*xtangm10*emass)
!                               xtkays   = 0.2194...d6
  real(8), parameter, public :: xtkays   = 1.0d-2*hbar/(ccm*2.0d0*pi*xtangm10**2*emass)
!                               xthz     = 0.6579...d16
  real(8), parameter, public :: xthz     = hbar/(2.0d0*pi*xtangm10*xtangm10*emass)
!                               xtev     = 0.3025...d2
  real(8), parameter, public :: xtev     = xtj/echarge
!                               xkjmol   = 2.9195...d3
  real(8), parameter, public :: xkjmol   = xtj*xfmol*1.0d-3
!                               xtkjml   = 2.9195...d3
  real(8), parameter, public :: xtkjml   = xkjmol
!                               xkcmol   = 697.7772...d0
  real(8), parameter, public :: xkcmol   = xkjmol/4.184d0
!                               xtkcml   = 697.7772...d0
  real(8), parameter, public :: xtkcml   = xkcmol
!                               xtnm     = 45.5788...d0
  real(8), parameter, public :: xtnm     = 1.0d7/xtkays
!                               xajoul   = 4.848...d0
  real(8), parameter, public :: xajoul   = 1.0d18*xtj
!                               autime   = 2.1752...d-17
  real(8), parameter, public :: autime   = hbar/xtj
!                               xfsec    = 2.1752...d-17
  real(8), parameter, public :: xfsec    = autime
!                               tesla    = 0.4254...d-5
  real(8), parameter, public :: tesla    = (xtang*xtang*echarge/hbar)*1.0d-20
!                               debye    = 2.5417...d0
  real(8), parameter, public :: debye    = echarge*xtang*ccm*1.0d11
!                               xtkmml   = 974.8801...d0
  real(8), parameter, public :: xtkmml   = (xfmol**2*pi*echarge**2/3.0d0)*1.0d-7
!                               xfamu    = 1822.8884...d0
  real(8), parameter, public :: xfamu    = umass/emass
!                               xfmp     = 1836.1526...d0
  real(8), parameter, public :: xfmp     = xfamu*pmass
!                               nmagn    = 0.5050...d-26
  real(8), parameter, public :: nmagn    = echarge*hbar/(2*xfmp*emass) ! nuclear magneton
!                               nmagnau  = 2.7230...d-4
  real(8), parameter, public :: nmagnau  = 1.0d0/(2*xfmp) ! nuclear magneton in au
!                               efaumksa = 5.1422...d0
  real(8), parameter, public :: efaumksa = hbar**2/(emass*echarge*xtangm10**3)*1.0d-11

  contains

  subroutine print_codata_reference
#include "priunit.h"
#if   defined (CODATA86)
   write(LUPRI,*)   "  The 1986 adjustment of the fundamental physical constants              "
   write(LUPRI,*)   "               E. Richard Cohen and Barry N. Taylor                      "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 59, No. 4, 1987                  "
#elif defined (CODATA98)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 1998  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999  "
#elif defined (CODATA02)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2002  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 77, No. 1, 2005                  "
#elif defined (CODATA06)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2006  "
   write(LUPRI,*)   "            Peter J. Mohr, Barry N. Taylor and David B. Nevell           "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 37, No. 3, 2008  "
#elif defined (CODATA10)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2010  "
   write(LUPRI,*)   "            Peter J. Mohr, Barry N. Taylor and David B. Nevell           "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 41, No. 4, 2012  "
#elif defined (CODATA14)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2014  "
   write(LUPRI,*)   "            Peter J. Mohr, David B. Newell and Barry N. Taylor           "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 45, No. 4, 2016  "
#elif defined (CODATA18)
   write(LUPRI,*)   "The 2018 CODATA Recommended Values of the Fundamental Physical Constants (Web Version 8.1)"
   write(LUPRI,*)   "       Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor                 "
   write(LUPRI,*)   "       Database developed by J. Baker, M. Douma, and S. Kotochigova                       "
   write(LUPRI,*)   " Available at http://physics.nist.gov/constants                                           "
   write(LUPRI,*)   " National Institute of Standards and Technology, Gaithersburg, MD 20899                   "
#else
   write(LUPRI,*)   "The 2018 CODATA Recommended Values of the Fundamental Physical Constants (Web Version 8.1)"
   write(LUPRI,*)   "       Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor                 "
   write(LUPRI,*)   "       Database developed by J. Baker, M. Douma, and S. Kotochigova                       "
   write(LUPRI,*)   " Available at http://physics.nist.gov/constants                                           "
   write(LUPRI,*)   " National Institute of Standards and Technology, Gaithersburg, MD 20899                   "
#endif

  end subroutine print_codata_reference

  subroutine print_codata
! control print out of some important constants...
#include "priunit.h"
    call print_codata_reference
    write(LUPRI,'(2x,a,f18.12,a)') '1 debye = ',debye,' a.u,'
    write(LUPRI,'(2x,a,f18.12,a)') 'Speed of light = ',cvel,' a.u,'
  end subroutine print_codata

end module
