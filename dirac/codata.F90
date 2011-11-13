!      Copyright (c) 2010 by the authors of DIRAC.
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

module codata

  implicit none

!  revised 1999/04/29 and combined with units.h
!
!  radovan: converted to f90 2009--09--23
!  miro: towards extension by older and newer constant sets, separated with preprocessor command; 2009--11--11
!

! TODO(miro): It would be good to have newer set of constants together with previous older sets of data.
! I'd prefer to have a keyword in DIRAC to switch between different CODATA sets. For this one would have to 'sacrifice' the parameter status
! of constants, however. Other choice is to use the preprocessor flag to compile DIRAC with chosen set of constants.

! TODO(miro):  replace constants accordingly to http://physics.nist.gov/cuu/Constants/index.html, now we have all from 1998

!---------------------------------------------------------------------------------------------------
!                               note that here we define pi here
!                               as soon as there is a pi module, use it here
!                               and remove this line
  real(8), parameter, public :: pi       = 3.14159265358979323846d0
!---------------------------------------------------------------------------------------------------

#if   defined (CODATA86)
! replace !!!!

!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant

! ------------ pure numerical constants ---------------------

  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: xfamu    = 1822.88848d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34

!                               xfmp     = m_p/m_e = proton mass in atomic units
  real(8), parameter, public :: xfmp     = 1836.152701d0

  real(8), parameter, public :: xfmol    = 6.02214199d23

!                               pmass in amu
  real(8), parameter, public :: pmass    = 1.007276470d0

!                               emass in kg
  real(8), parameter, public :: emass    = 9.10938188d-31

!                               ccm speed of light in m/s
  real(8), parameter, public :: ccm      = 2.997924580d8

  real(8), parameter, public :: xtkmml   = 974.864d0
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA98)

!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant

! ------------ pure numerical constants ---------------------

  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: xfamu    = 1822.88848d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34

!                               xfmp     = m_p/m_e = proton mass in atomic units
  real(8), parameter, public :: xfmp     = 1836.152701d0

  real(8), parameter, public :: xfmol    = 6.02214199d23

!                               pmass in amu
  real(8), parameter, public :: pmass    = 1.007276470d0

!                               emass in kg
  real(8), parameter, public :: emass    = 9.10938188d-31

!                               ccm speed of light in m/s
  real(8), parameter, public :: ccm      = 2.997924580d8

  real(8), parameter, public :: xtkmml   = 974.864d0
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA02)
! replace !!!!

!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant

! ------------ pure numerical constants ---------------------

  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: xfamu    = 1822.88848d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34

!                               xfmp     = m_p/m_e = proton mass in atomic units
  real(8), parameter, public :: xfmp     = 1836.152701d0

  real(8), parameter, public :: xfmol    = 6.02214199d23

!                               pmass in amu
  real(8), parameter, public :: pmass    = 1.007276470d0

!                               emass in kg
  real(8), parameter, public :: emass    = 9.10938188d-31

!                               ccm speed of light in m/s
  real(8), parameter, public :: ccm      = 2.997924580d8

  real(8), parameter, public :: xtkmml   = 974.864d0
  real(8), parameter, public :: autk     = 3.1577465d5

#elif defined (CODATA06)
! replace !!!!

!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant

! ------------ pure numerical constants ---------------------

  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: xfamu    = 1822.88848d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34

!                               xfmp     = m_p/m_e = proton mass in atomic units
  real(8), parameter, public :: xfmp     = 1836.152701d0

  real(8), parameter, public :: xfmol    = 6.02214199d23

!                               pmass in amu
  real(8), parameter, public :: pmass    = 1.007276470d0

!                               emass in kg
  real(8), parameter, public :: emass    = 9.10938188d-31

!                               ccm speed of light in m/s
  real(8), parameter, public :: ccm      = 2.997924580d8

  real(8), parameter, public :: xtkmml   = 974.864d0
  real(8), parameter, public :: autk     = 3.1577465d5

#else
! miro: this is the "default" for the moment, needs to be adjusted accordingly...

!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant

! ------------ pure numerical constants ---------------------

  real(8), parameter, public :: xtang    = 0.5291772083d0
  real(8), parameter, public :: xfamu    = 1822.88848d0
  real(8), parameter, public :: echarge  = 1.602176462d-19
  real(8), parameter, public :: hbar     = 1.054571596d-34

!                               xfmp     = m_p/m_e = proton mass in atomic units
  real(8), parameter, public :: xfmp     = 1836.152701d0

  real(8), parameter, public :: xfmol    = 6.02214199d23

!                               pmass in amu
  real(8), parameter, public :: pmass    = 1.007276470d0

!                               emass in kg
  real(8), parameter, public :: emass    = 9.10938188d-31

!                               ccm speed of light in m/s
  real(8), parameter, public :: ccm      = 2.997924580d8

  real(8), parameter, public :: xtkmml   = 974.864d0
  real(8), parameter, public :: autk     = 3.1577465d5

#endif

! -------------------- common derived constants --------------------------
  real(8), parameter, public :: xtangm10 = xtang*1.0d-10
  real(8), parameter, public :: cvel     = ccm*xtangm10*emass/hbar
  real(8), parameter, public :: alphac   = 1.0d0/cvel
  real(8), parameter, public :: alpha2   = alphac*alphac
  real(8), parameter, public :: xtj      = hbar**2/(xtangm10*xtangm10*emass)
  real(8), parameter, public :: xtkays   = 1.0d-2*hbar/(ccm*2.0d0*pi*xtangm10**2*emass)
  real(8), parameter, public :: xthz     = hbar/(2.0d0*pi*xtangm10*xtangm10*emass)
  real(8), parameter, public :: xtev     = xtj/echarge
  real(8), parameter, public :: xkjmol   = xtj*xfmol*1.0d-3
  real(8), parameter, public :: xkcmol   = xkjmol/4.184d0
  real(8), parameter, public :: xtkjml   = xkjmol
  real(8), parameter, public :: xtkcml   = xkcmol
  real(8), parameter, public :: xtnm     = 1.0d7/xtkays
  real(8), parameter, public :: xajoul   = 1.0d18*xtj
  real(8), parameter, public :: autime   = hbar/xtj
  real(8), parameter, public :: xfsec    = autime
  real(8), parameter, public :: tesla    = (xtang*xtang*echarge/hbar)*1.0d-20
  real(8), parameter, public :: debye    = echarge*xtang*ccm*1.0d11

  contains

  subroutine print_codata_reference
#include "priunit.h"
#if   defined (CODATA1986)
! add 
#elif defined (CODATA1998)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 1998  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999  "
#elif defined(CODATA2002)
! add
#elif defined (CODATA2006)
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2006  "
   write(LUPRI,*)   "            Peter J. Mohr, Barry N. Taylor and David B. Nevell        "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 37, No. 3, 2008  "
#else
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 1998  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999  "
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
