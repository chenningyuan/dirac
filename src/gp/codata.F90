!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module codata
  implicit none
  public
#include "codata.h"
#include "pi.h"
! ------------------------------------------------------------
! -------------------  Index  --------------------------------
! ------------------------------------------------------------
! xtang    : a_0 (bohr radius) in angstrom
! echarge  : e in Coulomb
! hbar     : hbar in J.s
! xfmol    : Avogadro constant
! umass    : m_u in kg (unified atomic mass unit)
! pmass    : m_p in amu
! emass    : m_e in kg
! ccm      : c (speed of light) in m/s
! autk     : Hartree-Kelvin relationship E_h/k
! fermicc  : Fermi coupling constant G_F/(hbar.c)^3 in GeV^-2
! ------------------------------------------------------------

!  "The 1986 adjustment of the fundamental physical constants"
!                     E. Richard Cohen and Barry N. Taylor
!  Reviews of Modern Physics, Vol. 59, No. 4, 1987
  real(8), parameter :: xtang1986    = 0.529177249d0
  real(8), parameter :: echarge1986  = 1.60217733d-19
  real(8), parameter :: hbar1986     = 1.05457266d-34
  real(8), parameter :: xfmol1986    = 6.0221367d23
  real(8), parameter :: umass1986    = 1.6605402d-27
  real(8), parameter :: pmass1986    = 1.007276470d0
  real(8), parameter :: emass1986    = 9.1093897d-31
  real(8), parameter :: ccm1986      = 2.99792458d8
  real(8), parameter :: autk1986     = 3.157733d5
!  real(8), parameter :: fermicc1986  = No values

!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
  real(8), parameter :: xtang1998    = 0.5291772083d0
  real(8), parameter :: echarge1998  = 1.602176462d-19
  real(8), parameter :: hbar1998     = 1.054571596d-34
  real(8), parameter :: xfmol1998    = 6.02214199d23
  real(8), parameter :: umass1998    = 1.66053873d-27
  real(8), parameter :: pmass1998    = 1.007276470d0
  real(8), parameter :: emass1998    = 9.10938188d-31
  real(8), parameter :: ccm1998      = 2.997924580d8
  real(8), parameter :: autk1998     = 3.1577465d5
  real(8), parameter :: fermicc1998  = 1.16639d-5

!  "CODATA Recommended Values of the Fundamental Physical Constants: 2002"
!                     Peter J. Mohr and Barry N. Taylor
!  Reviews of Modern Physics, Vol. 77, No. 1, 2005
  real(8), parameter :: xtang2002    = 0.5291772108d0
  real(8), parameter :: echarge2002  = 1.60217653d-19
  real(8), parameter :: hbar2002     = 1.05457168d-34
  real(8), parameter :: xfmol2002    = 6.0221415d23
  real(8), parameter :: umass2002    = 1.66053886d-27
  real(8), parameter :: pmass2002    = 1.00727646688d0
  real(8), parameter :: emass2002    = 9.1093826d-31
  real(8), parameter :: ccm2002      = 2.99792458d8
  real(8), parameter :: autk2002     = 3.1577465d5
  real(8), parameter :: fermicc2002  = 1.16639d-5

!  "CODATA Recommended Values of the Fundamental Physical Constants: 2006"
!                     Peter J. Mohr, Barry N. Taylor and David B. Newell
!  Journal of Physical and Chemical Reference Data, Vol. 37, No. 3, 2008
  real(8), parameter :: xtang2006    = 0.52917720859d0
  real(8), parameter :: echarge2006  = 1.602176487d-19
  real(8), parameter :: hbar2006     = 1.054571628d-34
  real(8), parameter :: xfmol2006    = 6.02214179d23
  real(8), parameter :: umass2006    = 1.660538782d-27
  real(8), parameter :: pmass2006    = 1.00727646677d0
  real(8), parameter :: emass2006    = 9.10938215d-31
  real(8), parameter :: ccm2006      = 2.99792458d8
  real(8), parameter :: autk2006     = 3.1577465d5
  real(8), parameter :: fermicc2006  = 1.16637d-5

!  "CODATA Recommended Values of the Fundamental Physical Constants: 2010"
!                     Peter J. Mohr, Barry N. Taylor and David B. Newell
!  Journal of Physical and Chemical Reference Data, Vol. 41, No. 4, 2012
  real(8), parameter :: xtang2010    = 0.52917721092d0
  real(8), parameter :: echarge2010  = 1.602176565d-19
  real(8), parameter :: hbar2010     = 1.054571726d-34
  real(8), parameter :: xfmol2010    = 6.02214129d23
  real(8), parameter :: umass2010    = 1.660538921d-27
  real(8), parameter :: pmass2010    = 1.007276466812d0
  real(8), parameter :: emass2010    = 9.10938291d-31
  real(8), parameter :: ccm2010      = 2.99792458d8
  real(8), parameter :: autk2010     = 3.1577504d5
  real(8), parameter :: fermicc2010  = 1.166364d-5

!  "CODATA Recommended Values of the Fundamental Physical Constants: 2014"
!                     Peter J. Mohr, David B. Newell and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 45, No. 4, 2016
  real(8), parameter :: xtang2014    = 0.52917721067d0
  real(8), parameter :: echarge2014  = 1.6021766208d-19
  real(8), parameter :: hbar2014     = 1.054571800d-34
  real(8), parameter :: xfmol2014    = 6.022140857d23
  real(8), parameter :: umass2014    = 1.660539040d-27
  real(8), parameter :: pmass2014    = 1.007276466879d0
  real(8), parameter :: emass2014    = 9.10938356d-31
  real(8), parameter :: ccm2014      = 2.99792458d8
  real(8), parameter :: autk2014     = 3.1577513d5
  real(8), parameter :: fermicc2014  = 1.1663787d-5

!  "The 2018 CODATA Recommended Values of the Fundamental Physical Constants (Web Version 8.1)"
!                     Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor
!                     Database developed by J. Baker, M. Douma, and S. Kotochigova
!  Available at http://physics.nist.gov/constants, National Institute of Standards and Technology, Gaithersburg, MD 20899
  real(8), parameter :: xtang2018    = 0.529177210903d0
  real(8), parameter :: echarge2018  = 1.602176634d-19
  real(8), parameter :: hbar2018     = 1.054571817d-34
  real(8), parameter :: xfmol2018    = 6.02214076d23
  real(8), parameter :: umass2018    = 1.66053906660d-27
  real(8), parameter :: pmass2018    = 1.007276466621d0
  real(8), parameter :: emass2018    = 9.1093837015d-31
  real(8), parameter :: ccm2018      = 2.99792458d8
  real(8), parameter :: autk2018     = 3.1577502480407d5
  real(8), parameter :: fermicc2018  = 1.1663787d-5


  contains

  subroutine set_codata_values(CODAT)
  implicit none
  character*8 CODAT

!  Define which data set is used in the calculation
   if(CODAT.eq.'CODATA86') then
    xtang    = xtang1986
    echarge  = echarge1986
    hbar     = hbar1986
    xfmol    = xfmol1986
    umass    = umass1986
    pmass    = pmass1986
    emass    = emass1986
    ccm      = ccm1986
    autk     = autk1986
! Warning! As the oldest value for fermicc is from 1998,
!          we use it also for the 1986 series.
    fermicc  = fermicc1998
   elseif(CODAT.eq.'CODATA98') then
    xtang    = xtang1998
    echarge  = echarge1998
    hbar     = hbar1998
    xfmol    = xfmol1998
    umass    = umass1998
    pmass    = pmass1998
    emass    = emass1998
    ccm      = ccm1998
    autk     = autk1998
    fermicc  = fermicc1998
   elseif(CODAT.eq.'CODATA02') then
    xtang    = xtang2002
    echarge  = echarge2002
    hbar     = hbar2002
    xfmol    = xfmol2002
    umass    = umass2002
    pmass    = pmass2002
    emass    = emass2002
    ccm      = ccm2002
    autk     = autk2002
    fermicc  = fermicc2002
   elseif(CODAT.eq.'CODATA06') then
    xtang    = xtang2006
    echarge  = echarge2006
    hbar     = hbar2006
    xfmol    = xfmol2006
    umass    = umass2006
    pmass    = pmass2006
    emass    = emass2006
    ccm      = ccm2006
    autk     = autk2006
    fermicc  = fermicc2006
   elseif(CODAT.eq.'CODATA10') then
    xtang    = xtang2010
    echarge  = echarge2010
    hbar     = hbar2010
    xfmol    = xfmol2010
    umass    = umass2010
    pmass    = pmass2010
    emass    = emass2010
    ccm      = ccm2010
    autk     = autk2010
    fermicc  = fermicc2010
   elseif(CODAT.eq.'CODATA14') then
    xtang    = xtang2014
    echarge  = echarge2014
    hbar     = hbar2014
    xfmol    = xfmol2014
    umass    = umass2014
    pmass    = pmass2014
    emass    = emass2014
    ccm      = ccm2014
    autk     = autk2014
    fermicc  = fermicc2014
   else
    xtang    = xtang2018
    echarge  = echarge2018
    hbar     = hbar2018
    xfmol    = xfmol2018
    umass    = umass2018
    pmass    = pmass2018
    emass    = emass2018
    ccm      = ccm2018
    autk     = autk2018
    fermicc  = fermicc2018
   end if

! -------------------- common derived constants --------------------------
! planck   = 6.6260...d-34
  planck   = hbar*2.0d0*pi
! xtangm10 = 0.5291...d-10
  xtangm10 = xtang*1.0d-10
! cvel     = 137.0359...d0
  cvel     = ccm*xtangm10*emass/hbar
! alphac   = 0.7297...d-2
  alphac   = 1.0d0/cvel
! alpha2   = 0.5325...d-4
  alpha2   = alphac*alphac
! xtj      = 0.4848...d-17
  xtj      = hbar**2/(xtangm10*xtangm10*emass)
! xtkays   = 0.2194...d6
  xtkays   = 1.0d-2*hbar/(ccm*2.0d0*pi*xtangm10**2*emass)
! xthz     = 0.6579...d16
  xthz     = hbar/(2.0d0*pi*xtangm10*xtangm10*emass)
! xtev     = 0.3025...d2
  xtev     = xtj/echarge
! xkjmol   = 2.9195...d3
  xkjmol   = xtj*xfmol*1.0d-3
! xkcmol   = 697.7772...d0
  xkcmol   = xkjmol/4.184d0
! xajoul   = 4.848...d0
  xajoul   = 1.0d18*xtj
! autime   = 2.1752...d-17
  autime   = hbar/xtj
! xfsec    = 2.1752...d-17
  xfsec    = autime
! tesla    = 0.4254...d-5
  tesla    = (xtang*xtang*echarge/hbar)*1.0d-20
! debye    = 2.5417...d0
  debye    = echarge*xtang*ccm*1.0d11
! xtkmml   = 974.8801...d0
  xtkmml   = (xfmol**2*pi*echarge**2/3.0d0)*1.0d-7
! xfamu    = 1822.8884...d0
  xfamu    = umass/emass
! xfmp     = 1836.1526...d0
  xfmp     = xfamu*pmass
! nmagn    = 0.5050...d-26
  nmagn    = echarge*hbar/(2.0d0*xfmp*emass) ! nuclear magneton
! nmagnau  = 2.7230...d-4
  nmagnau  = 1.0d0/(2*xfmp)  ! nuclear magneton in au
! efaumksa = 5.1422...d0
  efaumksa = hbar**2/(emass*echarge*xtangm10**3)*1.0d-11
! xtnm     = 45.5788...d0
  xtnm     = 1.0d7/xtkays
! cminv    = 4.55...d-6
  cminv    = 1.0d0/xtkays     !1 centimeter-to-minus-one in au
! nm       = 18.897...d0
  nm       = 1.0d1/xtang      !1 nanometer in au
! gfermi   = 2.222...d-14
  gfermi   = (1.0d-18*hbar*ccm**3*fermicc*emass)/(echarge**2*xtangm10)

  end subroutine set_codata_values


  subroutine print_codata_reference(CODAT)
  implicit none
#include "priunit.h"
  character*8 CODAT

   if(CODAT.eq.'CODATA86') then
   write(LUPRI,*)   "  The 1986 adjustment of the fundamental physical constants              "
   write(LUPRI,*)   "               E. Richard Cohen and Barry N. Taylor                      "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 59, 1121-1148 (1987)             "

   elseif(CODAT.eq.'CODATA98') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 1998  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 72, 351-495 (2000)               "

   elseif(CODAT.eq.'CODATA02') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2002  "
   write(LUPRI,*)   "               Peter J. Mohr and Barry N. Taylor                         "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 77, 1-107 (2005)                 "

   elseif(CODAT.eq.'CODATA06') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2006  "
   write(LUPRI,*)   "            Peter J. Mohr, Barry N. Taylor and David B. Newell           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 80, 633-730 (2008)               "

   elseif(CODAT.eq.'CODATA10') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2010  "
   write(LUPRI,*)   "            Peter J. Mohr, Barry N. Taylor and David B. Newell           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 84, 1527-1605 (2012)             "

   elseif(CODAT.eq.'CODATA14') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2014  "
   write(LUPRI,*)   "            Peter J. Mohr, David B. Newell and Barry N. Taylor           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 88, 035009 (2016)                "

   elseif(CODAT.eq.'CODATA18') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2018  "
   write(LUPRI,*)   "            Peter J. Mohr, David B. Newell and Barry N. Taylor           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 93, 025010 (2021)                "

!  Default CODATA values
   elseif(CODAT.eq.'NOCODATA') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2018  "
   write(LUPRI,*)   "            Peter J. Mohr, David B. Newell and Barry N. Taylor           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 93, 025010 (2021)                "

!  If the user asks for non existing data sets, we set default values and print a WARNING
   else
   write(LUPRI,*)   " *** WARNING: CODATA set ", CODAT, " does not exist. Default is used. ***"
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2018  "
   write(LUPRI,*)   "            Peter J. Mohr, David B. Newell and Barry N. Taylor           "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 93, 025010 (2021)                "
   end if

  end subroutine print_codata_reference
end module
