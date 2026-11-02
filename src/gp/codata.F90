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
! ----------------------------------------------------------------------------
! ---------------------------  Index  ----------------------------------------
! ----------------------------------------------------------------------------
! xtang    : a_0 (bohr radius) in angstrom
! echarge  : e in Coulomb
! planck   : Planck constant (h) in J.s
! xfmol    : Avogadro constant
! umass    : m_u in kg (unified atomic mass unit)
! pmass    : m_p in amu
! emass    : m_e in kg
! ccm      : c (speed of light in vacuum) in m/s
! autk     : Hartree-Kelvin relationship E_h/k
! cvel     : Inverse fine-structure constant
! gfac     : Electron g-factor
! fermicc  : Fermi coupling constant G_F/(hbar.c)^3 in GeV^-2
! s2thetaw : Weak mixing angle [sen^2(theta_W)], based on
!             particular variant of modified minimal subtraction (MS) scheme
!             and commented are also the on-shell scheme values
! ----------------------------------------------------------------------------

!  The 1986 CODATA Recommended Values of the Fundamental Physical Constants
!  E. Richard Cohen and Barry N. Taylor
!  J. Phys. Chem. Ref. Data 17 (4), 1795-1803 (1988)
  real(8), parameter :: xtang1986    = 0.529177249d0
  real(8), parameter :: echarge1986  = 1.60217733d-19
  real(8), parameter :: planck1986   = 6.6260755d-34
  real(8), parameter :: xfmol1986    = 6.0221367d23
  real(8), parameter :: umass1986    = 1.6605402d-27
  real(8), parameter :: pmass1986    = 1.007276470d0
  real(8), parameter :: emass1986    = 9.1093897d-31
  real(8), parameter :: ccm1986      = 2.99792458d8
  real(8), parameter :: autk1986     = 3.157733d5
  real(8), parameter :: cvel1986     = 137.0359895d0
  real(8), parameter :: gfac1986     = 2.002319304386d0
!  real(8), parameter :: fermicc1986  = Values not reported
!  real(8), parameter :: s2thetaw1986 = Values not reported

!  Review of Particle Properties
!  L. Montanet et al. (Particle Data Group)
!  Phys. Rev. D 50, 1173 (1994); Erratum: Phys. Rev. D 51, 3975 (1995)
  real(8), parameter :: fermicc1994  = 1.16639d-5
  real(8), parameter :: s2thetaw1994 = 0.2319d0    ! on-shell scheme value: not reported

!  CODATA Recommended Values of the Fundamental Physical Constants: 1998
!  Peter J. Mohr and Barry N. Taylor
!  J. Phys. Chem. Ref. Data 28 (6), 1713-1852 (1999)
  real(8), parameter :: xtang1998    = 0.5291772083d0
  real(8), parameter :: echarge1998  = 1.602176462d-19
  real(8), parameter :: planck1998   = 6.62606876d-34
  real(8), parameter :: xfmol1998    = 6.02214199d23
  real(8), parameter :: umass1998    = 1.66053873d-27
  real(8), parameter :: pmass1998    = 1.007276470d0
  real(8), parameter :: emass1998    = 9.10938188d-31
  real(8), parameter :: ccm1998      = 2.997924580d8
  real(8), parameter :: autk1998     = 3.1577465d5
  real(8), parameter :: cvel1998     = 137.0359997d0
  real(8), parameter :: gfac1998     = 2.0023193043737d0
  real(8), parameter :: fermicc1998  = 1.16639d-5
  real(8), parameter :: s2thetaw1998 = 0.23124d0    ! on-shell scheme value: 0.2224d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2002
!  Peter J. Mohr and Barry N. Taylor
!  Rev. Mod. Phys. 77 (1), 1-107 (2005)
  real(8), parameter :: xtang2002    = 0.5291772108d0
  real(8), parameter :: echarge2002  = 1.60217653d-19
  real(8), parameter :: planck2002   = 6.6260693d-34
  real(8), parameter :: xfmol2002    = 6.0221415d23
  real(8), parameter :: umass2002    = 1.66053886d-27
  real(8), parameter :: pmass2002    = 1.00727646688d0
  real(8), parameter :: emass2002    = 9.1093826d-31
  real(8), parameter :: ccm2002      = 2.99792458d8
  real(8), parameter :: autk2002     = 3.1577465d5
  real(8), parameter :: cvel2002     = 137.03599911d0
  real(8), parameter :: gfac2002     = 2.0023193043718d0
  real(8), parameter :: fermicc2002  = 1.16639d-5
  real(8), parameter :: s2thetaw2002 = 0.23124d0    ! on-shell scheme value: 0.22215d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2006
!  Peter J. Mohr, Barry N. Taylor and David B. Newell
!  J. Phys. Chem. Ref. Data 37 (3), 1187-1284 (2008)
  real(8), parameter :: xtang2006    = 0.52917720859d0
  real(8), parameter :: echarge2006  = 1.602176487d-19
  real(8), parameter :: planck2006   = 6.62606896d-34
  real(8), parameter :: xfmol2006    = 6.02214179d23
  real(8), parameter :: umass2006    = 1.660538782d-27
  real(8), parameter :: pmass2006    = 1.00727646677d0
  real(8), parameter :: emass2006    = 9.10938215d-31
  real(8), parameter :: ccm2006      = 2.99792458d8
  real(8), parameter :: autk2006     = 3.1577465d5
  real(8), parameter :: cvel2006     = 137.035999679d0
  real(8), parameter :: gfac2006     = 2.0023193043622d0
  real(8), parameter :: fermicc2006  = 1.16637d-5
  real(8), parameter :: s2thetaw2006 = 0.23122d0    ! on-shell scheme value: 0.22255d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2010
!  Peter J. Mohr, Barry N. Taylor and David B. Newell
!  J. Phys. Chem. Ref. Data 41 (4), 043109 (2012)
  real(8), parameter :: xtang2010    = 0.52917721092d0
  real(8), parameter :: echarge2010  = 1.602176565d-19
  real(8), parameter :: planck2010   = 6.62606957d-34
  real(8), parameter :: xfmol2010    = 6.02214129d23
  real(8), parameter :: umass2010    = 1.660538921d-27
  real(8), parameter :: pmass2010    = 1.007276466812d0
  real(8), parameter :: emass2010    = 9.10938291d-31
  real(8), parameter :: ccm2010      = 2.99792458d8
  real(8), parameter :: autk2010     = 3.1577504d5
  real(8), parameter :: cvel2010     = 137.035999074d0
  real(8), parameter :: gfac2010     = 2.00231930436153d0
  real(8), parameter :: fermicc2010  = 1.166364d-5
  real(8), parameter :: s2thetaw2010 = 0.23116d0    ! on-shell scheme value: 0.2223d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2014
!  Peter J. Mohr, David B. Newell and Barry N. Taylor
!  J. Phys. Chem. Ref. Data 45 (4), 043102 (2016)
  real(8), parameter :: xtang2014    = 0.52917721067d0
  real(8), parameter :: echarge2014  = 1.6021766208d-19
  real(8), parameter :: planck2014   = 6.626070040d-34
  real(8), parameter :: xfmol2014    = 6.022140857d23
  real(8), parameter :: umass2014    = 1.660539040d-27
  real(8), parameter :: pmass2014    = 1.007276466879d0
  real(8), parameter :: emass2014    = 9.10938356d-31
  real(8), parameter :: ccm2014      = 2.99792458d8
  real(8), parameter :: autk2014     = 3.1577513d5
  real(8), parameter :: cvel2014     = 137.035999139d0
  real(8), parameter :: gfac2014     = 2.00231930436182d0
  real(8), parameter :: fermicc2014  = 1.1663787d-5
  real(8), parameter :: s2thetaw2014 = 0.23126d0    ! on-shell scheme value: 0.2223d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2018
!  Eite Tiesinga, Peter J. Mohr, David B. Newell and Barry N. Taylor
!  J. Phys. Chem. Ref. Data 50 (3), 033105 (2021)
  real(8), parameter :: xtang2018    = 0.529177210903d0
  real(8), parameter :: echarge2018  = 1.602176634d-19
  real(8), parameter :: planck2018   = 6.62607015d-34
  real(8), parameter :: xfmol2018    = 6.02214076d23
  real(8), parameter :: umass2018    = 1.66053906660d-27
  real(8), parameter :: pmass2018    = 1.007276466621d0
  real(8), parameter :: emass2018    = 9.1093837015d-31
  real(8), parameter :: ccm2018      = 2.99792458d8
  real(8), parameter :: autk2018     = 3.1577502480407d5
  real(8), parameter :: cvel2018     = 137.035999084d0
  real(8), parameter :: gfac2018     = 2.00231930436256d0
  real(8), parameter :: fermicc2018  = 1.1663787d-5
  real(8), parameter :: s2thetaw2018 = 0.23122d0    ! on-shell scheme value: 0.22290d0

!  CODATA Recommended Values of the Fundamental Physical Constants: 2022
!  Peter J. Mohr, Eite Tiesinga, David B. Newell and Barry N. Taylor
!  J. Phys. Chem. Ref. Data 54 (3), 033105 (2025)
  real(8), parameter :: xtang2022    = 0.529177210544d0
  real(8), parameter :: echarge2022  = 1.602176634d-19
  real(8), parameter :: planck2022   = 6.62607015d-34
  real(8), parameter :: xfmol2022    = 6.02214076d23
  real(8), parameter :: umass2022    = 1.66053906892d-27
  real(8), parameter :: pmass2022    = 1.0072764665789d0
  real(8), parameter :: emass2022    = 9.1093837139d-31
  real(8), parameter :: ccm2022      = 2.99792458d8
  real(8), parameter :: autk2022     = 3.1577502480398d5
  real(8), parameter :: cvel2022     = 137.035999177d0
  real(8), parameter :: gfac2022     = 2.00231930436092d0
  real(8), parameter :: fermicc2022  = 1.1663787d-5
  real(8), parameter :: s2thetaw2022 = 0.23121d0    ! on-shell scheme value: 0.22305d0


  contains

  subroutine set_codata_values(CODAT)
  implicit none
  character*8 CODAT

!  Define CODATA data set to be used in the calculation
   if(CODAT.eq.'CODATA86') then
    xtang    = xtang1986
    echarge  = echarge1986
    planck   = planck1986
    xfmol    = xfmol1986
    umass    = umass1986
    pmass    = pmass1986
    emass    = emass1986
    ccm      = ccm1986
    autk     = autk1986
    cvel     = cvel1986
    gfac     = gfac1986
! Warning! As the oldest values of fermicc and s2thetaw are from 1994,
!          we use them also for the 1986's series.
    fermicc  = fermicc1994
    s2thetaw = s2thetaw1994
   elseif(CODAT.eq.'PDG94') then
    fermicc  = fermicc1994
    s2thetaw = s2thetaw1994
   elseif(CODAT.eq.'CODATA98') then
    xtang    = xtang1998
    echarge  = echarge1998
    planck   = planck1998
    xfmol    = xfmol1998
    umass    = umass1998
    pmass    = pmass1998
    emass    = emass1998
    ccm      = ccm1998
    autk     = autk1998
    cvel     = cvel1998
    gfac     = gfac1998
    fermicc  = fermicc1998
    s2thetaw = s2thetaw1998
   elseif(CODAT.eq.'CODATA02') then
    xtang    = xtang2002
    echarge  = echarge2002
    planck   = planck2002
    xfmol    = xfmol2002
    umass    = umass2002
    pmass    = pmass2002
    emass    = emass2002
    ccm      = ccm2002
    autk     = autk2002
    cvel     = cvel2002
    gfac     = gfac2002
    fermicc  = fermicc2002
    s2thetaw = s2thetaw2002
   elseif(CODAT.eq.'CODATA06') then
    xtang    = xtang2006
    echarge  = echarge2006
    planck   = planck2006
    xfmol    = xfmol2006
    umass    = umass2006
    pmass    = pmass2006
    emass    = emass2006
    ccm      = ccm2006
    autk     = autk2006
    cvel     = cvel2006
    gfac     = gfac2006
    fermicc  = fermicc2006
    s2thetaw = s2thetaw2006
   elseif(CODAT.eq.'CODATA10') then
    xtang    = xtang2010
    echarge  = echarge2010
    planck   = planck2010
    xfmol    = xfmol2010
    umass    = umass2010
    pmass    = pmass2010
    emass    = emass2010
    ccm      = ccm2010
    autk     = autk2010
    cvel     = cvel2010
    gfac     = gfac2010
    fermicc  = fermicc2010
    s2thetaw = s2thetaw2010
   elseif(CODAT.eq.'CODATA14') then
    xtang    = xtang2014
    echarge  = echarge2014
    planck   = planck2014
    xfmol    = xfmol2014
    umass    = umass2014
    pmass    = pmass2014
    emass    = emass2014
    ccm      = ccm2014
    autk     = autk2014
    cvel     = cvel2014
    gfac     = gfac2014
    fermicc  = fermicc2014
    s2thetaw = s2thetaw2014
   elseif(CODAT.eq.'CODATA18') then
    xtang    = xtang2018
    echarge  = echarge2018
    planck   = planck2018
    xfmol    = xfmol2018
    umass    = umass2018
    pmass    = pmass2018
    emass    = emass2018
    ccm      = ccm2018
    autk     = autk2018
    cvel     = cvel2018
    gfac     = gfac2018
    fermicc  = fermicc2018
    s2thetaw = s2thetaw2018
   elseif(CODAT.eq.'CODATA22') then
    xtang    = xtang2022
    echarge  = echarge2022
    planck   = planck2022
    xfmol    = xfmol2022
    umass    = umass2022
    pmass    = pmass2022
    emass    = emass2022
    ccm      = ccm2022
    autk     = autk2022
    cvel     = cvel2022
    gfac     = gfac2022
    fermicc  = fermicc2022
    s2thetaw = s2thetaw2022
   else
!   Set default fundamental constants, from the last CODATA set
    xtang    = xtang2022
    echarge  = echarge2022
    planck   = planck2022
    xfmol    = xfmol2022
    umass    = umass2022
    pmass    = pmass2022
    emass    = emass2022
    ccm      = ccm2022
    autk     = autk2022
    cvel     = cvel2022
    gfac     = gfac2022
    fermicc  = fermicc2022
    s2thetaw = s2thetaw2022
   end if

! -------------------- common derived constants --------------------------
! alphac   = 0.7297...d-2
  alphac   = 1.0d0/cvel
! alpha2   = 0.5325...d-4
  alpha2   = alphac*alphac
! hbar     = 1.0545...d-34
  hbar     = planck/(2.0d0*pi)
! xtangm10 = 0.5291...d-10
  xtangm10 = xtang*1.0d-10
! xtj      = 0.4848...d-17
  xtj      = hbar**2/(xtangm10*xtangm10*emass)
! autime   = 2.1752...d-17
  autime   = hbar/xtj
! xtkays   = 0.2194...d6
  xtkays   = 1.0d-2*hbar/(ccm*2.0d0*pi*xtangm10**2*emass)
! cminv    = 4.55...d-6
  cminv    = 1.0d0/xtkays    ! centimeter-to-minus-one in au
! debye    = 2.5417...d0
  debye    = echarge*xtang*ccm*1.0d11
! efaumksa = 5.1422...d0
  efaumksa = hbar**2/(emass*echarge*xtangm10**3)*1.0d-11
! gfermi   = 2.222...d-14
  gfermi   = (1.0d-18*hbar*ccm**3*fermicc*emass)/(echarge**2*xtangm10)
! g1fac    = -0.9976...
  g1fac    = gfac/2.0d0-1.0d0
! nm       = 18.897...d0
  nm       = 1.0d1/xtang     ! nanometer in au
! xfamu    = 1822.8884...d0
  xfamu    = umass/emass
! xfmp     = 1836.1526...d0
  xfmp     = xfamu*pmass
! nmagn    = 0.5050...d-26
  nmagn    = echarge*hbar/(2.0d0*xfmp*emass) ! nuclear magneton
! nmagnau  = 2.7230...d-4
  nmagnau  = 1.0d0/(2*xfmp)  ! nuclear magneton in au
! tesla    = 0.4254...d-5
  tesla    = (xtang*xtang*echarge/hbar)*1.0d-20
! xajoul   = 4.848...d0
  xajoul   = 1.0d18*xtj
! xkjmol   = 2.9195...d3
  xkjmol   = xtj*xfmol*1.0d-3
! xkcmol   = 697.7772...d0
  xkcmol   = xkjmol/4.184d0
! xtev     = 0.3025...d2
  xtev     = xtj/echarge
! xthz     = 0.6579...d16
  xthz     = hbar/(2.0d0*pi*xtangm10*xtangm10*emass)
! xtkmml   = 974.8801...d0
  xtkmml   = (xfmol**2*pi*echarge**2/3.0d0)*1.0d-7
! xtnm     = 45.5788...d0
  xtnm     = 1.0d7/xtkays

  end subroutine set_codata_values


  subroutine print_codata_reference(CODAT)
!      use priunit_mod, only: LUPRI
  implicit none
#include "priunit.h"
  character*8 CODAT

   if(CODAT.eq.'CODATA86') then
   write(LUPRI,*)   "      The 1986 Adjustment of the Fundamental Physical Constants          "
   write(LUPRI,*)   "               E. Richard Cohen and Barry N. Taylor                      "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 59, 1121-1148 (1987)             "

   elseif(CODAT.eq.'PDG94') then
   write(LUPRI,*)   " Only for the weak mixing angle and the Fermi coupling constant:         "
   write(LUPRI,*)   "                  Review of Particle Properties                          "
   write(LUPRI,*)   "               L. Montanet et al. (Particle Data Group)                  "
   write(LUPRI,*)   "        Physical Reviews D, Vol. 50, 1173-1814 (1994)                    "
   write(LUPRI,*)   "        Erratum: Physical Reviews D, Vol. 51, 3975-3977 (1995)           "

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
   write(LUPRI,*)   "    Eite Tiesinga, Peter J. Mohr, David B. Newell and Barry N. Taylor    "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 93, 025010 (2021)                "

   elseif(CODAT.eq.'CODATA22') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2022  "
   write(LUPRI,*)   "    Peter J. Mohr, Eite Tiesinga, David B. Newell and Barry N. Taylor    "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 97, 025002 (2025)                "

!  Default CODATA values
   elseif(CODAT.eq.'NOCODATA') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2022  "
   write(LUPRI,*)   "    Peter J. Mohr, Eite Tiesinga, David B. Newell and Barry N. Taylor    "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 97, 025002 (2025)                "

!  If the user asks for non existing data sets, we set default values and print a WARNING
   else
   write(LUPRI,*)   " *** WARNING: CODATA set ", CODAT, " does not exist. Default is used. ***"
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2022  "
   write(LUPRI,*)   "    Peter J. Mohr, Eite Tiesinga, David B. Newell and Barry N. Taylor    "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 97, 025002 (2025)                "
   end if
!  Constants other than s2thetaw and fermicc will be taken from the default CODATA set
   if(CODAT.eq.'PDG94') then
   write(LUPRI,*)   "  CODATA Recommended Values of the Fundamental Physical Constants: 2022  "
   write(LUPRI,*)   "    Peter J. Mohr, Eite Tiesinga, David B. Newell and Barry N. Taylor    "
   write(LUPRI,*)   "        Reviews of Modern Physics, Vol. 97, 025002 (2025)                "
   end if

  end subroutine print_codata_reference
end module
