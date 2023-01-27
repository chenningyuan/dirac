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

module polprp_cfg
     use iso_c_binding

! This module contains all input for the POLPRP/DAVIDSON module

  implicit none

  save    !  if values of variables are changed their values are retained

!
!  polarization propagator control variables
!
  integer, dimension(1:32), public :: polprp_statesym     = 0  
  logical, public                  :: polprp_doextended   = .false.
  real*8,  public                  :: polprp_writethr     = 0.0
  logical, public                  :: polprp_dodiag       = .true.
  logical, public                  :: polprp_dotrmo       = .false.
  integer, public                  :: polprp_printlev     = 0
  logical, public                  :: polprp_skipccseti   = .false.
!
!  davidson control variables
!
  integer, public                  :: polprp_davroots   = 5
  integer, public                  :: polprp_davmaxsp   = 200
  integer, public                  :: polprp_davmaxit   = 5
  real*8,  public                  :: polprp_davconv    = 1.0E-05
  logical, public                  :: polprp_davreort   = .false.
  
end module
