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

module adc_mat

! This module contains the matrix descriptor of the ADC matrix in
! the **current** run. All diagonalizers read the matrix information from
! this descriptor module.
! Additionally for the parallel run it contains local dimension descriptors

  implicit none

  save    !  if values of variables are changed they are still there after the next call 

  integer, public                  ::  reladc_md_iobase   = 0   !file handle for matrix handling
  integer, public                  ::  reladc_md_ionizl   = 0   !level od ionization for this matrix
  integer, public                  ::  reladc_md_ioldnew  = 0   !indicator for separate diag file
  integer, public                  ::  reladc_md_intbuf   = 0   !buffer length in matrix
  integer, public                  ::  reladc_md_desrep   = 0   !current symmetry
  integer, public                  ::  reladc_md_rcw      = 0   !real/complex matrix
  integer, public                  ::  reladc_md_lnzitr   = 0   !requested Krylov dimension
  integer, public                  ::  reladc_md_matdim   = 0   !matrix dimension
  integer, public                  ::  reladc_md_irecl    = 0   !record length of configuration data
  integer, public                  ::  reladc_md_nmain    = 0   !size of ADC matrix main block
  integer, public                  ::  reladc_md_nbufs    = 0   !number of buffers of len "intbuf"
  integer, public                  ::  reladc_md_neigenv  = 0   !number of requested full eigenvectors
  integer, public                  ::  reladc_md_davroots = 0   !number of requested Davidson roots
  integer, public                  ::  reladc_md_maxdavsp = 0   !maximum dimension of Davidson space
  integer, public                  ::  reladc_md_maxdavit = 0   !maximum number of Dav. iterations
  real*8,  public                  ::  reladc_md_convthr  = 0.0 !convergence of Davidson eigenvectors
  logical, public                  ::  reladc_md_davooc   = .true. !save memory in dav for large calculations
  character*6, public              ::  reladc_md_fileadc  = ''  !file name of ADC matrix
  character*6, public              ::  reladc_md_filediag = ''  !file name of diagonal of matrix, ioldnew=1
  character*6, public              ::  reladc_md_filecnf  = ''  !file name of configuration file
  character*5, public              ::  reladc_md_nmspec   = ''  !file name of spectral output file
  logical                          ::  reladc_md_isfano   = .false. !is this a fano run or not

  real, public, dimension(32,2)    ::  reladc_md_sip_eeigv= 0.0 !upper bound to eigenvalue for analysis
  real, public, dimension(32,2)    ::  reladc_md_dip_eeigv= 0.0 !upper bound to eigenvalue for analysis
  real, public                     ::  reladc_md_eeigv_lower = 0.0
  real, public                     ::  reladc_md_eeigv_upper = 0.0

end module

