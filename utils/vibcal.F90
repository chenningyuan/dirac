!dirac_copyright_start
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
!dirac_copyright_end

  program vibcal
!*******************************************************************************!
!                                                                               !
!     For a given normal mode this program will calculate the                   !
!     expectation value for a property in the five lowest                       !
!     vibrational states                                                        !
!                                                                               !
!     Written by T. Saue Feb 13 2007                                            !
!*******************************************************************************!
  use codata
  implicit none

#include "priunit.h"

  character :: potfil*30, prpfil*30, reply*1

  integer*4 :: i, j, n, iprp, ipot, ilogg, lnblnk, ier

  real*8            :: u, v2, v3, v4, p0, p1, p2, um, omega
  real*8            :: fac, a0, a1, a2, a3, v, w, har, anh, tot, dif, dif2

  real*8, parameter :: dp5=0.5d0, d1=1.0d0, dp16=0.0625d0, xtcm=2.19474625d+05

  real*8, dimension (0:4) :: b, h

  !-----------------------------------------------------------------------------!
  !                      Initialization of the code                             !
  !-----------------------------------------------------------------------------!

  call titler('VIBCAL for properties : Written by T. Saue ','*',110)

  ! read the name of potential file
  write(lupri,'(a)') 'Name input file (A30) with reduced mass (in Daltons)', &
                     'and force constants V2,V3,V4 (in au)'
  read(lucmd,'(a)') potfil

  ! read the name of property file
  write(lupri,'(a)') 'Name input file (A30) with property derivatives P0,P1,P2 (in au)'
  read(lucmd,'(a)') prpfil

  ! open potential file
  ipot = 1
  open(ipot,file=potfil,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(0,*) ' ERROR in opening file POTFIL:',potfil
    write(6,*) ' ERROR in opening file POTFIL:',potfil
    stop
  endif
  i=lnblnk(potfil)

  ! open property file
  iprp = 2
  open(iprp,file=prpfil,status='unknown',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(0,*) ' ERROR in opening file PRPFIL:',prpfil
    write(6,*) ' ERROR in opening file PRPFIL:',prpfil
    stop
  endif
  j=lnblnk(prpfil)

  ! open output file
  ilogg = 3
  open(ilogg,file=potfil(1:i)//'_'//prpfil(1:j)//'.VIB', &
       status='unknown',form='formatted',access='sequential')

  !-----------------------------------------------------------------------------!
  !                        Read data from the files                             !
  !-----------------------------------------------------------------------------!
  read(ipot,*,iostat=ier) u,v2,v3,v4
  if(ier /= 0)then
    write(0,*) ' ERROR in reading file POTFIL',potfil
    write(6,*) ' ERROR in reading file POTFIL',potfil
    stop
  endif

  write(ilogg,*) 'Reduced mass in Daltons: ',u
  write(ilogg,*) '2nd potential derivative: ', v2
  write(ilogg,*) '3rd potential derivative: ', v3
  write(ilogg,*) '4th potential derivative: ', v4

  read(iprp,*,iostat=ier) p0,p1,p2
  if(ier /= 0)then
    write(0,*) ' ERROR in reading file PRPFIL',prpfil
    write(6,*) ' ERROR in reading file PRPFIL',prpfil
    stop
  endif

  write(ilogg,*) '0th property derivative: ', p0
  write(ilogg,*) '1st property derivative: ', p1
  write(ilogg,*) '2nd property derivative: ', p2

  !-----------------------------------------------------------------------------!
  !                                Main code                                    !
  !-----------------------------------------------------------------------------!
  um = u*xfamu
  omega = sqrt(v2/um)
  fac=d1/um/omega
  a0 = p0
  a1 = dp5*p2*fac
  a2 = -dp5*p1*v3*fac*fac/omega
  a3 = -dp16*p2*v4*fac*fac*fac/omega

  write(ilogg,'(a,a)') '    n  Term 1 (P0)     Term 2 (P2)     Term 3 (P1 V3) ', &
                       ' Term 4 (P2 V4)  Harmonic cont.  Anharm. cont.   Total value'

  do n = 0,4
    v = n + dp5
    w = n*n+n+dp5
    har = a0 + a1*v
    anh = a2*v + a3*w
    tot = har + anh 
    b(n) = tot         
    h(n) = har
    write(ilogg,'(i5,7e16.8)') n,a0,a1*v,a2*v,a3*w,har,anh,tot
  enddo
  dif = b(1)-b(0)
  dif2 = dif*xthz*1000.0d0

  write(ilogg,'(a,e16.8,a,f10.1,a)') 'Difference :',dif, ' au = ',dif2, ' mHz'
  write(ilogg,'(a,f12.3,a)') 'PV shift:',2.0d0*dif2,' mHz'
  write(ilogg,'(a)') 'Harmonic model:'

  dif = h(1)-h(0)
  dif2 = dif*xthz*1000.0d0
  write(ilogg,'(a,e16.8,a,f10.1,a)') 'Difference :',DIF, ' au = ',dif2, ' mHz'
  write(ilogg,'(a,f12.3,a)') 'PV shift:',2.0d0*dif2,' mHz'

  !-----------------------------------------------------------------------------!
  end program
