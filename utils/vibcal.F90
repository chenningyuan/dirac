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
  use num_der
  implicit none

#include "priunit.h"

  character :: potfil*30, prpfil*30, reply*1, strkey*30

  integer*4 :: i, j, n, iprp, ipot, ilogg, lnblnk, ier
  integer*4 :: keymeth, npoints, ipoint

  integer*4, parameter :: nder=4, nprint=0

  real*8            :: u, v2, v3, v4, p0, p1, p2, um, omega, step
  real*8            :: fac, a0, a1, a2, a3, v, w, har, anh, tot, dif, dif2

  real*8, parameter :: dp5=0.5d0, d1=1.0d0, dp16=0.0625d0, xtcm=2.19474625d+05
  real*8, parameter :: engerr=1.0d-10

  real*8, dimension ( 0:4 ) :: b, h
  real*8, dimension(nder,3) :: derres

  real*8, allocatable, dimension(:) :: xcoor, ycoor

  !-----------------------------------------------------------------------------!
  !                      Initialization of the code                             !
  !-----------------------------------------------------------------------------!

  call titler('VIBCAL for properties : Written by T. Saue ','*',110)

  ! decide how to get derivatives
  write(lupri,'(a)') 'How to get force constants and property derivatives [default:1]:'
  write(lupri,'(a)') '      1. Read them from input file '
  write(lupri,'(a)') '      2. Calculate them numerically'
  read(lucmd,'(a)',iostat=ier) strkey; read(strkey,*,iostat=ier) keymeth
  if(keymeth /= 2) keymeth = 1

  if(keymeth == 1)then
    ! read the name of potential input file
    write(lupri,'(a)') 'Name input file (A30) with reduced mass (in Daltons)', &
                       'and force constants V2,V3,V4 (in au)'
    read(lucmd,'(a)') potfil
    ! read the name of property input file
    write(lupri,'(a)') 'Name input file (A30) with property derivatives P0,P1,P2 (in au)'
    read(lucmd,'(a)') prpfil
  elseif(keymeth == 2)then
    ! read reduced mass
    write(lupri,'(a)') 'Reduced mass (in Daltons)'
    read(lucmd,*) u
    ! read the name of potential file
    write(lupri,'(a)') 'Name of the input file (A30) with potential curve (in au.)'
    read(lucmd,'(a)') potfil
    ! read the name of property file
    write(lupri,'(a)') 'Name of the input file (A30) with property curve (in au.)'
    read(lucmd,'(a)') prpfil
  endif

  ! open potential input file
  ipot = 1
  open(ipot,file=potfil,status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(0,*) ' ERROR in opening file POTFIL:',potfil
    write(6,*) ' ERROR in opening file POTFIL:',potfil
    stop
  endif
  i=lnblnk(potfil)

  ! open property input file
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
  open(ilogg,file=potfil(1:i)//'_'//prpfil(1:j)//'.vib', &
       status='unknown',form='formatted',access='sequential')

  !-----------------------------------------------------------------------------!
  !                 Get force constants and property derivatives                !
  !-----------------------------------------------------------------------------!
  if(keymeth == 1)then
    read(ipot,*,iostat=ier) u,v2,v3,v4
    if(ier /= 0)then
      write(0,*) ' ERROR in reading file POTFIL',potfil
      write(ilogg,*) ' ERROR in reading file POTFIL',potfil
      stop
    endif

    read(iprp,*,iostat=ier) p0,p1,p2
    if(ier /= 0)then
      write(0,*) ' ERROR in reading file PRPFIL',prpfil
      write(ilogg,*) ' ERROR in reading file PRPFIL',prpfil
      stop
    endif
  endif

  if(keymeth == 2)then
    ! count number of points
    rewind ipot
    npoints = 0
    do
      read(ipot,*,iostat=ier)
      if(ier == -1)exit
      npoints = npoints + 1
    enddo
    if(nprint==1)write(ilogg,'(a,i4,a)') ' We have ',npoints,' points'

    ! allocate arrays
    allocate ( xcoor(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with xcoor  all.'
    allocate ( ycoor(npoints),stat=ier ); if(ier /= 0 ) stop 'Problem with ycoor  all.'

    !-----------------------!
    !  Get force constants  !
    !-----------------------!
    ! read points
    rewind ipot
    do i=1,npoints
      read(ipot,*,iostat=ier) xcoor(i), ycoor(i)
      if(ier /= 0)then
        write(ilogg,*) ' ERROR during reading POTFIL:',potfil
        write(ilogg,*) ' STOP in vibcal'
        write(0,*) ' STOP in vibcal (see output)'
        stop
      endif
    enddo

    ! the derivation will be calculated in this point, with this step
    ipoint = (npoints / 2) + mod(npoints,2)
    step = xcoor(2) - xcoor(1)
    if(nprint==1)then
      write(ilogg,'(a,i4,a)') ' We use  ',ipoint,'th point for calculating of derivatives'
      write(ilogg,'(a,es12.4)') ' step = ',step
      write(ilogg,*)
    endif

    ! calculate derivatives and estimate errors
    call nddrv(nder,'CENTRAL',ipoint,nprint,ilogg,npoints,xcoor,ycoor,engerr,derres)

    ! store results
    v2 = derres(2,1); v3 = derres(3,1); v4 = derres(4,1)

    !----------------------------!
    !  Get property derivatives  !
    !----------------------------!
    ! read points
    rewind iprp
    do i=1,npoints
      read(iprp,*,iostat=ier) xcoor(i), ycoor(i)
      if(ier /= 0)then
        write(ilogg,*) ' ERROR during reading PRPFIL:',prpfil
        write(ilogg,*) ' STOP in vibcal'
        write(0,*) ' STOP in vibcal (see output)'
        stop
      endif
    enddo

    ! the derivation will be calculated in this point
    ipoint = (npoints / 2) + mod(npoints,2)
    step = xcoor(2) - xcoor(1)
    if(nprint==1)then
      write(ilogg,'(a,i4,a)') ' We use  ',ipoint,'th point for calculating of derivatives'
      write(ilogg,'(a,es12.4)') ' step = ',step
      write(ilogg,*)
    endif

    ! calculate derivatives and estimate errors
    call nddrv(nder,'CENTRAL',ipoint,nprint,ilogg,npoints,xcoor,ycoor,engerr,derres)

    ! store results
    p0 = ycoor(ipoint); p1 = derres(1,1); p2 = derres(2,1)
  endif

  write(ilogg,*) 'Reduced mass in Daltons: ',u
  write(ilogg,*) '2nd potential derivative: ', v2
  write(ilogg,*) '3rd potential derivative: ', v3
  write(ilogg,*) '4th potential derivative: ', v4

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
