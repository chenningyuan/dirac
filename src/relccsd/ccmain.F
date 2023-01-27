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

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PROGRAM RELCCSD
C
      use interface_to_mpi
      implicit none
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted unrestricted coupled cluster code
C    CRAY version uses dynamic memory allocation and variable file buffers
C
C---------------Routines called-------------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "ccpar.inc"
#include "infpar.h"
C
C---------------Local variables--------------------------------------
C
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
      integer i
C
C---------------Executable code--------------------------------------
C
      CALL CPUUSED(SEC0)
      IN = 5
      IW = 6
      INTFCE = 'MOLFDIR   '
      CA = .FALSE.
#if defined (VAR_MPI)
C
C     MPI initialization
C
      call interface_mpi_INIT()
      call interface_mpi_COMM_RANK(global_communicator,MYPROC)
      call interface_mpi_COMM_SIZE(global_communicator,NMPROC)
      MASTER = 0
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
      CALL CCMAIN(IN,IW,CA,INTFCE)
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
#if defined (VAR_MPI)
C
C     MPI finalization
C
      call interface_mpi_FINALIZE()
#endif
C
 1000 FORMAT (//' CPU time (seconds) used in RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used before RELCCSD:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELCCSD Run ---')
C
      END
