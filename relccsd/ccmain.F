C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccmain.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PROGRAM RELCCSD
C
      IMPLICIT INTEGER (A-Z)
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
      PARAMETER (MAXCOR = 13 000 000)
#if !defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) 
      COMMON // CC(MAXCOR)
#endif
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC
      REAL*8 SEC
C
C---------------Executable code--------------------------------------
C
C     Define primary in and output units.
C     Other files units are defined in files include block.
C
      IN = 5
      IW = 6
C
      CALL CCHEADER
      MXCOR = MAXCOR
      CALL CCMAIN(CC,MXCOR)
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC
      CALL QUIT("NORMAL END OF PROGRAM")
 1000 FORMAT (//' Total amount of CPU time (seconds) used :',F14.4,
     &//'  --- Normal end of Run ---')
      END
