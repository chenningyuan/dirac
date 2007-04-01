!
! common block for handling CAP integrals
! to replace tedious reading from INPCAP,INPQD files
!
! PW 'theta'
!
! IPRCAP - print level
!
! ... maximum number of intervals for integration
!  also check the size MAXSPL when reading !

      PARAMETER (MAXSPL=100)

! ... flag for comparison between DIRAC-CAP overlap integrals
      LOGICAL CAPOVLP
      COMMON / DCCAPOVL / CAPOVLP

      REAL*8 PW,  X0(1:3),LAMDA(1:3),GINTPT(3,MAXSPL)
      COMMON / DCCAPR / PW, X0, LAMDA, GINTPT

      INTEGER*4  NGR(1:3), NGRS(3,MAXSPL), IPRCAP 
      COMMON / DCCAPI / NGR, NGRS, IPRCAP

