!
! common block CAP integration parameters
! read from the DIRAC input file are 
! replacing tedious reading from INPCAP,INPQD files
!
! PW       =  'theta'  x,y,z
! X0       =  'zeta'   x,y,z
! LAMDA    =  'lambda' x,y,z
!
! ... maximum number of intervals for integration
!     check the value MAXSPL when reading !
!
      PARAMETER (MAXSPL=30)

      REAL*8 PW, X0(1:3), LAMDA(1:3), GINTPT(1:3,MAXSPL)
      COMMON / DCCAPR / PW, X0, LAMDA, GINTPT

!
! ... DC field as part of the CAP operator
!
      REAL*8  X0DC(1:3), EPSDC(1:3)
      COMMON / CAP_DC / X0DC, EPSDC

      INTEGER*4  NGR(1:3)
      COMMON / DCCAP_INT1 / NGR

! ...  must be splitted 
      INTEGER*4  NGRS(1:3,MAXSPL)
      COMMON / DCCAP_INT2 / NGRS
