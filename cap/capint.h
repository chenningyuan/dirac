!
! common block CAP integration parameters
! replace tedious reading from INPCAP,INPQD files
!
! PW 'theta'
!
! ... maximum number of intervals for integration
!  also check the size MAXSPL when reading !
!
      PARAMETER (MAXSPL=100)

      REAL*8 PW, X0(1:3), LAMDA(1:3), GINTPT(1:3,MAXSPL)
      COMMON / DCCAPR / PW, X0, LAMDA, GINTPT

      INTEGER*4  NGR(1:3)
      COMMON / DCCAP_INT1 / NGR

! must be splitted 
      INTEGER*4  NGRS(1:3,MAXSPL)
      COMMON / DCCAP_INT2 / NGRS

