!
!  ... parameters for CCSD
!
!     parameter (nbig=90000000)
!     parameter (nbig=9000000)

! mi: fix leak complain : 
!  Array reference out of bounds for array 'igen', upper bound of dimension 1 exceeded (10098089 > 9000000)
      parameter (nbig=90000000)

      parameter (ndij=960000)

!     parameter (ndint=15000000)
      parameter (ndint=150000)

      parameter (ndt1=1000)
      parameter (ndt2=500000)

!     parameter (ndp3=7000000)
      parameter (ndp3=700000)

!     parameter (ndp4=20000000)
      parameter (ndp4=2000000)

      parameter (ndp5=325000)
      parameter (ndp6=ndp3+3*ndp5)

      parameter (nrle1=8*ndt1)
      parameter (nrle2=8*ndt2)

!miro: restored for CCSDT run
      parameter (nrle3=8*ndp3)
!      parameter (nrle3=1)

      parameter (nrle4=1000)

