MODULE RECP_INPUT_READ

  IMPLICIT NONE

CONTAINS

SUBROUTINE RECP_INPUT_READ_GENERAL(ECP_INPUT,ngen,ns,naords,ncons,ngcs,ncrs)
! Read general information, Y. C. Park

  IMPLICIT NONE
! variables : global
  INTEGER :: ECP_INPUT
  INTEGER :: ngen   !
  INTEGER :: ns     !
  INTEGER :: naords !
  INTEGER :: ncons  !
  INTEGER :: ngcs   !
  INTEGER :: ncrs   !

  READ (ECP_INPUT,*) ngen, ns, naords, ncons, ngcs, ncrs

! if ( ns .gt. msu ) then
!   call bummer('change msup (two places) to ',ns,faterr)
! elseif ( naords .gt. kaordp ) then
!   call bummer('change kaordp (one place) to ',naords,faterr)
! elseif ( ncons .gt. mconsu ) then
!   call bummer('change mconsp (two places) to ',ncons,faterr)
! elseif ( ngcs .gt. mgcsu ) then
!   call bummer('change mgcsup (one place) to ',ngcs,faterr)
! endif

END SUBROUTINE RECP_INPUT_READ_GENERAL


END MODULE RECP_INPUT_READ
