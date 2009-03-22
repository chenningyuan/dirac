!
! corresponds to different types of diagrams (pseudofiles)
!
! examples of file names: hhpp
!                         vphp
!                         gvhp (g-active hole), see README
!
! file names are stored in namefl(:)
!
! nstrt(1,:) - position of diagram in igen array
! nstrt(2,:) - position of diagram in gen  array
! 
! each diagram has own name, number, location on "gen" and "igen" 
! 
! nflcur - number of most recent created diagram (=total number of diagrams)
!
!
! diagram characteristics in iiia,iiva : iiia - 4 numbers (like 2233 - hhpp) iiva - 1000
!
!
!
!
!-------------------------------------------------------------------------
      integer  namefl(200),nstrt(2,0:200),nflcur,
     & iiia(200),iiva(200)

      common /cfiles/ namefl,nstrt,nflcur,iiia,iiva

!     common /cfiles/ namefl(200),nstrt(2,0:200),nflcur,
!    1 iiia(200),iiva(200)

