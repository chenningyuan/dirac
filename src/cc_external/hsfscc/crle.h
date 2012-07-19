!
! ... nrle2 is from dimension parameters
!

      integer isize(4) 
      real*8 trle1(nrle1),trle2(nrle2),trle3(nrle3),trle4(nrle4)    

      common /crle/ trle1,trle2,trle3,trle4,isize

! to replace this line:
!--------------------------------------
!      common /crle/ trle1(nrle1),trle2(nrle2),trle3(nrle3),
!    1 trle4(nrle4),isize(4)
!

