! This subroutine maps the roots of effective Hamiltonian
! of the current iteration to the roots of previous iteration.
! this is done by checking how much previous model space
! coefficients overlap with the current model space coefficients.

subroutine MapRoots

   use Parameters
   use HEFFModule
   use IOModule
   use ConvMod
   implicit none
   integer :: i, j
   double precision, dimension(NSTATE,NSTATE) :: OverLap, OverLap1
   double precision, dimension(NSTATE,NSTATE) :: TMPV
   integer, dimension(NSTATE) :: MapArray, MapArray1
   integer IDAMAX

   if (CurIter==0) return
   if ((CurIter==1).and.(.not.OPTION(RESTART))) return
   return

!  call DGEMM("T","N",NSTATE,NSTATE,NSTATE,1.0d00,BCL_PREV,NSTATE, &
!              BCR,NSTATE,0.0d00,OverLap,NSTATE)
!  call DGEMM("T","N",NSTATE,NSTATE,NSTATE,1.0d00,BCR_PREV,NSTATE, &
!              BCL,NSTATE,0.0d00,OverLap1,NSTATE)

   ! Now arrange the Model space vectors using these OverLap matrices

!  do i = 1, NSTATE
!     MapArray(i)  = IDAMAX(NSTATE,OverLap(1,i),1)
!     MapArray1(i) = IDAMAX(NSTATE,OverLap1(1,i),1)
!     write(UOUT,*) i,MapArray(i), MapArray1(i)
!  enddo

!  if (.not.CompareArrays(NSTATE,MapArray,MapArray1)) then
!     write(UOUT,*) "PROBLEM WITH OVERLAP MATRICES"
!     write(UOUT,*) "BOTH OVERLAP MATRICES PREDICT DIFFERENT ROOT MAPS"
!     write(UOUT,*) "UNABLE TO DETERMINE CORRECT MAPPING"
!     call flush(UOUT)
!     call abort
!  endif

   ! Print out there is a remap.

!  do i = 1, NSTATE
!     MapArray1(i) = i
!  enddo
!
!  if (.not.CompareArrays(NSTATE,MapArray,MapArray1)) then
!     write(UOUT,*) "There has been a remap at Iteration",CurIter
!  endif
!
!  ! Map Energy and Coefficients using the Mapping Array
!
!  do i = 1, NSTATE
!     TMPV(i,1) = EN(MapArray(i))
!  enddo
!  EN = TMPV(:,1)
!
!  do i = 1, NSTATE
!     TMPV(:,i) = BCL(:,MapArray(i))
!  enddo
!  BCL = TMPV
!
!  do i = 1, NSTATE
!     TMPV(:,i) = BCR(:,MapArray(i))
!  enddo
!  BCR = TMPV

   return

contains

! Returns true if A1 and A2 are same arrays

function CompareArrays(NSZ,A1,A2)

   implicit none
   integer, intent(in) :: NSZ
   integer :: A1(NSZ), A2(NSZ)
   logical :: CompareArrays
   integer :: i

   CompareArrays = .TRUE.
   do i = 1, NSZ
      if (A1(i)/=A2(i)) then
         CompareArrays = .FALSE.
         exit
      endif
   enddo

   return

end function CompareArrays

end subroutine MapRoots
