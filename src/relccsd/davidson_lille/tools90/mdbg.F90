module mdebug
!
!***********************************************************************
!    ----- this module contains the routines for debugging -----
!-----------------------------------------------------------------------
!    * subroutine subini initializes debugging parameters
!
!    * subroutine subent(subrou) and subret(subrou) control
!      the entry and exit of subroutine subrou
!
!    usage:
!    ======
!
!    program prog
!    implicit none
!    .....
!    .....
! !  first executable tatement
!    get somehow the value of SUB_FIRST
!    call subini(sub_first)
!    .....
!    .....
! !  end program prog
!    subroutine subxxx
!    implicit none
!    .....
!    .....
! !  first executable statement
!    call subent('subxxx')
!    .....
!    .....
!    call subret('subxxx',1)
!    return
!    .....
!    !....
!    call subret('subxxx',2)
!    return
!    end subroutine subxxx
!
!    to print a text
!    ----------------
!    call dbgprt(text)
!
!    to print a value
!    ----------------
!    call dbgprt(text,val)
!
!    to activate debugging at the beginning of spinor
!    ---------------------
!
!    spinor -debug <data >result
!
!    to activate debugging starting at subroutine toto
!    ---------------------
!
!    spinor -debug toto <data >result
!=======================================================================
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   use kindvar
   use caps, only: capson
   implicit none
!
!  ----- common variables -----
!
   character (len=60), save :: dmark=repeat(' ',60)
   character (len=10), save :: sdbug
   character (len=10), save, dimension (1000) :: caller
   integer, save :: lmark,n_call,lsubr,true_call
   logical, save :: ldbug,first
   private :: first,n_call,lsubr,true_call,dmark,lmark
   integer, private :: iw=6
!
!     ----- generic module for debug printing -----
!
   interface dbgprt
     module procedure dprt_a,dprt_sa,dprt_si4,dprt_sr8,dprt_ti4,dprt_tr8,  &
                      dprt_sl4,dprt_si2,dprt_ti2,dprt_sc8,dprt_tc8,dprt_tl4, &
                      dprt_tr4
   end interface dbgprt
!
   contains
   subroutine subini(debugs)
!
!***********************************************************************
!     ----- initialises the trace of subroutine calls -----
!-----------------------------------------------------------------------
!=======================================================================
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: debugs
!
!     +++++
!     start
!     +++++
!
      first = .TRUE.
      sdbug = debugs
      sdbug = adjustl(sdbug)
      lsubr = int(len_trim(sdbug))
      call capson(sdbug)
      ldbug = lsubr.eq.0
      lmark = 3
      n_call = 1
      true_call = 1
      caller(n_call) = '?   '
      return
   end subroutine subini
   subroutine subent(subrou)
!
!***********************************************************************
!  ----- this subroutine records the entrance in subroutine SUBROU -----
!  ----- it should be called on entry of subroutine  SUBROU        -----
!-----------------------------------------------------------------------
!   SUBROU(C??) = name of the subroutine (INPUT)
!=======================================================================
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: subrou
!
!     ----- local variable -----
!
      character (len=10) :: cdum
!
!     ----- formats -----
!
 9999 format(a,'*** ENTER ',a,'  CALLED BY ',a)
!
!     +++++
!     start
!     +++++
!
      cdum = subrou
      cdum = adjustl(cdum)
      call capson(cdum)
      if(first) then
         if(cdum.eq.sdbug) then
            ldbug = .TRUE.
            first = .FALSE.
         endif
      endif
      if(ldbug) then
         write(6,9999) dmark(1:lmark),subrou,caller(n_call)
         call flshfi(iw)
         if(cdum.ne.caller(n_call)) then
            lmark = lmark+2
            n_call = n_call+1
            caller(n_call) = cdum
         endif
         true_call = true_call+1
      endif
      return
   end subroutine subent
   subroutine subret(subrou,i)
!
!***********************************************************************
!  ----- this subroutine records the exit of subroutine SUBROU -----
!  ----- it should be called on exit of subroutine  SUBROU     -----
!-----------------------------------------------------------------------
!   SUBROU(C??) = name of the subroutine (INPUT)
!         I(I4) = an integer to identify multiple return points (INPUT)
!=======================================================================
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: subrou
      integer (kind=K_I4), intent(in) :: i
!
!     ----- local variable -----
!
      character (len=10) :: cdum
!
!     +++++
!     start
!     +++++
!
      if(ldbug) then
         cdum = subrou
         cdum = adjustl(cdum)
         if(n_call.eq.true_call) then
            lmark = lmark-2
            n_call = n_call-1
!           write(6,"(3x,a,'*** EXIT  ',a,i4)") dmark(1:lmark),subrou,i
!           call flshfi(iw)
         endif
         true_call = true_call-1
         write(6,"(a,'*** EXIT  ',a,i4)") dmark(1:lmark),subrou,i
         if(cdum.eq.sdbug) then
            ldbug = .FALSE.
            first = .TRUE.
         endif
         call flshfi(iw)
      endif
      return
    end subroutine subret
    subroutine dprt_a(text)
!
!***********************************************************************
!  ----- this subroutine prints a string -----
!-----------------------------------------------------------------------
!   text(C??) = text to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      call flshfi(iw)
      return
   end subroutine dprt_a
   subroutine dprt_sa(text,aval)
!
!***********************************************************************
!  ----- this subroutine prints a string and a character variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    AVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text,aval
!
!     +++++
!     start
!     +++++
!
      write(6,"(3a)") dmark(1:lmark+2),text,aval
      call flshfi(iw)
      return
   end subroutine dprt_sa
   subroutine dprt_sc8(text,ival)
!
!***********************************************************************
!  ----- this subroutine prints a string and an L4 variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      complex (kind=K_C8), intent(in) :: ival
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,'(',1pd23.15,',',d23.15,')')") dmark(1:lmark+2),text,ival
      call flshfi(iw)
      return
   end subroutine dprt_sc8
   subroutine dprt_si4(text,ival)
!
!***********************************************************************
!  ----- this subroutine prints a string and an I4 variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer (kind=K_I4), intent(in) :: ival
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,i10)") dmark(1:lmark+2),text,ival
      call flshfi(iw)
      return
   end subroutine dprt_si4
   subroutine dprt_si2(text,ival)
!
!***********************************************************************
!  ----- this subroutine prints a string and an I4 variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer (kind=K_I2), intent(in) :: ival
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,i10)") dmark(1:lmark+2),text,ival
      call flshfi(iw)
      return
   end subroutine dprt_si2
   subroutine dprt_sr8(text,r8val)
!
!***********************************************************************
!  ----- this subroutine prints a string and a R8 variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      double precision, intent(in) :: r8val
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,1pd23.15)") dmark(1:lmark+2),text,r8val
      call flshfi(iw)
      return
   end subroutine dprt_sr8
    subroutine dprt_ti4(text,ival,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and an I4 array of N elements -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      integer (kind=K_I4), intent(in), dimension (:) :: ival
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,' : ',i5)") dmark(1:lmark+2),i,ival(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_ti4
   subroutine dprt_ti2(text,ival,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and an I2 array -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      integer (kind=K_I2), intent(in), dimension (:) :: ival
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,i10)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,' : ',i5)") dmark(1:lmark+2),i,ival(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_ti2
   subroutine dprt_tr8(text,r8val,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and a R8 array of N elements -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      double precision, intent(in), dimension (:) :: r8val
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,1pd23.15)") dmark(1:lmark+2),i,r8val(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_tr8
   subroutine dprt_tr4(text,r4val,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and a R4 array of N elements -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      real, intent(in), dimension (:) :: r4val
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,1pd23.15)") dmark(1:lmark+2),i,r4val(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_tr4
   subroutine dprt_tc8(text,c8val,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and a C8 array of N elements -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    IVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      double complex, intent(in), dimension (:) :: c8val
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,'(',1pd23.15,',',d23.15,')')") dmark(1:lmark+2),i,c8val(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_tc8
   subroutine dprt_sl4(text,lval)
!
!***********************************************************************
!  ----- this subroutine prints a string and an L4 variable -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!    LVAL(I4) = value to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      logical (kind=K_L4), intent(in) :: lval
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a,l4)") dmark(1:lmark+2),text,lval
      call flshfi(iw)
      return
   end subroutine dprt_sl4
   subroutine dprt_tl4(text,l4val,n)
!
!***********************************************************************
!  ----- this subroutine prints a string and an L4 array of N elements -----
!-----------------------------------------------------------------------
!   TEXT(C??) = text to be printed (INPUT)
!   L4VAL(I4) = array to be printed (INPUT)
!=======================================================================
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
      implicit none
!
!     ----- arguments -----
!
      character (len=*), intent(in) :: text
      integer, intent(in) :: n
      logical, intent(in), dimension (:) :: l4val
!
!     ----- local variables -----
!
      integer :: i
!
!     +++++
!     start
!     +++++
!
      write(6,"(2a)") dmark(1:lmark+2),text
      do i=1,n
         write(6,"(a,i4,l4)") dmark(1:lmark+2),i,l4val(i)
      enddo
      call flshfi(iw)
      return
   end subroutine dprt_tl4
   subroutine flshfi(nft)
!
!***********************************************************************
!     ----- to flush the print buffer -----
!***********************************************************************
!
   use kindvar
   implicit none
!
!     ----- arguments -----
!
   integer, intent(in) :: nft
!
!     ----- local variables -----
!
#ifdef  _G95_64_
   integer (kind=K_I8) :: nft8
#else
   integer :: nft8
#endif
!
!     +++++
!     start
!     +++++
!
   nft8 = nft
#ifdef _UNDER_
   call flush_(nft8)
#else
   call flush(nft8)
#endif
   return
   end subroutine flshfi
end module mdebug
