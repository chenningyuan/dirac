module caps
!
!*******************************************************************************
!      ----- convert strings to uppercase or lower case        -----
!      ----- adjust strings to the left and compute the length -----
!===============================================================================
!   usage:
!           call capson(string)
!           call capsof(string)
!           call strl(length,string)
!-------------------------------------------------------------------------------
!     written by: J. P. Flament
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!*******************************************************************************
!
! implicit none
   character (len=26) , parameter :: minus = 'abcdefghijklmnopqrstuvwxyz', &
                                     majus = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
   contains
   subroutine capson(carte)
!
!***********************************************************************
!     convert the string CARTE in uppercase.
!***********************************************************************
!
      implicit none
!
      character (len=*), intent(inout) :: carte
!
      integer :: ll,i,j
!
      ll = len(carte)
      do i=1,ll
         j = index(minus,carte(i:i))
         if(j.ne.0) carte(i:i) = majus(j:j)
      enddo
      return
   end subroutine capson
   subroutine capsof(carte)
!
!***********************************************************************
!     convert the string CARTE in lowercase.
!***********************************************************************
!
      implicit none
!
      character (len=*), intent(inout) :: carte
!
      integer :: ll,i,j
!
      ll = len(carte)
      do i=1,ll
         j = index(majus,carte(i:i))
         if(j.ne.0) carte(i:i) = minus(j:j)
      enddo
      return
   end subroutine capsof
   subroutine strl(long,texte)
!
!***********************************************************************
!     ----- find the effective length of the string TEXTE -----
!     ----- and adjust the texte to the left              -----
!-----------------------------------------------------------------------
!    LONG = effective length of the string
!    TEXTE = contains the string
!***********************************************************************
!
!     ----- arguments -----
!
      character (len=*), intent(inout) :: texte
      integer, intent(out) :: long
!
!     +++++
!     start
!     +++++
!
      texte = adjustl(texte)
      long = len_trim(texte)
   return
   end subroutine strl
end module caps
