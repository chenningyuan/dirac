Module Indexed_Quantities

type indexed_term
! assigned unique ID 1 - Fij, 2 - Tij, 3 - Vpqrs, 4 - Tijab ... 
integer (kind=4) :: id = 0
character(len=20) :: it_name
! ... number of indexes (2,4,6...) - always even number
integer (kind=2) :: nindx 
end type indexed_term

contains

!function introduce_indexed_term() result(this)
!type(indexed_term) :: this
!end function introduce_indexed_term

End Module Indexed_Quantities

