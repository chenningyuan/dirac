MODULE Vectors

!use Memory
! deta type for strings
use mystring

!   printout parameter
integer(kind=4), private :: iprint_vect=0

type Vect
!private ! I had to release he private statement due to easy access to pointers !!!
! contains position on the PW (WORK) array 
  integer (kind=4) :: kstart 
!   the total size on the array
 ! integer (kind=8) :: array_size = 0
  integer :: array_size = 0
! ...important flag to check if vector is allocated !
  logical ::  is_alloc = .false.
! may be real or integer 
  logical ::  is_real = .false.
  logical ::  is_logi = .false.

! short text attached to the vector
  type(string), dimension(1) :: vect_text

! ... pointers - very effective tools 
  real(kind=8), dimension(:), pointer :: p_rA1D => NULL()
  integer,      dimension(:), pointer :: p_iA1D => NULL()
  logical,      dimension(:), pointer :: p_lA1D => NULL()

end type Vect

type Array2Dim
! ... 2-dimensional array based upon single vector
  integer(kind=4) :: dimens1=0, dimens2=0 ! 2nd dimension
  real(kind=8), dimension(:,:), pointer :: p_rA2D => NULL()
  integer,      dimension(:,:), pointer :: p_iA2D => NULL()
  logical,      dimension(:,:), pointer :: p_lA2D => NULL()
  logical  :: is_alloc = .false.
  type(Vect) :: Vect2D
end type Array2Dim

type Array3Dim
! ... 3-dimensional array based upon single vector
  integer(kind=4) :: dimens1=0, dimens2=0, dimens3=0 ! dimensions
  real(kind=8), dimension(:,:,:), pointer :: p_rA3D => NULL()
  integer,      dimension(:,:,:), pointer :: p_iA3D => NULL()
  logical,      dimension(:,:,:), pointer :: p_lA3D => NULL()
  logical  :: is_alloc = .false.
  type(Vect) :: Vect3D
end type Array3Dim

! yet TODO: program all 3Dim stuff 

! ... common routine names (methods overloading)
interface Fill2D_Vect
  module procedure Fill2D_Vect_re
  module procedure Fill2D_Vect_int
end interface Fill2D_Vect

interface ZeroObject
  module procedure ZeroVector
  module procedure Zero2DArray
end interface ZeroObject 

interface PrintObject
  module procedure print_vector 
  module procedure print_2Darray
end interface PrintObject 

interface DeAllocObject
  module procedure DeAllocVect
  module procedure DeAlloc2DArray
end interface DeAllocObject

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_isalloc_vect(self)
!-----------------------------------------------------------------------
!
! set the is_alloc attribute to .TRUE.
!
!-----------------------------------------------------------------------
type(Vect), intent(inout) :: self
call QENTER('set_isalloc_vect')
 self%is_alloc = .true.
call QEXIT('set_isalloc_vect')
end subroutine set_isalloc_vect

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine unset_isalloc_vect(self)
!-----------------------------------------------------------------------
!
!      set the is_alloc attribute to .FALSE.
!
!-----------------------------------------------------------------------
type(Vect), intent(inout) :: self
call QENTER('unset_isalloc_vect')
 self%is_alloc = .false.
call QEXIT('unset_isalloc_vect')
end subroutine unset_isalloc_vect

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_vect_kstart(self) result (kstart_out)
!----------------------------------------------------------------------
!
!    returns kstart value of entering vector data type "self"
!
!----------------------------------------------------------------------
#include "priunit.h"
type(Vect), intent(in) :: self
integer(kind=4) :: kstart_out
call QENTER('get_vect_kstart')
  if (.not.self%is_alloc) then
    call print_vect_attrib(self)
    call QUIT('get_vect_kstart: can not return "kstart" position'  &
            //' when memory for vector is NOT allocated !')
  else
    kstart_out = self%kstart
  endif
call QEXIT('get_vect_kstart')
end function get_vect_kstart

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_vect_iprint(iprint)
!==============================================================
!
! set up the print variable in the module: iprint_vect
!
!==============================================================
#include "priunit.h"
#include "iratdef.h"
integer (kind=8/IRAT), intent(in) :: iprint
CALL QENTER('set_vect_iprint')
 !write(LUPRI,*) 'set_vect_iprint: iprint=',iprint
 iprint_vect = iprint
 if (iprint_vect >= 5 ) then
  write(LUPRI,'(2X,A,I3)') 'set_vect_iprint: control printout variable is ',iprint_vect
 endif
CALL QEXIT('set_vect_iprint')
end subroutine set_vect_iprint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_vect_iprint() result(iprint)
!-------------------------------------------------------------------
! returns the printing parameter of the module
!-------------------------------------------------------------------
integer(kind=4) :: iprint
call QENTER('get_vect_iprint')
 iprint = iprint_vect
call QEXIT('get_vect_iprint')
end function get_vect_iprint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_vect_attrib(self)
!-------------------------------------------------------------------
!
!   print out only attributes of the vector data type
!
!-------------------------------------------------------------------
#include "priunit.h"
type(Vect), intent(in) :: self
logical :: is_pointer
call QENTER('print_vect_attrib')

 write(LUPRI,'(/2x,a)')  &
 '---------- Basic attributes of the "vector" data type -----------'

 write(LUPRI,'(2x,a,i3)') 'print level ',iprint_vect

 write(LUPRI,'(a,l1)')   &
 'Is "Vect" data type variable allocated (is_alloc) ?',self%is_alloc

 write(LUPRI,'(a,l1)') 'Is "Vect" data type variable real ?',self%is_real
 write(LUPRI,'(a,l1)') 'Is "Vect" data type variable logical ?',self%is_real

 if (self%is_alloc) then
  write(LUPRI,'(a,i7)') 'size of array:',self%array_size
 endif

 if (is_string(self%vect_text)) then
  write(LUPRI,'(4x,a)') '...attached text to the "Vect" type variable:'
  call print_string(self%vect_text)
 endif

 is_pointer = .false.

 if (associated(self%p_rA1D)) then
  write(LUPRI,'(a)') 'Pointer of "Vect" type variable for "real*8" 1D array.'
  is_pointer = .true.
 endif  

 if (associated(self%p_iA1D)) then
  write(LUPRI,'(a)') 'Pointer of "Vect" type variable for "integer" 1D array.'
  is_pointer = .true.
 endif  

 if (associated(self%p_lA1D)) then
  write(LUPRI,'(a)') 'Pointer of "Vect" type variable for "logical" 1D array.'
  is_pointer = .true.
 endif  

 if (self%is_alloc.and..not.is_pointer) then
   write(LUPRI,'(a)')   & 
   'print_vect_attrib: Vector is allocated but NO pointer is associated ! Wrong !'
   call QUIT('print_vect_attrib: Vector is allocated but NOT pointer is associated !')
 endif

call QEXIT('print_vect_attrib')
end subroutine print_vect_attrib

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_2Darray_attrib(self)
!------------------------------------------------------------------------
!
!   print out only attributes of the entering vector data type "self"
!
!------------------------------------------------------------------------
#include "priunit.h"
type(Array2Dim), intent(in) :: self
call QENTER('print_2Darray_attrib')

 write(LUPRI,'(/2x,a)')  &
 '---------- Basic attributes of the derived "Array2Dim" data type -----------'

 call print_vect_attrib(self%Vect2D)

 write(LUPRI,'(a,i3)') 'module print level:',iprint_vect

 write(LUPRI,'(a,l1)') 'Is "Array2Dim" allocated (is_alloc) ?',self%is_alloc
 write(LUPRI,'(a,2i4)') '"Array2Dim" dimensions :',self%dimens1,self%dimens2
 CALL FLSHFO(LUPRI)

 if (associated(self%p_rA2D)) then
   write(LUPRI,'(a)') 'Pointer of "Array2Dim" type variable for "real*8" 2D array.'
 else if (associated(self%p_iA2D)) then
   write(LUPRI,'(a)') 'Pointer of "Array2Dim"  type variable for "integer" 2D array.'
 else if (associated(self%p_lA2D)) then
   write(LUPRI,'(a)') 'Pointer of "Array2Dim" type variable for "logical" 2D array.'
 else 
   write(LUPRI,'(a)') 'No pointer of "Array2Dim" ("real*8"/"integer"/"logical") associated.'
   if (self%is_alloc) then
     write(LUPRI,'(a)')  &
     'Error ! No pointer of "Array2Dim" associated while self%is_alloc=.true.!'
     call QUIT('print_2Darray_attrib: Error ! No pointer of "Array2Dim" associated' & 
     //' while self%is_alloc=.true.!')
   endif
 endif  

 CALL FLSHFO(LUPRI)

call QEXIT('print_2Darray_attrib')
end subroutine print_2Darray_attrib

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_3Darray_attrib(self)
!-------------------------------------------------------------------
!
!   print out only attributes of the data type
!
!-------------------------------------------------------------------
#include "priunit.h"
type(Array3Dim), intent(in) :: self
call QENTER('print_3Darray_attrib')

 write(LUPRI,'(/2x,a)')  &
 '---------- Basic attributes of the derived "Array3Dim" data type -----------'
 call print_vect_attrib(self%Vect3D)

 write(LUPRI,'(a,l1)') 'Is "Array3Dim" allocated (is_alloc) ?',self%is_alloc
 write(LUPRI,'(a,3i4)') '"Array3Dim" dimensions :', &
              self%dimens1,self%dimens2,self%dimens3

 if (associated(self%p_rA3D)) then
   write(LUPRI,'(a)') 'Pointer of "Array3Dim" type variable for "real*8" 3D array.'
 else if (associated(self%p_iA3D)) then
   write(LUPRI,'(a)') 'Pointer of "Array3Dim"  type variable for "integer" 3D array.'
 else if (associated(self%p_lA3D)) then
   write(LUPRI,'(a)') 'Pointer of "Array3Dim" type variable for "logical" 3D array.'
 else 
   write(LUPRI,'(a)') & 
   'No pointer of "Array3Dim" ("real*8"/"integer"/"logical") associated.' &
  //' Not yet allocated.'
   if (self%is_alloc) then
     write(LUPRI,'(a)')  &
     'Error ! No pointer of "Array3Dim" associated while self%is_alloc=.true.!'
     call QUIT('print_3Darray_attrib: Error ! No pointer of "Array3Dim" associated' & 
     //' while self%is_alloc=.true.!')
   endif
 endif  

 CALL FLSHFO(LUPRI)

call QEXIT('print_3Darray_attrib')
end subroutine print_3Darray_attrib


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_vector(self)
!-------------------------------------------------------------------
!
!  print vector attributes and its elements as well
!
!-------------------------------------------------------------------
#include "priunit.h"
type(Vect), intent(in) :: self
call QENTER('print_vector')
call print_vect_attrib(self)
if (self%is_real.and.associated(self%p_rA1D)) then
    call print_real_vector(self%array_size,self%p_rA1D)
else if (.not.self%is_real.and.associated(self%p_iA1D)) then
  ! call print_int_vector(self%array_size,PW(self%kstart))
   call print_int_vector(self%array_size,self%p_iA1D)
else if (self%is_logi.and.associated(self%p_lA1D)) then
   call quit('not yet implemented !')
else
   write(LUPRI,*) 'print_vector: error branch !  '  
   write(LUPRI,*) 'self%is_real:',self%is_real
   write(LUPRI,*) 'self%is_logi:',self%is_logi
   write(LUPRI,*) 'associated(self%p_rA1D):',associated(self%p_rA1D)
   write(LUPRI,*) 'associated(self%p_iA1D):',associated(self%p_iA1D)
   write(LUPRI,*) 'associated(self%p_lA1D):',associated(self%p_lA1D)
   CALL FLSHFO(LUPRI)
   call quit('print_vector: error branch !')
endif
call QEXIT('print_vector')
end subroutine print_vector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_2Darray(self)
!---------------------------------------------------------------------
!
!       print out 2D array plus info on object
!
!---------------------------------------------------------------------
#include "priunit.h"
type(Array2Dim), intent(in) :: self
integer :: i,j
integer :: k
!integer(kind=4) :: i,j
call QENTER('print_2Darray')
call print_2Darray_attrib(self)
!write(LUPRI,*) 'self%dimens1,dimens2:',self%dimens1,self%dimens2
k = 0 
do i=1,self%dimens1
do j=1,self%dimens2
 k = k + 1
 if (associated(self%p_rA2D)) then
   write(LUPRI,'(1x,i3,1x,i3,A,D14.7,A,$)') i,j,'->',self%p_rA2D(i,j),'|'
   if (mod(k,4)==0) then
     write(LUPRI,*) 
     CALL FLSHFO(LUPRI)
     k = 0
   endif
 else if (associated(self%p_iA2D)) then
   write(LUPRI,'(1x,i3,1x,i3,A,I5,A,$)') i,j,'->',self%p_iA2D(i,j),'|'
   if (mod(k,5)==0) then
     write(LUPRI,*) 
     CALL FLSHFO(LUPRI)
     k = 0
   endif
 else if (associated(self%p_lA2D)) then
   call QUIT('this branch is not active...')
 else
   call QUIT('print_2Darray: no pointer associated !')
 endif
enddo
enddo

call QEXIT('print_2Darray')
end subroutine print_2Darray

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRINT_REAL_VECTOR(N,VEC)
!-----------------------------------------------------------------------
!
! Print out real vector
!
! Written by Miro Ilias, April 2007, Tel Aviv
!  Last change: MI, Octorer 2007, TA
!
!-----------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
      DIMENSION VEC(N)
      CALL QENTER('PRINT_REAL_VECTOR')

      DO I = 1, N
        WRITE(LUPRI,'(1X,I4,A,D14.6,A,$)') I,'->',VEC(I),'|'
        IF (MOD(I,5).EQ.0) THEN
          WRITE(LUPRI,*)
          CALL FLSHFO(LUPRI)
        ENDIF
      ENDDO

      CALL QEXIT('PRINT_REAL_VECTOR')
      RETURN
      END SUBROUTINE PRINT_REAL_VECTOR

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRINT_INT_VECTOR(N,IVEC)
!-----------------------------------------------------------------------
!
!  Print out the integer vector. 
!
! Written by Miro Ilias, April 2007, Tel Aviv
!  Last change: MI, Octorer 2007, TA
!
!-----------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
      DIMENSION IVEC(N)
      CALL QENTER('PRINT_INT_VECTOR')
      DO I = 1, N
        WRITE(LUPRI,'(1X,I3,A,I8,A,$)') I,'->',IVEC(I),'|'
        IF (MOD(I,5).EQ.0) THEN
           WRITE(LUPRI,*)
           CALL FLSHFO(LUPRI)
        ENDIF
      ENDDO
      CALL QEXIT('PRINT_INT_VECTOR')
      RETURN
      END  SUBROUTINE PRINT_INT_VECTOR

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function MakeVector(char_type,isize_in,text) result(self)
!--------------------------------------------------------------------------
!
!    basic constructor for s simple vector: real,integer or logical
!
!  On entry:
!===============
!  char_type  'REAL,'INTE','LOGI' for real, integer and logical variables
!
!  isize      length of the vector
!
!  text     string to be attached to the vector (optional)
!
!---------------------------------------------------------------------------
use Memory
use GeneralPointers
#include "iratdef.h"
#include "priunit.h"
character(*), intent(in) :: char_type
!integer(kind=8/IRAT),intent(in) :: isize
integer(kind=8/IRAT) :: isize
!integer(kind=8),intent(in) :: isize
integer,intent(in) :: isize_in
character(*), optional, intent(in) :: text
type(Vect) :: self
logical :: is_pointer
call QENTER('MakeVector')

isize = isize_in

! set string to zero by default
call set_no_string(self%vect_text)

self%is_real  = .false.;  self%is_logi  = .false.

! ... sometimes older compilers do not initialize variables by declaration
self%p_rA1D => NULL(); self%p_iA1D => NULL(); self%p_lA1D => NULL()

! add the text if is present
if (present(text)) then
   call make_string(self%vect_text, text)
endif

if (iprint_vect >= 10 ) then
 write(lupri,'(2X,A,I4)')  &
 'MakeVector: goingh to allocate vector in memory of lenght ',isize
 write(lupri,'(A,A)') 'MakeVector: entering descriptive text is :',text
 write(lupri,'(A,A)') 'MakeVector: entering data type :', char_type
 CALL FLSHFO(LUPRI)
endif

! ... allocate memmory in common DIRAC space
if (char_type(1:4) .EQ. 'REAL') then
  CALL mymemget('REAL',self%kstart,isize)
  if (iprint_vect >= 10 ) then
   write(LUPRI,*) 'MakeVector: mymemget(REAL...done'
   CALL FLSHFO(LUPRI)
  endif
  self%is_real  = .true.
  self%is_logi  = .false.
  call  AssignVect1D_real(PW(self%kstart), isize, iprint_vect )
  self%p_rA1D => p_real8; p_real8 => NULL()
else if (char_type(1:4) .EQ. 'INTE') then
  CALL mymemget('INTE',self%kstart,isize)
  if (iprint_vect >= 10 ) then
   write(LUPRI,*) 'MakeVector: mymemget(INTE...done'
   CALL FLSHFO(LUPRI)
  endif
  self%is_real  = .false.
  self%is_logi  = .false.
  call  AssignVect1D_inte(PW(self%kstart), isize, iprint_vect )
  self%p_iA1D => p_inte; p_inte => NULL()
else if (char_type(1:4) .EQ. 'LOGI') then
  self%is_real  = .false.
  self%is_logi  = .true.
  CALL mymemget('LOGI',self%kstart,isize)
  call quit('MakeVector - this branch is NOT yet ready!')
else
  write(LUPRI,*) '=== Output from MakeVector ==='
  write(LUPRI,*) 'Entering vector TYPE: ',char_type(1:4)
  call QUIT('MakeVector - neither REAL/INTE/LOGI ! Dead branch, check it ! ')
endif

self%array_size = isize
self%is_alloc = .true. ! array is allocated

Call ZeroVector(self) ! automatically make it zero

 is_pointer = .false.

 if (associated(self%p_rA1D)) then
  if (iprint_vect >= 10 ) then
    write(LUPRI,'(a)') 'MakeVector: pointer of "Vect" type variable is for "real*8" 1D array.'
  endif
  is_pointer = .true.
 endif  

 if (associated(self%p_iA1D)) then
  if (iprint_vect >= 10 ) then
    write(LUPRI,'(a)') 'MakeVector: pointer of "Vect" type variable is for "integer" 1D array.'
  endif
  is_pointer = .true.
 endif  

 if (associated(self%p_lA1D)) then
  if (iprint_vect >= 10 ) then
    write(LUPRI,'(a)') 'MakeVector: pointer of "Vect" type variable is for "logical" 1D array.'
    call Quit('MakeVector: logical vector is not yet programmed !')
  endif
  is_pointer = .true.
 endif  

 if (self%is_alloc.and..not.is_pointer) then
   write(LUPRI,'(a)')   & 
   'MakeVector: Vector is allocated but NO pointer is associated ! Wrong !'
   call QUIT('MakeVector: Vector is allocated but NO pointer is associated !')
 endif

if (iprint_vect >= 10 ) then
 write(lupri,'(2X,A,I4)')  &
 'MakeVector: vector was allocated in memory of lenght ',self%array_size
 CALL FLSHFO(LUPRI)
endif

call QEXIT('MakeVector')
end function MakeVector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function Make2DimArray(char_type,isize1,isize2,text) result(self)
!--------------------------------------------------------------------------
!
!    basic constructor for 2-dimensional array: real,integer or logical
!
!  On entry:
!===============
!
!  char_type  'REAL,'INTE','LOGI' for real, integer and logical variables
!
!  isize1, isize2    array dimensions
!
!  text      string to be attached to the 2DimArray2DimArray (optional)
!
!--------------------------------------------------------------------------
use Memory
#include "iratdef.h"
#include "priunit.h"
character(*), intent(in) :: char_type
!integer(kind=8/IRAT),intent(in) :: isize1,isize2
integer,intent(in) :: isize1,isize2
integer(kind=8/IRAT) :: isize
!integer(kind=8),intent(in) :: isize
character(*), optional, intent(in) :: text
type(Array2Dim) :: self
call QENTER('Make2DimArray')

isize = isize1*isize2

if (iprint_vect >= 10 ) then
 write(lupri,'(A,A)') 'Make2DimArray: entering text is :',text
 write(lupri,'(A,A)') 'Make2DimArray: entering character type :', char_type
 CALL FLSHFO(LUPRI)
endif

self%Vect2D = MakeVector(char_type,isize,text)
!  ...  assign array dimensions
self%dimens1 = isize1; self%dimens2 = isize2
self%is_alloc = .true.

if (iprint_vect >= 10 ) then
 write(lupri,'(A,3I4)') 'Make2DimArray: 2D array was allocated in memory'  &
                   //' with dimensions ',self%dimens1,self%dimens2
 CALL FLSHFO(LUPRI)
endif

! assign the pointer
! ... sometimes older compilers do not initialize variables by declaration
self%p_rA2D => NULL(); self%p_iA2D => NULL(); self%p_lA2D => NULL()
call Assign2DimArrayPointer(self)

call QEXIT('Make2DimArray')
end function Make2DimArray


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function Make3DimArray(char_type,isize1,isize2,isize3,text) result(self)
!--------------------------------------------------------------------------
!
!   basic constructor for 3-dimensional array: real,integer or logical
!
!  On entry:
!===============
!  char_type  'REAL,'INTE','LOGI' for real, integer and logical variables
!
!  isize1, isize2, isize3    array dimensions
!
!  text      string to be attached to the Array3Dim (optional)
!
!--------------------------------------------------------------------------
use Memory
#include "iratdef.h"
#include "priunit.h"
character(*), intent(in) :: char_type
!integer(kind=8/IRAT),intent(in) :: isize1,isize2
integer,intent(in) :: isize1,isize2,isize3
integer(kind=8/IRAT) :: isize
!integer(kind=8),intent(in) :: isize
character(*), optional, intent(in) :: text
type(Array3Dim) :: self
call QENTER('Make3DimArray')

isize = isize1*isize2*isize3

if (iprint_vect >= 10 ) then
 call HEADER('Output from Make3DimArray',-1)
 write(lupri,'(2X,A)') 'entering text is :',text
 write(lupri,'(2X,A)') 'entering character type :', char_type
 write(lupri,'(2X,A,3I4)') 'entering dimensions of 3D array :',isize1,isize2,isize3
 write(lupri,'(2X,A,3I4)') 'total size of 3D array (n1 x n2 x n3) :',isize
 CALL FLSHFO(LUPRI)
endif

self%Vect3D = MakeVector(char_type,isize,text)
!  ...  assign array dimensions
self%dimens1 = isize1; self%dimens2 = isize2; self%dimens3 = isize3
self%is_alloc = .true.

if (iprint_vect >= 10 ) then
 write(lupri,'(A,3I4)') 'Make3DimArray: 3D array was allocated in memory'  &
 //' (using MakeVector) with dimensions ',self%dimens1,self%dimens2,self%dimens3
 CALL FLSHFO(LUPRI)
endif

! assign the pointer
! ... sometimes older compilers do not initialize variables by declaration
self%p_rA3D => NULL(); self%p_iA3D => NULL(); self%p_lA3D => NULL()
call Assign3DimArrayPointer(self)

call QENTER('Make3DimArray')
end function Make3DimArray

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Assign2DimArrayPointer(self)
!--------------------------------------------------------------------------
!
!    assign pointer to the 2D array AFTER the memory allocation
!
!--------------------------------------------------------------------------
use Memory
use GeneralPointers
#include "priunit.h"
type(Array2Dim),intent(inout) :: self
call QENTER('Assign2DimArrayPointer')

if (.not.self%Vect2D%is_alloc) then
  call print_2Darray_attrib(self) ! ... control printout
  call QUIT('Assign2DimArrayPointer: vector NOT allocated !')
endif

if (self%Vect2D%is_real .AND. .not.self%Vect2D%is_logi) then
    call AssignVect2D_real(PW(self%Vect2D%kstart), iprint_vect, self%dimens1, self%dimens2 )
    self%p_rA2D => p2_real8;  p2_real8=> NULL()
    if (.not.associated(self%p_rA2D)) then
      call QUIT('Assign2DimArrayPointer: self%p_rA2D  not associated !')
    endif
else if (.not.self%Vect2D%is_real .AND. .not.self%Vect2D%is_logi) then
    call AssignVect2D_inte(PW(self%Vect2D%kstart), iprint_vect, self%dimens1, self%dimens2 )
    self%p_iA2D => p2_inte;  p2_inte => NULL()
    if (.not.associated(self%p_iA2D)) then
      call QUIT('Assign2DimArrayPointer: self%p_iA2D  not associated !')
    endif
else if (.not.self%Vect2D%is_real .AND. self%Vect2D%is_logi) then
    call quit('not ready yet')
else
    call QUIT('Assign2DimArrayPointer: vector NOT real, neither integer, nor logical !')
endif

call QEXIT('Assign2DimArrayPointer')
end subroutine Assign2DimArrayPointer

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Assign3DimArrayPointer(self)
!--------------------------------------------------------------------------
!
!    assign pointer to the 3D array AFTER the memory allocation
!
!--------------------------------------------------------------------------
use Memory
use GeneralPointers
#include "priunit.h"
type(Array3Dim),intent(inout) :: self
call QENTER('Assign3DimArrayPointer')

if (.not.self%Vect3D%is_alloc) then
  call print_3Darray_attrib(self) ! ... control printout
  call QUIT('Assign3DimArrayPointer: vector is NOT allocated !')
endif

if (self%Vect3D%is_real .AND. .not.self%Vect3D%is_logi) then
    call AssignVect3D_real(PW(self%Vect3D%kstart), iprint_vect, & 
                       self%dimens1, self%dimens2, self%dimens3 )
    self%p_rA3D => p3_real8;  p3_real8=> NULL()
    if (.not.associated(self%p_rA3D)) then
      call QUIT('Assign3DimArrayPointer: self%p_rA3D  not associated !')
    endif
else if (.not.self%Vect3D%is_real .AND. .not.self%Vect3D%is_logi) then
    call AssignVect3D_inte(PW(self%Vect3D%kstart), iprint_vect, &
                       self%dimens1, self%dimens2, self%dimens3 )
    self%p_iA3D => p3_inte;  p3_inte => NULL()
    if (.not.associated(self%p_iA3D)) then
      call QUIT('Assign3DimArrayPointer: self%p_iA3D not associated !')
    endif
else if (.not.self%Vect3D%is_real .AND. self%Vect3D%is_logi) then
    call quit('Assign3DimArrayPointer: branch not yet ready ...')
else
    call QUIT('Assign3DimArrayPointer: vector NOT real, neither integer, nor logical !')
endif

call QEXIT('Assign3DimArrayPointer')
end subroutine Assign3DimArrayPointer

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Fill2D_Vect_re(self,what)
!-----------------------------------------------------------------------
!
!   Fill real 2d array with specific number 'what'
!
!-----------------------------------------------------------------------
#include "priunit.h"
type(Array2Dim),intent(inout) :: self
real(kind=8), intent(in) :: what
integer i,j
call QENTER('Fill2D_Vect_re')

if (.not.self%is_alloc.or..not.associated(self%p_rA2D)) then
  call print_2Darray_attrib(self)
  call QUIT('Fill2D_Vect_re: 2D array is NOT allocated or pointer is not associated !')
endif

if (iprint_vect >= 10) then
  write(LUPRI,*) 'Fill2D_Vect_re: going to fill array  !'
  CALL FLSHFO(LUPRI)
endif

do i = 1, self%dimens1
do j = 1, self%dimens2
  self%p_rA2D(i,j) = what  ! store the number, using appropriate pointer  
  if (iprint_vect >= 30) then
    write(LUPRI,'(2X,A,I3,A,I3,A,D12.6)') 'self%p_rA2D(',i,';',j,')=',self%p_rA2D(i,j)
  endif
enddo
enddo

if (iprint_vect >= 10) then
 write(LUPRI,'(A,2I4,A,D12.6)') &
 'Fill2D_Vect_re: 2D array ',self%dimens1,self%dimens2,' filled with ',what
 CALL FLSHFO(LUPRI)
endif

call QEXIT('Fill2D_Vect_re')
end subroutine Fill2D_Vect_re

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Fill2D_Vect_int(self,what)
!-----------------------------------------------------------------------
!
!   Fill integer 2d array with specific number 'what'
!
!-----------------------------------------------------------------------
#include "priunit.h"
type(Array2Dim),intent(inout) :: self
integer, intent(in) :: what
integer i,j
call QENTER('Fill2D_Vect_int')

if (.not.self%is_alloc.or..not.associated(self%p_iA2D)) then
  call print_2Darray_attrib(self)
  call QUIT('Fill2D_Vect_int: 2D array is NOT allocated or pointer is not associated !')
endif

do i = 1, self%dimens1
do j = 1, self%dimens2
  self%p_iA2D(i,j) = what  ! store the number, using appropriate pointer  
  if (iprint_vect >= 30) then
    write(LUPRI,'(2X,A,I3,A,I3,A,I5)') 'self%p_rA2D(',i,';',j,')=',self%p_iA2D(i,j)
  endif
enddo
enddo

if (iprint_vect >= 10) then
 write(LUPRI,'(A,2I4,A,I6)') &
 'Fill2D_Vect_int: 2D array ',self%dimens1,self%dimens2,' filled with ',what
endif

call QEXIT('Fill2D_Vect_int')
end subroutine Fill2D_Vect_int


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine ZeroVector(self)
!-----------------------------------------------------------------------
!
!         Fill vector with zeroes !
!
! For logical vector set them to .FALSE.
!
!-----------------------------------------------------------------------
!use Memory
#include "priunit.h"
type(Vect), intent(inout) :: self
call QENTER('ZeroVector')
if (self%is_alloc) then
  if (self%is_real.and..not.self%is_logi) then
    ! call DZERO( PW(self%kstart), self%array_size )
     call DZERO( self%p_rA1D, self%array_size )
  else if (.not.self%is_real.and..not.self%is_logi) then
    ! call IZERO( PW(self%kstart), self%array_size )
     call IZERO( self%p_iA1D, self%array_size )
  else if (.not.self%is_real.and.self%is_logi) then
     !call LZERO( PW(self%kstart), self%array_size )
     call quit('ZeroVector: this branch is NOT implemented !')
     call LZERO( self%p_lA1D, self%array_size )
  else
     call print_vect_attrib(self)
     call QUIT('ZeroVector: dead branch !')
  endif
  if (iprint_vect >= 10 ) then
    write(lupri,'(2X,A)') 'ZeroVector: vector was filled with zeros '
  endif
endif
call QEXIT('ZeroVector')
end subroutine ZeroVector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine Zero2DArray(self)
!-----------------------------------------------------------------------
! 
! Fill 2D array with zeroes
!
!-----------------------------------------------------------------------
#include "priunit.h"
type(Array2Dim), intent(inout) :: self
call QENTER('Zero2DArray')
 call ZeroVector(self%Vect2D)
call QEXIT('Zero2DArray')
end subroutine Zero2DArray

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine DeAllocVect(self)
!-----------------------------------------------------------------------
!
!        deallocate the entering vector "self"
!
!-----------------------------------------------------------------------
use Memory
#include "priunit.h"
type(Vect),intent(inout) :: self
call QENTER('DeAllocVect')

if (iprint_vect >= 10) then
  write(LUPRI,'(/,2X,A)') 'DeAllocVect: going to deallocate the vector'
  call PrintObject(self)
  write(LUPRI,'(/2X,A)') 'DeAllocVect: related memory parameters to this vector:'
  call print_mem_param
endif

if (self%is_alloc) then
 call mymemrel('Vectors: DeAllocVect',self%kstart,self%kstart)
 self%is_alloc = .false.
 self%array_size = 0
 call cancel_string(self%vect_text) ! cancel the string text
else
 call QUIT('DeAllocVect: was NOT allocated !')
endif

if (iprint_vect >= 10) then
  write(LUPRI,'(/,2X,A)') 'DeAllocVect: vector was deallocated'
  call PrintObject(self)
  write(LUPRI,'(/2X,A)') 'DeAllocVect: related memory parameters of this released vector:'
  call print_mem_param
endif

call QEXIT('DeAllocVect')
end subroutine DeAllocVect


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine DeAlloc2DArray(self)
!-----------------------------------------------------------------------
!
!            deallocate the 2D array "self"
!
!-----------------------------------------------------------------------
use Memory
#include "priunit.h"
type(Array2Dim), intent(inout) :: self
call QENTER('DeAlloc2DArray')

if (iprint_vect >= 10 ) then
  write(LUPRI,'(/2X,A)') 'DeAlloc2DArray: going to deallocate the entering 2D-array:'
  call PrintObject(self)
  write(LUPRI,'(/2X,A)') 'DeAlloc2DArray: related memory parameters of this 2D-array:'
  call print_mem_param
endif

if (self%is_alloc) then
 call DeAllocVect(self%Vect2D)
 self%is_alloc = .false.
 self%dimens1 = 0;  self%dimens2 = 0
 self%p_rA2D => NULL(); self%p_iA2D => NULL(); self%p_lA2D => NULL()
else 
 call print_2Darray_attrib(self)
 call QUIT('DeAlloc2DArray: array was NOT allocated, thus cannot be deallocated !')
endif
if (iprint_vect >= 10) then
 write(LUPRI,*) 'DeAlloc2DArray: 2D array was deallocated in memory.'
 call print_mem_param
endif
call QEXIT('DeAlloc2DArray')
end subroutine DeAlloc2DArray

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function Access_p_iA1D(self) result(vect_p_iA1D)
!
! ... only test routine ...
!
! ... does not work for manipulating with indexes (reading/printing) 
! ... serves ONLY for passing POINTER
!
type(Vect), intent(in) :: self
!real(kind=8), dimension(:), pointer :: vect_p_rA1D => NULL()
integer, dimension(:), pointer :: vect_p_iA1D 
vect_p_iA1D => self%p_iA1D
end function Access_p_iA1D

END MODULE Vectors

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect1D_real(Vector,isize,iprint)
!-----------------------------------------------------------------------
!
!        Assign pointer to the real vector
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
#include "iratdef.h"
real(kind=8), target :: Vector(isize)
integer(kind=4), intent(in) :: iprint
integer(kind=8/IRAT),intent(in) :: isize
call QENTER('AssignVect1D_real')

if (iprint >= 10) then
  write(LUPRI,*) 'AssignVect1D_real: going to assing pointer....'
  CALL FLSHFO(LUPRI)
endif

p_real8 => Vector

if (iprint >= 10) then
  write(LUPRI,*) 'AssignVect1D_real: real 1dim pointer assigned'
  CALL FLSHFO(LUPRI)
endif

call QEXIT('AssignVect1D_real')
end subroutine AssignVect1D_real

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect1D_inte(Vector,isize,iprint)
!-----------------------------------------------------------------------
!
!          Assign pointer for integer vector
!
!   must be outside of the module !!!
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
#include "iratdef.h"
integer, target  :: Vector(isize)
integer(kind=4), intent(in) :: iprint
integer(kind=8/IRAT),intent(in) :: isize
call QENTER('AssignVect1D_inte')

if (iprint >= 10) then
  write(LUPRI,*) 'AssignVect1D_inte: going to assign pointer...'
  CALL FLSHFO(LUPRI)
endif

p_inte => Vector 

if (iprint >= 10) then
  write(LUPRI,*) 'AssignVect1D_inte: integer 1dim pointer assigned'
endif

call QEXIT('AssignVect1D_inte')
end subroutine AssignVect1D_inte

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect2D_real(Vector,iprint,n1,n2)
!-----------------------------------------------------------------------
!
! Assign pointer for real vector
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
!real(kind=8), dimension(:), pointer :: p_rA2D
!real(kind=8), dimension(:), intent(inout), target :: Vector
real(kind=8), target :: Vector(n1,n2)
integer(kind=4), intent(in) :: iprint,n1,n2
call QENTER('AssignVect2D_real')
p2_real8 => Vector 
if (iprint >= 10) then
 write(LUPRI,*) 'AssignVect2D_real: real 2dim pointer assigned'
endif
call QEXIT('AssignVect2D_real')
end subroutine AssignVect2D_real

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect3D_real(Vector,iprint,n1,n2,n3)
!-----------------------------------------------------------------------
!
!                 Assign pointer for 3D real array
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
integer(kind=4), intent(in) :: iprint,n1,n2,n3
real(kind=8), target :: Vector(n1,n2,n3)
call QENTER('AssignVect3D_real')
! ... assign pointer from GeneralPointers
p3_real8 => Vector 
!
if (iprint >= 10) then
 write(LUPRI,*) 'AssignVect3D_real: real 3dim pointer was assigned'
endif
call QEXIT('AssignVect3D_real')
end subroutine AssignVect3D_real

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect2D_inte(Vector,iprint,n1,n2)
!-----------------------------------------------------------------------
!
!     Assign pointer to the integer 2D-array
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
integer, target :: Vector(n1,n2)
integer(kind=4), intent(in) :: iprint,n1,n2
call QENTER('AssignVect2D_inte')
p2_inte => Vector 
if (iprint >= 10) then
 write(LUPRI,*) 'AssignVect2D_inte: integer 2dim pointer assigned.'
endif
call QEXIT('AssignVect2D_inte')
end subroutine AssignVect2D_inte

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignVect3D_inte(Vector,iprint,n1,n2,n3)
!-----------------------------------------------------------------------
!
!     Assign pointer to the integer 3D-array
!
!-----------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
integer(kind=4), intent(in) :: iprint,n1,n2,n3
integer, target :: Vector(n1,n2,n3)
call QENTER('AssignVect3D_inte')
p3_inte => Vector 
if (iprint >= 10) then
 write(LUPRI,*) 'AssignVect3D_inte: integer 2dim pointer was assigned.'
endif
call QEXIT('AssignVect3D_inte')
end subroutine AssignVect3D_inte

