Programmers' manual - Frequently asked questions
================================================


How can I get the AO density matrix?
------------------------------------

Use subroutine DENMAT dirac/dirden.F.


How should I program input reading and introduce new keywords?
--------------------------------------------------------------

Perhaps 98% of the input reading is done using tables and gotos.

Please do not introduce new keyword sections like that for at least three reasons:

1. It is a nightmare to read and maintain.
2. There is no type checking so the user
   gets useless error messages or unexpected results.
3. New scheme minimizes conflicts when merging branches.

Rather have a look in dirac/input_reader.F90.
This is the modern way to do it.


How can I get MO coefficients?
------------------------------

Use subroutine REACMO_NO_WORK
in dirac/dirgp.F.


How are density/Fock/... matrices blocked?
------------------------------------------

They are blocked according to NBBAS and IBBAS.


How should I allocate memory?
-----------------------------

Large parts of the DIRAC code use MEMGET/MEMREL. Please do not allocate
memory using MEMGET/MEMREL. Use routines alloc/dealloc.
We are trying to get rid of MEMGET/MEMREL.


Why do all routines end with "RETURN" before "END"?
---------------------------------------------------

Most of the routines end with::

  subroutine foo()
  ...
  return
  end

The "return" is redundant. Simply write "end" without "return".  To make it
better, we recommend programmers to stick to the Fortran90/95 coding standards.
So the structure of the routine shall be as::

  subroutine foo()
  ...
  end subroutine


What are all the "Decks" before the subroutines good for?
---------------------------------------------------------

Historical reasons.
They have no meaning today.


Can I use zgemm or do I have to use the module matrix_defop?
------------------------------------------------------------

There is nothing wrong with zgemm. You can use either one - whichever you prefer.


I want my new module to be completely modular. Can I then use DIRAC infrastructure routines like QUIT or READT?
---------------------------------------------------------------------------------------------------------------

If you want to be 100% independent of the DIRAC infrastructure then you cannot
use such routines.


I like to use implicit.h to avoid variables declaration
-------------------------------------------------------

The "implicit.h" is bad programing practice and I strongly recommend
to always use implicit none. 

Reasons:
1. With "implicit.h" it is easy to use variables that are undefined.

 example:

 call daxpy(size(array1), D1, array1, 1, array2, 2)

if one forgets to define the parameter D1, the code will compile
and D1 can be any number so the result is undefined.
With "implicit none" this will not compile.

Another example is a typo:
istead of ::

 NDIM = 12  ! correct code

you could type ::

 NDIN = 12  ! typo
 
With implicit none such mistake won't compile, with "implicit.h" the code will compile
and behave unexpectedly. In this case also the compiler flag "-Wuninitialized" won't help!

In the code we had (and unfortunatelly still have( many bugs like this and some of us have spent several long afternoons
chasing them. The insertions of "implicit.h"  is too dangerous and also experienced programmers
make typos.

2. The insertion of "implicit none" is nice for other people reading your code. They
can see which variables are local, which are global (from common blocks).
Without implicit none you have no idea whether this is a local or a global
variable if you don't know the common blocks by heart.

3. The  command "implicit none" makes it easier to identify include statements that are  
really 
used and includes that are just included but not used.
With "implicit.h"  you can comment out includes and the code may still  
compile.
Contraty, "implicit none" will warn you. If you remove an include file using "implicit none"
and the code compiles, then this include was useless. 

How about integers ?
--------------------
Please  in your codes always use *KIND free* declarations ::
 integer :: variable
instead of *KIND restricted* integer variables :: 
 integer(KIND=4) :: variable
 integer(KIND=8) :: variable

This is because we need flexibility for integer variables: compilation
of the DIRAC suite goes also with predefined integer*8 variables 
(through -i8 Fortranflag), 
for what all integer variables in the code must be declared without KIND 
restriction.

