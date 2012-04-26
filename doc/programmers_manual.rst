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
