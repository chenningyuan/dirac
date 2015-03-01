Good Fortran90 programming practices
====================================

Example module
--------------

.. literalinclude:: nice_module.F90
      :language: fortran

I/O control statements
----------------------

To ease the localization of input/output-type bugs during the code run please provide main I/O statements (OPEN, READ)
with the IOSTAT parameter evaluation or with other error-trapping programing scheme.

Examples: ::

  !... control reading if gaunt=true
  READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT,IGTINT
  IF (IOS.NE.0) THEN
    WRITE(LUPRI,'(/,2X,A)') 'Error in SCF INTFLG reading !'
  & //'4 parameters needed for Gaunt term !'
    CALL FLSHFO(LUPRI)
    CALL QUIT( 'DHFINP: Error in INTFLG reading for GAUNT=true!')
  ENDIF

or (assuming that LU is correct unit number:  ::

  WRITE(LU,desired_format,IOSTAT=IOS) SOMETHING
  IF (IOS.NE.0) THEN
  ! Yell, that some variables in "SOMETHING" are out of desired format-bounds ;
  ! this could indicate run-time error (getting variable(s) hurt before)
  ! might activate program stop (CALL QUIT)....
  ENDIF

F90 INCLUDE keyword
-------------------

Do not use it. Use #include preprocessor statements only.

SIXLTR SBRTNE NAMING OBSOLT
---------------------------

Because it is so hard to read and understand.
There is no need to to stick to 6 letter names anymore.
The limit is 31 characters. Use meaningful names so that other people (and you in a couple of months)
will understand your code.

Using dynamical strings
------------------------

The programmer may utilize short dynamical strings carrying descriptive labels, and attach this
data type to his own data structures.

The module is ::
 use mystring

Methods at hand are: 

i)  creating the string; ::

 type(string) :: this_string
 call make_string(this_string,'some description text')

ii) printing the string on the LUPRI channel ::

 call print_string(this_string)

iii) copying one string into another ::

  call copy_string(this_string, another_string)

iv) canceling (deallocating) the string ::

 call cancel_string(this_string)

For more details check the  *trunk/src/gp/dynamic_string.F90* source file.
