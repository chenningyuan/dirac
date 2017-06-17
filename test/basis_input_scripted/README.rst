Basis set files checking
------------------------

The  **test/basis_input_scripted** is dedicated for automated checking of (almost all) basis set sets files sitting in the DIRAC directory.

It does only the reading (.INPTEST) of preprocessed basis sets, 
and merely tests whether input parser does not choke on gives basis set. No reference outputs are involved.

When the typical string sequence is found in the basis set file, the proton number of a given element is extracted 
and used for the generic mol file. Afterward, the runtest library (https://runtest.readthedocs.io) based python script
:download:`test <../../../test/basis_input_scripted/test>`
launches Dirac run(s) with
:download:`input <../../../test/basis_input_scripted/input_test_basis.inp>`, 
:download:`mol <../../../test/basis_input_scripted/element.mol>` files and 
with further pam parameters containing selected basis set name and element's proton number as strings for replacements.

Default run
~~~~~~~~~~~

Test's default run checks all basis set files in given day once per month. 
This is achieved by running through basis set directories where each file with several basis sets is read. 
The total number of individual runs is more than 5000, and it takes more than 20 minutes.

On other days it makes empty run. 

User can specify by environment variables (BSFILE, ZELELE, RANDOM_BSF_Z) what basis files/elements to run.


Examples of usage
~~~~~~~~~~~~~~~~~

Run whole selected basis sets only:

::

 export BSFILE="basis_ecp/ECPDS28SDFSF"
 ./test -b <BUILDIR> 

Run only selected elements in selected basis sets:

::

 export BSFILE="basis_ecp/ECPDS28SDFSF"
 export ZELEM="30 35"
 ./test -b <BUILDIR> 

Run all elements in two randomly selected basis sets:

::

 export RANDOM_BSF_Z="2"
 ./test -b <BUILDIR> 

Run three elements in three randomly selected basis sets:

::

 export RANDOM_BSF_Z="3 3"
 ./test -b <BUILDIR> 


Nonworking cases
~~~~~~~~~~~~~~~~

Because of many basis set files (more than 300), one can not expect all working perfectly.

There are cases of wrong basis sets, let us mention few:

- does not work for basi set files like "basis_dalton/aug-pc-0_emsl" where the proton number can not be extracted.


