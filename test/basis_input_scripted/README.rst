Basis set files checking
------------------------

The DIRAC test **test/basis_input_scripted** is dedicated for automated checking of (almost all) 
basis set sets files sitting in the DIRAC directory.

How it works
~~~~~~~~~~~~

The tests performs reading of multiple basis set files, 
and merely checks that the input parser does not choke on picked basis set. No reference outputs are involved.

The principal Python script 
:download:`test <../../../test/basis_input_scripted/test>`
is employing the runtest library (https://runtest.readthedocs.io).

When the typical string sequence is found in the basis set file, the script extracts the proton number of a given element.
Script then launches run(s) with the
:download:`input <../../../test/basis_input_scripted/input_test_basis.inp>`, 
:download:`molecule <../../../test/basis_input_scripted/element.mol>` files and 
with pam parameters containing selected basis set name and element's proton number as strings for replacements by Python.

Default run
~~~~~~~~~~~

Test's default run checks all basis set files in given day once per month. 
This is achieved by running through basis set directories where each file with several basis sets inside is read. 
The total number of individual runs makes more than 5000, and it takes about 20 minutes.

On other days the test makes empty run. 

Examples of usage
~~~~~~~~~~~~~~~~~

User can specify by environment variables (BSFILE, ZELEM, RANDOM_BSF_Z) on what basis files/elements to run.

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

Because of many basis set files (more than 300), one can not expect all of them working.
There are cases of wrong basis sets, let us mention few:
- does not work for basi set files like "basis_dalton/aug-pc-0_emsl" where the proton number can not be extracted.
- ...
