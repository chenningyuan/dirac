Description of the test procedure
=================================


How to run tests
----------------

DIRAC uses the runtest script to run tests.
Go to the source root::

  $ cd /path/to/dirac/source

Then get a list of options::

  $ ./runtest --help

To run all tests, type::

  $ ./runtest --all

You can also run individual tests. Have a look at the options.

If you are in the subdirectory of an individual test (example: dft_ac), you can
do this instead to run just this one test::

  $ [cd test/dft_ac]
  $ ../../runtest


Some tests are broken, what now?
--------------------------------

You are a DIRAC user? Then please inform the DIRAC developers.

You are a developer? Perhaps it is not your fault. Have a look at the nightly
testing dashboard where you can check whether this test is marked as red
(broken) or green (passing): http://repo.ctcc.no/CDash/index.php?project=DIRAC
If it is green there but broken on your machine it can either mean that the
test has trouble with your architecture/compiler (in this case inform other
developers) or that you broke it with recent changes - then you have to go in
and fix it before you git push to origin/master.


Suddenly some red tests appear on the nightly testing dashboard, what happens now?
----------------------------------------------------------------------------------

Miro or Radovan or somebody else examining the nightly testing dashboard will
notice it and write to the developer who is likely to have introduced the
conflicting change.


How to add a new test
---------------------

1. Copy an existing test and choose a sensible name, e.g. not
   my_favorite_test or test2, but for instance symmetry_recognition to signal that
   we are testing the symmetry recognition procedure.

2. Put the inp and mol file(s) in the directory and modify the
   menu file to work with these filenames.

3. Put (a) correct output(s) in the result directory and modify the filter file to
   take out the part that is relevant for the functionality that you want to test.
   Don't test everything, e.g. the electronic Hartree-Fock energy is already
   tested at several places, no need to do this in every test.

4. Go to the source root directory and run the script runtest with your new test.  It
   is a good practice to test your test, try to introduce deliberately errors in
   your new code and check whether your new test catches the error - if not,
   modify the test. The script will complain if your filter does not filter out anything.

5. If the test works, look at how long it takes to run on a single PC to see where
   it should be quick or regular. A rule of thumb: < 10 seconds = quick, otherwise regular.

6. Add the test to pamadm (search for existing tests to see how).

7. Use git add and git commit to commit your work to the repository.  Make sure
   to include ALL files, don't forget the result directory that should contain the
   reference output and filter files. You can verify this with "git status".


How does the filtering work?
----------------------------

If one would use simply the unix command diff to compare the new output with
the reference output, one would always see differences, simply because dates
and execution times will be different, because insignificant last digits are
different, or because output formatting has changed.

To focus on relevant differences runscript uses the script crep to find the
relevant numbers that are to be checked. The idea is that you define a
characteristic text that is contained in the output just BEFORE the numbers
that are to be checked. Here you may either use a string or a regular
expressions, examples of both can be found when looking at the filter files.
This defines the start of the output section to be tested, but you also want to
define the end.

Here there are two possibilities, you either know in advance how many lines you
want to test, or you don't (keep in mind that others may later introduce
additional blank or text lines). In the first case you can simply specify the
number of lines to be test (using NLines), in the second case you have to
specify a characteristic text that closes the section.
