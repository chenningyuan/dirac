Reference test set
==================

It is *very important* that all functionality of DIRAC is thoroughly tested. This implies 
that whenever you add some new feature to DIRAC you should prepare an appropriate test for it as well.

How to add a new test
---------------------

1. Copy an existing test and choose a sensible name, e.g. not
   my_favorite_test or test2, but for instance symmetry_recognition to signal that
   we are testing the symmetry recognition procedure.

2. Put the input files in the directory and modify the
   menu file to work with these filenames.

3. Put (a) correct output(s) in the result directory and modify the filter file
   (see further below) to take out the part that is relevant for the functionality
   that you want to test.  Don't test everything, e.g. the electronic Hartree-Fock
   energy is already tested at several places, no need to do this in every test.

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

To focus on relevant differences runscript extracts the relevant numbers
according to keywords in result/filter.out.  Here you may either use strings or
regular expressions, examples of both can be found when looking at the filter
files.  There are two possibilities, you either know in advance how many lines
you want to test, or you don't (keep in mind that others may later introduce
additional blank or text lines). In the first case you can simply specify the
number of lines to be test (using nlines), in the second case you have to
specify a characteristic text that closes the section.

Note that runtest compares relative errors (DIRAC) and absolute errors (Dalton).
The default accepted relative/absolute error (DIRAC/Dalton) is 0.0001.
You can modify this using $ERROR. If you want to ignore sign, use $IGNORESIGN.
Also not that $ERROR and $IGNORESIGN are only valid for the next subtest.

Example 1: 
 ::
  $ERROR 0.000000001
  nlines 8:Electronic energy

  $ERROR 1.E-5
  from:Eigenvalues
  to:HOMO - LUMO

In this example we have two subtests. First we test Electronic energy and 7
following lines using tolerance 0.000000001, then we test Eigenvalues using
tolerance 1.E-5.  We need to specify 1.E-5 for the second subtest otherwise
runtest defaults to 0.0001.

Example 2:
 ::
 $ERROR 1.E-6
 $IGNORE_BELOW 1.E-20
 nlines 16:Expectation value for individual orbitals

In this case we compare expectation values for individual orbitals, and allow 
maximum differences of 1.E-6. We also specify "zero" limit for compared values of 1.E-20 (as $IGNORE_BELOW),
what means numbers below this threshold are simply ignored.

List of all possible filter keywords is here: ::

  $ERROR          # accepted error
  $IGNORE_BELOW   # ignore if numbers are below limit
  $IGNORESIGN     # ignore sign
  from:STRING     # start from string STRING
  fromre:STRING   # start from regex STRING
  to:STRING       # stop at string STRING, do not include end line
  toi:STRING      # stop at string STRING, include end line
  tore:STRING     # stop at regex STRING, do not include end line
  toire:STRING    # stop at regex STRING, include end line
  nlines N:STRING # consider line containing STRING and N-1 following lines
