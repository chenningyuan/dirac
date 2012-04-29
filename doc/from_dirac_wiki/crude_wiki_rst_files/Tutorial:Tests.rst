Tutorial:Tests
==============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

It is *very important* that you verify that your DIRAC installation
correctly reproduces the reference test set *before* running production
calculations.

To get documentation on the test set, run:

::

    ./runtest --help

As a starter try the quick test that covers most of the standard run
types:

::

    ./runtest --quick

If everything goes well the script will report calculations completed
and provide a summary:

::

     summary:   114 ok

     total time: 00h03m42s

If there is a problem the script will report that one or more tests
failed (wrong output) or crashed (execution error). To learn more go to
the directory **logs** where you will find a text file **\*.txt** with a
short summary of the test that failed. More importantly, there will be a
html-file **\*.htm** to be opened in a browser. Selected parts of your
output are compared to the reference output found in the **result**
subdirectory of your test and a coloring scheme in the html-file will
indicate where things go wrong.

Once you have identified the source of the problem, you may rerun
individual tests with:

::

    ./runtest --tests="test1 test2"

where test1 and test2 should be replaced by the name(s) of the test(s)
that you want to rerun: When everything seems fine try also the complete
test set:

::

    ./runtest --all

