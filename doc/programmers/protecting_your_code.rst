

General
-------

When you add your bleeding edge code to the master branch, there
are a number of things you try to make sure:
- That you do not break compilation for other developers because of missing dependencies.
- That your code does not get released in any way (official release or simply "make release").
- That your new tests are run iff your module is compiled and that they exit in a controlled way otherwise.


Conditional compilation
-----------------------

There are two complementary ways to introduce conditional compilation. Either
in pamadm (which sets CPP flags and selects directories to compile) and/or in
CMakeLists.txt where you can set CPP flags, too.  You can use both to your
maximum convenience but you should know that only pamadm controls the selection
of tests and the packaging of sources. So if you want to conditionally select
tests and avoid that your code gets released, you **have to** use pamadm for
this. Read the pamadm sources and look at other modules to understand how
directories are selected and tests packaged.  If you want to make sure that
your module does not end up in the release, verify it yourself. Do a "pamadm
--release; cd build; cmake ..; make release", extract the tarball and verify.
Do not rely on anybody to do that for you.


Should my module be part of "pamadm -a"?
----------------------------------------

Ideally yes. Ideally all modules should be part of "pamadm -a".  But in
practice it can be reasonable to exclude your module from "pamadm -a" if it
makes compilation for other developers very difficult or impossible.  Note that
most developers will compile with "pamadm -a". If you exclude it from "pamadm
-a", then you need to select your sources with "pamadm -a --yourmodule".
