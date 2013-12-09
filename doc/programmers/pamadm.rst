

Documentation for pamadm
========================


What are the pamadm modules and why do we have them?
----------------------------------------------------

First of all the pamadm modules have nothing to do with Fortran modules.  With
the pamadm modules it is possible to enable/disable the compilation and testing
of entire directories. A pamadm module will not only select a directory to
compile or ignore, but also define CPP definitions (MOD_FOO), and select tests.
Compilation depends on the selected pamadm modules (pamadm generates list of
sources to compile and sets compile definitions) It is important to realize
that when we release or export the code with "make release", then the result of
it depends on the pamadm flags. Code that is then not activated does not make
it to the release tarball. Same for tests.


Should I introduce a new pamadm module for my new functionality?
----------------------------------------------------------------

It can be useful to introduce a new pamadm module for your new functionality if
you want to have the functionality on master but you want to be able to disable
its release. Note that this can be achieved also with using MOD_UNRELEASED.

Personal opinion (Radovan): pamadm modules are an excuse to not use branches
properly. If we used feature branches properly and only merged to-be-released
code to master, we would not need pamadm modules.
