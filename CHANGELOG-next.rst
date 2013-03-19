

12.6
----

- Add basis_ecp to "make install" target.
- Make it possible to compile without XCFun.
- Enable creation of static executables with older linkers.
- Add --replace option to the pam script.
- introduce keyword .THRECI in KR-CI input to allow the user to change the default CI energy convergence criterion,
  present default::

  .THRECI
  1.0d-10
- Bugfix for MPI interface on Cray/MPICH2 systems (thanks to Juho for pointing out the problem)
