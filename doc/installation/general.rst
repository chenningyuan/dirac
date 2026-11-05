:orphan:
 

Basic installation
==================


Getting the Code
----------------

**We recommend users who do not plan to make any code changes
to download from** https://doi.org/10.5281/zenodo.18662049 (DIRAC26)

In this case the release tarballs are self-contained, and in addition we get useful download metrics which we can use
in future funding applications.


If you clone the code directly from GitLab, please clone with ``--recursive``, e.g. the latest master::

    $ git clone --recursive --branch master https://gitlab.com/dirac/dirac.git

By replacing ``master`` with a version number, for instance ``v26.0``, one can select a specific version.
In order to get the latest revision only::

    $ git clone --depth 1 --recursive git@gitlab.com:dirac/dirac.git


If you have downloaded a zip/tar.gz/tar.bz2/tar archive directly from GitLab,
you will require Git installed to build the code. After extracting the archive
you will need to run this additional step in order to fetch and update external
dependencies::

  $ git submodule update --init --recursive


Building
--------

The default installation (sequential) proceeds through four commands::

  $ cd dirac
  $ ./setup
  $ cd build
  $ make -j

DIRAC is configured using `CMake <https://www.cmake.org/>`_ , typically via the ``setup`` script,
and subsequently compiled using make (or gmake).
The ``setup`` script is a useful front-end to CMake.
You need python to run ``setup``. To see all options, run::

  $ ./setup --help

The ``setup`` script creates the directory "build" and
calls CMake with appropriate environment variables and flags.
By default CMake builds out of source. This means that all object files and the
final binary are generated outside of the source directory. Typically the build
directory is called "build", but you can change the name of the build directory
(e.g. "build_gfortran")::

  $ ./setup [--flags] build_gfortran
  $ cd build_gfortran
  $ make

You can compile the code on all available cores::

  $ make -j


Testing
-------

We strongly recommend that the installation is followed by testing your installation::

  $ ctest 

or ::

  $ make test

For more information about testing, see `here <testing.html>`_.

If you encounter any problems check the issues under `DIRAC <https://gitlab.com/dirac/dirac>`_ and create a new one if your problem is not mentioned.


Contributing
------------

- `Working with the public and private repositories <https://diracprogram.org/doc/master/development/multiple-remotes.html>`_
- `External contributions <https://diracprogram.org/doc/master/development/contributing-changes.html>`_
- `DIRAC testing dashboard <https://testboard.org/cdash/index.php?project=DIRAC>`_

Parallel build
--------------

Building DIRAC for parallel runs using `MPI <http://www.mpi-forum.org/docs/docs.html>`_ 
is in principle a very simple modification of the above procedure::

  $ ./setup [--flags] --mpi
  $ cd build
  $ make

Once again we recommend `testing the code <testing.html>`_ afterwards.


Typical examples
----------------

In order to get familiar with the configuration setup, let us demonstrate
some typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --mpi --fc=mpif90 --cc=mpicc --cxx=mpicxx

These compiler names are in fact default for MPI in Dirac, so there is a shortcut if these
are the compiler wrappers you want to use::

  $ ./setup --mpi

If the wrappers ``mpif90``, ``mpicc`` and ``mpicxx`` are not found in ``PATH``,
or if no ``--mpi`` flag is given, compilers that are not specified explicitly are
detected by CMake itself (respecting the environment variables ``FC``, ``CC``
and ``CXX``). On Cray systems, for instance, CMake picks the compiler wrappers
``ftn``, ``cc`` and ``CC``, so there it is sufficient to run::

  $ ./setup --mpi

Configure for sequential compilation using ifx/icx/icpx and link against parallel mkl::

  $ ./setup --fc=ifx --cc=icx --cxx=icpx --qmkl=parallel

Configure for sequential compilation using gfortran/gcc/g++::

  $ ./setup --fc=gfortran --cc=gcc --cxx=g++

You get the idea. The configuration is usually good at detecting math libraries
automatically, provided you export the proper environment variable ``MATH_ROOT``,
see :ref:`linking_to_math`.


What to do if CMake is not available or too old?
------------------------------------------------

If it is your machine and you have an Ubuntu or Debian-based distribution::

  $ sudo apt-get install cmake

On Fedora::

  $ sudo yum install cmake

Similar mechanisms exist for other distributions or
operating systems. Please consult Google.

If it is a cluster, please ask the Administrator to install/upgrade CMake.

If it is a cluster, but you prefer to install it yourself (it's easy):

1. `Download <https://www.cmake.org/cmake/resources/software.html>`_ the latest pre-compiled CMake tarball
2. Extract the tarball
3. Set correct PATH variable

What to do if hdf5 is not installed on your machine
---------------------------------------------------

You can download the latest version of hdf5 from `the hdfgroup website <https://www.hdfgroup.org/download-hdf5>`_

