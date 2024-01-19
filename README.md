# Quick start

## Homepage

More information can be found on the DIRAC homepage: https://www.diracprogram.org/doku.php

There is also a google group: https://groups.google.com/g/dirac-users

## Getting the Code

- **We recommend users who do not plan to make any code changes
to download from** https://zenodo.org/doi/10.5281/zenodo.7670748

  In this case the release tarballs are self-contained, and in addition we get useful download metrics which we can use
in future funding applications.

- If you clone the code directly from GitLab, please clone with ``--recursive``, e.g. the latest master::

  ```
  $ git clone --recursive --branch master https://gitlab.com/dirac/dirac.git
  ```

  By replacing ``master`` with a version number ``v23.0`` one can select a specific version.
  In order to get the latest revision only:
  ```
  $ git clone --depth 1 --recursive git@gitlab.com:dirac/dirac.git
  ```

- If you have downloaded a zip/tar.gz/tar.bz2/tar archive directly from GitLab,
you will require Git installed to build the code. After extracting the archive
you will need to run this additional step in order to fetch and update external
dependencies:
  ```
  $ git submodule update --init --recursive
  ```

## Building

The default installation (sequential) proceeds through four commands:
```
$ cd dirac
$ ./setup [--help]
$ cd build
$ make -j
```

DIRAC is configured using CMake https://cmake.org/ , typically via the ``setup`` script,
and subsequently compiled using make (or gmake).
The ``setup`` script is a useful front-end to CMake.
You need python to run ``setup``. To see all options, run
```
  $ ./setup --help
```

The ``setup`` script creates the directory "build" and
calls CMake with appropriate environment variables and flags.
By default CMake builds out of source. This means that all object files and the
final binary are generated outside of the source directory. Typically the build
directory is called "build", but you can change the name of the build directory
(e.g. ``build-gfortran``)
```
  $ ./setup [--flags] build-gfortran
  $ cd build-gfortran
  $ make
```

## Testing

We strongly recommend that the installation is followed by testing your installation:
```
$ ctest
```

or
```
$ make test
```

If you encounter any problems check the issues under https://gitlab.com/dirac/dirac and create a new one if your problem is not mentioned.

## Contributing

- [Working with the public and private repositories](https://diracprogram.org/doc/master/development/multiple-remotes.html)
- [External contributions](https://diracprogram.org/doc/master/development/contributing-changes.html)
- [DIRAC testing dashboard](https://testboard.org/cdash/index.php?project=DIRAC)
