# Quick start


## Contributing changes

- [Working with the public and private repositories](http://diracprogram.org/doc/master/development/multiple-remotes.html)
- [External contributions](http://diracprogram.org/doc/master/development/contributing-changes.html)


## Building and testing

Clone the complete repository:
```
$ git clone --recursive https://gitlab.com/dirac/dirac.git
```

or the latest revision only:
```
$ git clone --depth 1 --recursive https://gitlab.com/dirac/dirac.git
```

or a particular release:
```
$ git clone --recursive --branch v22.0 https://gitlab.com/dirac/dirac.git
```

If you have downloaded a zip/tar.gz/tar.bz2/tar archive directly from GitLab,
you will require Git installed to build the code. After extracting the archive
you will need to run this additional step in order to fetch and update external
dependencies:
```
$ git submodule update --init --recursive
```

**However, we recommend users who do not plan to make any code changes
to download from** https://doi.org/10.5281/zenodo.4836495.
In this case the above extra step is not needed as the release tarballs are
self-contained, and in addition we get useful download metrics which we can use
in future funding applications.

Build the code:
```
$ cd dirac
$ ./setup [--help]
$ cd build
$ make -j
```

Run the test set:
```
$ ctest
```

DIRAC testing dashboard: https://testboard.org/cdash/index.php?project=DIRAC
