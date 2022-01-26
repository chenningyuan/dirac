# Quick start


## Working with the public and private repositories

- http://diracprogram.org/doc/master/development/multiple-remotes.html


## Building and testing

Clone the complete repository
```
$ git clone --recursive git@gitlab.com:dirac/dirac.git
```

or the latest revision only
```
$ git clone --depth 1 --recursive git@gitlab.com:dirac/dirac.git
```

Build the code
```
$ cd dirac
$ ./setup [--help]
$ cd build
$ make -j
```

Run the test set
```
$ ctest
```

DIRAC testing dashboard: https://testboard.org/cdash/index.php?project=DIRAC
