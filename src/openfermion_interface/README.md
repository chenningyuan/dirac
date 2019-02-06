# Dirac interfaced with Openfermion

# Requirements

- Dirac 
- Openfermion
- all requirements for the aforementioned programs.

# Installation

INSTALL DIRAC: http://diracprogram.org/doku.php

INSTALL OPENFERMION : https://github.com/quantumlib/OpenFermion

INSTALL Openfermion-Dirac interface:

Clone the complete repository:
```
$ git clone https://github.com/bsenjean/Openfermion-Dirac.git
```

Build the code:
```
$ cd /path/to/Openfermion-Dirac/
$ pip install -e .
$ cd utils/
```

In /path/to/Openfermion-Dirac/openfermion_dirac/_run_dirac.py, set your own path to
- set/your/path/to/pam (pam is the run_script of the Dirac program), which is normally in /path/to/dirac/build.
- set/your/path/to/dirac_openfermion_mointegral_export.x, which is normally in /path/to/dirac/build.

Note that if you set the /path/to/dirac/build in your bash_profile, you can simply write pam and dirac_openfermion_mointegral_export.x instead of the whole path to them.

# Use

Different examples are furnished in the examples/ repository in python, as well as a tutorial in tutorial/. If one wants to play more with the tutorial, use jupyter notebook to do so.
