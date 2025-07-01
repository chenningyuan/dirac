This package offers a Python-based ASE (Atomic Simulation Environment) calculator specifically designed for DIRAC.
It is an integral component of the DIRAC package, designed to facilitate its use through an intuitive Python interface.


# INSTALLATION:

```
$> pip install -e .
```

Or, if you don't have permission:

```
$> pip install --user -e .
```

Also, if you would like to run the DIRAC-ASE tests
```
$> pip install --user -e ".[testing]"
$> pytest -v tests
```

or 

```
$> pip install --user -e ".[testing]"
$> pytest -v tests
```

Make sure that the script pam (generated in the build directory) is in your $PATH environment variable.

# RUN:

A simple demonstration that computes the MP2 energy of three molecules from the ASE molecules database is provided in the tests directory.
