This package offers a Python-based ASE (Atomic Simulation Environment) calculator specifically designed for DIRAC.
It is an integral component of the DIRAC package, designed to facilitate its use through an intuitive Python interface.


# INSTALLATION:

Before installing, make sure the `pam` script generated in the DIRAC build directory is available on your `$PATH`. Installation will abort if `pam` cannot be found.

For example:

```bash
$> cd ../build
$> export PATH=$(pwd):$PATH
$> cd ../ase_dirac
```

Install the package in the same Python environment that you will use to run the tests. The test dependency `pytest` is installed automatically.

```
$> pip install -e .
```

Or, if you don't have permission:

```
$> pip install --user -e .
```

To run the DIRAC-ASE tests:

```
$> python -m pytest -v tests
```

# TROUBLESHOOTING:

If `pytest` fails during import with

```text
NameError: name 'DIRAC' is not defined
```

you are likely using an older copy of `ase_dirac`. Reinstall the package from the current source tree:

```bash
$> pip install -e .
```

and rerun the tests with the same interpreter:

```bash
$> python -m pytest -v tests
```

If installation fails with an error saying that `pam` was not found, add the DIRAC build directory containing `pam` to your `$PATH` and rerun `pip install`.

# RUN:

A simple demonstration that computes the MP2 energy of three molecules from the ASE molecules database is provided in the tests directory.
