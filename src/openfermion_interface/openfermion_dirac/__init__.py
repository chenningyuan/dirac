"""
OpenFermion plugin to interface with Dirac
"""
from ._run_dirac import run_dirac
from ._molecular_data_Dirac import (
        MolecularData_Dirac,
        load_molecular_hamiltonian,
        periodic_table)
