name = "gen_mocoeff"

from ._run_dirac_mocoef import run_dirac_mocoef, average_of_configuration_openshell
from ._moldata import moldata
from ._AOMObases import BasisFunc, AOBasis, MO, MOBasis

from ._generate_adf_input import generate_adf_input

from ._utils import *

# Avoid psi4 import error in all tests if not installed
try:
 from ._run_psi4_mocoeff import run_psi4_mocoeff
 has_psi4 = True
except ImportError:
 has_psi4 = False

try:
 from ._run_pyscf_mocoeff import run_pyscf_mocoeff
 has_pyscf = True
except ImportError:
 has_pyscf = False

# is using run_pyscf_mocoeff and run_psi4_mocoeff:
from ._combine_MOs_B2 import *
