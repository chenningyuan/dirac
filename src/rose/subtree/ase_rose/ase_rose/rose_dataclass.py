"""
This module defines customized input options for Rose
"""
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Tuple, List, Union, Sequence
import numpy as np


@dataclass
class ROSEBaseTopology:
    """Abstract base class for atom/molecule topology in ROSE."""
    name: str = None
    atoms: Optional[
        List[Tuple[Union[str, np.ndarray], np.ndarray]]
    ] = field(default_factory=list)
    charge: int = 0
    multiplicity: int = 1
    basis: str = 'sto-3g'
    uncontract: bool = False

@dataclass
class ROSEFragment (ROSEBaseTopology):
    """A dataclass representing an atomic or molecular fragment in ROSE."""


@dataclass
class ROSETargetMolecule (ROSEBaseTopology):
    """A dataclass representing the target molecular system in ROSE."""


class ROSECalcType(Enum):
    """Enumerator class defining the type of calculation done by ROSE."""
    ATOM_FRAG = 'atom_frag'  # calc IAOs from atom-free AOs.
    MOL_FRAG = 'mol_frag'    # calc IFOs from free molecular fragment MOs.


class ROSEMOCalculators(Enum):
    """Enumerator class defining the program to which ROSE is interfaced."""
    PYSCF = 'pyscf'
    PSI4 = 'psi4'
    DIRAC = 'dirac'
    ADF = 'adf'
    GAUSSIAN = 'gaussian'


class ROSEIFOVersion(Enum):
    """Enumerator class defining the 'recipe' to obtain IFOs."""
    STNDRD_2013 = "Stndrd_2013"
    SIMPLE_2013 = "Simple_2013"
    SIMPLE_2014 = "Simple_2014"


class ROSEILMOExponent(Enum):
    """Enumerator class defining the exponent used to obtain ILMOs."""
    TWO = 2
    THREE = 3
    FOUR = 4


class ROSEDIRACHamiltonian(Enum):
    """Enumerator class defining the type of Hamiltonian used in DIRAC."""
    NONREL = "NONREL"
    X2C = "X2C"
    X2CMMF = "X2Cmmf"


@dataclass
class ROSEInputDataClass:
    """A dataclass representing input options for ROSE."""
    rose_target: ROSETargetMolecule
    """Target supermolecule => ASE Atoms object."""
    rose_frags: Union[ROSEFragment, Sequence[ROSEFragment]]
    """List of atomic or molecular fragments."""
    target_mo_file_exists: bool = False
    frags_mo_files_exist: bool = False
    """Are canonical mo files already present?"""
    rose_calc_type: ROSECalcType = ROSECalcType.ATOM_FRAG.value
    """Calculation type"""
    rose_mo_calculator: ROSEMOCalculators = (
        ROSEMOCalculators.PYSCF.value
    )
    """ROSE molecular orbitals calculator"""

    exponent: ROSEILMOExponent = ROSEILMOExponent.TWO.value
    """Exponent used in the localization procedure."""
    restricted: bool = True
    """Restricted calculation option"""
    openshell: bool = False
    """Open shell calculation option"""
    spatial_orbitals: bool = True
    """Spatial orbitals option - true for non-relativistic calculations."""
    test: bool = True
    """Test to calculate HF energy with the final localized orbitals."""
    include_core: bool = False
    """Frozen core option."""
    relativistic: bool = False
    """Relativistic calculation option."""
    version: ROSEIFOVersion = ROSEIFOVersion.STNDRD_2013.value
    """Version of the IAO construction."""
    spherical: bool = False
    "Spherical or cartesian coordinate GTOs."
    save_data: bool = True
    """Option to save final iaos/ibos in a chk file."""
    fcidump: bool = False
    """Option to save final iaos/ibos in a fcidump file"""
    additional_virtuals_cutoff: Optional[float] = None  # 2.0  # Eh
    """Add virtual orbitals with energies below this treshold."""
    frag_threshold: Optional[float] = None  # 10.0  # Eh
    """Set reference virtuals as those with energies below this treshold."""
    frag_valence: Optional[List[List[int]]] = field(default_factory=list)
    """# of valence electrons per fragment."""
    frag_core: Optional[List[List[int]]] = field(default_factory=list)
    """# of core electrons per fragment."""
    frag_bias: Optional[List[List[int]]] = field(default_factory=list)
    """Bias frags when assigning the loc orb (for non-polar bonding orb)."""

    avas_frag: Optional[List[int]] = field(default_factory=list)
    """Fragment IAO file to be extracted from ROSE."""
    nmo_avas: Optional[List[List[int]]] = field(default_factory=list)
    """List of spatial MOs (or spinors if restricted = False) to AVAS."""

    # DIRAC specific options
    dirac_hamiltonian: ROSEDIRACHamiltonian = ROSEDIRACHamiltonian.NONREL.value
    """DIRAC hamiltonian."""
    dirac_special_basis: Optional[List[str]] = None
    """Used to set basis for some atoms, e.g., ["STO-3G","H cc-pVDZ"]."""
    symmetry: bool = False
    """Use symmetry or not (GENIBO works without symmetry for now)."""
    dirac_spinfree: bool = False
    """Use Dyall's spin-free Hamiltonian."""
    dirac_point_nucleus: bool = True
    """Use the point-nucleus model."""
    dirac_speed_of_light: bool = False
    """Set the speed of light in a.u. (default is 137 a.u.)."""
    get_mrconee: bool = False
    """Extract the MRCONEE file containing the one-electron integrals."""

    def __post_init__(self):
        """This method is called after the instance is initialized."""
        if self.save_data:
            pass
