"""This module defines an ASE interface to ROSE.

Original paper & official release:
https://pubs.acs.org/doi/10.1021/acs.jctc.0c00964

GitLab page:
https://gitlab.com/quantum_rose/rose

Note: see also https://pubs.acs.org/doi/10.1021/ct400687b.
"""
import subprocess
import os
from typing import List
import copy
import numpy as np

from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, OldShellProfile

from .rose_dataclass import ROSEInputDataClass, ROSEMOCalculators, ROSECalcType
from .rose_io import (
    read_rose_out,
    write_rose_genibo_in,
    read_int,
    read_real_list
)


def _import_genibo(package):
    """Import required modules specific to qchem workflows."""
    try:
        from genibo_MOcoeff._moldata import moldata
        if package == 'pyscf_genibo':
            from genibo_MOcoeff._run_pyscf_mocoeff import run_pyscf_mocoeff, replace_pyscf_MOcoefficients
            from pyscf.scf.chkfile import dump_scf
            from pyscf.tools import fcidump
            return run_pyscf_mocoeff, replace_pyscf_MOcoefficients, moldata, dump_scf, fcidump
        if package == 'psi4_genibo':
            from genibo_MOcoeff._run_psi4_mocoeff import run_psi4_mocoeff, replace_psi4_MOcoefficients
            from psi4.driver import fcidump
            return run_psi4_mocoeff, replace_psi4_MOcoefficients, moldata, fcidump
        if package == 'dirac_genibo':
            from genibo_MOcoeff import run_dirac_mocoef
            return run_dirac_mocoef, moldata
        if package == 'adf_genibo':
            from genibo_MOcoeff import generate_adf_input
            return moldata, generate_adf_input
    except ImportError as error:
        raise ImportError(
            "This feature requires the installation of genibo_MOcoeff. "
            "It can be installed with: pip install python_scripts/genibo."
        ) from error


class ROSE(ROSEInputDataClass, FileIOCalculator):
    """A general ASE calculator for ROSE (Reduction of Orbital Space Extent).

    Args:
        ROSEInputDataClass (ROSEInputDataClass): Custom-made base dataclass
            for ROSE; see rose_dataclass.py.
        FileIOCalculator (FileIOCalculator): Base class for ASE calculators
            that write/read input/output files.
    """
    implemented_properties: List[str] = ['energy']
    command: str = ''
    restart: bool = None
    ignore_bad_restart_file: bool = FileIOCalculator._deprecated
    label: str = 'ROSE'
    atoms: Atoms = None
    _directory = '.'
    profile = OldShellProfile(command)

    def __init__(self, *args, **kwargs) -> None:
        """ASE-ROSE Class Constructor to initialize the object.

        Give an example here on how to use....
        """
        # initializing base classes
        ROSEInputDataClass.__init__(self, *args, **kwargs)
        FileIOCalculator.__init__(
            self,
            restart=self.restart,
            ignore_bad_restart_file=self.ignore_bad_restart_file,
            label=self.label,
            command=self.command,
            profile=self.profile,
            *args,
            **kwargs,
        )

        # define specific attributes to be used later
        self.mol_frags_filenames: List[str] = []
        self.rose_output_filename = "OUTPUT_ROSE"
        self.original_wfn = None
        self.mol_target = None

    def calculate(self, *args, **kwargs) -> None:
        """Execute ROSE workflow."""
        # FileIOCalculator.calculate(self, *args, **kwargs)
        self.get_input_genibo()

        self.get_rose_inputs()
        self.run_rose()

        if self.save_data:
            self.save_ibos()

    def get_input_genibo(self) -> None:
        """Generate fortran genibo file for ROSE.

        This method generates a file named "INPUT_GENIBO"
            containing ROSE input options.
        """
        input_genibo_filename = "INPUT_GENIBO"
        write_rose_genibo_in(
            input_genibo_filename,
            self
        )

    def get_rose_inputs(self) -> None:
        """Generate molecule/frags xyz's and initial MO coeffs"""
        call_selected_mo_calculator = {
            ROSEMOCalculators.PYSCF.value: self._get_pyscf_genibo_mo_coeff,
            ROSEMOCalculators.PSI4.value: self._get_psi4_genibo_mocoeff,
            ROSEMOCalculators.DIRAC.value: self._get_dirac_genibo_mocoeff,
            ROSEMOCalculators.ADF.value: self._get_adf_genibo_mocoeff
        }
        # run calculator defined in "mo_calculator" attribute.
        call_selected_mo_calculator[self.rose_mo_calculator]()

    def run_rose(self) -> None:
        """Run ROSE executable 'rose.x'."""
        self.profile.command = "rose.x > " + self.rose_output_filename
        FileIOCalculator.calculate(self)

    def read_results(self) -> None:
        """Read HF energy (IBO/MO) from ROSE output file."""
        energy = read_rose_out(self.rose_output_filename)
        if energy['energy'] is None:
            energy = {'energy': 0.0E+00}
        self.results = energy

    def save_ibos(self) -> None:
        """Generates a checkpoint file ('ibo.chk') with the final IBOs."""
        save_ibos_with_selected_calculator = {
            ROSEMOCalculators.PYSCF.value: self._save_pyscf_ibo,
            ROSEMOCalculators.PSI4.value: self._save_psi4_ibo,
            ROSEMOCalculators.DIRAC.value: self._save_dirac_ibo,
            ROSEMOCalculators.ADF.value: self._save_adf_ibo
        }
        save_ibos_with_selected_calculator[self.rose_mo_calculator]()

    def _get_pyscf_genibo_mo_coeff(self) -> None:
        """Helper function that calls genibo_MOcoeff pyscf."""
        # import module
        (run_pyscf_mocoeff, replace_pyscf_MOcoefficients,
         moldata, _, _) = _import_genibo('pyscf_genibo')

        # defining general calculation attributes.
        geometry = self.rose_target.atoms
        data_directory = os.getcwd()
        multiplicity = self.rose_target.multiplicity
        charge = self.rose_target.charge
        basis1 = self.rose_target.basis
        uncontract = self.rose_target.uncontract
        relativistic = self.relativistic
        spherical = self.spherical
        restricted = self.restricted
        openshell = self.openshell
        save = self.save_data
        fragments = self.rose_frags

        # generating xyz for the target molecular system.
        molecule = moldata(
            geometry=geometry,
            multiplicity=multiplicity,
            charge=charge,
            data_directory=data_directory
        )

        # calculating MOs
        molecule, mol, original_wfn = run_pyscf_mocoeff(
            molecule,
            basis=basis1,
            relativistic=relativistic,
            spherical=spherical,
            restricted=restricted,
            openshell=openshell,
            save=save,
            uncontract=uncontract
        )

        # renaming final xyz and orbital files to be consistent with ROSE.
        subprocess.check_call(
            f"mv -f {molecule.name}.h5 CHECKPOINT.h5",
            shell=True,
            cwd=data_directory
        )

        subprocess.check_call(
            f"mv -f {molecule.name}.xyz MOLECULE.XYZ",
            shell=True,
            cwd=data_directory
        )

        # saving final SCF wave function.
        if self.save_data:
            subprocess.check_call(
                f"mv -f {molecule.name}.chk MOLECULE.chk",
                shell=True,
                cwd=data_directory
            )

        # this will be latter used to save the final ibos.
        self.original_wfn = original_wfn
        self.mol_target = mol

        # Fragment calculation.
        for i, fragment_data in enumerate(fragments):
            frag_geo = fragment_data.atoms
            multiplicity_frag = fragment_data.multiplicity
            charge_frag = fragment_data.charge
            basis2 = fragment_data.basis

            # Generate xyz for the fragments.
            fragment = moldata(
                geometry=frag_geo,
                multiplicity=multiplicity_frag,
                charge=charge_frag,
                data_directory=data_directory
            )

            # Calculate atom-free AOs or free molecular fragment MOs.
            fragment, _, _ = run_pyscf_mocoeff(
                fragment,
                basis=basis2,
                relativistic=relativistic,
                spherical=spherical,
                openshell=openshell,
                restricted=restricted,
                uncontract=uncontract
            )

            # Generate ROSE fragment input files.
            frag_filename = self._get_fragment_filename(fragment, i)

            subprocess.check_call(
                f"mv -f {fragment.name}.h5 {frag_filename}.h5",
                shell=True,
                cwd=data_directory
            )

            subprocess.check_call(
                f"mv -f {fragment.name}.xyz {frag_filename}.xyz",
                shell=True,
                cwd=data_directory
            )

    def _save_pyscf_ibo(
            self,
            input_file: str = "ribo.h5",
            output_chk_file: str = "ibo.chk",
            output_fcidump_file: str = "ibo.fcidump"
    ) -> None:
        """Helper function to save ibos from pyscf."""
        # import module
        (run_pyscf_mocoeff, replace_pyscf_MOcoefficients,
         moldata, _, fcidump) = _import_genibo('pyscf_genibo')
        original_wfn = self.original_wfn
        mol = self.mol_target
        nmo, ibo_wfn = replace_pyscf_MOcoefficients(original_wfn,mol,'ribo')
        if self.fcidump:
            fcidump.from_scf(ibo_wfn, output_fcidump_file)

    def _get_psi4_genibo_mocoeff(self) -> None:
        """Helper function that calls genibo_MOcoeff psi4."""
        # import module
        (run_psi4_mocoeff, replace_psi4_MOcoefficients, moldata, _) = _import_genibo('psi4_genibo')

        # defining general calculation attributes.
        geometry = self.rose_target.atoms
        data_directory = os.getcwd()
        multiplicity = self.rose_target.multiplicity
        charge = self.rose_target.charge
        basis1 = self.rose_target.basis
        uncontract = self.rose_target.uncontract
        relativistic = self.relativistic
        spherical = self.spherical
        restricted = self.restricted
        openshell = self.openshell
        save = self.save_data
        fragments = self.rose_frags

        molecule = moldata(
            geometry=geometry,
            multiplicity=multiplicity,
            charge=charge,
            data_directory=data_directory
        )

        # calculating MOs
        molecule, scf_wfn = run_psi4_mocoeff(
            molecule,
            basis=basis1,
            relativistic=relativistic,
            spherical=spherical,
            restricted=restricted,
            openshell=openshell,
            fcidump=False,
            save=save,
            uncontract=uncontract
        )
        # this will be latter used to save the final ibos.
        self.original_wfn = scf_wfn

        # renaming final xyz and orbital files to be consistent with ROSE.
        subprocess.check_call(
            f"mv -f {molecule.name}.h5 CHECKPOINT.h5",
            shell=True,
            cwd=data_directory
        )

        subprocess.check_call(
            f"mv -f {molecule.name}.xyz MOLECULE.XYZ",
            shell=True,
            cwd=data_directory
        )

        # saving final SCF wave function.
        if self.save_data:
            subprocess.check_call(
                f"mv -f {molecule.name}.npy MOLECULE.npy",
                shell=True,
                cwd=data_directory
            )

        # Fragment calculation.
        for i, fragment_data in enumerate(fragments):
            frag_geo = fragment_data.atoms
            multiplicity_frag = fragment_data.multiplicity
            charge_frag = fragment_data.charge
            basis2 = fragment_data.basis

            # Generate xyz for the fragments.
            fragment = moldata(
                geometry=frag_geo,
                multiplicity=multiplicity_frag,
                charge=charge_frag,
                data_directory=data_directory
            )

            # Calculate atom-free AOs or free molecular fragment MOs.
            fragment, _ = run_psi4_mocoeff(
                fragment,
                basis=basis2,
                relativistic=relativistic,
                spherical=spherical,
                restricted=restricted,
                openshell=openshell,
                fcidump=False,
                save=False,
                uncontract=uncontract
            )

            # Generate ROSE fragment input files.
            frag_filename = self._get_fragment_filename(fragment, i)

            subprocess.check_call(
                f"mv -f {fragment.name}.h5 {frag_filename}.h5",
                shell=True,
                cwd=data_directory
            )

            subprocess.check_call(
                f"mv -f {fragment.name}.xyz {frag_filename}.xyz",
                shell=True,
                cwd=data_directory
            )

    def _save_psi4_ibo(
            self,
            new_basis: str = "ribo"
    ) -> None:
        """Helper function to save ibos from psi4."""
        # import module
        (_, replace_psi4_MOcoefficients, _, fcidump) = _import_genibo('psi4_genibo')

        scf_wfn = self.original_wfn
        nmo = replace_psi4_MOcoefficients(scf_wfn,new_basis)
        
        scf_wfn.to_file(new_basis+'.npy')

        if self.fcidump:
            fcidump(scf_wfn, fname=new_basis+'.fcidump')

    def _get_dirac_genibo_mocoeff(self) -> None:
        """Helper function that calls genibo_MOcoeff DIRAC."""
        raise NotImplementedError(
            "DIRAC not yet implemented as ROSE-ASE MO calculator."
        )

    def _save_dirac_ibo(self) -> None:
        """Helper function to save ibos from DIRAC."""
        raise NotImplementedError(
            "DIRAC not yet implemented as ROSE-ASE MO calculator."
        )

    def _get_adf_genibo_mocoeff(self) -> None:
        """Helper function that calls genibo_MOcoeff ADF."""
        raise NotImplementedError(
            "ADF not yet implemented as ROSE-ASE MO calculator."
        )

    def _save_adf_ibo(self) -> None:
        """Helper function to save ibos from ADF."""
        raise NotImplementedError(
            "ADF not yet implemented as ROSE-ASE MO calculator."
        )

    def _get_fragment_filename(self, fragment, frag_idx) -> List[str]:
        """Generates fragment filename."""
        if self.rose_calc_type == ROSECalcType.ATOM_FRAG.value:
            # use the Z number of the atom for filename
            frag_z = fragment.protons[0]
            frag_filename = f"{frag_z:03d}"
        else:
            # use the index of the fragment for filename
            frag_filename = f"frag{frag_idx+1}"
        return frag_filename
