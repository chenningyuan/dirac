"""Input/Output functions for ROSE."""
import re
from typing import Optional, List, Any, Union, TextIO, Dict

from ase.units import Ha

from .rose_dataclass import ROSEInputDataClass
from .rose_dataclass import ROSECalcType, ROSEIFOVersion, ROSEILMOExponent


def _write_line(file: TextIO, line: str) -> None:
    """
    Write a line of text to a file.

    Args:
        file (file object): The file object to write to.
        line (str): The line of text to write.
    """
    file.write(line + "\n")


def _write_section(
        file: TextIO, section_name: str, value: Optional[Any] = None
) -> None:
    """
    Write a section name and an optional value to a file.

    Args:
        file (file object): The file object to write to.
        section_name (str): The name of the section.
        value (optional[any]): The value associated with the
            section (default: None).
    """
    _write_line(file, f".{section_name.upper()}")
    if value is not None:
        _write_line(file, str(value))


def _write_key_value_pairs(
        file: TextIO, section_name: str, items: List[List[int]]
) -> None:
    """
    Write a section name and a list of key-value pairs to a file.

    Args:
        file (file object): The file object to write to.
        section_name (str): The name of the section.
        items (list[any]): A list of key-value pairs.
    """
    for item in items:
        _write_section(file, section_name)
        _write_line(file, str(item[0]))
        _write_line(file, str(item[1]))


def write_rose_genibo_in(filename: str, inp_data: ROSEInputDataClass) -> None:
    """
    Write ROSE genibo input datafile.

    Args:
        filename (str): The name of the file to write the input data to.
        inp_data (ROSEInputDataClass): The input data to write.
        **params: Additional parameters.

    Returns:
        None
    """
    with open(filename, "w", encoding='utf-8') as file:
        _write_line(file, "**ROSE")
        if inp_data.version != ROSEIFOVersion.STNDRD_2013.value:
            _write_section(file, "VERSION", inp_data.version)
        _write_section(file, "CHARGE", inp_data.rose_target.charge)
        if inp_data.exponent != ROSEILMOExponent.TWO.value:
            _write_section(file, "EXPONENT", inp_data.exponent)
        if inp_data.rose_mo_calculator == 'psi4' or inp_data.rose_mo_calculator == 'pyscf':
           file_format = 'h5'
        else:
           file_format = inp_data.rose_mo_calculator
        _write_section(file, "FILE_FORMAT", file_format)
        if inp_data.test:
            _write_section(file, "TEST")
        if not inp_data.restricted:
            _write_section(file, "UNRESTRICTED")
        if not inp_data.spatial_orbitals:
            _write_section(file, "SPINORS")
        if inp_data.include_core:
            _write_section(file, "INCLUDE_CORE")
        if inp_data.rose_calc_type == ROSECalcType.MOL_FRAG.value:
            _write_section(file, "NFRAGMENTS", len(inp_data.rose_frags))
            if inp_data.additional_virtuals_cutoff:
                _write_section(
                    file,
                    "ADDITIONAL_VIRTUALS_CUTOFF",
                    inp_data.additional_virtuals_cutoff
                )
            if inp_data.frag_threshold:
                _write_section(file, "FRAG_THRESHOLD", inp_data.frag_threshold)
            if inp_data.frag_valence:
                _write_key_value_pairs(
                    file,
                    "FRAG_VALENCE",
                    inp_data.frag_valence
                )
            if inp_data.frag_core:
                _write_key_value_pairs(file, "FRAG_CORE", inp_data.frag_core)
            if inp_data.frag_bias:
                _write_key_value_pairs(file, "FRAG_BIAS", inp_data.frag_bias)
        if len(inp_data.avas_frag) != 0:
            avas_threshold = '1.e-14' # hardwired to the default
            _write_line(file, ".AVAS_THRESHOLD")
            _write_line(file, avas_threshold)
            for frag, mos in zip(inp_data.avas_frag,inp_data.nmo_avas):
                 _write_line(file,".FRAG_AVAS")
                 _write_line(file,str(frag))
                 _write_line(file,str(len(mos)))
                 for mo in mos:
                     _write_line(file, str(mo))
        _write_line(file, "*END OF INPUT")


_rose_scf_energy_re = re.compile(r"HF energy \(IBO/MO\)\s+:\s+([-\d.]+)")


def read_rose_out(output_filename: str) -> Dict[str, Any]:
    """Reads ROSE output data from a file.

    Args:
        output_filename: The ROSE output file.

    Returns:
        A dictionary containing the extracted data from the output.
    """
    energy = None

    with open(output_filename, 'r', encoding='utf-8') as file:
        content = file.read()  # Read the entire content of the file

    ehfmatch = re.search(_rose_scf_energy_re, content)
    if ehfmatch is not None:
        energy = float(ehfmatch.group(1)) * Ha

    output = {'energy': energy}
    return output


def write_int(file: TextIO, text: str, var: int) -> None:
    """Write an integer value to a file in a specific format.

    Args:
        f (file object): The file object to write to.
        text (str): A string of text to precede the integer value.
        var (int): The integer value to be written to the file.

    Returns:
        None
    """
    file.write("{:43}I{:17d}\n".format(text, var))


def write_int_list(file: TextIO, text: str, var: List[int]) -> None:
    """Write a list of integers to a file in a specific format.

    Args:
        f (file object): The file object to write to.
        text (str): A string of text to precede the list of integers.
        var (list): The list of integers to be written to the file.

    Returns:
        None
    """
    file.write("{:43}{:3} N={:12d}\n".format(text, "I", len(var)))
    dim = 0
    buff = 6
    if len(var) < 6:
        buff = len(var)
    for i in range((len(var)-1)//6+1):
        for j in range(buff):
            file.write("{:12d}".format(var[dim+j]))
        file.write("\n")
        dim = dim + 6
        if (len(var) - dim) < 6:
            buff = len(var) - dim


def write_singlep_list(file: TextIO, text: str, var: List[float]) -> None:
    """Write a list of single-precision floating-point values to a file object.

    Args:
        file (file object): A file object to write to.
        text (str): The text to be written before the list.
        var (list): The list of single-precision floating-point
            values to write.

    Returns:
        None
    """
    file.write("{:43}{:3} N={:12d}\n".format(text, "R", len(var)))
    dim = 0
    buff = 5
    if len(var) < 5:
        buff = len(var)
    for i in range((len(var)-1)//5+1):
        for j in range(buff):
            file.write("{:16.8e}".format(var[dim+j]))
        file.write("\n")
        dim = dim + 5
        if (len(var) - dim) < 5:
            buff = len(var) - dim


def write_doublep_list(file: TextIO, text: str, var: List[float]) -> None:
    """Write a list of double precision floating point numbers to a file.

    Args:
        f (file object): the file to write the data to
        text (str): a label or description for the data
        var (list): a list of double precision floating point
            numbers to write to file

    Returns:
        None
    """
    file.write("{:43}{:3} N={:12d}\n".format(text, "R", len(var)))
    dim = 0
    buff = 5
    if len(var) < 5:
        buff = len(var)
    for _ in range((len(var)-1)//5+1):
        for j in range(buff):
            file.write("{:24.16e}".format(var[dim+j]))
        file.write("\n")
        dim = dim + 5
        if (len(var) - dim) < 5:
            buff = len(var) - dim


def read_int(text: str, file: TextIO) -> Union[int, None]:
    """Read an integer value from a text file.

    Args:
        text (str): The text to search for in the file.
        f (file object): The file object to read from.

    Returns:
        int: The integer value found in the file.
    """
    for line in file:
        if re.search(text, line):
            var = int(line.rsplit(None, 1)[-1])
            return var
    return None


def read_real(text: str, file: TextIO) -> Union[float, None]:
    """Read a floating-point value from a text file.

    Args:
        text (str): The text to search for in the file.
        f (file): The file object to read from.

    Returns:
        float: The floating-point value found in the file.
    """
    for line in file:
        if re.search(text, line):
            var = float(line.rsplit(None, 1)[-1])
            return var
    return None


def read_real_list(text: str, file: TextIO) -> Union[List[float], None]:
    """Read a list of floating-point values from a text file.

    Args:
        text (str): The text to search for in the file.
        file (file object): The file object to read from.

    Returns:
        list: A list of floating-point values found in the file.
    """
    for line in file:
        if re.search(text, line):
            n = int(line.rsplit(None, 1)[-1])
            var = []
            for _ in range((n-1)//5+1):
                line = next(file)
                for j in line.split():
                    var += [float(j)]
            return var
    return None
