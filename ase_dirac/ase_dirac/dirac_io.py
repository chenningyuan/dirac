"""Input/Output functions for Rose."""

from typing import Dict, Any, List
from copy import deepcopy
import re
from collections import OrderedDict

from ase.calculators.calculator import InputError
from ase.units import Ha


def _is_end_marker(line: str) -> bool:
    """Return True for DIRAC sentinel lines that should not become keys."""
    normalized = line.strip().upper()
    return normalized in {"*END OF", "*END OF INPUT"}


def _update_dict(
        dictionary: Dict[str, Any],
        key: str,
        value: Any
) -> Dict[str, Any]:
    """Updates a dict with a new key-value pair and put it at first position.

    Args:
        dictionary: The original dictionary to be updated.
        key: The key of the new element to be added.
        value: The value of the new element to be added.

    Returns:
        The updated dictionary with the new key-value pair at first position.
    """
    ordered_dict = OrderedDict(dictionary)
    ordered_dict[key] = value
    updated_dict = OrderedDict(
        [(key, ordered_dict[key]) for key in reversed(ordered_dict)]
    )
    return updated_dict


def _replace_underscores(dictionary: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively replaces underscores in dict keys and values with spaces.

    Args:
        dictionary: The dictionary to be processed.

    Returns:
        The dictionary with underscores replaced by spaces.
    """
    new_dict = {}
    for key, value in dictionary.items():
        if key.upper() not in ['VISUAL', '**VISUAL']:
            if isinstance(value, dict):
                new_dict[key.replace("_", " ")] = _replace_underscores(value)
            else:
                new_dict[key.replace("_", " ")] = value.replace("_", " ").strip()
        else:
            new_dict[key] = value
    return new_dict


def _format_value(arg: Any) -> str:
    """Formats the input value.

    Args:
        arg: The value to be formatted.

    Returns:
        The formatted value as a string.
    """
    if isinstance(arg, dict):
        formatted_values = []
        for key, val in arg.items():
            formatted_values.append(
                '{}\n{}'.format(str(key).upper(), str(val))
            )
        format_str = '\n'.join(formatted_values)
    elif isinstance(arg, (float, int, str)):
        if not arg:
            format_str = '{}'.format(str(arg))
        else:
            format_str = '\n{}'.format(str(arg))
    else:
        raise InputError('Format for', arg, 'not allowed.')
    return format_str


def _write_block(name: str, args: Dict[str, Any]) -> str:
    """Writes a block of formatted data.

    Args:
        name: The name of the block.
        args: The dictionary containing the data for the block.

    Returns:
        The formatted block as a string.
    """
    out = ['**{}'.format(name.upper())]
    for key, value in args.items():
        if key.startswith('.'):
            out.append('{}{}'.format(key.upper(), _format_value(value)))
        elif key.startswith('*'):
            out.append('{}\n{}'.format(key.upper(), _format_value(value)))
        else:
            raise InputError(
                'Check input! Options and subsections '
                'must start with . and *, respectively'
            )
    return "\n".join(out)


def _parse_dirac_inp(inp_str: str) -> Dict:
    """Internal helper to parse raw DIRAC text into a nested dictionary.
    
    Parameters
    ----------
    inp_str : str
        The raw DIRAC input as a string.
    Returns
    -------
    Dict
        A nested dictionary representing the structured content of the DIRAC input.
    """
    accepted_keywords: List[str] = [
        "dirac",
        "general",
        "molecule",
        "hamiltonian",
        "wave_function",
        "analyze",
        "properties",
        "visual",
        "integrals",
        "grid",
        "moltra",
    ]

    # Clean input and add empty lines where required
    inp_list = [
        line.strip()
        for line in inp_str.split("\n")
        if line.strip() != "" and not line.strip().startswith("!")
    ]

    to_insert = []
    for index in range(len(inp_list) - 1):
        line = inp_list[index]
        if line[0] == ".":
            next_line = inp_list[index + 1]
            if next_line[:2] == "**" or next_line[:1] in ["*", "."]:
                to_insert.append(index + 1)

    for index in reversed(to_insert):
        inp_list.insert(index, "")

    inp_list.append("")

    # Locate the different input sections
    double_stars = []
    single_stars = []
    points = []

    for index, line in enumerate(inp_list):
        if line[:2] == "**":
            double_stars.append(index)
        elif line[:1] == "*":
            single_stars.append(index)
        elif line[:1] == ".":
            points.append(index)

    # Extract the ** sections
    if len(double_stars) == 1:
        double_stars_slices = [slice(double_stars[0], len(inp_list))]
    elif len(double_stars) > 1:
        double_stars_slices = [
            slice(double_stars[i], double_stars[i + 1])
            for i in range(len(double_stars) - 1)
        ]
        double_stars_slices.append(slice(double_stars_slices[-1].stop, len(inp_list)))
    else:
        return {}

    # initialize dictionary
    input_dict = {}
    prev_double = None

    # Indicate type of each input line
    keyword_type = ["R" for _ in range(len(inp_list))]
    for index in double_stars:
        keyword_type[index] = "D"
    for index in single_stars:
        keyword_type[index] = "S"
    for index in points:
        keyword_type[index] = "P"

    # Parse by ** section
    for dslice in double_stars_slices:
        slice_list = inp_list[dslice]
        type_list = keyword_type[dslice]

        # read ** section keyword and initialize dictionary
        inp_keyword = slice_list[0].lower()
        inp_keyword = (
            "**wave function" if inp_keyword == "**wave functions" else inp_keyword
        )  # due to the keyword in the calculator
        dict_keyword = inp_keyword[2:].replace(" ", "_")

        # This overcomes cases where hierarchy is  ** -> ** -> . like **WAVE FUNCTION -> **RELCCSD -> .ENERGY
        if dict_keyword in accepted_keywords:
            input_dict.update({dict_keyword: {}})
            prev_double = dict_keyword
        else:
            type_list[0] = "S"

        # parent is defined to keep track of the ** -> * -> . hierarchy
        parent = None

        # Examine the hierarchy of the content of each ** slice
        for index, ttype in enumerate(type_list):
            # Define the ** keyword as parent of *
            if ttype == "S":
                if _is_end_marker(slice_list[index]):
                    parent = None
                    continue
                parent = slice_list[index]
                input_dict[prev_double].update({parent: {}})

            # The logic is that after any P, there must be at least an R due to construction
            if ttype == "P":
                start_p_index = index
                end_p_index = len(type_list)  # Assume longest possible slice

                # Check if the actual slice is smaller
                for inndex, tttype in enumerate(type_list[index + 1 :]):
                    if tttype != "R":
                        end_p_index = start_p_index + inndex + 1
                        break

                # once found, define slice, and content of P
                p_slice = slice(start_p_index, end_p_index)
                val_str = "\n".join(slice_list[p_slice][1:])

                # Use parent to determine if the P dictionary should be added as ** -> . directly or as ** -> * -> .
                # Within this logic, it is implicit that any ** that is not in the accepted keywords of the DIRAC class
                # will behave as ** -> ** -> . instead of trying to access a dictionary key that is not atmitted
                if parent is None:
                    input_dict[prev_double].update({slice_list[p_slice][0]: val_str})
                else:
                    input_dict[prev_double][parent].update(
                        {slice_list[p_slice][0]: val_str}
                    )

    return input_dict


def write_dirac_in(input_filename: str, **params: Dict[str, Any]) -> None:
    """Writes DIRAC input data to a file.

    Args:
        input_filename: input file to be generated.
        params: Additional parameters for the DIRAC input.

    Returns:
        None.
    """
    tmp_params = deepcopy(params)

    # replace any underscore by a white space
    # important for generating DIRAC input
    # Except VISUAL section (id grid != id_grid)
    tmp_params = _replace_underscores(tmp_params)

    out = []
    out += [_write_block(*item) for item in tmp_params.items()]
    out.append('*END OF INPUT')


    with open(input_filename, 'w') as file:
        to_write = '\n'.join(out)
        to_write = '\n'.join(line for line in to_write.splitlines() if line.strip()) + '\n'
        if 'id grid' in to_write:
            to_write.replace('id grid', 'id_grid')
        file.write(to_write)


_scf_energy_re = re.compile(r"Total energy\s+:\s+([-\d.]+)")
_mp_energy_re  = re.compile(r"Total MP2 energy\s+:\s+([-\d.]+)")
_cc_energy_re  = re.compile(r"Total CCSD energy\s+:\s+([-\d.]+)")
_cct_energy_re = re.compile(r"Total\s+CCSD\(T\)\s+energy\s*:\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
_ci_energy_re  = re.compile(r"CI Total Energy\s+:\s+([-\d.]+)")

def read_dirac_out(output_filename: str) -> Dict[str, Any]:
    """Reads DIRAC output data from a file.

    Args:
        output_filename: The DIRAC output file.

    Returns:
        A dictionary containing the extracted data from the DIRAC output.
    """
    output = {}
    energy = None

    with open(output_filename, 'r') as file:
        content = file.read()  # Read the entire content of the file

    ehfmatch = re.search(_scf_energy_re, content)
    if ehfmatch is not None:
        energy = float(ehfmatch.group(1)) * Ha
        output.update({'SCF_energy':energy})

    empmatch = re.search(_mp_energy_re, content)
    if empmatch is not None:
        energy = float(empmatch.group(1)) * Ha
        output.update({'MP2_energy':energy})

    eccmatch = re.search(_cc_energy_re, content)
    if eccmatch is not None:
        energy = float(eccmatch.group(1)) * Ha
        output.update({'CCSD_energy':energy})

    ecctmatch = re.search(_cct_energy_re, content)
    if ecctmatch is not None:
        energy = float(ecctmatch.group(1)) * Ha
        output.update({'CCSDT_energy':energy})

    ecimatch = re.search(_ci_energy_re, content)
    if ecimatch is not None:
        energy = float(ecimatch.group(1)) * Ha
        output.update({'CI_energy':energy})

    output.update({'energy': energy})

    return output
