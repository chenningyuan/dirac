"""
Setup file for DIRAC ASE calculator.
"""
import shutil

import setuptools


def validate_pam_available() -> None:
    """Abort installation when the DIRAC pam executable is not on PATH."""
    if shutil.which("pam") is None:
        raise RuntimeError(
            "The DIRAC executable 'pam' was not found on PATH.\n\n"
            "The ase_dirac package requires 'pam' to be available at "
            "installation time.\n"
            "Please add the DIRAC build directory containing 'pam' to PATH "
            "and rerun the installation.\n\n"
            "Example:\n"
            "  export PATH=/path/to/dirac/build:$PATH\n"
            "  pip install -e ."
        )


validate_pam_available()

with open("README.md", "r", encoding='utf-8') as file:
    long_description = file.read()
install_requires = [
    'ase',
    'h5py',
    'pytest',
]
setuptools.setup(
    name='ase_dirac',
    version='1.0.0',
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    author='Carlos M. R. Rocha',
    author_email='c.rocha@esciencecenter.nl',
    license="Apache-2.0",
    description='DIRAC ASE calculator',
    long_description=long_description,
    long_description_content_type='',
    url='',
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering",
    ],
    python_requires='>=3.8',
)
