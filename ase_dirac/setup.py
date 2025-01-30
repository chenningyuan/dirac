"""
Setup file for DIRAC ASE calculator.
"""
import setuptools
with open("README.md", "r", encoding='utf-8') as file:
    long_description = file.read()
install_requires = [
    'ase',
]
tests_require = [
    'pytest',
]
setuptools.setup(
    name='ase_dirac',
    version='1.0.0',
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    tests_require=tests_require,
    extras_require={
        'testing': tests_require,
    },
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