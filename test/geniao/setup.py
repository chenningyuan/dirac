import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
      
setuptools.setup(
    name="geniao_dfcoef_generator",
    version="0.0.1",
    author="bsenjean",
    author_email="bsenjean@gmail.com",
    description='Generate all the DFCOEF files for the construction of Intrinsic Fragment Orbitals',
    long_description=long_description,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
