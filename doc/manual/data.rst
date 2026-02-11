:orphan:
 

Data Handling in DIRAC
======================

Since release DIRAC23 all data suitable for restart and analysis purposes is stored on the file CHECKPOINT.h5

This file is in hdf5 format allowing for easy export of the generated data to other formats. The structure is
given by the data schema of which some important elements are highlighted below. This is usually so-called compound data,
also called data types, consisting of multiple individual elements. A full description can be found in the file `utils/DIRACschema.txt`,
we will merely summarize the most useful data here. For users familiar with Jupyter notebooks we recommend having a look
at the directory `demo_notebooks` that contains the Jupyter Notebook `DiracDataTutorial.ipynb` that demonstrates 
how data can be imported in Python for further processing. Similarly you may consult the notebook `TrexioConverter.ipynb`
to see how data can be converted to the TREXIO format. For developers we furthermore recommend consulting the
section on the CHECKPOINT file in the developer manual.

input
-----

This section contains of two sets of data: the `molecule` data type and the `aobasis` data type. 
 * The `molecule` data type contains the topology of the molecule with coordinates in xyz format, nuclear symbols and charges.
 * The `aobasis` data type contains all information regarding the employed atomic orbital basis: location of the expansion centers, exponential parameters, contraction coefficients, etc.

result
------

The amount of data in this section depends on the calculation that is carried out. Of interest for analysis is
typically the data type `mobasis` that contains the molecular orbitals and their energies. This data type is found for instance in
the `result/wavefunctions/scf` subsection that also contains the total energy of the SCF wavefunction.

Not all data is yet written to the CHECKPOINT file, for energies of wave functions other than SCF or MP2 it is still necessary to
consult the standard output file. In future releases the amount of data stored on CHECKPOINT.h5 will likely increase, but the aim
is to keep this file concise, large datasets like transformed 2-electron integrals will still be stored in a different file.

