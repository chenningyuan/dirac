"""Functions to prepare ADF input"""
from __future__ import absolute_import
import os
from ._utils import *

def rename(molecule,frags_list):
    xyz_file = molecule.filename + '.xyz'
    xyz = molecule.data_directory + "/MOLECULE.xyz"
    os.rename(xyz_file,xyz)
    for i in range(len(frags_list)):
       xyz_file = frags_list[i].filename + '.xyz'
       xyz = frags_list[i].data_directory + "/frag{:d}.xyz".format(i+1)
       os.rename(xyz_file,xyz)

def generate_adf_input(molecule,
                       basis,
                       frags_list,
                       frags_basis,
                       restricted,
                       openshell,
                       restart):

    # Creation of the Molecule XYZ files
    molecule_string = create_geometry_string(molecule.geometry)
    with open(molecule.filename + ".xyz", 'w') as f:
       f.write(str(len(molecule.geometry))+'\n')
       f.write(molecule.filename + ' # anything can be in this line\n')
       f.write(molecule_string)

    # Creation of the fragment XYZ files and count number of MOs in the minimal basis per fragment
    nmo_frag={}
    nfragments = len(frags_list)
    i = 0
    for fragment in frags_list:
       frag_geo = fragment.geometry
       molecule_string, nmo_frag[i] = create_geometry_string(frag_geo,minbas=True)
       with open(fragment.filename + ".xyz", 'w') as f:
          f.write(str(len(frag_geo))+'\n')
          f.write(fragment.filename + ' # anything can be in this line\n')
          f.write(molecule_string)
       i += 1
    
    # rename molecule and fragments xyz files:
    rename(molecule,frags_list)
    
    # Creation of the input to be read by ADF.py
    with open(molecule.data_directory + "/INPUT_ADF","w") as f:
       f.write(molecule.data_directory + "\n")
       f.write(str(restricted) + "\n")
       f.write(str(openshell) + "\n")
       f.write(basis + "\n")
       f.write(frags_basis + "\n")
       f.write(str(molecule.charge) + "\n")
       f.write(str(molecule.multiplicity) + "\n")
       f.write(str(nfragments) + "\n")
       for i in range(nfragments):
          tmpstr = frags_list[i].data_directory + "/frag{:d}.xyz".format(i+1)
          f.write("{:100s}{:10d}{:10d}{:10d}\n".format(tmpstr, nmo_frag[i],frags_list[i].charge,\
                                                        frags_list[i].multiplicity))
       f.write(str(restart) + "\n")
