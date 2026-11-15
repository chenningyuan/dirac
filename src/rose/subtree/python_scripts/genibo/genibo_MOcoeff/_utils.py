from __future__ import absolute_import
import re

periodic_table = [
        '?',
        'H', 'He',
        'Li', 'Be',
        'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg',
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca',
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
        'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
        'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
        'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
periodic_hash_table = {}
for atomic_number, atom in enumerate(periodic_table):
    periodic_hash_table[atom] = atomic_number

# Ghost atoms (centres that carry basis functions but no nucleus and no electrons), counterpoise corrections.
# A ghost is written as "Gh(O)" (psi4 notation; "@O" is accepted as well).
_ghost_re = re.compile(r"^\s*(?:[Gg][Hh]\(\s*([A-Za-z]{1,2})\s*\)|@([A-Za-z]{1,2}))\s*$")

def split_atom_label(label):
    """Return (element symbol, is_ghost). Real atoms are returned unchanged, so that a
    mistyped symbol still gives the KeyError of periodic_hash_table where it is used."""
    m = _ghost_re.match(label)
    if m is None:
        return label, False
    return (m.group(1) or m.group(2)).capitalize(), True

def basis_charge(label):
    """Charge of the element that carries the basis functions (also for ghosts)."""
    return periodic_hash_table[split_atom_label(label)[0]]

def nuclear_charge(label):
    """Nuclear charge: 0 for a ghost atom, otherwise the atomic number."""
    element, ghost = split_atom_label(label)
    return 0 if ghost else periodic_hash_table[element]

def is_ghost(label):
    return split_atom_label(label)[1]

# Spin polarization of atoms on period table.
periodic_polarization = [-1,
                         1, 0,
                         1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 3, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0,
                         1, 0, 1, 2, 5, 6, 5, 8, 9, 0, 1, 0, 1, 2, 3, 2, 1, 0]


def write_int(f,text,var):
  f.write("{:43}I{:17d}\n".format(text,var))

def write_int_list(f,text,var):
  f.write("{:43}{:3} N={:12d}\n".format(text,"I",len(var)))
  dim = 0
  buff = 6
  if (len(var) < 6): buff = len(var)
  for i in range((len(var)-1)//6+1):
      for j in range(buff):
         f.write("{:12d}".format(var[dim+j]))
      f.write("\n")
      dim = dim + 6
      if (len(var) - dim) < 6 : buff = len(var) - dim

def write_singlep_list(f,text,var):
  f.write("{:43}{:3} N={:12d}\n".format(text,"R",len(var)))
  dim = 0
  buff = 5
  if (len(var) < 5): buff = len(var)
  for i in range((len(var)-1)//5+1):
      for j in range(buff):
         f.write("{:16.8e}".format(var[dim+j]))
      f.write("\n")
      dim = dim + 5
      if (len(var) - dim) < 5 : buff = len(var) - dim

def write_doublep_list(f,text,var):
  f.write("{:43}{:3} N={:12d}\n".format(text,"R",len(var)))
  dim = 0
  buff = 5
  if (len(var) < 5): buff = len(var)
  for i in range((len(var)-1)//5+1):
      for j in range(buff):
         f.write("{:24.16e}".format(var[dim+j]))
      f.write("\n")
      dim = dim + 5
      if (len(var) - dim) < 5 : buff = len(var) - dim

def read_int(text,f):
   for line in f:
      if re.search(text, line):
         var=int(line.rsplit(None, 1)[-1])
         return var

def read_real(text,f):
   for line in f:
      if re.search(text, line):
         var=float(line.rsplit(None, 1)[-1])
         return var

def read_int_list(text,f):
   for line in f:
      if re.search(text, line):
         n=int(line.rsplit(None, 1)[-1])
         var=[]
         for i in range((n-1)//6+1):
             line = next(f)
             for j in line.split():
                 var += [int(j)]
         return var

def read_real_list(text,f):
   for line in f:
      if re.search(text, line):
         n=int(line.rsplit(None, 1)[-1])
         var=[]
         for i in range((n-1)//5+1):
             line = next(f)
             for j in line.split():
                 var += [float(j)]
         return var

def create_geometry_string(geometry,minbas=False):
    """This function converts MolecularData geometry to psi4 geometry.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances in
            angstrom. Use atomic symbols to specify atoms.

    Returns:
        geo_string: A string giving the geometry for each atom on a line, e.g.:
            H 0. 0. 0.
            H 0. 0. 0.7414
    """

    geo_string = ''
    nmo_minbas = 0
    for item in geometry:
        atom = item[0]
        coordinates = item[1]
        line = '{} {} {} {}'.format(atom,
                                    coordinates[0],
                                    coordinates[1],
                                    coordinates[2])
        if len(geo_string) > 0:
            geo_string += '\n'
        geo_string += line

        nmo_minbas += atom_minbas_size(nuclear_charge(atom))   # a ghost has no minimal-basis orbitals
    if minbas: return geo_string, nmo_minbas
    if not minbas: return geo_string

# returns size of minimal basis (NB: counting spinors, so there is a factor of 2 wrt orbitals)
def atom_minbas_size(atom_number):
    z = int(atom_number)
    if (z==0): # 0 is to be used for point charges that do not carry a basis set
        ncore = 0
        nvale = 0
    elif z in range(1,2+1): # H, He
        ncore = 0
        nvale = 2
    elif z in range(3,4+1): # Li, Be
        ncore = 2
        nvale = 8 # valence shell ns np (8 spinors): the np shell is needed for the valence virtuals
    elif z in range(5,10+1): # B..Ne
        ncore = 2
        nvale = 8
    elif z in range(11,12+1): # Na,Mg
        ncore = 10
        nvale = 8 # 3s3p
    elif z in range(13,18+1): # Al..Ar
        ncore = 10
        nvale = 8
    elif z in range(19,20+1): # K,Ca
        ncore = 18
        nvale = 8 # 4s4p
    elif z in range(21,30+1): # Sc..Zn
        ncore = 18
        nvale = 12
    elif z in range(31,36+1): # Ga..Kr
        ncore = 18
        nvale = 18
    elif z in range(37,38+1): # Rb,Sr
        ncore = 36
        nvale = 8 # 5s5p
    elif z in range(39,48+1): # Y..Cd
        ncore = 36
        nvale = 12
    elif z in range(49,54+1): # In..Xe
        ncore = 36
        nvale = 18
    elif z in range(55,56+1): # Cs,Ba
        ncore = 54
        nvale = 8 # 6s6p
    elif z in range(57,70+1): # La..Yb
        ncore = 54 
        nvale = 16 # assuming that we will enforce 4f rather than 5d occupation in the atomic scf
    elif z in range(71,80+1): # Lu..Hg
        ncore = 54 
        nvale = 26
    elif z in range(81,86+1): # Tl..Tn
        ncore = 68
        nvale = 18 # place the 4f in the core
    elif z in range(87,88+1): # Fr,Ra
        ncore = 86
        nvale = 8 # 7s7p
    elif z in range(89,102+1): # Ac,RNo
        ncore = 86
        nvale = 16 # assuming that we will enforce 5f rather than 6d occupation in the atomic scf
    elif z in range(103,112+1): # Lr..Cn
        ncore = 86
        nvale = 26
    elif z in range(113,118+1): # Nh..Og
        ncore = 100 # place the 5f in the core
        nvale = 18
    elif z > 119: # Venturing into the unknown superheavies, impossible to decide the valence
        ncore = 118
        nvale = atom_number - ncore
    nmo = (ncore+nvale)//2
    return nmo
