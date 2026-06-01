# This script reads a .fchk file and truncate it to keep only the important informations, that are then read in ROSE (using FORTRAN).

import sys
import re

filename = sys.argv[1]

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
def write_real_list(f,text,var):
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

with open(filename,"r") as f:
   head              = f.readline()
   natom             = read_int("Number of atoms", f)
   charge            = read_int("Charge", f)
   multiplicity      = read_int("Multiplicity", f)
   nelec             = read_int("Number of electrons", f)
   nalpha            = read_int("Number of alpha electrons", f)
   nbeta             = read_int("Number of beta electrons", f)
   nbasis            = read_int("Number of basis functions", f)
   atom_numbers      = read_int_list("Atomic numbers", f)
   nuc_charges       = read_real_list("Nuclear charges", f)
   cart_coord        = read_real_list("Current cartesian coordinates", f)
   nshell            = read_int("Number of primitive shells", f)
   pureamd           = read_int("Pure/Cartesian d shells", f)
   pureamf           = read_int("Pure/Cartesian f shells", f)
   orb_momentum      = read_int_list("Shell types", f)
   nprim_shell       = read_int_list("Number of primitives per shell", f)
   shell_atom_map    = read_int_list("Shell to atom map", f)
   prim_exp          = read_real_list("Primitive exponents", f)
   contraction_coeff = read_real_list("Contraction coefficients", f)
   coord_shell       = read_real_list("Coordinates of each shell", f)
   total_energy      = read_real("Total Energy", f)
   alpha_energies    = read_real_list("Alpha Orbital Energies", f)
   alpha_MO          = read_real_list("Alpha MO coefficients", f)
   beta_energies     = read_real_list("Beta Orbital Energies", f)
   beta_MO           = read_real_list("Beta MO coefficients", f)
f.close()

with open(filename,"w") as f:
     f.write(head)
     write_int(f,"Number of atoms",natom)
     write_int(f,"Charge",charge)
     write_int(f,"Multiplicity",multiplicity)
     write_int(f,"Number of electrons",nalpha+nbeta)
     write_int(f,"Number of alpha electrons",nalpha)
     write_int(f,"Number of beta electrons",nbeta)
     write_int(f,"Number of basis functions",nbasis)
     write_int_list(f,"Atomic numbers",atom_numbers)
     write_real_list(f,"Nuclear charges",nuc_charges)
     write_real_list(f,"Current cartesian coordinates",cart_coord)
     write_int(f,"Number of primitive shells",nshell)
     write_int(f,"Pure/Cartesian d shells",pureamd)
     write_int(f,"Pure/Cartesian f shells",pureamf)
     write_int_list(f,"Shell types",orb_momentum)
     write_int_list(f,"Number of primitives per shell",nprim_shell)
     write_int_list(f,"Shell to atom map",shell_atom_map)
     write_real_list(f,"Primitive exponents",prim_exp)
     write_real_list(f,"Contraction coefficients",[1]*len(prim_exp))
     write_real_list(f,"Coordinates of each shell",coord_shell)
     f.write("{:43}R{:27.15e}\n".format("Total Energy",total_energy))
     write_real_list(f,"Alpha Orbital Energies",alpha_energies)
     write_real_list(f,"Alpha MO coefficients",alpha_MO)
     if (beta_energies != None):
        write_real_list(f,"Beta Orbital Energies",beta_energies)
        write_real_list(f,"Beta MO coefficients",beta_MO)
