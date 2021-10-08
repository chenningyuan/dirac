#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 09:03:56 2021

@author: Lucas Visscher

Read/write data from DIRAC labeled files

Examples of usage.

1) Read DFCOEF (unique format)
file_name = 'DFCOEF'
contents = read_DFCOEF(file_name)
print('Energy = ',contents['energy'])
print('NAO    = ',contents['nao'])

2) Read AOPROPER (or other files with this format)
contents = read_DFILE(file_name)
print("Overlap Matrix")
np.set_printoptions(precision=3,linewidth=80)
print(contents['OVERLAP '])

"""

import numpy as np
from scipy.io import FortranFile

###### Functions to read and write labels and data

def read_label(f):
    """
    Read a DIRAC style label from an open file
    """
    record = f.read_record(np.dtype('4a8'))
    stars=record[0].decode("utf-8")
    if stars != '********':
        raise IOError('Error reading label, file may be corrupted')
    date_run=record[1].decode("utf-8")
    time_run=record[2].decode("utf-8")
    label=record[3].decode("utf-8")
    
    return label,date_run,time_run

def write_data(f,labels,data):
    """
    Write labeled data to open DIRAC file
    """
    record = ['********'.encode("utf-8")]
    for label in labels:
        record.append(label.encode("utf-8"))
    f.write_record(record)
    f.write_record(data)
    return


def read_DFILE(file_name):
    """
    Open DIRAC type file and return dictionary of its contents
    """
    f = FortranFile(file_name, 'r')
    # Make a dictionary to store the contents of this fixed format file
    contents = {}
    label =''
    while label != 'EOF':
        label, date_run, time_run = read_label(f)
        if label == 'EOFLABEL':
            break
        contents[label] = f.read_reals()
        
    return contents

def write_DFILE(file_name,contents):
    """
    Open DIRAC type file and store the dictionary contents
    """
    f = FortranFile(file_name, 'w')
    # Make a dictionary to store the contents of this fixed format file
    label =''
    for label in contents:
        labels = ['Molecule','xyzcoord',label]
        write_data(f,labels,contents[label])
    write_data(f,['EOFLABEL','EOFLABEL','EOFLABEL'],0.0)
    f.close()

    return

def read_DFCOEF(file_name,dtype_int='<i4'):
    """
    Open DFCOEF-type file and return dictionary of its contents
    """
    f = FortranFile(file_name, 'r')
    # Make a dictionary to store the contents of this fixed format file
    contents = {}
    # First record contains information about run time, info is a control string that can be ignored.
    info, date_run, time_run = read_label(f)
    contents['date'] = date_run
    contents['time'] = time_run
    # Second record contains a.o. information about dimensionality.
    try:
       record = f.read_record('a74', dtype_int, dtype_int, '3'+dtype_int,'f8')
    except:
       f.close()
       f = FortranFile(file_name, 'r')
       info, date_run, time_run = read_label(f)
       record = f.read_record('a74', dtype_int, dtype_int, '6'+dtype_int,'f8')
    contents['title']  = record[0]
    contents['nfsym']  = record[1]
    contents['nz']     = record[2]
    if (record[1] == 1):
        contents['npo']    = record[3][0]
        contents['nmo']    = record[3][1]
        contents['nto']    = record[3][0]+record[3][1]
        contents['nao']    = record[3][2]
    else:
        contents['npo']    = (record[3][0],record[3][3])
        contents['nmo']    = (record[3][1],record[3][4])
        contents['nto']    = (record[3][0]+record[3][1],record[3][3]+record[3][4])
        contents['nao']    = (record[3][2],record[3][5])
    contents['energy'] = record[4][0]
    label, date_run, time_run = read_label(f)
    contents['mo_coeff'] = f.read_reals()
    label, date_run, time_run = read_label(f)
    contents['eigenval'] = f.read_reals()
    label, date_run, time_run = read_label(f)
    contents['ibeig'] = f.read_ints()
    label, date_run, time_run = read_label(f)
    if label == 'KAPPA   ':
       # special case for atomic calculations with an extra record to process
       contents['kappa'] = f.read_ints()
       label, date_run, time_run = read_label(f)
    # 6 records with basis set information
    ints = f.read_ints(dtype=dtype_int)
    ao_nshells = ints[1]
    contents['ao_nbas']    = ints[0]   # size of the basis set
    contents['ao_npriexp'] = ints[2]   # number of primitive functions
    contents['ao_ncnt']    = ints[3]   # maximum contraction length
    contents['ao_nshells'] = ao_nshells # number of distint ao shells, then for each shell:
    ints = f.read_ints(dtype=dtype_int).reshape((3,ao_nshells))
    contents['ao_orbmom']  = ints[0,:] # orbital momentum (s=1,p=2, d=3, etc.)
    contents['ao_idcent']  = ints[1,:] # id of its center
    contents['ao_nprim']   = ints[2,:] # number of primitives
    ints = f.read_ints(dtype=dtype_int).reshape((2,ao_nshells))
    contents['ao_jstrt']   = ints[0,:] # start address of contraction coefficients
    contents['ao_numcf']   = ints[1,:] # column number of contraction coefficients
    contents['ao_cent']    = f.read_reals().reshape((3,ao_nshells)) # xyz coordinates of the center
    contents['ao_priexp']  = f.read_reals()    # exponential parameters
    contents['ao_priccf']  = f.read_reals()    # contraction coefficients
    f.close()
    return contents
