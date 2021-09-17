#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 09:03:56 2021

@author: Lucas Visscher

Read/write data from DIRAC labeled files
"""

import numpy as np
from scipy.io import FortranFile

###### Functions to read labels and data

def read_label(f):
    """
    Read a DIRAC style label from an open file
    """
    record = f.read_record(np.dtype('4a8'))
    stars=record[0].decode("utf-8")
    if stars != '********':
        raise IOError('Error reading label, file may be corrupted')
    time_run=record[1].decode("utf-8")
    date_run=record[2].decode("utf-8")
    label=record[3].decode("utf-8")
    
    return label,date_run,time_run

def read_AOPROPER(file_name):
    """
    Open AOPROPER type file and return dictionary of its contents
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


def read_DFCOEF(file_name):
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
    record = f.read_record('a74', '<i4','<i4','3<i4','f8')
    contents['title']  = record[0]
    contents['nfsym']  = record[1]
    contents['nz']     = record[2]
    contents['npo']    = record[3][0]
    contents['nmo']    = record[3][1]
    contents['nto']    = contents['nmo'] + contents['npo'] 
    contents['nao']    = record[3][2]
    contents['energy'] = record[4][0]
    label, date_run, time_run = read_label(f)
    contents['mo_coeff'] = f.read_reals()
    label, date_run, time_run = read_label(f)
    contents['eigenval'] = f.read_reals()
    label, date_run, time_run = read_label(f)
    contents['ibeig'] = f.read_ints()
    label, date_run, time_run = read_label(f)
    # 6 records with basis set information
    ao_nshells = f.read_ints()[1]
    contents['ao_nshells'] = ao_nshells    # number of distint ao shells, then for each shell:
    ints = f.read_ints().reshape((3,ao_nshells))
    contents['ao_orbmom']  = ints[0,:] # orbital momentum (s=1,p=2, d=3, etc.)
    contents['ao_idcent']  = ints[1,:] # id of its center
    contents['ao_nprim']   = ints[2,:] # number of primitives
    ints = f.read_ints().reshape((2,ao_nshells))
    contents['ao_jstrt']   = ints[0,:] # start address of contraction coefficients
    contents['ao_numcf']   = ints[1,:] # column number of contraction coefficients
    contents['ao_cent']    = f.read_reals().reshape((3,ao_nshells)) # xyz coordinates of the center
    contents['ao_priexp']  = f.read_reals()    # exponential parameters
    contents['ao_priccf']  = f.read_reals()    # contraction coefficients
    f.close()
    return contents

# Examples of usage. 
file_name = 'DFCOEF'
contents = read_DFCOEF(file_name)
print('Energy = ',contents['energy'])
print('NAO    = ',contents['nao'])
    
file_name = 'AOPROPER'
contents = read_AOPROPER(file_name) 
print("Overlap Matrix")
np.set_printoptions(precision=3,linewidth=80)
print(contents['OVERLAP '])


    
    

