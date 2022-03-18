#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 2022

@maintaining author: Stefan Knecht
# copied from
# Ref: https://stackoverflow.com/a/34401029/10462884
# with slight modifications

Combine different amfPCEC.h5 atom files into one for simple use in molecular run

Example of usage:
python merge_amf.py atom1/amfPCEC.h5 atom2/amfPCEC.h5 atom3/amfPCEC.h5
"""

import h5py
from sys import argv, exit

def h5py_dataset_iterator(g, prefix=''):
    for name, h5obj in g.items():
        path = '{}/{}'.format(prefix, name)
        if isinstance(h5obj, h5py.Dataset):  # test for dataset
            yield (h5obj, path)
        elif isinstance(h5obj, h5py.Group):  # test for group (go down)
            yield from h5py_dataset_iterator(h5obj, prefix=path)


with h5py.File('amfPCEC.h5', mode='w') as h5w:

    print("Total number of input files passed :", len(argv)-1)
    if(len(argv)-1 == 0):
        exit('no input .h5 files are given! I need at least two files.')
    print("Input files are:")
    for value in argv[1:]:
        print(value)

    for h5source in argv[1:]:
        print(f'\nWorking on file: {h5source}')
        with h5py.File(h5source, mode='r') as h5r:
            for (dset, path) in h5py_dataset_iterator(h5r):
                print(f'Copying dataset from: {path}')
                ds_obj = h5r[path]
                arr_dtype = ds_obj.dtype
                arr_shape = ds_obj.shape

                # If dataset doesn't exist, create new datset and copy data
                # Note: we can't use .copy() method b/c we are changing shape and maxshape
                if path not in h5w:
                    h5w.create_dataset(path, data=ds_obj,
                                       shape=arr_shape+(1,), maxshape=arr_shape+(None,))
                else:
                    # Check for compatiable dtype and shape
                    ds_dtype = h5w[path].dtype
                    ds_shape = h5w[path].shape
                    # If dataset exists and is compatibale, resize datset and copy data
                    if ds_dtype == arr_dtype and ds_shape[:-1] == arr_shape:
                        new_shape = ds_shape[0:-1] + (ds_shape[-1]+1,)
                        h5w[path].resize(new_shape)
                        h5w[path][...,-1] = ds_obj
                # Add attribute to dataset with source file name:
                h5w[path].attrs[f'File for index {h5w[path].shape[-1]-1} '] = h5source
