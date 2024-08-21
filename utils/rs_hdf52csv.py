#!/bin/env python
import sys
import h5py
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-fh5', help='HDF5 file on input')
parser.add_argument('-fcsv', help='CSV file on output')
parser.add_argument('-ftxt', help='TXT file on output (old DIRAC)')
parser.add_argument('-fcube', help='CUBE file on output')
parser.add_argument('-data', help='data to focus on: "grid" or "grid_function"')
parser.add_argument('-data_labels', help='select labels of grid functions')
args = parser.parse_args()

#
# functions copied from dirac_data.py
#
def read_hdf5(file_name):
    """
    Open hdf5-type file and return dictionary of its contents
    """
    import h5py
    data_dict = {}
    with h5py.File(file_name, 'r') as h5file:
        recursively_load_dict_contents_from_group(h5file, data_dict,'/')
    return data_dict

def recursively_load_dict_contents_from_group(h5file, data_dict, path):
    """
    Modified from code found at Stack Exchange to get flat dictionary
    """
    import h5py
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            data_dict[path+key] = {}
            data_dict[path+key]['value'] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
           recursively_load_dict_contents_from_group(h5file, data_dict, path + key + '/')


#
# other functions added for backward compatibility with older DIRAC
#

def write_cube():

data=read_hdf5(args.fh5)

if args.data == 'grid':
    labels=['coor_'+str(ind) for ind in range(4)]
    grid_xyz={}
    for k in labels:
        for x in data.keys():
            if k in x:
                grid_xyz[k] = data[x]["value"]
    
    
    df=pd.DataFrame.from_dict(grid_xyz)
    if args.fcsv is not None:
        df.to_csv(args.fcsv,index=False,float_format="%.8e")
    else:
        print('grid can only be converted to CSV file')

elif args.data == 'grid_function':
    if args.data_labels is not None:
        data_values={}
        for k in data_labels:
            for x in data.keys():
                if k in x:
                    data_values[k] = data[x]["value"]
    else:
        print('specify grid function labels')
        sys.exit()
else:
    print('unknown argument to "-data"')
    sys.exit()




