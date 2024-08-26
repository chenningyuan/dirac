#!/bin/env python
import sys
import h5py
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-fcsv', help='CSV file on input')
parser.add_argument('-csv_sep', default=',', help='separator on CSV file')
parser.add_argument('-fh5', help='HDF5 file on output')
parser.add_argument('-selected_columns', nargs='+', help='selected labels to be included in h5 file')
args = parser.parse_args()

data=pd.read_csv(args.fcsv, header=0, sep=args.csv_sep)
columns=data.columns
if args.selected_columns is not None:
    columns=args.selected_columns
    
selected_data=data[columns]
print(selected_data)
if args.fh5 is not None:
    selected_data.to_hdf(args.fh5, 'data', mode='w', format='table')
else:
    print('missing file name for HDF5 file')

# grid_function
#if args.data_labels is not None:
#    for k in data_labels:
#        for x in data.keys():
#            if k in x:
#                data_values[k] = data[x]["value"]
#    df2=pd.DataFrame.from_dict(data_values)
#    if args.fh5 is not None:
#        df2.to_hdf(args.fh5, 'data', append=True)
#    else:
#        print('missing file name for HDF5 file')
#
#else:
#    print('specify grid function labels')
#    sys.exit()
#
#
#
#
