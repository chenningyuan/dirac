#-----------------------------H5 File-------------------

def get_dset(path, dsetName = "dset", h5name = "amfPCEC.h5"):
    import h5py                       
    import numpy as np    
    direc = "/".join(str(x) for x in path)
    with h5py.File(h5name, "r") as hdf:
        group = hdf.get(direc)  
        dset = np.array(group[dsetName])
    return dset

def get_spinFree(path):
    import numpy as np
    z,x,y = get_dset(path).shape
    dset = np.zeros(z*x*y).reshape(z,x,y)
    dset[0] = get_dset(path)[0]
    return dset

def get_spinOrbit(path):
    import numpy as np
    z,x,y = get_dset(path).shape
    dset = np.zeros(z*x*y).reshape(z,x,y)
    dset[1:] = get_dset(path)[1:]
    return dset

def del_h5group(path, fname = "amfPCEC.h5"):     #takes path as list of folder names
    import h5py                              
    with h5py.File(fname, "a") as hdf:
        del hdf["/".join(path)]
    print(path[-1]+ " was deleted from " +"/".join(path[:-1]))

def modify_hdf5(path, dset, dsetName = "dset", fname = "amfPCEC.h5"):   #overwrites existing path and file
    import  h5py    
    direc = "/".join(path)
    with h5py.File(fname, "a") as hdf:    
        group = hdf.get(direc)
        del group[dsetName]
        dset = group.create_dataset(dsetName, data = dset)
        dset.attrs["ARRAY SIZE"] = matrix.shape
        dset.attrs["TYPE"] = "NUMPY ARRAY"
    print("/".join(path) + "/" + dsetName + " was modified")

def write_hdf5(path, dset, dsetName = "dset", h5name = "amfPCEC.h5"):     #writes hdf5 with default name "amfPCEC.h5"       
    import  h5py                                                          #atomic mean-field picture-change error corrections
    
    direc = "/".join(str(x) for x in path)    
    with h5py.File(h5name, "a") as hdf:            
        try:
            group = hdf.create_group(direc)     #creates directory and fill with data
            dataset = group.create_dataset(dsetName, data = dset)
            dataset.attrs["ARRAY SIZE"] = matrix.shape
            dataset.attrs["TYPE"] = "NUMPY ARRAY" 
            print("create " + direc + "/" + dsetName)
        except ValueError:                   #arises, if path and file exists 
            modify_hdf5(path, dset,dsetName)   

#visscher's code extended by z. 21 - 23: outputs attributes and empty file paths
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
            data_dict[path+key]['discription'] = list(item.attrs.values())
        elif len(list(item)) == 0:
            data_dict[path+key] = {"empty"}
        elif isinstance(item, h5py._hl.group.Group):
           recursively_load_dict_contents_from_group(h5file, data_dict, path + key + '/')
    return

#-----------------------------Fortranfile-------------------

def extract_path(record):     # returns list of entries     
    strings = [str(x) for x in record]
    splitted = str(strings).split("'")[1::2]
    ints = [x for x in "".join(splitted[4:]) if x.isdigit()]   
    
    o_number = list(splitted.pop(1))[1]
    method, bs, n_modell = splitted[:3]
    e_config = "C"+"".join(ints[:2])+"O"+"".join(ints[2:])    
    return [method.strip(),o_number, bs.strip(),n_modell.strip(),e_config]

def read_ff(fname):
    from scipy.io import FortranFile, FortranEOFError, FortranFormattingError
    
    recl = ["a20","i4", "a20", "a6", "a1", "i4", "i4","a1","i4"] #RECL = length in bytes of each record in the file. 
    i = 0                                                        # unknown length of the recl due to varying electron configuration
    while True:                                                  # it's added up until record can be read
        try:                                                      
            f = FortranFile(fname, "r")     # resets filepointer
            record = f.read_record(*recl)
            dim = f.read_ints("i4")     #nrow, ncol, nzc1  
            matrix = f.read_record("f8").reshape((dim[2],dim[0],dim[1]))     #reshapes matrix((nzc1,nrows,ncols))
            f.close()
            break
        except ValueError:
            recl += "i4","i4"
            f.close()                                            
        except FortranEOFError:
            break
        except FortranFormattingError:
            break
            
    return record, matrix
