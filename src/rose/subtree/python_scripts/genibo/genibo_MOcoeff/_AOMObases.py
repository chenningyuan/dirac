# Classes to define and read / write atomic and molecular orbital bases
# Used in ROSE and DIRAC, be sure to keep this consistent !
# Requires numpy and h5py
# Lucas Visscher - 2024
import numpy as np
import h5py

# The only dependence on h5py is contained in the following two functions

def recursively_load_dict_contents_from_group(h5file, data_dict, path):
    """
    Get a flattened dictionary with all data contained in a hdf5 file
    """
    import h5py
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            data_dict[path+key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            recursively_load_dict_contents_from_group(h5file, data_dict, path + key + '/')
    return

def read_hdf5(file_name):
    """
    Open hdf5-type file and return dictionary of its contents
    """
    import h5py
    data_dict = {}
    with h5py.File(file_name, 'r') as h5file:
        recursively_load_dict_contents_from_group(h5file, data_dict,'/')
    return data_dict

# the AO class

class BasisFunc():
    """
    Defines an AO basis function 
    """
    def __init__(self,coord=[0.,0.,0.],orb_momentum=1,\
                 exponent=[1.0],coefficient=[1.0]):
        self.coord    = coord
        self.orb_momentum = orb_momentum
        self.exponent = exponent
        self.coefficient = coefficient
        if len(exponent) != len(coefficient): raise ValueError
        self.n_primitives = len(exponent)
 
    def __repr__(self):
        return f"""<BasisFunc 
        coord:{self.coord} 
        orb_momentum:{self.orb_momentum}
        exponent:{self.exponent}
        coefficient:{self.coefficient}
        >"""

    def __str__(self):
        return f"""Basis function defined as:
                coordinates: {self.coord}
                orbital momentum: {self.orb_momentum}
                exponent(s): {self.exponent}
                coefficient(s): {self.coefficient}
                """
    def nfunctions(self,angular):
        if angular == 1: 
            return int(self.orb_momentum * (self.orb_momentum + 1) / 2) # Cartesian
        else:
            return 2*self.orb_momentum - 1                              # Spherical
        
class AOBasis():
    """
    Defines AO basis 
    """
    def __init__(self,angular=1,basis_funcs=[]):
        self.angular = angular # 1=cartesian, 2=spherical
        self.nshells = len(basis_funcs)
        self.basis_funcs = basis_funcs
        # finish by counting the number of individual functions and create a look-up table for shells
        self.nao = 0
        self.shell_indices = []
        for shell, basis_func in enumerate(self.basis_funcs):
            n_in_shell = basis_func.nfunctions(angular)
            self.shell_indices += n_in_shell * [shell+1]
            self.nao += n_in_shell
            
    def __str__(self):
        return f"""AO basis set characteristics:
                angular: {self.angular}
                number of shells: {self.nshells}
                number of basis functions: {self.nao}
                """
    
    def add_AOs(self,basis_funcs=[]):
        self.basis_funcs += basis_funcs
        # update the number of functions and the look-up table
        for shell, basis_func in enumerate(basis_funcs):
            n_in_shell = basis_func.nfunctions(angular)
            self.shell_indices += n_in_shell * [self.nshells+shell+1]
            self.nao += n_in_shell
        self.nshells = len(self.basis_funcs)
            
    def write_hdf5(self,file_name='AOMO.h5',section='/aobasis/'):
        """
        Write AO basis to hdf5-type file
        """
        h5file = h5py.File(file_name, 'a')
        # Arrays are written as extendable datasets, single strings or scalars are fixed size
        h5file.create_dataset(section+'angular',data=[self.angular])
        h5file.create_dataset(section+'n_shells',data=[self.nshells])
        h5file.create_dataset(section+'n_ao',data=[self.nao])
        orbmom = []
        n_prim = []
        n_cont = []
        center = []
        exponents = []
        contractions = []
        for basis_func in self.basis_funcs:
            orbmom.append(basis_func.orb_momentum)
            n_prim.append(basis_func.n_primitives)
            n_cont.append(1) # we don't allow general contraction to keep things simple
            center.append(basis_func.coord)
            exponents.append(basis_func.exponent)
            contractions.append(basis_func.coefficient)
        # we need to concatenate lists of numpy arrays to be able to write it as 1D array
        # lists of single numbers can be written directly
        center = np.hstack(center)
        exponents = np.hstack(exponents)
        contractions = np.hstack(contractions)
        h5file.create_dataset(section+'orbmom',data=orbmom, maxshape=(None))
        h5file.create_dataset(section+'n_prim',data=n_prim, maxshape=(None))
        h5file.create_dataset(section+'n_cont',data=n_cont, maxshape=(None))
        h5file.create_dataset(section+'center',data=center, maxshape=(None))
        h5file.create_dataset(section+'exponents',data=exponents, maxshape=(None))
        h5file.create_dataset(section+'contractions',data=contractions, maxshape=(None))
        h5file.close()

class AOBasisFromFile(AOBasis):
    """
    Defines AO basis by reading it from a hdf5-type file
    """
    def __init__(self,basis_file="CHECKPOINT.h5",section='/input/aobasis/1/'):
        h5file       = h5py.File(basis_file, 'r')
        angular      = h5file[section+'angular'][()][0]
        n_shells     = h5file[section+'n_shells'][()][0]
        orbmom       = h5file[section+'orbmom'][()]
        n_prim       = h5file[section+'n_prim'][()]
        n_cont       = h5file[section+'n_cont'][()]
        center       = h5file[section+'center'][()]
        exponents    = h5file[section+'exponents'][()]
        contractions = h5file[section+'contractions'][()]
        h5file.close()
        
        basis_funcs = []
        iPrimStart = 0
        for shell in range(n_shells):
            iPrimEnd = iPrimStart + n_prim[shell]
            basis_funcs.append( \
                BasisFunc( \
                coord        = center[shell*3:(shell+1)*3], \
                orb_momentum = orbmom[shell], \
                exponent     = exponents[iPrimStart:iPrimEnd],
                coefficient  = contractions[iPrimStart:iPrimEnd] ))
            iPrimStart += n_prim[shell]
        
        super().__init__(angular,basis_funcs) 

class MO():
    """
    Defines a molecular orbital 
    - nspin: 1 for spin-orbital, 2 for spinors that have both alpha and beta parts
    - algebra: real(r), complex(c) or quaternion(q)
    - eigenvalue: orbital energy
    - coeff: MO coefficients, 3-dimensional array [nspin,nmo,n_basis]
    - aobasis: pointer to the definition of the AO basis
    """
    def __init__(self,aobasis,coeff, nspin=1, algebra='r',eigenvalue=0.0):
        self.aobasis    = aobasis
        self.nspin      = nspin
        self.algebra    = algebra
        self.coeff      = coeff
        self.eigenvalue = eigenvalue
  
    def __repr__(self):
        return f"""<MO 
        algebra:{self.algebra} 
        eigenvalue:{self.eigenvalue}
        coeff:{self.coeff}
         >"""

    def __str__(self):
        return f"""Molecular Orbital defined as:
                algebra: {self.algebra}
                orbital energy: {self.eigenvalue}
                MO coeffcients: {self.coeff}
                """
    def remove_phase(self):
        """
        Remove the phase: convert quaternion MO into a real one.
        NB: will only work correctly for spinfree orbitals !
        """
        nzsp = np.shape(self.coeff)[0]
        nao  = self.aobasis.nao
        rijk_norm = np.zeros(nzsp)
        rcoeff = np.zeros((1,nao))
        for rijk in range(nzsp):
            rijk_norm[rijk] = np.linalg.norm(self.coeff[rijk,:])
        rcoeff[0,:]  = self.coeff[np.argmax(rijk_norm),:]
        self.coeff   = rcoeff
        self.nspin   = 1
        self.algebra = 'r'
        
class MOBasis():
    """
    MO basis (a collection of MOs)
    MOs can have real(r), complex(c) or quaternion(q) coefficient arrays
    MOs can be organized in one (restricted, nbar=1) or two sets (unrestricted, nbar=2) of orbitals
    """
    def __init__(self,mos=None,aobasis=None,eigenvalues=None,coeffs=None,algebra='r',restricted=True):
        self.aobasis = aobasis
        if restricted:
            self.nbar = 1
        else:
            self.nbar = 2
        if mos != None:
            # initialization via a list of predefined MO objects
            self.nmo = len(mos)
            self.MOs  = mos
        else:  
            # initialization using separate eigenvalue and coefficient arrays that refer to the same aobasis
            self.nmo   = int(len(eigenvalues)/self.nbar)
            if algebra == 'q': nzsp, nspin = 4, 2 # quaternion packed spinor pairs, factor 4 larger than real coefficients
            if algebra == 'c': nzsp, nspin = 4, 2 # complex spinors, also factor 4 larger (alpha/beta and real/imag)
            if algebra == 'r': nzsp, nspin = 1, 1 # real spin-orbitals  
            n_basis  = int(coeffs.size / (self.nmo * self.nbar * nzsp))
            shape = (self.nbar,nzsp,self.nmo,n_basis) # NB: reverse order compared to the Fortran definition
            coeffs.resize(shape)
            self.MOs = []
            for bar in range(self.nbar):
                for orbital in range(self.nmo):
                    self.MOs.append(MO(aobasis=aobasis,algebra=algebra,nspin=nspin,
                              eigenvalue=eigenvalues[orbital],
                              coeff=coeffs[bar,:,orbital,:]))
            
    def __str__(self):
        bartype = ['restricted','unrestricted']
        return f"""MO basis characteristics:
                type: {bartype[self.nbar-1]}
                algebra: {self.MOs[0].algebra}
                number of AOs: {self.MOs[0].aobasis.nao}
                number of MOs: {self.nmo}
                """
    
    def add_MOs(self,mos):
        self.nmo  += len(mos)
        self.MOs  += mos
        
    def remove_phases(self):
        for mo in self.MOs:
            mo.remove_phase()
   
    def write_hdf5(self,file_name='AOMO.h5',
                   aosection='/aobasis/',
                   aoset='',
                   mosection='/mobasis/',
                   mosym=''):
        """"
        Write MO basis to hdf5-type file
        aosection and mosection are used to specify the respective sections
        the aoset modifier is used for DIRAC that can have multiple ao basis sets
        the mosym modifier is used for DIRAC that stores mo coefficients with and without symmetry adaptation
        """
        # All MOs should refer to the same AO basis, we retrieve this from the first MO
        n_basis = self.MOs[0].aobasis.nao
        algebra = self.MOs[0].algebra
        if algebra == 'q': 
            nz = 4
        elif algebra == 'c': 
            nz = 2
        else:
            nz = 1

        # Write the aobasis as well to make sure this is consistent with the MOs
        self.MOs[0].aobasis.write_hdf5(file_name=file_name,section=aosection+aoset)
             
        # Write the MO section
        h5file = h5py.File(file_name, 'a')
        # Arrays are written as extendable datasets, single strings or scalars are fixed size
        h5file.create_dataset(mosection+'n_mo',data=[self.nmo])
        h5file.create_dataset(mosection+'nz'+mosym,data=[nz])
        h5file.create_dataset(mosection+'n_basis',data=[n_basis])
        
        eigenvalues = np.zeros((self.nbar*self.nmo))
        shape = (nz,self.nbar*self.nmo,n_basis) # NB: reverse order compared to the Fortran definition
        coeffs = np.zeros(shape)
        
        for orbital in range(self.nbar*self.nmo):
            eigenvalues[orbital]= self.MOs[orbital].eigenvalue
            coeffs[:,orbital,:] = self.MOs[orbital].coeff
   
        h5file.create_dataset(mosection+'eigenvalues'+mosym,data=eigenvalues.ravel(),chunks=True,maxshape=(None))
        h5file.create_dataset(mosection+'orbitals'+mosym,data=coeffs.ravel(),chunks=True,maxshape=(None))
        h5file.close()
        
    def make_mostring(self,mo_range=None,energy_tresholds=None):
        if mo_range != None:
            # make text string based on this list
            return '{:1}..{:1}\n'.format(*mo_range)
        else:
            return 'energy {:10.5f} {:10.5f} {:10.5f}\n'.format(*energy_tresholds)
        
    def get_eigenvalues(self,mostring):
        eigenvalues = []
        if 'energy' in mostring:
            lower_bound = float(mostring.split()[1])
            upper_bound = float(mostring.split()[2])
            for MO in self.MOs:
                if MO.eigenvalue >= lower_bound and MO.eigenvalue <= upper_bound:
                    eigenvalues.append(MO.eigenvalue)
                #todo: implement addition of eigenvalues just below or above this range 
            return np.array(eigenvalues)
        else:
            desired_MOs = mostring.split(',')
            for part in desired_MOs:
                if '..' in part: # range of MOs
                    start, end = mostring.split('..')
                    for imo in range(int(start)-1,int(end)):
                        eigenvalues.append(self.MOs[imo].eigenvalue)
                else: # single MO
                    imo = int(part)
                    eigenvalues.append(self.MOs[imo].eigenvalue)
   
        return np.array(eigenvalues)

    def get_mocoeffs(self,mostring):
        eigenvectors = []
        if 'energy' in mostring:
            lower_bound = float(mostring.split()[1])
            upper_bound = float(mostring.split()[2])
            for MO in self.MOs:
                if MO.eigenvalue >= lower_bound and MO.eigenvalue <= upper_bound:
                    eigenvectors.append(MO.coeff)
                #todo: implement addition of eigenvalues just below or above this range
        else:
            desired_MOs = mostring.split(',')
            for part in desired_MOs:
                if '..' in part: # range of MOs
                    start, end = mostring.split('..')
                    for imo in range(int(start)-1,int(end)):
                        eigenvectors.append(self.MOs[imo].coeff)
                else: # single MO
                    imo = int(part)
                    eigenvectors.append(self.MOs[imo].coeff)
        
        # convert to numpy array, also change ordering to nz last
        nmo_selected = len(eigenvectors)
        nz = eigenvectors[0].shape[0]
        nao = eigenvectors[0].shape[1]
        mocoeffs = np.zeros((nmo_selected,nao,nz))
        for ix, eigenvector in enumerate(eigenvectors):
            for iz in range(nz):
                mocoeffs[ix,:,iz] = eigenvector[iz,:]
            
        return mocoeffs

        
class MOBasisFromFile(MOBasis):
    """
    Initializes MO basis by reading it from a hdf5 file
    aosection and mosection are used to specify the respective sections
    the aoset modifier is used with DIRAC that can have multiple ao basis sets
    the mosym modifier is used with DIRAC that can store mo coefficients with and without symmetry adaptation
    """
    def __init__(self,
                 basis_file='AOMObasis.h5',
                 aosection='/aobasis/',
                 aoset='',
                 mosection='/mobasis/',
                 mosym=''):
        aobasis      = AOBasisFromFile(basis_file,section=aosection+aoset)
        h5file       = h5py.File(basis_file, 'r')
        eigenvalues  = h5file[mosection+'eigenvalues'+mosym][()]
        coeffs       = h5file[mosection+'orbitals'+mosym][()]
        nz           = h5file[mosection+'nz'+mosym][()][0]
        if nz == 4: 
            algebra = 'q'
        elif nz == 2:
            algebra = 'c'
        else:
            algebra = 'r'
        super().__init__(aobasis=aobasis,eigenvalues=eigenvalues,coeffs=coeffs,algebra=algebra,restricted=True)
        
# This class has become obsolete, remove after synchronization and adaptation within DIRAC
class MOBasisFromDIRAC(MOBasis):
    """
    Initializes MO basis by reading it from a DIRAC file
    Note that .ACMOUT should be specified in input as we cannot handle symmetry
    """
    def __init__(self,basis_file="CHECKPOINT.h5", \
                 aosection='/input/aobasis/3/', \
                 mosection='/result/wavefunctions/scf/mobasis/'):
        aobasis      = AOBasisFromFile(basis_file,section=aosection)
        h5file       = h5py.File(basis_file, 'r')
        eigenvalues  = h5file[mosection+'eigenvalues_C1'][()]
        coeffs       = h5file[mosection+'orbitals_C1'][()]
        nz           = h5file[mosection+'nz_C1'][()][0]
        if nz == 4: 
            algebra = 'q'
        elif nz == 2:
            algebra = 'c'
        else:
            algebra = 'r'
        super().__init__(aobasis=aobasis,eigenvalues=eigenvalues,coeffs=coeffs,algebra=algebra,restricted=True)
        
