from scm.plams import *
import numpy as np
import sys, os
import copy
import shutil, os
import re

##### Define functions ######

def read_input_options(input_file):
    """ Reads input from a file and returns this as a dictionary """
    with open(input_file,"r") as f:
       token = f.read().splitlines()
       data_directory = token[0]
       if token[1] == "True":
          restricted = True
       else:
          restricted = False
       if token[2] == "True":
          openshell = True
       else:
          openshell = False
       Basis1      = token[3]
       Basis2      = token[4]
       charge      = int(token[5])
       multiplicity= int(token[6])
       nfragments  = int(token[7])
       frags_list  = []
       nmo_frag    = []
       frag_charge = []
       frag_multiplicity = []
       for i in range(nfragments):
          frags_list.append(token[8+i].split()[0])
          nmo_frag.append(int(token[8+i].split()[1]))
          frag_charge.append(int(token[8+i].split()[2]))
          frag_multiplicity.append(int(token[8+i].split()[3]))
       if token[8+nfragments] == "False":
          restart = False
       else:
          restart = True
       input_options = { "data_directory":data_directory,"Basis1":Basis1, "Basis2":Basis2, "charge":charge, "multiplicity":multiplicity, 
                         "restricted":restricted, "openshell":openshell, "restart":restart,
                         "nfragments":nfragments,"frags_list":frags_list, "nmo_frag":nmo_frag, "frag_charge":frag_charge,"frag_multiplicity":frag_multiplicity}
       return input_options

def add_fragments(frags_result):
    """ Takes a list of PLAMS result objects and combines them into a supermolecule """
    m_tot = Molecule()
    for i, frags in enumerate(frags_result):
        mol = frags.get_input_molecule().copy()
        for a in mol:
            a.properties.suffix = 'adf.f=frag'+str(i)

        m_tot += mol

    return m_tot

def extract_MOs(results,spin):
    """ Extracts MO coefficients from a PLAMS result object """
    NAO      = int(results.readrkf('Basis','naos',file='adf'))
    NMO      = int(results.readrkf('A','nmo_'+spin,file='adf'))
    MO       = results.readrkf('A','Eigen-Bas_'+spin,file='adf')
    e_orb    = np.array(results.readrkf('A','eps_'+spin,file='adf')).tolist()
    if isinstance(e_orb,float): e_orb = [e_orb]

    ##### Reorder the MO coefficient matrices ########
    order         = results.readrkf('A','npart',file='adf')
    if isinstance(order,int): order = [order]
    mocoef      = np.array(MO, dtype="float64").reshape(NMO, NAO)
    mocoef_sort = np.zeros((NMO, NAO), dtype="float64")
    for i in range(NMO):
        for j in range(NAO):
            mocoef_sort[i][order[j]-1] = mocoef[i][j]

    return NAO, NMO, e_orb, order, mocoef_sort

def extract_other_data(results):
    """ Extracts essential other data from a PLAMS result object """
    nel         = int(results.readrkf('General','electrons',file='adf'))
    natoms      = int(results.readrkf('Geometry','nr of atoms',file='adf'))
    nspin       = int(results.readrkf('General','nspin',file='adf'))
    geom        = np.array(results.readrkf('Geometry','xyz InputOrder',file='adf')).reshape(natoms,3)
    tot_enrgy   = float(results.readrkf('Energy','Bond Energy',file='adf'))

    return nel, natoms, nspin, geom, tot_enrgy

##### Execute the desired type of run ####

# First define some general AMS/ADF settings
settings = Settings()
settings.input.ams.task = "SinglePoint"
settings.input.adf.symmetry = "NOSYM"
settings.input.adf.basis.core = "None"
settings.input.adf.xc.gga = "PBE"
settings.input.adf.noprint = "Logfile"
settings.input.adf.save = "TAPE15"
#settings.input.adf.Relativity.Level = "None"    # CAREFUL : IN ADF2020-21 ZORA is default 
# Then get the user-defined options into a dictionary
input_options = read_input_options("INPUT_ADF")

# Put some generic options also in variables for convenience
restart    = input_options["restart"]
restricted = input_options["restricted"]
openshell  = input_options["openshell"]
nfragments = input_options["nfragments"]
nmo_frag   = input_options["nmo_frag"]
newInterface = True

if not restart:
###### Molecule in Basis B1###############
########## Running ADF for Molecule in Basis B1########
    B1_settings = settings.copy()
    B1_settings.input.adf.basis.type = input_options["Basis1"]
    B1_settings.input.ams.system.charge = str(input_options["charge"])
    B1_settings.input.adf.spinpolarization = str(input_options["multiplicity"] - 1)
    if ((not restricted) or openshell):
        B1_settings.input.adf.Unrestricted = " "
    mol = Molecule(input_options["data_directory"] + "/MOLECULE.xyz")
    print("Starting ADF calculation for Molecule in Basis B1...")
    results_B1 = AMSJob(molecule=mol, settings=B1_settings, name='Molecule_B1_MO1').run()

    adf_file   = results_B1["adf.rkf"]
    tape15file = results_B1["TAPE15"]

    # copy the required rkf files (supermolecule coefficients and sm/fragments overlap)
    shutil.copy(adf_file,  os.getcwd()+"/CHECKPOINT.rkf")
    shutil.copy(adf_file,  os.getcwd()+"/Molecule_B1_MO1.rkf")
    shutil.copy(tape15file,os.getcwd()+"/Molecule_B1_MO1.t15")
    # Write path to SCF calculation B1 in RESTART_AMS
    with open("RESTART_AMS","w") as f:
         f.write(os.getcwd()+"/Molecule_B1_MO1_ribo.rkf")

    print("Molecule Basis B1 Successful")


########## fragment calculation ##################
########## Running ADF for fragment calculation ###################

    mols_list = [Molecule(m) for m in input_options["frags_list"]]

    frag_settings = settings.copy()
    frag_settings.input.adf.Basis.Type = input_options["Basis2"]

    if not restricted or openshell:
        frag_settings.input.adf.Unrestricted = " "

    frag_TAPE21_original = []
    frag_TAPE21_NMO2 = []
    for index, mol in enumerate(mols_list):
        frag_settings.input.ams.system.Charge = str(input_options["frag_charge"][index])
        frag_settings.input.adf.SpinPolarization = str(input_options["frag_multiplicity"][index] - 1)
        name = 'frag' + str(index+1)
        print("Starting ADF calculation for frag"+str(index+1) + "...")
        frag_B2 = AMSJob(molecule=mol, settings=frag_settings, name=name).run()
        source = frag_B2["adf.rkf"]
        destination = os.getcwd()+"/"+ name + ".rkf"
        shutil.copy(source, destination)   

        frag_TAPE21_original.append(destination)
        frag_TAPE21_NMO2.append(frag_B2)

######## Running ADF calculation for overlap matrix ############

    tot_list = [results_B1["adf.rkf"]] + frag_TAPE21_original
    m_tot = add_fragments(([results_B1] + frag_TAPE21_NMO2 ))
    ovrlp_settings = settings.copy()
    s_frag = ovrlp_settings.input.adf.fragments

    for i, frags in enumerate(tot_list):
        s_frag['frag'+str(i)] = frags

    ovrlp_settings.input.adf.CALCOVERLAPONLY = ""
    ovrlp_settings.input.ams.system.CHARGE = str(input_options["charge"])
    ovrlp_settings.input.ams.system.AllowCloseAtoms = ""
    ovrlp_settings.input.adf.SpinPolarization = str(input_options["multiplicity"]-1)

    if not restricted or openshell:
        ovrlp_settings.input.adf.UNRESTRICTEDFRAGMENTS = " "
        ovrlp_settings.input.adf.UNRESTRICTED = " "

    print("Starting Calculation of overlap matrix")

    ovrlp_job = AMSJob(molecule=m_tot, settings=ovrlp_settings, name="Molecule_B2_S12").run()

    source = ovrlp_job["TAPE15"]
    destination = os.getcwd()+"/Molecule_B2_S12.t15"
    shutil.copy(source, destination)

    print("ADF Calculation completed successfully!!!")

##### ADF calculations completed ################

else: # RESTART ADF with new MO coefficients

    with open("RESTART_AMS","r") as f:
         destination_ibo = f.read()

    os.remove(os.getcwd()+"/Molecule_B1_MO1.t15")
    os.remove(os.getcwd()+"/Molecule_B2_S12.t15")
    for index in range(nfragments):
        name = 'frag' + str(index+1)
        os.remove(os.getcwd()+"/"+ name + ".rkf")

    # Restart calculation
    restart_settings = settings.copy()
    restart_settings.input.ams.EngineRestart = destination_ibo
    restart_settings.input.adf.Basis.Type = input_options["Basis1"]
    restart_settings.input.ams.system.CHARGE = str(input_options["charge"])
    restart_settings.input.adf.SpinPolarization = str(input_options["multiplicity"]-1)
    if openshell or not restricted:
        restart_settings.input.adf.UNRESTRICTED = " "
    mol = Molecule(input_options["data_directory"] + "/MOLECULE.xyz")
    print("Starting ADF calculation for Molecule in the new IBO basis...")
    results_new = AMSJob(molecule=mol, settings=restart_settings, name='Molecule_B1_new').run()
    print("ADF calculation in the new basis ended successfully")
