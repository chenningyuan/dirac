from __future__ import absolute_import

import sys
import numpy as np

try:
    from ._run_pyscf_mocoeff import *
except ImportError:
    pass

try:
    from ._run_psi4_mocoeff import *
except ImportError:
    pass


def construct_B2_from_IFOs(molecule,
                           basis,
                           restricted,
                           nao,
                           wfn_object,
                           avas_frag,
                           mo_avas_frag,
                           data_directory,
                           extension):
    """
    This function combine orbitals (expressed in the same basis of nao basis functions)
    from different fragments to generate a B2 basis, used in AVAS for instance
    """

    B2_dict = {}
    nmo_avas = sum([len(mos) for mos in mo_avas_frag])
    alpha_energy = np.zeros(nmo_avas)
    alpha_MO = np.zeros((nao,nmo_avas))
    if not restricted: 
      beta_energy = np.zeros(nmo_avas)
      beta_MO = np.zeros((nao,nmo_avas))

    offset = 0
    # We'll take the orbitals from mo_avas_frag for each considered fragment.
    for frag in avas_frag:
      with open(data_directory+"/frag{}_iao.".format(frag)+extension,"r") as f:
         alpha_energies = read_real_list("Alpha Orbital Energies", f)
         alpha_coeff    = read_real_list("Alpha MO coefficients", f)
         beta_energies  = read_real_list("Beta Orbital Energies", f)
         beta_coeff     = read_real_list("Beta MO coefficients", f)
      f.close()

      nmo = len(alpha_energies)
      alpha_frag_coeff = np.zeros((nao,nmo),dtype=float)
      beta_frag_coeff = np.zeros((nao,nmo),dtype=float)

      ij = 0
      for i in range(nmo):
         for j in range(nao):
              alpha_frag_coeff[j,i] = alpha_coeff[ij]
              if not restricted:
                 beta_frag_coeff[j,i] = beta_coeff[ij]
              else:
                 beta_frag_coeff[j,i] = alpha_coeff[ij]
              ij += 1

      i = 0
      for mo in mo_avas_frag[frag]:
         if restricted:
            alpha_MO[:,offset+i] = alpha_frag_coeff[:,mo-1]
            alpha_energy[offset+i] = alpha_energies[mo-1]
         else:
            alpha_MO[:,offset+i] = alpha_frag_coeff[:,mo-1]
            alpha_energy[offset+i] = alpha_energies[mo-1]
            beta_MO[:,offset+i] = beta_frag_coeff[:,mo-1]
            beta_energies[offset+i] = beta_energies[mo-1]
         i += 1
      offset += len(mo_avas_frag[frag])

      alpha_MO_list = []
      for i in range(alpha_MO.shape[1]):
        for j in range(alpha_MO.shape[0]):
          alpha_MO_list.append(alpha_MO[j][i])
     
      if not restricted:
        beta_MO_list = [] 
        for i in range(beta_MO.shape[1]):
          for j in range(beta_MO.shape[0]):
            beta_MO_list.append(beta_MO[j][i])

    if restricted:
      B2_dict = {'wfn_object':wfn_object,'alpha_MO':alpha_MO_list,'alpha_energies':alpha_energy}
    else:
      B2_dict = {'wfn_object':wfn_object,'alpha_MO':alpha_MO_list,'alpha_energies':alpha_energy,
                 'beta_MO':beta_MO_list,'beta_energies':beta_energy}

    if extension == "psi4":
       run_psi4_mocoeff(molecule,
                        basis=basis,
                        restricted=restricted,
                        just_write=B2_dict)
    elif extension == "pyscf":
       run_pyscf_mocoeff(molecule,
                        basis=basis,
                        restricted=restricted,
                        just_write=B2_dict)
    else:
       sys.exit("extension should be psi4 or pyscf for now")
