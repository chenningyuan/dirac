      parameter (Lmax=4)   ! max. angular momentum of basis functions  DO NOT INCREASE
cbs                          TO MORE THAN four  !!!!!!!!!!!!!!!!!!!
cbs                          if you do, you will have to edit the ixyzpow array by hand...............
cbs                          (in datapow.h)
cbs  you will also have to add the correct normalizations when using 
CBS  the old DALTON Normalisation..
      parameter (Lmax_occ=3)   ! highest L-value for occupied orbitals                     
cmi   parameter (MxprimL=30) ! max. of primitives per angular momentum
      parameter (MxprimL=35) ! max. of primitives per angular momentum
cmi   parameter (MxcontL=30) ! max. of contracted functions per angular momentum
      parameter (MxcontL=35) ! max. of contracted functions per angular momentum
cmi   parameter (MxCart=360) ! max. number of contracted functions in the atom
      parameter (MxCart=320) ! max. number of contracted functions in the atom
cmi   parameter (MxCart=300) ! max. number of contracted functions in the atom

      parameter (ndfmx=4*Lmax+4) ! dimension of precomputed double factorials     
