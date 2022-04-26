import numpy as np

Hbar_binary_relccsd = np.fromfile('M_small__.bin', dtype=np.complex, count=-1)

# hardcoded dimensions, but can get those from the number of elements in the 1d array read from file
Hbar_numpy_2d = np.reshape(Hbar_binary_relccsd, (8,8), order='F')

(evalues,evectors) = np.linalg.eig(Hbar_numpy_2d)

# sort eigenvalues and eigenvectors in ascending order
idx = evalues.argsort()[::1] 
evalues = evalues[idx]
evectors = evectors[:,idx]

counter = range(0,8)

print("RELCCSD formatted Hbar")
print("Eigenvalues:")
for i in counter :
   print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(evalues[i]),end=" ")
print("")

print("Eigenvectors")
for i in counter :
   for j in counter :
      print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(evectors[i][j]),end=" ")
   print("")

# now, we want to expand the HBar from RELCCSD to Exacorr, with the appropriate dimensionalities 
nocc=2
nvirt=4
