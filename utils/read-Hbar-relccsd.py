import numpy as np

print("read in RELCCSD formatted Hbar")
Hbar_binary_relccsd = np.fromfile('M_small__.bin', dtype=np.complex, count=-1)

number_elements = Hbar_binary_relccsd.shape[0]
number_rows = int(np.sqrt(number_elements))
print("Number of rows:",number_rows)
# hardcoded dimensions, but can get those from the number of elements in the 1d array read from file
Hbar_numpy_2d = np.reshape(Hbar_binary_relccsd, (number_rows, number_rows), order='F')

(evalues,evectors) = np.linalg.eig(Hbar_numpy_2d)

# sort eigenvalues and eigenvectors in ascending order
idx = evalues.argsort()[::1] 
evalues = evalues[idx]
evectors = evectors[:,idx]

counter = range(0,number_rows)

print("read in RELCCSD formatted Hbar")

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
print("read in EXACC formatted Hbar")
Hbar_binary_exacc = np.fromfile('M_small_EXACC.bin', dtype=np.complex, count=-1)
number_elements_exacc = Hbar_binary_exacc.shape[0]
number_rows_exacc = int(np.sqrt(number_elements_exacc))
print("Number of rows:",number_rows_exacc)
# hardcoded dimensions, but can get those from the number of elements in the 1d array read from file
Hbar_exacc_numpy_2d = np.reshape(Hbar_binary_exacc, (number_rows_exacc, number_rows_exacc), order='F')


(evalues_exacc,evectors_exacc) = np.linalg.eig(Hbar_exacc_numpy_2d)

# sort eigenvalues and eigenvectors in ascending order
idx = evalues_exacc.argsort()[::1]
evalues_exacc = evalues_exacc[idx]
evectors_exacc = evectors_exacc[:,idx]

counter = range(0,number_rows_exacc)

print("read in EXACC formatted Hbar")

print("Eigenvalues:")
for i in counter :
   print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(evalues_exacc[i]),end=" ")
print("")

print("Eigenvectors")
for i in counter :
   for j in counter :
      print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(evectors_exacc[i][j]),end=" ")
   print("")


