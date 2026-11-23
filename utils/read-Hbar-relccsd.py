import numpy as np



def main ():
   print("read in RELCCSD formatted Hbar")
   Hbar_binary_relccsd = np.fromfile('M_small__.bin', dtype=np.complex, count=-1)

   number_elements = Hbar_binary_relccsd.shape[0]
   number_rows = int(np.sqrt(number_elements))
   print("Number of rows:",number_rows)
# hardcoded dimensions, but can get those from the number of elements in the 1d array read from file
   Hbar_numpy_2d = np.reshape(Hbar_binary_relccsd, (number_rows, number_rows), order='F')

   diagonalize_and_print_values_vectors(Hbar_numpy_2d,number_rows,"Hbar RELCCSD")

   Hbar_numpy_2d_nodd = Hbar_numpy_2d
   counter_remove_dd = range(2,number_rows)
   for i in counter_remove_dd :
      for j in counter_remove_dd :
#     if i == j :
         Hbar_numpy_2d_nodd[i][j] = 0
#     else : 
#        Hbar_numpy_2d_nodd[i][j] = 0.0001

   diagonalize_and_print_values_vectors(Hbar_numpy_2d_nodd,number_rows,"Hbar RELCCSD with zeroed out HDD")


# now, we want to expand the HBar from RELCCSD to Exacorr, with the appropriate dimensionalities 
   nocc=2
   nvirt=6
   print("\nRead in EXACC formatted Hbar")
   Hbar_binary_exacc = np.fromfile('M_small_EXACC.bin', dtype=np.complex, count=-1)
   number_elements_exacc = Hbar_binary_exacc.shape[0]
   number_rows_exacc = int(np.sqrt(number_elements_exacc))
   print("Number of rows:",number_rows_exacc)
# hardcoded dimensions, but can get those from the number of elements in the 1d array read from file
   Hbar_exacc_numpy_2d = np.reshape(Hbar_binary_exacc, (number_rows_exacc, number_rows_exacc), order='F')

   diagonalize_and_print_values_vectors(Hbar_exacc_numpy_2d,number_rows_exacc,"Hbar EXACC")

   Hbar_exacc_numpy_2d_nodd = Hbar_exacc_numpy_2d
   counter_remove_dd = range(nocc,number_rows_exacc)
   for i in counter_remove_dd :
      for j in counter_remove_dd :
         Hbar_exacc_numpy_2d_nodd[i][j] = 0

   diagonalize_and_print_values_vectors(Hbar_exacc_numpy_2d,number_rows_exacc,"Hbar EXACC with zeroed out HDD")

def print_matrix(m,nrows,name="matrix") :
   counter = range(0,nrows)
   print(name)
   for i in counter :
      for j in counter :
        print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(m[i][j]),end=" ")
      print("")

def diagonalize_and_print_values_vectors(m,nrows,name="matrix"):
   print_matrix(m,nrows,name)
   (evalues,evectors) = np.linalg.eig(m)
# sort eigenvalues and eigenvectors in ascending order
   idx = evalues.argsort()[::1]
   evalues = evalues[idx]
   evectors= evectors[:,idx]
   print_vector(evalues,nrows,name="\n Eigenvalues:")
   print_matrix(evectors,nrows,name=" Eigenvectors")

def print_vector(v,nelements,name="vector"):
   counter = range(0,nelements)
   print(name)
   for i in counter :
      print('({0.real:8.4f} + {0.imag:8.4f}i)'.format(v[i]),end=" ")
   print("")

main()
