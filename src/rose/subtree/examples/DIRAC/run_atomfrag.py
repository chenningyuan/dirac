from genibo_MOcoeff import periodic_table
import subprocess

# This script is not working, nice idea but needs definition of xyz files.
for atom in periodic_table:
  try:
    subprocess.check_call("python $ROSEDIR/python_scripts/atomfrag_dirac.py " + atom + "_atom 0 1> " + atom + "_test 2> " + atom + "_error", shell = True)
    print("Success for atom ", atom)
  except:
    print("Failure for atom ", atom)
