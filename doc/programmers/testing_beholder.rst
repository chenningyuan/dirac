

.. _testing_beholder:


beholder_dirac.py
-----------------

The beholder_dirac.py is an extension of the beholder library.  The beholder
library is a low-level program-independent library that provides infrastructure
for running calculations and extracting and comparing numbers against reference
outputs. For the documentation of the beholder library, please consult
https://github.com/rbast/beholder.

You have the liberty to not use beholder_dirac.py and use beholder directly but
this would mean more typing when creating new tests. Therefore we have
introduced beholder_dirac.py to make it easier to add tests that are typical
for DIRAC.


Easy example
------------

Let us look at an easy example.

First we load some standard modules and import functionality from the library
(highlighted lines).  This part is generic to all DIRAC tests.

.. code-block:: python
   :emphasize-lines: 3,4,5,6,7,8,9,10

   #!/usr/bin/env python
   
   import os
   import sys
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = '@   Elements of the electric dipole',
         to_string   = '@   anisotropy')
   f.add(from_string = ' *************************** Expectation values',
         to_string   = '    s0 = T : Expectation value')
   
   test_run.execute(['PBE0gracLB94.inp', 'GLLBsaopLBalpha.inp'], ['Ne.mol'], f)
   
   sys.exit(0)

Then we construct the filter object. It consists of two filter tasks. We can construct
as many filter objects as we like and each can consist of as many tasks as we like.

.. code-block:: python
   :emphasize-lines: 12,13,14,15,16

   #!/usr/bin/env python
   
   import os
   import sys
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = '@   Elements of the electric dipole',
         to_string   = '@   anisotropy')
   f.add(from_string = ' *************************** Expectation values',
         to_string   = '    s0 = T : Expectation value')
   
   test_run.execute(['PBE0gracLB94.inp', 'GLLBsaopLBalpha.inp'], ['Ne.mol'], f)
   
   sys.exit(0)

With test_run.execute() we really run the job and run the filter through it.
Note how we pass the filter object. If we omit to pass it, then the
calculations will be run but not verified. This is useful for multi-step jobs.
Also observe how we give a list of input files and molecule files (in this case
two input files and one molecule file). The test library will run and test all
input/molecule file combinations (in this case two).  We could have executed
them separately in two lines. This would make no difference, just more typing.

.. code-block:: python
   :emphasize-lines: 18

   #!/usr/bin/env python
   
   import os
   import sys
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = '@   Elements of the electric dipole',
         to_string   = '@   anisotropy')
   f.add(from_string = ' *************************** Expectation values',
         to_string   = '    s0 = T : Expectation value')
   
   test_run.execute(['PBE0gracLB94.inp', 'GLLBsaopLBalpha.inp'], ['Ne.mol'], f)
   
   sys.exit(0)

Finally on the last line we exit with zero (success). This is important.


Filter keywords
---------------

The filter keywords are parsed directly by the low-level test library and
therefore are documented here:
https://github.com/rbast/beholder/blob/master/doc/filter.md.


Multi-step tests
----------------

Here is an example for a multi-step test.
Note how only every second run is actually verified by passing the filter object.

.. code-block:: python
   :emphasize-lines: 20,24,28

   #!/usr/bin/env python
   
   import os
   import sys
   import shutil
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = 'Energy at final geometry is',
         num_lines   = 3,
         tolerance   = 1.0e-4)
   
   test_run.execute(['O.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc_noamfi.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc_noamfi.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   # cleanup
   os.unlink('DFCOEF')
   os.unlink('DFPROJ')
   
   sys.exit(0)

The other runs only serve to prepare files and are not checked (no filter passed as argument).

.. code-block:: python
   :emphasize-lines: 18,22,26

   #!/usr/bin/env python
   
   import os
   import sys
   import shutil
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = 'Energy at final geometry is',
         num_lines   = 3,
         tolerance   = 1.0e-4)
   
   test_run.execute(['O.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc_noamfi.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc_noamfi.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   # cleanup
   os.unlink('DFCOEF')
   os.unlink('DFPROJ')
   
   sys.exit(0)

Always clean up after yourself!

.. code-block:: python
   :emphasize-lines: 30,31,32

   #!/usr/bin/env python
   
   import os
   import sys
   import shutil
   
   sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
   from beholder import Filter
   from beholder_dirac import DiracTestRun
   
   test_run = DiracTestRun(__file__, sys.argv)
   
   f = Filter()
   f.add(from_string = 'Energy at final geometry is',
         num_lines   = 3,
         tolerance   = 1.0e-4)
   
   test_run.execute(['O.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   test_run.execute(['O.2c_iotc_noamfi.inp'], ['O.mol'], args='--get=DFCOEF')
   shutil.copy('DFCOEF', 'DFPROJ')
   test_run.execute(['H2O.2c_iotc_noamfi.inp'], ['H2O.mol'], f, args='--copy=DFPROJ')
   
   # cleanup
   os.unlink('DFCOEF')
   os.unlink('DFPROJ')
   
   sys.exit(0)
