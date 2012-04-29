Calculating and plotting the magnetically induced current density
=================================================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

`|image0| </dirac/index.php/Image:Dirac.png>`_

**The functionality described on this page works with DIRAC10 but not
with the latest development version and will not work with DIRAC11.**

In this tutorial we wish to calculate and plot the magnetically induced
current density in cyclobutadiene using Hartree-Fock, Dirac-Coulomb
Hamiltonian, and a common gauge origin. The calculation can be performed
on a laptop in a minute or two. At the same time we will learn how to
restart response calculations from converged SCF wave function
calculations.

We will use the following molecule file and call it **cbd.mol**:

::

    INTGRL


    C   2
           6.     4
    C     -1.36000 -1.36000  0.0
    C     -1.36000  1.36000  0.0
    C      1.36000 -1.36000  0.0
    C      1.36000  1.36000  0.0
    LARGE BASIS cc-pVDZ
           1.     4
    H     -2.83000 -2.83000  0.0
    H     -2.83000  2.83000  0.0
    H      2.83000 -2.83000  0.0
    H      2.83000  2.83000  0.0
    LARGE BASIS cc-pVDZ
    FINISH

The molecule is oriented in the xy-plane and we will place the magnetic
field vector along the perpendicular z-axis. Note that the here employed
(square) structure of cyclobutadiene does not correspond to the
equilibrium structure and that the basis set (cc-pVDZ) is not ideal
either. In real life calculations the basis set needs to be well
calibrated.

First we need to optimize the SCF wave function. We will use the
following input file and call it **scf.inp**:

::

    **DIRAC
    .WAVE FUNCTION
    **HAMILTONIAN
    .LVCORR
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    *END OF

And we can run DIRAC with the following run script:

::

    #!/bin/bash

      pam -nobackup -noarch        -outcmo                               -mem 100 cbd.mol scf.inp
    # pam -nobackup -noarch -incmo -outcmo                -get "PAMXVC"  -mem 100 cbd.mol response.inp
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_total.inp
    # mv PLOT_2D j_total.2d
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_para.inp
    # mv PLOT_2D j_para.2d
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_dia.inp
    # mv PLOT_2D j_dia.2d

    exit 0

Make sure that the **pam** script is in your $PATH. Note that we save
the optimized MO coefficients file (DFCOEF) using the **-outcmo** flag.

The calculation should finish in about a minute and at the end you
should see an output file (scf\_cbd.out) and the DFCOEF file in your
directory.

Now we are ready for the second step - the response calculation with the
following input file that we can call **response.inp**:

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    **HAMILTONIAN
    .LVCORR
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    **PROPERTIES
    *LINEAR
    .A OPERATOR
     'B_z oper'
     ZAVECTOR
     'YDIPLEN'
     'XDIPLEN'
     COMFACTOR
     -68.517999904721
    .B OPERATOR
     'B_z oper'
     ZAVECTOR
     'YDIPLEN'
     'XDIPLEN'
     COMFACTOR
     -68.517999904721
    .ANALYZE
    *END OF

We will copy the DFCOEF file to the working directory using the
**-incmo** flag and at the end save the optimized response vector
(PAMXVC) using **-get "PAMXVC":**

::

    #!/bin/bash

    # pam -nobackup -noarch        -outcmo                               -mem 100 cbd.mol scf.inp
      pam -nobackup -noarch -incmo -outcmo                -get "PAMXVC"  -mem 100 cbd.mol response.inp
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_total.inp
    # mv PLOT_2D j_total.2d
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_para.inp
    # mv PLOT_2D j_para.2d
    #
    # pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_dia.inp
    # mv PLOT_2D j_dia.2d

    exit 0

You should check that you obtain the following response functions:

::

    grep "at freq" response_cbd.out
    response_cbd.out:<< 1, 1>>: -0.33026893E+02 a.u.,     E--contribution at frequency    0.0000000000 a.u.
    response_cbd.out:<< 1, 1>>:  0.37735752E+02 a.u.,     P--contribution at frequency    0.0000000000 a.u.
    response_cbd.out:    4.708858878995 a.u. at frequency      0.00000000 a.u.    5.72E-06   (converged)

These are the paramagnetic, diamagnetic, and the total zz
magnetizability tensor element, respectively.

Fine. Now we have both the DFCOEF and the PAMXVC file and are ready to
calculate the corresponding induced current density. We will do that in
2D planes (xy, xz, and yz) as well as in 3 chosen points and we will use
3 input files to calculate the total, the paramagnetic, and the
diamagnetic current density.

**current\_total.inp**:

::

    **DIRAC
    # it is important to comment
    # out .WAVE FUNCTION
    # otherwise the MO coefficients
    # may change during SCF step
    #.WAVE FUNCTION
    **HAMILTONIAN
    .LVCORR
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    **VISUAL
    .2D
    .QUANTITIES
     2
    PAMXVC 1 5 0
    PAMXVC 2 5 0
    .2DADD
     4.0
    .2DSTEP
     0.2
    .2DORIGIN
     0.0 0.0 0.0
    .LIST
     3
     0.5 0.5 0.0
     1.0 1.0 0.0
     2.0 2.0 0.0
    *END OF

**current\_para.inp**:

::

    **DIRAC
    # it is important to comment
    # out .WAVE FUNCTION
    # otherwise the MO coefficients
    # may change during SCF step
    #.WAVE FUNCTION
    **HAMILTONIAN
    .LVCORR
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    **VISUAL
    .2D
    .QUANTITIES
     1
    PAMXVC 1 5 0
    .2DADD
     4.0
    .2DSTEP
     0.2
    .2DORIGIN
     0.0 0.0 0.0
    .LIST
     3
     0.5 0.5 0.0
     1.0 1.0 0.0
     2.0 2.0 0.0
    *END OF

**current\_dia.inp**:

::

    **DIRAC
    # it is important to comment
    # out .WAVE FUNCTION
    # otherwise the MO coefficients
    # may change during SCF step
    #.WAVE FUNCTION
    **HAMILTONIAN
    .LVCORR
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    **VISUAL
    .2D
    .QUANTITIES
     1
    PAMXVC 2 5 0
    .2DADD
     4.0
    .2DSTEP
     0.2
    .2DORIGIN
     0.0 0.0 0.0
    .LIST
     3
     0.5 0.5 0.0
     1.0 1.0 0.0
     2.0 2.0 0.0
    *END OF

The 3 files differ in the **.QUANTITIES** keyword. The differing numbers
1 and/or 2 correspond to records on the PAMXVC file. The first record is
the paramagnetic response (more correctly the positive-positive energy
response), the second record the diamagnetic response (more correctly
the positive-negative energy response).

Let us run the calculations with the following run script:

::

    #!/bin/bash

    # pam -nobackup -noarch        -outcmo                               -mem 100 cbd.mol scf.inp
    # pam -nobackup -noarch -incmo -outcmo                -get "PAMXVC"  -mem 100 cbd.mol response.inp
     
      pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_total.inp
      mv PLOT_2D j_total.2d
     
      pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_para.inp
      mv PLOT_2D j_para.2d
     
      pam -nobackup -noarch -incmo         -copy "PAMXVC" -get "PLOT_2D" -mem 100 cbd.mol current_dia.inp
      mv PLOT_2D j_dia.2d

    exit 0

Now we have j\_total.2d, j\_para.2d, and j\_dia.2d in our directory. In
tutorial we are only interested in the xy-plane which we will cut out
using sed:

::

    sed '1,/xy/d;/plane/,/$/d' j_total.2d > j_total.2d_xy
    sed '1,/xy/d;/plane/,/$/d' j_para.2d  > j_para.2d_xy
    sed '1,/xy/d;/plane/,/$/d' j_dia.2d   > j_dia.2d_xy

Finally we can plot the current densities using gnuplot for instance:

::

    #!/usr/bin/gnuplot

    set term postscript enhanced

    set label 'C' at -1.92, 0.00
    set label 'C' at  1.92, 0.00
    set label 'C' at  0.00,-1.92
    set label 'C' at  0.00, 1.92
    set label 'H' at -4.00, 0.00
    set label 'H' at  4.00, 0.00
    set label 'H' at  0.00,-4.00
    set label 'H' at  0.00, 4.00

    set size square
    set xrange [-5:5]
    set yrange [-5:5]

    scale=5.0

    set output "j_total.ps"
    plot "j_total.2d_xy" using 1:2:($5*scale):($6*scale) not w vector

    set output "j_para.ps"
    plot "j_para.2d_xy" using 1:2:($5*scale):($6*scale) not w vector

    set output "j_dia.ps"
    plot "j_dia.2d_xy" using 1:2:($5*scale):($6*scale) not w vector

The resulting figures are not exactly publication quality but they give
us a already a quick idea of the induced current density. For more
advanced plots we recommend for instance
`PyNGL <http://www.pyngl.ucar.edu/>`_, but there are many nice
alternatives.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Calculating\_and\_plotting\_the\_magnetically\_induced\_current\_density <http://wiki.chem.vu.nl/dirac/index.php/Calculating_and_plotting_the_magnetically_induced_current_density>`_"

Views
'''''

-  `Article </dirac/index.php/Calculating_and_plotting_the_magnetically_induced_current_density>`_
-  `Discussion </dirac/index.php?title=Talk:Calculating_and_plotting_the_magnetically_induced_current_density&action=edit>`_
-  `View
   source </dirac/index.php?title=Calculating_and_plotting_the_magnetically_induced_current_density&action=edit>`_
-  `History </dirac/index.php?title=Calculating_and_plotting_the_magnetically_induced_current_density&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Calculating_and_plotting_the_magnetically_induced_current_density>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Calculating_and_plotting_the_magnetically_induced_current_density>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Calculating_and_plotting_the_magnetically_induced_current_density>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Calculating_and_plotting_the_magnetically_induced_current_density&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Calculating_and_plotting_the_magnetically_induced_current_density&oldid=26311>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 22:25, 20 January 2011.
-  This page has been accessed 961 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
