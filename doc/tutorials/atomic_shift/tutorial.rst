:orphan:

Atomic shift operator
=====================

DIRAC allows the definition of an atomic shift operator :ref:`HAMILTONIAN_.ASHIFT` to be added to the Hamiltonian

.. math::

   V = \sum_{A}\sum_{I\in A}\sum_{i\in I}|\varphi^A_i\rangle \varepsilon^A_I \langle \varphi^A_i|

where :math:`A` is a sum over atomic types. For each atomic types a group :math:`I` of orbitals may be
specified with an associated shift parameter :math:`\varepsilon^A_I`. In this manual we shall see
an example of its use.

Degeneracy-driven covalency: the case of CO
-------------------------------------------

We consider carbon monoxide. We run a straightforward HF calculation using the input files *CO.xyz* 

.. literalinclude:: CO.xyz

and *CO.inp* 

.. literalinclude:: CO.inp

using the command::		    

   pam --inp=CO --mol=CO.xyz

Looking at the Mulliken population analysis in the output, we find for the lowest two occupied molecular orbitals::

  * Electronic eigenvalue no.   1: -20.684805978278       (Occupation : f = 1.0000)  m_j=  1/2
  ============================================================================================
  
  * Gross populations greater than 0.01000
  
  Gross     Total   |    L A1 O  s   
  --------------------------------------
   alpha    0.9995  |      0.9992
   beta     0.0005  |      0.0000
  
  * Electronic eigenvalue no.   2: -11.369805624916       (Occupation : f = 1.0000)  m_j=  1/2
  ============================================================================================
  
  * Gross populations greater than 0.01000
  
  Gross     Total   |    L A1 C  s   
  --------------------------------------
   alpha    0.9997  |      0.9990
   beta     0.0003  |      0.0000

We see that the first and second lowest occupied MO are :math:`O1s` and :math:`C1s`, respectively, separated by a sizable energy of :math:`9.315\ E_h`.
     
Let us now see what happens if we bring the two atomic core orbitals into resonance. First we
generate the carbon atomic orbitals, using input files *C.inp* 

.. literalinclude:: C.xyz

and *C.inp* 

.. literalinclude:: C.inp

and the command::		    

   pam --inp=C --mol=C.xyz 

The checkpoint file *C_C.h5* contains the :math:`C_1` coefficients::

  h5dump -n C_C.h5 | grep C1
   dataset    /result/wavefunctions/scf/mobasis/eigenvalues_C1
   dataset    /result/wavefunctions/scf/mobasis/orbitals_C1
   
We next prepare the input file *COshift.inp*

.. literalinclude:: COshift.inp

which notably defines an atomic shift operator for carbon which builds a projector from the :math:`C1s` orbital and shifts it down by :math:`9.315\ E_h`. We run our calculation using the command::

  pam 	--inp=COshift --mol=CO.xyz --put "ac.C=AFCXXX"

The Mulliken population analysis of the two lower orbitals now read::

   * Electronic eigenvalue no.  1: -20.687541416409       (Occupation : f = 1.0000)  m_j=  1/2
   ==========================================================================================
   
   * Gross populations greater than 0.00010
   
   Gross     Total   |    L A1 C  s      L A1 O  s      A1 O  _small   B1 O  _small   B2 O  _small
   --------------------------------------------------------------------------------------------------
    alpha    0.9996  |      0.4790         0.5203         0.0001         0.0000         0.0000
    beta     0.0004  |      0.0000         0.0000         0.0000         0.0001         0.0001
   
   * Electronic eigenvalue no.  2: -20.681935936276       (Occupation : f = 1.0000)  m_j=  1/2
   ==========================================================================================
   
   * Gross populations greater than 0.00010
   
   Gross     Total   |    L A1 C  s      L A1 O  s      A1 O  _small   B1 O  _small   B2 O  _small
   --------------------------------------------------------------------------------------------------
    alpha    0.9996  |      0.5202         0.4790         0.0001         0.0000         0.0000
    beta     0.0004  |      0.0000         0.0000         0.0000         0.0001         0.0001

and we indeed see that the two lowest molecular orbitals now are almost degenerate with concomitant strong mixing of the :math:`O1s` and :math:`C1s` orbitals.

By scanning shift values in the resonance region we obtain the following plot of orbital energies of the two lowest occupied molecular orbitals


.. image:: COshift_eig.png
    :width: 600px

as well as the mixing of :math:`O1s` and :math:`C1s` orbitals in the LOMO.

.. image:: COshift.png
    :width: 600px

At this point it may be useful to consider the diagonalization of a :math:`2\times2` matrix

.. math::

   \left[\begin{array}{cc}\varepsilon_1&\delta\\\delta&\varepsilon_2\end{array}\right],

where eigenvalues are given by

.. math::

   \lambda_\pm = \frac{1}{2}\left[(\varepsilon_1+\varepsilon_2)\pm\sqrt{(\varepsilon_1-\varepsilon_2)^2 +4\delta^2}\right];\quad\Delta\lambda=\sqrt{(\varepsilon_1-\varepsilon_2)^2 +4\delta^2}

We notably see that as we bring the two starting energies into resonance, that is :math:`\varepsilon_1=\varepsilon_2`, the difference :math:`\Delta\lambda` becomes equal to the    
twice the interaction :math:`2|\delta|`. In this particular case :math:`2|\delta|=5.6\ mE_h`.

What can learn from this ? Neidig, Clark and Martin :cite:`Neidig_CoorChemRev2013` have introduced
the concept of *near-degeneracy driven covalency*, where bonding interaction arises from near-degeneracy, even in the absence of overlap. In the above example we clearly see this mechanism at work: bringing the :math:`O1s` and :math:`C1s` orbitals into resonance causes strong orbital mixing.
*However*, it would be absurd to call this bonding. The two atomic core orbitals can be disentangled by localization and we also note that their interaction is extremely weak.

QED-effects on molecular properties
-----------------------------------

As a second example we will use the atomic shift operator to provide an estimate on QED effects on
atomic properties. We shall use plumbane :math:`PbH_4` as our toy molecule.

First we need an estimate of QED-effects on the atom and therefore run a GRASP calculation 
using the input

.. literalinclude:: PbQED.ginp

From the output we extract the following information::

   ORBITAL  V.P. INTEGRAL    S.E. Flambaum   Total
   
     1S    -1.8464162D+00    8.1221155D+00   6.275699E00
     2S    -2.6815009D-01    1.2549909D+00   9.868408E−01
     2P-   -3.2180190D-02    1.7973416D-01   1.475540E−01
     2P    -1.8664568D-03    2.0271463D-01   2.008482E−01
     3S    -6.1459602D-02    2.9244620D-01   2.309866E−01
     3P-   -8.2611417D-03    3.8450656D-02   3.018951E−02
     3P    -5.0526653D-04    4.3839003D-02   4.333374E−02
     4S    -1.5604120D-02    7.4811895D-02   5.920778E−02
     3D-   -7.5017064D-06   -1.1977737D-03  −1.205275E−03
     3D    -1.4544609D-06    7.8807517D-03   7.879297E−03
     4P-   -2.0742359D-03    9.4211937D-03   7.346958E−03
     4P    -1.2756741D-04    1.0569246D-02   1.044168E−02
     5S    -3.1433537D-03    1.5156024D-02   1.201267E−02
     4D-   -2.0711814D-06   -9.3972075D-05  −9.604326E−05
     4D    -4.0004592D-07    1.8413906D-03   1.840991E−03
     5P-   -3.7559795D-04    1.7022380D-03   1.326640E−03
     5P    -2.2197734D-05    1.8210527D-03   1.798855E−03
     4F-   -6.9186631D-10   -3.5875200D-04  −3.587527E−04
     4F    -1.8748766D-10    2.8402916D-04   2.840290E−04
     5D-   -2.4901405D-07   -8.8467308D-06  −9.095745E−06
     5D    -4.6301885D-08    2.0818817D-04   2.081419E−04
     6S    -3.5592555D-04    1.7220539D-03   1.366128E−03
     6P-   -2.6181266D-05    1.1878363D-04   9.260236E−05
     6P    -1.3212900D-06    1.0828944D-04   1.069681E−04
		      
