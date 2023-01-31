:orphan:

X-ray absorption spectroscopy beyond the electric-dipole approximation
======================================================================

Introduction
------------

The validity of the electric dipole approximation comes into question in the X-ray regime, because at these energy scales, the wave length becomes comparable to the spatial extent of the molecule. Therefore, we will assess its validity by comparing the valence and core transitions of the barium atom. In the following, we will use two schemes to simulate effects beyond the electric dipole approximation: truncating the multipole expansion beyond zeroth order and applying the full semi--classical light--matter interaction operator. The expansion of the oscillator strength with respect to the wave vector assumes the form

.. math::

	f_{fi}=f_{fi}^{[0]}+f_{fi}^{[2]}+f_{fi}^{[4]}+\cdots,

where :math:`f_{fi}^{[0]}` corresponds to the electric dipole approximation. At arbitrary order, we have

.. math::
	f_{fi}^{[2n]}=\frac{2\omega_{fi}}{\hbar e^{2}}\sum_{m=0}^{n}(2-\delta_{m0})\epsilon_{p}\epsilon_{q}k_{j_{1}}k_{j_{2}}\ldots k_{j_{2n}}\text{Re}\bigg\{\langle f\lvert\hat{X}_{j_{1}\ldots j_{n+m};p}^{[n+m]}(\omega_{fi})\lvert i\rangle\langle f\lvert\hat{X}_{j_{n+m+1}\ldots j_{2n};q}^{[n-m]}(\omega_{fi})\lvert i\rangle^{*}\bigg\},

in which appears the multipole operator

.. math::

	\hat{X}_{j_{1}\ldots j_{n};p}^{[n]}(\omega)=\frac{1}{(n+1)!}\hat{Q}_{j_{1}\ldots j_{n};p}^{[n+1]}-\frac{i}{\omega}\frac{1}{n!}\hat{m}_{j_{1}\ldots j_{n-1};r}^{[n]}\varepsilon_{rj_{n}p},

constructed from the electric and magnetic multipole moment. The multipole operator can be expressed in two distinct forms that become equivalent in the complete basis set limit: the length and velocity representation. In the length representation, the multipole moments are represented in their conventional form


.. math::

	&\hat{Q}_{j_{1}\ldots j_{n}p}^{[n+1]}=-er_{j_{1}}\ldots r_{j_{n}}r_{p}\\&\hat{m}_{j_{1}\ldots j_{n-1};r}^{[n]}=\frac{n}{n+1}r_{j_{1}}\ldots r_{j_{n-1}}(\mathbf{r}\times\hat{\mathbf{j}})_{r};\quad\hat{\mathbf{j}}=-ec\boldsymbol{\alpha}.

In the velocity representation, the magnetic multipole remains unchanged, while the electric multipole operator reads

.. math::

	\hat{\mathcal{Q}}_{j_{1}\ldots j_{n};p}^{[n+1]}=\frac{ie}{\omega}r_{j_{1}}\ldots r_{j_{n-1}}(c\alpha_{p}r_{j_{n}}+nc\alpha_{j_{n}}r_{p}).

A nice feature of this formalism is that the isotropic averaging can be carried out analytically.	

For the full interaction, the oscillator strength can be written as

.. math::

	f_{fi}=\frac{2\omega_{fi}}{\hbar e^{2}}\big\lvert\langle f\lvert\hat{T}(\omega_{fi})\lvert i\rangle\big\lvert^{2};\quad\hat{T}(\omega)=\frac{e}{\omega}(c\boldsymbol{\alpha}\cdot\boldsymbol{\epsilon})e^{i\mathbf{k}\cdot\mathbf{r}}.

This expression contains all corrections from zero to infinity. However, it comes at a somewhat higher computational cost because the integrals depend on the frequency and the isotropic averaging needs to be carried out on a Lebedev grid. For more information, see ref. :cite:`List_JCP2020`.


Valence Transitions
-------------------

We start of by calculating the valence transitions in the barium atom, using the input `valence.inp`

.. literalinclude:: valence.inp

the xyz file of the barium atom

.. literalinclude:: Ba.xyz

and the command::

	pam --inp=valence.inp --mol=Ba.xyz --mw=700


In this file, the keyword :ref:`EXCITATION_ENERGIES_.BED` activates the calculations done with the full interaction operator. The keyword :ref:`EXCITATION_ENERGIES_.INTENS` is used for the truncated oscillator strenght, where the number below denotes the order. For each boson symmetry we will calculate one excitation. In the electric dipole approximation, this would be rather redundant, because the electric dipole operator only allows excitations that transform as either one of the Cartesian coordinates :math:`(B_{1u}, B_{2u}, B_{3u})`. The full and truncated oscillator strengths, allow more transition symmetries. We thus calculate them in this example. 

Below you can find the results for the truncated oscillator strengths::

	* Isotropic oscillator strengths (generalized velocity gauge) above threshold :1.00E-08
   - [f( 0) corresponds to the electric dipole approximation]
	================================================================================================================================

	Level    Frequency (eV) Symmetry    f(total)        f( 0)           f( 2)           f( 4)           f( 6)           f( 8) 
	--------------------------------------------------------------------------------------------------------------------------------
	    6          1.19696  B1u  1u     6.426967E-04    6.426598E-04    3.688481E-08    5.141472E-13   -4.299514E-19    1.801579E-25
	                                           Q-Q      6.426598E-04   -1.656356E-10    3.127739E-17   -5.019542E-24    7.179734E-31
	                                           M-Q      0.000000E+00    3.705044E-08   -1.988940E-14    5.751024E-21   -1.224538E-27
	                                           M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
	    7          1.19696  B2u  0u-    6.426967E-04    6.426598E-04    3.688481E-08    5.141472E-13   -4.299514E-19    1.801579E-25
	                                           Q-Q      6.426598E-04   -1.656356E-10    3.127739E-17   -5.019542E-24    7.179733E-31
	                                           M-Q      0.000000E+00    3.705044E-08   -1.988940E-14    5.751024E-21   -1.224538E-27
	                                           M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
	    8          1.19696  B3u  0u-    6.426967E-04    6.426598E-04    3.688481E-08    5.141472E-13   -4.299514E-19    1.801579E-25
	                                           Q-Q      6.426598E-04   -1.656356E-10    3.127739E-17   -5.019542E-24    7.179733E-31
	                                           M-Q      0.000000E+00    3.705044E-08   -1.988940E-14    5.751024E-21   -1.224538E-27
	                                           M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
	--------------------------------------------------------------------------------------------------------------------------------
	Sum of oscillator strengths (general length) :  1.92809E-03

	* Isotropic oscillator strengths (generalized length gauge) above threshold :1.00E-08
	   - [f( 0) corresponds to the electric dipole approximation]
	================================================================================================================================

	Level    Frequency (eV) Symmetry    f(total)        f( 0)           f( 2)           f( 4)           f( 6)           f( 8) 
	--------------------------------------------------------------------------------------------------------------------------------
	    6          1.19696  B1u  1u     6.330127E-04    6.329761E-04    3.660752E-08    5.143085E-13   -4.300219E-19    1.801837E-25
	                                           Q-Q      6.329761E-04   -1.627207E-10    2.997773E-17   -4.615767E-24    6.312181E-31
	                                           M-Q      0.000000E+00    3.677024E-08   -1.972681E-14    5.680071E-21   -1.198602E-27
 	                                          M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
	    7          1.19696  B2u  0u-    6.330127E-04    6.329761E-04    3.660752E-08    5.143085E-13   -4.300219E-19    1.801837E-25
	                                           Q-Q      6.329761E-04   -1.627207E-10    2.997772E-17   -4.615767E-24    6.312180E-31
	                                           M-Q      0.000000E+00    3.677024E-08   -1.972681E-14    5.680071E-21   -1.198602E-27
 	                                          M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
 	   8          1.19696  B3u  0u-    6.330127E-04    6.329761E-04    3.660752E-08    5.143085E-13   -4.300219E-19    1.801837E-25
	                                           Q-Q      6.329761E-04   -1.627207E-10    2.997772E-17   -4.615767E-24    6.312180E-31
 	                                          M-Q      0.000000E+00    3.677024E-08   -1.972681E-14    5.680071E-21   -1.198602E-27
	                                           M-M      0.000000E+00    0.000000E+00    5.340053E-13   -4.356974E-19    1.813817E-25
	--------------------------------------------------------------------------------------------------------------------------------
	Sum of oscillator strengths (general length) :  1.89904E-03


First of all, note that only the dipole allowed transitions are visible. The code did calculate the other transitions, but they fell below the threshold of :math:`10^{-8}` and are hence not displayed. Furthermore, the numbers clearly indicate that the multipole expansion is already converged at zeroth order, i.e. the electric dipole approximation. This was to be expected, since it is known that the electric dipole approximation produces good results in the UV-vis regime. 

The results from the full interaction confirm the validity of the electric dipole approximation as well: only the dipole allowed transitions were large enough to be printed and their values agree with :math:`f^{[0]}`::

	* Isotropic oscillator strengths (full light-matter interaction) above threshold :1.00E-08
	==========================================================================================

	Level    Frequency (eV) Symmetry    f(total)        
	------------------------------------------------------------------------------------------
	    6          1.19696  B1u  1u     6.426967E-04
	    7          1.19696  B2u  0u-    6.426967E-04
	    8          1.19696  B3u  0u-    6.426967E-04
	------------------------------------------------------------------------------------------
	Sum of oscillator strengths (full light-matter interaction) :  1.92809E-03


Core Transitions
----------------

We will now proceed to calculate the core transitions using the input `core.inp`

.. literalinclude:: core.inp

the same xyz file and the command::

	pam --inp=core.inp --mol=Ba.xyz

Note that this input is nearly identical to the previous one. The main difference is the :ref:`EXCITATION_ENERGIES_.OCCUP` keyword. Once activated, this keyword invokes the restricted-excitation-window approximation :ref:`stener2003time` :ref:`south20164`. Without this approximation, core excitations can only be calculated by first computing all higher excitations, since the diagonalization routine in DIRAC computes the excitation energies from lowest to highest. The restricted-excitation-window approximation allows us to fix the occupied orbital involved in the excitation. In this example, we set this orbital to the :math:`1s_{1/2}` orbital, hence calculating the Barium K-edge.

Below you can find the truncated oscillator strengths::

	* Isotropic oscillator strengths (generalized velocity gauge) above threshold :1.00E-08
	   - [f( 0) corresponds to the electric dipole approximation]
	================================================================================================================================

	Level    Frequency (eV) Symmetry    f(total)        f( 0)           f( 2)           f( 4)           f( 6)           f( 8) 
	--------------------------------------------------------------------------------------------------------------------------------
	    1      37334.75964  B3g  1g    -3.968737E-06    0.000000E+00    1.745343E-09   -2.506329E-09    1.488143E-08   -3.982858E-06
	                                           Q-Q      0.000000E+00    1.390509E-21    1.553802E-23    7.343148E-21   -1.065761E-18
	                                           M-Q      0.000000E+00    0.000000E+00    4.873376E-22   -1.212851E-19    1.136316E-17
	                                           M-M      0.000000E+00    1.745343E-09   -2.506329E-09    1.488143E-08   -3.982858E-06
	    2      37334.75964  B2g  1g    -3.968712E-06    0.000000E+00    1.745343E-09   -2.506326E-09    1.488107E-08   -3.982832E-06
	                                           Q-Q      0.000000E+00    5.887396E-22    2.138116E-22    2.441124E-21   -9.801645E-19
	                                           M-Q      0.000000E+00    0.000000E+00    1.203034E-21   -2.030584E-19    1.563374E-17
	                                           M-M      0.000000E+00    1.745343E-09   -2.506326E-09    1.488107E-08   -3.982832E-06
	    3      37334.75964  B1g  1g    -3.968687E-06    0.000000E+00    1.745343E-09   -2.506324E-09    1.488075E-08   -3.982806E-06
	                                           Q-Q      0.000000E+00    3.309524E-27   -4.389399E-25    2.141423E-23    1.557558E-22
	                                           M-Q      0.000000E+00    0.000000E+00    3.374693E-22    2.777935E-20    6.783862E-19
	                                           M-M      0.000000E+00    1.745343E-09   -2.506324E-09    1.488075E-08   -3.982806E-06
 	   4      37334.76073  Ag   2g     1.262812E-04    0.000000E+00    9.737751E-07   -2.670012E-08   -1.008622E-06    1.263427E-04
	                                           Q-Q      0.000000E+00    9.737751E-07   -6.744802E-08   -7.083675E-07    3.108814E-05
	                                           M-Q      0.000000E+00    0.000000E+00    4.074789E-08   -3.006805E-07    9.526082E-05
	                                           M-M      0.000000E+00    0.000000E+00    0.000000E+00    4.262768E-10   -6.261504E-09
	    6      37336.26040  B2u  0u-    1.589508E+00    2.163748E-05    1.184342E-06    5.310252E-06   -4.229244E-03    1.593709E+00
	                                           Q-Q      2.163748E-05   -2.459370E-07    8.602505E-07   -9.749223E-05   -2.545372E-02
	                                           M-Q      0.000000E+00    1.430279E-06    4.426366E-06   -4.131899E-03    1.619299E+00
	                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465643E-07   -1.363362E-04
	    7      37336.26040  B1u  1u     1.589507E+00    2.163748E-05    1.184342E-06    5.310240E-06   -4.229240E-03    1.593708E+00
	                                           Q-Q      2.163748E-05   -2.459370E-07    8.602463E-07   -9.749113E-05   -2.545392E-02
	                                           M-Q      0.000000E+00    1.430279E-06    4.426358E-06   -4.131895E-03    1.619298E+00
	                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465641E-07   -1.363360E-04
	    8      37336.26040  B3u  0u-    1.589505E+00    2.163748E-05    1.184342E-06    5.310230E-06   -4.229234E-03    1.593706E+00
	                                           Q-Q      2.163748E-05   -2.459370E-07    8.602443E-07   -9.748983E-05   -2.545428E-02
	                                           M-Q      0.000000E+00    1.430280E-06    4.426350E-06   -4.131891E-03    1.619297E+00
	                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465638E-07   -1.363359E-04
	--------------------------------------------------------------------------------------------------------------------------------
	Sum of oscillator strengths (general length) :  4.76863E+00

	* Isotropic oscillator strengths (generalized length gauge) above threshold :1.00E-08
   - [f( 0) corresponds to the electric dipole approximation]
	================================================================================================================================

	Level    Frequency (eV) Symmetry    f(total)        f( 0)           f( 2)           f( 4)           f( 6)           f( 8) 
	--------------------------------------------------------------------------------------------------------------------------------
    1      37334.75964  B3g  1g    -3.968737E-06    0.000000E+00    1.745343E-09   -2.506329E-09    1.488143E-08   -3.982858E-06
                                           Q-Q      0.000000E+00    1.310439E-21   -6.202397E-21    1.536998E-18   -1.002294E-16
                                           M-Q      0.000000E+00    0.000000E+00    7.330683E-22   -9.185172E-20    9.616797E-18
                                           M-M      0.000000E+00    1.745343E-09   -2.506329E-09    1.488143E-08   -3.982858E-06
    2      37334.75964  B2g  1g    -3.968712E-06    0.000000E+00    1.745343E-09   -2.506326E-09    1.488107E-08   -3.982832E-06
                                           Q-Q      0.000000E+00    5.601652E-22   -2.836271E-21    6.654874E-19   -4.324261E-17
                                           M-Q      0.000000E+00    7.722475E-32    1.199779E-21   -1.747560E-19    1.533937E-17
                                           M-M      0.000000E+00    1.745343E-09   -2.506326E-09    1.488107E-08   -3.982832E-06
    3      37334.75964  B1g  1g    -3.968687E-06    0.000000E+00    1.745343E-09   -2.506324E-09    1.488075E-08   -3.982806E-06
                                           Q-Q      0.000000E+00    2.098863E-27   -4.366621E-25    5.233625E-23   -4.489705E-21
                                           M-Q      0.000000E+00    0.000000E+00    1.446503E-24   -3.713084E-22    4.904970E-20
                                           M-M      0.000000E+00    1.745343E-09   -2.506324E-09    1.488075E-08   -3.982806E-06
    4      37334.76073  Ag   2g    -7.052296E-02    0.000000E+00    9.179494E-07   -4.608196E-06    1.103012E-03   -7.162228E-02
                                           Q-Q      0.000000E+00    9.179494E-07   -4.647758E-06    1.103403E-03   -7.173915E-02
                                           M-Q      0.000000E+00    0.000000E+00    3.956263E-08   -3.907210E-07    1.168808E-04
                                           M-M      0.000000E+00    0.000000E+00    0.000000E+00    4.262768E-10   -6.261504E-09
    6      37336.26040  B2u  0u-   -5.065597E+02    2.165894E-05    3.399533E-05   -1.473260E-02    3.376557E+00   -5.099215E+02
                                           Q-Q      2.165894E-05    3.256434E-05   -1.473814E-02    3.381175E+00   -5.116489E+02
                                           M-Q      0.000000E+00    1.430989E-06    5.512443E-06   -4.617888E-03    1.727525E+00
                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465643E-07   -1.363362E-04
    7      37336.26040  B1u  1u    -5.065597E+02    2.165894E-05    3.399533E-05   -1.473260E-02    3.376557E+00   -5.099216E+02
                                           Q-Q      2.165894E-05    3.256434E-05   -1.473814E-02    3.381175E+00   -5.116490E+02
                                           M-Q      0.000000E+00    1.430989E-06    5.512434E-06   -4.617885E-03    1.727525E+00
                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465641E-07   -1.363360E-04
    8      37336.26040  B3u  0u-   -5.065597E+02    2.165894E-05    3.399533E-05   -1.473260E-02    3.376557E+00   -5.099215E+02
                                           Q-Q      2.165894E-05    3.256434E-05   -1.473814E-02    3.381175E+00   -5.116489E+02
                                           M-Q      0.000000E+00    1.430989E-06    5.512427E-06   -4.617880E-03    1.727523E+00
                                           M-M      0.000000E+00    0.000000E+00    2.363606E-08    1.465638E-07   -1.363359E-04
	--------------------------------------------------------------------------------------------------------------------------------
	Sum of oscillator strengths (general length) : -1.51975E+03

Contrary to previous example, these oscillator strengths do not appear to converge to a certain value. It appears that the magnitude of the oscillator strengths monotonically increases, which implies that the accumulated oscillator strength diverges. In ref. :cite:`List_JCP2020`, it is shown that the oscillator strengths do not diverge, but rather converge extremely slowly. 

Due to these artifacts, it is much more preferable to apply the full interaction::

	* Isotropic oscillator strengths (full light-matter interaction) above threshold :1.00E-08
	==========================================================================================

	Level    Frequency (eV) Symmetry    f(total)        
	------------------------------------------------------------------------------------------
	    4      37334.76073  Ag   2g     9.407196E-07
	    6      37336.26040  B2u  0u-    2.280511E-05
	    7      37336.26040  B1u  1u     2.280511E-05
	    8      37336.26040  B3u  0u-    2.280511E-05
	------------------------------------------------------------------------------------------
	Sum of oscillator strengths (full light-matter interaction) :  6.93566E-05

The numbers calculated with this method do not suffer from aforementioned problems. However, if we compare the full oscillator strengths with the lowest non-zero truncated oscillator strength, we find that the two merely differ by 5\%. This small difference can be understood by the compactness of the :math:`1s_{1/2}` orbital. The electric dipole approximation breaks down if the limit :math:`kr \ll 1` holds. Even though the wave length is short in this example, the spatial extent of the :math:`1s_{1/2}` orbital is small, hence reducing non-dipolar effects. Therefore, compared to the electric dipole approximation, the full interaction gives modest corrections to dipole-allowed transitions, but it is essential to include dipole-forbidden transitions. However, it remains an open question whether non-dipolar effects play a more important role in extended systems. I leave it to you to find it out. 
