:orphan:

===============================
Nuclear spin-rotation constants
===============================

Introduction
------------

In this tutorial we introduce the calculation of nuclear spin-rotation tensors as given in the DIRAC code, based on the theoretical
developments by I. Agustín Aucar et al. For details, you are welcome to see :cite:`Aucar_JCP2012`.

It is important to highlight that the current implementation gives theoretically correct results only for molecules in equilibrium.


The nuclear spin-rotation (SR) tensor elements of a nucleus :math:`N` are given by

.. math::

   M_{N,\alpha\beta} = M_{N,\alpha\beta}^{NU} + M_{N,\alpha\beta}^{EV} + M_{N,\alpha\beta}^{LR}

They have three terms: the first of them is independent on the electronic variables, the second one depends on an expectation value (EV) on the electronic ground state, and the latter is given by a linear response (LR) calculation.

For molecules in their equilibrium geometry, it is found that the sum of the first two terms, :math:`M_{\alpha\beta}^{nuc} + M_{\alpha\beta}^{EV}`, is equal to a new tensor element which is independent of the electronic dynamics (see Eq. 60 of :cite:`Aucar_JCP2012`),

.. math::

   M_{N,\alpha\beta}^{nuc} = M_{N,\alpha\beta}^{NU} + M_{N,\alpha\beta}^{EV}


Therefore, the electronic dependence of the nuclear spin-rotation tensor of a nucleus :math:`N` in a molecule in equilibrium is completely given by its linear response term, :math:`M_{N,\alpha\beta}^{LR}` (see Eq. 59 of :cite:`Aucar_JCP2012`).


Application to the HF molecule
------------------------------

As an example, we show a calculation of the SR constant of the fluorine nucleus at the Hydrogen fluoride molecule. The input file `spinrot.inp` is given by

.. literalinclude:: spinrot.inp

whereas the molecular input file `HF_acv4z.mol` is

.. literalinclude:: HF_acv4z.mol

The calculation is run using::

   pam --inp=spinrot --mol=HF_acv4z


As a result, SR constants are obtained at the Random Phase Approximation (RPA) level of approach. Other approaches (as DFT) are also available.

It is also interesting to note that as :ref:`HAMILTONIAN_.URKBAL` is requested in the present calculation, the results will be very close to those obtained using the RKB prescription, because the diamagnetic-like (e-p) contributions to :math:`M_{N,\alpha\beta}^{LR}` are almost zero (see Eq. 60 of :cite:`Aucar_JCP2012`).


Reading the output file
-----------------------

As the :ref:`PROPERTIES_.PRINT` value choose in the input file was 3, the results are fully detailed; they are also given in both, kHz and ppm units.

The SR constant of the fluorine nucleus, in kHz, will look like:

.. literalinclude:: SR-HF-khz

One should recall that in the case of linear molecules, as the present one, only one tensor element is printed out (the spin-rotation constant) because for these molecules the SR tensor has only two equal and non-zero elements.

As it can be seen, the total SR constant of the fluorine nucleus is given, and then it is splitted out in its three terms (NU, EV and LR). The latter is also splitted out in their (e-e) and (e-p) linear response contributions, as well as in their :math:`\mathbf{L}` and :math:`\mathbf{S}` parts. It is seen how the (e-p) contribution to the linear response term of the SR constants is almost zero.


Finally, as one should be interested in the relationship between SR constants and :ref:`PROPERTIES_.SHIELDING` (see :cite:`Aucar_JPCL2016` and :cite:`AucarChap2019`), we print also the SR constants in ppm, by multiplying each SR tensor element :math:`M_{N,\alpha\beta}` by :math:`\frac{m_p I_{\beta\gamma}}{g_N}`, where :math:`m_p` is the mass of a proton and :math:`g_N` is the g-factor of the nucleus `N`.


In this way, one obtains the following output:

.. literalinclude:: SR-HF-ppm

