VIBCAL
======

VIBCAL is a utility program allowing the calculation of vibrationally averaged properties.

Theory
------

Suppose that we have optimized the geometry of our molecule and have carried out a harmonic
vibrational analysis. For a selected normal mode we expand the electronic energy :math:`E` 
as function of the associated normal coordinate :math:`q` about the equilibrium position 
:math:`q=0`. 

.. math::

   E(q) = E^{[0]} + E^{[1]}q + \frac{1}{2}E^{[2]}q^2 + \frac{1}{6}E^{[3]}q^3 + \ldots;\quad
          E^{[n]} = \left.\frac{\partial^n E}{\partial q^n}\right|_{q=0}

We note that :math:`q=0`  implies :math:`E^{[1]}=0` and that the second derivative :math:`E^{[2]}` corresponds to the force constant :math:`k`. We furthermore expand some property :math:`P\equiv P(q)` in the same manner. 

From perturbation theory we now obtain an expression for the expectation value of the property :math:`P` for vibrational level :math:`\nu` of the 
selected normal mode

.. math::

   P_{\nu}=\left<\nu\left|P(q)\right|\nu\right>_q = P^{[0]} + \frac{1}{2}\left(\frac{\hbar}{\mu\omega_e}\right)\left(\nu+\frac{1}{2}\right)
      -\frac{1}{2}P^{[1]}E^{[3]}\left(\frac{\hbar}{\mu\omega_e}\right)^2\frac{\left(\nu+\frac{1}{2}\right)}{\hbar\omega_e}+\ldots

..      -\frac{1}{16}P^{[2]}E^{[4]}\left(\frac{\hbar}{\mu\omega_e}\right)^3\frac{n^2+n+\frac{1}{2}}{\hbar\omega_e}

where appears the reduced mass :math:`\mu` and harmonic frequency :math:`\omega_e=\sqrt{\frac{k}{\mu}}` of the selected normal mode.

Experimentally one can observe the change in this expectation value between two vibrational levels, given by

.. math::

   P_{\nu'} - P_{\nu} = \frac{1}{2}\frac{\hbar}{\mu\omega_e}\left[P^{[2]}-\frac{1}{\mu\omega_e^2}P^{[1]}E^{[3]}\right]\Delta\nu;\quad \Delta\nu=\nu' - \nu

We note that the purely electronic contribution :math:`P^{[0]}` does not contribute to the vibrational shift. Furthermore, we may note that for a purely
harmonic mode (:math:`E^{[3]}=0`) and a property purely linear in the normal coordinate (:math:`P^{[2]}=0`) the vibrational shift is strictly zero.
