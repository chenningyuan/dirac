Basis sets for relativistic calculations
========================================

Gaussian Type Orbitals (GTOs)
-----------------------------

*Cartesian Gaussians* are defined as

.. math::

   G^{\alpha}_{ijk}=G^{\alpha}_i(x)G^{\alpha}_j(y)G^{\alpha}_k(z)

with

.. math::

   G^{\alpha}(x)=\left(\frac{2\alpha}{\pi}\right)^{1/4}\sqrt{\frac{\left(4\alpha\right)^i}{\left(2i-1\right)!!}}
                  x^i\exp\left[-\alpha x^2\right]

Alternatively one may express a Cartesian Gaussian as

.. math::

   G^{\alpha}_{ijk}=N^{\alpha}_{ijk}x^iy^jz^k\exp\left[-\alpha r^2\right];\quad N^{\alpha}_{ijk}=\left(\frac{2\alpha}{\pi}\right)^{3/4}\sqrt{\frac{2^l}{F){ijk}}}\left(\sqrt{2\alpha}\right)^l

The factor 

.. math::

   F_{ijk}=\left(2i-1\right)!!\left(2j-1\right)!!\left(2k-1\right)!!

shows that the Cartesian components of a given shell may have different normalization constants. In the HERMIT 
integral code a single normalization is chosen for each shell by ignoring :math:`F_{ijk}` meaning that Cartesian
Gaussians are normalized to

.. math::

   \left<G^{\alpha}_{ijk}|G^{\alpha}_{ijk}\right>=F_{ijk}

In practice this means that s- and p-functions are normalized to one. So are :math:`d110`, :math:`d101` and 
:math:`d011`, whereas :math:`d200`, :math:`d020` and :math:`d002` are normalized to three. :math:`f111` is
normalized to one, :math:`f210`, :math:`f201`, :math:`f120`, :math:`f102` and :math:`f012` are normalized to 
three, and :math:`f300`, :math:`f030` and :math:`f003` are normalized to 15.

   
