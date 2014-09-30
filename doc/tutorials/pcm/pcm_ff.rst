Calculation of polarizabilities in solution: finite field approach
==================================================================

`Five-point stencil <http://en.wikipedia.org/wiki/Five-point_stencil>`_:

.. math::
   
   f'(x) &\simeq \frac{-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)}{12h} \\
   f''(x) &\simeq \frac{-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) + f(x-2h)}{12h^2}
