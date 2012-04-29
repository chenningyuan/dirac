How to build math libraries for 64-bit integers
===============================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Before you continue please read the
`FAQ </dirac/index.php/FAQ#I_am_about_to_install_DIRAC_on_a_64-bit_machine.2C_should_I_compile_using_64-bit_integers_or_32-bit_integers.3F>`_
and verify whether you really need 64-bit integers.

It is extra work and perhaps not needed.

If you decide to build DIRAC for 64-bit integers, you have to link to
64-bit integer aware math library/libraries providing both BLAS and
LAPACK.

The following things can happen:

-  Your math library is compiled for 64-bit integers and DIRAC will link
   and run properly.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   not link.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   link but stop at runtime.
-  Your math library is not compiled for 64-bit integers and DIRAC will
   link and not stop at runtime because you deactivated the self-test
   and may produce wrong numbers.

For cases 2-4 you may have to either:

-  Verify linking.
-  Compile your own math library.
-  Use DIRAC's internal math implementation (slow).
-  Go back to 32-bit integers.

MKL
---

MKL provides bindings for 64-bit integers and normally DIRAC should
correctly link if you have specified MATH\_ROOT correctly (see the
installation guide).

Atlas
-----

It is not completely trivial but doable to compile your own Atlas
library with 64-bit integer support.

The difficulty is that you need to compile both LAPACK and Atlas.

To do this we recommend to follow the nice example: `Installing ATLAS
with full LAPACK on
Linux/AMD64 <http://math-atlas.sourceforge.net/atlas_install/atlas_install.html#SECTION00090000000000000000>`_.

