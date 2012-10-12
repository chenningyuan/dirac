

Pick the right basis for your calculation
=========================================

The development of basis sets suitable to use in relativistic calculations
reflects the relatively lateness of the field's development. Because the more
consistent efforts in method development started at about the mid 1980's, it
wasn't until well into the late 1990's that the pioneering works of the early
and mid 1990's were substantially complemented and improved upon.

This situations has dramatically improved in recent years, notably with the
work of K. G. Dyall, and Dirac users are strongly advised to use Dyall's basis
sets whenever they are available. These sets follow roughly the
"correlation-consisten" philosophy introduced by Dunning and coworkers
[Dunning1989]_, so they already contain polarization functions.

At the time of writing, these are avaiable from Ga onwards, and a detailed
description can be found at the dirac website, under the "Basis Set Repository"
link`[1] <http://dirac.chem.sdu.dk/basisarchives/dyall/index.html>`_.  For
convenience, Dirac collects all of Dyall's available basis sets (those
published and areadly in distributable form, as the ones listed in the
repository below, as well as those in preparation) in basis set files. One must
keep in mind that currently these files support only the use of fully
uncontracted basis sets.

The following files are present in Dirac's distribution:

dyall.v2z contains the double-zeta basis dyall.v3z contains the triple-zeta
basis dyall.v4z contains the quadruple-zeta basis dyall.av2z contains the
double-zeta basis (as in dyall.v2z) plus diffuse functions (p-block only)
dyall.av3z contains the triple-zeta basis (as in dyall.v3z) plus diffuse
functions (p-block only) dyall.av4z contains the quadruple-zeta basis (as in
dyall.v4z) plus diffuse functions (p-block only)

For light elements (say up to Argon) it is advisable to continue using standard
non-relativistic basis set, such as the correlation-consistent sets of Dunning
and coworkers. It is advisable that, in order to have a balanced description
when light and heavy elements are present, that one uses either contracted or
uncontracted sets thoughout.

Apart from Dyall's sets, one can choose several different basis sets based upon
geometric progressions of exponents. One such set is that of K. Faegri, also
available in the basis set library, but with the drawback that the user needs
to extended it by adding polarization functions.

The DIRAC distribution shares a large library of standard
non-relativistic and scalar-relativistic basis sets with the
`Dalton <http://www.kjemi.uio.no/software/dalton>`_ program.
These basis sets can be found
in the directory **basis\_dalton** of the DIRAC distribution.

These basis sets are not all suitable for relativistic calculations,
especially not for the heavier elements.
Basis sets developed for full relativistic calculations (including
spin-orbit coupling) can be found in the directory **basis**.

As Ken Dyall extends the archive actively, the latest of these basis
sets may not be available in the **basis** directory of your DIRAC
distribution.

See the `Dyall basis set
repository <http://dirac.chem.sdu.dk/basisarchives/dyall/index.html>`_
for the latest updates by Ken Dyall and the appropriate basis set
references.

We recommend that you use Dyall basis sets
(see, for instance http://dirac.chem.sdu.dk/basisarchives/dyall/index.html)
whenever they are available for the
elements of interest. In order to make that usage as convenient as possible,
the following files, containing all sets currently available at the
URL above (published or to be published), are made available:

+------------------+-------------+----------------+
|     quality      |  valence    |  core-valence  |
+==================+=============+================+
|     double-zeta  |  dyall.v2z  |  dyall.cv2z    |
+------------------+-------------+----------------+
|     triple-zeta  |  dyall.v3z  |  dyall.cv3z    |
+------------------+-------------+----------------+
|  quadruple-zeta  |  dyall.v4z  |  dyall.cv4z    |
+------------------+-------------+----------------+

While the division in "valence" and "core-valence" can be at times
not so clear-cut as for lighter elements, the option was made to
stick to the usual jargon of non-relativistic theory, particularly
in relation to the "correlation-consistent" family of basis sets.

Therefore, one should keep in mind that the major difference between
the "valence" and "core-valence" sets above is that for the latter
all the exponents contained in the archives found on the URL above
are present, whereas for "valence" basis usually all but
the exponents related to core correlation are included. Users are
encouraged to look in these files and in the original archives
published at Theor. Chem. Acc. (as, for instance, for the 5p TZ
basis, http://dirac.chem.sdu.dk/basisarchives/dyall/5p_tz_archive)
to get a feel for what is included in each case.

In the case of p-block elements, additional files containing the
corresponding diffuse functions, in addition to the valence and
core-valence basis sets mentioned above, are also provided.

+------------------+--------------+-----------------+
|     quality      |  valence     |  core-valence   |
+==================+==============+=================+
|     double-zeta  |  dyall.av2z  |  dyall.acv2z    |
+------------------+--------------+-----------------+
|     triple-zeta  |  dyall.av3z  |  dyall.acv3z    |
+------------------+--------------+-----------------+
|  quadruple-zeta  |  dyall.av4z  |  dyall.acv4z    |
+------------------+--------------+-----------------+

In case of errors or omissions on any of the files in this directory,
users are kindly asked to contact the authors of DIRAC.
