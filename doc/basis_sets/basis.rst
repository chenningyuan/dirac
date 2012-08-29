

Pick the right basis for your calculation
=========================================

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
