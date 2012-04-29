Notes:Screening
===============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

We have some number of old density matrices *D*\ :sub:`*k* −
1`\ ,*D*\ :sub:`*k* − 2`\  etc, and want to write the new density matrix
as

.. figure:: /dirac/images/math/6/c/8/6c86ef3df163c5e325edbc447abb31d9.png
   :align: center
   :alt: D\_k = \\Delta D + \\sum\_{i=1}^N c\_i D\_{k-i},

   D\_k = \\Delta D + \\sum\_{i=1}^N c\_i D\_{k-i},
where Δ*D* will be sent to ``TWOFCK`` and should give efficient
screening. Currently we are minimizing the Frobenius norm (sum of
squares of elements) of Δ*D*, using one old density matrix. In the
general case we have to solve *S**c* = *y*, where

*S*\ :sub:`*i**j*`\  = (*D*\ :sub:`*k* − *i*`\ ,*D*\ :sub:`*k* − *j*`\ )
and *y*\ :sub:`*i*`\  = (*D*\ :sub:`*k*`\ ,*D*\ :sub:`*k* − *i*`\ ).

Here |(\\cdot,\\cdot)| is the inner product between the *D* "vectors".
However, the use of the Frobenius norm gives the larges weight to the
large elements of Δ*D*. This is not really what we want, because we care
more about producing as many small elements in Δ*D* as possible.
Therefore we can imagine to introduce a metric matrix *M*, so that

.. figure:: /dirac/images/math/8/7/7/877cd0b8aba2f7fcefddaa52ff414fdd.png
   :align: center
   :alt: (A,B) = \\sum\_{ij} A^\*\_{ij}M\_{ij}B\_{ij},

   (A,B) = \\sum\_{ij} A^\*\_{ij}M\_{ij}B\_{ij},
and set *M* to emphasize the elements that have a chance of being made
small. In practice one could imagine

-  Calculate Δ*D* using a unit (*M*\ :sub:`*i**j*`\  = 1) metric.
-  Find the largest elements (with some definition of largest) of Δ*D*
   and set the corresponding elements in *M* to zero.
-  Recalculate *S* using the new metric.
-  Solve for the new Δ*D*

Doing this for the water molecule SCF in test 1 shrinks the elements in
the final iteration Δ*D* by four orders of magnitude. In this case I
used five old density matrices. Just using a large number of old
densities without modifying the metric did very little. A simple test
inside Dirac shows that it is not enough to use one historical density
matrix with a modified metric. The reason is that at the end of the SCF
the change in density is almost orthogonal to the last density matrix.

It would also be possible to give larger weight to certain elements of
Δ*D*, for example those that come from different integral blocks or
centers.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Notes:Screening <http://wiki.chem.vu.nl/dirac/index.php/Notes:Screening>`_"

Views
'''''

-  `Article </dirac/index.php/Notes:Screening>`_
-  `Discussion </dirac/index.php?title=Talk:Notes:Screening&action=edit>`_
-  `View source </dirac/index.php?title=Notes:Screening&action=edit>`_
-  `History </dirac/index.php?title=Notes:Screening&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Notes:Screening>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Notes:Screening>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Notes:Screening>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Notes:Screening&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Notes:Screening&oldid=13661>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 13:11, 13 October 2008.
-  This page has been accessed 532 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |(\\cdot,\\cdot)| image:: /dirac/images/math/8/c/6/8c642c478f6aead536b9edcb76af1f42.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
