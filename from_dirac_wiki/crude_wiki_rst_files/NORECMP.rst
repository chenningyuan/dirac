NORECMP
=======

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSORT -- Sorting module <#.26CCSORT_--_Sorting_module>`_

   -  `1.1 NORECMP <#NORECMP>`_

&CCSORT -- Sorting module
=========================

NORECMP
-------

For the open shell calculation, turn off this flag. Use it with USEOE
option.

::

    &RELCCSD NELEC=4,3 TIMING=T, IPRNT=1, DOSORT=T, DOENER=T &END
    &CCSORT USEOE=F, NORECMP=F &END
    &CCENER DOMP2=T, DOCCSD=T, DOCCSDT=F &END

