Meeting2008
===========

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Dirac meeting Odense June 9-13
   2008 <#Dirac_meeting_Odense_June_9-13_2008>`_

   -  `1.1 Schedule <#Schedule>`_

      -  `1.1.1 Monday 9/6 <#Monday_9.2F6>`_

         -  `1.1.1.1 Theme: Tools for distributed
            development <#Theme:_Tools_for_distributed_development>`_

      -  `1.1.2 Tuesday 10/6 <#Tuesday_10.2F6>`_

         -  `1.1.2.1 Theme: Modularity <#Theme:_Modularity>`_

      -  `1.1.3 Wednesday 11/6 <#Wednesday_11.2F6>`_

         -  `1.1.3.1 Theme: Release <#Theme:_Release>`_

      -  `1.1.4 Thursday 12/6 <#Thursday_12.2F6>`_

         -  `1.1.4.1 Theme:
            Interoperability <#Theme:_Interoperability>`_

   -  `1.2 Discussion and work topics <#Discussion_and_work_topics>`_
   -  `1.3 Speakers and tentative
      topics <#Speakers_and_tentative_topics>`_
   -  `1.4 Discussion Board <#Discussion_Board>`_

Dirac meeting Odense June 9-13 2008
===================================

Schedule
--------

Please add speakers to the schedule.

It is better that we already define possible tasks to be performed
during the meeting before we get to Odense -andre

Monday 9/6
~~~~~~~~~~

Theme: Tools for distributed development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  9:00 -- 10:45 The Fathers + Node leaders

   -  HJJ's initial words
   -  The Fathers speak
   -  Node leaders on the respective developments

-  10:45 -- 11:00 Break
-  11:00 -- 12:30 A walkthrough on the latest developments

   -  11:00: `Subversion </dirac/index.php/SVN>`_,
      `Wiki </dirac/index.php/Wiki_guidelines>`_ and
      `other <http://radon.chem.vu.nl/cgi-bin/gnatsweb.pl>`_ project
      management-related topics
      (`test <http://www.pgbuildfarm.org/cgi-bin/show_status.pl>`_/`compile
      farm <http://en.wikipedia.org/wiki/Compile_farm>`_ etc). Q/A and
      feedback. (UE & HJJ)
   -  11:45: Overhauling Dirac's build system with
      `Foray <http://code.google.com/p/foraytool>`_ (documentation
      written under its former name
      `tcbuild <http://macresearch.org/tcbuild-new-build-tool-fortran>`_)
      (ASPG)
      `Presentation <http://wiki.chem.vu.nl/dirac/images/f/f0/DiracMeetingOdense2008_Foray.pdf>`_

-  12:30 -- 13:30 Lunch
-  13:30 -- 15:45 Time to work together on programming or other tasks.
-  15:45 -- 16:00 Break
-  16:00 -- 17:00 "Content" talks

   -  16:00 -- 16:15 Overview of F90 functionality (ASPG)
      `Presentation <http://wiki.chem.vu.nl/dirac/images/0/0d/DiracMeetingOdense2008_F90inDirac.pdf>`_
   -  16:15 -- 16:30 Getting rid of MEMGET(): dynamic memory allocation
      in RELCCSD (ASPG)
      `Presentation <http://wiki.chem.vu.nl/dirac/images/b/bc/DiracMeetingOdense2008_allocation.pdf>`_
      `Draft
      documentation <http://wiki.chem.vu.nl/dirac/images/0/07/Allocator_reference.pdf>`_
   -  16:30 -- 17:00 My recent technical developments of F90 in DIRAC,
      second quantization operators strings. (MI)

Tuesday 10/6
~~~~~~~~~~~~

Theme: Modularity
^^^^^^^^^^^^^^^^^

-  9:00 -- 9:30 Advantages of increasing modularity by object oriented
   programming`[1] <http://www.macanics.net/downloads/ProgrammingParadigms.zip>`_.
   (LV)
-  9.30 -- 10.15 Dirac system requirements, writing clean 64 bit code
   `[2] <http://www.unix.org/whitepapers/64bit.html>`_
   `[3] <http://developer.apple.com/documentation/Darwin/Conceptual/64bitPorting/intro/chapter_1_section_1.html>`_
   `[4] <http://www.acmqueue.org/modules.php?name=Content&pa=showpage&pid=421&page=4>`_
   (and testing it). (UE & others)
-  10:15 -- 10:45 Break

-  11:00 -- 11:45 Current status of the code: unnecessary entanglements
   (UE)
-  11:45 -- 12:15 Possibilities for using object-oriented
   design/implementation in Dirac
   `[5] <http://en.wikipedia.org/wiki/Object_oriented>`_
   `[6] <ftp://ftp.liv.ac.uk/pub/F90Course/CourseNotes.ps.gz>`_ (JHS)
-  12:15 -- 12.30 Interoperability with other codes or scripting
   languages (LV)
-  15:45 -- 16:00 Break
-  16:00 -- 17:00 "Content" talks

   -  16:00 -- 16:30 Recent features developments in Dirac : Complex
      Absorption Potentials, proposal for interfacing Hilbert Space CC
      (MI)
   -  16:30 -- 17:00 Open-ended response behind the scenes: bringing
      DALTON and DIRAC closer together with F90 operator overloading
      (RB)

Wednesday 11/6
~~~~~~~~~~~~~~

Theme: Release
^^^^^^^^^^^^^^

-  9:00 -- ?? Roundtable discussion about the release of Dirac08.
   Everyone brings some small things that can be improved in the
   afternoon (bug fixes etc).
-  10:30 -- 10:45
-  10:45 -- 11:00 Break
-  12:30 -- 13:30 Lunch
-  15:45 -- 16:00 Short commercial for Jonas Juselius' libgetkw to show
   how input parsing can look like (RB)
-  16:00 -- 17:00 "Content" talks

   -  ??:?? Soccer match

Thursday 12/6
~~~~~~~~~~~~~

Theme: Interoperability
^^^^^^^^^^^^^^^^^^^^^^^

-  9:00 -- 9:30 Data communication between Dirac and other codes using
   `XML <http://www.w3schools.com/xml/default.asp>`_(JHS)
-  9:30 -- 10:30 Defining a new interface file format (DFCOEF+more)
-  10:30 -- 11:00 Break
-  12:30 -- 13:30 Lunch
-  15:45 -- 16:00 Break

Discussion and work topics
--------------------------

Manual and user friendliness

-  DiracWiki and web in general
-  Input interface (GUI but also input parsing)
-  Routines for release and documentation
-  Machine readable output (XML?)
-  Should we have here (on Diracwiki) the DIRAC programmers manual
   accessible only for dirac developers ?

   -  There is a very good integrated wiki/bug tracker/svn interface
      called Trac that is better used for this
      It is very difficult to restrict access to some pages of this
      wiki. --Ulf

Development (technical)

-  F90: modules, generating module dependencies for proper subsequent
   compilation (sfmakedepend.pl)
-  F90: how to "nicely" access common block variables (like ntbas,
   norbt, ...)
-  F90 derived data types for encapsulating quaternion matrices
   (operators)
-  pamadm and the build system (configure/makefiles etc) in general.
-  SVN server
-  Build/test farm (and selection of the "standard" F90 - gfortran?),
   automatic testing, validation, unit tests, ...
-  Memory allocation/deallocation/management/statistics
-  Modularization (Should Dirac remain one single monolithic program?)
-  File formats (extensions and documentation).
-  Standard input parser (now we have a lot of redundancy here)
-  Requirements for Interoperability (a fancy word to indicate that
   Dirac is used together with other codes, e.g. ADF)
-  What should be the system requirements for running Dirac? Can we
   assume users have python, Fortran 2003 etc?
-  Removal of redundant functions, general clean-up.
-  When will we have proper arrays dimensions in abacus (it means
   possibility to run DIRAC compiled with an 'arrays dimension check'
   flag)
-  Since we are using only ONE F90 compiler for the whole DIRAC, do we
   need both "F90\_OPTFLAGS" and "OPTFLAGS" in the "Makefile.config" ?

Development (features)

-  Fast DFT (we need a f90 XC integrator for several projects and it is
   a good excuse to re-discuss the fast DFT)
-  Future direction for DIRAC (functionality, efficiency, manpower)
-  Tensor contraction engine (Luuk?; I'm very interested --Miro)
-  Structure optimization and frequencies
-  What to do ourselves and where to use interfaces to other codes ?

Speakers and tentative topics
-----------------------------

-  Ulf Ekstrom - Modern collaborative program development: SVN/Wiki/Trac
   and all that.
-  Radovan Bast - F90 derived type for quaternion matrices and operator
   overloading - some ideas to start a discussion.
-  Radovan Bast - Short commercial for Jonas Juselius' libgetkw, "a
   simple input file parser with a C-ish syntax, with language bindings
   for C, Fortran 90 and F77" - even if we don't use it we can learn
   from it - there are some nice ideas for input parsing.
-  Miroslav Ilias - Recent features developments in Dirac : Complex
   Absorption Potentials, proposal for Hilbert Space CC.
-  Miroslav Ilias - My recent technical developments of F90 in DIRAC,
   second quantization operators strings.
-  Andre Gomes - General overview of the F90 work so far; memory
   allocation and its use in RELCCSD.
-  Jetze Sikkema - I could speak about XML and / or about modules.
-  Luuk Visscher - The need and possibilities for modularization and
   clean interfaces

Discussion Board
----------------

-  Jetze, Andre, Miro, ..., if you like you can talk about the
   quaternion derived data type instead of me. it's not that i know much
   about it. it's only that i need it and that it would be nice to agree
   on a standard definition --radovan
-  Radovan, I was thinking of presenting an overview of what's going on,
   and then talking about what I've done in more detail. I think then
   the people working on a particular section should be given the floor.
   For that, we should coordinate a bit before the meeting -- Andre
-  I agree with Andre, let each of us presents what he has done so far.
   For example, F90 modules for vectors and matrices are from me and
   from Radovan as well; we should discuss a common direction --Miro
-  OK, my aim is to show how we do matrix operations "the lazy way" with
   overloaded operators and how we can hide matrix flags like point
   group symmetry, packing time-reversal symmetry in the matrix type. I
   agree that there should be a general overview first --radovan
-  On the meeting set up: what about talks only in the morning, say
   9-11.30, and leave the afternoon to actually work together ? -- Luuk
-  Maybe we could try to move important decisions (about pamadm, foray,
   DFCOEF) to the first day so that we have time to work on the topics
   after the decisions are made? --radovan

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Meeting2008 <http://wiki.chem.vu.nl/dirac/index.php/Meeting2008>`_"

Views
'''''

-  `Article </dirac/index.php/Meeting2008>`_
-  `Discussion </dirac/index.php/Talk:Meeting2008>`_
-  `View source </dirac/index.php?title=Meeting2008&action=edit>`_
-  `History </dirac/index.php?title=Meeting2008&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Meeting2008>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Meeting2008>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Meeting2008>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Meeting2008&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Meeting2008&oldid=10641>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 07:20, 10 June 2008.
-  This page has been accessed 2,471 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
