Wiki guidelines
===============

From Diracwiki
~~~~~~~~~~~~~~

(Redirected from
`Guidelines </dirac/index.php?title=Guidelines&redirect=no>`_)

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

These guidelines are meant to help, not to restrict. They are not
definitive and all contributions/changes/improvements/ideas are most
welcome.

Contents
--------

-  `1 General guidelines <#General_guidelines>`_
-  `2 Mediawiki syntax <#Mediawiki_syntax>`_

   -  `2.1 Formulas <#Formulas>`_
   -  `2.2 File uploads <#File_uploads>`_

-  `3 Hints and suggestions <#Hints_and_suggestions>`_
-  `4 References <#References>`_
-  `5 Style templates <#Style_templates>`_
-  `6 DIRAC amboxes <#DIRAC_amboxes>`_
-  `7 Internal and external links <#Internal_and_external_links>`_
-  `8 Units <#Units>`_

General guidelines
------------------

The user manual is for DIRAC developers but also for general users. It
should be as simple, fast, and convenient as possible to find the
necessary information. Please mark difficult to understand content with
{{clarifyme}} (looks like\ :sup:`[*please
`clarify </dirac/index.php/Clarify>`_*]`\ ) or {{example}} (looks
like\ :sup:`[*please provide an
`example </dirac/index.php/ExampleNeeded>`_*]`\ ). Such content can
later be found automatically through `this
page </dirac/index.php/Special:Categories>`_.

Try to avoid internal variable names and too formal Fortran format
specifications since they may confuse the general user. Try to use
explicit examples which can be copy-pasted into the input file and
adapted. Be careful with the blank characters in front of the input
examples.

*Slightly more useful for the developer:*

::

    == {{keyword|.PRINT}} ==
    Set the print level.

    ''Arguments:'' Integer IPRLOC
    ''Default:'' IPRLOC = 1

*More useful for the user (favored):*

::

    == {{keyword|.PRINT}} ==
    Set the print level.

    ''Default:''
     .PRINT
      1

We should try to arrange keywords in a logical order. The user should be
able to find related keywords either close to each other, or linked to
each other.

Mediawiki syntax
----------------

A quick way to learn the syntax is to use "perturbation theory". Browse
through this wiki or other Mediawiki-based wikis, find a part that looks
almost like what you have in mind, click on edit (top of the page or on
the right), then copy-paste the source to the new wiki page, adapt, and
improve it. See also the `Mediawiki user
help <http://www.mediawiki.org/wiki/Help:Contents>`_, the `Mediawiki
magic words <http://en.wikipedia.org/wiki/Help:Magic_words>`_, and the
`Wikipedia how-to
edit <http://en.wikipedia.org/wiki/Wikipedia:How_to_edit_a_page>`_.

Formulas
~~~~~~~~

Formulas can be written using a large subset of all TeX and amsmath
commands. Tell the wiki that you are writing formulas by enclosing the
TeX in math tags, like

::

    <math>f(x) = \int_0^x dx</math> 

which renders as |f(x) = \\int\_0^x dx|. Note that you should not
enclose your formulas in $$ or \\begin{equation}.

File uploads
~~~~~~~~~~~~

It is possible to upload files to a global file area of the wiki, and
then link these files from particular pages. Files are uploaded through
the link in the left hand side toolbox, and can be linked either as
images

::

    [[Image:yourimage.png]]

or, for non-image file types, as

::

    [[Media:yourfile.zip]]

**NB: All file names have to be unique. This means that you should name
your files after the context where they appear. Names such as
example.inp are strictly *verboten*.**

Hints and suggestions
---------------------

-  If you use Firefox 2 and prefer to use an external text editor, check
   out a great tool called `It's All
   Text <https://addons.mozilla.org/en-US/firefox/addon/4125>`_ (see
   also
   `[1] <http://en.wikipedia.org/wiki/Wikipedia:Text_editor_support#How_to_open_articles_in_your_text_editor>`_).

-  Server runs on UTC time, update your wiki preferences to set your own
   time zone.

References
----------

Use the style templates {{article}}, {{book\_chapter}}, {{book}}, and
{{unpublished}} for references.

Example for an article:

::

    <ref name="Saue:2002:jcomputchem">
    {{article
    | authors = T. Saue and T. Helgaker
    | title   = Four-component relativistic Kohn-Sham theory
    | journal = J. Comput. Chem.
    | volume  = 23
    | page    = 814
    | year    = 2002
    | doi     = 10.1002/jcc.10066
    }}
    </ref>

The tag name="Saue:2002:jcomputchem" makes it possible to use the same
reference several times within one wiki page. For a second (third ...)
citation within the same page, it would be sufficient to type <ref
name="Saue:2002:jcomputchem" /> (mind the closing XML slash character).
It is important to provide the DOI to have a working link to the
electronic version of the reference. If the DOI is not known (e.g. no
electronic version exists), leave the entire doi line out.

Example for a book chapter:

::

    <ref name="visscher:post-dhf">
    {{book_chapter
    | authors   = L. Visscher
    | title     = Post Dirac-Hartree-Fock methods - Electron Correlation
    | editors   = P. Schwerdtfeger
    | booktitle = Relativistic Electronic Structure Theory - Part 1: Fundamentals
    | publisher = Elsevier
    | address   = Amsterdam
    | year      = 2002
    }}
    </ref>

Style templates
---------------

Style templates make it easier to have a common layout and to easily
change the functionality for a larger number of objects inside the text.

Please use these style templates, improve them, and add new ones (see
e.g. `this
guide <http://strategywiki.org/wiki/StrategyWiki:Guide/Templates>`_).

+------------------------------------------------------------------------------------+-----------------------------+
| **how you type it**                                                                | **how it looks**            |
+------------------------------------------------------------------------------------+-----------------------------+
| {{dirac}}                                                                          | DIRAC                       |
+------------------------------------------------------------------------------------+-----------------------------+
| {{dirac04}}                                                                        | DIRAC04                     |
+------------------------------------------------------------------------------------+-----------------------------+
| {{dirac08}}                                                                        | DIRAC08                     |
+------------------------------------------------------------------------------------+-----------------------------+
| {{dalton}}                                                                         | DALTON                      |
+------------------------------------------------------------------------------------+-----------------------------+
| {{grasp}}                                                                          | GRASP                       |
+------------------------------------------------------------------------------------+-----------------------------+
| {{molfdir}}                                                                        | MOLFDIR                     |
+------------------------------------------------------------------------------------+-----------------------------+
| {{relccsd}}                                                                        | RELCCSD                     |
+------------------------------------------------------------------------------------+-----------------------------+
| {{lucita}}                                                                         | LUCITA                      |
+------------------------------------------------------------------------------------+-----------------------------+
| {{keyword\|.KEYWORD}}                                                              | .KEYWORD                    |
+------------------------------------------------------------------------------------+-----------------------------+
| {{unreleased\_keyword\|.KEYWORD}}                                                  | .KEYWORD                    |
+------------------------------------------------------------------------------------+-----------------------------+
| {{highlight\|this}}                                                                | this                        |
+------------------------------------------------------------------------------------+-----------------------------+
| {{pink\|INPUT}}{{green\|INPUT}}{{blue\|INPUT}}{{yellow\|INPUT}}{{orange\|INPUT}}   | INPUTINPUTINPUTINPUTINPUT   |
+------------------------------------------------------------------------------------+-----------------------------+
| {{fix\_me\|obsolete}}                                                              | obsolete                    |
+------------------------------------------------------------------------------------+-----------------------------+
| {{news\|2007-02-07\|Good news.}}                                                   | 2007-02-07 Good news.       |
+------------------------------------------------------------------------------------+-----------------------------+

DIRAC amboxes
-------------

This section lists various DIRAC ambox (article message box) templates.
Feel free to change them and to add new ones. For changing them edit the
page and scroll to the bottom of the page.

**how you type it**

**how it looks**

{{warning\_unreleased\_material}}

`|image1| </dirac/index.php/Image:Dirac.png>`_

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

{{advanced\_options}}

`|image2| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

{{programmers\_options}}

`|image3| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

Internal and external links
---------------------------

::

    Link to an html page like this [http://www.povray.org/ POV-Ray] or like that [http://www.povray.org/].

*Gives:* Link to a html page like this
`POV-Ray <http://www.povray.org/>`_ or like that
`[2] <http://www.povray.org/>`_.

::

    Link to the keyword .LVNEW on the same wiki page like this [[#.LVNEW|.LVNEW]].

*Gives:* Link to the keyword .LVNEW on the same wiki page like this
`.LVNEW <#.LVNEW>`_.

::

    Link to the keyword .EVCCNV on a different wiki page like this [[Manual:WAVE_FUNCTION:DHFCAL#.EVCCNV|.EVCCNV]].

*Gives:* Link to the keyword .EVCCNV on a different wiki page like this
`.EVCCNV </dirac/index.php/Manual:WAVE_FUNCTION:DHFCAL#.EVCCNV>`_.

::

    Link to a [[Manual:One-electron operators|different wiki page]] like this.

*Gives:* Link to a `different wiki
page </dirac/index.php/Manual:One-electron_operators>`_ like this.

Units
-----

It is not correct to write out units capitalized or pluralised, the
correct way is to write: Å, *a*\ :sub:`0`\ , Hz, ångström, bohr, hertz,
... (see for instance `[3] <http://en.wikipedia.org/wiki/SI>`_).

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Wiki\_guidelines <http://wiki.chem.vu.nl/dirac/index.php/Wiki_guidelines>`_"

`Categories </dirac/index.php/Special:Categories>`_: `Articles needing
clarification </dirac/index.php/Category:Articles_needing_clarification>`_
\| `Articles needing
examples </dirac/index.php?title=Category:Articles_needing_examples&action=edit>`_

Views
'''''

-  `Article </dirac/index.php/Wiki_guidelines>`_
-  `Discussion </dirac/index.php?title=Talk:Wiki_guidelines&action=edit>`_
-  `View source </dirac/index.php?title=Wiki_guidelines&action=edit>`_
-  `History </dirac/index.php?title=Wiki_guidelines&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Wiki_guidelines>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Wiki_guidelines>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Wiki_guidelines>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Wiki_guidelines&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Wiki_guidelines&oldid=9311>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 13:15, 30 May 2008.
-  This page has been accessed 4,150 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |f(x) = \\int\_0^x
dx| image:: /dirac/images/math/2/b/d/2bd9372cd1b075fbcab938c6ffc8c5f0.png
.. |image1| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image2| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image3| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
