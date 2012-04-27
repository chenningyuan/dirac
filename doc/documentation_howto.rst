

How this documentation works and how you can contribute
=======================================================

This documentation is generated using sphinx http://sphinx.pocoo.org/ which
translates RST (ReStructured Text) markup to html and latex.  Sphinx and RST
are extremely widely used in the python community.  The nice thing about RST
markup that it looks nice by itself in the terminal.  RST is very simple and
intuitive and for first steps use these guides:

 * http://docutils.sourceforge.net/docs/user/rst/quickstart.html
 * http://docutils.sourceforge.net/docs/user/rst/quickref.html
 * http://docutils.sourceforge.net/docs/user/rst/cheatsheet.txt

Typically RST markups is auto-colored in the terminal so you can
see when the markup is correct or wrong.

Installation of the documentation tool
--------------------------------------

In order to create documentation you need to install "python-sphinx"
and "python-matplotlib" with related packages on your Linux machine.


Generation of html pages
------------------------

The central documentation web pages are generated automatically
from the RST files every hour.

But you can also generate the web pages yourself. For this install
the python-sphinx package, then::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make html

Now point your browser to build/_build/html/index.html


How to add new sections and text
--------------------------------

Modify existing RST files or add new ones and reference them
in index.rst.
