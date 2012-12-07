

How this documentation works and how you can contribute
=======================================================

This documentation is generated using `sphinx <http://sphinx.pocoo.org/>`_
which translates RST (reStructured Text) markup to html and latex.  Sphinx and
RST are extremely widely used in the python community.  The nice thing about
RST markup that it looks nice by itself in the terminal.  RST and Sphinx are
very simple and intuitive and for first steps consult the `reStructuredText
Primer <http://sphinx.pocoo.org/rest.html>`_ and the `Sphinx Markup Constructs
<http://sphinx.pocoo.org/markup/index.html>`_.  Typically RST markups is
auto-colored in the terminal so you can see when the markup is correct or
wrong.

Generation of html pages
------------------------

The central documentation web pages are generated automatically
from the RST files every hour.
But you can also generate the web pages yourself. For this install
the packages "python-sphinx" and "python-matplotlib", then::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make html

Now point your browser to build/_html/index.html


How to add new sections and text
--------------------------------

Modify existing RST files or add new ones and reference them
in index.rst.
