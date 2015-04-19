.. _documentation_howto:

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


Sphinx installation
-------------------
The sphinx generator is part of standard packages, and on most systems
can be installed via their official repositories.
Common names of packages are "python-sphinx" and "python-matplotlib", but there might be
additional packages needed depending on the operating system.

Linux users, for example on Ubuntu, can get the documentation generator as::

 sudo apt-get install python-sphinx python-matplotlib

For  *Mac* users using the macport ::

 sudo port install python-sphinx
 sudo port select --set sphinx py27-sphinx
 sudo port install py27-matplotlib


Generation of html pages
------------------------

The central documentation web pages are generated automatically
from the RST files every hour.
But you can also generate the web pages by yourself. For this install
the sphinx-packages first (see above), then::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make -k html > build_html.log 2> build_html.log

Now point your browser to build/html/index.html. We recommend the developer to check the *build_html.log* file for possible errors
in generating the documentation.


How to add new sections and text
--------------------------------

Modify existing RST files or add new ones and reference them
in index.rst.
