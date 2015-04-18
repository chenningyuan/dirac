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

For  *Mac* users using the macport::

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

New theme
---------

The `sphinx_rtd_theme <https://github.com/snide/sphinx_rtd_theme>`_ has been customized and integrated as the new theme for Dirac. Theme extension can be installed using pip:

.. code-block:: bash

    pip install sphinx_rtd_theme

Here is example configuration with the fallback to the old theme:

.. code-block:: python

    try:
        import sphinx_rtd_theme
        html_theme = 'sphinx_rtd_theme'

        # set paths for loading local and also theme static resources
        html_theme_path = ['_themes', sphinx_rtd_theme.get_html_theme_path()]
        html_static_path = ['_static']

        # styles can be customized in this file
        html_context = {
            'css_files': ['_static/theme-custom.css'],
        }
    except ImportError:
        # fallback to old theme
        html_theme = 'dirac'
        html_theme_path = ['_themes']

If you want to customize the HTML templates, just copy the original ones from the theme and place them
inside the ``_templates`` folder inside the project. Sphinx will automatically prefer the ones in the ``_templates`` folder
over the global ones from the theme itself.

.. code-block:: bash

    # layouts can be customized in these paths
    _templates/layout.html
    _templates/footer.html

To override/customize the CSS styles, there is a css file placed under the ``_static`` folder called ``theme-custom.css``.

New format of citations
-----------------------

Citations have been updated to use the BibTeX format. This format is more maintainable and standardized. Support for this format is provided by the `sphinxcontrib-bibtex <http://sphinxcontrib-bibtex.readthedocs.org>`_ extension and it can be installed using pip:

.. code-block:: bash

    pip install sphinxcontrib-bibtex

Citations can be added using the new :rst:dir:`cite` directive:

.. code-block:: rest

    # old format
    [Dubillard2006]_

    # new format
    :cite:`Dubillard2006`

The new format also enables the option to download the citations in BibTeX. This code will generate a link to static BibTeX file with all citations:

.. code-block:: rest

    :download:`BibTeX <zreferences.bib>`

The list of references can be added to the document using the :rst:dir:`bibliography` directive:

.. code-block:: rest

    .. bibliography:: references.bib
      :enumtype: upperroman
      :style: diracstyle

.. note::

    Ensure that the bibliography directive is processed after all cites. The simplest way is to
    name the references page so that it starts with the "z" letter in order to be alphabeticaly last.
    The issue is being tracked `here <http://sphinxcontrib-bibtex.readthedocs.org/en/latest/usage.html#unresolved-citations-across-documents>`_.

Download mathjax formulas in LaTeX format:
------------------------------------------

The mathjax extension has been extended and renamed to **mathjaxplus**. The new **mathjaxplus** extension provides a way to download formulas in LaTeX format. There are also configuration parameters exposed to customize the generation process:

.. code-block:: python

    # disable or enable
    mathjax_generate_latex = True

    # style of save as link
    mathjax_latex_style = 'font-size:14px;text-align:right;'

    # title of the save as link
    mathjax_latex_title = 'Save as Latex'

Gitlog - display the GIT commit history for pages
-------------------------------------------------

A new extension called **gitlog** has been developed and integrated to the documentation. The motivation for this extension was to display information about the last git commit for each document in the documentation. There are some parameters exposed to configure the process:

.. code-block:: python

    # path to the git repository
    gitlog_source_path = '../doc/'

    # date format used to format git commit date in templates
    gitlog_date_format = '%d %b %Y, %H:%M:%S'

The extension adds following data available to the context of given document:

.. code-block:: bash

    git_author_name  # name of the commit author
    git_author_email # email address of the commit author
    git_commit_sha   # sha1 hash of the commit
    git_commit_date  # date of the commit

This code shows a simple usage in Jinja2 template:

.. code-block:: jinja

    {% if git_author_name %}
        Committer name: {{ git_author_name }}
        Committer email: {{ git_author_email }}
        Commit date: {{ git_commit_date }}
        Commit sha1 hash: {{ git_commit_sha }}
    {%- endif %}


Google Analytics
----------------

Google Analytics support for both slides and main documentation has been added. In order to support this
feature a new ``_templates/shared/macros.html`` file has been added to enable Jinja 2 macros sharing between
two documentation systems.

Setting Google Account for main documentation (in file ``conf.py``):

.. code-block:: python

    html_context = {
        'google_analytics_account': 'UA-XXXXX',
    }

Setting Google Account for slides (in file ``workshop-slides/conf.py``):

.. code-block:: python

    slide_theme_options = {
        'google_analytics_account': 'UA-XXXXX',
    }

Example usage of the ``google_analytics`` macro in Jinja 2 template:

.. code-block:: jinja

    {% import "shared/macros.html" as macros %}
    {{ macros.google_analytics('UA-XXXXX') }}
