How to safely move your uncommited files to SVN
===============================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

If you decide/need to move your modified files to SVN **without**
committing them to CVS first, this is a safe way:

-  Update your CVS trunk working copy or branch working copy
-  Check out a fresh SVN working copy (either trunk or a specific
   branch)
-  Move modified files by hand to the new SVN working copy
-  **Important**: Check with ``svn status`` whether SVN marks the same
   files as modified (the list of files given by ``svn status`` should
   look like the list of files given by ``cvs update``). If this is not
   the case and/or you have doubts, don't delete the CVS working copy
   too early.

