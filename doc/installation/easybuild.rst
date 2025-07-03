:orphan:
 

Installation using EasyBuild
============================

We recommend using EasyBuild to build DIRAC on HCP systems.

If you or your cluster administrator are unfamiliar with it, we recommend following the `installation <https://docs.easybuild.io/installation/>`_ and `usage <https://docs.easybuild.io/using-easybuild/>`_ instructions provided on the `EasyBuild website <https://docs.easybuild.io/>`_. Additionally, there is an instructive `tutorial <https://easybuild.io/tutorial/>`_ that can be used to configure it.

Once EasyBuild is installed on your system, as a user, you can install DIRAC employing, for example, the following commands::

  $ module load EasyBuild/5.1.0
  $ module use /path-to-your-modules-directory/software/modules/all
  $ eb -r --buildpath=/your-build-directory/ --installpath=/path-to-your-modules-directory/software --repositorypath=/path-to-your-modules-directory/software/repository --sourcepath=/path-to-your-modules-directory/software/sources DIRAC-23.0-intel-2023a.eb

Once installed, you can use DIRAC simply by typing in your runfile, for example: ::

  $ module load DIRAC/23.0-intel-2023a
  
and to run DIRAC you must use: ::

  $ pam-dirac --inp=... --mol=... [...]

Current versions of DIRAC supported by EasyBuild can be found `here <https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/d/DIRAC/>`_ and `here <https://docs.easybuild.io/version-specific/supported-software/d/DIRAC/>`_.

However, if you still want to compile the program manually, you can follow the instructions provided in the :doc:`Installation instructions for system administrators <sysadmins>` section.
