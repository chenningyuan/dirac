Localizing bugs in DIRAC run using debugger
===========================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Since the DIRAC program system is still evolving, it gets spoiled with
bugs during its perpetual development. We would like to share knowledge
on how to localize crashes/bugs in the code.

We believe that this short tutorial is useful not only for DIRAC
developers, but also for ordinary users, who likewise could send us
precise localization and specification of the bug.

Contents
--------

-  `1 Introduction <#Introduction>`_
-  `2 Debugging through the "pam"
   script <#Debugging_through_the_.22pam.22_script>`_

   -  `2.1 Stack printing <#Stack_printing>`_

-  `3 Debugging in the plain mode with the GNU/gdb
   debugger <#Debugging_in_the_plain_mode_with_the_GNU.2Fgdb_debugger>`_

   -  `3.1 Demonstration of GNU/gdb debugging
      commands <#Demonstration_of_GNU.2Fgdb_debugging_commands>`_

-  `4 Intel Fortran
   compiler/debugger <#Intel_Fortran_compiler.2Fdebugger>`_
-  `5 Debugging DIRAC in the parallel
   mode <#Debugging_DIRAC_in_the_parallel_mode>`_

Introduction
------------

We assume a serial compiled DIRAC program suite.

Make sure that DIRAC is compiled with "-g" flag for both Fortran and C
compilers. This flag is the be the default part of compiler flags.

When you would like to run multiple debugging runs, please choose *safe
flags* for the compilation (it means omitting optimization flags, see
the "Makefile.config" file) in order to speed up repeating DIRAC
compilation and improve ease of debugging. Note that some bugs (crashes)
might "disappear" when the optimization is turned off, because not all
compilers have bug-free optimization procedure.

Ensure that proper debugger is placed into the "pam" script (and is
likewise installed). For the GNU set (gfortran,gcc) this is GNU/gdb.
Some commercial compilers are providing their own debuggers: ifort,icc -
idbg, for pgf90/pgcc - pgdbg etc.

Debugging through the "pam" script
----------------------------------

Running DIRAC with "pam -debug ..." is suitable for one-shoot debugging.
Pam script calls the debugger, which afterwards shows up its command
line.

Then the user can set up breakpoints, stops for the code and run
'dirac.x' inside the debugger.

Stack printing
~~~~~~~~~~~~~~

There is (in the development version of DIRAC) the possibility to print
out the stack content (programmed by Andre).

Any 'QUIT' crash of the code causes a SegFault termination you can
investigate with the debugger - you can 'visit' individual procedures
within the stack sequence and examine their variables.

For example, beeing the "gdb" debugger session:

::

    > run
    (Dirac runs and crashes somewhere)
    > backtrace
    (show the stack)
    > up {n}
    (go to the n-the previous function in the call stack; see help up)
    > print <VARIABLE>
    (whatever you want to inspect).

Oposite way can be achiewed using 'down' command.

Debugging in the plain mode with the GNU/gdb debugger
-----------------------------------------------------

We give a simple demo-session of an thorough and effective debugging.
Recent version of the widely used the GNU/gdb debugger is recommended.

Make a working directory with enough disk space. Provide a link of the
$DIRAC/dirac.x executable file there. This is because you will repeat
modifications of DIRAC source files and subsequent compilation, so it is
easier to have linked executable rather than copying it.

Make sure that the 'DIRAC.INP', 'MOLECULE.INP', *gdb\_demo* and of other
necessary files (e.g. 'DFCOEF' if you are restarting SCF iterations) are
in the working directory as well.

Demonstration of GNU/gdb debugging commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This *gdb\_demo* file (create&copy from this web-page) contains some
commands (and macros) of the gdb debugger:

::

    ###  load the executable file with all symbols
    file dirac.x
     
    ### This settings is necessary for complex breakpoint conditions in fortran...
    #set language fortran
    set language c
    set language auto
     
    ### Write where I am...
    echo In directory:
    shell echo pwd
    echo \n Files in dir:
    shell ls -a
     
    ##################################################################
    ### Delete  DIRAC files in the scratch directory
    ### - needed for restarting and for further catching of bugs..
    ##################################################################
    shell /bin/rm -r DF*  AOPROPER BSSMAT
    echo \n After deletion remaining files in dir:
    shell ls -a
     
    ### clean the desk: delete all previous breakpoints
    d b
     
    ###################################
    # Example of macros
    ###################################
    define my_macro
    printf "in GMOTRA begin NZ=%d N2BBASXQ=%d SPINFR=%d\n",nz,n2bbasxq,spinfr
    end
     
    define my_macro1
    printf "in GMOTRA begin SPINFR2=%d\n",spinfr2
    end
     
    ### demo of a conditional breakpoint
    b gmotra_ if (nz==1 | (n2bbasxq==0 && spinfr==0))
    command
     my_macro
     my_macro1
     echo type continue or sipmly c...\n
    # c
    end
     
    ### example of the unconditional breakpoint
    b dirone.F:1399
    command
     if ( (nfsym.eq.1 && n2bbasxq.gt.10) | nz.eq.1 | spinfr.eq.0)
      printf "NFSYM=%d\n", NFSYM
      if (nfsym.eq.2)
        printf "NTMO(1)=%d NTMO(2)=%d\n",NTMO(1),NTMO(2)
      else
        printf "One symmetry...NTMO(1)=%d \n",NTMO(1)
      end
     else
      printf "... else branch...."
     end
      
     ##### while cycle demo #####
     printf" Few values of the WORK(KTMAT) array, KTMAT=%d... \n",KTMAT
     set $indx=0
     while ($indx.lt.5)
      printf "WORK(%d)=%lf\n",ktmat+$indx,(*WORK)(ktmat+$indx) 
      #  p ktmat+$indx
      set $indx=$indx+1
     end
    end

Being the working directory type 'gdb dirac.x' to run the 'dirac.x' in
the debugging mode. When the gdb command line appears, insert 'source
gdb\_demo' to load the debugger commands source file.

Each time when you modify and recompile DIRAC please retype again 'file
dirac.x' in the gdb mode to reload the fresh executable.

If you modify the *gdb\_demo* debugger source file apply again the gdb
command 'source gdb\_demo' to reload it.

More GNU/gdb know-how is given at the web-page
`http://sources.redhat.com/gdb/current/onlinedocs/gdb.html <http://sources.redhat.com/gdb/current/onlinedocs/gdb.html>`_.

Intel Fortran compiler/debugger
-------------------------------

For debugging programs compiled with ifort you need to use the Intel
Debugger idb, which is the same set of commands as the GNU/gdb.

Debugging DIRAC in the parallel mode
------------------------------------

For parallel debugging we recommend some commercial software, which can
follow individual threads.

...totalview debugger ?

