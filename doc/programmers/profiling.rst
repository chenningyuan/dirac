
Profiling how to
================================================

This is a quick and dirty guide to profiling. You are **strongly** recommended to read
the profiling software documentation to get deeper insight in the commands suggested 
in the following.

Profiling with Intel VTune Amplifier XE  
------------------------------------------------

This software is a nice tool for profiling and has a GUI working also on Linux.
Moreover, you will not need to compile your program with profiling flags on.
The only bad thing is that it is commercial software.
Profiling takes two steps:
1. run your program under VTune;
2. analyze the collected data.

Assuming your program is called ``foobar``, and further assuming that you only want to 
collect the *hotspots*, for step 1 you will launch::

	amplxe-cl -collect hotspot -r /results/directory ./foobar

the results of the sampling will be put in ``/results/directory``. If the
results directory is not specified VTune will put everything in a subdirectory
of the current directory.
For step 2::
	
	amplxe-cl -report report_name -r /results/directory

this will produce a hotspots report called ``report_name``.

Profiling DIRAC with Intel VTune Amplifier XE
------------------------------------------------

The situation is a bit different if you want to profile DIRAC, because we usually
launch the executable through a script. 
If you are using the ``wrapper.py`` script, you will specify as launcher::

	--launcher="amplxe-cl -collect hotspot -r /results/directory"

the rest of the procedure remains the same.
What if we want to profile during an MPI run? We modify the launcher as follows::

	--launcher="mpirun -np 12 amplxe-cl -collect hotspots -follow-child -mrte-mode=auto -target-duration-type=medium -no-allow-multiple-runs -no-analyze-system -data-limit=100 -slow-frames-threshold=40 -fast-frames-threshold=100 -r /results/directory"

Some more words of comment (shamelessly copied from ``amplxe-cl -help collect``):

- ``follow-child``, collects data on processes launched by the target process;

- ``mrte-mode``, selects the profiling mode;

- ``target-duration-type``, estimates the application duration time. This value affects the size of collected data;

- ``no-allow-multiple-runs``, disables multiple runs to achieve more precise results for hardware event-based collections;

- ``no-analyze-system``, disables analyzing all processes running on the system;

- ``data-limit``, limits the amount of raw data to be collected;
 
- ``slow-frames-threshold``, specifies a threshold to separate slow and good frames. It must be smaller than the threshold for fast frames; 

- ``fast-frames-threshold``, specifies the threshold to separate good and fast frames.



	
	
