=====================================================================
Davidson size-consistency corrections for truncated CI wave functions
=====================================================================

----------
Background
----------

A serious deficiency of truncated CI expansions, including the most popular ones, namely (multi-reference; MR)CI-SD is
the lack of size-consistency; i.e., the energy of the system does not scale properly with the system size. Besides
corrections to the CI equations which will improve the energy and the wave function, *a posteriori* corrections (most often
referred to as *Davidson-type* corrections or *+Q*) to the energy are very popular. A comprehensive review on this topic has been
given in a recent `review paper on multiconfigurational and multirefernce methods
<http://dx.doi.org/10.1021/cr200137a>`_. A number of different flavors of the latter *a posteriori* corrections to the
energy have been implemented for the modules KRCI and LUCITA and for the nomenclature and detailed equations we refer to
the above review article (chapter 2.1.3.1). 

--------------
Example: Water
--------------

Following the example for water in the `review article <http://dx.doi.org/10.1021/cr200137a>`_ we need to run two
calculations. 

Step 1: Get reference wave function and configuration(s)
--------------------------------------------------------

We start first a multiconfigurational SCF calculation with a CAS(4,4) space using a cc-pVDZ basis set. This will provide
us a suitable reference wave function for the subsequent CI. The MCSCF module will in addition automatically create a file named **refvec.luci** which
contains the reference configurations needed for the *+Q* corrections. Our CASSCF input ``scf-casscf.inp`` is

.. literalinclude:: scf-casscf.inp

The molecular structure is located in ``h2o.xyz`` 

.. literalinclude:: h2o.xyz

Suppose that we have run the CASSCF calculation using::

   pam --inp=scf-casscf.inp  --mol=h2o.xyz --outkrmc --get="refvec.luci"

where we keep the MCSCF coefficients (in KRMCSCF) and the reference vectors stored in refvec.luci. 

Step 2: Run the MRCISD calculation and compute the +Q corrections
-----------------------------------------------------------------

In the second step we now run a MRCISD calculation using the reference wave function coefficients from the previous step
as starting point and compute the *+Q* corrections to the MRCISD energy. 
Using the same molecular input structure we run the MRCISD+Q correction using the input file ``krci-q.inp``

.. literalinclude:: krci-q.inp

with the following start command::

   pam --inp=scf-casscf.inp  --mol=h2o.xyz --inkrmc --put="refvec.luci"

We will find s summary of different *+Q* corrections printed at the end of the output: 

.. literalinclude:: summary_q.txt

A detailed account of the general *performance* of the different *+Q* corrections (and which other QC programs have
implemented which *+Q* correction by default) is given in chapter 2.1.3.1 of the above review article. We refer the user to carefully read this section. 
