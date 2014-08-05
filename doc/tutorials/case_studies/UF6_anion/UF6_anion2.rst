=========================================
The :math:`UF_6` anion - another approach
=========================================

Different and more robust solution for obtaining the SCF convergence of the :math:`UF_6` anion is provided here.

It is based on the damping of the Fock matrix from iteration to iteration.

In the first step we follow the previous approach of starting from the closed-shell UF6 system.
We also used a shift in the virtual space to reduce mixing between occupied and virtual spinors.

 pam --noarch --inp=UF6neutral.inp  --mol=UF6.mol --get "DFCOEF=DFCOEF.UF6neutral"

We compare two calculations starting from neutral molecular spinors of UF6. 

One calculation uses the DIIS method and the other uses the damping with the damping factor of 0.90.

 pam --noarch --inp=UF6_anion_scfdamp.inp --mol=UF6.mol --put "DFCOEF.UF6neutral=DFCOEF" --get "DFCOEF=DFCOEF.UF6anion_damp_pass1" --replace dampfactor=0.9
 pam --noarch --inp=UF6_anion_scfdamp.inp --mol=UF6.mol --put "DFCOEF.UF6anion_damp_pass1=DFCOEF" --get "DFCOEF=DFCOEF.UF6anion_damp_pass2" --replace dampfactor=0.5

With damping convergence speedup the calculation converges to 1,4 10^-8 after 50 iterations.

The damping scheme uses more of the input Fock matrix that the output Fock matrix.
Concerning the damping factor, the lower than 1 it is, we should get faster convergence. When the factor is getting too small,
iterations diverge. The closer we are to the converged states, the smaller the damping factor can be.

On the other hand the DIIS iterations initially converge to about 5*10^-5 but then began to diverge
and after 150 iterations the error gradient of the Fock matrix is cca 2.10^-5 .

However, if we began from the converged solutions from the damping, the DIIS iterations converge 
rapidly after few iterations. But with deactivated active-active shells interaction:

 pam --noarch --inp=UF6_anion2.inp --mol=UF6.mol --put "DFCOEF.UF6anion_damp_pass2=DFCOEF" --get "DFCOEF=DFCOEF.UF6anion_diis"

Also with OPENFAC=1 we get the convergence, not  with the default DIIS, but with damping factor
  pam --noarch --inp=UF6anion3.inp --mol=UF6.mol --put "DFCOEF.UF6anion_diis=DFCOEF"  
With this setting we did not get the convergence.

We try to get wacac convergence with damping
  oarsub -n basista -l nodes=1/cpu=4 ./runUF6m 


TODO: automatic adjustement o the damping factor, suggestions for keyword
.CNVDMP
1.D-4 1.D-5 1.0D-6
.DAMPDEL ! decrease level of damp. factor
0.2

