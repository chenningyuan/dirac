:orphan:

Properties in solution: TDDFT calculations of excitation energies and solvent shifts
===================================================================================

In this tutorial we investigate property calculations in solution through 
polarizable embedding (PE). The tutorial builds on the implementation reported in 
:cite:`Hedegaard2017`. This implementation allows calculation of electric properties 
up to linear response, including an explicit environment (this corresponds to electric 
properties up to second order). Examples include electric dipole moments and 
polarizabilities. In addition, the linear response module allows calculation of 
excitation energies and transition momements, including the effect of a surrounding 
environment. Here we focus on TD-DFT includinga solvent. A few comments regarding 
calculation the solvent effect on dipole moments and polarizabilities are given in 
the end of this tutorial. ::  


TDDFT including an explicit environment with PE
-----------------------------------------------
The generalized eigenvalue problems that is solved when solving the response equations 

.. math::

   \left(E_{\text{tot}}^{[2]} -\hbar\lambda S^{[2]}\right)X_m=0

is modified by the pressence of an environment. The modification is manifested in an 
additional term, `E_{\text{pe}}^{[2]}`, to the electronic Hessian, `E_0^{[2]}` 
(cf. Eqs. 36-38) in :cite:`Hedegaard2017`

.. math::
 E^{[2]} =  E_0^{[2]} + E_{\text{pe}}^{[2]} = \left(\begin{array}{cc}
                                 \mathbf{A} & \mathbf{B} \\
                                 \mathbf{B}^* & \mathbf{A}^{*} 
                                 \end{array} \right) , 
with the matrix elements

.. math::

 A_{ai;bj} = \langle 0 \vert [\hat{q}^{\dagger}_{ai},[\hat{q}_{bj},\hat{f}_0 + \hat{v}^{\mathrm{PE}}]] \vert 0 \rangle + 
   \langle 0 \vert [\hat{q}^{\dagger}_{ai},\hat{v}^{\mathrm{ind}}]] \vert 0 \rangle \label{Hessian_2a} \\
 B_{ai;bj} = \langle 0 \vert [\hat{q}_{ai},[\hat{q}_{bj},\hat{f}_0 + \hat{v}^{\mathrm{PE}}]] \vert 0 \rangle + 
   \langle 0 \vert [\hat{q}_{ai},\hat{v}^{\mathrm{ind}}] \vert 0 \rangle . 

A derivation and more thorough explaination of the `E_{\text{pe}}^{[2]}` term can be 
found in cite:`Hedegaard2017`. cf. Eqs. 36-38) in :cite:`Hedegaard2017`) 

To demonstrate how to calculate the effect of the environment with PE, we use the same 
system as in the `PE-HF <pe_scf.rst>` tutorial. The employed h2o.mol file and potential 
3_h2o.pot can be found in this tutorial. The input pe_response.inp is similar to a 
regular TD-DFT calculation but with `.PEQM` under the `**HAMILTONIAN` flag, just as for 
`PE-HF <pe_scf.rst>` ::

   **DIRAC
   .TITLE
   H2O
   .WAVE FUNCTION
   .PROPERTIES
   **HAMILTONIAN
   .DFT
   B3LYP
   .PEQM     ! enable PElib 
   *PEQM     ! PE-related keywords 
   .POTENT   ! Specify what PE potential to use 
   3_h2o.pot ! name of the used PE potential 
   **INTEGRALS
   *TWOINT
   .SCREEN
   1.0D-12
   *READIN
   .UNCONT
   **WAVE FUNCTIONS
   .SCF
   *SCF
   .CLOSED SHELL
   10
   **PROPERTIES
   *EXCITATION ENERGIES
   .EXCITA
   1 5
   .INTENS
   0
   **END OF

As for the case of PE-HF, specification of `.PEQM` will thus make the program calculate the 
PE-specific terms; in case of linear response that is the `E_{\text{pe}}^{[2]}` term. 
Note that use of the `.GSPOL` keyword under the `*PEQM` flag will result in a calculation 
where only the term involving `\hat{v}^{\text{pe}}` is included (physcally, this corresponds 
to only allow polarization of the ground-state during the excitation process). Now we have 
all ingredients for starting the calculation::


   pam --inp=pe-hf.inp --mol=h2o.mol --put="3_h2o.pot" 


The summary will look identical to a vacuum calculation ::

   * Isotropic DL-DL non-zero oscillator strengths (f)                                                      
   ===================================================                                                      
   DL = dipole length                                                                                       
   Rate = Dipole radiation rate (s-1)   
   Lifetime = corresponding radiation lifetime (s)
       
   Level  Frequency (eV)    f           Rate         Lifetime    Symmetry                                   
   ------------------------------------------------------------------------                                 
       1      7.00316    0.000001    1.26019E+03    7.93528E-04                                             
       2      7.00317    0.000000    3.70356E+02    2.70010E-03       
       3      7.00317    0.000002    3.30199E+03    3.02848E-04                                             
       4      7.30170    0.049543    1.14615E+08    8.72484E-09                                             
       5      8.53008    0.000000    5.32729E+01    1.87713E-02                                             
   ------------------------------------------------------------------------                                 
   Sum of oscillator strengths:      0.04955

Note that the first three excitations have very low oscillator strength. These excitations 
corresponds to triplet excitations, which are automatically included in a relativistic 
framework. The fourth excitation is the observed one. To obtain the solvent shift, a 
calculation without PE must also be conducted. We collect the result of such a calculation 
in the table below together with the calculated and experimental solvent shifts. 




Electric dipole moments and polarizabilities
--------------------------------------------

Calculation of dipole moments and polarizabilities with PE is quite straightforward. We just need to add the PCM section to
a dipole moment calculation input file, e.g., for a ::

   **DIRAC
   .WAVE FUNCTION
   .PROPERTIES
   **HAMILTONIAN
   .PCM
   **WAVE FUNCTION
   .SCF
   **PROPERTIES
   .DIPOLE
   *END OF

This PCMSolver input file requires the CPCM model with cyclohexane as solvent. The solute is enclosed
in a cavity made from a single sphere, centered at the origin and with a 20 angstrom radius.

The output for the dipole moment calculation looks just like usual::

Electric dipole moment polarizabilities
---------------------------------------

We recall that the electric dipole polarizability is a second order property as is defined as
the second derivative, at zero field strength, of the energy:

Again, we just need to add the PCM section to an electric dipole polarizability calculation input file::

   **DIRAC
   .PROPERTIES
   **HAMILTONIAN
   .PCM
   **WAVE FUNCTION
   .SCF
   **PROPERTIES
   .POLARIZABILITY
   *END OF

the same .mol file as above for water is used in this example. The following PCMSolver input file is used::   

The output for the dipole moment calculation looks just like usual::

                         +--------------------------------+
                         ! Electric dipole polarizability !
                         +--------------------------------+
      
      
      
     @   Elements of the electric dipole polarizability tensor
     
     @   xx            7.76376447 a.u.   (converged)      
     @   yy            3.21397401 a.u.   (converged)      
     @   zz            5.57773363 a.u.   (converged)      
     
     @   average       5.51849070 a.u.
     @   anisotropy    3.941      a.u.
     
     @   xx            1.15047118 angstrom**3
     @   yy            0.47626180 angstrom**3
     @   zz            0.82653484 angstrom**3
     
     @   average       0.81775594 angstrom**3
     @   anisotropy    0.584      angstrom**3

----------------------------------------------


