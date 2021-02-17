:orphan:


Converging atoms   
=======

DIRAC presently does not perform an open-shell Hartree--Fock in a strict sense,
rather an average-of-configuration calculation, which amounts to the
optimization of the average energy of a set of configurations (or determinants)
generated from the specification of a given number of open shells and their
electron occupations. In the case of atoms, a convergence problem can occur when
the inner open-shell orbitals are more stable than outer "closed" shell.
From DIRAC21, we can specify the occupation number of each :math:`\kappa` value. by using the keyword :ref:`SCF_.KPSELE`. 

As an example of such a calculation, we consider the Nb and Np atom.

-----------------------------------------------------
Example 1:  Nb [Kr]4d\ :sup:`4`\ 5s\ :sup:`1`\
-----------------------------------------------------
::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   18 18
  .OPEN SHELL
   2
   4/10,0
   1/2,0
  .KPSELE                                                                                             
   5                                                                                                   
   -1  1 -2  2 -3                                                                                      
    8  6 12  4  6                                                                                      
    0  0  0  4  6                                                                                      
    2  0  0  0  0   
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

--------------------------------
Example 2:  Np [Rn]5f\ :sup:`4`\6d\ :sup:`1`\7s\ :sup:`2`\
--------------------------------
::

  **DIRAC                                                                                             
  .WAVE FUNCTION                                                                                      
  .ANALYZE                                                                                            
  **HAMILTONIAN                                                                                       
  .X2C                                                                                                
  **INTEGRALS                                                                                     
  *READIN                                                                                             
  .UNCONTRACT                                                                                         
  **WAVE FUNCTION                                                                                     
  .SCF                                                                                                
  *SCF                                                                                                
  .CLOSED SHELL                                                                                       
   44 44                                                                                              
  .OPEN SHELL                                                                                         
  2                                                                                                   
  4/0,14                                                                                              
  1/10,0                                                                                              
  .KPSELE                                                                                             
  7                                                                                                   
  -1  1  -2  2 -3  3 -4                                                                               
  14 10  20 12 18  6  8                                                                               
   0  0   0  0  0  6  8                                                                               
   0  0   0  4  6  0  0                                                                               
  **ANALYZE                                                                                           
  .MULPOP                                                                                             
  *MULPOP                                                                                             
  .VECPOP                                                                                             
   1..oo                                                                                              
   1..oo                                                                                              
  *END OF

--------------------------------
Example 3:  Sorting the orbitals for Np
--------------------------------
The above approach is better, but the below one is interesting as a tutorial.
First, we remove the valence electrons and to compute the [Rn] configuration and to keep the DFCOEF file::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   42 44
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

This converges nicely.  Now we could be tempted to restart from DFCOEF,
reintroduce the  :math:`7s^{2}` and the open shells :math:`5f^{4}` and :math:`6d^{1}` and impose overlap
selection like this::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   44 44
  .OPEN SHELL
   2
   1/10,0
   4/0,14
  .OVLSEL
  .NODYNSEL
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

However this does not converge (please try). Why not? To understand this have a
look at the Mulliken population analysis of the [Rn] configuration (first
calculation).  In there you will see that gerade orbitals 22--26 are *d*, and
27 is an *s*.  It is important to know that DIRAC assumes that open-shell
orbitals come at the end, after closed-shell orbitals. In this case DIRAC would
occupy the *d* orbital 22 fully and then distribute one electron over 4 *d* and
1 *s* orbital. But this is not what we want. What we will do is to restart from
[Rn] DFCOEF and reorder the virtual gerade orbitals. We move the empty *d*
orbitals after the empty *s* orbital::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  .REORDER
  1..21,27,22,23,24,25,26

  *SCF
  .CLOSED SHELL
   44 44
  .OPEN SHELL
   2
   1/10,0
   4/0,14
  .OVLSEL
  .NODYNSEL
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

So in the above input we ask to keep gerade orbitals 1--21
as they are, then to move 22-26 after 27.
The empty line that is following means that ungerade
orbitals should be kept as they are.

Alright, we run this input, we keep the new DFCOEF
and in a final run we restart again, we remove reordering
(because orbitals are already reordered) and we relax
the overlap selection::

  **DIRAC
  .WAVE FUNCTION
  .ANALYZE
  **HAMILTONIAN
  .X2C
  **INTEGRALS
  *READIN
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *SCF
  .CLOSED SHELL
   44 44
  .OPEN SHELL
   2
   1/10,0
   4/0,14
  **ANALYZE
  .MULPOP
  *MULPOP
  .VECPOP
   1..oo
   1..oo
  *END OF

This converges nicely.
We can verify that the final valence orbitals are as they should be.
First the gerade::

  * Electronic eigenvalue no. 22: -0.2053405574542       (Occupation : f = 1.0000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np s
  --------------------------------------
   alpha    1.0000  |      1.0000
   beta     0.0000  |      0.0000

  * Electronic eigenvalue no. 23: -0.1931234363014       (Occupation : f = 0.1000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np dxx    L Ag Np dyy    L B1gNp dxy    L B2gNp dxz    L B3gNp dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.8000  |      0.2000         0.2000         0.4000         0.0000         0.0000
   beta     0.2000  |      0.0000         0.0000         0.0000         0.1000         0.1000

  * Electronic eigenvalue no. 24: -0.1931234361060       (Occupation : f = 0.1000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np dxx    L Ag Np dyy    L Ag Np dzz    L B2gNp dxz    L B3gNp dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.4000  |      0.0667         0.0667         0.2667         0.0000         0.0000
   beta     0.6000  |      0.0000         0.0000         0.0000         0.3000         0.3000

  * Electronic eigenvalue no. 25: -0.1829590945173       (Occupation : f = 0.1000)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np dxx    L Ag Np dyy    L Ag Np dzz    L B2gNp dxz    L B3gNp dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.6000  |      0.1000         0.1000         0.4000         0.0000         0.0000
   beta     0.4000  |      0.0000         0.0000         0.0000         0.2000         0.2000

  * Electronic eigenvalue no. 26: -0.1829590945233       (Occupation : f = 0.1000)  m_j=  5/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np dxx    L Ag Np dyy    L B1gNp dxy
  --------------------------------------------------------------------
   alpha    1.0000  |      0.2500         0.2500         0.5000
   beta     0.0000  |      0.0000         0.0000         0.0000

  * Electronic eigenvalue no. 27: -0.1829590945048       (Occupation : f = 0.1000)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L Ag Np dxx    L Ag Np dyy    L B1gNp dxy    L B2gNp dxz    L B3gNp dyz
  --------------------------------------------------------------------------------------------------
   alpha    0.2000  |      0.0500         0.0500         0.1000         0.0000         0.0000
   beta     0.8000  |      0.0000         0.0000         0.0000         0.4000         0.4000

And here is the 5f4 shell::

  * Electronic eigenvalue no. 23: -0.3986217183147       (Occupation : f = 0.2857)  m_j=  5/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B2uNp fxxy   L B2uNp fyyy   L B1uNp fxxz   L B1uNp fyyz   L Au Np fxyz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.1429  |      0.0000         0.0000         0.0000         0.0000         0.0357         0.0357         0.0714
   beta     0.8571  |      0.1071         0.3214         0.3214         0.1071         0.0000         0.0000         0.0000

  * Electronic eigenvalue no. 24: -0.3986217184250       (Occupation : f = 0.2857)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B3uNp fxzz   L B2uNp fxxy   L B2uNp fyyy   L B2uNp fyzz   L B1uNp fxxz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.7143  |      0.0000         0.0000         0.0000         0.0000         0.0000         0.0000         0.1786
   beta     0.2857  |      0.0214         0.0071         0.1143         0.0071         0.0214         0.1143         0.0000

  Gross  | L B1uNp fyyz   L Au Np fxyz
  --------------------------------------
   alpha |   0.1786         0.3571
   beta  |   0.0000         0.0000

  * Electronic eigenvalue no. 25: -0.3986217184451       (Occupation : f = 0.2857)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B3uNp fxzz   L B2uNp fxxy   L B2uNp fyyy   L B2uNp fyzz   L B1uNp fxxz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.4286  |      0.0000         0.0000         0.0000         0.0000         0.0000         0.0000         0.1286
   beta     0.5714  |      0.0429         0.0143         0.2286         0.0143         0.0429         0.2286         0.0000

  Gross  | L B1uNp fyyz   L B1uNp fzzz
  --------------------------------------
   alpha |   0.1286         0.1714
   beta  |   0.0000         0.0000

  * Electronic eigenvalue no. 26: -0.3652679300791       (Occupation : f = 0.2857)  m_j= -7/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B2uNp fxxy   L B2uNp fyyy
  -----------------------------------------------------------------------------------
   alpha    0.0000  |      0.0000         0.0000         0.0000         0.0000
   beta     1.0000  |      0.1250         0.3750         0.3750         0.1250

  * Electronic eigenvalue no. 27: -0.3652679301320       (Occupation : f = 0.2857)  m_j=  1/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B3uNp fxzz   L B2uNp fxxy   L B2uNp fyyy   L B2uNp fyzz   L B1uNp fxxz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.5714  |      0.0000         0.0000         0.0000         0.0000         0.0000         0.0000         0.1714
   beta     0.4286  |      0.0321         0.0107         0.1714         0.0107         0.0321         0.1714         0.0000

  Gross  | L B1uNp fyyz   L B1uNp fzzz
  --------------------------------------
   alpha |   0.1714         0.2286
   beta  |   0.0000         0.0000

  * Electronic eigenvalue no. 28: -0.3652679300598       (Occupation : f = 0.2857)  m_j= -3/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B3uNp fxzz   L B2uNp fxxy   L B2uNp fyyy   L B2uNp fyzz   L B1uNp fxxz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.2857  |      0.0000         0.0000         0.0000         0.0000         0.0000         0.0000         0.0714
   beta     0.7143  |      0.0536         0.0179         0.2857         0.0179         0.0536         0.2857         0.0000

  Gross  | L B1uNp fyyz   L Au Np fxyz
  --------------------------------------
   alpha |   0.0714         0.1429
   beta  |   0.0000         0.0000

  * Electronic eigenvalue no. 29: -0.3652679300957       (Occupation : f = 0.2857)  m_j=  5/2
  ==========================================================================================

  * Gross populations greater than 0.00010

  Gross     Total   |    L B3uNp fxxx   L B3uNp fxyy   L B2uNp fxxy   L B2uNp fyyy   L B1uNp fxxz   L B1uNp fyyz   L Au Np fxyz
  --------------------------------------------------------------------------------------------------------------------------------
   alpha    0.8571  |      0.0000         0.0000         0.0000         0.0000         0.2143         0.2143         0.4286
   beta     0.1429  |      0.0179         0.0536         0.0536         0.0179         0.0000         0.0000         0.0000
