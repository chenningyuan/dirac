Weak interactions in atoms and molecules
========================================


P-odd interaction constants
---------------------------

Demonstration of DC-SCF and DC-DFT calculations of the :math:`W_{A}` parameter of the P-odd spin-rotational Hamiltonian
for small diatomic molecule.

For more details see Ref. [Borschevsky2012]_ .


Other interactions
-------------------

Also experimental tests of the e_EDM property.


B atom 
~~~~~~~

4s2p STO-2G decontracted basis set

seeking <1s_1/2|H_e_EDM|2p_12> elements

D2h:
Au, x, (1,9)             :  0.0003366642        0.0000000000
Au, y, (1,9)             :  0.0003366642        0.0000000000
Au, z, (1,9)             :  0.0003366642        0.0000000000

C2v:
B2, x, (1,13) -> (1a,7a) :  0.0000000000       -0.0003366642
A1, y, (1,3)             : -0.0003366642        0.0000000000
B1, z, (1,13)            :  0.0003366642        0.0000000000

C2:
x, (1,2)                 :  0.0003354219        0.0000288954 => sqrt(re*re+im*im)=0.0003366642
y, (1,13)                :  0.0000288954       -0.0003354219
z, (1,3)                 : -0.0000288954        0.0003354219

C1: 
x, (1,3)                 : -0.0001381189       -0.0002682541  => sqrt(re*re+im*im)=.000301723 ?!
y, (1,3)                 : -0.0001381189       -0.0002682541
z, (1,3)                 : -0.0001381189       -0.0002682541

RPA,D2h, nonzero terms:


Na atom
~~~~~~~
SCF energy:  -155.866369658577725

RPA, D2h, nonzero terms:
   Propagator            Frequency             Real part        Imaginary part
@ <<e_EDM x ;ZDIPLEN >> (   0.00000)        0.0000000000000    -0.0026007587066
@ <<e_EDM y ;ZDIPLEN >> (   0.00000)        0.0000000000000    -0.0026007587066
@ <<e_EDM z ;ZDIPLEN >> (   0.00000)        0.0000000000000    -0.0026329560412

RPA,C2v, nonzero terms:
@ <<e_EDM x ;XDIPLEN >> (   0.00000)       -0.0000000000000    -0.0026003598269
@ <<e_EDM z ;YDIPLEN >> (   0.00000)        0.0000000000000     0.0025681624894

RPA, C1, nonzero terms:
@ <<e_EDM x ;ZDIPLEN >> (   0.00000)       -0.0000000000000     0.0026007587060
@ <<e_EDM y ;ZDIPLEN >> (   0.00000)       -0.0000000000000     0.0026007587060
@ <<e_EDM z ;ZDIPLEN >> (   0.00000)        0.0000000000000     0.0026329560407

Li atom
~~~~~~~
SCF energy : -7.109472981377731

RPA, D2h
   Propagator            Frequency             Real part        Imaginary part
@ <<e_EDM x ;ZDIPLEN >> (   0.00000)       -0.0000000000000    -0.0000980446380
@ <<e_EDM y ;ZDIPLEN >> (   0.00000)       -0.0000000000000    -0.0000980446380
@ <<e_EDM z ;ZDIPLEN >> (   0.00000)       -0.0000000000000    -0.0000981969272 

C2v
@ <<e_EDM x ;XDIPLEN >> (   0.00000)       -0.0000000000000    -0.0000980434369
@ <<e_EDM z ;YDIPLEN >> (   0.00000)        0.0000000000000     0.0000978911478

C1
@ <<e_EDM x ;ZDIPLEN >> (   0.00000)       -0.0000000000000     0.0000980446380
@ <<e_EDM y ;ZDIPLEN >> (   0.00000)       -0.0000000000000     0.0000980446380
@ <<e_EDM z ;ZDIPLEN >> (   0.00000)       -0.0000000000000     0.0000981969271
