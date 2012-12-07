

Calculation of NMR shieldings using simple magnetic balance
===========================================================

In this tutorial we look at 4-component relativistic calculations of NMR shielding constants. For best results we recommend to use the scheme of

simple magnetic balance (sMB) in conjunction with London orbitals <ref name="gosia:smb">.

{{article

| authors = Ma&#322;gorzata Olejniczak, Radovan Bast, Trond Saue, and Magdalena Pecul

| title   = A simple scheme for magnetic balance in four-component relativistic Kohn-Sham calculations of nuclear magnetic resonance shielding constants in a Gaussian basis

| journal = J. Chem. Phys.

| volume  = 136

| page    = 014108

| year    = 2012

| doi     = 10.1063/1.3671390

}}

</ref>


Basic theory of magnetic balance
--------------------------------

We encourage the reader to consult the original paper <ref name="gosia:smb"> </ref> for full details; here we just provide some key points.



In the absence of magnetic fields the coupling between the large and small components of the Dirac equation is

.. math::

  c \psi^S = \frac{1}{2m} \left[1 + \frac{E-V}{2mc^2} \right]^{-1} (\vec \sigma \cdot \vec p) \psi^L

The exact coupling can not be implemented ''per se'' at the basis set level due to the energy dependence of the coupling.

One therefore employs the non-relativistic limit of this coupling

.. math::

  \lim_{c \to \infty} c \psi^S = \frac{1}{2m} (\vec \sigma \cdot \vec p) \psi^L,

valid for positive-energy solutions (and extended nuclear charge distributions). Construction of the small component basis in this

manner assures the correct representation of the kinetic energy operator in the non-relativistic limit, and the coupling is therefore

referred to as '''kinetic balance''' <ref name="stanton:havriliak">

{{article

|  authors = R. E. Stanton and S. Havriliak

|  title   = Kinetic balance: A Partial solution to the problem of variational safety in Dirac calculations

|  journal = J. Chem. Phys.

|  volume  = 81

|  page    = 1910

|  year    = 1984

|  doi     = 10.1063/1.447865

}}

</ref>

. In a 2-spinor basis kinetic balance can be implemented in such a manner that there is a 1:1 ratio between the sizes of the large and small

component basis; this is referred to as '''restricted kinetic balance''' (RKB).

.. math::

  \chi_{\mu}^{2c;S} = (\vec \sigma \cdot \vec p) \chi_{\mu}^{2c;L}

In a scalar basis the straightforward implementation of RKB is not possible and one typically resorts to ''unrestricted kinetic balance'' (UKB)

.. math::

  \{ \chi_{\mu}^{1c;S} \} = \{ \vec p \chi_{\mu}^{1c;L} \},

basically generating the small component basis functions as derivative of the large component ones.


Magnetic fields are introduced through minimal substitution

.. math::

  \vec p \to \vec \pi = \vec p + e \vec A

which manifestly modifies the coupling between the large and the small components. The non-relativistic limit now reads

.. math::

  \lim_{c \to \infty} c \psi^S = \frac{1}{2m} (\vec \sigma \cdot \vec \pi) \psi^L

and is referred to as ''magnetic balance'' <ref name="aucar:dso">

{{article

|  authors = G. Aucar and T. Saue and H. J. Aa. Jensen and L. Visscher

|  title   = On the origin and contribution of the diamagnetic term in four-component relativistic calculations of magnetic properties

|  journal = J. Chem. Phys.

|  volume  = 110

|  page    = 6208

|  year    = 1999

|  doi     = 10.1063/1.479181

}}

</ref>.



Magnetic balance needs to be taken into account when constructing the first-order correction to the wave function upon inclusion of an

external magnetic field and is not assured in a RKB basis <ref name="pecul:shpol">

{{article

|  authors = M. Pecul, T. Saue, K. Ruud and A. Rizzo

|  title   = Electric field effects on the shielding constants of noble gases: A 4-component relativistic Hartree-Fock study

|  journal = J. Chem. Phys.

|  volume  = 121

|  page    = 3051

|  year    = 2004

|  doi     = 10.1063/1.1771635

}}

</ref>.



In the cited paper <ref name="gosia:smb"> </ref> we show that magnetic balance can be obtained in a Gaussian basis combining UKB with London orbitals. However, UKB generally leads to a small component basis larger than the large component basis set from

which it was generated and may increase computational cost as well as introduce linear dependencies in the basis set as compared to RKB. We therefore propose calculations of NMR shielding parameters by a simple two-step procedure described in the next section.



== A sample calculation ==

We will consider a 4-component relativistic density functional calculation of the NMR shielding constant of hydrogen fluoride using the KT2 functional

<ref name="KT2 functional">

{{article

|  authors = T. W. Keal and D. J. Tozer

|  title   = The exchange-correlation potential in Kohn--Sham nuclear magnetic resonance shielding calculation

|  journal = J. Chem. Phys.

|  volume  = 119

|  page    = 3015

|  year    = 2003

|  doi     = 10.1063/1.1590634

}}

</ref>.



The molecular input file ''HF.mol'' reads

.. literalinclude:: hf.mol

===Step 1: SCF calculation===



The first step is the generation of the Kohn-Sham wave function in a RKB basis. Although DIRAC employs scalar Gaussian basis functions restricted kinetic balance is by default imposed when transforming to orthonormal basis<ref name="visscher-saue-00">

{{article

| authors = L. Visscher and T. Saue

| title   = Approximate relativistic electronic structure methods based on the quaternion modified Dirac equation

| journal = J. Chem. Phys.

| volume  = 113

| page    = 3996

| year    = 2000

| doi     = 10.1063/1.1288371

}}

</ref>. The menu file ''step1.inp'' reads

.. literalinclude:: scf.inp

and we make sure to recover the coefficient file ''DFCOEF'' when running the calculation

.. literalinclude:: scf.run

Note that the above command refers to a parallel calculation with eight nodes.



=== Step2 : Response calculation ===



In the second step we start from the RKB coefficient obtained in the SCF calculation, but DIRAC will add the UKB complement

(keyword [[Manual:GENERAL#.RKBIMP|.RKBIMP]])



and

calculate the NMR shielding constant using unrestricted kinetic balance (keyword [[Manual:HAMILTONIAN#.URKBAL|.URKBAL]])

combined with London orbitals (keyword [[Manual:PROPERTIES:NMR#.LONDON|.LONDON]]). The menu file ''step2.inp'' accordingly reads

.. literalinclude:: response.inp

and the ''pam'' command is

.. literalinclude:: response.run

=== Comparison with RKB calculations ===



With the above computational procedure we obtain the following NMR shielding constants



            isotropic shielding

        ----------------------------

 atom          total          dia         para

 ---------------------------------------------------

 H            29.9502       -1.9992       31.9494

 F           415.8316       14.5348      401.2967

 ---------------------------------------------------



By not adding the UKB complement (deleting .RKBIMP and .URKBAL), that is, using RKB, gives



             isotropic shielding

         ----------------------------

 atom          total          dia         para

 --------------------------------------------------------

 H            28.8028       -1.9992       30.8020

 F           323.2672       14.5348      308.7323

 --------------------------------------------------------



If we in addition switch off the use of London-orbitals (deleting .LONDON) we obtain



             isotropic shielding

         ----------------------------

 atom          total          dia         para

 ------------------------------------------------------

 H            28.2920       20.2563        8.0357

 F           322.9065      387.1375      -64.2310

 ------------------------------------------------------



The experimental value for the absolute isotropic shielding available for H in HF is 28.72 ppm.<ref name="Chesnut-1994">.

{{article

| authors = D. B. Chesnut

| title   = Ab Initio Calculations of NMR Chemical Shielding

| journal = Annual Reports on NMR Spectroscopy

| volume  = 29

| page    = 71

| year    = 1994

| doi     = 10.1016/S0066-4103(08)60131-3

}}

</ref>.



== References ==



<references/>

.



In the cited paper
