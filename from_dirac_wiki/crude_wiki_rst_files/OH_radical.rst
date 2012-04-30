OH radical
==========

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Hydroxyl radical <#Hydroxyl_radical>`_

   -  `1.1 Settig the SCF
      configuration <#Settig_the_SCF_configuration>`_

      -  `1.1.1 Spin-free approach <#Spin-free_approach>`_

   -  `1.2 Mj selection <#Mj_selection>`_
   -  `1.3 Correlated approach <#Correlated_approach>`_

-  `2 Hydroxyl radical cation OH+ <#Hydroxyl_radical_cation_OH.2B>`_

Hydroxyl radical
----------------

Here we present open-shell SCF calculations on the OH.(known as the
hydroxyl radical) system, followed by complete small open shell CI
method (cosci) for getting individual states.

The ground state of the OH readical is X 2\_Pi\_3/2 resulting from the
electron configuration of

::

    1s_sig(2) 2s_sig(2) 2pz_sig(2) 2pxy_pi_1/2(2) 2pxy_pi_3/2(1)

Settig the SCF configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

What does not converge is the 1 electron in 1 open-shell:

::

    .CLOSED SHELL
     8
    .OPEN S
    1
    1/2

Therefore for this case one has to resort to the
average-of-configurations (3 electrons over 4 pi shells):

::

    .CLOSED SHELL
     6
    .OPEN S
    1
    3/4

The Dirac-Coulomb Hamiltonian gives noncorrelated fine structure
splitting of 143.842458 cm-1 as the result of complete open-shell CI
over pi-shells.

::

    pam OH.ccpVDZ.lsym.mol OH.dc_rkb.scf_os_res.inp
    pam OH.ccpVDZ.C2v.mol OH.dc_rkb.scf_os_res.inp

COSCI energies are :

::

     -75.446610400724 (   2 * )
     -75.445955006263 (   2 * )

The average SCF energy is

::

    -75.446282703494589 

and is equal to the average values of COSCI energies

::

    ((-75.446610400724-75.445955006263)/2=-75.446282703493

Spin-free approach
^^^^^^^^^^^^^^^^^^

Utilizing the Dirac-Coulomb spin-free Hamiltonian confirms the four-fold
(pi\_x=pi\_y) degeneracy of the valence open shells:

::

    pam OH.ccpVDZ.C2v.mol OH.dc_sf.scf_os_res.inp

With the resulting spin-free SCF (and COSCI) energy of

::

     -75.446280480624 (   4 * ).

Mj selection
~~~~~~~~~~~~

For this heteronuclear diatomic molecule one has the advantage of using
the omega-number occupation of shells. One can even place one electron
in one shell, because the ".MJSELECTION" keyword ensures the convergence
of both "X 2Pi\_3/2"

::

    .CLOSED SHELL
     8
    .OPEN S
    1
    1/2
    # 2Pi_3/2
    .MJSELECTION
    3
    4 0 0
    0 1 0

and of the "A 2Pi\_1/2" first excited state:

::

    .CLOSED SHELL
     8
    .OPEN S
    1
    1/2
    # 2Pi_1/2
    .MJSELECTION
    3
    3 1 0
    1 0 0

::

    pam OH.ccpVDZ.lsym.mol OH.dc.scf_mj_2Pi32.inp
    pam OH.ccpVDZ.lsym.mol OH.dc.scf_mj_2Pi12.inp

giving total SCF energies of the ground and the excited states:

::

    X 2Pi_3/2 :  -75.446717125968732
    A 2Pi_1/2 :  -75.446063726708829

The fine-structure splitting then makes

::

    FSS : (-75.446717125968732+75.446063726708829)*219474.631280634 = 143.404561 cm-1

which is very close to the above given COSCI value of 143.842458 cm-1.

Energy order of spinors for the OH "X 2Pi\_3/2" ground state is:

::

    1: -20.623437096718       (Occupation : f = 1.0000)  m_j=  1/2 ; 1s_sig(2)
    2: -1.2970494235574       (Occupation : f = 1.0000)  m_j=  1/2 ; 2s_sig(2)
    3: -0.6455362075645       (Occupation : f = 1.0000)  m_j=  1/2 ; 2pz_sig(2)
    4: -0.5439526090110       (Occupation : f = 1.0000)  m_j=  1/2 ; 2pxy_pi_1/2(2) 
    5: -0.5863405891346       (Occupation : f = 0.5000)  m_j= -3/2 ; 2pxy_pi_3/2(1)

While for the first excited "A 2Pi\_1/2" state of the OH radical we have
this order:

::

    1: -20.623901392971       (Occupation : f = 1.0000)  m_j=  1/2  ;  1s_sig(2)
    2: -1.2971646385507       (Occupation : f = 1.0000)  m_j=  1/2  ;  2s_sig(2)
    3: -0.6456515101399       (Occupation : f = 1.0000)  m_j=  1/2  ;  2pz_sig(2)
    4: -0.5431435116553       (Occupation : f = 1.0000)  m_j= -3/2  ;  2pxy_pi_3/2(2)
    5: -0.5874674338778       (Occupation : f = 0.5000)  m_j=  1/2  ;  2pxy_pi_1/2(1)

Correlated approach
~~~~~~~~~~~~~~~~~~~

The user can try to approach experimental
`data <http://webbook.nist.gov/cgi/cbook.cgi?ID=C3352576&Units=SI&Mask=1000#Diatomic>`_.

Hydroxyl radical cation OH+
---------------------------

Again one employs the average-of-configuration:

::

    .CLOSED SHELL
     6
    .OPEN S
    1
    2/4

Number of electronic states is higher.

Electronic states can be found
`here <http://webbook.nist.gov/cgi/cbook.cgi?ID=C12259299&Units=SI&Mask=1000#Diatomic>`_.

