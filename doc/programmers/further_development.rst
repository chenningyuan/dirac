Past and current works on DIRAC development
==========================================

Web-page dedicated to present and discuss various areas for the development of the DIRAC software.
This space is mainly for finished and on-going bachelor and master thesis of IT-students, plus general wishlist.


Master thesis
-------------

Comparison of CPU and GPU mathematical libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aim of this thesis is to compare the speed of some numerical procedures from CUBLAS and CULA libraries with their Fotran’s analogues BLAS and LAPACK libraries. 
Using these libraries we wrote a short test program (utils/gpu_test.F90) in programming languages Fortran and C, 
which compares the time complexity for various input parameters. 
Notice of this thesis is in future applying, which we tested, CUBLAS / CULA numerical library procedures in the entire system of quantum-chemical program DIRAC. 
(Finished in summer 2012, D.Kuzma, supervised by Miro Ilias.)


Bachelor thesis
---------------

CMake platform universal system for the software buildup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
First, learn the CMake buildup system for simple testing programs written in C,C++,Fortran.
Later, try to program more complex buildup CMake prescriptions, including macros.
Assigned student I.Hrasko. Miro Ilias supervising.

Documentation generator Sphinx and its usage for DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dokumentačný systém Sphinx (http://sphinx.pocoo.org), postavený na báze interpretačného jazyka Python, má v súčasnosti široké využitie pre tvorbu on-line manuálov rôznych softvérových projektov. Jedným z nich, kde Sphinx našiel svoje uplatnenie, je on-line dokumentácia vedeckého softvéru DIRAC (diracprogram.org). Existuje široký priestor na formálne  (i obsahové) vylepšenie dokumentácie k DIRAC-u, čo dáva náplň tejto bakalárskej práci.
Ciele: Popísať dokumentačný generátor Sphinx a s využitím rôznych jeho pokročilých funkcií vylepšiť dokumentáciu vedeckého programu DIRAC.


Adpatation of DIRAC for Microsoft Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adaptácia vedeckého softvéru DIRAC pre operačný systém Microsoft Windows.
Vďaka systému CMake (www.cmake.org) má vedecký softvér DIRAC (diracprogram.org), používaný v operačných systémoch Linux,Unix a Mac,  potenciál fungovať aj pod komerčným operačným systémom Microsoft Windows. Efektívne adaptovanie DIRACu pre MS Windows by rozšírilo možnosti  vedeckého počítania aj pod týmto sytémom.
Ciele: Skonfigurovať, skompilovať a spúštať DIRAC v prostredí MS Windows (od 7 vyššie). Upraviť rôzne parametre operačného systému MS Windows pre optimálny priebeh kvantovochemických výpočtov programu DIRAC.

IT-aspects of grid-computations with DIRAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Informačno-technologické aspekty výpočtov s vedeckým softvérom DIRAC v gridovom prostredí.
V súčasnosti gridové prostredie ponúka rozsiahly priestor zosieťovaných počítačov, kde možno spúšťať rôzne výpočtové úlohy. Nemej zaujímavé je vykonávanie vedeckých výpočtov v gride, v našom prípade z oblasti relativistickej kvantovej chémie, a to softvérom DIRAC (diracprogram.org). Aby testovacie a produkčné výpočty s DIRACom v gride boli vykonávané čo najefektívnejšie, treba poznať a správne ponastavovať viacero parametrov, či už na strane softvéru, a či v spúšťacích skriptoch pre grid.
Ciele: Nájsť najlepší spôsob pre spúšťanie programu DIRAC  v heterogénnom  gridovom prostredí.

Testing various mathematical libraries for DIRAC scientific software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Vplyv matematických knižníc na rýchlosť vedeckého softvéru DIRAC.
Ako mnohé vedecké programy, tak aj systém kvantovochemických programov DIRAC využíva pre svoj chod matematické knižnice BLAS a LAPACK. 
Nakoľko existuje viacero "producentov" týchto knižníc (uvádzame voľne dostupné GNU BLAS+LAPACK, GNU-ATLAS, Intel MKL), je vhodné odmerať rýchlosť vedeckých výpočtov v závislosti od druhu prilinkovaných knižníc a ešte od iných kompilačných parametrov.
Ciele: Zistiť, aký druh matematických knižníc a aké kompilačné parametre pre softvér DIRAC sú najvhodnejšie pre produkčné výpočty.

Git version control system and its usage in DIRAC software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Verzionovací systém Git a jeho uplatnenie vo vedeckom softvéri DIRAC.
Moderný spôsob vývoja softvéru zahrňuje aj ukladanie jeho jednotlivých verzií, k čomu slúži niekoľko verzionovacích nástrojov. Jedným z najviac používaných je  kontrolný systém verzií (angl. "version control system") Git, viď http://www.git-scm.com/. Git našiel svoje pevné miesto vo vývoji vedeckého softvéru DIRAC (diracprogram.org) a stále existuje priestor na zlepšenie používania Git-u v praktickom programovaní.
Ciele: Popísať verzionovací systém Git; ukázať jeho uplatnenie vo vedeckom projekte DIRAC.

Other TODO wishes
------------------

New input reader with documenting of keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let us write/use an keyword input reader which requires keywords to have documentation.
Otherwise the code won't even compile. Then the "make html" command would go through the source code and generate
the RST which would exactly reflect the state of the code (unless of course somebody changes the keyword without
changing the text 2 lines below but we could introduce some punishment for it). Radovan Bast.
