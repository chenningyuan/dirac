
Connection to Sham's modified Gamess program for making MO-transformed integrals.

Note on compilation: I had to lower some array sizes to reduce the size of common blocks below 2GB
Otherwise it would not be possible to use DIRAC makefile/dependencies system and one would have resort
to write own dependencies of each executable. With the so called medium memory model, flag 'mcmodel=medium'.

Also I had to rename some routines in some (not all) independent executable files because of overlaping their names.
Now all object files are used for linking executable files, therefore all function names must be unique.

These GAMESS files (produced due to Sham's changes) are read :
 HCRMO  $SCR/$JOB.F281
 MULMO  $SCR/$JOB.F282

I use the "ARCH32BIT" preprocessor flag to specify 32-bit architectures.

Efforts with old ACESII were not fruitfull - despite succesfull compilation there is some memory crash
at the beginning of the run.

Finally code has been adapted under spin-free DIRAC machinery.

Replacing static allocation with dynamic within modules:
========================================================

subroutine msub1 (ccsd1.F) - extra dynamic allocation is needed !
subroutine msub2 (ccsd1.F) - extra dynamic allocation is needed (separate module)
...  dimension ia(8),ib(10000),ibe(0:10000) ....

replace ipack4(i,j,k,l)=ior(ior(ish... with some module

I. Fock-space cc programs:
========================================
csd - closed shell      (0h,0p)-sector

csv - add 1 electron    (0,1)-sector
csu - add 2nd electron  (0,2)-sector
csw - add 3nd electron  (0,3)-sector
csx - add 4nd electron  (0,4)-sector
csy - add 5nd electron  (0,5)-sector
csz - add 6nd electron  (0,6)-sector

chv - ionize 1 electron    (1,0)-sector
chu - ionize 2nd electron  (2,0)-sector
chw - ionize 3nd electron  (3,0)-sector
chx - ionize 4nd electron  (4,0)-sector
chy - ionize 5nd electron  (5,0)-sector
chz - ionize 6nd electron  (6,0)-sector

chp,chq - excite 1 electron (1,1) sector; first run (0,1) and (1,0) sectors 

 miro&ephraim: not iterative scheme, some approximations, to be explored

    In (1,2) sector do "cst.f" then "m21b.f"
    In (2,1) sector for CCSD or CCSDT-1 do "cht.f" then "m12b.f"
    In (2,1) sector for CCSDT do "cht.f" only

TODO(miro): (2,2) FSCCSD 

II.    CCSD    available for (2,1)
       CCSDT   available for (0,0) (1,1) (1,2) 
                             (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) 
                             (1,0) (2,0) (3,0) (4,0) (5,0) (6,0) 
     
III. Files needed for compilation:
for all programs - c??.f csdsub.f
for csv-z, chv-z  - also eneru.f (look at the Makefile)


Miro: all compilations were moved under the DIRAC makefile system

IV. Formatted 'input' file (c**.in):
line 1 - title
line 2 - iopt(1:20) (20i3)
 iopt(1) - print level (normally 0) - use it more abudantly for control printouts !
 iopt(2) - max number of iterations.
 iopt(3) - convergence criterion. convergence reached when all
           t amplitudes agree to 10**(-iopt(3))

 iopt(4) .ne.0 - do rle (convergence accelerator) every iopt(4) iterations.
          good number is 4, must not be 1

 iopt(5) - t amplitudes (and other quantities, depending on iopt(1))
           printed only if elements gt 10**(-iopt(5))

 iopt(6)  = XYY - controls damping. iopt(6)=100*i61+i62, where i62 (or YY)
           is the iter. number on which damping starts (0 = no damping),
           i61 (or X) is damping factor:
            for i61 = 1,..,6 :  fnew=1/i61 ;
                    t(n,in)=[t(n-1,out)*1 + t(n-1,in)*i61]/(1+i61) ;
            for i61 > 6 :  fnew=i61-5 ;
                    t(n,in)=[t(n-1,out)*(i61-5) + t(n-1,in)*1]/(1+i61-5) ;
            default: i61=1 .

 iopt(7) - if .ne.0, continuation run. ???

 iopt(8) .eq.0 - for Hartree-Fock orbitals only: read only orbital energies
         .gt.0 - for arbitrary orbitals: read Fock matrix and orbital energies
         .eq.2 - (usual practice) - orb.energies are taken from diagonal
                 of Fock matrix.
         .le.0 - shifts for orb. energies (doesn't work properly)

 iopt(9) .ne.0 - incomplete model space, implemented fully in csv, csu;
                 partially in chp, chq (to delete whole symmetry blocks
                 of 1-electron excitations).

 iopt(10) must be = 0
 
    "sixv" and "sixh" do CCSD+T for the (0,n) and (n,0) sectors respectvely.

     iopt (11) - desired value of n (only needed for n.ge.2)

   "sixp" does CCSD+T for the (1,1) sector and must be preceded by "sixv" and
   "sixh", run in that order.
  
 iopt (16) = 1 : sixp does not follow
               = 2 : sixp does follow

 iopt (17) = 1 : CCSDT-1a               (noniterative ?)
           = 2 : CCSDT-1b (= CCSDT-1)
           = 3 : CCSDT
           = 0 : CCSD

 iopt (19) = 1 : ACES II (miro: not working now)
 iopt (19) = 0 : GAMESS (miro: very limited adaptation thanks to Sham's changed in GAMESS)
 iopt (19) = 7 : svn DIRAC (miro's working adaptation)


line 3 - norb = # of (spacial) orbitals in transformed integral file. (i3)
line 4 - cc(1:norb) = type of orbitals. following types are allowed:
                    c - inactive (applies to inactive orbitals appearing
                        in the integral file)
                    h - hole
                    p - particle
                    g - valence hole
                    v - valence particle
             note: the same input may be used for all programs, as programs
             which do not handle valence holes or particles interpret
             g as h or v as p.
line 5 - F(i3) NNNN = number of symmetries
line 6 - F(24i3) NORB(I), I=1,NNNN - number of orbitals in each symmetry
line 7 - F(24i3) L(I), I=1,NNNN - L=0 or 1 depend on symmetry group
(one line for every symmetry operation. indicators for orb groups, 0 if
group symmetric under operation, 1 if group antisymmetric (20i3).)
last line - number '9' 

next lines (if are presented):

incomplete model space info (for iopt(9) ne 0):

in (chv) csv (subroutine rdrsp1):
      read (5,*) nrsp1,((it(j,i),j=1,2),i=1,nrsp1)
   where a transition from it(1,i) to it(2,i) (both valence orbitals)
   takes you from p space to q space.

in (chu) csu (subroutine rdrsp2):
      read (5,*) nrsp2,((it(j,i),j=1,2),i=1,nrsp2)
   where (it(1,i),it(2,i)) defines an r space determinant (should be in p
   but isn't)

Additional info: incomplete model space input, if (iopt(9).ne.0). csv and csu only.
r-space is the part of the would-be complete model space that is now i
q space. e.g., if 3 and 4 are valence orbs, and p space includes 33 and
34 only, 44 is r space. in csv, input number of 1-electron transitions
which can move from p to r, and then the transition (free format). in
the example, input is
1, 3,4
because 3 to 4 can take 34 to 44. in csu, put number of r-space determinants,
then determinants. in the example,
1, 4,4.
c ------------------------------------------------------------------------------------

=========================================================================
V. SOME INPUT EXAMPLES AND EXPLANATIONS:
=========================================================================

 ***  csd  *** file isdnh4p3
 -------------
nh43. nh4+3. r=3.
  0 30  5  0  2
 41
hhvvpppppppppppppppppppppppphpppppppppppp
  2
 28 13
  0  1
9
 
 
 ***  csv  *** file isvnh4p3
 -------------
nh43. nh4+3. r=3.
  0 30  5  0  2  0  0  0  1
 41
hhvvpppppppppppppppppppppppphpppppppppppp
  2
 28 13
  0  1
9
1, 3,4
 
 
 ***  csu  *** file isunh4p3
 -------------
nh43. nh4+3. r=3.
  0 30  5  0  2  0  0  0  1
 41
hhvvpppppppppppppppppppppppphpppppppppppp
  2
 28 13
  0  1
9
1, 4,4
 
---------------------------------------------------
another example. li atom, li+ reference state, 2s and 2p are valence orbs,
2p**2 not in p space (incomplete model space)
 
  **  file isdli
li+.
  0 50  6
 26
hvppppppvppppvpppppvpppppp
  7
  8  5  5  1  5  1  1
  0  1  0  1  0  1  0
  0  0  1  1  0  0  1
  0  0  0  0  1  1  1
9
 
  **  file isvli
li+.
  0 50  6  0  0  3
 26
hvppppppvppppvpppppvpppppp
  7
  8  5  5  1  5  1  1
  0  1  0  1  0  1  0
  0  0  1  1  0  0  1
  0  0  0  0  1  1  1
9
  --  --  --  --
note: while this is an incomplete model space, the p-to-r 1-electron
transitions are 2s to 2p, which are symmetry forbidden. isvli looks
therefore like a complete space run.
 
  **  file isuli
li+.
  0 50  6  0  0  3  0  0  1
 26
hvppppppvppppvpppppvpppppp
  7
  8  5  5  1  5  1  1
  0  1  0  1  0  1  0
  0  0  1  1  0  0  1
  0  0  0  0  1  1  1
9
9, 9,9, 9,14, 9,20, 14,9, 14,14, 14,20, 20,9, 20,14, 20,20


Mg2+ testing system to mimic spinfree DIRAC FSCCSD 00,01,11,10,22 results
(basis set from 19.test)

 Mg test with DIRAC-HSFS
  0 99 08  4  2  0  0  0  0  0  0  0  0  0  0  0  0  0  7
 20
hhvvppppgvpgvppgvppp
  8
  8  3  3  1  3  1  1  0
  0  0  0  0  0  0  0  0
9

run as "hs_fscc.x csd csv csu chv chu chp chq > out"

CCSD 02: -0.13113488967157
relcc  : -0.131134892087296

01:   -0.52345057E+00. degen  2.
      -0.523450570249
      -0.22231511E+00. degen  2
      -0.222315110263
      0.60482483E+00. degen  6.
      0.604824830845

02:   -0.77584492E+00. degen  1.
      -0.775844923290
      -0.60277004E+00. degen  3.
      -0.602770039654
      -0.58644739E+00. degen  1.
      -0.586447396947
      -0.32064907E+00. degen  1.
      -0.320649068442

20:    0.68347932E+01. degen  9.
       6.834793186609
       0.70102683E+01. degen  5.
       7.010268354366
       0.71926553E+01. degen  1.
       7.192655334165

11:    0.19901840E+01. degen  9.
       1.989659151079
       0.20075755E+01. degen  3.
       2.007207787789
       0.25112430E+01. degen  9.
       2.511157022314 

__________________________________________________________________
__________________________________________________________________
 
VI. Reading the output:
******************
tmax and tmax(rle), largest t2 amplitudes, are printed for each iteration,
together with their identity (hexadecimal). in csd, csv, and chv, largest
t1 also appears. largest t differences relative to preceding iteration
 are shown. message after end of iterations may be "converged after n
iterations" or "iteration limit". last t amplitudes and energies are then
given.
what is the energy?
******************
the closed-shell energy is e(scf) + correlation energy from csd.
to get open-shell energies, add to that the values from the open-shell
programs: csv for one additional e, csu for two, chv for ionizing one e,
chu for two, chq for exciting one e.
chq gives the dominant determinant for each eigenvalue. this is described
by the orbital which lost an electron and the one which gained an electron.
minus indicates beta spin. in csv, csu, chv, chu, the determinants are given
by the valence spinorbs. 2n-1 is orb n with spin alpha, 2n - with spin beta.

__________________________________________________________________
VII. Some parameters info:
C intrc = 600 - old size of the integral buffer
c intrc = 9600 - new size of the integral buffer (from MOLCAS)
C nort - number of spatial MO
C norb - number of molecular spin-orbitals (SO-ls)
C eps(io) - orbital energies
c
c "io" - spin-orbitals (SO-ls) have continuous numeration
c        through all symmetry classes 
c
c    Sectors: t(m,n), (m<=6 - valence holes, n<=6 - valence particles)
c ic(io) = 1 - "c" (inactive)
c        = 2 - "h" (hole), "g" (valence hole)
c        = 3 - "p" (particle), "v" (valence particle)
c icv(io) = 0,1; so that:
c 2*ic(io)-1+icv(io) = 1 - "c" (inactive)
c                    = 2 - "c" (inactive)
c                    = 3 - "h" (hole)
c                    = 4 - "g" (valence hole)
c                    = 5 - "p" (particle)
c                    = 6 - "v" (valence particle)
c icm(io) = icl (# of irrep)
c
c nsym - number of basic symmetry operations
C nclass - number of the symmetry classes (irreps)
c iclass(icl) - number of the orbitals in "icl"-irrep
c isym(icl,i) - multiplication table for irreps (sym. classes)
c ibias(icl) + 1 = io - the number of the first spin-orbital in "icl"-irrep
c
c    *** Be carefull!!! ***


===============================================================================

 VIII. SJOB.trle1, $JOB.trle2, $JOB.vrle, $JOB.rleps - for RLE acceleration
 of convergency (Devidson-like procedure, see: Purvis III & Bartlett,
 JCP 75(3), 1284 (1981). These files are temporary and they can be
 removed before running the next program.

========================================
IX. "param.f" file  parameters description: 
 (miro:  renamed to "param_ccsd.h")
========================================


ndis (=      90 000) - number of symmetry blocks in PseudoFile (PF) [PF - array
                       of <ij|kl> integrals of some definite structure,
                       like  HHPP or PPPP  etc, where H - hole & P - particle].
                     
ndij  (=  1 500 000) - number of ij-blocks in PF , where ij - the first two
                       indexes in <ij|kl>.   <See:  common /cij1/  higher>

nbig (= 70 000 000) - the dimension of the array igen(), where all integral
                       indexes of all PF are placed.  <Used in: common /i/
                       igen(nbigi), subr. dimx( ... ,nbig, ...)>

nbig  (= 65 000 000) - the dimension of the array  gen(), where all integral
                       values of all PF are placed.  <Used in: common //
                       gen(nbig), subr. dimx( ... ,nbig, ...)>

ndp4  (= 12 000 000) - max size of HPPP file.

ndt2  (=  3 000 000) - max size of HHPP file.

ndp41 (=     50 000) - max size of PP file.

ndp44 (= 39 000 000) - max size of PPPP file.

ndint (= 20 000 000) - max size of original ordered 2-el. integrals, which was
                       created by the Hess program (with spin-orbitals!)
                       Int-ls >= 1.d-8, for any nbl (in fact, Coulomb int-ls).

ndt1  (=     50 000) - max size of original ordered 1-el. integrals, which was
                       created by the Hess program.

