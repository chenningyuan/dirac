# Design note: AO centers as the primitive for atomic IAO fragments

**Status:** implemented (steps 2a and 2b); see sections 7 and 8
**Scope:** prerequisite for labelled-XYZ fragment input; behaviour-neutral for conventional systems
**Relates to:** step 1 (fragment/group disentanglement, already implemented)

---

## 1. Motivation

The step-1 change separated the two roles that `n_fragments` used to play, giving us
an IAO-fragment layer and a recanonisation-group layer connected by `frag_to_group(:)`.
The next question is how the user supplies the grouping.

A labelled XYZ file (`C  0.0 0.0 0.0  methane`) is far friendlier than a list of
integers, but it introduces a second copy of the coordinates. Before adding it we need
to settle what a *center* is in ROSE, because the answer determines whether that second
copy is harmless or a source of silent inconsistency.

The conclusion of this note is that ROSE's real primitive is the **AO center**, not the
nucleus; that enumerating atomic fragments over AO centers removes the ambiguity, fixes
an existing bug, and lets the coordinates be read exactly once; and that the labelled
XYZ can then be introduced safely as a pure lookup table.

---

## 2. Analysis of the current code

### 2.1 The fragment mapping is purely geometric

Fragment-to-supermolecule matching is done by analytic Gaussian overlap, using nothing
but the shell coordinates. From `get_overlap.F90` (lines 169-180):

```fortran
xl = ao_basis_ket%gtos(lsh)%coord(1)
yl = ao_basis_ket%gtos(lsh)%coord(2)
zl = ao_basis_ket%gtos(lsh)%coord(3)
...
xk = ao_basis_bra%gtos(ksh)%coord(1)
```

There is no index join and no atom-number lookup. The field `gtos%atom_number` is
written only by `fchk_interface.F90:418` and is never read back for fragment mapping --
only by fchk's own writer (lines 564, 602, 676). It is absent or meaningless for the
other interfaces.

**Therefore: the coordinate of an AO shell is the only identity a center has, as far as
the IAO construction is concerned.**

### 2.2 What the nucleus list actually does

In the atomic-fragment path the nucleus list (`molecule%firstatom`) does four separable
jobs. Only the last is genuinely about nuclei.

| Job | Site | What it is really asking |
|---|---|---|
| Enumeration | `import_fragments.F90:26` -- `n_fragments = molecule%natoms` | "how many centers carry a basis?" |
| Placement | `import_fragments.F90:123,127` -- `gtos(j)%coord = fragments(i)%firstatom%x` | correction, because atomic library files store the atom at the origin |
| Identity | `import_fragments.F90:117` -- `nint(fragments(i)%firstatom%bascharge)` | which atomic library file to read |
| Charge bookkeeping | `properties.F90:45,126`; `set_coreval_info` | genuinely about nuclei |

Note that the identity job already uses `bascharge` rather than `charge`: the code is
half-aware that *basis* identity is the relevant quantity. And the analysis routines in
`properties.F90` touch only `fragments(i)%Znumber` and `fragments(i)%name` -- a scalar
and a string. No coordinate reaches the analysis.

### 2.3 The information doubling already exists

`hdf5_interface.F90` reads two independent datasets from the same checkpoint:

```fortran
call lab_read(datafile, '/input/molecule/geometry', rdata=geometry)   ! line  91  -> nuclei
call lab_read(datafile, ao_label//'center',         rdata=coord)      ! line 311  -> AO centers
```

ROSE then assumes, without ever checking, that nucleus *i* owns the shells sitting at
nucleus *i*'s position. The two lists are never reconciled.

So the doubling that a labelled XYZ would appear to introduce is already present. An
external file would be a *third* copy, and would make an existing latent inconsistency
visible rather than create a new one.

### 2.4 Ghosts are handled; point nuclei are not

`nuc_charge` and `bas_charge` (`hdf5_interface.F90:84-89`) already encode the
distinction the user cares about. The handling is asymmetric.

**Ghost center** (nuclear charge 0, basis charge Z) -- correct throughout. Named
`Gh(X)` at `import_fragments.F90:43`, and `set_coreval_info` assigns it
`ncore = 0, nvale = 0`. That is right: a ghost has basis functions but no electrons, so
it contributes no occupied reference orbitals to the IAO minimal basis, while remaining
available as a projection target.

**Point nucleus** (nuclear charge Z, basis charge 0) -- **broken**. It is still counted
at line 26, so it receives a fragment slot. `set_coreval_info` falls through to
`atom_minbas_size(nint(0.0), ...)`, and line 117 then writes `'000'` as the library
filename for a file that does not exist. Nothing filters it out, because the enumeration
runs over nuclei rather than over centers that carry a basis.

---

## 3. Proposed design

### 3.1 Enumerate atomic IAO fragments over distinct AO centers

Replace `n_fragments = molecule%natoms` with an enumeration derived from the AO basis
itself: collapse `ao_basis%gtos(1:nshells)%coord` onto distinct centers within a
tolerance, producing

- `n_centers` -- the number of distinct AO centers,
- `center_coord(3, n_centers)` -- their coordinates,
- `shell_to_center(1:nshells)` -- the map from shell to center.

`shell_to_center` is the important product. The IAO block structure is currently built
on an *implicit* assumption about shell ordering; making the map explicit is what allows
the block structure to be derived rather than assumed.

Consequences, all of which fall out without further work:

- **Coordinates are read once.** They arrive as an intrinsic part of the orbital
  definition. There is nowhere else they could come from.
- **The placement overwrite at lines 123/127 keeps its job but changes its source.**
  (Corrected during implementation -- an earlier draft of this note wrongly said the
  overwrite becomes unnecessary.) `ao_basis_2(i)` is the *atomic library* basis, which
  genuinely does store the atom at the origin and must still be translated. What changes
  is only the target coordinate: the AO center instead of the nucleus. In practice this
  is done by setting `fragments(i)%firstatom%x` from the center list, so line 123 is
  untouched and simply picks up the better coordinate.
- **Point nuclei drop out automatically** -- no basis means no center means no fragment.
- **Ghosts stay in automatically** -- they have a basis, so they are centers.
- **The nucleus list is demoted** to charge bookkeeping in the analysis, plus supplying
  `bascharge` for the library lookup via a coordinate match against the center list.

### 3.2 New tolerance parameters

Two module parameters in `read_xyzfile.F90`, alongside the existing `ghost_threshold`
(line 20). Both are declared directly in the internal units, so no conversion happens
anywhere:

```fortran
REAL(KIND=8), PARAMETER :: ghost_threshold  = 1.0D-10   ! unchanged
REAL(KIND=8), PARAMETER :: center_match_tol = 2.0D-2    ! bohr
REAL(KIND=8), PARAMETER :: charge_match_tol = 1.0D-2    ! electron
```

Declaring `center_match_tol` in bohr rather than Angstrom removes the `/xtang`
conversion and the risk of it being lost at a call site. Internal coordinates are
already in bohr -- both `read_xyz_atoms` and `hdf5_interface` apply `xyz = xyz / xtang`
on input -- so the comparison is direct.

**Why 0.02 and not 0.01.** This is a same-center identity test, not a chemical
tolerance, so it wants to be tight -- but it also has to survive hand-typed coordinates
in the labelled XYZ of section 4. A coordinate written with only two decimals in Angstrom
carries a rounding error of up to 0.005 Angstrom per component, or 0.0094 bohr; under a
Euclidean test the worst-case three-component error is sqrt(3) x 0.0094 = 0.0163 bohr.
At 0.01 bohr that case fails to match; at 0.02 bohr (about 0.0106 Angstrom) it matches
with room to spare, while still being far tighter than any real interatomic distance.

The metric is Euclidean distance, implemented as a squared-distance comparison in
`match_center`. Both boundary cases are covered by the unit test described in section 7.

### 3.3 Relationship to `ghost_threshold` -- decided

`ghost_threshold = 1.0D-10` is an "is this exactly zero" test, used at
`import_fragments.F90:43` and `:502` and `rose_structure.F90:57`.

`charge_match_tol = 0.01` is a different notion: a physically meaningful tolerance on
charge comparison, which would also accommodate fractional or effective charges.

**Decided:** `ghost_threshold` stays unchanged, and `charge_match_tol` is used only for
the new comparisons introduced by this design. Redefining `ghost_threshold` to 0.01 would
change the classification of any center with a small but nonzero charge; that is a
behaviour change to take deliberately and separately, not as a side effect of this work.

## 4. The labelled XYZ, made safe

Once AO centers are ground truth, the labelled file stops being a geometry source and
becomes a **lookup table from coordinate to label**. The governing rule:

> A coordinate read from the label file is used only as a join key. It is matched
> against the AO center list within `center_match_tol`, the label is taken, and the
> coordinate is then discarded. It is never copied into any data structure.

This removes the "which coordinate wins" question entirely, because the file's
coordinates never enter the calculation.

Error handling follows directly:

- a center matched by no XYZ line -> error naming the center's coordinate;
- an XYZ line matching no center -> error naming the line;
- an XYZ line matching more than one center -> error; the tolerance is too loose for
  this geometry.

This is redundancy as a *checked join key*, not redundancy as a second source of truth.
It is also what makes the format friendly: the user already wrote that XYZ to run the
SCF, so they add one column to a file they already have, and ROSE tells them immediately
if it does not line up with the orbitals.

---

## 5. Sequencing

**Step 2a -- center-derived enumeration (behaviour-neutral).**
Introduce the center-collapse helper, the tolerance parameters, and `shell_to_center`.
Switch the atomic-fragment enumeration onto the center list. This also requires reordering
the main driver so the AO basis is read before the fragments are defined (see section 7). For any system in which
every nucleus carries a basis and no ghosts are present, this produces an identical list
in an identical order, so it is a no-op -- while fixing the point-nucleus case and
removing the placement overwrite.

**Step 2b -- labelled XYZ input.**
Add the extended-XYZ reader, the coordinate-to-label join, and the keyword. This lands
on a foundation where "center" already means one unambiguous thing, so the file does not
have to carry the disambiguation itself.

Splitting the work this way keeps step 2a reviewable against existing test cases, since
any change in its output is by definition a bug.

---

## 6. Decisions and remaining questions

### 6.1 Settled at review

1. **`ghost_threshold` stays at `1.0D-10`.** See 3.3.

2. **Tolerances are declared in bohr and electrons directly** -- `center_match_tol =
   1.0D-2` bohr, `charge_match_tol = 1.0D-2` e. No unit conversion anywhere. See 3.2 for
   the precision consequence this has for hand-written labelled XYZ files.

3. **Centers are numbered by first appearance** in the shell list -- but that numbering
   is *not* what orders the fragments. See section 7: the premise behind this decision
   turned out to be wrong, and the fragment ordering is taken from the nucleus list
   instead. `collapse_ao_centers` still numbers centers by first appearance, since that
   ordering is internal to the map and never reaches the fragment numbering.

4. **Small-component consistency is checked, at no extra cost.** See 6.2.

### 6.2 Small-component check -- the check is the map

`gtos_small` has its own shell coordinates when `has_small_component` is true. Rather
than collapsing them into a second center list and comparing the two, map each small
shell onto the *large-component* center list using the same matcher:

```fortran
if (ao_basis%has_small_component) then
   do ishell = 1, ao_basis%nshells_small
      icen = match_center(ao_basis%gtos_small(ishell)%coord, center_coord, center_match_tol)
      if (icen == 0) stop 'small-component shell has no matching large-component center'
      shell_to_center_small(ishell) = icen
   end do
end if
```

This is the minimum possible amount of code because the consistency check and the
construction of `shell_to_center_small` are the same loop -- and that map is needed
anyway for the small-component IAO blocks. Checking costs nothing beyond what building
the map already requires, and it is strictly stronger than comparing center counts.

### 6.3 Undeclared basis centers -- clarification and proposal

This is the item that needed clarifying. The phrase "a basis center that is not a
declared ghost" was doing too much work in one clause.

Every AO center that becomes an atomic fragment needs a `bascharge`, because
`import_fragments.F90:117` uses it to choose the atomic library file. That value lives in
the nucleus list, so each center must be matched back to a nucleus-list entry. There are
three outcomes:

| Case | Nucleus-list entry | `bascharge` available? |
|---|---|---|
| (a) Ordinary atom | present, nuc > 0, bas > 0 | yes |
| (b) **Declared** ghost | present, nuc approx 0, bas > 0 | yes |
| (c) **Undeclared** basis center | **absent entirely** | no |

The word *declared* refers to occupying a row in `/input/molecule/geometry`. A ghost is
declared: it is listed in the molecule block, simply carrying zero nuclear charge. That
is why case (b) already works today -- the element is known, so the library file can be
found. Case (c) is a center that exists **only** in the AO basis, with no counterpart in
the molecule block at all.

Case (c) can arise from a floating or bond-centered basis function placed away from any
listed center, from an interface that writes AO centers but a shorter molecule list, or
from exactly the dataset inconsistency that section 2.3 warns about.

It needs a decision because such a center has no determinable element, therefore no
reference atomic basis, therefore cannot become an atomic IAO fragment -- yet its basis
functions still span part of the AO space. Silently dropping it would leave those
functions belonging to no fragment, making the IAO projection incomplete and changing
results without warning.

**Decided: hard error, naming the coordinate.** The user can then either add the center
to the molecule list as a declared ghost, which makes it work immediately, or switch to
molecular fragments. Silent skipping is not an option, because the failure mode is a
quietly wrong answer rather than a crash.

Implemented as a claim-and-verify pass: each nucleus that matches a center claims it, and
any center left unclaimed at the end triggers the error with its coordinate printed. The
same pass catches the converse inconsistency -- a molecule-list entry declaring no basis
(`bascharge <= charge_match_tol`) that nonetheless sits on a center -- and a tolerance so
loose that two entries claim the same center.


---

## 7. Implementation notes (step 2a)

Step 2a is implemented. Two claims made earlier in this note turned out to be wrong and
were corrected in the code rather than worked around; both are recorded above at the
point where they were originally asserted.

### 7.1 The ordering premise was wrong

Section 6.1 originally justified numbering centers by first appearance on the grounds that
it "preserves the current ordering for conventional systems". It does not.
`add_atom` in `read_xyzfile.F90` does not append: it inserts each new atom immediately
after the last atom **of the same name**, so the nucleus list is grouped by element. For
an input order H, C, H the list is built as H, H, C. The AO shell ordering follows the
checkpoint and does not reproduce that grouping, so numbering centers by first appearance
in the shell list would in general have permuted the fragments.

That permutation would not have been harmless: besides relabelling the partial-charge
output, `frag_bias(1) = 0.1` biases whichever fragment is first, so the tie-breaking for
non-polar bonds would have changed.

The implementation therefore splits the two concerns:

- **Ordering comes from the nucleus list.** The construction loop walks
  `molecule%firstatom` exactly as before, so the fragment numbering is unchanged.
- **Membership comes from the AO basis.** An entry with no matching center (a point
  nucleus) is skipped; every center must be claimed by exactly one entry.

This delivers everything section 3.1 promised while being genuinely behaviour-neutral.

### 7.2 The AO basis must be read before the fragments are defined

`rose_define_fragments` did not have access to the supermolecular AO basis: `rose_main`
called it at line 63 and only read `ao_basis_1` at line 64. The two calls are now swapped
and the basis is passed in.

This is safe in both directions: `rose_read_mos` was checked and touches no fragment state
at all (it depends only on `molecule%name`, `file_for_MOs`, and the `restricted` /
`spatial_orbitals` flags), and `rose_define_fragments` gains a dependency on exactly the
object the swap makes available.

**Untested risk:** the reorder applies to every interface, including ADF/AMS. The ADF/AMS
fragment path itself is deliberately untouched -- `atom_fragments` is false there, so those
runs take the unchanged molecular branch -- but the reordering of the two calls has not
been exercised against an ADF/AMS checkpoint.

### 7.3 Output may differ in the last digits

The atomic library basis is now translated onto the AO center rather than onto the nucleus.
For a conventional molecule these are the same point, but they need not be bitwise equal:
the geometry is stored in Angstrom and divided by `xtang` on input, while the shell centers
are read directly. Any difference is at the level of coordinate round-off and feeds into
the overlap integrals, so regression comparisons should expect agreement to within
round-off rather than bit-for-bit identity. This is the intended direction -- the AO center
is the coordinate the supermolecular orbitals were actually built on.

### 7.4 What is computed but not yet used

`shell_to_center` and `shell_to_center_small` are built and immediately discarded. They are
the input to deriving the IAO block structure rather than assuming it, which is step 2b
work; using them now would change behaviour. The small-component *check* is not wasted,
however, since it happens inside `collapse_ao_centers` while the map is being built.

### 7.5 Verification

The development machine has no Fortran compiler reachable from the automation
environment (the shell available there is a Linux VM with only the project folder
mounted; the Homebrew toolchain lives on the macOS side). Verification was therefore
done in a scratch Linux container with gfortran 13.3.

**Full build of the tree: clean.** CMake configured (BLAS/LAPACK found, HDF5 not
present so `ENABLE_HDF5=OFF`) and the whole project compiled and linked to `rose.x`
with **zero errors**. 40 Fortran objects plus the static library and the executable.

**No new warnings.** Under `-Wall -Wextra`, `rose_ao.F90` -- which is entirely new code
-- is completely clean. The only warnings in the other touched files are pre-existing and
unrelated: an unused label 563 and an unused private `ios` in `read_xyzfile.F90`, and an
unused `ierr` at line 174 of `import_fragments.F90`, which is inside
`rose_read_fragments`, not the subroutine changed here.

**Unit tests of the new helpers, all passing:**

- `collapse_ao_centers` on six shells over three centers in a deliberately non-grouped
  order (A B A C B A) gives `n_centers = 3` and `shell_to_center = 1 2 1 3 2 1`.
- Tolerance boundaries: 0.015 bohr matches, 0.030 bohr does not, and the worst-case
  two-decimal-Angstrom diagonal of 0.0163 bohr matches -- the case that motivated raising
  the tolerance from 0.01 to 0.02.
- The small-component path, with small-component shells given in a different order from
  the large component and offset by round-off, builds `shell_to_center_small` correctly.

**What this build did not cover.** With `ENABLE_HDF5=OFF` the low-level bindings
`hdf5_util/mh5.F90` were replaced by `mh5_nohdf5.F90`; note that `hdf5_interface.F90`
itself -- the file holding `read_from_hdf5`, which populates the nucleus list and the AO
basis -- *was* compiled. Also not built in this configuration: `rose_ao_ADF.F90`,
`correlating_basisfunctions_dirac.F90`, and the real AMS bindings. None of these is
touched by this change, but it does mean the ADF/AMS build variant is unverified, which
compounds the untested reorder noted in 7.2.

**Still outstanding:** an end-to-end run. The natural next check is the existing test
cases (`dirac_atomfrags`, `psi4_atomfrags`, `pyscf_atomfrags` in particular), confirming
the output matches to round-off rather than bit-for-bit, for the reason given in 7.3.


---

## 8. Implementation notes (step 2b -- labelled xyz input)

Step 2b is implemented. Step 2a passed the full regression suite on all four
interfaces (DIRAC, PSI4, PySCF, ADF), so the center enumeration underneath this is
confirmed behaviour-neutral in practice, not just by argument.

### 8.1 What was added

**`read_xyz_labels` in `read_xyzfile.F90`.** Reads an extended xyz file whose atom
lines carry a fifth field:

```
22
biphenyl, split into two phenyl groups at the 3-7 ring-ring bond
C       -0.44664000       1.44424000       0.93720000   phenyl1
```

Coordinates are converted Angstrom -> bohr on read, matching `read_xyz_atoms`, so the
later comparison happens in internal units. A missing label is a hard error naming the
offending line, rather than a silent fallback.

**`.FRAGMENT_XYZ` keyword** (`fragment_xyz_file` in `rose_input_data`). It does not
require `.NFRAGMENTS` to precede it, which is the whole point: the file is matched
against the AO centers after those are known, so it works for atomic fragments, where
the fragment count is not available at input-reading time. The earlier `.GROUPS`
keyword could not do this, and has since been removed -- see section 10.

**Label-to-group mapping in `rose_define_fragments`.** Each line is matched to an AO
center with `match_center`, mapped through `center_owner` to a fragment index, and
assigned a group number by first appearance of its label. `n_groups` and
`frag_to_group` are then set directly, so `rose_init_groups` validates them and leaves
them alone rather than installing the identity.

The section 4 rule is honoured exactly: the coordinates in the labelled file are
consumed by the match and then discarded. Nothing is ever placed using them.

### 8.2 Errors it raises

- a labelled line matching no AO center (geometries differ, or the tolerance is too
  tight), printing the offending coordinate;
- two lines claiming the same center;
- a fragment left without a label, printing the fragment name and coordinate;
- `.FRAGMENT_XYZ` with molecular fragments, which is not implemented.

The grouping is also printed to `luout` under a `Recanonization Groups` banner, so it
is visible in the output and can be filtered by a regression test.

### 8.3 The test case

Placed under `test/`, modelled on the existing runtest cases. It was first built on
biphenyl; that molecule made the run too slow for a regression test and was replaced
by fluoroethane, which is the same exercise at 8 atoms instead of 22 (see 8.6).

- `test/xyz_files/fluoroethane.xyz` -- the experimental geometry from NIST CCCBDB
  (Hellwege 1976, CAS 353-36-6).
- `test/xyz_files/fluoroethane_groups.xyz` -- the same geometry labelled
  `fluoromethyl` / `methyl`, 4 atoms each. The split was derived by connectivity
  rather than by eye: the C-C bond is the only bond whose removal disconnects the
  graph, atoms 1-2 at 1.5050 A.
- `python_scripts/groupfrag_psi4.py` -- launcher, a patched copy of
  `atomfrag_psi4.py` that writes `.FRAGMENT_XYZ` and stages the labelled file into the
  run directory.
- `test/psi4_groupfrags/test` -- runtest script.
- `test/runtest_config.py` -- additive selector: a `mol` field of the form
  `geometry:labels` picks `groupfrag_*.py`. Existing tests contain no colon and are
  unaffected (verified by simulating the selector on all three shapes).
- `CMakeLists.txt` -- globs `groupfrag_*.py`, copies it to the build directory, and
  registers `rose_test(psi4_groupfrags ...)`.

**Which group gets the severed bond.** The tie is broken by `frag_bias`, which
defaults to 0.1 on the *first* fragment. The nucleus list starts with atom 1, the
fluorinated carbon, so the C-C bond lands in `fluoromethyl`. That asymmetry is
deliberate and is what the example exists to exercise.

### 8.4 Verification

- Full build from scratch in the gfortran container: clean, `rose.x` links, no new
  warnings (the three pre-existing ones from 7.5 are unchanged).
- `read_xyz_labels` unit-tested on the real labelled file: correct line count,
  correct Angstrom-to-bohr conversion, correct first and last labels.
- The label-to-group assignment was exercised against a synthetic AO basis whose
  centers were supplied in **reversed** order relative to the labelled file, giving 2
  groups of 11 with the right names -- confirming the join is by coordinate and not by
  line order.
- `groupfrag_psi4.py`, `runtest_config.py` and the test script all byte-compile, and
  the selector was checked to route all three `mol` shapes correctly.

**Not verified, and needing a run on real hardware:** there is no reference output in
`test/psi4_groupfrags/result/` yet, because generating it requires PSI4. That has to be
produced by running the case once and checking the result is sensible before it is
frozen as the regression reference.

### 8.5 First run: fragment multiplicities (fixed)

The first attempt crashed before ROSE was reached, inside PSI4's molecule parser:

```
ValidationError: Inconsistent or unspecified chg/mult:
   sys chg: None, frag chg: [0.0], sys mult: None, frag mult: [1]
```

The fault was in the test input, not in the labelled-xyz code. The second bracket of
the runtest input is a list of multiplicities: the first entry is the molecule, and the
rest are one per **distinct element**, in order of first appearance in the xyz file.
The test had been written with `[1,1,1]`, which asks for a hydrogen atom -- one electron
-- in a singlet. That is impossible, so qcelemental rejected it. Carbon at multiplicity
1 passed validation, since six electrons can pair, but was equally unintended.

The convention was already established by `pyscf_atomfrags`, which runs ethylene as
`[0,0,0] [1,3,2]` and HCN as `[0,0,0,0] [1,4,3,2]` -- atoms take their own ground-state
multiplicity: H is 2S, C is 3P, N is 4S. Biphenyl is 12 carbons followed by 10
hydrogens, the same element order as ethylene, so it takes the same `[0,0,0] [1,3,2]`.
The test now carries that string with a comment explaining the ordering rule.

Two guards were added to `groupfrag_psi4.py` so this fails legibly next time:

- the number of charges and multiplicities supplied is checked against the number of
  distinct elements before any SCF is attempted;
- each fragment's electron count is checked against its requested multiplicity for
  parity, naming the element and the counts. This rejects only *impossible*
  combinations, not merely non-ground-state ones, so it does not second-guess a
  deliberate choice of electronic state.

`atomfrag_psi4.py` has the same gap but was left untouched, since it is covered by
passing regression tests and the change would be gratuitous there.


---

## 9. The additional-virtuals segfault (step 1 bug, found by the biphenyl run)

The first successful biphenyl run reached the IAO stage cleanly and then segfaulted where
the additional virtuals are merged. The cause was a latent bug in **step 1**, not in the
labelled-xyz code, and it could only ever appear once `n_groups /= n_fragments` -- which
the biphenyl case is the first to do. That is why the whole regression suite passed.

### 9.1 What was wrong

Step 1 made `nocc_frag`, `nvir_frag` and `nlo_frag` **group-sized**, and I recorded at
the time that "all downstream routines are already size-agnostic". That was wrong for
`rose_merge_virtuals`, which was hard-coded to `n_fragments` in five places.

The immediate crash: `nvir_frag` arrives with size `n_groups` (2), `nadd_frag` with size
`n_fragments` (22, one block per atomic fragment, since the hard virtuals come from each
fragment basis). The routine then did

```fortran
allocate(nvir_all(n_fragments))
nvir_all = nvir_frag + nadd_frag     ! (2) + (22) -- not conformable
do i = 1, n_fragments
   ... Lvir(:,off2+1:off2+nvir_frag(i),:)   ! nvir_frag(i) read past the end
```

so the merge loop indexed a 2-element array up to 22, got garbage extents, and walked off
the end of `Lvirt`.

### 9.2 Why aggregating the counts is not enough

`Lvir` is ordered by **group** (that is what `sort_on_fragments` produces), while `Ladd`
is built one block per **fragment**. Summing `nadd_frag` into a group-sized array fixes
conformability but still assumes each group's additional virtuals are contiguous in
`Ladd`, and they are not.

Biphenyl shows this plainly. The nucleus list is grouped by element, so the fragment
order is C1..C12 then H13..H22, while the ring membership is

```
fragments  1- 6  -> phenyl1      fragments 13-17 -> phenyl1
fragments  7-12  -> phenyl2      fragments 18-22 -> phenyl2
```

Each group's fragments come in two disjoint runs. The additional virtuals therefore have
to be *gathered* by group, not merely counted.

### 9.3 The fix

`rose_merge_virtuals` now records the offset of every fragment block inside `Ladd`,
aggregates `nadd_frag` into `nadd_group`, and walks groups: for each group it copies that
group's valence virtuals (contiguous in `Lvir`) and then gathers the additional-virtual
block of every fragment belonging to it, wherever those sit in `Ladd`.

`copy_original_occ` and `copy_original_vir` were allocating their block-count arrays with
`n_fragments` as well; both now use `n_groups`, since what they describe is a
recanonization block, not a fragment.

The heading of the per-block energy listing follows suit, but conditionally: it reads
`Group` only when `n_groups /= n_fragments`, and stays `Fragment` otherwise. An
unconditional rename was the first attempt and it broke both ase tests, whose
`MO_ENERGIES_REGEX` matches that header literally against stored reference output (see
section 12).

`nadd_frag` itself stays fragment-sized everywhere it is produced -- in
`generate_fragvirs`, `generate_corrfuns/build_fragments`, and the "no additional
virtuals" branch of `rose_main`. That is correct: the reference virtuals genuinely do
come from the individual fragment basis sets, and the group view is needed only at the
merge.

### 9.4 Verification

The merge loop was extracted verbatim into a standalone program and run on two cases:

- **Interleaved, biphenyl-like** (6 fragments, groups `1 1 2 2 1 2`): column order comes
  out `V1 V2 | A1 A2 A3 A6 A7 | V3 V4 V5 | A4 A5 A8`, i.e. group 1 correctly collects the
  additional virtuals of fragments 1, 2 and 5 from their scattered positions, and
  `nvir_all = [7,6]` sums to the 13 columns available.
- **Identity** (`frag_to_group(i) = i`): reproduces the old per-fragment interleaving
  exactly, `V1 A1 V2 A2 ...`, confirming the change is behaviour-neutral for every
  existing test.

In both cases every source column is placed exactly once and the running offset ends on
the last column. The full tree also rebuilds clean, with no new warnings (the two
`-Wconversion` notices in `generate_ibo.F90` are pre-existing, both at
`nocc_active = molecule%znumber - charge`).

### 9.5 Audit of the remaining consumers

Checked for the same class of mismatch, and clean:

- `recanonize` takes `size(nlo_frag)` and is genuinely size-agnostic.
- `rose_write_mos` uses `size(Locc,2)` / `size(Lvir,2)` and passes the block counts
  straight through.
- `write_data_rose` writes `size(ribo_occ_per_frag)` entries, so it records group counts
  now, which is the right description of the recanonized blocks.
- `write_ROSE_results` is AMS-only and cannot be reached with groups, because
  `.FRAGMENT_XYZ` is rejected for the ADF/AMS path where `atom_fragments` is false.
- `build_fragments` in `project_corrfuns.F90` allocates `nadd_frag(n_fragments)`, which
  is correct and is handled by the aggregation at the merge.


---

## 10. Removal of the `.GROUPS` keyword

`.GROUPS` -- the original grouping input, a list of one integer per fragment on a
single line -- has been removed. It was never used in anger, and `.FRAGMENT_XYZ`
supersedes it on both counts that matter.

It was **unfriendly**: a bare row of integers whose meaning depended entirely on
knowing the fragment ordering, which, as section 7.1 shows, is not the order anyone
would guess. And it was **structurally limited**: it read its integers at
input-parsing time, so it required `.NFRAGMENTS` to have been seen already, which
confined it to molecular fragments and made it unusable for exactly the atomic-fragment
case this work was for. That limitation is what prompted the labelled-xyz design in the
first place.

What was removed:

- the keyword handler in `input_reader_sections.F90`;
- the "both keywords given" guard in `rose_define_fragments`, now unreachable because
  nothing else allocates `frag_to_group` before that point.

What stays: the whole group layer -- `n_groups`, `frag_to_group`, `rose_init_groups`
and every consumer. Only the input route changed.

One consequence worth stating plainly. `.GROUPS` was the only grouping input that
worked for *molecular* fragments, so with it gone, grouping is available for atomic
fragments only. Nothing is lost in practice, since the keyword was never used, but the
molecular case is now genuinely unimplemented rather than merely untested. The error
raised for `.FRAGMENT_XYZ` with molecular fragments says so, and points at extending
that keyword as the intended route rather than at a keyword that no longer exists. The
`frag_to_group` layer is indifferent to which route fills it, so that extension remains
a matter of input handling, not of structure.


---

## 11. Replacing biphenyl with fluoroethane

Biphenyl was a poor regression case: 22 atoms and 210 basis functions made the run
too slow to sit in a suite that is meant to be run often. It was replaced by
fluoroethane, 8 atoms, cut at the C-C bond into `fluoromethyl` (CH2F) and `methyl`
(CH3).

**The geometry is the experimental one** from NIST CCCBDB (Hellwege 1976, CAS
353-36-6), not an optimised or invented structure. Its Cartesian coordinates were
checked against the internal coordinates quoted on the same page and reproduce them:
C-C 1.5050 (1.505), C-F 1.3980 (1.398), C-H 1.0950 and 1.0910 (1.095, 1.091), C-C-F
109.70 (109.7), H-C-F 106.04 (106.1), H-C-C 112.90 on the fluorinated carbon (112.9),
and all three H-C-H values 108.81 / 108.90 / 108.70 (108.8 / 108.9 / 108.7).

**The coverage is preserved**, which is the part that mattered. Section 9 turns on
each group's fragments not being contiguous in fragment order, so that the additional
virtuals have to be gathered rather than merely counted. Because `add_atom` groups the
nucleus list by element, fluoroethane gives

```
frag_to_group = [1, 2, 1, 1, 1, 2, 2, 2]
   group 1 (fluoromethyl) owns fragments 1, 3, 4, 5
   group 2 (methyl)       owns fragments 2, 6, 7, 8
```

Both groups are non-contiguous, so the gather path is exercised exactly as biphenyl
exercised it. Fragment 1 is the fluorinated carbon, so `fluoromethyl` takes the
severed C-C bond through `frag_bias`, keeping the two groups inequivalent.

**Input changes.** Three distinct elements now, in first-appearance order C, F, H, so
the runtest input becomes `[0,0,0,0] [1,3,2,2]` -- ground-state multiplicities C 3P,
F 2P, H 2S. The filter line counts were reduced with the shorter output (the
recanonization-group block from 30 to 16 lines, the partial-charge block from 26 to
14). The output prefix follows the molecule name, so the reference file to freeze is
now `psi4_fluoroethane.stdout`.

Sections 8.5 and 9 still describe the biphenyl runs, since that is where those
problems were found; the history is left as it happened.


---

## 12. The ase test failures, and why the heading is conditional

Renaming the per-block energy heading from `Fragment` to `Group` was a one-line change
for accuracy, and it broke both ase test targets.

`ase_rose/tests/test_*.py` define

```python
MO_ENERGIES_REGEX = 'Fragment          MO              Energy.*\n...'
```

and apply it to **both** sides of the comparison: the stored reference under
`ase_rose/tests/data/` and the output of the fresh run. The references still contain
`Fragment          MO              Energy`, so the expected side matched while the
actual side no longer did, and `extract_number` raised
`ValueError("No pattern found in text.")`. Only `test_virtual_energies` failed in each
target; `test_partial_charges` and `test_avas` use different patterns.

Auditing the whole diff for output-visible text showed that this single header was the
*only* change any pre-existing test could see. Everything else added is either an error
path, reached only on failure, or inside the `.FRAGMENT_XYZ` block, reached only when
grouping is requested.

The heading is therefore now chosen at run time: `Fragment` when
`n_groups == n_fragments`, `Group` otherwise. The identity-case line was checked
byte-for-byte against the stored references and is identical, so no reference output and
no regex needed changing.

This is the same invariant the rest of the design holds to -- a run that requests no
grouping should be indistinguishable from one before this work -- and it is worth
extending to output, not just to results. The alternative, updating four reference files
and four test files to chase a cosmetic rename, would have put churn into unrelated tests
for no gain.
