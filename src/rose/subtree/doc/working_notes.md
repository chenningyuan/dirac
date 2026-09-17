# ROSE working notes

Practical knowledge accumulated while implementing the fragment/group separation.
The *design* reasoning lives in `atomic_fragment_centers.md`; this file holds the
things that are easy to rediscover the hard way -- conventions, traps, the build
and test workflow, and what is still open.

---

## 1. Conventions that will bite you

These are all verified against the code or against real data, not inferred.

### 1.1 The nucleus list is ordered by element, not by input order

`add_atom` in `read_xyzfile.F90` does not append. It inserts each new atom
immediately after the last atom **of the same name**, so an xyz written H, C, H
produces a list ordered H, H, C.

Anything that assumes the nucleus list follows the file is wrong. This is why the
fragment ordering must be driven by walking `molecule%firstatom` rather than by the
AO shell order, which follows the checkpoint instead (see
`atomic_fragment_centers.md` section 7.1 -- getting this wrong silently permutes the
fragments and moves `frag_bias(1)`).

### 1.2 AO centers are the real primitive, not nuclei

Fragment matching is done by analytic Gaussian overlap on shell coordinates alone
(`get_overlap.F90`); `gtos%atom_number` is written only by the fchk interface and
never read back for mapping. A ghost center carries a basis and no charge; a point
nucleus carries a charge and no basis. Full argument in
`atomic_fragment_centers.md` sections 2 and 3.

### 1.3 `orbmom` in the ROSE h5 is **1-based**

`l=1` is s, `l=2` is p. Reading it as 0-based gave 522 AOs instead of 210 for the
biphenyl/STO-3G case this was first checked on. Always reconcile the implied AO count against `aobasis/n_ao` before
trusting an interpretation.

### 1.4 `contractions` holds the primitive normalisation constant

Across all 138 shells of the biphenyl/STO-3G case this was measured on,
`c/N = 1.00000000` exactly,
with

```
N = (2a/pi)^(3/4) * (4a)^(l/2) / sqrt((2l-1)!!)
```

So ROSE stores unnormalised primitives with the normalisation folded into the
coefficient. Molden expects coefficients **for normalised primitives**, so a writer
must divide by `N` -- which for a decontracted basis makes every coefficient 1.0.
Skip this and the orbitals come out scaled per shell: plausible pictures, wrong
densities.

### 1.5 Group-sized versus fragment-sized arrays

Since the step-1 separation these two lengths differ whenever `.FRAGMENT_XYZ` is
used, and mixing them is a segfault rather than an error.

| length `n_groups` | length `n_fragments` |
|---|---|
| `nocc_frag`, `nvir_frag`, `nlo_frag`, `nvir_all`, `nadd_group` | `nmo_frag`, `nadd_frag`, `frag_to_group`, `frag_bias`, `frag_ncor`, `frag_nval`, `fragments(:)` |

`nadd_frag` is deliberately fragment-sized: the reference virtuals genuinely come
from the individual fragment basis sets. It is aggregated to group level only at the
merge in `rose_merge_virtuals`.

Watch for two separate mistakes here. Aggregating the *counts* is not enough -- the
`Ladd` columns also have to be **gathered** by group, because fragments belonging to
one group need not be adjacent in fragment order. In the regression case the nucleus list is
element-grouped, so phenyl1 owns fragments 1-6 and 13-17: two disjoint runs. See
`atomic_fragment_centers.md` section 9.

### 1.6 Tolerances

`center_match_tol = 2.0D-2` bohr and `charge_match_tol = 1.0D-2` e, declared in
`read_xyzfile.F90` next to `ghost_threshold`, which stays at `1.0D-10`. The center
tolerance is in bohr deliberately, so no conversion happens at any call site. 0.02
bohr is about 0.0106 Angstrom, chosen so that a coordinate hand-typed with two
decimals in Angstrom still matches under a Euclidean test (worst case
sqrt(3) x 0.0094 = 0.0163 bohr).

---

## 2. Test harness

### 2.1 Launcher selection

`test/runtest_config.py` picks the launcher from the shape of the `mol` field:

| `mol` looks like | launcher |
|---|---|
| `geometry:labels` | `groupfrag_<program>.py` -- atomic IAOs grouped by a labelled xyz |
| `[a,b,c]` | `molfrag_<program>.py` -- molecular fragments |
| `name` | `atomfrag_<program>.py` -- atomic fragments |

The colon branch was added for `psi4_groupfrags`; existing tests contain no colon and
are unaffected.

### 2.2 Charges and multiplicities

The second bracket is `[charges] [multiplicities]`. The **first** entry is the
molecule; the rest are one per **distinct element**, in order of first appearance in
the xyz file -- *not* one per atom, and not one per fragment.

Atomic fragments take their own ground-state multiplicity: H = 2, C = 3, N = 4,
O = 3. Established by `pyscf_atomfrags`, which runs ethylene as `[0,0,0] [1,3,2]`
and HCN as `[0,0,0,0] [1,4,3,2]` (HCN is ordered N, C, H).

Getting this wrong fails deep inside qcelemental with an opaque
`Inconsistent or unspecified chg/mult` error. `groupfrag_psi4.py` now guards against
it; `atomfrag_psi4.py` does not.

### 2.3 Reference outputs

`runtest` diffs filtered stdout against `test/<case>/result/<program>_<mol>.stdout`.

**Every `get_filter` needs a tolerance** -- `rel_tolerance` or `abs_tolerance`. Omitting
it does not fall back to a default: the tolerance field of runtest's filter namedtuple
stays unset and the comparison fails deep inside the library with

```
TypeError: '<=' not supported between instances of 'float' and 'collections._tuplegetter'
```

which names neither the filter nor the test. For a block of structural integers (indices,
counts) use `abs_tolerance=1.0e-8`, which is an exact match in practice.
**`test/psi4_groupfrags/result/` is still empty** -- the reference has never been
frozen. Run the case, check the result is sensible, and only then save it. The file
it wants is `psi4_fluoroethane.stdout`, the prefix following the molecule name.

---

## 3. Building and verifying without a local compiler

The automation environment attached to this repo is a **Linux aarch64 VM** with only
the connected folders mounted. It has `gcc` and `make` but **no Fortran compiler**,
no `cmake`, no `h5py`, and cannot install any (uid 1001, `no new privileges` blocks
sudo, and the distro repos are blocked by the egress proxy). The Homebrew toolchain
lives on the macOS side and is not reachable from it -- `/opt/homebrew` and even
`/Users` do not exist there.

What works instead: stage the sources to a cloud container and build there.

```
# in the container, once
apt-get install -y gfortran cmake liblapack-dev libblas-dev

# stage src/*.F90 + src/CMakeLists.txt + CMakeLists.txt
#   + src/hdf5_util/* + src/interest/*      (49 files, two staging calls)

# the top-level CMakeLists installs python packages from ase_rose/ and
# python_scripts/, which are not staged -- strip that block from the COPY,
# between "# Now, install and/or upgrade runtest and genibo" and
# "#### final of Updates from the original ROSE CMakeLists.txt"

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build -j4
```

HDF5 is optional (`ENABLE_HDF5=OFF` if not found); BLAS and LAPACK are required.
This builds and links `rose.x` cleanly.

**What that build does not cover:** with HDF5 off, `hdf5_util/mh5.F90` is replaced by
`mh5_nohdf5.F90` -- though `hdf5_interface.F90` itself, which holds `read_from_hdf5`,
*is* compiled. Also unbuilt: `rose_ao_ADF.F90`, `correlating_basisfunctions_dirac.F90`
and the real AMS bindings. So the ADF/AMS variant is never exercised this way.

**Pre-existing warnings**, so they are not alarming when they appear: unused label 563
and unused private `ios` in `read_xyzfile.F90`; unused `ierr` at line 174 of
`import_fragments.F90` (inside `rose_read_fragments`); two `-Wconversion` notices in
`generate_ibo.F90` at `nocc_active = molecule%znumber - charge`.

A useful trick for checking one subroutine without the whole tree: extract it into a
harness module against the real `read_xyzfile` and `rose_ao` with a stub
`rose_input_data`. That is enough to typecheck interfaces and catch structural errors.

---

## 4. Utilities

### `test/npy2molden.py`

Writes a Molden file from a psi4 wavefunction archive (`ribo.npy`). Verified working
on a real run.

Three things worth knowing:

- The archive carries the basis **name** only, never the primitives, so psi4 has to
  rebuild the basis. A numpy-only converter from `.npy` is impossible.
- A wavefunction restored by `Wavefunction.from_file` is a bare
  `psi4.core.Wavefunction` with **no `occupation_a()`**. `wfn.write_molden` calls it
  unconditionally and raises `AttributeError`; `psi4.molden()` only warns and writes
  zeros. The script tries the former and falls back to the latter, then fills the
  occupations in itself from the archive's `nalpha`/`nbeta`.
- ROSE's psi4 runs use an **ROHF** reference, so psi4 emits an alpha block and an
  identical beta block -- twice as many orbitals as the calculation has. The script
  drops the beta block when the archive says both spins share orbitals.

`psi4.molden` is deprecated and expected to stop working around psi4 1.5, at which
point this script needs occupation vectors attached to the wavefunction instead.

### Wanted: an h5-to-molden converter

`ribo.h5` already holds everything required -- geometry, `nuc_charge`, per-shell
`center`/`exponents`/`orbmom`, `mobasis/orbitals` as `(nmo, nao)` row-major, and
`eigenvalues`. Occupations follow from `sum(nuc_charge) - charge`, which sidesteps the
`occupation_a` problem entirely. Mind sections 1.3 and 1.4 above. Unlike the psi4
route this would also serve DIRAC and ADF.

---

## 5. Open work

1. **Freeze the `psi4_groupfrags` reference** (section 2.3).
2. **Installation testing inside AMS and DIRAC** -- not done since the group work
   landed.
3. **Mixing atomic and molecular fragments** in one calculation. The `frag_to_group`
   layer was deliberately built to make this possible without further structural
   change.
4. **Grouping for molecular fragments.** Currently rejected with a clear error. The
   old `.GROUPS` keyword used to cover this case and was removed (see
   `atomic_fragment_centers.md` section 10), so it is now unimplemented rather than
   just untested. The route is to extend `.FRAGMENT_XYZ` to that case; the
   `frag_to_group` layer needs no change.
5. **`ghost_threshold` unification** -- deferred decision, see
   `atomic_fragment_centers.md` section 3.3.
6. **h5-to-molden converter** (section 4).

---

## 6. State of the tree

Implemented and regression-tested across DIRAC, PSI4, PySCF and ADF: the step-1
fragment/group separation and the step-2a center-derived enumeration.

Implemented and exercised once end to end: step 2b, the labelled-xyz input, plus the
`rose_merge_virtuals` fix that the run exposed. That fix has not been through the full
regression suite since it changed -- it is behaviour-neutral by construction for
`n_groups == n_fragments`, and was unit-tested against the identity case, but a suite
run would confirm it.

The end-to-end run was done on biphenyl. The committed regression case is
**fluoroethane**, which is much faster and preserves the coverage that matters: both
groups are non-contiguous in fragment order (`frag_to_group = [1,2,1,1,1,2,2,2]`), so
the gather path in `rose_merge_virtuals` is still exercised. Fluoroethane itself has
not yet been run end to end. See `atomic_fragment_centers.md` section 11.
