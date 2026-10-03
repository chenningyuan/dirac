#!/usr/bin/env amspython
"""
rose_kf_bridge -- move data between an ADF binary (KF) file and the HDF5
checkpoint format that ROSE reads.

This is a developer tool, not part of the test suite. It is run by hand under
amspython, and the files it produces are checked in as fixtures. If it is ever
called from a test, the AMS dependency that step 3 of the restructuring exists
to remove has come back through the back door.

    amspython rose_kf_bridge.py inspect adf.rkf
    amspython rose_kf_bridge.py forward adf.rkf -o CHECKPOINT.h5 [--t15 x.t15]

inspect   prints the sections and variables the conversion depends on, with
          shapes and leading values. Use it on a new case before converting.
forward   writes CHECKPOINT.h5 from adf.rkf. With --t15 it also stores ADF's
          own AO overlap in /validation/, which ROSE then checks its Slater
          overlaps against; see the gate in hdf5_interface.F90.

    amspython rose_kf_bridge.py backward ribo.h5 -o Molecule_B1_MO1_ribo.rkf \
              --template Molecule_B1_MO1.rkf
    amspython rose_kf_bridge.py fixture <plams_workdir> -o <dir> [--input ...]

backward  writes the Rose Data section ADF's TDDFT consumes. It takes over from
          write_mos_to_ADF and write_data_rose, removed with the file-based ADF
          route and recoverable from the git history. It needs nfrag, the per-group
          counts and the RIBO Fock matrix, which replace_mobasis now stores
          beside the MO dimensions as /mobasis/ribo_*; an h5 file written
          before that cannot be converted back and says so.
fixture   converts a whole plams_workdir in one go: the supermolecule to
          CHECKPOINT.h5 and each fragment to NNN.h5, named by atomic number,
          because that is how rose_read_fragments builds the name.

WHAT THE ADF SLATER BASIS LOOKS LIKE
------------------------------------
ADF stores Cartesian-power Slater functions,

    chi = N x^kx y^ky z^kz r^kr exp(-alpha r),

with kx+ky+kz = l and kr = n-1-l, so a d shell has six components and not five.
In ROSE terms that is basis_angular = 1 (Cartesian) and orb_momentum = l+1.
The component order ADF uses, xx xy xz yy yz zz, is the order ROSE already
generates in get_spherical_projection, so no reordering is needed inside a
shell. ADF's Slater functions are uncontracted, one primitive per shell.

The Basis section is held per atom TYPE at two levels:

    nbaspt(1:ntypes+1)  points into the per-shell arrays
                        nqbas (n), lqbas (l), alfbas (zeta), basnrm
    nbptr (1:ntypes+1)  points into the per-component arrays
                        kx, ky, kz, kr, alf, bnorm

so type t owns shells nbaspt(t) .. nbaspt(t+1)-1 and components
nbptr(t) .. nbptr(t+1)-1. The types are then replicated over the atoms that
carry them to give naos functions in total, atom by atom.

Normalisation is per Cartesian component and is reproduced exactly by
cart_angular_norm of slater_overlap, so basnrm and bnorm are read only to be
checked, never to be transported.

AO ORDERING
-----------
A%Eigen-Bas_A is NOT in the Basis order: the two p shells of an atom are
interleaved so that the px functions run together, then the py, then the pz.
A%npart holds the permutation. Its direction was established from ROSE's own
former ADF reader, which formed

    S11_AO(i,j) = smat_ao_full(order_ao(i),order_ao(j))

with S11_AO in MO-coefficient order and smat in Basis order. So npart(i) is the
Basis-order index of the i-th function in Eigen-Bas order. That reader has since
been deleted along with the rest of the file-based ADF route, so this comment is
the record of it. The overlap gate re-checks the conclusion on every converted
file, which is why it is worth keeping the reference overlap in the fixtures.

Since that permutation cuts across shell boundaries it cannot be expressed as a
reordering of whole shells, and /aobasis/ is shell structured. The conversion
therefore writes the basis in the Basis order and permutes the COEFFICIENTS
into it, c_basis[npart[i]-1] = c_adf[i].
"""

import argparse
import sys

import numpy as np

# UNITS. The two coordinate arrays in the h5 file are in DIFFERENT units, and
# this is not a mistake of the format but what ROSE's reader expects:
#
#   /input/molecule/geometry   ANGSTROM.  read_hdf5 divides it by xtang, the
#                              same way an xyz file is treated, so writing bohr
#                              here puts every nucleus 1.89 times too far out.
#   /aobasis/center            BOHR.      read_ao_basis takes it unconverted.
#
# Nothing fails loudly when they disagree: the AO overlap is built from the
# centres alone and comes out right, while collapse_ao_centers quietly fails to
# match any AO to any nucleus, because it compares the two against a 0.02 bohr
# tolerance. KF holds both in bohr, so the geometry is converted on the way out
# and the centres are not.
ANGSTROM_PER_BOHR = 0.5291772083        # xtang of read_xyzfile.F90, to the digit


def kffile_class():
    """Import KFFile on use, so that the parts of this module that are plain
    arithmetic can be exercised without AMS on the path."""
    try:
        from scm.base import KFFile
    except ImportError:                                 # older AMS layouts
        try:
            from scm.libbase import KFFile
        except ImportError:
            raise SystemExit('rose_kf_bridge: KFFile not found. Run this '
                             'under amspython.')
    return KFFile


def open_kf(path, writable=False):
    """Open an existing KF file for reading, and refuse anything else.

    KFFile's default open mode is OpenMode.Any, which CREATES the file when it
    is not there. Opening a mistyped path therefore leaves an empty KF file
    behind and then reports it as having no sections, which looks exactly like
    a file whose layout this tool does not understand. Hence Old, and the
    checks below: everything here is read-only and nothing should ever be
    created.
    """
    import os

    KFFile = kffile_class()
    if not os.path.exists(path):
        raise SystemExit('rose_kf_bridge: %s does not exist. Nothing has been '
                         'created.' % path)

    is_kffile = getattr(KFFile, 'is_kffile', None)
    if is_kffile is not None:
        try:
            if not is_kffile(path):
                raise SystemExit('rose_kf_bridge: %s is not a KF file.' % path)
        except SystemExit:
            raise
        except Exception:                               # noqa: BLE001
            pass

    mode = getattr(KFFile, 'OpenMode', None)
    old = getattr(mode, 'Old', None) if mode is not None else None
    if writable:
        # autosave, so that a section created here is flushed even if a later
        # write raises. Old still means the file must exist.
        try:
            kf = KFFile(path, True, old) if old is not None else KFFile(path, True)
        except TypeError:
            kf = KFFile(path, old) if old is not None else KFFile(path)
    else:
        kf = KFFile(path, old) if old is not None else KFFile(path)

    num_sections = getattr(kf, 'num_sections', None)
    if num_sections is not None:
        try:
            count = num_sections() if callable(num_sections) else num_sections
        except Exception:                               # noqa: BLE001
            count = None
        if count == 0:
            raise SystemExit(
                'rose_kf_bridge: %s has no sections -- it is an empty KF file. '
                'An earlier run of this tool, or of anything else that opened '
                'it with the default OpenMode.Any, will have created it. Delete '
                'it and point this at the real results file.' % path)
    return kf


# ----------------------------------------------------------------------------
#  KF access
# ----------------------------------------------------------------------------

def kf_read(kf, section, variable, default=None):
    """Read one KF variable, returning default when it is absent.

    KFFile has carried more than one accessor across AMS versions, so the forms
    are tried in turn and the one that works is remembered. Without this a
    signature mismatch looks exactly like an empty file.
    """
    exists = getattr(kf, 'var_exists', None)
    if exists is not None:
        try:
            if not exists(section, variable):
                return default
        except Exception:                               # noqa: BLE001
            pass

    forms = (lambda: kf.read(section, variable),
             lambda: kf.read_reals_np(section, variable),
             lambda: kf.read_ints_np(section, variable),
             lambda: kf[(section, variable)],
             lambda: kf.read('%s%%%s' % (section, variable)))
    order = _READ_FORM[0]
    if order is not None:
        forms = (forms[order],) + tuple(f for i, f in enumerate(forms)
                                        if i != order)
    for i, form in enumerate(forms):
        try:
            value = form()
        except Exception:                               # noqa: BLE001
            continue
        if _READ_FORM[0] is None:
            _READ_FORM[0] = i
        return value
    return default


_READ_FORM = [None]


def kf_require(kf, section, variable):
    value = kf_read(kf, section, variable)
    if value is None:
        raise KeyError('%s%%%s not on file' % (section, variable))
    return value


def as_array(value, dtype):
    return np.atleast_1d(np.asarray(value, dtype=dtype))


def as_scalar(value):
    return int(np.atleast_1d(np.asarray(value)).ravel()[0])


# ----------------------------------------------------------------------------
#  Geometry: internal atom order, types and coordinates
# ----------------------------------------------------------------------------

#  The atoms carrying the basis are in ADF's internal order, which groups them
#  by type and so is generally not the input order.
#
#  Two coordinate arrays are therefore written, and they are not interchangeable:
#  /input/molecule/geometry takes Molecule%Coords, in INPUT order, because that
#  is what every other ROSE front end writes and what the user's xyz file is
#  numbered by; the AO centres take Geometry%xyz, in INTERNAL order, because the
#  shell list has to follow the AO order that npart and Eigen-Bas assume. The
#  two describe the same nuclei, and ROSE associates AOs with atoms by
#  coordinate rather than by index, so this is consistent -- but confusing the
#  two arrays is silent, hence the permutation check below.
#
#  The atom-to-type map comes from Geometry%nqptr, which is ntyp+1 long and
#  gives the first atom of each type in internal order. That is unambiguous.
#  Geometry%fragment and atomtype index is kept as a fallback for files without
#  nqptr: it is 2*natoms long, the fragment index per atom followed by the type
#  index per atom, in two blocks rather than interleaved -- the same layout as
#  Geometry%atom order index.

COORD_INTERNAL = ['xyz', 'Atom coordinates', 'coordinates']
COORD_INPUT = [('Molecule', 'Coords'), ('Molecule', 'coords'),
               ('Geometry', 'xyz InputOrder')]
TYPE_CANDIDATES = ['fragment and atomtype index', 'atomtype index', 'atomtype']


def read_geometry(kf, ntypes):
    """Return (coord_internal[natoms,3], atom_type[natoms] 1-based)."""
    natoms = as_scalar(kf_require(kf, 'Geometry', 'nr of atoms'))

    coord = None
    for name in COORD_INTERNAL:
        value = kf_read(kf, 'Geometry', name)
        if value is not None:
            value = as_array(value, float).ravel()
            if value.size == 3 * natoms:
                coord = value.reshape((natoms, 3))
                break
    if coord is None:
        raise KeyError('no Geometry coordinate array of length 3*natoms found; '
                       'tried %s. Run inspect and add the right name to '
                       'COORD_INTERNAL.' % COORD_INTERNAL)

    ntyp_kf = kf_read(kf, 'Geometry', 'nr of atomtypes')
    if ntyp_kf is not None and as_scalar(ntyp_kf) != ntypes:
        raise ValueError('Basis implies %d atom types but Geometry reports %d'
                         % (ntypes, as_scalar(ntyp_kf)))

    atype = _atom_types(kf, natoms, ntypes)
    return coord, atype


def _atom_types(kf, natoms, ntypes):
    """Type index per atom, in internal order, 1 based."""
    nqptr = kf_read(kf, 'Geometry', 'nqptr')
    if nqptr is not None:
        nqptr = as_array(nqptr, int).ravel()
        if nqptr.size == ntypes + 1 and nqptr[-1] - 1 == natoms:
            atype = np.zeros(natoms, dtype=int)
            for t in range(ntypes):
                atype[nqptr[t] - 1:nqptr[t + 1] - 1] = t + 1
            if atype.min() >= 1:
                return atype
        raise ValueError('Geometry%%nqptr is %s, which does not describe %d '
                         'atoms in %d types' % (nqptr.tolist(), natoms, ntypes))

    for name in TYPE_CANDIDATES:
        value = kf_read(kf, 'Geometry', name)
        if value is None:
            continue
        value = as_array(value, int).ravel()
        if value.size == 2 * natoms:
            atype = value[natoms:].copy()        # two blocks, not interleaved
        elif value.size == natoms:
            atype = value.copy()
        else:
            continue
        if atype.min() >= 1 and atype.max() <= ntypes:
            return atype
    raise KeyError('no usable Geometry atom-type source; tried nqptr and %s'
                   % TYPE_CANDIDATES)


# ----------------------------------------------------------------------------
#  Basis
# ----------------------------------------------------------------------------

class AdfBasis(object):
    """The Slater basis of one ADF calculation, expanded over the atoms.

    After construction:
        shells   list of dicts with l, n, zeta, coord
        nao      number of Cartesian functions, which must equal Basis%naos
    """

    def __init__(self, kf):
        self.naos = as_scalar(kf_require(kf, 'Basis', 'naos'))
        nbaspt = as_array(kf_require(kf, 'Basis', 'nbaspt'), int).ravel()
        nbptr = as_array(kf_require(kf, 'Basis', 'nbptr'), int).ravel()
        self.ntypes = nbaspt.size - 1
        if nbptr.size - 1 != self.ntypes:
            raise ValueError('nbaspt and nbptr disagree on the number of '
                             'atom types (%d vs %d)'
                             % (nbaspt.size - 1, nbptr.size - 1))

        nqbas = as_array(kf_require(kf, 'Basis', 'nqbas'), int).ravel()
        lqbas = as_array(kf_require(kf, 'Basis', 'lqbas'), int).ravel()
        alfbas = as_array(kf_require(kf, 'Basis', 'alfbas'), float).ravel()
        kx = as_array(kf_require(kf, 'Basis', 'kx'), int).ravel()
        ky = as_array(kf_require(kf, 'Basis', 'ky'), int).ravel()
        kz = as_array(kf_require(kf, 'Basis', 'kz'), int).ravel()
        kr = as_array(kf_require(kf, 'Basis', 'kr'), int).ravel()

        self._check_core(kf)
        self._check_shape(nbaspt, nbptr, nqbas, lqbas, kx, ky, kz, kr)

        self.nbaspt = nbaspt
        self.nqbas = nqbas
        self.lqbas = lqbas
        self.alfbas = alfbas
        self.basnrm = as_array(kf_read(kf, 'Basis', 'basnrm', default=[]),
                               float).ravel()
        self.bnorm = as_array(kf_read(kf, 'Basis', 'bnorm', default=[]),
                              float).ravel()
        self.kxyz = (kx, ky, kz, kr)

        self.coord, self.atom_type = read_geometry(kf, self.ntypes)
        self.shells = self._expand()
        self.nao = sum(ncomp(s['l']) for s in self.shells)
        if self.nao != self.naos:
            raise ValueError('expanded basis has %d functions but Basis%%naos '
                             'is %d. The atom-to-type assignment or the atom '
                             'order is wrong.' % (self.nao, self.naos))

    @staticmethod
    def _check_core(kf):
        ncset = kf_read(kf, 'Basis', 'ncset')
        if ncset is None:
            ncset = kf_read(kf, 'ActiveFrag', 'ncset')
        if ncset is not None and as_scalar(ncset) > 0:
            raise NotImplementedError(
                'this calculation has %d frozen-core basis sets. The core '
                'functions are not part of the valence set written here, so '
                'the converted basis would not span what Eigen-Bas is '
                'expressed in. Rerun with core None.' % as_scalar(ncset))

    def _check_shape(self, nbaspt, nbptr, nqbas, lqbas, kx, ky, kz, kr):
        """Verify the two documented identities before trusting anything."""
        nshell = nbaspt[-1] - 1
        ncomp_total = nbptr[-1] - 1
        if nqbas.size < nshell or lqbas.size < nshell:
            raise ValueError('nqbas/lqbas shorter than nbaspt implies')
        for arr, name in ((kx, 'kx'), (ky, 'ky'), (kz, 'kz'), (kr, 'kr')):
            if arr.size < ncomp_total:
                raise ValueError('%s shorter than nbptr implies' % name)

        # component by component: kx+ky+kz = l and kr = n-1-l
        pos = 0
        for t in range(self.ntypes):
            for ish in range(nbaspt[t] - 1, nbaspt[t + 1] - 1):
                l = int(lqbas[ish])
                n = int(nqbas[ish])
                for _ in range(ncomp(l)):
                    if kx[pos] + ky[pos] + kz[pos] != l:
                        raise ValueError(
                            'component %d has kx+ky+kz = %d but lqbas says %d'
                            % (pos + 1, kx[pos] + ky[pos] + kz[pos], l))
                    if kr[pos] != n - 1 - l:
                        raise ValueError(
                            'component %d has kr = %d but n-1-l = %d'
                            % (pos + 1, kr[pos], n - 1 - l))
                    pos += 1
        if pos != ncomp_total:
            raise ValueError('shells of all types expand to %d components but '
                             'nbptr says %d. The basis is probably contracted '
                             'or carries functions this reader does not know.'
                             % (pos, ncomp_total))

    def _expand(self):
        """Shells atom by atom, in the Basis (not Eigen-Bas) order."""
        shells = []
        for iatom, t in enumerate(self.atom_type):
            for ish in range(self.nbaspt[t - 1] - 1, self.nbaspt[t] - 1):
                shells.append(dict(l=int(self.lqbas[ish]),
                                   n=int(self.nqbas[ish]),
                                   zeta=float(self.alfbas[ish]),
                                   coord=self.coord[iatom].copy(),
                                   atom=iatom + 1))
        return shells

    def spherical_deficit(self):
        """How many fewer MOs than AOs the Cartesian-to-spherical step costs.

        ADF defines the basis in Cartesian powers but diagonalises in the solid
        harmonics, as DIRAC does, so each shell contributes
        (l+1)(l+2)/2 - (2l+1) = l(l-1)/2 functions that never carry an orbital:
        one per d shell, three per f, six per g. For a d shell the discarded
        direction is the contaminant x^2+y^2+z^2 = r^2 exp(-alpha r), an s-type
        function of principal quantum number n.

        So nao - nmo is structural and predictable, and is NOT a sign of linear
        dependence: the Cartesian overlap of a real ADF basis is well
        conditioned (smallest eigenvalue 7.2e-3 for nitrogen DZP), because that
        contaminant is not in the s set. A discrepancy here means something else
        is going on -- a frozen core, symmetry blocking, or a genuine dependency
        removal -- and the converted file should not be trusted until it is
        understood.
        """
        return sum(s['l'] * (s['l'] - 1) // 2 for s in self.shells)

    def check_normalisation(self, tol=1.0e-3):
        """Compare ADF's own bnorm against the convention slater_overlap uses.

        A mismatch here means the converted basis spans the right functions
        with the wrong scale, which the overlap gate would also catch but less
        legibly. Returns the largest relative deviation, or None when ADF did
        not store bnorm.
        """
        if self.bnorm.size == 0:
            return None
        kx, ky, kz, _ = self.kxyz
        worst = 0.0
        pos = 0
        for t in range(self.ntypes):
            for ish in range(self.nbaspt[t] - 1, self.nbaspt[t + 1] - 1):
                n = int(self.nqbas[ish])
                zeta = float(self.alfbas[ish])
                for _ in range(ncomp(int(self.lqbas[ish]))):
                    mine = sto_norm(n, zeta) * cart_angular_norm(
                        int(kx[pos]), int(ky[pos]), int(kz[pos]))
                    theirs = float(self.bnorm[pos])
                    if theirs != 0.0:
                        worst = max(worst, abs(mine / theirs - 1.0))
                    pos += 1
        return worst


def ncomp(l):
    return (l + 1) * (l + 2) // 2


def dfact(n):
    f = 1.0
    while n >= 2:
        f *= n
        n -= 2
    return f


def sto_norm(n, zeta):
    import math
    return (2.0 * zeta) ** (n + 0.5) / math.sqrt(math.factorial(2 * n))


def cart_angular_norm(kx, ky, kz):
    import math
    l = kx + ky + kz
    return math.sqrt(dfact(2 * l + 1) /
                     (4.0 * math.pi * dfact(2 * kx - 1) *
                      dfact(2 * ky - 1) * dfact(2 * kz - 1)))


# ----------------------------------------------------------------------------
#  Orbitals
# ----------------------------------------------------------------------------

def read_orbitals(kf, naos, spin='A'):
    """MO coefficients in BASIS order, eigenvalues, from the A section."""
    nmo = as_scalar(kf_require(kf, 'A', 'nmo_' + spin))
    flat = as_array(kf_require(kf, 'A', 'Eigen-Bas_' + spin), float).ravel()
    if flat.size != naos * nmo:
        raise ValueError('Eigen-Bas_%s has %d elements, expected naos*nmo = %d'
                         % (spin, flat.size, naos * nmo))
    # The former ADF reader did reshape(flat,(naos,nmo)) in Fortran order, so
    # the AO index runs fastest.
    coeff = flat.reshape((nmo, naos))

    npart = kf_read(kf, spin, 'npart')
    if npart is None:
        npart = kf_read(kf, 'A', 'npart')
    if npart is None:
        raise KeyError('A%npart not on file; without it the coefficients '
                       'cannot be put in Basis order')
    npart = as_array(npart, int).ravel()
    if npart.size != naos:
        raise ValueError('npart has %d entries, expected naos = %d'
                         % (npart.size, naos))
    if sorted(npart.tolist()) != list(range(1, naos + 1)):
        raise ValueError('npart is not a permutation of 1..naos')

    permuted = np.zeros_like(coeff)
    permuted[:, npart - 1] = coeff

    eps = as_array(kf_require(kf, 'A', 'eps_' + spin), float).ravel()
    if eps.size != nmo:
        raise ValueError('eps_%s has %d entries, expected %d'
                         % (spin, eps.size, nmo))
    return permuted, eps, npart


def read_overlap(kf, naos):
    """ADF's AO overlap, lower triangle packed row by row, in Basis order."""
    tri = kf_read(kf, 'Matrices', 'Smat')
    if tri is None:
        return None
    tri = as_array(tri, float).ravel()
    if tri.size != naos * (naos + 1) // 2:
        raise ValueError('Matrices%%Smat has %d elements, expected %d'
                         % (tri.size, naos * (naos + 1) // 2))
    return tri


# ----------------------------------------------------------------------------
#  Writing
# ----------------------------------------------------------------------------

def write_h5(path, basis, coeff, eps, natoms, nuc_charge, geometry,
             overlap_tri=None):
    import h5py

    nmo, nao = coeff.shape
    with h5py.File(path, 'w') as f:
        f.create_dataset('/input/molecule/n_atoms', data=[natoms])
        f.create_dataset('/input/molecule/nuc_charge', data=nuc_charge)
        # Angstrom here, bohr in /aobasis/center; see the note at the top.
        f.create_dataset('/input/molecule/geometry',
                         data=(geometry * ANGSTROM_PER_BOHR).ravel())

        orbmom, n_prim, n_cont = [], [], []
        center, exponents, contractions, prinqn = [], [], [], []
        for s in basis.shells:
            orbmom.append(s['l'] + 1)        # ROSE stores l+1
            n_prim.append(1)                 # ADF Slater shells are uncontracted
            n_cont.append(1)
            center.extend(s['coord'])
            exponents.append(s['zeta'])
            contractions.append(1.0)
            prinqn.append(s['n'])

        ao = '/aobasis/'
        f.create_dataset(ao + 'angular', data=[1])       # Cartesian
        f.create_dataset(ao + 'basis_type', data=[2])    # Slater
        f.create_dataset(ao + 'n_shells', data=[len(basis.shells)])
        f.create_dataset(ao + 'n_ao', data=[basis.nao])
        f.create_dataset(ao + 'orbmom', data=orbmom, maxshape=(None,))
        f.create_dataset(ao + 'n_prim', data=n_prim, maxshape=(None,))
        f.create_dataset(ao + 'n_cont', data=n_cont, maxshape=(None,))
        f.create_dataset(ao + 'n', data=prinqn, maxshape=(None,))
        f.create_dataset(ao + 'center', data=np.asarray(center, float),
                         maxshape=(None,))
        f.create_dataset(ao + 'exponents', data=np.asarray(exponents, float),
                         maxshape=(None,))
        f.create_dataset(ao + 'contractions',
                         data=np.asarray(contractions, float), maxshape=(None,))

        mo = '/mobasis/'
        f.create_dataset(mo + 'n_basis', data=[nao])
        f.create_dataset(mo + 'n_mo', data=[nmo])
        f.create_dataset(mo + 'nz', data=[1])
        f.create_dataset(mo + 'eigenvalues', data=eps, chunks=True,
                         maxshape=(None,))
        f.create_dataset(mo + 'orbitals', data=coeff.ravel(), chunks=True,
                         maxshape=(None,))

        if overlap_tri is not None:
            # ADF's own AO overlap, kept so that the Slater overlap routine can
            # be checked against it; see the gate in hdf5_interface.F90. Lower
            # triangle, row by row: element (i,j) with j <= i sits at
            # i*(i-1)/2 + j, one based.
            f.create_dataset('/validation/ao_overlap_tri', data=overlap_tri,
                             chunks=True, maxshape=(None,))


# ----------------------------------------------------------------------------
#  Commands
# ----------------------------------------------------------------------------

INSPECT_SECTIONS = ['General', 'Geometry', 'Molecule', 'Basis', 'A',
                    'Symmetry', 'Matrices', 'ActiveFrag']


def kf_sections(kf):
    """Section name -> list of variable names, however this KFFile reports it.

    scm.base._internal.KFFile exposes sections() and variables(section), which
    is the explicit pair and is tried first. get_skeleton() exists on that class
    too but has come back empty, and older builds offer neither, so the other
    shapes are kept as fallbacks: a dict of lists, a dict of dicts, or a flat
    sequence of 'Section%var' strings. Returning an empty map is a real answer,
    and the caller says so rather than printing nothing.
    """
    # 1. sections() + variables(section)
    lister = getattr(kf, 'sections', None)
    varlister = getattr(kf, 'variables', None)
    if lister is not None and varlister is not None:
        try:
            names = lister() if callable(lister) else lister
            out = {}
            for raw in names:
                sec = str(raw).strip()
                vars_ = None
                for key in (raw, sec):
                    try:
                        vars_ = [str(v).strip() for v in varlister(key)]
                        break
                    except Exception:                   # noqa: BLE001
                        continue
                out[sec] = vars_ if vars_ is not None else []
            if out:
                return out, 'sections()/variables()'
        except Exception as exc:                        # noqa: BLE001
            print('  (sections() failed: %s)' % exc)

    # 2. a single call returning the whole layout
    skeleton, how = None, None
    for name in ('get_skeleton', 'get_contents', '_sections'):
        getter = getattr(kf, name, None)
        if getter is None:
            continue
        try:
            skeleton = getter() if callable(getter) else getter
            how = name
            if skeleton:
                break
        except Exception as exc:                        # noqa: BLE001
            print('  (%s() failed: %s)' % (name, exc))
    if not skeleton:
        return {}, how

    out = {}
    try:
        items = list(skeleton.items())
    except AttributeError:
        for entry in skeleton:
            if isinstance(entry, (tuple, list)) and len(entry) == 2:
                sec, var = entry
            elif isinstance(entry, str) and '%' in entry:
                sec, var = entry.split('%', 1)
            else:
                continue
            out.setdefault(str(sec).strip(), []).append(str(var).strip())
        return out, how

    for sec, vars_ in items:
        try:
            names = list(vars_)
        except TypeError:
            names = []
        out[str(sec).strip()] = [str(v).strip() for v in names]
    return out, how


def cmd_probe(args):
    """Print the KFFile binding's own signatures and try to open the file.

    These are pybind11 bindings, so every method carries its signature in its
    docstring. Reading them beats guessing at the call, which is how this
    command came to exist.
    """
    KFFile = kffile_class()
    print('class: %s.%s' % (KFFile.__module__, KFFile.__name__))
    print('NOTE: the open attempts below use OpenMode.Any among others, which '
          'CREATES a missing file. Only run probe on a file that exists.')

    for name in ('__init__', 'sections', 'variables', 'read', 'get_skeleton',
                 'num_sections', 'num_variables', 'section_exists',
                 'var_length', 'var_type', 'read_reals_np', 'read_ints_np'):
        member = getattr(KFFile, name, None)
        doc = getattr(member, '__doc__', None)
        print('\n--- %s' % name)
        if doc:
            for line in str(doc).strip().split('\n')[:8]:
                print('    %s' % line)
        else:
            print('    (no docstring)')

    mode = getattr(KFFile, 'OpenMode', None)
    if mode is not None:
        print('\n--- OpenMode')
        members = getattr(mode, '__members__', None)
        if members:
            print('    %s' % ', '.join(members))
        else:
            print('    %s' % [a for a in dir(mode) if not a.startswith('_')])

    print('\n--- opening %s' % args.rkf)
    attempts = [('KFFile(path)', lambda: KFFile(args.rkf))]
    if mode is not None:
        for member in getattr(mode, '__members__', {}):
            attempts.append(('KFFile(path, OpenMode.%s)' % member,
                             lambda m=member: KFFile(args.rkf,
                                                     getattr(mode, m))))
    for label, build in attempts:
        try:
            kf = build()
        except Exception as exc:                        # noqa: BLE001
            print('    %-34s construction failed: %s' % (label, exc))
            continue
        bits = []
        for probe in ('is_kffile', 'num_sections', 'num_variables'):
            value = getattr(kf, probe, None)
            try:
                bits.append('%s=%s' % (probe,
                                       value() if callable(value) else value))
            except Exception as exc:                    # noqa: BLE001
                bits.append('%s=<%s>' % (probe, exc))
        try:
            secs = list(kf.sections())
            bits.append('sections=%d %s' % (len(secs), secs[:6]))
        except Exception as exc:                        # noqa: BLE001
            bits.append('sections=<%s>' % exc)
        print('    %-34s %s' % (label, '; '.join(bits)))
    return 0


def cmd_inspect(args):
    kf = open_kf(args.rkf)
    print('file: %s' % args.rkf)

    sections, how = kf_sections(kf)
    if not sections:
        print('\nCould not list the sections of this file.')
        print('  KFFile class : %s.%s' % (type(kf).__module__,
                                          type(kf).__name__))
        print('  accessor used: %s' % how)
        public = [a for a in dir(kf) if not a.startswith('__')]
        print('  it offers    : %s' % ', '.join(public))
        for probe in ('is_kffile', 'num_sections', 'num_variables'):
            value = getattr(kf, probe, None)
            try:
                print('  %-13s: %s' % (probe,
                                       value() if callable(value) else value))
            except Exception as exc:                    # noqa: BLE001
                print('  %-13s: <%s>' % (probe, exc))
        print('\nRun the probe subcommand on the same file: it prints the '
              'binding signatures and tries every constructor form.')
        return 1
    print('listed by %s: %d sections\n' % (how, len(sections)))

    wanted = [s for s in INSPECT_SECTIONS if s in sections]
    missing = [s for s in INSPECT_SECTIONS if s not in sections]
    others = sorted(s for s in sections if s not in INSPECT_SECTIONS)

    for section in wanted:
        print('[%s]' % section)
        for name in sorted(sections[section]):
            value = kf_read(kf, section, name)
            if value is None:
                print('  %-34s (unreadable)' % name)
                continue
            arr = np.atleast_1d(np.asarray(value).ravel())
            head = ', '.join(str(v) for v in arr[:10])
            if arr.size > 10:
                head += ', ...'
            print('  %-34s n=%-6d %s' % (name, arr.size, head))
        print('')

    if missing:
        print('sections the conversion wants but this file does not have: %s'
              % ', '.join(missing))
    if others:
        print('other sections on the file: %s' % ', '.join(others))

    print('\nThe conversion needs: Geometry nr of atoms, a coordinate array of '
          'length 3*natoms in INTERNAL order, and an atom-type index per atom; '
          'Basis naos/nbaspt/nbptr/nqbas/lqbas/alfbas/kx/ky/kz/kr; A '
          'nmo_A/Eigen-Bas_A/eps_A/npart; Molecule AtomicNumbers and coords in '
          'INPUT order.')
    return 0


def cmd_forward(args):
    kf = open_kf(args.rkf)

    nspin = kf_read(kf, 'General', 'nspin')
    if nspin is not None and as_scalar(nspin) != 1:
        raise NotImplementedError('only restricted (nspin = 1) calculations '
                                  'are converted; this one has nspin = %d'
                                  % as_scalar(nspin))
    symlab = kf_read(kf, 'Symmetry', 'symlab')
    if symlab is not None and str(symlab).strip().upper() not in ('A', 'NOSYM'):
        raise NotImplementedError('only NOSYM calculations are converted; this '
                                  'one reports symmetry %r' % symlab)

    basis = AdfBasis(kf)
    coeff, eps, npart = read_orbitals(kf, basis.naos)

    natoms = as_scalar(kf_require(kf, 'Geometry', 'nr of atoms'))
    nuc_charge = as_array(kf_require(kf, 'Molecule', 'AtomicNumbers'),
                          float).ravel()
    geometry = None
    for section, name in COORD_INPUT:
        value = kf_read(kf, section, name)
        if value is not None:
            geometry = as_array(value, float).ravel()
            break
    if geometry is None or geometry.size != 3 * natoms:
        raise KeyError('no Molecule coordinate array of length 3*natoms found')
    geometry = geometry.reshape((natoms, 3))

    # The AO centres are in internal order and this array is in input order.
    # They must at least be the same set of nuclei; if they are not, one of the
    # two names resolved to the wrong variable and everything downstream is
    # quietly wrong.
    def _sorted_rows(a):
        a = np.round(np.asarray(a, float), 8)
        return a[np.lexsort((a[:, 2], a[:, 1], a[:, 0]))]

    if not np.allclose(_sorted_rows(geometry), _sorted_rows(basis.coord),
                       atol=1.0e-8):
        raise ValueError('the input-order and internal-order coordinate arrays '
                         'are not permutations of each other; check which '
                         'Geometry/Molecule variables were picked up')

    overlap_tri = None
    if args.t15:
        kf15 = open_kf(args.t15)
        naos15 = as_scalar(kf_require(kf15, 'Basis', 'naos'))
        if naos15 != basis.naos:
            raise ValueError('%s has naos = %d, the rkf has %d'
                             % (args.t15, naos15, basis.naos))
        overlap_tri = read_overlap(kf15, basis.naos)
        if overlap_tri is None:
            raise KeyError('Matrices%Smat not on %s' % args.t15)

    write_h5(args.output, basis, coeff, eps, natoms, nuc_charge, geometry,
             overlap_tri=overlap_tri)

    print('%s -> %s' % (args.rkf, args.output))
    print('  atoms          %d in %d types' % (natoms, basis.ntypes))
    print('  shells         %d' % len(basis.shells))
    print('  AO functions   %d (Cartesian)' % basis.nao)
    print('  MOs            %d' % coeff.shape[0])
    deficit = basis.nao - coeff.shape[0]
    expected = basis.spherical_deficit()
    if deficit == expected:
        print('                 %d fewer than the AOs, which is exactly the '
              'Cartesian-to-spherical' % deficit)
        print('                 reduction for this basis (sum of l(l-1)/2 over '
              'shells)')
    else:
        print('                 **WARNING** %d fewer than the AOs, but the '
              'Cartesian-to-spherical' % deficit)
        print('                 reduction accounts for only %d. The remaining '
              '%d are unexplained:' % (expected, deficit - expected))
        print('                 a frozen core, symmetry blocking or a genuine '
              'dependency removal.')
        print('                 Do not trust this file until that is '
              'understood.')
    reordered = int((npart != np.arange(1, basis.naos + 1)).sum())
    print('  npart          %d of %d functions moved' % (reordered, basis.naos))
    print('  units          geometry written in Angstrom, AO centres in bohr')
    worst = basis.check_normalisation()
    if worst is None:
        print('  normalisation  bnorm not on file, not checked')
    else:
        print('  normalisation  agrees with ADF bnorm to %.2e (relative)'
              % worst)
    if overlap_tri is None:
        print('  validation     no --t15 given, no reference overlap stored')
    else:
        print('  validation     ADF overlap stored; ROSE will check its own '
              'Slater overlaps against it')
    return 0


def cmd_backward(args):
    """ribo.h5 -> a Rose Data section on an rkf.

    Takes over from write_mos_to_ADF and write_data_rose, which did this from
    inside ROSE until the file-based ADF route was removed. Everything is taken
    from the h5 file except froc, derived as they derived it from
    General%electrons and Molecule%Charge on the rkf itself, and npart, which is
    needed to put the coefficients back into ADF's own AO order.
    """
    import os
    import shutil

    import h5py

    with h5py.File(args.ribo, 'r') as f:
        def need(path):
            if path not in f:
                raise SystemExit('rose_kf_bridge: %s has no %s. An h5 file '
                                 'written before the RIBO metadata was carried '
                                 'cannot be converted back.' % (args.ribo, path))
            return f[path][()]

        nmo = as_scalar(need('/mobasis/n_mo'))
        nao = as_scalar(need('/mobasis/n_basis'))
        nz = as_scalar(need('/mobasis/nz'))
        if nz != 1:
            raise NotImplementedError('only real restricted orbitals are '
                                      'written back; this file has nz = %d' % nz)
        coeff = np.asarray(need('/mobasis/orbitals'), float).reshape((nmo, nao))
        eps = np.asarray(need('/mobasis/eigenvalues'), float).ravel()
        # The RIBO metadata sits beside the MO dimensions rather than in a
        # group of its own: ROSE's writer cannot create a top-level group.
        nfrag = as_scalar(need('/mobasis/ribo_nfrag'))
        occ_per_frag = as_array(need('/mobasis/ribo_occ_per_frag'), int).ravel()
        vir_per_frag = as_array(need('/mobasis/ribo_vir_per_frag'), int).ravel()
        ncore = as_scalar(need('/mobasis/ribo_n_core'))
        nribo = as_scalar(need('/mobasis/ribo_n_active'))
        fock = np.asarray(need('/mobasis/ribo_fock'), float)
        # The molecule, to check that the rkf really describes this system.
        # ribo.h5 inherits /input/molecule/ from CHECKPOINT.h5, in Angstrom.
        h5_z = as_array(need('/input/molecule/nuc_charge'), float).ravel()
        h5_xyz = (np.asarray(need('/input/molecule/geometry'), float).ravel()
                  / ANGSTROM_PER_BOHR).reshape((-1, 3))

    if occ_per_frag.size != nfrag or vir_per_frag.size != nfrag:
        raise ValueError('nfrag is %d but the group counts have %d and %d '
                         'entries' % (nfrag, occ_per_frag.size,
                                      vir_per_frag.size))
    if nribo != nmo - ncore:
        raise ValueError('fock_ribo is %d wide but nmo - n_core is %d'
                         % (nribo, nmo - ncore))
    if fock.size != nribo * nribo:
        raise ValueError('fock_ribo has %d elements, expected %d'
                         % (fock.size, nribo * nribo))

    # Settle the container before writing a single byte.
    #
    # When a template is given it is ALWAYS copied, overwriting whatever the
    # output path held. The alternative, reusing an existing file, is precisely
    # the trap write_mos_to_ADF fell into: it copied only when the target was
    # absent, so a second case running in the same directory wrote its
    # Rose Data into the first case's rkf and produces a file whose Basis
    # section describes one molecule and whose Rose Data describes another.
    # Nothing downstream notices.
    if args.template:
        if not os.path.exists(args.template):
            raise SystemExit('rose_kf_bridge: template %s does not exist'
                             % args.template)
        if (os.path.exists(args.output)
                and os.path.samefile(args.template, args.output)):
            raise SystemExit('rose_kf_bridge: --template and --output are the '
                             'same file; the template must stay untouched')
        replaced = os.path.exists(args.output)
        shutil.copyfile(args.template, args.output)
        print('%s %s -> %s' % ('replaced' if replaced else 'copied',
                               args.template, args.output))
    elif not os.path.exists(args.output):
        raise SystemExit('rose_kf_bridge: %s does not exist and no --template '
                         'was given to copy from' % args.output)
    else:
        print('WARNING writing into the existing %s, which this run did not '
              'create;' % args.output)
        print('        pass --template to start from a known rkf instead')

    kf = open_kf(args.output, writable=True)

    naos = as_scalar(kf_require(kf, 'Basis', 'naos'))
    if naos != nao:
        raise ValueError('%s has naos = %d, the h5 file has %d'
                         % (args.output, naos, nao))

    # naos agreeing is necessary but not sufficient: two runs on different
    # geometries of the same molecule, or on different molecules with the same
    # basis size, would pass it. Compare the nuclei as well.
    rkf_z = as_array(kf_require(kf, 'Molecule', 'AtomicNumbers'), float).ravel()
    rkf_xyz = None
    for section, name in COORD_INPUT:
        value = kf_read(kf, section, name)
        if value is not None:
            value = as_array(value, float).ravel()
            if value.size == 3 * rkf_z.size:
                rkf_xyz = value.reshape((-1, 3))
                break
    if rkf_xyz is None:
        raise KeyError('no coordinate array of length 3*natoms on %s'
                       % args.output)
    if rkf_z.size != h5_z.size:
        raise ValueError('%s has %d atoms, the h5 file has %d'
                         % (args.output, rkf_z.size, h5_z.size))

    def _fingerprint(z, xyz):
        a = np.column_stack((np.round(z, 6), np.round(xyz, 6)))
        return a[np.lexsort((a[:, 3], a[:, 2], a[:, 1], a[:, 0]))]

    if not np.allclose(_fingerprint(rkf_z, rkf_xyz),
                       _fingerprint(h5_z, h5_xyz), atol=1.0e-6):
        raise ValueError('%s and %s do not describe the same nuclei; the rkf '
                         'is from a different calculation' % (args.output,
                                                              args.ribo))

    # Back into ADF's AO order. Forward did c_basis[npart[i]-1] = c_adf[i];
    # this is that read backwards.
    npart = as_array(kf_require(kf, 'A', 'npart'), int).ravel()
    if sorted(npart.tolist()) != list(range(1, nao + 1)):
        raise ValueError('A%npart is not a permutation of 1..naos')
    coeff_adf = coeff[:, npart - 1]

    # froc, exactly as write_mos_to_ADF derived it.
    nel = int(as_array(kf_require(kf, 'General', 'electrons'), float).ravel()[0])
    chg = int(as_array(kf_require(kf, 'Molecule', 'Charge'), float).ravel()[0])
    nocc = (nel - chg) // 2
    froc = np.zeros(nmo)
    froc[:nocc] = 2.0

    section = 'Rose Data'
    exists = getattr(kf, 'section_exists', None)
    had_section = bool(exists(section)) if exists is not None else False
    if had_section:
        kf.delete_section(section)
    kf.create_section(section)

    # Eigen-Bas_A is nao-fastest, which is what ravel() of an (nmo,nao) array
    # gives; see the forward direction, where the same layout is read.
    kf.write(section, 'nmo_A', int(nmo))
    kf.write(section, 'Eigen-Bas_A', coeff_adf.ravel().tolist())
    kf.write(section, 'eps_A', eps.tolist())
    kf.write(section, 'froc_A', froc.tolist())
    kf.write(section, 'nfrag', int(nfrag))
    kf.write(section, 'ribo_occ_per_frag_A', occ_per_frag.tolist())
    kf.write(section, 'ribo_vir_per_frag_A', vir_per_frag.tolist())
    kf.write(section, 'Fock_RIBO_A', fock.ravel().tolist())
    save = getattr(kf, 'save', None)
    if save is not None:
        save()

    print('%s -> %s' % (args.ribo, args.output))
    print('  section        Rose Data (%s)'
          % ('replaced' if had_section else 'created'))
    print('  nmo_A          %d, of which %d frozen core' % (nmo, ncore))
    print('  Eigen-Bas_A    %d x %d, put back in ADF AO order' % (nmo, nao))
    print('  froc_A         %d electrons, charge %d, %d doubly occupied'
          % (nel, chg, nocc))
    print('  groups         %d, occ %s, vir %s'
          % (nfrag, occ_per_frag.tolist(), vir_per_frag.tolist()))
    print('  Fock_RIBO_A    %d x %d' % (nribo, nribo))
    return 0


def cmd_fixture(args):
    """Convert a whole plams_workdir into the fixture set ROSE reads.

    The supermolecule becomes CHECKPOINT.h5 and each fragment becomes NNN.h5,
    named by atomic number, because rose_read_fragments builds that name as
    nint(bascharge) in i3.3 format. Everything is converted while the -2 route
    still works; afterwards these inputs cannot be regenerated.
    """
    import os

    if not os.path.isdir(args.workdir):
        raise SystemExit('rose_kf_bridge: %s is not a directory' % args.workdir)
    os.makedirs(args.output, exist_ok=True)

    jobs = []
    molecule = os.path.join(args.workdir, 'Molecule_B1_MO1')
    if not os.path.isdir(molecule):
        raise SystemExit('rose_kf_bridge: no Molecule_B1_MO1 in %s'
                         % args.workdir)
    jobs.append((molecule, 'CHECKPOINT.h5'))

    for name in sorted(os.listdir(args.workdir)):
        if not name.startswith('frag'):
            continue
        path = os.path.join(args.workdir, name)
        rkf = os.path.join(path, 'adf.rkf')
        if not os.path.isfile(rkf):
            continue
        kf = open_kf(rkf)
        z = as_array(kf_require(kf, 'Molecule', 'AtomicNumbers'), int).ravel()
        if z.size != 1:
            raise SystemExit('%s holds %d atoms; atomic fragment files are '
                             'named by a single atomic number' % (path, z.size))
        jobs.append((path, '%03d.h5' % int(z[0])))

    for path, out in jobs:
        rkf = os.path.join(path, 'adf.rkf')
        t15 = os.path.join(path, 'TAPE15')
        sub = argparse.Namespace(rkf=rkf,
                                 output=os.path.join(args.output, out),
                                 t15=t15 if os.path.isfile(t15) else None)
        print('\n=== %s -> %s' % (os.path.basename(path), out))
        cmd_forward(sub)

    genibo = os.path.join(args.output, 'INPUT_GENIBO')
    if args.input and os.path.isfile(args.input):
        with open(args.input) as fh:
            text = fh.read()
        text = text.replace('\nadf\n', '\nh5\n').replace('\nams\n', '\nh5\n')
        with open(genibo, 'w') as fh:
            fh.write(text)
        print('\nINPUT_GENIBO taken from %s with .FILE_FORMAT set to h5'
              % args.input)
    else:
        print('\nNo --input given, so no INPUT_GENIBO was written. Copy the one '
              'the original run used and change .FILE_FORMAT to h5.')

    print('%d files in %s' % (len(jobs), args.output))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    sub = parser.add_subparsers(dest='command')

    p = sub.add_parser('probe', help='print the KFFile binding signatures and '
                                     'try to open a file; use when inspect '
                                     'cannot list the sections')
    p.add_argument('rkf')
    p.set_defaults(func=cmd_probe)

    p = sub.add_parser('inspect', help='print the KF data the conversion uses')
    p.add_argument('rkf')
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser('forward', help='adf.rkf -> CHECKPOINT.h5')
    p.add_argument('rkf')
    p.add_argument('-o', '--output', default='CHECKPOINT.h5')
    p.add_argument('--t15', help='TAPE15 holding Matrices%%Smat, stored as the '
                                 'reference for the overlap gate')
    p.set_defaults(func=cmd_forward)

    p = sub.add_parser('backward', help='ribo.h5 -> a Rose Data section on an rkf')
    p.add_argument('ribo', help='the h5 file ROSE wrote, normally ribo.h5')
    p.add_argument('-o', '--output', required=True,
                   help='the rkf to write into, normally '
                        'Molecule_B1_MO1_ribo.rkf')
    p.add_argument('--template', help='rkf to copy from when the output does '
                                      'not exist yet, normally '
                                      'Molecule_B1_MO1.rkf')
    p.set_defaults(func=cmd_backward)

    p = sub.add_parser('fixture', help='a whole plams_workdir -> a fixture set')
    p.add_argument('workdir', help='the plams_workdir of an ADF run')
    p.add_argument('-o', '--output', required=True, help='where to write them')
    p.add_argument('--input', help='the INPUT_GENIBO of the original run, '
                                   'copied with .FILE_FORMAT set to h5')
    p.set_defaults(func=cmd_fixture)

    args = parser.parse_args(argv)
    if not getattr(args, 'func', None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
