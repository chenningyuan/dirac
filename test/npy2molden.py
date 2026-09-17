#!/usr/bin/env python3
"""
Write a Molden file from a psi4 wavefunction archive.

The ROSE psi4 launchers leave a psi4 wavefunction archive behind (ribo.npy), whose
MO coefficients are the RIBO orbitals copied out of ribo.h5 -- the two are bitwise
identical, the .npy simply being a psi4 container around the same numbers. This
script turns that archive into a Molden file for visual inspection.

psi4 does the part that has to be right: a Molden file must state the AO basis in
exactly the order its MO coefficient columns assume, and the archive carries only
the basis *name* (STO-3G-DECON, say), not the primitives. psi4 rebuilds the basis
and writes the [GTO] and [MO] blocks from one wavefunction object, so they cannot
disagree.

Two things psi4 cannot do from an archive, which this script fixes afterwards:

  * Occupation numbers. A wavefunction restored with Wavefunction.from_file has no
    occupation vectors, so psi4 warns and writes zeros for every orbital. They are
    filled in here from the archive's own nalpha/nbeta.

  * Redundant beta orbitals. ROSE's psi4 runs use an ROHF reference, so psi4 emits
    an alpha block and an identical beta block -- twice as many orbitals as the
    calculation has. When the archive says both spins share their orbitals and the
    counts match, the beta block is dropped and the alpha occupations become 2.0.

Usage
-----
    python npy2molden.py ribo.npy
    python npy2molden.py ribo.npy -o my_orbitals.molden
    python npy2molden.py ribo.npy --occupied-only
    python npy2molden.py ribo.npy --quiet        # psi4's own output to a file
    python npy2molden.py ribo.npy --raw          # skip the fixes, keep psi4's file

psi4 prints as it loads the wavefunction; that is left alone by default. Note that
psi4's set_output_file also opens a log at <name>.log, so it cannot be pointed at
os.devnull.

Orbital order
-------------
ROSE writes the occupied orbitals first and then the virtuals, each ordered by
recanonisation group and recanonised within that group. Orbital energies therefore
restart at every group boundary and are not monotonic across the file. That is
expected; a viewer that re-sorts by energy will scramble the group structure.
"""

import argparse
import os
import sys

HEADER_KEYS = ("sym=", "ene=", "spin=", "occup=")


# --------------------------------------------------------------------------- #
#  Molden text handling. No psi4 here, so all of it is directly testable.
# --------------------------------------------------------------------------- #

def split_molden(path):
    """Split a Molden file into everything up to [MO], and the MO records after it."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("[mo]"):
            start = i
            break
    if start is None:
        return lines, []
    preamble = lines[:start + 1]
    records, current = [], None
    for line in lines[start + 1:]:
        if line.strip().lower().startswith("sym="):
            if current is not None:
                records.append(current)
            current = [line]
        elif current is not None:
            current.append(line)
        else:
            preamble.append(line)
    if current is not None:
        records.append(current)
    return preamble, records


def record_spin(record):
    for line in record:
        stripped = line.strip().lower()
        if stripped.startswith("spin="):
            return stripped.split("=", 1)[1].strip()
    return "alpha"


def record_occupation(record):
    for line in record:
        stripped = line.strip().lower()
        if stripped.startswith("occup="):
            try:
                return float(stripped.split("=", 1)[1])
            except ValueError:
                return 0.0
    return 0.0


def set_occupation(record, value):
    """Return a copy of the record with Occup= set, inserting the line if absent."""
    out, seen, last_header = [], False, -1
    for i, line in enumerate(record):
        stripped = line.strip().lower()
        if stripped.startswith("occup="):
            out.append(" Occup= %.6f" % value)
            seen = True
            last_header = i
        else:
            if any(stripped.startswith(k) for k in HEADER_KEYS):
                last_header = i
            out.append(line)
    if not seen:
        out.insert(last_header + 1, " Occup= %.6f" % value)
    return out


def repair_molden(path, nalpha, nbeta, same_orbs, occupied_only=False):
    """Fill in occupations, and drop the beta block when it merely repeats alpha.

    Returns a dict describing what was done.
    """
    preamble, records = split_molden(path)
    if not records:
        return {"records": 0, "dropped_beta": False, "occupied_only": False}

    alpha = [r for r in records if record_spin(r).startswith("a")]
    beta = [r for r in records if record_spin(r).startswith("b")]

    drop_beta = bool(beta) and same_orbs and nalpha == nbeta and len(alpha) == len(beta)
    single_block = drop_beta or not beta

    if single_block:
        # One block of spatial orbitals: doubly occupied up to nbeta, singly
        # occupied from there to nalpha. For a closed shell that is just 2.0
        # across the occupied set; for a restricted open shell it puts 1.0 on
        # the singly occupied orbitals, which is what Molden expects.
        ndouble, nsingle = min(nalpha, nbeta), max(nalpha, nbeta)
        fixed = []
        for i, record in enumerate(alpha):
            if i < ndouble:
                occ = 2.0
            elif i < nsingle:
                occ = 1.0
            else:
                occ = 0.0
            fixed.append(set_occupation(record, occ))
    else:
        # Genuinely unrestricted: both spins kept, one electron each.
        fixed = []
        for spin_records, nocc in ((alpha, nalpha), (beta, nbeta)):
            for i, record in enumerate(spin_records):
                fixed.append(set_occupation(record, 1.0 if i < nocc else 0.0))

    if occupied_only:
        fixed = [r for r in fixed if record_occupation(r) > 0.0]

    body = [line for record in fixed for line in record]
    with open(path, "w") as fh:
        fh.write("\n".join(preamble + body) + "\n")
    return {"records": len(fixed), "dropped_beta": drop_beta,
            "occupied_only": occupied_only}


def summarize_molden(path):
    """Read a Molden file back and report what it contains."""
    info = {"natoms": 0, "nmo": 0, "nocc": 0, "nelec": 0.0, "angular_tags": []}
    section = None
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("["):
                tag = stripped.split("]")[0].lstrip("[").strip().upper()
                if tag in ("5D", "7F", "9G", "5D7F", "5D10F"):
                    if tag not in info["angular_tags"]:
                        info["angular_tags"].append(tag)
                section = tag
                continue
            if not stripped:
                continue
            if section == "ATOMS":
                info["natoms"] += 1
            elif section == "MO":
                low = stripped.lower()
                if low.startswith("ene"):
                    info["nmo"] += 1
                elif low.startswith("occup"):
                    try:
                        occ = float(stripped.split("=", 1)[1])
                    except (IndexError, ValueError):
                        continue
                    info["nelec"] += occ
                    if occ > 0.0:
                        info["nocc"] += 1
    return info


def read_archive(stem):
    import numpy as np
    data = np.load(stem + ".npy", allow_pickle=True).item()
    return {"nmo": int(data["int"]["nmo"]),
            "nalpha": int(data["int"]["nalpha"]),
            "nbeta": int(data["int"]["nbeta"]),
            "same_orbs": bool(data["boolean"]["same_a_b_orbs"]),
            "puream": bool(data["boolean"]["basispuream"]),
            "reference": str(data["string"]["name"])}


# --------------------------------------------------------------------------- #

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Write a Molden file from a psi4 wavefunction archive (.npy).")
    parser.add_argument("archive",
                        help="the psi4 archive, e.g. ribo.npy (.npy suffix optional)")
    parser.add_argument("-o", "--output", default=None,
                        help="output file (default: the archive name with .molden)")
    parser.add_argument("--occupied-only", action="store_true",
                        help="keep only the occupied orbitals")
    parser.add_argument("--raw", action="store_true",
                        help="leave psi4's file exactly as written, occupations and "
                             "duplicated beta block included")
    parser.add_argument("--quiet", action="store_true",
                        help="send psi4's own output to <output>_psi4.out")
    args = parser.parse_args(argv)

    stem = args.archive[:-4] if args.archive.endswith(".npy") else args.archive
    if not os.path.isfile(stem + ".npy"):
        sys.exit("no such archive: " + stem + ".npy")
    output = args.output or (stem + ".molden")

    try:
        import psi4
    except ImportError:
        sys.exit("psi4 is not importable. Run this inside the environment the "
                 "calculation was run in.")

    meta = read_archive(stem)

    # Redirecting psi4's output is a convenience only, and must never cost the
    # conversion: set_output_file also opens a log at <name>.log, which is why
    # os.devnull cannot be used (it becomes /dev/null.log).
    if args.quiet:
        try:
            psi4.set_output_file(os.path.splitext(output)[0] + "_psi4.out", False)
        except Exception as exc:                              # noqa: BLE001
            print("note: could not redirect psi4 output (%s)" % exc)

    wfn = psi4.core.Wavefunction.from_file(stem)

    # Try the supported wfn.write_molden first, then fall back to psi4.molden.
    #
    # A wavefunction restored with Wavefunction.from_file is a bare
    # psi4.core.Wavefunction, and write_molden calls self.occupation_a(), which such
    # an object does not have. psi4.molden() is deprecated but more forgiving: it
    # warns and writes zeros, and the occupations are filled in below in either
    # case. Testing hasattr is not enough - the method exists, it just fails inside,
    # so the call itself has to be attempted.
    writer_used, errors = None, []
    for name, call in (("wfn.write_molden", lambda: wfn.write_molden(output)),
                       ("psi4.molden", lambda: psi4.molden(wfn, output))):
        try:
            call()
            writer_used = name
            break
        except Exception as exc:                              # noqa: BLE001
            errors.append("%s -> %s: %s" % (name, type(exc).__name__, exc))
    if writer_used is None:
        sys.exit("psi4 could not write a Molden file:\n  " + "\n  ".join(errors))

    if args.raw:
        repaired = {"records": 0, "dropped_beta": False, "occupied_only": False}
    else:
        repaired = repair_molden(output, meta["nalpha"], meta["nbeta"],
                                 meta["same_orbs"], args.occupied_only)

    info = summarize_molden(output)
    print()
    print("wrote %s   (%s reference, via %s)"
          % (output, meta["reference"], writer_used))
    if errors:
        print("  note             : %s" % errors[0])
    print("  atoms            : %d" % info["natoms"])
    print("  orbitals written : %d" % info["nmo"])
    print("  with occupation  : %d   (%.1f electrons)" % (info["nocc"], info["nelec"]))
    print("  angular tags     : %s"
          % (", ".join(info["angular_tags"]) if info["angular_tags"]
             else "none -> cartesian d/f"))
    if repaired["dropped_beta"]:
        print("  note             : dropped the duplicate beta block and set the "
              "alpha occupations to 2.0")

    expected_mo = meta["nalpha"] if args.occupied_only else meta["nmo"]
    expected_occ = meta["nalpha"]
    if not args.raw:
        if info["nmo"] != expected_mo:
            print("  WARNING: expected %d orbitals, file has %d"
                  % (expected_mo, info["nmo"]))
        if info["nocc"] != expected_occ:
            print("  WARNING: archive says %d occupied orbitals, file marks %d"
                  % (expected_occ, info["nocc"]))
    if meta["puream"] and not info["angular_tags"]:
        print("  WARNING: the archive uses a spherical basis but the file carries "
              "no [5D]/[7F] tag")
    return 0


if __name__ == "__main__":
    sys.exit(main())
