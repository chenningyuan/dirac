"""
Script that updates copyright headers. (c) Radovan Bast.

Usage (run this in the Git repository root folder):
    $ python maintenance/update_copyright.py

The script modifies files in place.
"""

import os
import re
import io
from collections import defaultdict


def convert_to_c_comments(text):
    c_text = "\n".join([" *" + line[2:] for line in text.splitlines()])
    return f"\n/*\n{c_text}\n */\n"


copyright_header = {}

copyright_header[
    ("dirac", "fortran")
] = """
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
"""
copyright_header[("dirac", "c")] = convert_to_c_comments(
    copyright_header[("dirac", "fortran")]
)

copyright_header[
    ("dalton", "fortran")
] = """
!      Dalton, a molecular electronic structure program
!      Copyright (c) by the authors of Dalton.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
"""

copyright_header[("dalton", "c")] = convert_to_c_comments(
    copyright_header[("dalton", "fortran")]
)

files = defaultdict(list)

for path, _, _files in os.walk(os.getcwd()):
    for f in _files:
        if f.endswith(".F90") or f.endswith(".f90") or f.endswith(".F"):
            files["fortran"].append(os.path.join(path, f))
        if f.endswith(".c"):
            files["c"].append(os.path.join(path, f))

for lang in ["fortran", "c"]:
    for code in ["dirac", "dalton"]:
        if lang == "fortran":
            p = re.compile(
                f"(?<=!{code}_copyright_start).*(?=!{code}_copyright_end)",
                re.DOTALL,
            )
        else:
            p = re.compile(
                f"(?<=/\* {code}_copyright_start \*/).*(?=/\* {code}_copyright_end \*/)",
                re.DOTALL,
            )
        for file_name in files[lang]:
            with io.open(file_name, "r", encoding="ascii", errors="ignore") as f:
                s = f.read()
            if "_copyright_start" in s:
                s = re.sub(p, copyright_header[(code, lang)], s)
                with io.open(file_name, "w", encoding="ascii", errors="ignore") as f:
                    f.write(s)
