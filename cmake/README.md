The CMake infrastructure for this project was originally generated using [Autocmake]
by Radovan Bast, Roberto Di Remigio, Jonas Juselius and contributors.
The `update.py` Python script and the contents of the directories `autocmake` and `downloaded` are licensed
under the terms of the [BSD-3-Clause license], unless otherwise stated.

**Note: `update.py` is obsolete and should not be run.**
The top-level `setup` script and `CMakeLists.txt`, as well as modules in `downloaded`,
have since been edited by hand, and `autocmake.yml` no longer describes them completely.
Running `update.py` would overwrite these files with incomplete versions
(e.g. dropping `add_subdirectory(src)` and several `setup` options) and replace
modules in `downloaded` with different upstream versions.
Edit `setup`, `CMakeLists.txt` and the CMake modules directly instead.
The directory `autocmake` is still required, since `setup` imports from it.

[Autocmake]: http://autocmake.org
[BSD-3-Clause license]: https://tldrlegal.com/license/bsd-3-clause-license-(revised)
