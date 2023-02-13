#Miro: remove problematic linking library
list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES  "pthread")

# Ivan Hrasko's determination of compiler versions - very good tool !
include(FindCompilerVersion)

# Miro Ilias restored CMake printout, associated with build_info.h generation
include(build_info)

find_package(HDF5 COMPONENTS C HL)
