# - Find a LAPACK library
#
# This module will first look in LAPACK_ROOT before considering the default
# system pahts.
# The linker language can be defined by setting the varable LAPACK_LANG
#
# This module defines:
#
#  LAPACK_INCLUDE_DIRS Where to find lapack.h (or equivalent)
#  LAPACK_LIBRARIES Libraries to link against to use LAPACK
#  LAPACK_FOUND Defined if LAPACK is available
#  HAVE_LAPACK To be used in #ifdefs
#  LAPACK_H Name of LAPACK header file
#
# None of the above will be defined unless LAPACK can be found.
#
#=============================================================================
# Copyright 2011 Jonas Juselius <jonas.juselius@uit.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

include(MathLibFunctions)

if (LAPACK_LANG STREQUAL "C")
# fixme
endif()

set(path_suffixes lib lib/atlas)
if (LAPACK_LANG STREQUAL "C")
# fixme
else()
	set(lapack_libs lapack)
endif()

find_math_header("lapack")
find_math_libs("lapack")

unset(lapack_libs)
unset(path_suffixes)

cache_math_result(Atlas "lapack")

