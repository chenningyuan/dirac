# set minimum cmake version
cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

#Needed on macosx so that it finds dynamically linked libraries
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# project name
project(rose LANGUAGES Fortran C)  # C support is required by a dependency of LAPACK detection

# if CMAKE_BUILD_TYPE undefined, we set it to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(ROSE_LINALG_LIBRARIES)
if(TARGET LAPACK::LAPACK)
  list(APPEND ROSE_LINALG_LIBRARIES LAPACK::LAPACK)
else()
  list(APPEND ROSE_LINALG_LIBRARIES ${LAPACK_LIBRARIES})
endif()
if(TARGET BLAS::BLAS)
  list(APPEND ROSE_LINALG_LIBRARIES BLAS::BLAS)
else()
  list(APPEND ROSE_LINALG_LIBRARIES ${BLAS_LIBRARIES})
endif()

set(HDF5_FIND_DEBUG TRUE)
find_package(HDF5 COMPONENTS C HL)
if(HDF5_FOUND)
   option(ENABLE_HDF5 "Enable HDF5" ON)
   message(STATUS "HFD5 found, enabling HDF5 interface ")
   message(STATUS "HDF5 include directory location: ${HDF5_INCLUDE_DIRS} ")
   include_directories(${HDF5_INCLUDE_DIRS})
else()
   option(ENABLE_HDF5 "Enable HDF5" OFF)
   message(STATUS "HDF5 interface can not be used, no HDF5 installation found")
endif()

add_subdirectory(src)


#### Updates from the original ROSE CMakeLists.txt

# Find the Python interpreter. Avoiding here the use of FindPython [e.g, find_package(Python COMPONENTS Interpreter REQUIRED)].
find_program(PYTHON_INTERPRETER python3)
if(NOT PYTHON_INTERPRETER)
    message(FATAL_ERROR "Could not find python3 interpreter.")
endif()

# Now, install and/or upgrade runtest and genibo
message("Installing and/or upgrading required python packages...\n")

# Setting up relevant variables
set(ROSE_PYTHON_SCRIPTS_PATH ${PROJECT_SOURCE_DIR}/python_scripts)
set(ROSE_ASE_PATH ${PROJECT_SOURCE_DIR}/ase_rose)
# Manipulating wildcards in cmake
file(GLOB atomfrag_files "${ROSE_PYTHON_SCRIPTS_PATH}/atomfrag_*.py")
file(GLOB molfrag_files "${ROSE_PYTHON_SCRIPTS_PATH}/molfrag_*.py")
file(GLOB groupfrag_files "${ROSE_PYTHON_SCRIPTS_PATH}/groupfrag_*.py")

# Execute pip if python is found.
#
# Each package gets its own execute_process call. Several COMMAND entries in a
# single call are not run in sequence: CMake chains them into a pipeline, so the
# output of all but the last disappears into a pipe and RESULT_VARIABLE reports
# only the exit status of the last one. A failing install then goes unnoticed
# until something tries to import the package.
function(rose_pip_install)
  set(pip_args ${ARGN})
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} -m pip install --user ${pip_args} --upgrade
    RESULT_VARIABLE pip_user_result
  )
  if(NOT pip_user_result EQUAL 0)
    message(STATUS "pip install --user ${pip_args} failed, retrying without --user")
    execute_process(
      COMMAND ${PYTHON_INTERPRETER} -m pip install ${pip_args} --upgrade
      RESULT_VARIABLE pip_plain_result
    )
    if(NOT pip_plain_result EQUAL 0)
      message(FATAL_ERROR
        "Failed to pip install ${pip_args}. Check the python interpreter and pip.")
    endif()
  endif()
endfunction()

rose_pip_install(runtest)
rose_pip_install(-e ${ROSE_PYTHON_SCRIPTS_PATH}/genibo/)
rose_pip_install(-e ${ROSE_ASE_PATH}[testing])

message("\n Runtest and genibo (python_scripts/genibo/) python packages installed/upgraded successfully!\n")

# Copying python files to build directory
message("Copying atomfrag and molfrag python scripts to ${CMAKE_CURRENT_BINARY_DIR} directory...\n")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy ${atomfrag_files} ${molfrag_files} ${groupfrag_files} ${CMAKE_CURRENT_BINARY_DIR}
  RESULT_VARIABLE FILES_COPY_CHECK_RESULT
)
if(NOT FILES_COPY_CHECK_RESULT)
  message("Files copied successfully!\n")
else()
  message(FATAL_ERROR "Failed to copy atomfrag_*.py and molfrag_*.py files from ${ROSE_PYTHON_SCRIPTS_PATH}.")
endif()

#### final of Updates from the original ROSE CMakeLists.txt


# Tests are launched with the bare name python3, not with ${PYTHON_INTERPRETER}.
# The latter is an absolute path resolved at configure time and baked into
# CTestTestfile.cmake, which pins the tests to whatever environment happened to
# be active when cmake ran. Resolving python3 through PATH at test time instead
# means the tests follow the environment that is active when they are run, and
# runtest_config.py then hands that same interpreter to the launcher it spawns.
macro(rose_test _name _label)
    add_test(${_name} python3 ${PROJECT_SOURCE_DIR}/test/${_name}/test --binary-dir=${PROJECT_BINARY_DIR} --work-dir=${PROJECT_BINARY_DIR}/test/${_name} --verbose)
    if(NOT "${_label}" STREQUAL "")
      set_tests_properties(${_name} PROPERTIES LABELS "${_label}")
    endif()
endmacro()

# activate ctest
include(CTest)
enable_testing()
rose_test(dirac_atomfrags "dirac" "")
rose_test(dirac_molfrags "dirac" "")
rose_test(psi4_molfrags "psi4" "")
rose_test(psi4_atomfrags "psi4" "")
rose_test(psi4_groupfrags "psi4" "")
rose_test(pyscf_molfrags "pyscf" "")
rose_test(pyscf_atomfrags "pyscf" "")
rose_test(pyscf_libfrags "pyscf" "")
rose_test(ams_atomfrags "ams" "")
rose_test(ams_molfrags "ams" "")

# adding ase_rose tests
set(ASE_WORKDIR ${PROJECT_SOURCE_DIR}/ase_rose/tests/)
add_test(NAME rose_ase_pyscf_tests
         COMMAND python3 -m pytest -v
          ${PROJECT_SOURCE_DIR}/ase_rose/tests/test_pyscf_water.py
          ${PROJECT_SOURCE_DIR}/ase_rose/tests/test_pyscf_water-ammonia.py
)
add_test(NAME rose_ase_psi4_tests
         COMMAND python3 -m pytest -v
          ${PROJECT_SOURCE_DIR}/ase_rose/tests/test_psi4_water-ammonia.py
)

# set the working as rose_ase/test
set_tests_properties(rose_ase_pyscf_tests rose_ase_psi4_tests PROPERTIES WORKING_DIRECTORY ${ASE_WORKDIR})
