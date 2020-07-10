#Deactivate ExaTensor if MPI is not available
if(ENABLE_MPI AND (NOT ENABLE_64BIT_INTEGERS))
    option(ENABLE_EXATENSOR "Enable ExaTensor library" ON)
else()
    option(ENABLE_EXATENSOR "Enable ExaTensor library" OFF)
endif()

# This is configurable to make it possible to configure on machines
# on which the repository resides locally but where network cannot be accessed
# in which case EXATENSOR_GIT_REPO_LOCATION can be set to point to a path on the hard disk
# ./setup --cmake-options="-DEXATENSOR_GIT_REPO_LOCATION='/path/to/exatensor/'"
set(EXATENSOR_GIT_REPO_LOCATION "https://gitlab.com/DmitryLyakh/ExaTensor.git" CACHE STRING "ExaTENSOR Git repository location")
set(EXATENSOR_GIT_HASH c0a96b2fca285b4994cae7635ea9cb0ed3e470cd CACHE STRING "ExaTENSOR Git repository hash in use")

if(ENABLE_EXATENSOR)
    ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY ${EXATENSOR_GIT_REPO_LOCATION}
        GIT_TAG ${EXATENSOR_GIT_HASH}
        CONFIGURE_COMMAND true  # currently no configure command
        BUILD_COMMAND cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/ && make
        INSTALL_COMMAND true  # currently no install command
        INSTALL_DIR "${PROJECT_BINARY_DIR}/exatensor/install"
        )

    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib/libexatensor.a)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib/libtalsh.a)
    #Make sure to manually add all ExaTENSOR dependencies here (BLAS, OpenMP, CUDA, C++):
    #set(EXTERNAL_LIBS ${EXTERNAL_LIBS} -L/sw/summit/essl/6.1.0-2/essl/6.1/lib64 -L/sw/summit/xl/16.1.1-5/xlC/16.1.1/lib -L/sw/summit/xl/16.1.1-5/xlf/16.1.1/lib -lessl -lxlf90_r -lxlfmath -L/sw/summit/cuda/10.1.243/lib64 -lcublas -lcudart -lnvToolsExt -lstdc++ -lgomp)

    include_directories(${PROJECT_BINARY_DIR}/exatensor/src/exatensor/include)
endif()

message(STATUS "Enable ExaTENSOR library: ${ENABLE_EXATENSOR}")
