#deactivate exatensor if mpi is not available
if(ENABLE_MPI AND (NOT ENABLE_64BIT_INTEGERS))
    option(ENABLE_EXATENSOR "Enable ExaTensor library" ON)
else()
    option(ENABLE_EXATENSOR "Enable ExaTensor library" OFF)
endif()
# configure with (e.g.):
# $ export WRAP=NOWRAP && export TOOLKIT=GNU && export BUILD_TYPE=DEV && export MPILIB=OPENMPI && export BLASLIB=ATLAS && export GPU_CUDA=NOCUDA && export PATH_OPENMPI=/home/bast/openmpi/1.10.2 && ./setup
# in future we could use cmake discovery to infer this information, now it has to be given by the user

# this is configurable to make it possible to configure on machines
# on which the repository is located but where network cannot be accessed
# in this case EXATENSOR_GIT_REPO_LOCATION can be set to point to a path on the hard disk
# ./setup --cmake-options="-DEXATENSOR_GIT_REPO_LOCATION='/path/to/exatensor/'"
set(EXATENSOR_GIT_REPO_LOCATION "https://github.com/ornl-qci/exatensor.git" CACHE STRING "ExaTENSOR Git repository location")
set(EXATENSOR_GIT_HASH 2ab28a45499aa0e0dea4a3efc29aa15a4012a905 CACHE STRING "ExaTENSOR Git repository hash in use")

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

    include_directories(${PROJECT_BINARY_DIR}/exatensor/src/exatensor/include)
endif()

message(STATUS "Enable ExaTensor library: ${ENABLE_EXATENSOR}")
