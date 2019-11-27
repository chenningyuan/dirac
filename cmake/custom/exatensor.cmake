#deactivate exatensor if mpi is not available
if(ENABLE_MPI AND (NOT ENABLE_64BIT_INTEGERS))
    option(ENABLE_EXATENSOR "Enable ExaTensor library" ON)
else()
    option(ENABLE_EXATENSOR "Enable ExaTensor library" OFF)
endif()
# configure with (e.g.):
# $ export WRAP=NOWRAP && export TOOLKIT=GNU && export BUILD_TYPE=DEV && export MPILIB=OPENMPI && export BLASLIB=ATLAS && export GPU_CUDA=NOCUDA && export PATH_OPENMPI=/home/bast/openmpi/1.10.2 && ./setup
# in future we could use cmake discovery to infer this information, now it has to be given by the user


if(ENABLE_EXATENSOR)
    ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY https://github.com/ornl-qci/exatensor.git
        GIT_TAG a3a550e67d2beba1750719cc1d9882a3c308f5b2
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
