option(ENABLE_EXATENSOR "Enable ExaTensor library" ON)

# configure with (e.g.):
# $ export WRAP=NOWRAP && export TOOLKIT=GNU && export BUILD_TYPE=DEV && export MPILIB=OPENMPI && export BLASLIB=ATLAS && export GPU_CUDA=NOCUDA && export PATH_OPENMPI=/home/bast/openmpi/1.10.2 && ./setup
# in future we could use cmake discovery to infer this information, now it has to be given by the user

if(ENABLE_EXATENSOR)
    ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY git@gitlab.com:DmitryLyakh/ExaTensor.git
        GIT_TAG dirac
        CONFIGURE_COMMAND true  # currently no configure command
        BUILD_COMMAND cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/ && make
        INSTALL_COMMAND true  # currently no install command
        INSTALL_DIR "${PROJECT_BINARY_DIR}/exatensor/install"
        )

    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/TALSH/libtalsh.a -lgomp -lcudart -lcublas)

    include_directories(${PROJECT_BINARY_DIR}/exatensor/src/exatensor/TALSH)
endif()

message(STATUS "Enable ExaTensor library: ${ENABLE_EXATENSOR}")
