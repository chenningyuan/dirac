option(ENABLE_EXATENSOR "Enable ExaTensor library" ON)

# configure with (e.g.):
# $ export WRAP=NOWRAP && export TOOLKIT=GNU && export BUILD_TYPE=DEV && export MPILIB=OPENMPI && export BLASLIB=ATLAS && export GPU_CUDA=NOCUDA && export PATH_OPENMPI=/home/bast/openmpi/1.10.2 && ./setup
# in future we could use cmake discovery to infer this information, now it has to be given by the user

if(ENABLE_EXATENSOR)
    ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY https://gitlab.com/DmitryLyakh/ExaTensor.git
        GIT_TAG devel
        CONFIGURE_COMMAND true  # currently no configure command
        BUILD_COMMAND cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/ && make
        INSTALL_COMMAND true  # currently no install command
        INSTALL_DIR "${PROJECT_BINARY_DIR}/exatensor/install"
        )

    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib/libexatensor.a)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/lib/libtalsh.a -L/sw/summit/essl/6.1.0-2/essl/6.1/lib64 -L/sw/summit/xl/16.1.1-1/xlC/16.1.1/lib -L/sw/summit/xl/16.1.1-1/xlf/16.1.1/lib -lessl -lxlf90_r -lxlfmath -L/sw/summit/cuda/9.2.148/lib64 -lgomp -lcublas -lcudart -lnvToolsExt -lstdc++)

    include_directories(${PROJECT_BINARY_DIR}/exatensor/src/exatensor/include)
endif()

message(STATUS "Enable ExaTensor library: ${ENABLE_EXATENSOR}")
