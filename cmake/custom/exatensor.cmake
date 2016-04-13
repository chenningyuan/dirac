option(ENABLE_EXATENSOR "Enable ExaTensor library" OFF)

# configure with (e.g.):
# $ ./setup --cmake-options="-DENABLE_EXATENSOR=ON -DEXATENSOR_CONFIGURATION='WRAP=NOWRAP;TOOLKIT=GNU;BUILD_TYPE=DEV;MPILIB=OPENMPI;BLASLIB=ATLAS;GPU_CUDA=NOCUDA;PATH_OPENMPI=/home/bast/openmpi/1.10.2'"
# in future we could use cmake discovery to infer this information, now it has to be given by the user

if(ENABLE_EXATENSOR)
    set(EXATENSOR_CONFIGURATION)  # by default empty

    set(_build_command
           cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/DDSS
        && ${EXATENSOR_CONFIGURATION} make
        && cd ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/GFC
        && ${EXATENSOR_CONFIGURATION} make)

    ExternalProject_Add(exatensor
        PREFIX "${PROJECT_BINARY_DIR}/exatensor"
        GIT_REPOSITORY git@gitlab.com:DmitryLyakh/ExaTensor.git
        GIT_TAG bde1c836  # later we will follow proper tag
        CONFIGURE_COMMAND true  # currently no configure command
        BUILD_COMMAND ${_build_command}
        INSTALL_COMMAND true  # currently no install command
        INSTALL_DIR "${PROJECT_BINARY_DIR}/exatensor/install"
        )

    set(EXATENSOR_LIBRARIES
        ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/DDSS/libddss.a
        ${PROJECT_BINARY_DIR}/exatensor/src/exatensor/GFC/libgfc.a
        )

    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${EXATENSOR_LIBRARIES})
endif()

message(STATUS "Enable ExaTensor library: ${ENABLE_EXATENSOR}")
