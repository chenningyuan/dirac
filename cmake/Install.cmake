message("-- DIRAC: CMAKE_INSTALL_PREFIX is \"${CMAKE_INSTALL_PREFIX}\" ")

install(
    DIRECTORY
    DESTINATION
    dirac
    )

install(
    DIRECTORY
    DESTINATION
    dirac/utils
    )

install(
    TARGETS dirac.x
    DESTINATION dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    FILES
    ${PROJECT_SOURCE_DIR}/pam
    DESTINATION
    dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    FILES
    ${PROJECT_SOURCE_DIR}/utils/common.py
    DESTINATION
    dirac/utils
    )

install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_SOURCE_DIR}/basis_dalton
    DESTINATION
    dirac
    PATTERN .svn EXCLUDE
    )
