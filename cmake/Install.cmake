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

foreach(
    BINARY
    ${executables}
    )
    install(
        TARGETS ${BINARY}
        DESTINATION dirac
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

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
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_SOURCE_DIR}/basis_dalton
    DESTINATION
    dirac
    PATTERN .svn EXCLUDE
    )
