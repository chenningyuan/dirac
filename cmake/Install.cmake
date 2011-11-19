foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    install(
        TARGETS ${EXECUTABLE}
        DESTINATION share/dirac
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

install(
    FILES ${PROJECT_BINARY_DIR}/pam-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    FILES ${PROJECT_SOURCE_DIR}/utils/wrapper.py
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/basis ${PROJECT_SOURCE_DIR}/basis_dalton
    DESTINATION share/dirac
    PATTERN .git EXCLUDE
    )
