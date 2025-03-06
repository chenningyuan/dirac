file(COPY ${PROJECT_SOURCE_DIR}/basis        DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/basis_dalton DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/basis_ecp    DESTINATION ${PROJECT_BINARY_DIR})
message(STATUS "Copied DIRAC basis set directories into the build directory")

file(COPY ${PROJECT_SOURCE_DIR}/utils/DIRACschema.txt    DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/utils/dirac_data.py      DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/utils/process_schema.py  DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/utils/merge_amf.py       DESTINATION ${PROJECT_BINARY_DIR})
message(STATUS "Copied data schema and python utilities into the build directory")

foreach(
    EXECUTABLE
    ${_list_of_executables}
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
    FILES ${PROJECT_BINARY_DIR}/pam
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

# workaround to install pam-dirac symlink:
# 1) copy pam to pam-dirac
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/pam ${PROJECT_BINARY_DIR}/pam-dirac)")
# 2) install real file pam-dirac
install(
    FILES ${PROJECT_BINARY_DIR}/pam-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )
# 3) remove real file
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")
# 4) create symlink
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/pam ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")

# workaround to install the utilities executable symlink:
# 1) copy *.x to *-dirac
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/diag.x ${PROJECT_BINARY_DIR}/diag-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/dirac_mointegral_export.x ${PROJECT_BINARY_DIR}/dirac_mointegral_export-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/mx2fit.x ${PROJECT_BINARY_DIR}/mx2fit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/pcmo_addlabels.x ${PROJECT_BINARY_DIR}/pcmo_addlabels-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/polfit.x ${PROJECT_BINARY_DIR}/polfit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/twofit.x ${PROJECT_BINARY_DIR}/twofit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/vibcal.x ${PROJECT_BINARY_DIR}/vibcal-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/dirac_data.py ${PROJECT_BINARY_DIR}/dirac_data-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/process_schema.py ${PROJECT_BINARY_DIR}/process_schema-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/merge_amf.py ${PROJECT_BINARY_DIR}/merge_amf-dirac)")
# 2) install real files *-dirac
install(
    FILES
    ${PROJECT_BINARY_DIR}/diag-dirac
    ${PROJECT_BINARY_DIR}/dirac_mointegral_export-dirac
    ${PROJECT_BINARY_DIR}/mx2fit-dirac
    ${PROJECT_BINARY_DIR}/pcmo_addlabels-dirac
    ${PROJECT_BINARY_DIR}/polfit-dirac
    ${PROJECT_BINARY_DIR}/twofit-dirac
    ${PROJECT_BINARY_DIR}/vibcal-dirac
    ${PROJECT_BINARY_DIR}/dirac_data-dirac
    ${PROJECT_BINARY_DIR}/process_schema-dirac
    ${PROJECT_BINARY_DIR}/merge_amf-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )
# 3) remove real files
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/diag-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/dirac_mointegral_export-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/mx2fit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/pcmo_addlabels-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/polfit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/twofit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/vibcal-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/dirac_data-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/process_schema-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/merge_amf-dirac)")
# 4) create symlink
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/diag.x ${CMAKE_INSTALL_PREFIX}/bin/diag-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/dirac_mointegral_export.x ${CMAKE_INSTALL_PREFIX}/bin/dirac_mointegral_export-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/mx2fit.x ${CMAKE_INSTALL_PREFIX}/bin/mx2fit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/pcmo_addlabels.x ${CMAKE_INSTALL_PREFIX}/bin/pcmo_addlabels-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/polfit.x ${CMAKE_INSTALL_PREFIX}/bin/polfit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/twofit.x ${CMAKE_INSTALL_PREFIX}/bin/twofit-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/vibcal.x ${CMAKE_INSTALL_PREFIX}/bin/vibcal-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/dirac_data.py ${CMAKE_INSTALL_PREFIX}/bin/dirac_data-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/process_schema.py ${CMAKE_INSTALL_PREFIX}/bin/process_schema-dirac)")
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/merge_amf.py ${CMAKE_INSTALL_PREFIX}/bin/merge_amf-dirac)")

# write git hash to build dir
file(WRITE ${PROJECT_BINARY_DIR}/GIT_HASH "${_git_last_commit_hash}")
# copy version info to install dir
install(
    FILES ${PROJECT_BINARY_DIR}/GIT_HASH ${PROJECT_SOURCE_DIR}/VERSION
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_SOURCE_DIR}/basis_dalton
    ${PROJECT_SOURCE_DIR}/basis_ecp
    DESTINATION share/dirac
    PATTERN .git EXCLUDE
    )

install(
    FILES
    ${PROJECT_SOURCE_DIR}/utils/dirac_data.py
    ${PROJECT_SOURCE_DIR}/utils/process_schema.py
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    FILES
    ${PROJECT_SOURCE_DIR}/utils/DIRACschema.txt
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )
