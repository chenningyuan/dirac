cmake_dependent_option(ENABLE_XML "XML input/output module" ON
                       "ENABLE_UNRELEASED" OFF)

if(ENABLE_XML)
    add_definitions(-DMOD_XML)
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        src/gp/xml/xml_file.F90
        src/gp/xml/xml_structure.F90
        src/gp/xml/xml_parser.F90
        )
else()
    set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_XML)
endif()

message(STATUS "XML input/output module: ${ENABLE_XML}")
