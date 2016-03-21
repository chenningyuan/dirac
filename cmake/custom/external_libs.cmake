
#message("cmake/custom/external_libs.cmake: from src/CMakeFiles.txt... obtained EXTERNAL_LIBS=${EXTERNAL_LIBS}")

set(EXTERNAL_LIBS
    ${EXTERNAL_LIBS}
    ${MATH_LIBS}
    )

if(NOT EXPLICIT_LIBS STREQUAL "off")
    set(EXTERNAL_LIBS
        ${EXTERNAL_LIBS}
        ${EXPLICIT_LIBS}
        )
endif()

#message("cmake/custom/external_libs.cmake: added mathlibs/explicitlibs EXTERNAL_LIBS=${EXTERNAL_LIBS}")

