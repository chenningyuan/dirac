set(COMPILATION_INFO_FILE ${CMAKE_SOURCE_DIR}/gp/compilation_info.F90)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/compilation_info.py.in
    ${CMAKE_BINARY_DIR}/compilation_info.py
    )
execute_process(
    COMMAND python ${CMAKE_BINARY_DIR}/compilation_info.py
    TIMEOUT 1
    OUTPUT_VARIABLE COMPILER_VERSIONS
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# split output variable into two EQUAL substrings with Fortran and C compiler versions (solution for cmake 2.6)
string(LENGTH ${COMPILER_VERSIONS} COMP_VERS_STRLEN)
math(EXPR  halflength  "${COMP_VERS_STRLEN}/2")
string(SUBSTRING "${COMPILER_VERSIONS}"  "0"               "${halflength}"  FORTRAN_COMPILER_VERSION  )
string(SUBSTRING "${COMPILER_VERSIONS}"  "${halflength}"   "${halflength}"  C_COMPILER_VERSION        )

# find out git revision info, both for Fortran subroutine and CMake-build output
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/git_revision_info.py.in
    ${CMAKE_BINARY_DIR}/git_revision_info.py
    )
# one can get stdout either into file or into variable, can not have both 
execute_process(
    COMMAND python ${CMAKE_BINARY_DIR}/git_revision_info.py
    TIMEOUT 1
    OUTPUT_FILE ${CMAKE_SOURCE_DIR}/gp/git_revision_info.F90
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# copy content of file "git_revision_info.F90" into GIT_REVISION variable
file(READ ${CMAKE_SOURCE_DIR}/gp/git_revision_info.F90 GIT_REVISION)
# clean the string a little
string(REPLACE "subroutine print_git_revision_info()" "" GIT_REVISION  "${GIT_REVISION}")
string(REPLACE "end subroutine"                       "" GIT_REVISION  "${GIT_REVISION}")
string(REPLACE "write(6, '(a)') "                     "" GIT_REVISION  "${GIT_REVISION}")
string(REPLACE "'"                                    "" GIT_REVISION  "${GIT_REVISION}")
#message("GIT_REVISION_MOD=${GIT_REVISION_MOD}")



