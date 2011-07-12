#
#  add the list of standalone_executables and their single-file sources 
#
#   input: ${standalone_executables} pairs of executable, source code
#   output:  ${executables}
#
if(VERBOSE_OUTPUT)
    message("-- List of standalone executables and their source codes:")
endif()
set( i "0") # initialize counter
foreach(code ${standalone_executables})
 # radovan: what is this math stuff and why is it so complicated?
 math( EXPR i "${i} + 1 ")     # i:=i+1
 math( EXPR i1 "2*(${i}/2)")   # i1=2*(i/2)
 set(this_src ${code})
 set(this_exec ${this_exec1})
 if ("${i1}" STREQUAL "${i}")
    # if i is even, retrieve executable and its source file
    if(VERBOSE_OUTPUT)
        message ("      ${this_exec}  ${this_src}")
    endif()
    # add standalone executable and its source file
    add_executable(
    ${this_exec}
    ${this_src}
    )
    # add next executable to the list
    set(executables ${executables}  ${this_exec})
 endif()
 # shift variable in the cycle ...
 set(this_exec1 ${this_src})
endforeach(code)
