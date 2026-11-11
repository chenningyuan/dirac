find_package(Git QUIET)

set(_git_last_commit_hash "unknown")
set(_git_last_commit_author "unknown")
set(_git_last_commit_date "unknown")
set(_git_branch "unknown")

if(GIT_FOUND)
  execute_process(
    # abbreviated hash, with "-dirty" appended if tracked files have local modifications
    # (the non-matching --match pattern makes describe ignore tags and fall back to the hash)
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty --match=NeVeRmAtCh
    OUTPUT_VARIABLE _git_last_commit_hash
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )

  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%aN -n 1
    OUTPUT_VARIABLE _git_last_commit_author
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )

  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%ad -n 1
    OUTPUT_VARIABLE _git_last_commit_date
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )

  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE _git_branch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
endif()

configure_file(
  ${INPUT_DIR}/git_info.h.in
  ${OUTPUT_DIR}/git_info.h
  @ONLY
  )

unset(_git_last_commit_hash)
unset(_git_last_commit_author)
unset(_git_last_commit_date)
unset(_git_branch)
