# ——— Helper function to add & register tests —————————————————————————
include_guard()

function(ppc_add_test test_name test_src USE_FLAG)
  if(${USE_FLAG})
    add_executable(${test_name} "${PROJECT_SOURCE_DIR}/${test_src}")
    enable_testing()
    add_test(NAME ${test_name} COMMAND ${test_name})
    install(TARGETS ${test_name} RUNTIME DESTINATION bin)
  endif()
endfunction()

# Collects implementations from settings.json in one pass: - all keys from
# PPC_IMPLEMENTATIONS marked "disabled" - keys marked "enabled" that also have
# an existing implementation directory
function(ppc_collect_implementations_from_settings SETTINGS_PATH SUBDIR
         OUT_ENABLED_IMPLEMENTATIONS OUT_DISABLED_IMPLEMENTATIONS)
  set(ENABLED_IMPLEMENTATIONS "")
  set(DISABLED_IMPLEMENTATIONS "")

  if(EXISTS "${SETTINGS_PATH}")
    file(READ "${SETTINGS_PATH}" SETTINGS_JSON_CONTENT)

    foreach(IMPL IN LISTS PPC_IMPLEMENTATIONS)
      string(
        JSON
        IMPL_STATUS
        ERROR_VARIABLE
        IMPL_STATUS_ERROR
        GET
        "${SETTINGS_JSON_CONTENT}"
        tasks
        "${IMPL}")
      if(IMPL_STATUS_ERROR)
        continue()
      endif()

      string(TOLOWER "${IMPL_STATUS}" IMPL_STATUS)
      if(IMPL_STATUS STREQUAL "disabled")
        list(APPEND DISABLED_IMPLEMENTATIONS "${IMPL}")
      elseif(IMPL_STATUS STREQUAL "enabled"
             AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/${IMPL}")
        list(APPEND ENABLED_IMPLEMENTATIONS "${IMPL}")
      endif()
    endforeach()
  endif()

  list(REMOVE_DUPLICATES ENABLED_IMPLEMENTATIONS)
  list(REMOVE_DUPLICATES DISABLED_IMPLEMENTATIONS)
  set(${OUT_ENABLED_IMPLEMENTATIONS}
      "${ENABLED_IMPLEMENTATIONS}"
      PARENT_SCOPE)
  set(${OUT_DISABLED_IMPLEMENTATIONS}
      "${DISABLED_IMPLEMENTATIONS}"
      PARENT_SCOPE)
endfunction()

# Function to configure tests
function(add_tests test_flag exec_target subdir)
  if(${test_flag})
    # Gather all source files under tests/<subdir>
    file(GLOB_RECURSE src_files "${TEST_DIR}/${subdir}/*.cpp"
         "${TEST_DIR}/${subdir}/*.cxx" "${TEST_DIR}/${subdir}/*.cc")
    target_sources(${exec_target} PRIVATE ${src_files})
    list(APPEND TEST_EXECUTABLES ${exec_target})
    set(TEST_EXECUTABLES
        "${TEST_EXECUTABLES}"
        PARENT_SCOPE)
  endif()
endfunction()

# ============================================================================
# Function: setup_implementation - NAME:       implementation sub‐directory name
# (e.g. “mpi”) - PROJ_NAME:  project base name - BASE_DIR:   root source
# directory - TESTS:      list of test executables to link against
# ============================================================================
function(setup_implementation)
  # parse named args: NAME, PROJ_NAME, BASE_DIR; multi‐value: TESTS
  cmake_parse_arguments(SETUP "" # no plain options
                        "NAME;PROJ_NAME;BASE_DIR" "TESTS" ${ARGN})

  # skip if impl dir doesn't exist
  set(IMP_DIR "${SETUP_BASE_DIR}/${SETUP_NAME}")
  if(NOT EXISTS "${IMP_DIR}")
    return()
  endif()
  message(STATUS "  -- ${SETUP_NAME}")

  # Collect sources and create library: STATIC if implementation has .cpp files,
  # otherwise INTERFACE.
  file(GLOB_RECURSE IMPL_CPP_SOURCES "${IMP_DIR}/src/*.cpp")
  file(GLOB_RECURSE IMPL_SOURCES "${IMP_DIR}/include/*.h"
       "${IMP_DIR}/include/*.hpp" "${IMP_DIR}/src/*.cpp")

  set(LIB_NAME "${SETUP_PROJ_NAME}_${SETUP_NAME}")
  if(IMPL_CPP_SOURCES)
    add_library(${LIB_NAME} STATIC ${IMPL_SOURCES})
    set(LIB_LINK_SCOPE PUBLIC)
  else()
    add_library(${LIB_NAME} INTERFACE)
    target_sources(${LIB_NAME} INTERFACE ${IMPL_SOURCES})
    set(LIB_LINK_SCOPE INTERFACE)
  endif()

  # link core module
  target_link_libraries(${LIB_NAME} ${LIB_LINK_SCOPE} core_module_lib)

  # and link into each enabled test executable
  foreach(test_exec ${SETUP_TESTS})
    target_link_libraries(${test_exec} PUBLIC ${LIB_NAME})
  endforeach()
endfunction()

# Function to configure each subproject
function(ppc_configure_subproject SUBDIR)
  set(SETTINGS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/settings.json")

  # Keep per-task settings/id macros available even when the task is skipped.
  # Some tests reference settings from other tasks.
  add_compile_definitions(PPC_SETTINGS_${SUBDIR}="${SETTINGS_PATH}"
                          PPC_ID_${SUBDIR}="${SUBDIR}")

  ppc_collect_implementations_from_settings(
    "${SETTINGS_PATH}" "${SUBDIR}" ENABLED_IMPLEMENTATIONS
    DISABLED_IMPLEMENTATIONS)
  if(DISABLED_IMPLEMENTATIONS)
    list(JOIN DISABLED_IMPLEMENTATIONS ", " DISABLED_IMPLEMENTATIONS_STR)
    message(
      STATUS
        "${SUBDIR} (skipped: disabled in settings.json -> ${DISABLED_IMPLEMENTATIONS_STR})"
    )
    return()
  endif()

  if(NOT ENABLED_IMPLEMENTATIONS)
    message(
      STATUS "${SUBDIR} (skipped: no enabled implementations in settings.json)")
    return()
  endif()

  # Switch project context to the subproject
  project(${SUBDIR})

  # Directory with tests and list of test executables (populated by
  # setup_implementation)
  set(TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/tests")
  set(TEST_EXECUTABLES "")

  # Register functional and performance test runners
  add_tests(USE_FUNC_TESTS ${FUNC_TEST_EXEC} functional)
  add_tests(USE_PERF_TESTS ${PERF_TEST_EXEC} performance)

  message(STATUS "${SUBDIR}")

  # List of implementations to configure
  foreach(IMPL IN LISTS ENABLED_IMPLEMENTATIONS)
    setup_implementation(
      NAME
      ${IMPL}
      PROJ_NAME
      ${SUBDIR}
      TESTS
      "${TEST_EXECUTABLES}"
      BASE_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}")
  endforeach()
endfunction()
