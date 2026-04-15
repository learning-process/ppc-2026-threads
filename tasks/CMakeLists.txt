project(parallel_programming_course LANGUAGES C CXX)

message(STATUS "Student's tasks")

# Test runner executables
set(FUNC_TEST_EXEC ppc_func_tests)
set(PERF_TEST_EXEC ppc_perf_tests)

# ——— Include helper scripts ——————————————————————————————————————
include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

# ——— Initialize test executables —————————————————————————————————————
ppc_add_test(${FUNC_TEST_EXEC} common/runners/functional.cpp USE_FUNC_TESTS)
ppc_add_test(${PERF_TEST_EXEC} common/runners/performance.cpp USE_PERF_TESTS)

# ——— List of implementations ————————————————————————————————————————
set(PPC_IMPLEMENTATIONS "all;mpi;omp;seq;stl;tbb" CACHE STRING "Implementations to build (semicolon-separated)")

# ——— Configure ONLY MY subproject ——————————————————————————————————
# Мы закомментировали автоматический поиск всех папок (GLOB)
# и прописали вашу папку вручную.

set(MY_TASK "votincev_d_radixmerge_sort")

if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${MY_TASK}")
    message(STATUS "Configuring only: ${MY_TASK}")
    ppc_configure_subproject(${MY_TASK})
else()
    message(FATAL_ERROR "Directory ${MY_TASK} not found!")
endif()