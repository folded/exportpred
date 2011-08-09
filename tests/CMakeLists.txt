include_directories("${exportpred_SOURCE_DIR}/include")

PROJECT(ghmm_tests)

SET(GTEST_DIR "${exportpred_SOURCE_DIR}/external/gtest-1.6.0")

include_directories(${GTEST_DIR})
include_directories(${GTEST_DIR}/include)

function(cxx_test_with_flags name cxx_flags libs)
  add_executable(${name} ${ARGN})
  set_target_properties(${name}
    PROPERTIES
    COMPILE_FLAGS "${cxx_flags}")
  # To support mixing linking in static and dynamic libraries, link each
  # library in with an extra call to target_link_libraries.
  foreach (lib "${libs}")
    target_link_libraries(${name} ${lib})
  endforeach()
  add_test(${name} ${name})
endfunction()

# cxx_test(name libs srcs...)
#
# creates a named test target that depends on the given libs and is
# built from the given source files.  Unlike cxx_test_with_flags,
# test/name.cc is already implicitly included in the source file list.
function(cxx_test name libs)
  cxx_test_with_flags("${name}" "${CMAKE_CXX_FLAGS}" "${libs}" "${name}.cpp" ${ARGN})
endfunction()
