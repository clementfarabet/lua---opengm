FILE(GLOB headers ../../include/opengm/* ../../include/opengm/inference/*)

add_executable(test-opengm test.cxx ${headers})
add_test(test-opengm test-opengm)

if(Boost_FOUND)
  add_executable(test-graphcut test-graphcut.cxx ${headers})
  add_test(test-graphcut test-graphcut)
endif()
