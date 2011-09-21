FILE(GLOB headers ../../include/opengm/* ../../include/opengm/inference/*)

if(HDF5_FOUND)
  add_executable(tutorial1 tutorial1.cxx ${headers})
  target_link_libraries(tutorial1 ${HDF5_CORE_LIBRARY})

  add_executable(tutorial2 tutorial2.cxx ${headers})
  target_link_libraries(tutorial2 ${HDF5_CORE_LIBRARY})
endif()