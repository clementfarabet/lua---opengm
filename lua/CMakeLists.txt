
# find Torch7
string (REGEX REPLACE "(.*)lib/luarocks/rocks.*" "\\1" TORCH_PREFIX "${CMAKE_INSTALL_PREFIX}" )
message (STATUS "Found Torch7, installed in: " ${TORCH_PREFIX})
find_package (Torch REQUIRED)

include_directories (${TORCH_INCLUDE_DIR})
add_library (lopengm SHARED lopengm.cxx ${headers})
link_directories (${TORCH_LIBRARY_DIR})
target_link_libraries (lopengm ${TORCH_LIBRARIES})

install_files(/lua/opengm init.lua)
install_files(/lua/opengm factors.lua)
install_targets(/lib lopengm)
