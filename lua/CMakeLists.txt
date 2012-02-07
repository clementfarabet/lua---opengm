
# find Torch7
find_package (Torch REQUIRED)

SET(src lopengm.cxx ${headers})
SET(luasrc)
ADD_TORCH_PACKAGE(lopengm "${src}" "${luasrc}" "Graphs")
TARGET_LINK_LIBRARIES(lopengm luaT TH)

SET(src)
SET(luasrc init.lua factors.lua help.lua)
ADD_TORCH_PACKAGE(opengm "${src}" "${luasrc}" "Graphs")
