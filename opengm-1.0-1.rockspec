
package = "opengm"
version = "1.0-1"

source = {
   url = "opengm-1.0-1.tgz"
}

description = {
   summary = "Lua bindings for OpenGM, a Graphical Model package.",
   detailed = [[
            Allows easy graphical model descriptions in Lua, and 
            efficient inference using OpenGM's A* and Belief Propagation
            implementations.
   ]],
   homepage = "http://www.andres.sc/opengm",
   license = "MIT/X11"
}

dependencies = {
   "lua >= 5.1",
   "torch",
   "xlua",
   "luagraph"
}

build = {
   type = "cmake",
   variables = {
      CMAKE_INSTALL_PREFIX = "$(PREFIX)"
   }
}
