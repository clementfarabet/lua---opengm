# LOpenGM: Lua bindings for OpenGM

[OpenGM](http://www.andres.sc/opengm) is a C++ library for graphical 
modeling, and inference. The Lua
bindings provide a simple way of describing graphs, from Lua, and then
optimizing them with OpenGM.

## License

LOpenGM Copyright (c) 2011 Clement Farabet (Lua Bindings)

OpenGM  Copyright (c) 2010 by Bjoern Andres and Joerg Hendrik Kappes.

This software was developed by Bjoern Andres and Joerg Hendrik Kappes.
Enquiries shall be directed to:

bjoern.andres@iwr.uni-heidelberg.de, kappes@math.uni-heidelberg.de

All advertising materials mentioning features or use of this software must
display the following acknowledgement: ``This product includes the OpenGM
library developed by Bjoern Andres and Joerg Hendrik Kappes. Please direct
enquiries concerning OpenGM to bjoern.andres@iwr.uni-heidelberg.de,
kappes@math.uni-heidelberg.de''.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- All advertising materials mentioning features or use of this software must
  display the following acknowledgement: ``This product includes the OpenGM
  library developed by Bjoern Andres and Joerg Hendrik Kappes. Please direct
  enquiries concerning OpenGM to bjoern.andres@iwr.uni-heidelberg.de,
  kappes@math.uni-heidelberg.de''.
- The names of the authors must not be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Install dependencies 

1/ third-party libraries:

On Linux (Ubuntu > 9.04):

``` sh
$ apt-get install gcc g++ git libreadline5-dev cmake graphviz
```

On Mac OS (Leopard, or more), using [Homebrew](http://mxcl.github.com/homebrew/):

``` sh
$ brew install git readline cmake graphviz
```

2/ Lua 5.1 + Luarocks + xLua:

``` sh
$ git clone https://github.com/clementfarabet/lua4torch
$ cd lua4torch
$ make install PREFIX=/usr/local
```

3/ opengm:

clone this repo and then:

``` sh
$ luarocks make
```

or, without the repo:

``` sh
$ luarocks install opengm
```

(for info: this will first install Torch7, which is used to efficiently
represent N-dim arrays, for variables, factors, ...)

## Use the library

### API, in very short:

Load/start up package:

``` lua
require 'opengm'
```

Construct a graph:

``` lua
g = opengm.Graph(...)
```

Optimize a graph:

``` lua
g:optimize{}
```

Display a graph, using Graphviz:

``` lua
g:show{}
```

A simple complete example:

```lua
-- load opengm
require 'opengm'

-- args
op = xlua.OptionParser('%prog [options]')
op:option{'-dp', '--display', action='store_true', dest='display',
          help='display optimized graph (energies + states)'}
opt = op:parse()

-- standard factors
f = opengm.factors

-- define variables
variables = {'car', 'person', 'building', 'street'}

-- define factors
factors = {-- unary factors (prior probabilities of each class):
           {f.prior(0.1), {1}},
           {f.prior(0.1), {2}},
           {f.prior(0.51),{3}},
           {f.prior(0.6), {4}},
           -- Potts factors (joint probabilities):
           {f.joint(0),   {1, 2}},
           {f.joint(0),   {3, 4}}}

-- create graph
g = opengm.Graph(variables, factors)

-- optimize graph
g:optimize{verbose=true}

-- show graph
if opt.display then
   g:show{}
else
   print(g)
end
```