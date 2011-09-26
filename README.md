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

-- standard factors
f = opengm.factors

-- define variables
variables = {'car', 'person', 'building', 'street', 'vehicle'}

-- define factors
factors = {-- unary factors (prior probabilities of each class):
           {f.prior(0.9),  {'car'}},
           {f.prior(0.01), {'person'}},
           {f.prior(0.7),  {'building'}},
           {f.prior(0.8),  {'street'}},
           {f.prior(0.4),  {'vehicle'}},
           -- Potts factors (joint probabilities):
           {f.band(0),     {'car',      'person'}},
           {f.band(0),     {'person',   'building'}},
           {f.band(0),     {'building', 'street'}},
           {f.band(0),     {'car',      'building'}},
           {f.band(0),     {'building', 'vehicle'}},
           {f.band(0),     {'street',   'vehicle'}},
           {f.band(0),     {'person',   'vehicle'}},
           {f.bimplies(1), {'car',      'vehicle'}}}

-- create graph
g = opengm.Graph(variables, factors)

-- optimize graph
g:optimize{method='a*', verbose=true}

-- print graph
print(g)
```

Running the script above outputs:

```
<opengm> optimizing... 
step 1: E=3.99758, c=0
step 2: E=3.63212, c=2.19722
step 3: E=3.63212, c=2.19722
<opengm.Graph>
  + nb of variables: 4
  + nb of factors: 6
  + graph is acyclic
  + current (optimized) variable states: 
    - car [1]
    - person [0]
    - building [0]
    - street [0]
    - vehicle [1]
```
