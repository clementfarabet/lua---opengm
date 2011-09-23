----------------------------------------------------------------------
--
-- Copyright (c) 2011 Clement Farabet
-- 
-- Permission is hereby granted, free of charge, to any person obtaining
-- a copy of this software and associated documentation files (the
-- "Software"), to deal in the Software without restriction, including
-- without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to
-- permit persons to whom the Software is furnished to do so, subject to
-- the following conditions:
-- 
-- The above copyright notice and this permission notice shall be
-- included in all copies or substantial portions of the Software.
-- 
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
-- NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
-- LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
-- OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
-- WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
-- 
----------------------------------------------------------------------
-- description:
--     opengm - a wrapper around the OpenGM library.
--              this wrapper allows easy creation of graphs, and
--              efficient inference.
--
-- history: 
--     September 23, 2011, 1:43AM - first skeleton - Clement Farabet
----------------------------------------------------------------------

require 'torch'
require 'image'

-- create global opengm table:
opengm = {}

-- c lib:
require 'liblopengm'

-- help stirngs:
local help = {
   new = [[
This function creates a new graphical model. Here is a simple
example. Say we want to create a graph with 4 variables
x1, x2, x3 and x4. The first 3 variables are binary, and the
fourth is tertiary. We first define a state table:

> states = {2, 2, 2, 3}

Then we want to define 6 factors, 4 unary factors, and 2 Potts
factors. The table goes like this:

> factors = {{1}, {2}, {3}, {4}, {1,2}, {1,4}}

Each entry of this table defines a new factor, which takes for
arguments the variables listed. The 5th factor of the list is
noted f(x1,x2).

Lastly, we need to specify which all possible energies for each
factor. Energies are real-valued, an energy of 0 represents a
likely configuration, whereas a large energy represents an 
unlikely configuration.

> energies = {{0,0}, {0,100}, {100,0}, {0,40}, {0,100,0,0}, {100,100,0,0}}

In that example, x2 is likely to be 0 (low energy) and very unlikely
to be 1 (high energy). On the other hand, x1 is agnostic, having a
low energy for both 0 and 1.
For multi-term factors, like f(x1,x2), the energies are listed in the
same order they would be on a boolean table:

x1 | x2 | f(x1,x2)
 0    0    0
 0    1    100
 1    0    0
 1    1    0

Finally, the graph can be obtained like this:
> graph = opengm.Graph(states, factors, energies)    ]]
}

----------------------------------------------------------------------
-- visualize a graph
--
local function showgraph(graph, ...)
   -- usage
   local _, save, show = xlua.unpack(
      {...},
      'opengm.Graph:show()', 'render the graph using luagraph/graphviz',
      {arg='save', type='string', help='save to file (png)'},
      {arg='show', type='boolean', help='show', default=true}
   )

   -- simple print
   print(self)

   -- get structure from graph
   local states = graph.states
   local factors = graph.factors
   local energies = graph.energies
   local names = graph.names or {}

   -- then try to load luagraph package
   gr = xrequire 'graph'
   if not gr then print('cant show graph, luagraph required (luarocks install luagraph)') return end

   -- local symbols
   local node, edge, subgraph, cluster, digraph, strictdigraph =
      gr.node, gr.edge, gr.subgraph, gr.cluster, gr.digraph, gr.strictdigraph

   -- build graph
   local gargs = {'G',
                  compound = '1',
                  rankdir = 'LR',
                  size='1000,1000',
                  comment = 'Graphical Model'}

   -- names ?
   for i in ipairs(states) do
      names[i] = names[i] or ('x' .. i)
   end

   -- insert variable nodes
   for i,states in ipairs(states) do
      table.insert(gargs, node{'v'..i, label=names[i]})
   end

   -- insert factor nodes and their edges
   for i,factor in ipairs(factors) do
      local str = 'f('..names[factor[1]]
      if factor[2] then str = str .. ',' .. names[factor[2]] end
      str = str .. ')'
      table.insert(gargs, node{'f'..i, label=str, shape='square'})
      for _,k in ipairs(factor) do
         table.insert(gargs, edge{'f'..i, 'v'..k, color='red'})
      end
   end

   -- render graph
   local g = strictdigraph(gargs)
   if show then
      g:showdotty()
   end
   if save then
      g:layout() g:render('png', save)
   end
   g:close()
end

----------------------------------------------------------------------
-- a nice constructor to build graphs
--
local function newgraph(self,...)
   -- usage
   local _, states, factors, energies, names = xlua.unpack(
      {...},
      'opengm.Graph', help.new,
      {arg='states', type='table', help='a list of possible states, for each variable node', req=true},
      {arg='factors', type='table', help='a list of all factors\' arguments', req=true},
      {arg='energies', type='table', help='a list of all possible energies, for each factor', req=true},
      {arg='names', type='table', help='a list of names, one per variable'}
   )

   -- create graph, straight from arguments
   local graph = opengm.Graph.new(states, factors, energies)

   -- create larger container that retains structure tables as well
   local complete = {graph=graph, states=states, factors=factors, 
                     energies=energies, names=names}
   complete.show = function(self,...) showgraph(self,...) end
   complete.optimize = function(self,...) self.graph:optimize(...) end
   complete.variables = function(self,...) self.states:states(...) end
   setmetatable(complete, {__tostring = function(self) return tostring(self.graph) end})

   -- return complete container
   return complete
end
setmetatable(opengm.Graph, {__call = newgraph})
