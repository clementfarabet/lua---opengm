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
--     September 24, 2011, 9:07PM - first usable state - Clement Farabet
--     September 23, 2011, 1:43AM - first skeleton - Clement Farabet
----------------------------------------------------------------------

require 'torch'
require 'xlua'

-- create global opengm table:
opengm = {}

-- c lib:
require 'liblopengm'

-- submodules:
dofile(sys.concat(sys.fpath(),'factors.lua'))

-- help strings:
local help = {
   new = [[
This function creates a new graphical model. Here is a simple
example. Say we want to create a graph with 4 variables
x1, x2, x3 and x4. The first 3 variables are binary, and the
fourth is tertiary. We define the variables like this:

> variables = {{'x1',2}, {'x2',2}, {'x3',2}, {'x4',3}}

(note: binary variables are the default, so {'x1',2} can be 
 replaced by 'x1', for simplicity)

Then we need to define the factors of the graph.
Say we want to define 6 factors, 4 unary factors, and 2 Potts
factors. The table goes like this:

> factors = {{func1, {1}}, {func2, {2}}, {func3, {3}}, {func4, {4}}, 
             {func5, {1,2}}, {func6, {1,4}}}

Each entry of this table defines a new factor, which takes for
arguments the variables listed. The first part of each tuple
is a lua function that was defined apriori.

Such a function is simple to define: it must return an energy for each 
configuration of its inputs. For instance, func5 could be defined as:

func5 = function(x1,x2) 
           if x1 == 1 and x2 == 0 then return 1e3
           else return 0 end
        end

That function assigns a high cost to (x1,x2) == (1,0), and 0-costs to
all other configurations.

(note: the table opengm.factors provides a couple of standard factor
 functions)

The graph can be obtained like this:
> graph = opengm.Graph(variables, factors)    ]],

optimize = [[
This function optimizes the current graphical model (inference). Given
a graph g, a call to g:optimize{} fills the states vector with the
optimal states obtained. These can be retrieved like this:
states = g:getstates().

Several inference algorithms are available:
- A* (a*)
- Belief Propagation (bp)
- Tree-Reweighted Belief Propagation (trbp)
- Iterated Conditional Modes (icm)
- Lazy Flipper (lf) 

Once optimized, printing a graph reveals its optimal states:
> g:optimize{}
> print(g)
<opengm.Graph>
  + nb of variables: 4
  + nb of factors: 6
  + graph is acyclic
  + current (optimized) variable states: 
    - x1 [0]
    - x2 [0]
    - x3 [0]
    - x4 [1]
]]
}

----------------------------------------------------------------------
-- visualize a graph
--
local function showgraph(graph, ...)
   -- usage
   local _, save, show = xlua.unpack(
      {...},
      'opengm.Graph:show', 'render the graph using luagraph/graphviz',
      {arg='save', type='string', help='save to file (png)'},
      {arg='show', type='boolean', help='show', default=true}
   )

   -- simple print
   print(graph)

   -- get structure from graph
   local states = graph.states
   local factors = graph.factors
   local energies = graph.energies
   local functions = graph.functions
   local names = graph.variables or {}

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

   -- are states available ?
   local optstates = g:getstates() or {}

   -- insert variable nodes
   for i,states in ipairs(states) do
      local name = names[i]
      if optstates[i] then name = name .. ' [=' .. optstates[i] .. ']' end
      table.insert(gargs, node{'v'..i, label=name})
   end

   -- insert factor nodes and their edges
   for i,factor in ipairs(factors) do
      local str = 'f('..names[factor[1]]
      if factor[2] then str = str .. ',' .. names[factor[2]] end
      str = str .. ')'
      if optstates[1] then
         local func = functions[i][1]
         local args = functions[i][2]
         local val
         if #factor == 1 then
            val = func(optstates[args[1]])
         elseif #factor == 2 then
            val = func(optstates[args[1]], optstates[args[2]])
         end
         str = str .. ' [=' .. tostring(val) .. ']'
      end
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
-- optimize a graph
--
local function optimize(graph, ...)
   -- usage
   local _, method, iterations, verbose = xlua.unpack(
      {...},
      'opengm.Graph:optimize', help.optimize,
      {arg='method', type='string', help='optimization method: a* | bp | trbp | lf | icm', default='bp'},
      {arg='verbose', type='boolean', help='verbose'},
      {arg='iterations', type='number', help='max number of iterations (bp & trbp)', default=100},
      {arg='damping', type='number', help='damping coefficient (bp & trbp)', default=0.0}
   )

   -- optimize
   graph.graph:optimize(method, verbose, iterations, damping)
end

----------------------------------------------------------------------
-- optimize a graph
--
local function graphtostring(graph)
   local str = tostring(graph.graph)
   local optstates = graph:getstates()
   if optstates and graph.variables then
      str = str .. '\n  + current (optimized) variable states: '
      for i,name in ipairs(graph.variables) do
         str = str .. '\n    - ' .. name .. ' [' .. optstates[i] .. ']'
      end
   end
   return str
end

----------------------------------------------------------------------
-- a nice constructor to build graphs
--
local function newgraph(self,...)
   -- usage
   local _, variables, factors = xlua.unpack(
      {...},
      'opengm.Graph', help.new,
      {arg='variables', type='table', help='a list of variables: {{varname1,nbstates}, {varname2,nbstates}, ...}', req=true},
      {arg='factors', type='table', help='a list of factors: {{factor1,arg1[,arg2]}, {factor2,arg1[,arg2]}, ...}', req=true}
   )

   -- variables may contain their possible states
   local states = {}
   for i,var in ipairs(variables) do
      if type(var) == 'table' then
         variables[i] = var[1]
         states[i] = var[2]
      else
         -- default nb of states: binary
         states[i] = 2
      end
   end

   -- evaluate energies for all factors
   local energies = {}
   local functions = factors
   local factors = {}
   for i,factor in ipairs(functions) do
      local en = factor[1]
      local args = factor[2]
      energies[i] = {}
      local tab = energies[i]
      if #args == 1 then -- unary factor
         -- eval energy for all states
         for s = 0,states[args[1]]-1 do
            table.insert(tab, en(s))
         end
      elseif #args == 2 then -- 2-term factor
         -- eval energy for all states
         for s1 = 0,states[args[1]]-1 do
            for s2 = 0,states[args[2]]-1 do
               table.insert(tab, en(s1,s2))
            end
         end
      end
      factors[i] = args
   end

   -- create graph, straight from arguments
   local graph = opengm.Graph.new(states, factors, energies)

   -- create larger container that retains structure tables as well
   local complete = {graph=graph, states=states, factors=factors, 
                     functions = functions, energies=energies, variables=variables}
   complete.show = function(self,...) showgraph(self,...) end
   complete.optimize = function(self,...) return optimize(self,...) end
   complete.getstates = function(self,...) return self.graph:states(...) end
   setmetatable(complete, {__tostring = function(self) return graphtostring(self) end})

   -- return complete container
   return complete
end
setmetatable(opengm.Graph, {__call = newgraph})
