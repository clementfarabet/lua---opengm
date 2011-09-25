
-- help strings:
opengm.help = {
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
