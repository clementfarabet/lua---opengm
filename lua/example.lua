
-- load opengm
require 'opengm'

-- args
op = xlua.OptionParser('%prog [options]')
op:option{'-d', '--display', action='store_true', dest='display',
          help='display optimized graph (energies + states)'}
op:option{'-m', '--method', action='store', dest='method',
          help='optimization method: a* | bp | trbp | lf | icm', default='bp'}
opt = op:parse()

-- standard factors
f = opengm.factors

-- define variables
variables = {'car', 'person', 'building', 'street'}

-- define factors
factors = {-- unary factors (prior probabilities of each class):
           {f.prior(0.01), {1}},
           {f.prior(0.01), {2}},
           {f.prior(0.7),  {3}},
           {f.prior(0.8),  {4}},
           -- Potts factors (joint probabilities):
           {f.joint(0),   {1, 2}},
           {f.joint(0),   {3, 4}}}

-- create graph
g = opengm.Graph(variables, factors)

-- optimize graph
g:optimize{method=opt.method, verbose=true}

-- show graph
if opt.display then
   g:show{}
else
   print(g)
end
