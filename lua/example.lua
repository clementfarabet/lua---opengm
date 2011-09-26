
-- load opengm
require 'opengm'

-- args
op = xlua.OptionParser('%prog [options]')
op:option{'-d', '--display', action='store_true', dest='display',
          help='display optimized graph (energies + states)'}
op:option{'-m', '--method', action='store', dest='method',
          help='optimization method: a* | bp | trbp | lf | icm', default='a*'}
opt = op:parse()

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
g:optimize{method=opt.method, verbose=true}

-- show graph
if opt.display then
   g:show{}
else
   print(g)
end
