
-- load opengm
require 'opengm'

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
print(g)
--g:show{}
