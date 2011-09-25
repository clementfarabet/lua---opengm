
-- factors is table that defines standard factor functions
opengm.factors = {}

-- epsilon: the minimum probabilty that we want to be able to 
-- represent (a proba of 0 is an infinite energy)
local eps = 1e-20

-- neutral energy: no prior
opengm.factors.noprior = function()
                            return function()
                                      return 0
                                   end
                         end

-- standard unary energy/prior
opengm.factors.prior = function(proba)
                          proba = math.max(proba,eps)
                          return function(x)
                                    if x == 1 then return -math.log(proba)
                                    else return -math.log(1-proba) end
                                 end
                       end

-- standard joint probability btwn two variables
opengm.factors.joint = function(proba)
                          proba = math.max(proba,eps)
                          return function(x1,x2)
                                    if x1 == x2 then return -math.log(proba)
                                    else return -math.log(1-proba) end
                                 end
                       end
