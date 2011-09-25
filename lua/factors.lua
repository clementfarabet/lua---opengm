
-- factors is table that defines standard factor functions
opengm.factors = {}

-- epsilon: the minimum probabilty that we want to be able to 
-- represent (a proba of 0 is an infinite energy)
local eps = 1e-20

-- neutral energy: no prior
opengm.factors.noprior = function()
                            return function(x)
                                      return -math.log(1/2)
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

-- probabilty of x1 and x2 both active:
opengm.factors.band = function(proba)
                        proba = math.max(proba,eps)
                        return function(x1,x2)
                                  if ((x1 == x2) and (x1 == 1)) then return -math.log(proba)
                                  else return -math.log((1-proba)/3) end
                               end
                     end

-- probabilty of neither x1 nor x2 active:
opengm.factors.bnor = function(proba)
                         proba = math.max(proba,eps)
                         return function(x1,x2)
                                   if ((x1 == x2) and (x1 == 0)) then return -math.log(proba)
                                   else return -math.log((1-proba)/3) end
                                end
                      end

-- probabilty of either x1 or x2, but not both:
opengm.factors.bneq = function(proba)
                         proba = math.max(proba,eps)
                         return function(x1,x2)
                                   if (x1 ~= x2) then return -math.log(proba/2)
                                   else return -math.log((1-proba)/2) end
                                end
                      end

-- probabilty of x1 == x2:
opengm.factors.beq = function(proba)
                        proba = math.max(proba,eps)
                        return function(x1,x2)
                                  if (x1 == x2) then return -math.log(proba/2)
                                  else return -math.log((1-proba)/2) end
                               end
                     end

-- probabilty of x1 => x2 (implication):
opengm.factors.bimplies = function(proba)
                             proba = math.max(proba,eps)
                             return function(x1,x2)
                                       if (x1 == 0) or (x2 == 1) then return -math.log(proba/3)
                                       else return -math.log((1-proba)) end
                                    end
                          end
