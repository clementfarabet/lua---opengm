function s = opengm_info(gm)
%
s.number_of_variables = opengm_number_of_variables(gm);
s.number_of_factors = opengm_number_of_factors(gm);

s.numbers_of_states = unique(gm.numbers_of_states);
s.numbers_of_states_abundance = zeros(numel(s.numbers_of_states), 1);
for j = 1:length(s.numbers_of_states)
    s.numbers_of_states_abundance(j) = nnz(gm.numbers_of_states == s.numbers_of_states(j));
end

factor_orders = zeros(s.number_of_factors, 1);
for j = 1:s.number_of_factors
    vi = opengm_factor(gm, j);
    factor_orders(j) = length(vi);
end
s.factor_orders = unique(factor_orders);
s.factor_orders_abundance = zeros(numel(s.factor_orders), 1);
for j = 1:length(s.factor_orders)
    s.factor_orders_abundance(j) = nnz(factor_orders == s.factor_orders(j));
end

fov = opengm_factors_of_variables(gm);
factors_per_variable = zeros(s.number_of_variables, 1);
for j = 1:s.number_of_variables
    factors_per_variable(j) = numel(fov{j});
end
s.factors_per_variable = unique(factors_per_variable);
s.factors_per_variable_abundance = zeros(numel(s.factors_per_variable), 1);
for j = 1:length(s.factors_per_variable)
    s.factors_per_variable_abundance(j) = nnz(factors_per_variable == s.factors_per_variable(j));
end
%
end