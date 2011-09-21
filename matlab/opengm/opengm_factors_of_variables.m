function fov = opengm_factors_of_variables(gm)
%
fov = cell(opengm_number_of_variables(gm), 1);
for j = 1:opengm_number_of_factors(gm)
    variable_indices = opengm_factor(gm, j);
    for k = 1:length(variable_indices)
        vi = variable_indices(k);
        fov{vi} = [fov{vi} j];
    end
end
%
end