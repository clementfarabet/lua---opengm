function index = opengm_add_factor(gm, variable_indices, table)
%
if length(variable_indices) == 1
    table = table(:);
else
    if ndims(table) ~= length(variable_indices)
        error('ndims(table) ~= length(variable_indices).');
    end
end
for j = 1:length(variable_indices)
    if variable_indices(j) > length(gm.numbers_of_states)
        error('variable index %d at index %d is too large.', variable_indices(j), j);
    end
    if size(table, j) ~= gm.numbers_of_states(variable_indices(j))
        error('table does not have the right shape.');
    end
end
gm.variable_indices.addElement(variable_indices(:));
gm.tables.addElement(table);
index = gm.tables.size();
%
end