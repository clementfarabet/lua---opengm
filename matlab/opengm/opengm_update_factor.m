function opengm_update_factor(gm, index, variable_indices, table)
%
if index > gm.variable_indices.size()
    error('index exceeds number of factors.');
end
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
gm.variable_indices.setElementAt(variable_indices(:), index-1);
gm.tables.setElementAt(table, index-1);
%
end