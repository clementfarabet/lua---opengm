function [vi, table] = opengm_factor(gm, index)
%
if index > gm.variable_indices.size()
    error('index out of bounds.');
end

vi = gm.variable_indices.get(index-1);
table = gm.tables.get(index-1);
%
end