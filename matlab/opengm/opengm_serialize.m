function ser = opengm_serialize(gm)
%
% Warning: variable indices in the serialization 
% start from zero
%

% pre-compute size
size = 2 + opengm_number_of_variables(gm);
for j = 1:opengm_number_of_factors(gm)
    [vi, t] = opengm_factor(gm, j);
    size = size + 1 + numel(vi) + numel(t);
end
ser = zeros(size, 1);

% write output
ser(1) = opengm_number_of_variables(gm);
ser(2) = opengm_number_of_factors(gm);
n = 2;
for j = 1:opengm_number_of_variables(gm)
    n = n + 1;
    ser(n) = gm.numbers_of_states(j);
end
for j = 1:opengm_number_of_factors(gm)
    [vi, t] = opengm_factor(gm, j);
    n = n + 1;
    ser(n) = length(vi);
    for k = 1:length(vi)
        n = n + 1;
        ser(n) = vi(k)-1;
    end
    ser( n+1 : n+numel(t) ) = t;
    n = n + numel(t);
end

end