function gm = opengm_deserialize(ser)
%
number_of_variables = ser(1);
number_of_potentials = ser(2);
numbers_of_states = zeros(number_of_variables, 1);
j = 2;
for k = 1:number_of_variables
    j = j + 1;
    numbers_of_states(k) = ser(j);
end
gm = opengm_create(numbers_of_states);
for k = 1:number_of_potentials
    j = j + 1;
    potential_number_of_variables = ser(j);
    potential_variable_indices = zeros(potential_number_of_variables, 1);
    for m = 1:potential_number_of_variables
        j = j + 1;
        potential_variable_indices(m) = ser(j) + 1; % MATLAB indexing
    end
    potential_table_dimension = numel(potential_variable_indices);
    potential_table_size = numbers_of_states(potential_variable_indices);
    if potential_table_dimension == 1
        potential_table = zeros(potential_table_size, 1);
    else
        potential_table = zeros(potential_table_size');
    end
    for m = 1:numel(potential_table)
        j = j + 1;
        potential_table(m) = ser(j);
    end
    opengm_add_factor(gm, potential_variable_indices, potential_table);
end
%
end